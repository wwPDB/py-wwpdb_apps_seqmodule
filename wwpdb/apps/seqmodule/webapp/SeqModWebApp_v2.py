##
# File:  SeqModWebApp.py
# Date:  16-Dec-2009
# Updates:
#
# 14-Jan-2010  Refactor to SeqToolWebApp and remove remaining webserver dependencies.
# 01-Feb-2010  Revise json packaging.
# 10-Feb-2010  Menu edit features added.
# 14-Feb-2010  Save feature added to export results.
# 03-Mar-2010  Make url mapping with functor dictionary.
# 03-Mar-2010  Add polling methods -
# 12-Mar-2010  Provide top-level adjustment for rcsb path and example path details
# 09-Apr-2010  Use form provided file type and unrestrict file name convention.
# 20-Apr-2010  Renamed to SeqModWebApp and ported to seqmodule package
# 24-Apr-2010  Streamline data import operations with DataImporter() class
# 28-Apr-2010  Add classId and siteId
# 28-Apr-2010  Add new methods for new_sequence and new_taxid
#              Add status database updates
# 02-May-2010  Update SequenceDataExport -- Add file path handling for export files.
# 05-May-2010  Add identifiers to alignment template.
#              Fixed file upload isssue
# 27-Jul-2010  RPS: Added support for accommodating different ordering of sequence types as per user preferences
# 02-Aug-2010  RPS: splitting SeqModWebAppWorker's "save/done" processing into separate calls
# 12-Aug-2010  RPS: Updated to accommodate new strategy for Save vs. Back to WFM buttons.
# 22-Sep-2011  RPS: topPath now being derived from ConfigInfoData based on siteId value.
# 19-Dec-2011  RPS: __isFileUpload() updated to return False in cases where type of 'file' is types.UnicodeType.
# 10-Oct-2012  RPS: Now deriving path of sessions directory from ConfigInfoData.
#
# 20-Feb-2013  jdw Refactoring - update to current ConfigInfoData, common request model, common session model.
# 03-Feb-2013  jdw common file upload handling
# 07-Mar-2013  jdw move all of the child process management to DetachUtils()
# 22-Mar-2013  jdw add download option for assignments and updated models.
# 25-Mar-2013  jdw add reload and reload checking op
# 27-Mar-2013  jdw add data reload worker method implementation.
# 16-Sep-2013  jdw add services for entityreview  -
# 15-Nov-2013  jdw update status tracking option
# 12-Dec-2013  jdw emit warning message on export if sample/xyz conflicts are found
# 07-Feb-2014  jdw use Wftracking in the api.status
# 10-Feb-2014  jdw update closing/finish processing protocol.
# 20-Feb-2014  jdw fix seq_search_op overwrite -
# 22-Feb-2014  jdw update archive/workflow files only on completion
# 22-Mar-2014  jdw refactor edit forms operations --
# 20-Apr-2014  jdw update input export file copy operations -
# 28-Apr-2014  jdw retire DataImporter.updateData*  modules
# 14-May-2014  jdw assign any missing aligntag to the current entity group
# 22-May-2014  jdw new alignment protocol -- push entry details in summary update response
#  8-Jul-2014  jdw export group id list in the summary view return object
#  2-Dec-2014  jdw capture updates in selection and alignment id lists after alignment edits --
# 24-Aug-2017  zf  add '/service/sequence_editor/rerun_blast/start' & _reRunBlastOp for rerun blast search
# 19-Oct-2022  zf  add __getSummaryPageContent() method
#  6-Jan 2024  zf  add /service/sequence_editor/respond_form/seqbuilder and _respondSeqBuilderOp() method
##
#    WF Testing entry points -
#
#   /service/sequence_editor/new_session/wf?classID=AnnMod&identifier=D_1000000000&filesource=wf-archive&instance=W_000
#   /service/sequence_editor/new_session/wf?classID=AnnMod&identifier=D_1000000001&filesource=wf-archive&instance=W_000
##
"""
Sequence editor tool web request and response processing modules.

This software was developed as part of the World Wide Protein Data Bank
Common Deposition and Annotation System Project

Copyright (c) wwPDB

This software is provided under a Creative Commons Attribution 3.0 Unported
License described at http://creativecommons.org/licenses/by/3.0/.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.07"

try:
    import cPickle as pickle
except ImportError:
    import pickle as pickle
#
import os
import shutil
import sys
import time
import traceback

from wwpdb.utils.config.ConfigInfo import ConfigInfo

from wwpdb.utils.wf.dbapi.WfTracking import WfTracking

from wwpdb.apps.seqmodule.align.AlignmentDepictionTools import AlignmentDepictionTools
from wwpdb.apps.seqmodule.align.AlignmentFrontEndUIEditor import AlignmentFrontEndUIEditor

from wwpdb.apps.seqmodule.control.DataImporter import DataImporter
from wwpdb.apps.seqmodule.control.SequenceDataAssemble_v2 import SequenceDataAssemble
from wwpdb.apps.seqmodule.control.SummaryView_v2 import SummaryView
from wwpdb.apps.seqmodule.control.SummaryViewDepiction import SummaryViewDepiction

from wwpdb.apps.seqmodule.io.SequenceDataExport_v2 import SequenceDataExport

from wwpdb.apps.seqmodule.update.UpdatePolymerEntityPartitions import UpdatePolymerEntityPartitions
from wwpdb.apps.seqmodule.update.UpdatePolymerEntitySourceDetails import UpdatePolymerEntitySourceDetails
from wwpdb.apps.seqmodule.update.UpdatePolymerEntityReference_v2 import UpdatePolymerEntityReference

from wwpdb.apps.seqmodule.view3d.ModelViewer3D import ModelViewer3D

from wwpdb.apps.seqmodule.webapp.SeqModWebRequest import SeqModInputRequest, SeqModResponseContent

from wwpdb.utils.db.DBLoadUtil import DBLoadUtil
from wwpdb.utils.detach.DetachUtils import DetachUtils
from wwpdb.utils.session.WebUploadUtils import WebUploadUtils
from wwpdb.utils.session.WebDownloadUtils import WebDownloadUtils
from wwpdb.utils.session.UtilDataStore import UtilDataStore
from wwpdb.io.locator.PathInfo import PathInfo


class SeqModWebApp(object):
    """Handle request and response object processing for sequence editor tool application."""

    def __init__(self, parameterDict=None, verbose=False, log=sys.stderr, siteId="WWPDB_DEPLOY_TEST"):
        """
        Create an instance of `SeqModWebApp` to manage a sequence editor web request.

         :param `parameterDict`: dictionary storing parameter information from the web request.
             Storage model for GET and POST parameter data is a dictionary of lists.
         :param `verbose`:  boolean flag to activate verbose logging.
         :param `log`:      stream for logging.

        """
        if parameterDict is None:
            parameterDict = {}
        self.__verbose = verbose
        self.__lfh = log
        self.__debug = False
        self.__siteId = siteId
        self.__cI = ConfigInfo(self.__siteId)
        #
        # the <site-specificit...>/webapps  directory path --  containing htdocs & fcgi subdirs
        self.__topPath = self.__cI.get("SITE_WEB_APPS_TOP_PATH")

        # The path containing the sessions directory (i.e. <top_sessions_path>/sessions )
        self.__topSessionPath = self.__cI.get("SITE_WEB_APPS_TOP_SESSIONS_PATH")

        if isinstance(parameterDict, dict):
            self.__myParameterDict = parameterDict
        else:
            self.__myParameterDict = {}

        if self.__debug:
            self.__lfh.write("+SeqModWebApp.__init() - SERVICE REQUEST STARTING with input parameter dictionary \n")
            self.__lfh.write("%s" % ("".join(self.__dump())))

        self.__reqObj = SeqModInputRequest(self.__myParameterDict, verbose=self.__verbose, log=self.__lfh)
        self.__templatePath = os.path.join(self.__topPath, "htdocs", "seqmodule")
        self.__reqObj.setValue("TopSessionPath", self.__topSessionPath)
        self.__reqObj.setValue("TemplatePath", self.__templatePath)
        self.__reqObj.setValue("TopPath", self.__topPath)
        self.__reqObj.setValue("WWPDB_SITE_ID", self.__siteId)
        os.environ["WWPDB_SITE_ID"] = self.__siteId
        #
        #
        # Example data path details -- for internally stored example data ---
        # self.__examplePath      =os.path.join(self.__topPath,"scripts","wwpdb","apps","seqmodule","examples")
        # self.__examplePathData  =os.path.join(self.__examplePath,"rcsb-data")
        # self.__examplePathSeq   =os.path.join(self.__examplePath,"rcsb-sequence")
        #
        # RCSB Specific paths
        #
        # self.__reqObj.setValue("RcsbDataPath","/annotation")
        # self.__reqObj.setValue("RcsbReferenceSequencePath","/www-rcsb/supertool/blast/rcsb")
        # self.__reqObj.setValue("RcsbDataPathExample",self.__examplePathData)
        # self.__reqObj.setValue("RcsbReferenceSequencePathExample",self.__examplePathSeq)
        #
        #

    def doOp(self):
        """
        Execute request and package results in response dictionary.

        :returns:
             A dictionary containing response data for the input request.
             Minimally, the content of this dictionary will include the
             keys: CONTENT_TYPE and REQUEST_STRING.

        """
        #
        if self.__verbose:
            rqp = self.__reqObj.getRequestPath()
            self.__lfh.write("\n#########################################################################################################\n")
            self.__lfh.write("+SeqModWebApp.__doOP() - Beginning request:  %s\n" % rqp)
            self.__reqObj.printIt(ofh=self.__lfh)

        stw = SeqModWebAppWorker(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        # doOp returns a response content object --
        rC = stw.doOp()
        if self.__debug:
            self.__lfh.write("+SeqModWebApp.doOp() working method completed with response format %s\n" % self.__reqObj.getReturnFormat())
            if rC is not None:
                self.__lfh.write("%s" % ("".join(rC.dump())))
            else:
                self.__lfh.write("+SeqModWebApp.doOp() return object is empty\n")

        self.__lfh.write("+SeqModWebApp.doOp() request completed \n\n")
        return rC.get()

    def __dump(self):
        """Utility method to format the contents of the internal parameter dictionary
        containing data from the input web request.

        :returns:
            ``list`` of formatted text lines

        """
        retL = []
        retL.append("\n +SeqModWebAppInput.__dump()  request dictionary length = %d\n" % len(self.__myParameterDict))
        for k, vL in self.__myParameterDict.items():
            retL.append("  - Key:  %-35s" % k)
            for v in vL:
                retL.append(" ->  %s\n" % v)
        retL.append("    -------------------------------------------------------------\n")
        return retL


class SeqModWebAppWorker(object):
    def __init__(self, reqObj=None, verbose=False, log=sys.stderr):
        """
        Worker methods for the sequence editor application

        Performs URL - application mapping and application launching
        for sequence editor tool.

        All operations can be driven from this interface which can
        supplied with control information from web application request
        or from a testing application.


        """
        self.__verbose = verbose
        self.__lfh = log
        #
        self.__reqObj = reqObj
        self.__siteId = str(self.__reqObj.getValue("WWPDB_SITE_ID"))
        self.__doWfTracking = True
        self.__sObj = None
        self.__sessionPath = None
        #
        self.__appPathD = {
            "/service/environment/dump": "_dumpOp",
            "/service/sequence_editor/store_alignment/start": "_storeAlignmentStartOp",
            "/service/sequence_editor/store_alignment/check": "_storeAlignmentCheckOp",
            "/service/sequence_editor/load_summary": "_loadSummaryOp",
            "/service/sequence_editor/reload_summary": "_reloadSummaryOp",
            "/service/sequence_editor/save": "_saveOp",
            "/service/sequence_editor/close_unfinished": "_closeUnfinishedOp",
            "/service/sequence_editor/close_completed": "_closeCompletedOp",
            "/service/sequence_editor/save_partial_assignment": "_savePartialAssignmentOp",
            "/service/sequence_editor/remove_partial_assignment": "_removePartialAssignmentOp",
            "/service/sequence_editor/new_session": "_newSessionOp",
            "/service/sequence_editor/new_session/wf": "_newSessionWfOp",
            "/service/sequence_editor/edit": "_editOp",
            "/service/sequence_editor/global_edit": "_globalEditOp",
            "/service/sequence_editor/global_auth_seq_edit": "_globalAuthSeqEditOp",
            "/service/sequence_editor/global_edit_menu": "_globalEditMenuOp",
            "/service/sequence_editor/move": "_moveEditOp",
            "/service/sequence_editor/undo_edit": "_undoEditOp",
            "/service/sequence_editor/delete": "_deleteOp",
            "/service/sequence_editor/molviewer/jmol": "_launchJmolViewerOp",
            "/service/sequence_editor/molviewer/astexviewer": "_launchAstexViewerOp",
            #
            "/service/sequence_editor/load_data/start/rcsb": "_loadDataStartOp",
            "/service/sequence_editor/load_data/check/rcsb": "_loadDataCheckOp",
            "/service/sequence_editor/load_data/check": "_loadDataCheckOp",
            #
            "/service/sequence_editor/load_data/start/wf": "_loadDataWfStartOp",
            "/service/sequence_editor/load_data/check/wf": "_loadDataCheckOp",
            #
            "/service/sequence_editor/rerun_blast/start": "_reRunBlastOp",
            #
            "/service/sequence_editor/load_form/taxonomy": "_loadTaxonomyFormOp",
            "/service/sequence_editor/load_form/seqdbref": "_loadSeqDbRefFormOp",
            "/service/sequence_editor/load_form/entityreview": "_loadEntitySourceDetailsFormOp",
            "/service/sequence_editor/respond_form/taxonomy": "_respondTaxonomyFormOp",
            "/service/sequence_editor/respond_form/seqbuilder": "_respondSeqBuilderOp",
            "/service/sequence_editor/respond_form/seqdbref": "_respondSeqDbRefFormOp",
            "/service/sequence_editor/respond_form/entityreview": "_respondEntitySourceDetailsFormOp",
            "/service/sequence_editor/download": "_respondDownloadOp",
            "/service/sequence_editor/reload_data/start": "_reloadDataStartOp",
            "/service/sequence_editor/reload_data/check": "_loadDataCheckOp",
            #
            "/service/sequence_editor/align_view": "_alignViewOp",
            "/service/sequence_editor/align_view_frame": "_alignViewFrameOp",
            "/service/sequence_editor/polymer_linkage_table": "_getPolymerLinkageTableOp",
        }

    def doOp(self):
        #
        # Map path to operation -
        #
        try:
            reqPath = self.__reqObj.getRequestPath()
            if reqPath not in self.__appPathD:
                # bail out if operation is unknown -
                self.__reqObj.setReturnFormat(return_format="json")
                rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
                rC.setError(errMsg="Unknown operation")
            else:
                mth = getattr(self, self.__appPathD[reqPath], None)
                rC = mth()
            return rC
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.__reqObj.setReturnFormat(return_format="json")
            rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            rC.setError(errMsg="Operation failure")
            return rC

    def setLogHandle(self, log=sys.stderr):
        """Reset the stream for logging output."""
        try:
            self.__lfh = log
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            return False

    # ------------------------------------------------------------------------------------------------------------
    #      Top-level REST method resolvers ---
    #
    #
    def _dumpOp(self):
        self.__reqObj.setReturnFormat(return_format="html")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        oL = self.__reqObj.dump(format="html")
        rC.setHtmlList(oL)
        return rC

    def _alignViewOp(self):
        """ """
        return self.__alignView()

    def _alignViewFrameOp(self):
        """ """
        return self.__alignView(frameOpt=True)

    def _saveOp(self):
        return self.__saveSelection()

    def _closeCompletedOp(self):
        """RPS, 2010-Aug-02: created independent function for
        closing sequence editor app independent of saving selection
        RPS, 2010-Aug-12: renamed for explicit use for completion of seq edit process
        """
        return self.__closeSeqEditor(mode="completed")

    def _closeUnfinishedOp(self):
        """RPS, 2010-Aug-12: created independent function for
        closing sequence editor app independent of saving selection
        *BUT* with intention of returning to complete editing activity
        """
        return self.__closeSeqEditor(mode="unfinished")

    def _savePartialAssignmentOp(self):
        """Save partial assignment information into data archive directory"""
        self.__getSession()
        #
        # Remove previous saved partial assignment information before new saving operation
        #
        self.__remove_partial_assignment_file()
        #
        identifier = self.__reqObj.getValue("identifier")
        #
        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        #
        sEx = SequenceDataExport(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        # ok, numConflicts, conflictList, warningMsg = sEx.exportAssignments()
        entityIdList = sEx.getAllEntityIdList()
        #
        pI = PathInfo(siteId=self.__siteId, sessionPath=self.__sessionPath, verbose=self.__verbose, log=self.__lfh)
        #
        partialFilePath = pI.getFilePath(dataSetId=identifier, contentType="partial-seq-annotate", formatType="txt", fileSource="session")
        fph = open(partialFilePath, "w")
        #
        copyFileList = []
        for entityId in entityIdList:
            for cType in ("seqdb-match", "seq-align-data", "mismatch-warning"):
                sourceFilePath = pI.getFilePath(dataSetId=identifier, contentType=cType, formatType="pic", fileSource="session", partNumber=entityId)
                if (not sourceFilePath) or (not os.access(sourceFilePath, os.R_OK)):
                    continue
                #
                exportFilePath = pI.getFilePath(dataSetId=identifier, contentType=cType, formatType="pic", partNumber=entityId, versionId="next")
                copyFileList.append((sourceFilePath, exportFilePath))
                fph.write("%s pic %s %s\n" % (cType, str(entityId), exportFilePath))
            #
        #
        foundFlag = True
        for (cType, fType, version) in (("model", "pdbx", "next"), ("seq-data-stats", "pic", "next"), ("partial-seq-annotate", "txt", "latest")):
            sourceFilePath = pI.getFilePath(dataSetId=identifier, contentType=cType, formatType=fType, fileSource="session")
            if (not sourceFilePath) or (not os.access(sourceFilePath, os.R_OK)):
                foundFlag = False
                continue
            #
            exportFilePath = pI.getFilePath(dataSetId=identifier, contentType=cType, formatType=fType, versionId=version)
            copyFileList.append((sourceFilePath, exportFilePath))
            if cType != "partial-seq-annotate":
                fph.write("%s %s 1 %s\n" % (cType, fType, exportFilePath))
            #
        #
        fph.close()
        #
        if (not foundFlag) or (len(copyFileList) < 5):
            rC.setError(errMsg="Saving partial assignments failed.")
            return rC
        #
        for ftuple in copyFileList:
            shutil.copyfile(ftuple[0], ftuple[1])
        #
        rC.setStatusCode("OK")
        return rC

    def _removePartialAssignmentOp(self):
        """Remove saved partial assignment files"""
        self.__getSession()
        self.__remove_partial_assignment_file()
        #
        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        rC.setStatusCode("OK")
        return rC

    def _storeAlignmentStartOp(self):
        return self.__storeAlignmentStart()

    def _storeAlignmentCheckOp(self):
        return self.__storeAlignmentCheck()

    def _editOp(self):
        return self.__editAlignment()

    def _moveEditOp(self):
        return self.__editAlignment()

    def _globalEditOp(self):
        self.__reqObj.setValue("operation", "global_edit")
        return self.__editAlignment()

    def _globalAuthSeqEditOp(self):
        self.__reqObj.setValue("operation", "global_edit_auth_seq")
        return self.__editAlignment()

    def _globalEditMenuOp(self):
        return self.__editAlignment()

    def _undoEditOp(self):
        return self.__editAlignment()

    def _deleteOp(self):
        return self.__editAlignment()

    def _respondDownloadOp(self):
        """Respond to download request for sequence assignment file - type = [model,assignment]"""
        self.__getSession()
        identifier = self.__reqObj.getValue("identifier")
        cType = self.__reqObj.getValue("type")
        #

        if cType in ["model"]:
            self.__reqObj.setValue("content_type", "model")
        elif cType in ["assignment"]:
            self.__reqObj.setValue("content_type", "seq-assign")
        self.__reqObj.setValue("data_set_id", identifier)
        self.__reqObj.setValue("file_source", "session")
        #
        wdu = WebDownloadUtils(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        rC = wdu.makeDownloadResponse()
        #
        if rC.isError():
            self.__reqObj.setReturnFormat(return_format="json")
        else:
            rC.setStatusCode("ok")
            self.__reqObj.setReturnFormat(return_format="binary")
        return rC

    def _respondTaxonomyFormOp(self):

        # Do this before updating the session information to avoid overwrite  -
        uds = UtilDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        seq_search_op = self.__reqObj.getValue("seq_search_op")
        uds.set("seq_search_op", seq_search_op)
        withref_info = self.__reqObj.getValue("withref_info")
        uds.set("withref_info", withref_info)
        uds.serialize()
        if self.__verbose:
            self.__lfh.write("+SeqModWebApp._respondTaxonomyFormOp() - starting with seq_search_op %s\n" % seq_search_op)
            self.__lfh.write("+SeqModWebApp._respondTaxonomyFormOp() - starting with withref_info %s\n" % withref_info)
        #
        self.__getSession()
        #
        sdU = UpdatePolymerEntityPartitions(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        status = sdU.polymerEntityPartEditFormResponder()
        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        if status:
            rC.setStatusCode("ok")
        else:
            rC.setError(errMsg="Taxonomy form update processing failure")
        return rC

    def _loadTaxonomyFormOp(self):
        self.__getSession()
        authId = self.__reqObj.getValue("auth_id")
        identifier = self.__reqObj.getValue("identifier")
        if authId is not None and len(authId) > 1:
            tL = str(authId).split("_")
            entityId = tL[2]
        else:
            entityId = "???"  # keep pylint happy
        cD = self.__makeTaxonomyEditForm(entityId=entityId, entryId=identifier)
        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        rC.addDictionaryItems(cD=cD)
        return rC

    def __makeTaxonomyEditForm(self, entityId="1", entryId=""):
        rD = {}
        sdU = UpdatePolymerEntityPartitions(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        rD["htmlcontent"] = sdU.makePolymerEntityPartEditForm(entityId=entityId, entryId=entryId)
        return rD

    def _respondSeqBuilderOp(self):
        self.__getSession()
        sdU = UpdatePolymerEntityPartitions(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        seqFetchError, hasRefSeq, returnSeq, returnSeqInfoList = sdU.seqBuilderResponder()
        #
        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        if seqFetchError:
            rC.setError(errMsg=seqFetchError)
        elif (returnSeq is not None) and len(returnSeq) > 0:
            rC.setStatusCode("ok")
            rC.setText(text=returnSeq)
            rC.set("seqpartinfolist", returnSeqInfoList)
            if hasRefSeq:
                rC.set("withref", "yes")
            else:
                rC.set("withref", "no")
            #
        else:
            rC.setError(errMsg="Building sequence processing failure")
        #
        return rC

    def _respondSeqDbRefFormOp(self):
        self.__getSession()
        sdU = UpdatePolymerEntityReference(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        seqFetchError, packedRefSeqLabel = sdU.seqDbRefFormResponder()
        #
        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        if seqFetchError:
            rC.setError(errMsg=seqFetchError)
        elif packedRefSeqLabel is not None:
            rC.setStatusCode("ok")
        else:
            rC.setError(errMsg="Reference sequence form update processing failure")
        #
        return rC

    def _loadSeqDbRefFormOp(self):
        self.__getSession()
        refId = self.__reqObj.getValue("ref_id")
        identifier = self.__reqObj.getValue("identifier")
        if refId is not None and len(refId) > 1:
            tL = str(refId).split("_")
            entityId = str(tL[1])
            partId = str(tL[3])
        else:
            entityId = ""
            partId = ""
        #
        sdU = UpdatePolymerEntityReference(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        cD = sdU.makeSeqdbrefEditForm(entityId, partId, entryId=identifier)

        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        rC.addDictionaryItems(cD=cD)
        return rC

    def _loadEntitySourceDetailsFormOp(self):
        self.__getSession()
        identifier = self.__reqObj.getValue("identifier")
        entityId = self.__reqObj.getValue("groupid")

        authId = self.__reqObj.getValue("auth_id")
        if authId is not None and len(authId) > 1:
            tL = str(authId).split("_")
            entityId = tL[2]
        #

        sdU = UpdatePolymerEntitySourceDetails(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        cD = sdU.makeEntitySourceDetailsForm(entityId, entryId=identifier)

        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        rC.addDictionaryItems(cD=cD)
        return rC

    def _respondEntitySourceDetailsFormOp(self):
        self.__getSession()
        sdU = UpdatePolymerEntitySourceDetails(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        ok = sdU.entitySourceDetailsFormResponder()
        #
        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        if ok:
            rC.setStatusCode("ok")
        else:
            rC.setError(errMsg="Entity review form update processing failure")

        return rC

    def _launchJmolViewerOp(self):
        self.__getSession()
        if self.__verbose:
            self.__lfh.write("+SeqModWebApp._launchJmolViewerOp() - starting with file source\n")

        viewer = ModelViewer3D(reqObj=self.__reqObj, verbose=self.__verbose)
        cD = viewer.launchJmol()
        #
        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        rC.addDictionaryItems(cD=cD)
        return rC

    def _launchAstexViewerOp(self):
        self.__getSession()
        viewer = ModelViewer3D(reqObj=self.__reqObj, verbose=self.__verbose)
        cD = viewer.launchAstex()
        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        rC.addDictionaryItems(cD=cD)
        return rC

    def _loadSummaryOp(self):
        """Loads pre-calculated HTML pages containing alignment summary --

        Performs a default selection of the best matching reference sequences.
        """
        return self.__loadSummary(op="load")

    def _reloadSummaryOp(self):
        """Reloads pre-calculated HTML pages containing alignment summary -"""
        return self.__loadSummary(op="reload")

    def _newSessionOp(self):
        """Entry point for new sessions when invoked as the sequence editor tool."""
        return self.__newSession()

    def _newSessionWfOp(self):
        """Entry point for new sessions when launched as a module from the wf environment."""
        return self.__newSessionWf()

    def _getPolymerLinkageTableOp(self):
        """Return the polymer linkage table"""
        self.__getSession()
        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)

        sessionId = self.__reqObj.getSessionId()
        # sessionPath    =self.__reqObj.getSessionPath()
        identifier = self.__reqObj.getValue("identifier")
        pI = PathInfo(siteId=self.__siteId, sessionPath=self.__sessionPath, verbose=self.__verbose, log=self.__lfh)
        polyLinkHtmlPath = pI.getPolyLinkReportFilePath(identifier, fileSource="session")
        if self.__verbose:
            self.__lfh.write("+SeqModWebApp.__getPolymerLinkageTableOp() report file path %s\n" % polyLinkHtmlPath)
        #
        htmlList = []
        if os.access(polyLinkHtmlPath, os.R_OK):
            ifh = open(polyLinkHtmlPath, "r")
            htmlList = ifh.readlines()
            ifh.close()
        else:
            htmlList = []

        (_head, tail) = os.path.split(polyLinkHtmlPath)
        PathRel = os.path.join("/sessions", sessionId, tail)
        rC.setHtmlContentPath(PathRel)
        rC.setHtmlList(htmlList)
        return rC

    # --------------------------------------------------------------------------------------------------
    #  WF data loading and initialization entry point
    #
    def _loadDataWfStartOp(self):
        """Load sequence data from model and sequence database matching results.  Compute preliminary
        alignment statistics.   Data is taken from the workflow storage system.
        """
        self.__getSession()
        fileSource = self.__reqObj.getValue("filesource")
        identifier = self.__reqObj.getValue("identifier")
        instance = self.__reqObj.getValue("instance")
        if fileSource in ["archive", "wf-archive"]:
            fileSource = "archive"

        if self.__verbose:
            self.__lfh.write("+SeqModWebApp._loadDataWfStartOp() - starting with file source %s identifier %s instance %s\n" % (fileSource, identifier, instance))
        rC = self.__loadDataStart(fileSource=fileSource)
        return rC

    # --------------------------------------------------------------------------------------------------
    #   Local data loading and initialization entry point methods ---
    #
    def _loadDataStartOp(self):
        """Load sequence data from model and perform sequence database matching results.  Compute preliminary
        alignment statistics.   Model data can be uploaded (fileSource=local-upload) or copied from
        data processing repository (fileSource=local-repository)
        """
        self.__getSession()
        if self.__verbose:
            self.__lfh.write("+SeqModWebApp._loadLocalDataStartOp() - starting with session %s\n" % self.__sessionPath)

        wuu = WebUploadUtils(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        if wuu.isFileUpload():
            fileSource = "local-upload"
            fileType = self.__reqObj.getValue("filetype")
            if self.__verbose:
                self.__lfh.write("+SeqModWebApp._loadLocalDataStartOp() - file source %s type %s\n" % (fileSource, fileType))

            modelFileName = wuu.copyToSession(fileTag="file")
            #
            if modelFileName is None:
                rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
                rC.setError(errMsg="Upload processing failure")
                self.__reqObj.setReturnFormat(return_format="json")
                return rC
            else:
                #
                fId, _idType = wuu.perceiveIdentifier(modelFileName)
                self.__reqObj.setValue("identifier", fId)
                self.__reqObj.setValue("UploadFileName", modelFileName)
                self.__reqObj.setValue("UploadFileType", fileType)

                uds = UtilDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
                uds.set("identifier", fId)
                uds.set("UploadFileName", modelFileName)
                uds.set("UploadFileType", fileType)
                uds.serialize()
            #
        else:
            fileSource = "local-repository"
            identifier = self.__reqObj.getValue("identifier")
            #
            uds = UtilDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            uds.set("identifier", identifier)
            uds.serialize()

            if self.__verbose:
                self.__lfh.write("+SeqModWebApp._loadLocalDataStartOp() - file source %s identifier %s\n" % (fileSource, identifier))
            #
        #
        rC = self.__loadDataStart(fileSource=fileSource)
        return rC

    def _reloadDataStartOp(self):
        """Launch a child process to handle reloading operations required by taxonomy/part edit operations.
        This method supports both local and wf data/file sources.
        """
        #
        self.__getSession()
        if self.__verbose:
            self.__lfh.write("+SeqModWebApp._reloadDataStartOp() - STARTING at site %s\n" % self.__siteId)
        #
        dU = DetachUtils(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        assembleUtil = SequenceDataAssemble(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        dU.set(workerObj=assembleUtil, workerMethod="doUpdateSelections")
        dU.runDetach()

        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        identifier = self.__reqObj.getValue("identifier")
        rC.setIdentifier(identifier)
        rC.setStatusCode("running")

        return rC

    def _loadDataCheckOp(self):
        """Shared checking method for data initialization and loading/uploading operations."""
        return self.__loadDataCheck()

    # --------------------------------------------------------------------------------------------------
    #   Re-run blast search entry point methods ---
    #
    def _reRunBlastOp(self):
        """Re-run blast search, perform sequence database matching results, compute preliminary alignment statistics."""
        self.__getSession()
        if self.__verbose:
            self.__lfh.write("+SeqModWebApp._reRunBlastOp() - starting with session %s\n" % self.__sessionPath)
        #

        dU = DetachUtils(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        assembleUtil = SequenceDataAssemble(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        dU.set(workerObj=assembleUtil, workerMethod="reRunSearch")
        dU.runDetach()

        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        identifier = self.__reqObj.getValue("identifier")
        rC.setIdentifier(identifier)
        rC.setStatusCode("running")

        return rC

    # --------------------------------------------------------------------------------------------
    #                              Supporting methods -
    #
    def __newSession(self):
        #
        if self.__verbose:
            self.__lfh.write("+SeqModWebApp.__newSession() - starting\n")

        self.__getSession(forceNew=True)
        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        return rC

    def __getSession(self, forceNew=False, useContext=True):
        """Join existing session or create new session as required.

        Import any saved parameters into the current request object from UtilDataStore().
        """
        #
        sessionId = self.__reqObj.getSessionId()
        if self.__verbose:
            self.__lfh.write("+SeqModWebApp.__getSession() - starting with session %s\n" % sessionId)
        #
        self.__sObj = self.__reqObj.newSessionObj(forceNew=forceNew)

        # self.__sessionId = self.__sObj.getId()
        self.__sessionPath = self.__sObj.getPath()
        # self.__rltvSessionPath = self.__sObj.getRelativePath()
        #
        if useContext:
            uds = UtilDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            dd = uds.getDictionary()
            self.__lfh.write("+SeqModWebApp.__getSession() - importing persisted session parameters: %r\n" % dd.items())
            self.__reqObj.setDictionary(dd, overWrite=True)
        #
        if self.__verbose:
            self.__lfh.write("+SeqModWebApp.__getSession() - Leaving with session path %s\n" % self.__sessionPath)

    def __loadDataStart(self, fileSource="local-repository"):
        """Launch a child process to handle data initialization and loading operations.

        This method supports both local and wf data/file sources.
        """
        #
        if self.__verbose:
            self.__lfh.write("+SeqModWebApp.__loadDataStart() - STARTING using file source %s at site %s\n" % (fileSource, self.__siteId))

        dU = DetachUtils(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)

        dI = DataImporter(reqObj=self.__reqObj, fileSource=fileSource, verbose=self.__verbose, log=self.__lfh)

        dU.set(workerObj=dI, workerMethod="loadData")

        dU.runDetach()

        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        identifier = self.__reqObj.getValue("identifier")
        rC.setIdentifier(identifier)
        rC.setStatusCode("running")

        return rC

    def __loadDataCheck(self):
        """Performs a check on the contents of a semaphore file and returns the associated status.

        This method currently supports both rcsb and wf filesources.
        """
        #
        self.__getSession(useContext=True)
        if self.__verbose:
            self.__lfh.write("+SeqModWebApp.__loadDataCheck() - starting\n")

        sph = self.__reqObj.getSemaphore()
        delayValue = self.__reqObj.getValue("delay")
        if self.__verbose:
            self.__lfh.write("+SeqModWebApp.__loadDataCheck() Checking status of semaphore %s with delay %s\n" % (sph, delayValue))

        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        identifier = self.__reqObj.getValue("identifier")
        rC.setIdentifier(identifier)

        dU = DetachUtils(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        if dU.semaphoreExists(sph):
            status = dU.getSemaphore(sph)
            if self.__verbose:
                self.__lfh.write("+SeqModWebApp.__loadDataCheck() status value for semaphore %s is %s\n" % (sph, str(status)))
            if status == "OK":
                rC.setStatusCode("completed")
            else:
                rC.setStatusCode("failed")
        else:
            if self.__verbose:
                self.__lfh.write("+SeqModWebApp.__loadDataCheck() semaphore %s pending - waiting %s\n" % (sph, delayValue))
            time.sleep(int(delayValue))
            rC.setStatusCode("running")

        return rC

    def __loadSummary(self, op="load"):
        #
        if self.__verbose:
            self.__lfh.write("+SeqModWebApp.__loadSummary() - starting\n")

        self.__getSession(useContext=True)
        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        try:
            sV = SummaryView(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            sumObj, warningMsg = sV.loadSummary(operation=op)
            if warningMsg:
                rC.setWarning(warnMsg=warningMsg)
            #
            eD = sV.getEntryDetails()
            gIdList = sV.getGroupIdList()
            activeGroupId = self.__reqObj.getValue("activegroupid")
            sVD = SummaryViewDepiction(verbose=self.__verbose)
            oL = sVD.buildSummaryView(sumObj, activeGroupId=activeGroupId)
            #
            sessionId = self.__reqObj.getSessionId()
            sPath = self.__reqObj.getSessionPath()
            htmlPathAbs = os.path.join(sPath, sessionId, "current-alignment-summary.html")
            fp = open(htmlPathAbs, "w")
            fp.write("%s" % "".join(oL))
            fp.close()
            rC.setHtmlContentPath(os.path.join("/sessions", sessionId, "current-alignment-summary.html"))
            #
            rC.setIdentifier(self.__reqObj.getValue("identifier"))
            rC.setStructTitle(eD["STRUCT_TITLE"])
            rC.setCitationTitle(eD["CITATION_TITLE"])
            rC.setPdbCode(eD["PDB_ID"])
            rC.setGroupIdList(gIdList)
            #
            pI = PathInfo(siteId=self.__siteId, sessionPath=self.__sessionPath, verbose=self.__verbose, log=self.__lfh)
            sourceFilePath = pI.getFilePath(dataSetId=self.__reqObj.getValue("identifier"), contentType="partial-seq-annotate", formatType="txt")
            if sourceFilePath and os.access(sourceFilePath, os.R_OK):
                rC.set("haspartialassignment", True)
            #
            # Set warning message
            #
            picklePath = pI.getFilePath(self.__reqObj.getValue("identifier"), contentType="mismatch-warning", formatType="pic", fileSource="session")
            if picklePath and os.access(picklePath, os.F_OK):
                fb = open(picklePath, "rb")
                warningD = pickle.load(fb)
                fb.close()
                #
                warningMsg = ""
                if ("seq_warning_info" in warningD) and len(warningD["seq_warning_info"]) > 0:
                    warningMsg = self.__reformatText(warningD["seq_warning_info"])
                #
                if ("mismatch" in warningD) and len(warningD["mismatch"]) > 0:
                    if warningMsg:
                        warningMsg += "</br>\n"
                    #
                    if len(warningD["mismatch"]) > 1:
                        warningMsg += self.__reformatText("Sequence/coordinates mismatch in entities: " + ",".join(warningD["mismatch"]))
                    else:
                        warningMsg += "Sequence/coordinates mismatch in entity: " + warningD["mismatch"][0]
                    #
                #
                if ("not_found_existing_match" in warningD) and len(warningD["not_found_existing_match"]) > 0:
                    if warningMsg:
                        warningMsg += "</br>\n"
                    #
                    if len(warningD["not_found_existing_match"]) > 1:
                        warningMsg += self.__reformatText("Sequences not matched to existing entries in entities: " + ",".join(warningD["not_found_existing_match"]))
                    else:
                        warningMsg += "Sequences not matched to existing entries in entity: " + warningD["not_found_existing_match"][0]
                    #
                #
                for warningType in ("missing_residue", "mix_mse_met"):
                    if (warningType in warningD) and len(warningD[warningType]) > 0:
                        for textMsg in warningD[warningType]:
                            if warningMsg:
                                warningMsg += "</br>\n"
                            #
                            warningMsg += textMsg
                        #
                    #
                #
                if warningMsg:
                    rC.set("entrywarningmessage", warningMsg)
                #
            #
            # Set summary page content
            #
            rC.set("summaryinfo", self.__getSummaryPageContent(sV.getSummaryPageObj()))
            #
            #  Reset --
            #
            self.__reqObj.setNewRefId("")
        except:  # noqa: E722 pylint: disable=bare-except
            rC.setError(errMsg="Sequence alignment summary preparation has failed.")
            rC.setStatusCode("failed")
            if self.__verbose:
                self.__lfh.write("+SeqModWebApp.__loadSummary() - failed with exception\n")
                traceback.print_exc(file=self.__lfh)
        #
        return rC

    def __reformatText(self, inputText):
        """ """
        length = 150
        if len(inputText) < length:
            return inputText
        #
        break_index = []
        for i in range(0, len(inputText)):
            if (inputText[i] == " ") or (inputText[i] == ","):
                break_index.append(i)
            #
        #
        if len(break_index) == 0:
            return inputText
        #
        rangeList = []
        insertList = []
        start = 0
        prev_i = break_index[0]
        for i in break_index:
            if (i - start) > length:
                rangeList.append((start, prev_i + 1))
                insertList.append(prev_i)
                start = prev_i + 1
                continue
            #
            prev_i = i
        #
        if (prev_i > start) and (prev_i < len(inputText)) and ((len(insertList) == 0) or (prev_i != insertList[-1])) and ((len(inputText) - start) > length):
            rangeList.append((start, prev_i + 1))
            start = prev_i + 1
        #
        if start < len(inputText):
            rangeList.append((start, len(inputText)))
        #
        outputText = ""
        for rangeTup in rangeList:
            if outputText:
                outputText += "</br>\n"
            #
            outputText += inputText[rangeTup[0] : rangeTup[1]]
        #
        return outputText

    def __getSummaryPageContent(self, spObj):
        """ Generate summary page content
        """
        columnNameKeyList = (("Entity", "entity_id"), ("Chain", "chain_ids"), ("Type", "entity_types"),
                             ("Molecule name<br/>(Depositor-provided name)", "mol_names"),
                             ("Source name/strain<br/>(Depositor-provided source name)", "source_names"),
                             ("Source TaxID<br/>(Depositor-provided TaxID)", "tax_ids"),
                             ("Scores", "identity_scores"),
                             ("Top hit ID [TaxID]<br/>(Depositor-provided ID)", "ref_db_ids"),
                             ("Mismatch and differences", "warn_err_msgs"))
        #
        oL = []
        oL.append('<table class="summary_table">\n')
        oL.append("<thead>\n")
        oL.append("<tr>")
        for columnNameKeTuple in columnNameKeyList:
            oL.append("<th>%s</th>" % columnNameKeTuple[0])
        #
        oL.append("</tr>\n")
        oL.append("</thead>\n")
        oL.append("<tbody>\n")
        #
        gIdList = list(spObj.keys())
        gIdList.sort(key=int)
        oddEven = True
        for gId in gIdList:
            gObj = spObj[gId]
            if oddEven:
                oL.append('<tr class="odd">')
                oddEven = False
            else:
                oL.append('<tr class="even">')
                oddEven = True
            #
            for columnNameKeTuple in columnNameKeyList:
                if columnNameKeTuple[1] in gObj:
                    oL.append("<td>%s</td>" % gObj[columnNameKeTuple[1]])
                else:
                    oL.append("<td>-</td>")
                #
            #
            oL.append("</tr>\n")
        #
        oL.append("</tbody>\n")
        oL.append("</table>\n")
        return "".join(oL)

    def __storeAlignmentStart(self):
        """Launch a subprocess to compute, store and render the aligned input sequence list."""
        #
        self.__getSession(useContext=True)
        _sessionId = self.__reqObj.getSessionId()  # noqa: F841
        alignTag = self.__reqObj.getValue("aligntag")
        operation = self.__reqObj.getValue("operation")
        viewalign_order = self.__reqObj.getValue("viewalign_order")
        alignids = self.__reqObj.getValue("alignids")
        activeGroupId = self.__reqObj.getValue("activegroupid")

        if self.__verbose:
            self.__lfh.write("+SeqModWebApp.__storeAlignmentStart() - starting\n")
            self.__lfh.write("+SeqModWebApp.__storeAlignmentStart() operation       %s\n" % operation)
            self.__lfh.write("+SeqModWebApp.__storeAlignmentStart() aligntag        %s\n" % alignTag)
            self.__lfh.write("+SeqModWebApp.__storeAlignmentStart() activegroupid   %s\n" % activeGroupId)
            self.__lfh.write("+SeqModWebApp.__storeAlignmentStart() alignids        %s\n" % alignids)
            self.__lfh.write("+SeqModWebApp.__storeAlignmentStart() viewalign_order %s\n" % viewalign_order)

        if len(alignTag) < 1:
            alignTag = activeGroupId
            if self.__verbose:
                self.__lfh.write("+SeqModWebApp.__storeAlignmentStart() creating new alignment tag %s\n" % alignTag)
            #
        #
        self.__reqObj.setValue("aligntag", alignTag)

        dU = DetachUtils(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        aD = AlignmentDepictionTools(reqObj=self.__reqObj, entityId=alignTag, verbose=self.__verbose, log=self.__lfh)
        dU.set(workerObj=aD, workerMethod="doRender")
        dU.runDetach()

        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        rC.setStatusCode("running")
        rC.setAlignTag(alignTag)
        rC.setViewAlignOrder(viewalign_order)
        # reset edit operation -
        rC.setEditOp(0)
        self.__lfh.write("+SeqModWebApp.__storeAlignmentStart() Parent process completed\n")
        return rC

    def __storeAlignmentCheck(self):
        #
        self.__getSession()
        if self.__verbose:
            self.__lfh.write("\n+SeqModWebApp.__storeAlignmentCheck() - starting\n")
        #
        sph = self.__reqObj.getSemaphore()
        sessionId = self.__reqObj.getSessionId()
        alignTag = self.__reqObj.getValue("aligntag")
        delayValue = self.__reqObj.getValue("delay")
        #
        if self.__verbose:
            self.__lfh.write("+SeqModWebApp.__storeAlignmentCheck() Checking status of semaphore %s with delay %s\n" % (sph, delayValue))
        #
        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        dU = DetachUtils(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)

        if dU.semaphoreExists(sph):
            status = dU.getSemaphore(sph)
            if self.__verbose:
                self.__lfh.write("+SeqModWebApp.__storeAlignmentCheck() status value for semaphore %s is %s\n" % (sph, str(status)))
            #
            rC.setIdentifier(self.__reqObj.getValue("identifier"))
            rC.setStructTitle(self.__reqObj.getValue("title"))
            rC.setPdbCode(self.__reqObj.getValue("pdbid"))
            #
            rC.setConflictReportFlag(False)
            rC.setAnnotationReportFlag(False)
            #
            sPath = self.__reqObj.getSessionPath()
            rptPathAbs = os.path.join(sPath, sessionId, "current-alignment-annotation-" + alignTag + ".html")
            if os.access(rptPathAbs, os.F_OK):
                pth = os.path.join("/sessions", sessionId, "current-alignment-annotation-" + alignTag + ".html")
                rC.setAnnotationReportPath(pth)
                rC.setAnnotationReportFlag(True)
            #
            miscD = {}
            rptPathAbs = os.path.join(sPath, sessionId, "alignment-" + alignTag + "-misc.pic")
            if os.access(rptPathAbs, os.F_OK):
                fb = open(rptPathAbs, "rb")
                miscD = pickle.load(fb)
                fb.close()
            #
            if "warning" in miscD:
                rC.setWarning(warnMsg=miscD["warning"])
            #
            if status == "OK":
                # reset edit operation -
                rC.setEditOp(0)
                rC.setAlignTag(alignTag)
                rC.setStatusCode("completed")
                #
                rptPathAbs = os.path.join(sPath, sessionId, "current-alignment-" + alignTag + ".html")
                if os.access(rptPathAbs, os.F_OK):
                    rC.setHtmlContentPath(os.path.join("/sessions", sessionId, "current-alignment-" + alignTag + ".html"))
                #
                rptPathAbs = os.path.join(sPath, sessionId, "conflict-report-" + alignTag + ".html")
                if os.access(rptPathAbs, os.F_OK):
                    pth = os.path.join("/sessions", sessionId, "conflict-report-" + alignTag + ".html")
                    rC.setConflictReportPath(pth)
                    rC.setConflictReportFlag(True)
                    if self.__verbose:
                        self.__lfh.write("+SeqModWebApp.__storeAlignmentCheck() conflict report %s\n" % pth)
                    #
                #
                if "gedittype" in miscD:
                    rC.addDictionaryItems({"gedittype": miscD["gedittype"]})
                else:
                    rC.addDictionaryItems({"gedittype": "no-mismatch"})
                #
                for item in ("alignids", "selectids", "alignmentblocklist", "missingauthseqmap", "blockedithtml", "repdelhtml"):
                    if item in miscD:
                        rC.addDictionaryItems({item: miscD[item]})
                    #
                #
            else:
                rC.setStatusCode("failed")
            #
        else:
            if self.__verbose:
                self.__lfh.write("+SeqModWebApp.__storeAlignmentCheck() semaphore %s pending - waiting %s\n" % (sph, delayValue))
            #
            time.sleep(int(delayValue))
            rC.setStatusCode("running")
        #
        return rC

    def __editAlignment(self):
        """ """
        sessionId = self.__reqObj.getSessionId()
        alignTag = self.__reqObj.getValue("aligntag")
        #
        operation = self.__reqObj.getValue("operation")
        if self.__verbose:
            self.__lfh.write("+SeqModWebApp.__editAlignment() - starting\n")
            self.__lfh.write("+SeqModWebApp.__editAlignment() - session id %s\n" % sessionId)
            self.__lfh.write("+SeqModWebApp.__editAlignment() - operation %s\n" % operation)
            self.__lfh.write("+SeqModWebApp.__editAlignment() - alignment tag %s\n" % alignTag)

        #
        aE = AlignmentFrontEndUIEditor(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        aE.setAlignmentTag(alignTag)
        #
        cD = aE.edit()
        if self.__verbose:
            self.__lfh.write("+SeqModWebApp.__editAlignment() edit completed for session %s\n" % sessionId)
        #
        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        rC.addDictionaryItems(cD=cD)
        return rC

    def __saveSelection(self):
        """Save assignment and annotation file"""
        self.__getSession()
        if self.__verbose:
            self.__lfh.write("+SeqModWebApp.__saveSelection() - starting\n")

        # selectIdList = self.__reqObj.getSummarySelectList()
        selectIdList = []
        sEx = SequenceDataExport(reqObj=self.__reqObj, exportList=selectIdList, verbose=self.__verbose, log=self.__lfh)
        ok, numConflicts, conflictList, warningMsg = sEx.exportAssignments()
        if self.__verbose:
            self.__lfh.write("+SeqModWebApp.__saveSelection() - exporting assignments status = %r and conflict count  %d\n" % (ok, numConflicts))
        # if (ok and (numConflicts < 3)):
        if ok and (numConflicts == 0):
            ok = sEx.applyAssignmentsToModel()
            if self.__verbose:
                self.__lfh.write("+SeqModWebApp.__saveSelection() - exporting model status = %r\n" % ok)
        #
        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        if warningMsg:
            rC.setWarning(warnMsg=warningMsg)
        #
        if not ok:
            rC.setError(errMsg="Sequence data export failed")
        elif numConflicts > 0:
            oL = []
            for ctup in conflictList:
                oL.append("Entity %s(%d)" % (ctup[0], ctup[1]))
            tS = " , ".join(oL)
            rC.setError(errMsg="Error: Total sample conflict count = %d  - %s</br>\nSequence data export failed" % (numConflicts, tS))
        else:
            self.__remove_partial_assignment_file()
            rC.set("removepartialassignment", True)
        #
        return rC

    def __closeSeqEditor(self, mode):
        """RPS, 2010-Aug-02: providing dedicated function to accommodate user request to end sesssion when
        done with using sequence editor interface
        RPS, 2010-Aug-12: refactored to support different 'modes' = ('completed' | 'unfinished')

        """
        self.__getSession()
        if self.__verbose:
            self.__lfh.write("\n\n+SeqModWebApp.__closeSeqEditor() - starting with operation %s\n" % mode)

        state = "???"  # pylint noted could be unset
        if mode == "completed":
            state = "closed(0)"
        elif mode == "unfinished":
            state = "waiting"

        depId = self.__reqObj.getValue("identifier")
        instId = self.__reqObj.getValue("instance")
        classId = self.__reqObj.getValue("classid")
        sessionId = self.__reqObj.getSessionId()
        fileSource = str(self.__reqObj.getValue("filesource")).strip().lower()
        #
        skipStatus = str(self.__reqObj.getValue("skipstatus")).strip()
        if skipStatus.lower() in ["yes", "y"]:
            self.__doWfTracking = False
        elif skipStatus.lower() in ["no", "n"]:
            self.__doWfTracking = True

        if self.__verbose:
            self.__lfh.write("+SeqModWebApp.__closeSeqEditor() - depId   %s\n" % depId)
            self.__lfh.write("+SeqModWebApp.__closeSeqEditor() - instId  %s\n" % instId)
            self.__lfh.write("+SeqModWebApp.__closeSeqEditor() - classID %s\n" % classId)
            self.__lfh.write("+SeqModWebApp.__closeSeqEditor() - sessionID %s\n" % sessionId)
            self.__lfh.write("+SeqModWebApp.__closeSeqEditor() - filesource %r\n" % fileSource)
            self.__lfh.write("+SeqModWebApp.__closeSeqEditor() - doWfTracking  %r\n" % self.__doWfTracking)
            self.__lfh.write("+SeqModWebApp.__closeSeqEditor() - close mode  %r\n" % mode)

        #
        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)

        if fileSource in ["archive", "wf-archive"]:
            # copy session files back to the archive directory -- independent of close mode --
            dI = DataImporter(reqObj=self.__reqObj, fileSource=fileSource, verbose=self.__verbose, log=self.__lfh)
            ok = dI.copyFilesOnClose(includeSeqAssignFile=True)

        if fileSource in ["wf-instance"]:
            # Copy session files back to the invoking instance directory -
            dI = DataImporter(reqObj=self.__reqObj, fileSource=fileSource, verbose=self.__verbose, log=self.__lfh)
            ok = dI.copyFilesOnClose(includePolyLinkFile=True, includeSeqAssignFile=True)
            # And if we are done then also copy back to the archive --
            ok1 = True
            if mode in ["completed"]:
                ok1 = dI.copyFiles(
                    inputFileSource="session",
                    outputFileSource="wf-archive",
                    versionIndex=4,
                    includeModelFile=True,
                    includeSeqAssignFile=True,
                    messageHead="DataImporter.copyFilesOnClose()",
                )
                #
                self.__getSession()
                pI = PathInfo(siteId=self.__siteId, sessionPath=self.__sessionPath, verbose=self.__verbose, log=self.__lfh)
                sourceFilePath = pI.getFilePath(dataSetId=depId, contentType="model", formatType="pdbx", fileSource="session")
                if os.access(sourceFilePath, os.R_OK):
                    dbLoader = DBLoadUtil(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
                    dbLoader.doLoading([sourceFilePath])
                #
            #
            ok = ok and ok1

        okDb = True
        if self.__doWfTracking:
            try:
                okDb = self.__updateWfStatus(depId=depId, instId=instId, classId=classId, statusCode=state)
            except:  # noqa: E722 pylint: disable=bare-except
                okDb = False

        ok = ok and okDb

        if not ok:
            self.__lfh.write("+SeqModWebApp.__closeSeqEditor() - Quit/Finish operation failed\n")
            rC.setError(errMsg="Quit/Finish file or database update operation has failed")
        else:
            self.__lfh.write("+SeqModWebApp.__closeSeqEditor() - status update completed\n")
            rC.setHtmlText("Quit/Finish operation successful!")

        return rC

    def __newSessionWf(self):

        if self.__verbose:
            self.__lfh.write("+SeqModWebApp.__newSessionWf() - starting\n")

        sessionId = self.__reqObj.getSessionId()
        if len(sessionId) == 0:
            sObj = self.__reqObj.newSessionObj()
            sessionId = sObj.getId()

        #
        #  Attibutes that will fill the HTML template...
        #
        identifier = self.__reqObj.getValue("identifier")
        instance = self.__reqObj.getValue("instance")
        fileSource = self.__reqObj.getValue("filesource")
        skipStatus = str(self.__reqObj.getValue("skipstatus")).strip()
        if skipStatus.lower() in ["yes", "y"]:
            self.__doWfTracking = False
        elif skipStatus.lower() in ["no", "n"]:
            self.__doWfTracking = True

        #
        myD = {}
        myD["sessionid"] = sessionId
        myD["identifier"] = identifier
        myD["instance"] = instance
        myD["classid"] = self.__reqObj.getValue("classID")
        myD["filesource"] = fileSource

        myD["height"] = "100%"
        # myD['height']='1000px'
        myD["width"] = "100%"
        #
        uds = UtilDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        uds.set("identifier", identifier)
        uds.set("filesource", fileSource)
        uds.set("instance", instance)
        uds.set("skipstatus", skipStatus)
        uds.serialize()
        try:
            _ok = self.__updateWfStatus(depId=myD["identifier"], instId=myD["instance"], classId=myD["classid"], statusCode="open")  # noqa: F841
        except:  # noqa: E722 pylint: disable=bare-except
            _ok = False  # noqa: F841
        #
        self.__reqObj.setReturnFormat(return_format="html")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        templateFilePath = os.path.join(self.__reqObj.getValue("TemplatePath"), "summary_template.html")
        webIncludePath = os.path.join(self.__reqObj.getValue("TopPath"), "htdocs")
        rC.setHtmlTextFromTemplate(templateFilePath=templateFilePath, webIncludePath=webIncludePath, parameterDict=myD)
        return rC

    def __alignView(self, frameOpt=True):
        """Send a skeleton HTML template with placeholders for the sequence alignment and conflict report."""
        viewalign_order = self.__reqObj.getValue("viewalign_order")
        if self.__verbose:
            self.__lfh.write("+SeqModWebApp.__alignView() - starting with %s \n" % viewalign_order)
        #
        self.__getSession(forceNew=False, useContext=True)

        sessionId = self.__reqObj.getSessionId()
        alignTag = self.__reqObj.getValue("aligntag")
        viewalign_order = self.__reqObj.getValue("viewalign_order")
        activeEntityGroupId = self.__reqObj.getValue("activegroupid")

        if self.__verbose:
            self.__lfh.write("+SeqModWebApp.__alignView() Loading alignment template for sessionId %s alignment tag %s order %s\n" % (sessionId, alignTag, viewalign_order))
        myD = {}
        myD["sessionid"] = self.__reqObj.getSessionId()
        myD["seqview"] = "fromlist"
        myD["alignids"] = ",".join(self.__reqObj.getAlignList())
        # myD["selectids"] = ",".join(self.__reqObj.getSummarySelectList())
        myD["selectids"] = self.__reqObj.getValue("selectids")
        myD["aligntag"] = activeEntityGroupId
        myD["activegroupid"] = activeEntityGroupId

        myD["identifier"] = self.__reqObj.getValue("identifier")
        myD["instance"] = self.__reqObj.getValue("instance")
        myD["viewalign_order"] = self.__reqObj.getValue("viewalign_order")
        myD["modelfilename"] = self.__reqObj.getValue("modelfilename")

        if self.__verbose:
            self.__lfh.write("+SeqModWebApp.__alignView() parameter dictionary %s\n" % myD.items())

        self.__reqObj.setReturnFormat(return_format="html")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        templateFilePath = os.path.join(self.__reqObj.getValue("TemplatePath"), "alignment_template.html")
        webIncludePath = os.path.join(self.__reqObj.getValue("TopPath"), "htdocs")
        if frameOpt:
            # no difference here always go to frame --
            rC.setHtmlTextFromTemplate(templateFilePath=templateFilePath, webIncludePath=webIncludePath, parameterDict=myD)
        else:
            rC.setHtmlTextFromTemplate(templateFilePath=templateFilePath, webIncludePath=webIncludePath, parameterDict=myD)
        rC.setViewAlignOrder(viewalign_order)
        if self.__verbose:
            self.__lfh.write("+SeqModWebApp.__alignView() Just set viewalign_order in ResponseContent as %s\n" % viewalign_order)

        return rC

    def __updateWfStatus(self, depId, instId, classId, statusCode):
        """Update progress and tracking status --"""
        if not self.__doWfTracking:
            if self.__verbose:
                self.__lfh.write("+SeqModWebApp.__updateWfStatus() - skipped setting depId %s instId %s classId %s status %s\n" % (depId, instId, classId, statusCode))
            return True
        try:
            if self.__verbose:
                self.__lfh.write("+SeqModWebApp.__updateWfStatus() - setting depId %s instId %s classId %s status %s\n" % (depId, instId, classId, statusCode))
            wft = WfTracking(verbose=self.__verbose, log=self.__lfh)
            ok = wft.setInstanceStatus(depId=depId, instId=instId, classId=classId, status=statusCode)
            return ok
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+SeqModWebApp.__updateWfStatus() - failed depId %s instId %s classId %s status %s\n" % (depId, instId, classId, statusCode))
            traceback.print_exc(file=self.__lfh)
        return False

    def __remove_partial_assignment_file(self):
        """Remove saved partial assignment files"""
        identifier = self.__reqObj.getValue("identifier")
        #
        pI = PathInfo(siteId=self.__siteId, sessionPath=self.__sessionPath, verbose=self.__verbose, log=self.__lfh)
        sourceFilePath = pI.getFilePath(dataSetId=identifier, contentType="partial-seq-annotate", formatType="txt")
        if sourceFilePath and os.access(sourceFilePath, os.R_OK):
            fph = open(sourceFilePath, "r")
            inData = fph.read()
            fph.close()
            #
            for line in inData.split("\n"):
                tupL = line.split(" ")
                if len(tupL) == 4 and os.access(tupL[3], os.R_OK):
                    os.remove(tupL[3])
                #
            #
            os.remove(sourceFilePath)
        #


def main():
    sTool = SeqModWebApp()
    d = sTool.doOp()
    for k, v in d.items():
        sys.stdout.write("Key - %s  value - %r\n" % (k, v))


if __name__ == "__main__":
    main()
