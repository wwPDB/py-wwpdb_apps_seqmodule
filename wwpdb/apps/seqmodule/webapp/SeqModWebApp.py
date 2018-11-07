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

import filecmp, os, shutil, string, sys, time, traceback, types
from json import loads, dumps

from wwpdb.apps.seqmodule.webapp.SeqModWebRequest import SeqModInputRequest, SeqModResponseContent
from wwpdb.apps.seqmodule.view3d.ModelViewer3D import ModelViewer3D
from wwpdb.apps.seqmodule.control.SummaryView import SummaryView
from wwpdb.apps.seqmodule.control.SummaryViewDepiction import SummaryViewDepiction
from wwpdb.apps.seqmodule.control.DataImporter import DataImporter
#
from wwpdb.apps.seqmodule.update.UpdatePolymerEntityPartitions import UpdatePolymerEntityPartitions
from wwpdb.apps.seqmodule.update.UpdatePolymerEntitySourceDetails import UpdatePolymerEntitySourceDetails
from wwpdb.apps.seqmodule.update.UpdatePolymerEntityReference import UpdatePolymerEntityReference
#
from wwpdb.apps.seqmodule.control.ReferenceSequenceDataUpdate import ReferenceSequenceDataUpdate

from wwpdb.apps.seqmodule.align.MultiAlignPseudo import MultiAlignPseudo
from wwpdb.apps.seqmodule.align.AlignmentViewDepiction import AlignmentViewDepiction
from wwpdb.apps.seqmodule.align.AlignmentEdit import AlignmentEdit
from wwpdb.apps.seqmodule.align.AlignmentUtils import AlignmentUtils
from wwpdb.apps.seqmodule.util.SequenceLabel import SequenceLabelUtils

from wwpdb.apps.seqmodule.io.SequenceDataExport import SequenceDataExport
from wwpdb.api.status.dbapi.WfTracking import WfTracking

from wwpdb.api.facade.ConfigInfo import ConfigInfo
from wwpdb.utils.rcsb.DetachUtils import DetachUtils
from wwpdb.utils.rcsb.WebUploadUtils import WebUploadUtils
from wwpdb.utils.rcsb.WebDownloadUtils import WebDownloadUtils
from wwpdb.utils.rcsb.UtilDataStore import UtilDataStore
from wwpdb.utils.rcsb.PathInfo import PathInfo


class SeqModWebApp(object):
    """Handle request and response object processing for sequence editor tool application.

    """

    def __init__(self, parameterDict={}, verbose=False, log=sys.stderr, siteId="WWPDB_DEPLOY_TEST"):
        """
        Create an instance of `SeqModWebApp` to manage a sequence editor web request.

         :param `parameterDict`: dictionary storing parameter information from the web request.
             Storage model for GET and POST parameter data is a dictionary of lists.
         :param `verbose`:  boolean flag to activate verbose logging.
         :param `log`:      stream for logging.

        """
        self.__verbose = verbose
        self.__lfh = log
        self.__debug = False
        self.__siteId = siteId
        self.__cI = ConfigInfo(self.__siteId)
        #
        # the <site-specificit...>/webapps  directory path --  containing htdocs & fcgi subdirs
        self.__topPath = self.__cI.get('SITE_WEB_APPS_TOP_PATH')

        # The path containing the sessions directory (i.e. <top_sessions_path>/sessions )
        self.__topSessionPath = self.__cI.get('SITE_WEB_APPS_TOP_SESSIONS_PATH')

        if isinstance(parameterDict, types.DictType):
            self.__myParameterDict = parameterDict
        else:
            self.__myParameterDict = {}

        if (self.__debug):
            self.__lfh.write("+SeqModWebApp.__init() - SERVICE REQUEST STARTING with input parameter dictionary \n")
            self.__lfh.write("%s" % (''.join(self.__dump())))

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
        #self.__examplePath      =os.path.join(self.__topPath,"scripts","wwpdb","apps","seqmodule","examples")
        #self.__examplePathData  =os.path.join(self.__examplePath,"rcsb-data")
        #self.__examplePathSeq   =os.path.join(self.__examplePath,"rcsb-sequence")
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
        if (self.__verbose):
            rqp = self.__reqObj.getRequestPath()
            self.__lfh.write("\n#########################################################################################################\n")
            self.__lfh.write("+SeqModWebApp.__doOP() - Beginning request:  %s\n" % rqp)
            self.__reqObj.printIt(ofh=self.__lfh)

        stw = SeqModWebAppWorker(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        # doOp returns a response content object --
        rC = stw.doOp()
        if (self.__debug):
            self.__lfh.write("+SeqModWebApp.doOp() working method completed with response format %s\n" % self.__reqObj.getReturnFormat())
            if rC is not None:
                self.__lfh.write("%s" % (''.join(rC.dump())))
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
        self.__reqObj = reqObj
        self.__doWfTracking = True
        #
        self.__appPathD = {'/service/environment/dump': '_dumpOp',
                           '/service/sequence_editor/store_alignment/start': '_storeAlignmentStartOp',
                           '/service/sequence_editor/store_alignment/check': '_storeAlignmentCheckOp',
                           '/service/sequence_editor/load_summary': '_loadSummaryOp',
                           '/service/sequence_editor/reload_summary': '_reloadSummaryOp',
                           '/service/sequence_editor/save': '_saveOp',
                           '/service/sequence_editor/close_unfinished': '_closeUnfinishedOp',
                           '/service/sequence_editor/close_completed': '_closeCompletedOp',
                           '/service/sequence_editor/save_session_info': '_saveSessionInfoOp',
                           '/service/sequence_editor/new_session': '_newSessionOp',
                           '/service/sequence_editor/new_session/wf': '_newSessionWfOp',
                           '/service/sequence_editor/edit': '_editOp',
                           '/service/sequence_editor/global_edit': '_globalEditOp',
                           '/service/sequence_editor/global_edit_form': '_globalEditFormOp',
                           '/service/sequence_editor/global_edit_menu': '_globalEditMenuOp',
                           '/service/sequence_editor/move': '_moveEditOp',
                           '/service/sequence_editor/undo_edit': '_undoEditOp',
                           '/service/sequence_editor/delete': '_deleteOp',
                           '/service/sequence_editor/molviewer/jmol': '_launchJmolViewerOp',
                           '/service/sequence_editor/molviewer/astexviewer': '_launchAstexViewerOp',
                           #
                           '/service/sequence_editor/load_data/start/rcsb': '_loadDataStartOp',
                           '/service/sequence_editor/load_data/check/rcsb': '_loadDataCheckOp',
                           '/service/sequence_editor/load_data/check': '_loadDataCheckOp',
                           #
                           '/service/sequence_editor/load_data/start/wf': '_loadDataWfStartOp',
                           '/service/sequence_editor/load_data/check/wf': '_loadDataCheckOp',
                           #
                           '/service/sequence_editor/rerun_blast/start': '_reRunBlastOp',
                           #
                           '/service/sequence_editor/load_form/taxonomy': '_loadTaxonomyFormOp',
                           '/service/sequence_editor/load_form/seqdbref': '_loadSeqDbRefFormOp',
                           '/service/sequence_editor/load_form/entityreview': '_loadEntitySourceDetailsFormOp',

                           '/service/sequence_editor/respond_form/taxonomy': '_respondTaxonomyFormOp',
                           '/service/sequence_editor/respond_form/seqdbref': '_respondSeqDbRefFormOp',
                           '/service/sequence_editor/respond_form/entityreview': '_respondEntitySourceDetailsFormOp',

                           '/service/sequence_editor/download': '_respondDownloadOp',

                           '/service/sequence_editor/reload_data/start': '_reloadDataStartOp',
                           '/service/sequence_editor/reload_data/check': '_reloadDataCheckOp',
                           #
                           #
                           '/service/sequence_editor/load_new_sequence/start': '_loadDataSeqDbStartOp',
                           '/service/sequence_editor/load_new_sequence/check': '_loadDataCheckOp',
                           '/service/sequence_editor/load_new_taxid/start': '_loadDataTaxidStartOp',
                           '/service/sequence_editor/load_new_taxid/check': '_loadDataCheckOp',
                           #
                           '/service/sequence_editor/get_edit_timestamp': '_getEditTimeStampOp',
                           '/service/sequence_editor/get_alignment_timestamp': '_getAlignmentTimeStampOp',
                           '/service/sequence_editor/align_view': '_alignViewOp',
                           '/service/sequence_editor/align_view_frame': '_alignViewFrameOp',
                           '/service/sequence_editor/polymer_linkage_table': '_getPolymerLinkageTableOp'
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
                rC.setError(errMsg='Unknown operation')
            else:
                mth = getattr(self, self.__appPathD[reqPath], None)
                rC = mth()
            return rC
        except:
            traceback.print_exc(file=self.__lfh)
            self.__reqObj.setReturnFormat(return_format="json")
            rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            rC.setError(errMsg='Operation failure')
            return rC

    def setLogHandle(self, log=sys.stderr):
        """  Reset the stream for logging output.
        """
        try:
            self.__lfh = log
            return True
        except:
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
        """
        """
        return self.__alignView()

    def _alignViewFrameOp(self):
        """
        """
        return self.__alignView(frameOpt=True)

    def _saveOp(self):
        return self.__saveSelection()

    def _closeCompletedOp(self):
        """ RPS, 2010-Aug-02: created independent function for
            closing sequence editor app independent of saving selection
            RPS, 2010-Aug-12: renamed for explicit use for completion of seq edit process
        """
        return self.__closeSeqEditor(mode='completed')

    def _closeUnfinishedOp(self):
        """ RPS, 2010-Aug-12: created independent function for
            closing sequence editor app independent of saving selection
            *BUT* with intention of returning to complete editing activity
        """
        return self.__closeSeqEditor(mode='unfinished')

    def _saveSessionInfoOp(self):
        """ Save session information into data archive directory for debugging
        """
        self.__getSession()
        htmlTxt = self.__saveSessionInfo()
        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        rC.setHtmlText(htmlTxt)
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

    def _globalEditMenuOp(self):
        return self.__editAlignment()

    def _globalEditFormOp(self):
        return self.__editAlignment()

    def _undoEditOp(self):
        return self.__editAlignment()

    def _deleteOp(self):
        return self.__editAlignment()

    def _respondDownloadOp(self):
        """  Respond to download request for sequence assignment file - type = [model,assignment]
        """
        self.__getSession()
        identifier = self.__reqObj.getValue("identifier")
        cType = self.__reqObj.getValue("type")
        #

        if cType in ['model']:
            self.__reqObj.setValue('content_type', 'model')
        elif cType in ['assignment']:
            self.__reqObj.setValue('content_type', 'seq-assign')
        self.__reqObj.setValue('data_set_id', identifier)
        self.__reqObj.setValue('file_source', 'session')
        #
        wdu = WebDownloadUtils(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        rC = wdu.makeDownloadResponse()
        #
        if rC.isError():
            self.__reqObj.setReturnFormat(return_format="json")
        else:
            rC.setStatusCode('ok')
            self.__reqObj.setReturnFormat(return_format="binary")
        return rC

    def _respondTaxonomyFormOp(self):

        # Do this before updating the session information to avoid overwrite  -
        uds = UtilDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        seq_search_op = self.__reqObj.getValue("seq_search_op")
        uds.set('seq_search_op', seq_search_op)
        uds.serialize()
        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp._respondTaxonomyFormOp() - starting with seq_search_op %s\n" % seq_search_op)
        #
        self.__getSession()
        #
        sdU = UpdatePolymerEntityPartitions(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        status = sdU.polymerEntityPartEditFormResponder()
        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        if (status):
            rC.setStatusCode('ok')
        else:
            rC.setError(errMsg='Taxonomy form update processing failure')
        return rC

    def _loadTaxonomyFormOp(self):
        self.__getSession()
        authId = self.__reqObj.getValue("auth_id")
        identifier = self.__reqObj.getValue("identifier")
        if authId is not None and len(authId) > 1:
            tL = str(authId).split('_')
            entityId = tL[2]
        cD = self.__makeTaxonomyEditForm(entityId=entityId, entryId=identifier)
        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        rC.addDictionaryItems(cD=cD)
        return rC

    def __makeTaxonomyEditForm(self, entityId='1', entryId=''):
        rD = {}
        sdU = UpdatePolymerEntityPartitions(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        rD['htmlcontent'] = sdU.makePolymerEntityPartEditForm(entityId=entityId, entryId=entryId)
        return rD

    def _respondSeqDbRefFormOp(self):
        self.__getSession()
        sdU = UpdatePolymerEntityReference(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        packedRefSeqLabel = sdU.seqDbRefFormResponder()
        seqFetchError = sdU.getSeqFetchError()
        #
        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        if seqFetchError:
            rC.setError(errMsg=seqFetchError)
        elif (packedRefSeqLabel is not None):
            rC.setStatusCode('ok')
        else:
            rC.setError(errMsg='Reference sequence form update processing failure')

        return rC

    def _loadSeqDbRefFormOp(self):
        self.__getSession()
        refId = self.__reqObj.getValue("ref_id")
        identifier = self.__reqObj.getValue("identifier")
        if refId is not None and len(refId) > 1:
            tL = str(refId).split('_')
            entityId = str(tL[1])
            partId = str(tL[3])

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
            tL = str(authId).split('_')
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
        if (ok):
            rC.setStatusCode('ok')
        else:
            rC.setError(errMsg='Entity review form update processing failure')

        return rC

    def _launchJmolViewerOp(self):
        self.__getSession()
        if (self.__verbose):
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
        """ Loads pre-calculated HTML pages containing alignment summary --

            Performs a default selection of the best matching reference sequences.
        """
        return self.__loadSummary(op='load')

    def _reloadSummaryOp(self):
        """ Reloads pre-calculated HTML pages containing alignment summary -
        """
        return self.__loadSummary(op='reload')

    def _newSessionOp(self):
        """ Entry point for new sessions when invoked as the sequence editor tool.
        """
        return self.__newSession()

    def _getEditTimeStampOp(self):
        """Return the timestamp for the last edit.
        """
        self.__getSession()
        tS = self.__getTimeStamp(fileName='sequenceEditStore.pic')
        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        rC.setEditTimeStamp(tS)
        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp.__getEditTimeStamp() time stamp %s\n" % tS)
        return rC

    def _getAlignmentTimeStampOp(self):
        """Return the timestamp for the alignment store.
        """
        self.__getSession()
        tS = self.__getTimeStamp(fileName='alignmentDataStore.pic')
        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        rC.setAlignmentTimeStamp(tS)
        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp.__getAlignmentTimeStamp() time stamp %s\n" % tS)
        return rC

    def _newSessionWfOp(self):
        """ Entry point for new sessions when launched as a module from the wf environment.
        """
        return self.__newSessionWf()

    def _getPolymerLinkageTableOp(self):
        """Return the polymer linkage table
        """
        self.__getSession()
        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)

        sessionId = self.__reqObj.getSessionId()
        #sessionPath    =self.__reqObj.getSessionPath()
        siteId = self.__reqObj.getValue("WWPDB_SITE_ID")
        identifier = self.__reqObj.getValue("identifier")
        pI = PathInfo(siteId=siteId, sessionPath=self.__sessionPath, verbose=self.__verbose, log=self.__lfh)
        polyLinkHtmlPath = pI.getPolyLinkReportFilePath(identifier, fileSource='session')
        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp.__getPolymerLinkageTableOp() report file path %s\n" % polyLinkHtmlPath)
        #
        htmlList = []
        if os.access(polyLinkHtmlPath, os.R_OK):
            ifh = open(polyLinkHtmlPath, 'r')
            htmlList = ifh.readlines()
            ifh.close()
        else:
            htmlList = []

        (head, tail) = os.path.split(polyLinkHtmlPath)
        PathRel = os.path.join('/sessions', sessionId, tail)
        rC.setHtmlContentPath(PathRel)
        rC.setHtmlList(htmlList)
        return rC

    # --------------------------------------------------------------------------------------------------
    #  WF data loading and initialization entry point
    #
    def _loadDataWfStartOp(self):
        """ Load sequence data from model and sequence database matching results.  Compute preliminary
            alignment statistics.   Data is taken from the workflow storage system.
        """
        self.__getSession()
        fileSource = self.__reqObj.getValue("filesource")
        identifier = self.__reqObj.getValue("identifier")
        instance = self.__reqObj.getValue("instance")
        if fileSource in ['archive', 'wf-archive']:
            fileSource = 'archive'

        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp._loadDataWfStartOp() - starting with file source %s identifier %s instance %s\n" %
                             (fileSource, identifier, instance))
        rC = self.__loadDataStart(fileSource=fileSource)
        return rC

    # --------------------------------------------------------------------------------------------------
    #   Local data loading and initialization entry point methods ---
    #
    def _loadDataStartOp(self):
        """ Load sequence data from model and perform sequence database matching results.  Compute preliminary
            alignment statistics.   Model data can be uploaded (fileSource=local-upload) or copied from
            data processing repository (fileSource=local-repository)
        """
        self.__getSession()
        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp._loadLocalDataStartOp() - starting with session %s\n" % self.__sessionPath)

        wuu = WebUploadUtils(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        if wuu.isFileUpload():
            fileSource = 'local-upload'
            fileType = self.__reqObj.getValue("filetype")
            if (self.__verbose):
                self.__lfh.write("+SeqModWebApp._loadLocalDataStartOp() - file source %s type %s\n" % (fileSource, fileType))

            modelFileName = wuu.copyToSession(fileTag='file')
            #
            if modelFileName is None:
                rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
                rC.setError(errMsg='Upload processing failure')
                self.__reqObj.setReturnFormat(return_format="json")
                return rC
            else:
                #
                fId, idType = wuu.perceiveIdentifier(modelFileName)
                self.__reqObj.setValue('identifier', fId)
                self.__reqObj.setValue('UploadFileName', modelFileName)
                self.__reqObj.setValue('UploadFileType', fileType)

                uds = UtilDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
                uds.set('identifier', fId)
                uds.set('UploadFileName', modelFileName)
                uds.set('UploadFileType', fileType)
                uds.serialize()
        else:
            fileSource = 'local-repository'
            identifier = self.__reqObj.getValue("identifier")
            saved_session_id = self.__reqObj.getValue("saved_session_id")
            ok = self.__checkSavedSessionInfo(identifier, saved_session_id)
            if ok:
                self.__reqObj.setValue('activegroupid', '1')
                rC = self.__loadSummary()
                rC.setIdentifier(identifier)
                rC.setStatusCode('ok')
                return rC
            #
            uds = UtilDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            uds.set('identifier', identifier)
            uds.serialize()

            if self.__verbose:
                self.__lfh.write("+SeqModWebApp._loadLocalDataStartOp() - file source %s identifier %s\n" % (fileSource, identifier))

        rC = self.__loadDataStart(fileSource=fileSource)
        return rC

    # --------------------------------------------------------------------------------------------------
    #   Re-run blast search entry point methods ---
    #
    def _reRunBlastOp(self):
        """ Re-run blast search, perform sequence database matching results, compute preliminary alignment statistics.
        """
        self.__getSession()
        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp._reRunBlastOp() - starting with session %s\n" % self.__sessionPath)
        #

        dU = DetachUtils(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        dI = DataImporter(reqObj=self.__reqObj, fileSource="session", verbose=self.__verbose, log=self.__lfh)
        dU.set(workerObj=dI, workerMethod="runBlastSearch")
        dU.runDetach()

        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        identifier = self.__reqObj.getValue('identifier')
        rC.setIdentifier(identifier)
        rC.setStatusCode('running')

        return rC

    # --------------------------------------------------------------------------------------------------
    #  Add a single reference sequence from an sequence database.
    #
    def _loadDataSeqDbStartOp(self):
        """ Load a single reference database sequence.  Compute preliminary
            alignment statistics.   Data is loaded via a webservice from the sequence data.
        """

        fileSource = "sequence-database"

        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp._loadDataSeqDb() - starting with file source %s\n" % fileSource)

        return self.__loadDataStart(fileSource=fileSource)

    # --------------------------------------------------------------------------------------------------
    #  Change the Taxid for reference sequence database search/sorting
    #
    def _loadDataTaxidStartOp(self):
        """ Reassign the Taxid for the selected entity and re-sort the reference sequences  Compute preliminary
            alignment statistics.   No data is loaded here.
        """

        fileSource = "input-taxid"

        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp._loadDataTaxidOp() - starting with file source %s\n" % fileSource)

        return self.__loadDataStart(fileSource=fileSource)

    def _loadDataCheckOp(self):
        """ Shared checking method for data initialization and loading/uploading operations.
        """
        return self.__loadDataCheck()

    # --------------------------------------------------------------------------------------------
    #                              Supporting methods -
    #
    def __getTimeStamp(self, fileName='alignmentDataStore.pic'):

        sessionId = self.__reqObj.getSessionId()
        sessionPath = self.__reqObj.getSessionPath()
        fPathAbs = os.path.join(sessionPath, sessionId, fileName)
        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp.__getTimeStamp() - checking %s in path %s\n" % (fileName, fPathAbs))
        try:
            mTime = os.path.getmtime(fPathAbs)
        except:
            mTime = 0.0
        return mTime

    def __newSession(self):
        #
        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp.__newSession() - starting\n")

        self.__getSession(forceNew=True)
        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        return rC

    def __getSession(self, forceNew=False, useContext=True):
        """ Join existing session or create new session as required.

            Import any saved parameters into the current request object from UtilDataStore().
        """
        #
        sessionId = self.__reqObj.getSessionId()
        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp.__getSession() - starting with session %s\n" % sessionId)
        #
        self.__sObj = self.__reqObj.newSessionObj(forceNew=forceNew)

        self.__sessionId = self.__sObj.getId()
        self.__sessionPath = self.__sObj.getPath()
        self.__rltvSessionPath = self.__sObj.getRelativePath()
        #
        if (useContext):
            uds = UtilDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            dd = uds.getDictionary()
            self.__lfh.write("+SeqModWebApp.__getSession() - importing persisted session parameters: %r\n" % dd.items())
            self.__reqObj.setDictionary(dd, overWrite=True)
        #
        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp.__getSession() - Leaving with session path %s\n" % self.__sessionPath)

    def __saveSessionInfo(self):
        """ Copy session directory under entry data archive directory
        """
        identifier = self.__reqObj.getValue("identifier")
        fileList = self.__getSessionFileList(identifier, self.__sessionPath)
        #
        htmlTxt = "Saving entry identifier=" + identifier + " session ID=" + self.__sessionId + " information failed."
        if not fileList:
            return htmlTxt
        #
        siteId = self.__reqObj.getValue("WWPDB_SITE_ID")
        pI = PathInfo(siteId=siteId, sessionPath=self.__sessionPath, verbose=self.__verbose, log=self.__lfh)
        modelFile = pI.getFilePath(dataSetId=identifier, contentType='model', formatType='pdbx')
        if modelFile:
            dirPath = os.path.dirname(modelFile)
            saveSessionPath = os.path.join(dirPath, self.__sessionId)
            if os.access(saveSessionPath, os.F_OK):
                shutil.rmtree(saveSessionPath)
            #
            os.makedirs(saveSessionPath)
            if os.access(saveSessionPath, os.F_OK):
                successList = []
                for filename in fileList:
                    sourceFile = os.path.join(self.__sessionPath, filename)
                    targetFile = os.path.join(saveSessionPath, filename)
                    shutil.copyfile(sourceFile, targetFile)
                    if os.access(targetFile, os.F_OK) and filecmp.cmp(sourceFile, targetFile):
                        successList.append(filename)
                    #
                #
                if len(fileList) == len(successList):
                    htmlTxt = "Saved entry identifier=" + identifier + " session ID=" + self.__sessionId + " information for debugging."
                #
            #
        #
        return htmlTxt

    def __checkSavedSessionInfo(self, identifier, saved_session_id):
        """ copy saved session directory information to current session directory
        """
        if (not identifier) or (not saved_session_id):
            return False
        #
        siteId = self.__reqObj.getValue("WWPDB_SITE_ID")
        pI = PathInfo(siteId=siteId, sessionPath=self.__sessionPath, verbose=self.__verbose, log=self.__lfh)
        modelFile = pI.getFilePath(dataSetId=identifier, contentType='model', formatType='pdbx')
        if not modelFile:
            return False
        #
        dirPath = os.path.dirname(modelFile)
        saveSessionPath = os.path.join(dirPath, saved_session_id)
        if not os.access(saveSessionPath, os.F_OK):
            return False
        #
        fileList = self.__getSessionFileList(identifier, saveSessionPath)
        if not fileList:
            return False
        #
        for filename in fileList:
            sourceFile = os.path.join(saveSessionPath, filename)
            targetFile = os.path.join(self.__sessionPath, filename)
            shutil.copyfile(sourceFile, targetFile)
        #
        return True

    def __getSessionFileList(self, identifier, sessionPath):
        """
        """
        generalFlag = False
        currentFlag = False
        entryFlag = False
        fileList = []
        os.chdir(sessionPath)
        for filename in os.listdir('.'):
            if filename.startswith(identifier):
                entryFlag = True
                fileList.append(filename)
            elif filename == 'general-util-session.pic':
                generalFlag = True
                fileList.append(filename)
            elif filename == 'current-alignment-summary.html':
                currentFlag = True
                fileList.append(filename)
            #
        #
        if (not generalFlag) or (not currentFlag) or (not entryFlag):
            return []
        #
        return fileList

    def __loadDataStart(self, fileSource='local-repository'):
        """ Launch a child process to handle data initialization and loading operations.

            This method supports both local and wf data/file sources.
        """
        #
        siteId = self.__reqObj.getValue('WWPDB_SITE_ID')
        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp.__loadDataStart() - STARTING using file source %s at site %s\n" % (fileSource, siteId))

        dU = DetachUtils(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)

        dI = DataImporter(reqObj=self.__reqObj, fileSource=fileSource, verbose=self.__verbose, log=self.__lfh)

        dU.set(workerObj=dI, workerMethod="loadData")

        dU.runDetach()

        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        identifier = self.__reqObj.getValue('identifier')
        rC.setIdentifier(identifier)
        rC.setStatusCode('running')

        return rC

    def __loadDataCheck(self):
        """Performs a check on the contents of a semaphore file and returns the associated status.

           This method currently supports both rcsb and wf filesources.
        """
        #
        self.__getSession(useContext=True)
        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp.__loadDataCheck() - starting\n")

        sph = self.__reqObj.getSemaphore()
        delayValue = self.__reqObj.getValue("delay")
        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp.__loadDataCheck() Checking status of semaphore %s with delay %s\n" % (sph, delayValue))

        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        identifier = self.__reqObj.getValue('identifier')
        rC.setIdentifier(identifier)

        dU = DetachUtils(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        if (dU.semaphoreExists(sph)):
            status = dU.getSemaphore(sph)
            if (self.__verbose):
                self.__lfh.write("+SeqModWebApp.__loadDataCheck() status value for semaphore %s is %s\n" % (sph, str(status)))
            if (status == "OK"):
                rC.setStatusCode('completed')
            else:
                rC.setStatusCode('failed')
        else:
            if (self.__verbose):
                self.__lfh.write("+SeqModWebApp.__loadDataCheck() semaphore %s pending - waiting %s\n" % (sph, delayValue))
            time.sleep(int(delayValue))
            rC.setStatusCode('running')

        return rC

    def __loadSummary(self, op='load'):
        #
        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp.__loadSummary() - starting\n")

        self.__getSession(useContext=True)
        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        try:
            sV = SummaryView(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            sumObj,warningMsg = sV.loadSummary(operation=op)
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
            fp = open(htmlPathAbs, 'w')
            fp.write("%s" % ''.join(oL))
            fp.close()
            rC.setHtmlContentPath(os.path.join("/sessions", sessionId, "current-alignment-summary.html"))
            #
            rC.setIdentifier(self.__reqObj.getValue("identifier"))
            rC.setStructTitle(eD['STRUCT_TITLE'])
            rC.setCitationTitle(eD['CITATION_TITLE'])
            rC.setPdbCode(eD['PDB_ID'])
            rC.setGroupIdList(gIdList)
            #
            #  Reset --
            #
            self.__reqObj.setNewRefId('')
        except:
            rC.setError(errMsg="Sequence alignment summary preparation has failed.")
            rC.setStatusCode('failed')
            if (self.__verbose):
                self.__lfh.write("+SeqModWebApp.__loadSummary() - failed with exception\n")
                traceback.print_exc(file=self.__lfh)
        #
        return rC

    def __storeAlignmentStart(self):
        """ Launch a subprocess to compute, store and render the aligned input sequence list.
        """
        #
        self.__getSession(useContext=True)
        sessionId = self.__reqObj.getSessionId()
        alignTag = self.__reqObj.getValue("aligntag")
        operation = self.__reqObj.getValue("operation")
        viewalign_order = self.__reqObj.getValue("viewalign_order")
        alignids = self.__reqObj.getValue("alignids")
        activeGroupId = self.__reqObj.getValue("activegroupid")

        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp.__storeAlignmentStart() - starting\n")
            self.__lfh.write("+SeqModWebApp.__storeAlignmentStart() operation       %s\n" % operation)
            self.__lfh.write("+SeqModWebApp.__storeAlignmentStart() aligntag        %s\n" % alignTag)
            self.__lfh.write("+SeqModWebApp.__storeAlignmentStart() activegroupid   %s\n" % activeGroupId)
            self.__lfh.write("+SeqModWebApp.__storeAlignmentStart() alignids        %s\n" % alignids)
            self.__lfh.write("+SeqModWebApp.__storeAlignmentStart() viewalign_order %s\n" % viewalign_order)

        if (len(alignTag) < 1):
            #sLabU = SequenceLabelUtils()
            #alignTag =  sLabU.getAlignGroupId(self.__reqObj.getAlignIdList()  )
            #alignTag = str(time.strftime("%Y%m%d%H%M%S", time.localtime()))
            alignTag = activeGroupId
            if (self.__verbose):
                self.__lfh.write("+SeqModWebApp.__storeAlignmentStart() creating new alignment tag %s\n" % alignTag)

        self.__reqObj.setValue('aligntag', alignTag)
        #

        #

        dU = DetachUtils(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        dU.set(workerObj=self, workerMethod="_storeAlignmentFiles")
        dU.runDetach()

        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        rC.setStatusCode('running')
        rC.setAlignTag(alignTag)
        rC.setViewAlignOrder(viewalign_order)
        # reset edit operation -
        rC.setEditOp(0)
        self.__lfh.write("+SeqModWebApp.__storeAlignmentStart() Parent process completed\n")
        return rC

    def __storeAlignmentCheck(self):
        #
        self.__getSession()
        if (self.__verbose):
            self.__lfh.write("\n+SeqModWebApp.__storeAlignmentCheck() - starting\n")
        #
        sph = self.__reqObj.getSemaphore()
        sessionId = self.__reqObj.getSessionId()
        alignTag = self.__reqObj.getValue("aligntag")
        delayValue = self.__reqObj.getValue("delay")
        activeGroupId = self.__reqObj.getValue("activegroupid")
        #
        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp.__storeAlignmentCheck() Checking status of semaphore %s with delay %s\n" % (sph, delayValue))
        #
        #
        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        dU = DetachUtils(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)

        if (dU.semaphoreExists(sph)):
            status = dU.getSemaphore(sph)
            if (self.__verbose):
                self.__lfh.write("+SeqModWebApp.__storeAlignmentCheck() status value for semaphore %s is %s\n" % (sph, str(status)))

            if (status == "OK"):
                # reset edit operation -
                rC.setEditOp(0)
                rC.setAlignTag(alignTag)
                rC.setIdentifier(self.__reqObj.getValue("identifier"))
                rC.setStructTitle(self.__reqObj.getValue("title"))
                rC.setPdbCode(self.__reqObj.getValue("pdbid"))
                #
                rC.setHtmlContentPath(os.path.join("/sessions", sessionId, "current-alignment-" + alignTag + ".html"))
                rC.setStatusCode('completed')
                #
                sPath = self.__reqObj.getSessionPath()
                rptPathAbs = os.path.join(sPath, sessionId, "conflict-report-" + alignTag + ".html")
                if (os.access(rptPathAbs, os.F_OK)):
                    pth = os.path.join("/sessions", sessionId, "conflict-report-" + alignTag + ".html")
                    rC.setConflictReportPath(pth)
                    rC.setConflictReportFlag(True)
                    if (self.__verbose):
                        self.__lfh.write("+SeqModWebApp.__storeAlignmentCheck() conflict report %s\n" % pth)
                else:
                    rC.setConflictReportFlag(False)

                rptPathAbs = os.path.join(sPath, sessionId, "current-alignment-annotation-" + alignTag + ".html")
                if (os.access(rptPathAbs, os.F_OK)):
                    pth = os.path.join("/sessions", sessionId, "current-alignment-annotation-" + alignTag + ".html")
                    rC.setAnnotationReportPath(pth)
                    rC.setAnnotationReportFlag(True)
                else:
                    rC.setAnnotationReportFlag(False)

                rptPathAbs = os.path.join(sPath, sessionId, "warning-alignment-" + alignTag + ".html")
                if (os.access(rptPathAbs, os.F_OK)):
                    fp = open(rptPathAbs, 'r')
                    warningMsg = fp.read()
                    fp.close()
                    rC.setWarning(warnMsg=warningMsg)
                #

                misMatchTypeFile = os.path.join(sPath, sessionId, "mismatch-type.txt")
                if (os.access(misMatchTypeFile, os.F_OK)):
                    fp = open(misMatchTypeFile, 'r')
                    misMatchType = fp.read()
                    fp.close()
                    rC.addDictionaryItems({ 'gedittype' : misMatchType })
                else:
                    rC.addDictionaryItems({ 'gedittype' : 'no-mismatch' })
                #

                # --------------------------------------------
                if False:
                    aU = AlignmentUtils(self.__reqObj, verbose=self.__verbose, log=self.__lfh)
                    idList = aU.readAlignmentIdList()
                    rC.setAlignIdList(alignIdList=idList)
                # ---------------------------------------------

            else:
                rC.setStatusCode('failed')
        else:
            if (self.__verbose):
                self.__lfh.write("+SeqModWebApp.__storeAlignmentCheck() semaphore %s pending - waiting %s\n" % (sph, delayValue))
            time.sleep(int(delayValue))
            rC.setStatusCode('running')
        #
        self.__lfh.write("%s" % (''.join(rC.dump())))
        return rC

    def _storeAlignmentFiles(self):
        """Compute and store the rendered alignment for the input sequence list.
        """
        sessionId = self.__reqObj.getSessionId()
        alignTag = self.__reqObj.getValue("aligntag")
        alignViewOrder = self.__reqObj.getAlignmentOrdering()
        operation = self.__reqObj.getValue("operation")
        identifier = self.__reqObj.getValue("identifier")
        #
        alignIdList = self.__reqObj.getAlignList()
        # alignIdList=self.__reqObj.getSummaryAlignList()
        selectedIdList = self.__reqObj.getSummarySelectList()
        alignGroupId = self.__reqObj.getValue("activegroupid")
        #
        if (self.__verbose):
            self.__lfh.write("\n\n+SeqModWebApp._storeAlignmentFiles() - starting\n")
            self.__lfh.write("+SeqModWebApp._storeAlignmentFiles() - sessionId %s\n" % sessionId)
            self.__lfh.write("+SeqModWebApp._storeAlignmentFiles() - using alignment tag %s\n" % alignTag)
            self.__lfh.write("+SeqModWebApp._storeAlignmentFiles() - operation %s\n" % operation)
            self.__lfh.write("+SeqModWebApp._storeAlignmentFiles() - identifier %s\n" % identifier)
            self.__lfh.write("+SeqModWebApp._storeAlignmentFiles() - active/alignGroupId  %s\n" % alignGroupId)
            self.__lfh.write("+SeqModWebApp._storeAlignmentFiles() - alignIdList %s\n" % alignIdList)
            self.__lfh.write("+SeqModWebApp._storeAlignmentFiles() - selectIdList %s\n" % selectedIdList)
        #
        map = MultiAlignPseudo(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        alignSeqList = map.run(operation=operation, identifier=identifier, alignIdList=alignIdList, alignGroupId=alignGroupId, selectedIdList=selectedIdList)
        aD = map.getSelectedAnnotations()
        eD = map.getEntryDetails()
        #
        # Update stored select and align id list with any revisions from previous alignment.
        #
        alignIdList = map.updateSequenceIdList(alignIdList)
        selectedIdList = map.updateSequenceIdList(selectedIdList)
        #
        if self.__verbose:
            self.__lfh.write("+SeqModWebApp._storeAlignmentFiles() - revised alignIdList %s\n" % alignIdList)
            self.__lfh.write("+SeqModWebApp._storeAlignmentFiles() - revised selectIdList %s\n" % selectedIdList)

        #
        uds = UtilDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        uds.set('title', eD['STRUCT_TITLE'])
        uds.set('pdbid', eD['PDB_ID'])
        uds.set('rev-allalignids', alignIdList)
        uds.set('rev-selectids', selectedIdList)
        uds.serialize()
        #
        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp.__storeAlignmentFiles() - align entity group %s length align sequence list %d\n" %
                             (alignGroupId, len(alignSeqList)))

        avd = AlignmentViewDepiction(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        oL = avd.renderAlignment(alignSeqList=alignSeqList, alignGroupId=alignGroupId, type='original', viewOrderCode=alignViewOrder)
        cL = avd.renderConflictTable(alignSeqList=alignSeqList, alignGroupId=alignGroupId, type='original')
        aL = avd.renderAnnotations(aD=aD)
        #
        sPath = self.__reqObj.getSessionPath()
        htmlPathAbs = os.path.join(sPath, sessionId, "warning-alignment-" + alignTag + ".html")
        if (os.access(htmlPathAbs, os.F_OK)):
             os.remove(htmlPathAbs)
        #
        misMatchTypeFile = os.path.join(sPath, sessionId, "mismatch-type.txt")
        fp = open(misMatchTypeFile, 'w')
        fp.write("%s" % avd.getMisMatchType())
        fp.close()
        #
        for aTup in alignSeqList:
            seqLabel = aTup[1]
            (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = seqLabel.get()
            if seqType == 'auth' and aD['source_method'].upper() == 'NAT':
                eelCommentL = []
                refSeq = aTup[2]
                for idx in range(0, len(refSeq)):
                    if (len(refSeq[idx][5]) > 0):
                        comment = str(refSeq[idx][5])
                        for cType in ( 'engineered mutation', 'expression tag', 'linker' ):
                            if (comment.find(cType) != -1) and (cType not in eelCommentL):
                                eelCommentL.append(cType)
                            #
                        #
                    #
                #
                if len(eelCommentL) > 0:
                    warningMsg = "Entity '" + alignGroupId + "' with natural source has '" + "', '".join(eelCommentL) + "' sequence annotation information.<br />"
                    fp = open(htmlPathAbs, 'w')
                    fp.write("%s\n" % warningMsg)
                    fp.close()
                #
            #
        #
        htmlPathAbs = os.path.join(sPath, sessionId, "current-alignment-" + alignTag + ".html")
        fp = open(htmlPathAbs, 'w')
        fp.write("%s" % ''.join(oL))
        fp.close()
        #
        rptPathAbs = os.path.join(sPath, sessionId, "conflict-report-" + alignTag + ".html")
        if (len(cL) > 0):
            fp = open(rptPathAbs, 'w')
            fp.write("%s" % ''.join(cL))
            fp.close()
        else:
            try:
                os.remove(rptPathAbs)
            except:
                pass
            #
        #
        htmlPathAbs = os.path.join(sPath, sessionId, "current-alignment-annotation-" + alignTag + ".html")
        if len(aL) > 0:
            fp = open(htmlPathAbs, 'w')
            fp.write("%s" % ''.join(aL))
            fp.close()
        else:
            try:
                os.remove(htmlPathAbs)
            except:
                pass
            #
        #
        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp._storeAlignmentFiles() - completed \n")
        #
        return True

    def _reloadDataWorker(self):
        seqSearchOp = self.__reqObj.getValue('seq_search_op')
        entityId = self.__reqObj.getValue('activegroupid')

        if seqSearchOp == 'on':
            if (self.__verbose):
                self.__lfh.write("\n\n+SeqModWebApp._reloadDataWorker() - starting with entity %s and searchOp %s\n" % (entityId, seqSearchOp))
            try:
                rsdU = ReferenceSequenceDataUpdate(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
                ok1 = rsdU.doUpdate(entityId=entityId)
                ok2 = rsdU.doUpdateSelections(entityId=entityId)
            except:
                if (self.__verbose):
                    self.__lfh.write("+SeqModWebApp._reloadDataWorker() - failing for entity %s\n" % entityId)
                    traceback.print_exc(file=self.__lfh)
                return False
        else:
            if (self.__verbose):
                self.__lfh.write("+SeqModWebApp._reloadDataWorker() - skipping search with entity %s and searchOp %s\n" % (entityId, seqSearchOp))
            try:
                rsdU = ReferenceSequenceDataUpdate(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
                ok = rsdU.doUpdateSelections(entityId=entityId)
            except:
                if (self.__verbose):
                    self.__lfh.write("+SeqModWebApp._reloadDataWorker() - failing for entity %s\n" % entityId)
                    traceback.print_exc(file=self.__lfh)
                return False

        return True

    def _reloadDataStartOp(self):
        """ Launch a child process to handle reloading operations required by taxonomy/part edit operations.

            This method supports both local and wf data/file sources.
        """
        #
        self.__getSession()
        siteId = self.__reqObj.getValue('WWPDB_SITE_ID')
        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp._reloadDataStartOp() - STARTING at site %s\n" % siteId)

        dU = DetachUtils(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)

        dU.set(workerObj=self, workerMethod="_reloadDataWorker")

        dU.runDetach()

        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        identifier = self.__reqObj.getValue('identifier')
        rC.setIdentifier(identifier)
        rC.setStatusCode('running')

        return rC

    def _reloadDataCheckOp(self):
        """Performs a check on the contents of a semaphore file and returns the associated status.

           This method currently supports both rcsb and wf filesources.
        """
        #
        self.__getSession(useContext=True)
        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp._reloadDataCheckOp() - starting\n")

        sph = self.__reqObj.getSemaphore()
        delayValue = self.__reqObj.getValue("delay")
        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp._reloadDataCheckOp() Checking status of semaphore %s with delay %s\n" % (sph, delayValue))

        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        identifier = self.__reqObj.getValue('identifier')
        rC.setIdentifier(identifier)

        dU = DetachUtils(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        if (dU.semaphoreExists(sph)):
            status = dU.getSemaphore(sph)
            if (self.__verbose):
                self.__lfh.write("+SeqModWebApp._reloadDataCheckOp() status value for semaphore %s is %s\n" % (sph, str(status)))
            if (status == "OK"):
                # rC=self.__loadSummary(op='reload')
                rC.setStatusCode('completed')
            else:
                rC.setStatusCode('failed')
        else:
            if (self.__verbose):
                self.__lfh.write("+SeqModWebApp._reloadDataCheckOp() semaphore %s pending - waiting %s\n" % (sph, delayValue))
            time.sleep(int(delayValue))
            rC.setStatusCode('running')

        return rC

    def __editAlignment(self):

        sessionId = self.__reqObj.getSessionId()
        alignTag = self.__reqObj.getValue("aligntag")
        #
        operation = self.__reqObj.getValue("operation")
        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp.__editAlignment() - starting\n")
            self.__lfh.write("+SeqModWebApp.__editAlignment() - session id %s\n" % sessionId)
            self.__lfh.write("+SeqModWebApp.__editAlignment() - operation %s\n" % operation)
            self.__lfh.write("+SeqModWebApp.__editAlignment() - alignment tag %s\n" % alignTag)

        #
        aE = AlignmentEdit(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        aE.setAlignmentTag(alignTag)
        #
        cD = aE.edit()
        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp.__editAlignment() edit completed for session %s\n" % sessionId)
        #
        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        rC.addDictionaryItems(cD=cD)
        return rC

    def __saveSelection(self):
        """ Save assignment and annotation file
        """
        self.__getSession()
        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp.__saveSelection() - starting\n")

        selectIdList = self.__reqObj.getSummarySelectList()
        sEx = SequenceDataExport(reqObj=self.__reqObj, exportList=selectIdList, verbose=self.__verbose, log=self.__lfh)
        ok, numConflicts, conflictList, warningMsg = sEx.exportAssignments()
        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp.__saveSelection() - exporting assignments status = %r and conflict count  %d\n" % (ok, numConflicts))
        # if (ok and (numConflicts==0)):
        if (ok and (numConflicts < 3)):
            ok = sEx.applyAssignmentsToModel()
            if (self.__verbose):
                self.__lfh.write("+SeqModWebApp.__saveSelection() - exporting model status = %r\n" % ok)

        self.__reqObj.setReturnFormat(return_format="json")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        if warningMsg:
            rC.setWarning(warnMsg=warningMsg)
        #
        if (not ok):
            rC.setError(errMsg='Sequence data export failed')
        elif numConflicts > 0:
            oL = []
            for ctup in conflictList:
                oL.append("Entity %s(%d)" % (ctup[0], ctup[1]))
            tS = ' , '.join(oL)
            rC.setError(errMsg='Warning total sample conflict count = %d  - %s' % (numConflicts, tS))
        return rC

    def __closeSeqEditor(self, mode):
        """ RPS, 2010-Aug-02: providing dedicated function to accommodate user request to end sesssion when
            done with using sequence editor interface
            RPS, 2010-Aug-12: refactored to support different 'modes' = ('completed' | 'unfinished')

        """
        self.__getSession()
        if (self.__verbose):
            self.__lfh.write("\n\n+SeqModWebApp.__closeSeqEditor() - starting with operation %s\n" % mode)

        if (mode == 'completed'):
            state = "closed(0)"
        elif (mode == 'unfinished'):
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

        if (self.__verbose):
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

        if fileSource in ['archive', 'wf-archive']:
            # copy session files back to the archive directory -- independent of close mode --
            dI = DataImporter(reqObj=self.__reqObj, fileSource=fileSource, verbose=self.__verbose, log=self.__lfh)
            ok = dI.copyFilesOnClose()

        if fileSource in ['wf-instance']:
            # Copy session files back to the invoking instance directory -
            dI = DataImporter(reqObj=self.__reqObj, fileSource=fileSource, verbose=self.__verbose, log=self.__lfh)
            ok = dI.copyFilesOnClose(includePolyLinkFile=True)
            # And if we are done then also copy back to the archive --
            ok1 = True
            if mode in ['completed']:
                ok1 = dI.copyFiles(inputFileSource="session",outputFileSource='wf-archive',includeModelFile=True,versionIndex=4, \
                                   messageHead="DataImporter.copyFilesOnClose()")
            #
            ok = ok and ok1

        okDb = True
        if self.__doWfTracking:
            try:
                okDb = self.__updateWfStatus(depId=depId, instId=instId, classId=classId, statusCode=state)
            except:
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

        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp.__newSessionWf() - starting\n")

        sessionId = self.__reqObj.getSessionId()
        if (len(sessionId) == 0):
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
        myD['sessionid'] = sessionId
        myD['identifier'] = identifier
        myD['instance'] = instance
        myD['classid'] = self.__reqObj.getValue("classID")
        myD['filesource'] = fileSource

        myD['height'] = '100%'
        # myD['height']='1000px'
        myD['width'] = '100%'
        #
        uds = UtilDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        uds.set('identifier', identifier)
        uds.set('filesource', fileSource)
        uds.set('instance', instance)
        uds.set('skipstatus', skipStatus)
        uds.serialize()
        try:
            ok = self.__updateWfStatus(depId=myD['identifier'], instId=myD['instance'], classId=myD['classid'], statusCode="open")
        except:
            ok = False
        #
        self.__reqObj.setReturnFormat(return_format="html")
        rC = SeqModResponseContent(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        templateFilePath = os.path.join(self.__reqObj.getValue("TemplatePath"), "summary_template.html")
        webIncludePath = os.path.join(self.__reqObj.getValue("TopPath"), "htdocs")
        rC.setHtmlTextFromTemplate(templateFilePath=templateFilePath, webIncludePath=webIncludePath, parameterDict=myD)
        return rC

    def __alignView(self, frameOpt=True):
        """ Send a skeleton HTML template with placeholders for the sequence alignment and conflict report.
        """
        viewalign_order = self.__reqObj.getValue("viewalign_order")
        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp.__alignView() - starting with %s \n" % viewalign_order)
        #
        self.__getSession(forceNew=False, useContext=True)

        sessionId = self.__reqObj.getSessionId()
        alignTag = self.__reqObj.getValue("aligntag")
        viewalign_order = self.__reqObj.getValue("viewalign_order")
        activeEntityGroupId = self.__reqObj.getValue("activegroupid")

        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp.__alignView() Loading alignment template for sessionId %s alignment tag %s order %s\n"
                             % (sessionId, alignTag, viewalign_order))
        myD = {}
        myD['sessionid'] = self.__reqObj.getSessionId()
        myD['seqview'] = 'fromlist'
        myD['alignids'] = ','.join(self.__reqObj.getAlignList())
        myD['selectids'] = ','.join(self.__reqObj.getSummarySelectList())
        #myD['aligntag']  = str(time.strftime("%Y%m%d%H%M%S", time.localtime()))
        myD['aligntag'] = activeEntityGroupId
        myD['activegroupid'] = activeEntityGroupId

        myD['identifier'] = self.__reqObj.getValue("identifier")
        myD['instance'] = self.__reqObj.getValue("instance")
        myD['viewalign_order'] = self.__reqObj.getValue("viewalign_order")
        # initial view -
        myD['modelfilename'] = self.__reqObj.getValue("modelfilename")

        if (self.__verbose):
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
        if (self.__verbose):
            self.__lfh.write("+SeqModWebApp.__alignView() Just set viewalign_order in ResponseContent as %s\n" % viewalign_order)

        return rC

    def __updateWfStatus(self, depId, instId, classId, statusCode):
        """ Update progress and tracking status --
        """
        if not self.__doWfTracking:
            if (self.__verbose):
                self.__lfh.write("+SeqModWebApp.__updateWfStatus() - skipped setting depId %s instId %s classId %s status %s\n"
                                 % (depId, instId, classId, statusCode))
            return True
        try:
            if (self.__verbose):
                self.__lfh.write("+SeqModWebApp.__updateWfStatus() - setting depId %s instId %s classId %s status %s\n"
                                 % (depId, instId, classId, statusCode))
            wft = WfTracking(verbose=self.__verbose, log=self.__lfh)
            ok = wft.setInstanceStatus(depId=depId, instId=instId, classId=classId, status=statusCode)
            return ok
        except:
            if (self.__verbose):
                self.__lfh.write("+SeqModWebApp.__updateWfStatus() - failed depId %s instId %s classId %s status %s\n"
                                 % (depId, instId, classId, statusCode))
            traceback.print_exc(file=self.__lfh)
        return False


if __name__ == '__main__':
    sTool = SeqModWebApp()
    d = sTool.doOp()
    for k, v in d.items():
        sys.stdout.write("Key - %s  value - %r\n" % (k, v))
