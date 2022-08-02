##
# File:    UpdatePolymerEntityReference.py
# Date:    22-Mar-2014
#
# Updates:
# 04-July-2014  jdw Fix GB reference sequence update
# 29-July-2014  jdw Extend for sequence isoforms
# 21-Nov.-2018  zf  Re-organize the code: move fetching reference DB sequence process into FetchReferenceSequenceUtils class
#                   and aligning author/reference sequences process into AlignmentTools class
##
"""
Utilities for updating reference sequence of matching residue ranges

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.08"

import sys
import traceback

from wwpdb.apps.seqmodule.align.AlignmentTools import AlignmentTools
from wwpdb.apps.seqmodule.util.FetchReferenceSequenceUtils import FetchReferenceSequenceUtils
from wwpdb.apps.seqmodule.util.SequenceLabel import SequenceLabel
from wwpdb.apps.seqmodule.util.SequenceReferenceData import SequenceReferenceData
from wwpdb.apps.seqmodule.util.UpdateSequenceDataStoreUtils import UpdateSequenceDataStoreUtils


class UpdatePolymerEntityReference(UpdateSequenceDataStoreUtils):
    """Utilities for updating reference sequence of matching residue ranges"""

    def __init__(self, reqObj=None, verbose=False, log=sys.stderr):
        super(UpdatePolymerEntityReference, self).__init__(reqObj=reqObj, seqDataStore=None, verbose=verbose, log=log)
        self._reqObj = reqObj
        self._verbose = verbose
        self._lfh = log
        #
        self.__srd = SequenceReferenceData(verbose=self._verbose, log=self._lfh)
        self.__setup()

    def __setup(self):
        try:
            self.__placeHolderValue = "click-to-edit"
            self.__siteId = self._reqObj.getValue("WWPDB_SITE_ID")
            self.__sessionId = self._reqObj.getSessionId()
            self.__sessionObj = self._reqObj.getSessionObj()
            # self.__sessionPath = self.__sessionObj.getPath()
            #
            self.__selectIdList = self._reqObj.getSummarySelectList()
            #
        except:  # noqa: E722 pylint: disable=bare-except
            if self._verbose:
                self._lfh.write("+UpdatePolymerEntityReference.__setup() sessionId %s failed\n" % (self.__sessionObj.getId()))

    def makeSeqdbrefEditForm(self, entityId, partId=1, entryId=""):
        """Return a preliminary form to input sequence database references."""
        #
        form_template = """
        <div id="sectseqdbref">
        <h3>Reference Sequence Database Data Form for Entry %(entryid)s Entity %(entityid)s Part %(partid)s</h3>
        <form name="formseqdbref_%(partid)s" id="formseqdbref_%(partid)s" action="/service/sequence_editor/respond_form/seqdbref" method="post" class="auth_ajaxform">
            <input type="hidden" name="sessionid" value="%(sessionid)s" />
            <input type="hidden" name="identifier" value="%(entryid)s" />
            <input type="hidden" name="selected_auth_id" value="%(selected_auth_id)s" />
            <input type="hidden" name="partid" value="%(partid)s" />
            <input type="hidden" name="entityid" value="%(entityid)s" />
            <table>
            <tr>
               <th>Reference Resource</th>
               <th>Accession Code</th>
               <th>Seq Begin</th>
               <th>Seq End</th>
             </tr>
            <tr>
            <td><span id="dbname"      class="ief %(dbname_css)s" data-ief-edittype="select" data-ief-selectvalues='[{"value":"UNP","label":"UNP","selected":false},{"value":"GB","label":"GB","selected":false}]'>%(dbname)s</span></td>
            <td><span id="dbaccession" class="ief %(dbaccession_css)s">%(dbaccession)s</span></td>
            <td><span id="dbseqbegin"  class="ief %(dbseqbegin_css)s">%(dbseqbegin)s</span></td>
            <td><span id="dbseqend"    class="ief %(dbseqend_css)s">%(dbseqend)s</span></td>
            </tr>
            </table>
            <input type="submit" name="submit" value="Submit" class="disableonclick submitparentform" />
            <!-- <input type="reset" name="reset" value="Reset" /> -->
        </form>
        </div>
        """
        #
        dbName, _dbCode, dbAccession, dbIsoform, dbSeqBegin, dbSeqEnd = self.__getCurrentRefDetails(entityId, partId=int(partId))
        selectedAuthId = ""
        for sId in self.__selectIdList:
            if sId.startswith("auth"):
                tL = sId.split("_")
                if (len(tL) > 3) and tL[1] == entityId:
                    selectedAuthId = sId
                    break
                #
            #
        #
        pD = {}
        pD["entryid"] = entryId
        pD["sessionid"] = self.__sessionId
        pD["selected_auth_id"] = selectedAuthId
        pD["partid"] = partId
        pD["entityid"] = entityId
        #
        pD["dbname"] = self.__srd.convertDbNameToResource(dbName)
        #
        if dbIsoform is not None and len(dbIsoform) > 0:
            pD["dbaccession"] = dbIsoform
        else:
            pD["dbaccession"] = dbAccession
        #
        pD["dbseqbegin"] = dbSeqBegin
        pD["dbseqend"] = dbSeqEnd
        #
        for item in ("dbname", "dbaccession", "dbseqbegin", "dbseqend"):
            kyCss = item + "_css"
            if pD[item] is None or len(str(pD[item])) < 1:
                pD[item] = self.__placeHolderValue
                pD[kyCss] = "greyedout"
            else:
                pD[kyCss] = ""
            #
        #
        rD = {}
        rD["htmlcontent"] = form_template % pD
        return rD

    def seqDbRefFormResponder(self):
        """Update the sequence data store using data sequence database reference
        Form data encoded in the input request object --
        Return the packed sequence label of the added reference sequence or None
        """
        try:
            selectedAuthId = self._reqObj.getValue("selected_auth_id")
            partId = int(str(self._reqObj.getValue("partid")))
            entityId = self._reqObj.getValue("entityid")
            dbName = self._reqObj.getValue("dbname")
            dbAccession = self._reqObj.getValue("dbaccession")
            if (dbName == self.__placeHolderValue) or (dbAccession == self.__placeHolderValue):
                # This will cause a crash - until we can figure out what is supposed to happen
                return packedSeqLabel  # noqa: F821 pylint: disable=undefined-variable
            #
            dbSeqBegin = self._reqObj.getValue("dbseqbegin")
            dbSeqEnd = self._reqObj.getValue("dbseqend")
            if (dbSeqBegin == self.__placeHolderValue) or (dbSeqEnd == self.__placeHolderValue):
                dbSeqBegin = None
                dbSeqEnd = None
            else:
                try:
                    dbSeqBegin = int(str(dbSeqBegin))
                except:  # noqa: E722 pylint: disable=bare-except
                    dbSeqBegin = None
                    dbSeqEnd = None
                #
                try:
                    dbSeqEnd = int(str(dbSeqEnd))
                except:  # noqa: E722 pylint: disable=bare-except
                    dbSeqBegin = None
                    dbSeqEnd = None
                #
            #
            dbIsoform = ""
            if dbName in ["UNP", "SP", "TR"]:
                tL = dbAccession.split("-")
                if len(tL) > 1:
                    dbIsoform = dbAccession
                    dbAccession = tL[0]
                #
            #
            selectedAuthId, authFeatureDict = self.__getAuthFeatures(selectedAuthId, entityId, partId)
            if not authFeatureDict:
                return "Can not find auth sequence information", ""
            #
            polyTypeCode = "AA"
            if "POLYMER_TYPE" in authFeatureDict:
                polyTypeCode = authFeatureDict["POLYMER_TYPE"]
            #
            fetchUtil = FetchReferenceSequenceUtils(siteId=self.__siteId, seqReferenceData=self.__srd, verbose=self._verbose, log=self._lfh)
            errMsg, refFeatureDict, refSeqList = fetchUtil.fetchReferenceSequence(
                dbName, dbAccession, dbIsoform, polyTypeCode=polyTypeCode, refSeqBeg=dbSeqBegin, refSeqEnd=dbSeqEnd
            )
            #
            if errMsg:
                return errMsg, ""
            #
            alignTool = AlignmentTools(reqObj=self._reqObj, entityId=entityId, seqDataStore=self.getSequenceDataStoreObj(), verbose=self._verbose, log=self._lfh)
            #
            if dbSeqBegin is None:
                dbSeqBegin, dbSeqEnd = alignTool.getAlignRefSeqRange(authSeqId=selectedAuthId, partId=partId, refSeqs=refSeqList)
                errMsg, refFeatureDict, refSeqList = fetchUtil.fetchReferenceSequence(
                    dbName, dbAccession, dbIsoform, polyTypeCode=polyTypeCode, refSeqBeg=dbSeqBegin, refSeqEnd=dbSeqEnd
                )
                #
                if errMsg:
                    return errMsg, ""
                #
            #
            nextAltId = self.getNextAlternativeNumber("ref", entityId, partId)
            self.setSequence(refSeqList, "ref", entityId, seqPartId=partId, seqAltId=nextAltId)
            self.setFeature(
                self.getRefFeatureObj(
                    polyTypeCode, int(partId), authFeatureDict["AUTH_SEQ_NUM_BEGIN"], authFeatureDict["AUTH_SEQ_NUM_END"], authFeatureDict["AUTH_SEQ_PART_TYPE"], refFeatureDict
                ).get(),
                "ref",
                entityId,
                seqPartId=partId,
                seqAltId=nextAltId,
            )
            #
            sL = SequenceLabel()
            sL.set(seqType="ref", seqInstId=entityId, seqPartId=partId, seqAltId=nextAltId, seqVersion=1)
            nextRefLabel = sL.pack()
            self._reqObj.setNewRefId(nextRefLabel)
            #
            alignTool.addRefAlignIndices(authSeqId=selectedAuthId, refSeqId=nextRefLabel, partId=partId)
            #
            self.saveSequenceDataStore()
            #
            return "", nextRefLabel
        except:  # noqa: E722 pylint: disable=bare-except
            self._lfh.write("+UpdatePolymerEntityReference.seqDbRefResponder() failing\n")
            traceback.print_exc(file=self._lfh)
            return "Fetch reference sequence [ dbName=" + dbName + ", Accession=" + dbAccession + "] failed.", ""
        #

    def __getCurrentRefDetails(self, entityId, partId=1):
        self._lfh.write("+UpdatePolymerEntityReference.__getCurrentRefDetails() entityId %r partId %r\n" % (entityId, partId))
        self._lfh.write("+UpdatePolymerEntityReference.__getCurrentRefDetails() selectIdList %r\n" % (self.__selectIdList))
        sL = SequenceLabel()
        for sId in self.__selectIdList:
            if sId.startswith("ref"):
                sL.unpack(sId)
                _seqType, seqInstId, seqPartId, seqAltId, seqVersion = sL.get()
                self._lfh.write("+UpdatePolymerEntityReference.__getCurrentRefDetails() testing seqInstId %r seqPartId %r\n" % (seqInstId, seqPartId))
                if entityId == seqInstId and partId == seqPartId:
                    fD = self.getFeature("ref", seqInstId, seqPartId, seqAltId, seqVersion)
                    return fD["DB_NAME"], fD["DB_CODE"], fD["DB_ACCESSION"], fD["DB_ISOFORM"], fD["REF_MATCH_BEGIN"], fD["REF_MATCH_END"]
                #
            #
        #
        return None, None, None, None, None, None

    def __getAuthFeatures(self, authId, entityId, partId):
        """Get feature data for the author sequence entity sequence.
        Returns a dictionary of features for the more recent sequence version.
        """
        sL = SequenceLabel()
        if sL.unpack(authId):
            seqType, seqInstId, seqPartId, seqAltId, seqVersion = sL.get()
            featureD = self.getFeature(seqType, seqInstId, seqPartId, seqAltId, seqVersion)
            if featureD:
                return authId, featureD
            #
        #
        verList = self.getVersionIdList("auth", entityId, partId, 1)
        if len(verList) == 0:
            return authId, {}
        #
        sL.set(seqType="auth", seqInstId=entityId, seqPartId=partId, seqAltId=1, seqVersion=verList[0])
        selectedAuthId = sL.pack()
        return selectedAuthId, self.getFeature("auth", entityId, partId, 1, verList[0])
