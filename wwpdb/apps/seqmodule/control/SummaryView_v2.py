##
# File:    SummaryView.py
# Date:    19-Jan-2010
# Updates:
# 09-Feb-2010  jdw Adopt SequenceFeature() to access to feature dictionaries.
# 12-Feb-2010  jdw Move example data loader to the SequenceDataImportExample class.
# 14-Feb-2010  jdw Add statistics updater after primary sequence data load.
# 20-Apr-2010  jdw Ported to module seqmodule
# 27-Apr-2010  jdw Add row-level data dictionary to the summary data object.
# 02-May-2010  jdw Add SequenceSelection()
# 09-Aug-2010  rps __loadSummary() --> Highest numbered version of coordinate sequence to be selected by default in "SELECT" column of UI
# 08-Jan-2013  jdw Select higher version of author sequence in reload operations
# 03-Mar-2013  jdw Refactor -
# 05-Mar-2013  jdw Move default sequence selection from SequenceSelection()
#                  Add default self-reference selections
#                  Implement stored sorting order for reference sequeneces --
#                  Limit reference sequence to maxRefAlign --
# 12-Mar-2013  jdw adjust the assembly of author sequence details -
# 13-Apr-2013  jdw change self reference identifiers
# 15-Nov-2013  jdw major overhaul to support expanded author content display and entity review
# 12-Dec-2013  jdw make default alignment selections.
# 19-Dec-2013  jdw select any added or edited reference sequence by default
# 15-May-2014  jdw add status field for instance id with xyz sequence groups
# 22-May-2014  jdw add method to provide entry details --
# 30-Aug-2017  zf  change self.__reqObj.getSummaryAlignList & self.__reqObj.getSummarySelectList to use latest UI input values
# 29-Jun-2021  zf  added self.__checkPolyAlaAssignment()
# 19-Oct-2022  zf  added section for generation summary landing page
##
"""
Controlling class for the production of data for the summary sequence view.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.08"

try:
    import cPickle as pickle
except ImportError:
    import pickle as pickle
#

import os
import sys
import traceback

from wwpdb.apps.seqmodule.io.AlignmentDataStore import AlignmentDataStore
from wwpdb.apps.seqmodule.io.FetchSeqInfoUtils import FetchSeqInfoUtils
from wwpdb.apps.seqmodule.io.SequenceDataStore import SequenceDataStore
from wwpdb.apps.seqmodule.util.SequenceLabel import SequenceLabel, SequenceFeature
from wwpdb.apps.seqmodule.util.SequenceFeatureDepict import SequenceFeatureDepict

from wwpdb.apps.seqmodule.util.SequenceAssign import SequenceAssignDepositor
from wwpdb.apps.seqmodule.update.UpdatePolymerEntitySourceDetails import UpdatePolymerEntitySourceDetails

from wwpdb.io.misc.FormatOut import FormatOut


class SummaryView(object):
    """Controlling class for the production of data for the summary sequence view.

    Supported operations:

     load             = load summary from current sequence data store using default sequence selection.
     reload           = reload summary using current sequence data and current user selections.

    """

    def __init__(self, reqObj=None, maxRefAlign=100, verbose=False, log=sys.stderr):
        self.__verbose = verbose
        self.__debug = False
        self.__reqObj = reqObj
        self.__lfh = log
        self.__maxRefAlign = maxRefAlign
        #
        # placeholders for sequence identifiers picked on the summary page as selected and/or to be aligned.
        self.__summarySeqAlignList = []
        self.__summarySeqSelectList = ""
        #
        self.__natureSourceTaxIds = {}
        self.__partRangeErrorMsg = ""
        self.__sds = None
        #
        self.__summaryPageInfoMap = {}
        #
        self.__setup()

    def __setup(self):
        try:
            self.__sessionObj = self.__reqObj.getSessionObj()
            self.__sessionPath = self.__sessionObj.getPath()  # pylint: disable=unused-private-member
            self.__sds = SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            self.__sdu = UpdatePolymerEntitySourceDetails(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            #
            # Selections coming from the the web request --
            #
            # self.__summarySeqAlignList = self.__reqObj.getSummaryAlignList(usingRevAllAlignIds=False)
            # self.__summarySeqSelectList = self.__reqObj.getSummarySelectList(usingRevSlectIds=False)
            self.__summarySeqAlignList = str(self.__reqObj.getValue("allalignids")).split(",")
            self.__summarySeqSelectList = str(self.__reqObj.getValue("selectids")).split(",")
            self.__updateSelectList()
            if self.__verbose:
                self.__lfh.write("+SummaryView.__setup() request input align list  %r\n" % (self.__summarySeqAlignList))
                self.__lfh.write("+SummaryView._setup()  request input select list %r\n" % (self.__summarySeqSelectList))
            #
            # Reset any newly loaded reference sequence id as it is included above --
            self.__reqObj.setNewRefId("")
        except Exception as _e:  # noqa: F841
            if self.__verbose:
                self.__lfh.write("+SummaryView.__setup() sessionId %s failed\n" % (self.__sessionObj.getId()))
                traceback.print_exc(file=self.__lfh)
            #
        #

    def __updateSelectList(self):
        """ """
        updatedids = str(self.__reqObj.getValue("updatedids"))
        if not updatedids:
            return
        #
        updatedList = updatedids.split(",")
        #
        entityId = str(self.__reqObj.getValue("activegroupid"))
        chainIdList = self.__sds.getGroup(entityId)
        for seqId in self.__summarySeqSelectList:
            tL = str(seqId).strip().split("_")
            if len(tL) < 3:
                continue
            #
            elif (tL[0] in ("selfref", "auth", "ref")) and (tL[1] == entityId):
                continue
            elif (tL[0] == "xyz") and (tL[1] in chainIdList):
                continue
            #
            updatedList.append(seqId)
        #
        self.__summarySeqSelectList = updatedList

    def __finish(self):
        try:
            # Persist the current summary selection list --
            self.__reqObj.setSummarySelectList(self.__summarySeqSelectList)
            self.__reqObj.setSummaryAlignList(self.__summarySeqAlignList)
            self.__sds.setSelectedIds(idList=self.__summarySeqSelectList)
            self.__sds.serialize()
        except Exception as _e:  # noqa: F841
            if self.__verbose:
                self.__lfh.write("+SummaryView.__finish() sessionId %s failed\n" % (self.__sessionObj.getId()))
                traceback.print_exc(file=self.__lfh)

    def getEntryDetails(self, kyList=None):
        if kyList is None:
            kyList = ["STRUCT_TITLE", "CITATION_TITLE", "PDB_ID"]

        eD = {}
        for ky in kyList:
            eD[ky] = self.__sds.getEntryDetail(detailKey=ky)
        #
        return eD

    def loadSummary(self, operation="load"):
        """Create the data structure to populate the HTML pages containing alignment summary --

        Options for operation:
        + load   - Performs a default selection of the best matching reference sequences then
                   creates the summary data structure.
        + reload - reloads current data from the sequence and alignment data stores.
        """
        if self.__verbose:
            self.__lfh.write("+SummaryView.loadSummary() operation %s : sessionId %s\n" % (operation, self.__sessionObj.getId()))
            self.__lfh.write("+SummaryView.loadSummary() request input align list  %r\n" % (self.__summarySeqAlignList))
            self.__lfh.write("+SummaryView.loadSummary() request input select list %r\n" % (self.__summarySeqSelectList))
        if operation == "load":
            # Using coming from persisted data only here
            #  -- ignore user selections from web context as these may not exist or may need to be overridden --
            self.__summarySeqSelectList = self.__sds.getSelectedIds()
        #
        self.__summarySeqAlignList = self.__summarySeqSelectList

        self.__sdu.updateAuthEntityDetails(selectIdList=self.__summarySeqSelectList)
        self.__sds = SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        return self.__loadSummary(op=operation)

    def getGroupIdList(self):
        return self.__sds.getGroupIds()

    def getSummaryPageObj(self):
        return self.__summaryPageInfoMap

    def __loadSummary(self, op="load"):
        """Assemble the data for sequence summary view using current contents of the session sequence data store.

        Returns:  summaryDataObj [<entity/groupId>]  ['auth'|'xyz'] ->  {data dict}
                                                     ['ref']        -> [(partId,{data dict}),(partId,{data dict}),...]
                               {data dict}  with keys

                               dT['ROW_IDS']       =rowIdList
                               dT['ROW_STATUS']    =rowStatusList
                               dT['ROW_DATA_DICT'] =rowDataDictList
                               dT['SELF_REFERENCE_FLAG']=Boolean
        """
        #
        if self.__verbose:
            self.__lfh.write("+SummaryView.__loadSummary() with sessionId %s\n" % self.__sessionObj.getId())
        #
        summaryDataObj = {}
        #
        gIdList = self.__sds.getGroupIds()
        #
        if self.__verbose:
            self.__lfh.write("+SummaryView.__loadSummary() group list is %r\n" % gIdList)
            self.__lfh.write("+SummaryView.__loadSummary() starting selection list %r\n" % self.__summarySeqSelectList)
            self.__lfh.write("+SummaryView.__loadSummary() starting alignment list %r\n" % self.__summarySeqAlignList)
        #
        polyAlaCaseList = []
        warningMsg = ""
        for gId in gIdList:
            #
            summaryDataObj[gId] = {}
            #
            dT, polyALA_assignment, authSummaryPageD, method, polyType, dbAccessionList, seqLength = self.__buildAuthSection(groupId=gId, op=op)
            if polyALA_assignment > 0:
                polyAlaCaseList.append((gId, polyALA_assignment))
            #
            summaryDataObj[gId]["auth"] = dT
            #
            rL, withGapScoreList, withoutGapScoreList, hitDbInfoList, hitDbIdList, authTaxIdList = self.__buildReferenceSection(groupId=gId)
            summaryDataObj[gId]["ref"] = rL
            #
            summaryDataObj[gId]["xyz"], missingResidueInfoList = self.__buildCoordinateSection(groupId=gId, authSeqLength=seqLength)
            if len(missingResidueInfoList) > 0:
                for missingInfo in missingResidueInfoList:
                    warningMsg += missingInfo + "<br />\n"
                #
            #
            decoration_start = ""
            decoration_end = ""
            if len(self.__sds.getGroup(gId)) == 0:
                decoration_start = '<p style="text-decoration: line-through;">'
                decoration_end = "</p>"
            #
            self.__summaryPageInfoMap[gId] = {}
            if decoration_start:
                self.__summaryPageInfoMap[gId]["entity_id"] = decoration_start + '<span class="width20px"> &nbsp;' + gId + '&nbsp; </span>' + decoration_end
            else:
                # self.__summaryPageInfoMap[gId]["entity_id"] = '<a id="page_' + gId + '" href="#" class="page_control"><span class="width20px"> &nbsp;' \
                self.__summaryPageInfoMap[gId]["entity_id"] = '<a id="page_' + gId + '" href="#closecompleted" class="page_control"><span class="width20px">' \
                    + ' &nbsp;' + gId + '&nbsp; </span></a>'
            #
            self.__summaryPageInfoMap[gId]["chain_ids"] = ",".join(self.__sds.getGroup(gId))
            for tupL in (("mol_names", "ENTITY_DESCRIPTION"), ("source_names", "SOURCE_ORGANISM"), ("tax_ids", "SOURCE_TAXID")):
                self.__summaryPageInfoMap[gId][tupL[0]] = decoration_start + self.__buildNameSourceSummaryPage(tupL[1], authSummaryPageD) + decoration_end
            #
            if (len(withGapScoreList) > 0) and (len(withoutGapScoreList) > 0):
                self.__summaryPageInfoMap[gId]["identity_scores"] = decoration_start + '<span class="detailkey">w/ gaps: </span><span class="detailvalue">' + \
                    ",".join(withGapScoreList) + '</span><br/><span class="detailkey">w/o gaps: </span><span class="detailvalue">' + \
                    ",".join(withoutGapScoreList) + '</span>' + decoration_end
            #
            if len(hitDbInfoList) > 0:
                self.__summaryPageInfoMap[gId]["ref_db_ids"] = decoration_start + "/".join(hitDbInfoList) + "<br/>"
            else:
                self.__summaryPageInfoMap[gId]["ref_db_ids"] = decoration_start + "-<br/>"
            #
            if len(dbAccessionList) > 0:
                displayList = []
                for valTup in dbAccessionList:
                    accCode = valTup[0]
                    if hitDbIdList and (valTup[0] not in hitDbIdList):
                        accCode = '<span style="color:red">' + valTup[0] + '</span>'
                    #
                    taxId = ""
                    if len(valTup[1]) > 0:
                        taxId = "[" + valTup[1] + "]"
                        if valTup[1] not in authTaxIdList:
                            taxId = '[<span style="color:red">' + valTup[1] + '</span>]'
                        #
                    #
                    if taxId:
                        accCode += " " + taxId
                    #
                    displayList.append(accCode)
                #
                self.__summaryPageInfoMap[gId]["ref_db_ids"] += "(" + "/".join(displayList) + ")" + decoration_end
            else:
                self.__summaryPageInfoMap[gId]["ref_db_ids"] += "(-)" + decoration_end
            #
            displayMethod = method.upper()
            warningErrorMsgs = ""
            alignData = AlignmentDataStore(reqObj=self.__reqObj, entityId=gId, verbose=self.__verbose, log=self.__lfh)
            warningErrorD = alignData.getSummaryPageInfo()
            if warningErrorD:
                for item in ("acetylation", "amidation", "cloning artifact", "conflict", "chromophore", "deletion", "engineered mutation", "expression tag",
                             "initiating methionine", "insertion", "linker", "microheterogeneity", "microheterogeneity/modified residue", "mismatch",
                             "modified residue", "variant"):
                    if item not in warningErrorD:
                        continue
                    #
                    redFlag = False
                    if item == "mismatch":
                        redFlag = True
                    elif item == "expression tag":
                        if displayMethod == "NAT":
                            displayMethod = '<span style="color:red">' + method.upper() + '</span>'
                            redFlag = True
                        elif warningErrorD[item] > 19:
                            redFlag = True
                        #
                    #
                    if redFlag:
                        val = '<span style="color:red">' + item.capitalize() + '(' + str(warningErrorD[item]) + ')</span>'
                    else:
                        val = item.capitalize() + "(" + str(warningErrorD[item]) + ")"
                    #
                    if warningErrorMsgs != "":
                        warningErrorMsgs += "<br/>"
                    #
                    warningErrorMsgs += val
                #
            #
            if not warningErrorMsgs:
                warningErrorMsgs = "-"
            #
            self.__summaryPageInfoMap[gId]["entity_types"] = decoration_start + displayMethod + "<br/>" + polyType + decoration_end
            self.__summaryPageInfoMap[gId]["warn_err_msgs"] = decoration_start + warningErrorMsgs + decoration_end
        #
        self.__finish()
        #
        if len(self.__natureSourceTaxIds) > 1:
            if warningMsg:
                warningMsg += "<br />\n"
            #
            warningMsg += "Entry contains multiple natural sources:<br />\n"
            for k, v in self.__natureSourceTaxIds.items():
                if len(v) > 1:
                    warningMsg += "Entities '" + "', '".join(v) + "' have "
                else:
                    warningMsg += "Entity '" + "', '".join(v) + "' has "
                #
                warningMsg += "source taxonomy Id '" + k + "'.<br />\n"
            #
        #
        if self.__partRangeErrorMsg != "":
            if warningMsg:
                warningMsg += "<br />\n"
            #
            warningMsg += "Part range errors:<br />\n" + self.__partRangeErrorMsg
        #
        if polyAlaCaseList:
            for polyAlaCaseTupL in polyAlaCaseList:
                if warningMsg:
                    warningMsg += "<br />\n"
                #
                warningMsg += "Entity [" + polyAlaCaseTupL[0] + "]"
                if polyAlaCaseTupL[1] == 1:
                    warningMsg += " is composed only of poly-ALA.<br />\n"
                else:
                    warningMsg += " has stretches of poly-ALA.<br />\n"
                #
            #
        #
        return summaryDataObj, warningMsg

    def __buildAuthSection(self, groupId=None, op="load"):
        """Assemble the data content for the author sequence summary view.

        Returns: summaryDataObject)

        ** Always show a single full sequence for each entity and store multi-part information as details

        """
        #
        # Each author entity sequence is shown once and is labeled by its first instance in the group.
        #
        # seqId0 = groupId
        #
        partIdList = self.__sds.getPartIds(groupId, dataType="sequence", seqType="auth")
        if len(partIdList) == 0:
            return {}, 0, {}, "", "", [], 0
        #
        altId = 1
        verList = self.__sds.getVersionIds(groupId, partId=partIdList[0], altId=altId, dataType="sequence", seqType="auth")
        if len(verList) == 0:
            return {}, 0, {}, "", "", [], 0
        #
        if self.__verbose:
            self.__lfh.write("SummaryView.__buildAuthSection() groupId %r op %s \n" % (groupId, op))
        #
        detailsD, dbAccessionList = self.__getAuthFeaturesAll(groupId, groupId, partIdList)
        #
        dT = self.__getAuthSection(groupId, partIdList[0], altId, verList, detailsD, len(partIdList))
        #
        # Get polyALA_assignment
        #
        polyALA_assignment, method, polyType, seqLength = self.__checkPolyAlaAssignment(seqId=groupId, partId=partIdList[0], altId=altId, version=verList[0])
        #
        summaryPageD = self.__getSummaryPageInfo(seqId=groupId, partIdList=partIdList, altId=altId)
        #
        return dT, polyALA_assignment, summaryPageD, method, polyType, dbAccessionList, seqLength

    def __buildCoordinateSection(self, groupId=None, authSeqLength=0):  # pylint: disable=unused-argument
        """Assemble the data content for the coordinate sequence summary view.

        Returns: summaryDataObject

        """
        seqLabel = SequenceLabel()
        seqFeature = SequenceFeature()
        #
        # For each coordinate sequence -
        #
        seqIdList = self.__sds.getGroup(groupId)
        #
        rowIdList = []
        rowStatusList = []
        rowDataDictList = []
        missingResidueInfoList = []
        #
        partId = 1
        altId = 1
        for seqId in seqIdList:
            verList = self.__sds.getVersionIds(seqId=seqId, partId=partId, altId=altId, dataType="sequence", seqType="xyz")
            if len(verList) < 1:
                continue
            #
            maxVrsnNum = verList[0]
            seqLabel.set(seqType="xyz", seqInstId=seqId, seqPartId=partId, seqAltId=altId, seqVersion=maxVrsnNum)
            maxSeqXyzId = seqLabel.pack()
            #
            # select the highest version sequence
            if not (maxSeqXyzId in self.__summarySeqSelectList):
                self.__summarySeqSelectList.append(maxSeqXyzId)
                self.__summarySeqAlignList.append(maxSeqXyzId)
            #
            for ver in verList:
                seqXyz = self.__sds.getSequence(seqId, "xyz", partId=partId, altId=altId, version=ver)
                seqXyzFD = self.__sds.getFeature(seqId, "xyz", partId=partId, altId=altId, version=ver)
                seqLabel.set(seqType="xyz", seqInstId=seqId, seqPartId=partId, seqAltId=altId, seqVersion=ver)
                seqXyzId = seqLabel.pack()
                #
                seqFeature.set(seqXyzFD)
                sfd = SequenceFeatureDepict(sfObj=seqFeature, verbose=self.__verbose, log=self.__lfh)
                detailString = sfd.markupXyzFeatures()
                #
                rowLabel = seqId
                rowIdList.append(seqXyzId)
                #
                # unselect any selected lower version for this sequence -
                if (ver < maxVrsnNum) and (seqXyzId in self.__summarySeqSelectList):
                    self.__summarySeqSelectList.remove(seqXyzId)
                    if seqXyzId in self.__summarySeqAlignList:
                        self.__summarySeqAlignList.remove(seqXyzId)
                #
                isSelected = seqXyzId in self.__summarySeqSelectList
                isAligned = seqXyzId in self.__summarySeqAlignList
                # JDW add instance here to provide a selection group
                if ver == maxVrsnNum:
                    rowStatusList.append((isSelected, isAligned, seqId, "maxver"))
#                   if (authSeqLength > 0) and ((authSeqLength * 3) >= (10 * len(seqXyz))):
#                       percent = float((authSeqLength - len(seqXyz)) * 100) / float(authSeqLength)
#                       missingResidueInfoList.append("%.1f" % percent + "% residues of chain '" + seqId + "' (%d/%d residues) are missing in coordinates." %
#                                                     (authSeqLength - len(seqXyz), authSeqLength))
#                   #
                else:
                    rowStatusList.append((isSelected, isAligned, seqId, ""))
                #
                rowDataDict = {}
                rowDataDict["ROW_ID_CODE"] = rowLabel
                rowDataDict["ROW_VERSION"] = ver
                rowDataDict["ROW_SEQ_LENGTH"] = len(seqXyz)
                rowDataDict["ROW_DETAIL_STRING"] = detailString
                rowDataDict.update(seqXyzFD)
                rowDataDictList.append(rowDataDict)
            #
        #
        dT = {}
        dT["ROW_IDS"] = rowIdList
        dT["ROW_STATUS"] = rowStatusList
        dT["ROW_DATA_DICT"] = rowDataDictList
        return dT, missingResidueInfoList

    def __buildReferenceSection(self, groupId=None):
        """Assemble the data content for the reference sequence summary view.

        Returns: [(partId,summaryDataObject),]

               where:
                  partId is the integer identifier for pieces of a multipart entity sequence
                  summaryDataObject contains the 'ref' sequence summary data for the input groupId.
        """
        rL = []
        withGapScoreList = []
        withoutGapScoreList = []
        hitDbInfoList = []
        hitDbIdList = []
        authTaxIdList = []
        seqLabel = SequenceLabel()
        seqFeature = SequenceFeature()
        #
        # use the leading sequence Id in the group as the id for the reference -
        #
        seqIdRef = groupId
        #
        partInfoList = []
        partIdList = self.__sds.getPartIds(seqIdRef, dataType="sequence", seqType="auth")
        for partId in partIdList:
            #
            authVerList = self.__sds.getVersionIds(seqId=seqIdRef, partId=partId, altId=1, dataType="feature", seqType="auth")
            # Changed to original tax ID ( DAOTHER-6126 )
            OrigTaxId = ""
            skipPart = False
            if len(authVerList) > 0:
                authFObj = self.__sds.getFeatureObj(seqIdRef, "auth", partId=partId, altId=1, version=authVerList[0])
                sourceInfo = authFObj.getEntitySourceMethod().upper()
                OrigTaxId = authFObj.getSourceTaxIdOrig()
                if OrigTaxId and (OrigTaxId not in authTaxIdList):
                    authTaxIdList.append(OrigTaxId)
                #
                TaxId = authFObj.getSourceTaxId()
                if (sourceInfo == "NAT") and TaxId:
                    if TaxId in self.__natureSourceTaxIds:
                        self.__natureSourceTaxIds[TaxId].append(groupId)
                    else:
                        self.__natureSourceTaxIds[TaxId] = [groupId]
                    #
                #
                authPartId, authSeqBeg, authSeqEnd, authSeqPartType = authFObj.getAuthPartDetails()
                partInfoList.append((authPartId, authSeqBeg, authSeqEnd))
                skipPart = len(authSeqPartType) > 1 and authSeqPartType.lower() != "biological sequence"
            #
            if skipPart:
                self.__lfh.write("+SummaryView.__buildReferenceSection() entity %s skipping part %s type %s\n" % (groupId, partId, authSeqPartType))
                continue
            #
            rowIdList = []
            rowStatusList = []
            rowDataDictList = []
            # List of reference sequences for this group only for the leading sequence -
            #
            altIdList = self.__sds.getAlternativeIds(seqIdRef, dataType="sequence", seqType="ref", partId=partId)

            srId = "selfref_" + str(groupId) + "_" + str(partId)
            if srId in self.__summarySeqSelectList:
                if self.__verbose:
                    self.__lfh.write("+SummaryView.__buildReferenceSection() using self-refereence for entity %s part %s\n" % (groupId, partId))
                selfRefFlag = True
            else:
                selfRefFlag = False
            #
            for altId in altIdList[: self.__maxRefAlign]:
                #
                verList = self.__sds.getVersionIds(seqId=seqIdRef, partId=partId, altId=altId, dataType="sequence", seqType="ref")
                for ver in verList:
                    #
                    seqRefFD = self.__sds.getFeature(seqIdRef, "ref", partId=partId, altId=altId, version=ver)
                    seqFeature.set(seqRefFD)
                    seqLabel.set(seqType="ref", seqInstId=seqIdRef, seqPartId=partId, seqAltId=altId, seqVersion=ver)
                    seqRefId = seqLabel.pack()
                    #
                    taxIdWarningFlag = False
                    RefTaxId = seqFeature.getSourceTaxId()
                    if OrigTaxId and RefTaxId and (OrigTaxId != RefTaxId):
                        taxIdWarningFlag = True
                    #
                    isSelected = seqRefId in self.__summarySeqSelectList
                    if isSelected:
                        if ("AUTH_REF_SEQ_SIM_WITH_GAPS" in seqRefFD) and (float(seqRefFD["AUTH_REF_SEQ_SIM_WITH_GAPS"]) > 0.001):
                            withGapScoreList.append("%6.3f" % float(seqRefFD["AUTH_REF_SEQ_SIM_WITH_GAPS"]))
                        #
                        if ("AUTH_REF_SEQ_SIM" in seqRefFD) and (float(seqRefFD["AUTH_REF_SEQ_SIM"]) > 0.001):
                            withoutGapScoreList.append("%6.3f" % float(seqRefFD["AUTH_REF_SEQ_SIM"]))
                        #
                        dbName = ""  # noqa: F841 pylint: disable=unused-variable
                        if ("DB_NAME" in seqRefFD) and seqRefFD["DB_NAME"]:
                            dbName = seqRefFD["DB_NAME"]  # noqa: F841
                        #
                        dbAccession = ""
                        if ("DB_ACCESSION" in seqRefFD) and seqRefFD["DB_ACCESSION"]:
                            dbAccession = seqRefFD["DB_ACCESSION"]
                        #
                        if ("DB_ISOFORM" in seqRefFD) and seqRefFD["DB_ISOFORM"]:
                            dbAccession = seqRefFD["DB_ISOFORM"]
                        #
                        if dbAccession:
                            refTaxId = RefTaxId
                            if taxIdWarningFlag:
                                refTaxId = '<span style="color:red">' + RefTaxId + "</span>"
                            #
                            if refTaxId:
                                hitDbInfoList.append(dbAccession + " [" + refTaxId + "]")
                            else:
                                hitDbInfoList.append(dbAccession)
                            #
                            if dbAccession not in hitDbIdList:
                                hitDbIdList.append(dbAccession)
                            #
                        #
                    #
                    isAligned = seqRefId in self.__summarySeqAlignList
                    rowStatusList.append((isSelected, isAligned))
                    #
                    sfd = SequenceFeatureDepict(sfObj=seqFeature, verbose=self.__verbose, log=self.__lfh)
                    detailString = sfd.markupReferenceAlignmentFeatures()
                    featureString = sfd.markupReferenceFeatures()
                    rowLabel = sfd.markupDatabaseReferenceWithUrl(altId)
                    rowIdList.append(seqRefId)

                    lengthRefSeq = seqFeature.getMatchLength()
                    # authSimWithGaps = seqFeature.getAuthRefSimWithGaps()
                    #
                    #
                    rowDataDict = {}
                    rowDataDict["ROW_ID_CODE"] = rowLabel
                    rowDataDict["ROW_VERSION"] = ver
                    rowDataDict["ROW_SEQ_LENGTH"] = lengthRefSeq
                    # rowDataDict['ROW_AUTH_REF_SIM'] = "%.3f" % authSimWithGaps
                    rowDataDict["ROW_AUTH_REF_SIM"] = sfd.markupReferenceSimilarttFeatures()
                    rowDataDict["ROW_DETAIL_STRING"] = detailString
                    rowDataDict["ROW_FEATURE_STRING"] = featureString
                    rowDataDict["ROW_IS_SELECTED"] = isSelected
                    rowDataDict["ROW_IS_ALIGNED"] = isAligned
                    # rowDataDict.update(seqRefFD)
                    for key, _val in seqRefFD.items():
                        if (key == "SOURCE_TAXID") and taxIdWarningFlag:
                            rowDataDict[key] = '<span style="color:red">' + seqRefFD[key] + "</span>"
                        else:
                            rowDataDict[key] = seqRefFD[key]
                        #
                    #
                    rowDataDictList.append(rowDataDict)
            #
            dT = {}
            dT["ROW_IDS"] = rowIdList
            dT["ROW_STATUS"] = rowStatusList
            dT["ROW_DATA_DICT"] = rowDataDictList
            dT["SELF_REFERENCE_FLAG"] = selfRefFlag
            rL.append((partId, dT))
        #
        self.__checkPartRange(groupId, partInfoList)
        #
        return rL, withGapScoreList, withoutGapScoreList, hitDbInfoList, hitDbIdList, authTaxIdList

    def __buildNameSourceSummaryPage(self, key, summaryInfoD):
        """
        """
        currList = []
        origList = []
        if (key in summaryInfoD) and summaryInfoD[key]:
            for val in summaryInfoD[key]:
                if val:
                    currList.append(val)
                else:
                    currList.append("-")
                #
            #
        else:
            currList.append("-")
        #
        origKey = key + "_ORIG"
        if (origKey in summaryInfoD) and summaryInfoD[origKey]:
            for val in summaryInfoD[origKey]:
                if val:
                    origList.append(val)
                else:
                    origList.append("-")
                #
            #
        else:
            origList.append("-")
        #
        chimera = ""
        if (key == "ENTITY_DESCRIPTION") and (len(currList) > 1):
            chimera = " chimera"
        #
        return "/".join(currList) + chimera + "<br/>(" + "/".join(origList) + chimera + ")"

    def __getAuthFeaturesAll(self, entityId, seqId0, partIdList):
        """Consolidate the feature data for the author sequences entity sequence respecting
        potential mutiple source details.

        Original and current features are assembled.

        Returns a dictionaries of features for each sequence version.

        Details[ver][partId]['current']={'partdetails','sourceAndStrain','description','hostorg'}
        Details[ver][partId]['author']={'partdetails','sourceAndStrain','description','hostorg'}
        """
        if self.__verbose:
            self.__lfh.write("+SummaryView.__getAuthFeaturesAll() entityId %r seqId0 %r partIdList %r\n" % (entityId, seqId0, partIdList))

        #
        authRefAssignText, dbAccessionList = self.__markupAuthorAssignments(entityId)
        #
        #
        altId = 1
        detailsD = {}
        for partId in partIdList:
            verList = self.__sds.getVersionIds(seqId0, partId=partId, altId=altId, dataType="sequence", seqType="auth")
            for ver in verList:
                if ver not in detailsD:
                    detailsD[ver] = {}

                sfObj = self.__sds.getFeatureObj(seqId=seqId0, seqType="auth", partId=partId, altId=altId, version=ver)
                sfd = SequenceFeatureDepict(sfObj=sfObj, verbose=self.__verbose, log=self.__lfh)

                if partId not in detailsD[ver]:
                    detailsD[ver][partId] = {}

                detailsD[ver][partId]["current"] = sfd.markupCurrentEntityDetails()
                detailsD[ver][partId]["author"] = sfd.markupAuthorEntityDetails()

                if len(authRefAssignText) > 0:
                    detailsD[ver][partId]["author"]["description"] += authRefAssignText

        for ver in detailsD.keys():
            np = 0
            for partId in partIdList:
                if partId in detailsD[ver]:
                    np += 1
            detailsD[ver]["numparts"] = np
        if self.__debug:
            self.__lfh.write("+SummaryView.__getAuthFeaturesAll() detailsD entity %r seqId %r partIdList %r\n" % (entityId, seqId0, partIdList))
            fOut = FormatOut()
            fOut.autoFormat("Summary Entity Details", detailsD, 3, 3)
            fOut.writeStream(self.__lfh)
            fOut.clear()
        #
        return detailsD, dbAccessionList

    def __checkPolyAlaAssignment(self, seqId="", partId="1", altId=1, version="1"):
        """Check if the sequence contains 10 or more consecutive ALA residues
        """
        authFObj = self.__sds.getFeatureObj(seqId, "auth", partId=partId, altId=altId, version=version)
        if authFObj.getPolymerType() != "AA":
            return 0, authFObj.getEntitySourceMethod(), authFObj.getPolymerLinkingType(), 0
        #
        seqAuth = self.__sds.getSequence(seqId=seqId, seqType="auth", partId=partId, altId=altId, version=version)
        has_consecutive_ALA = False
        count = 0
        for seqTupL in seqAuth:
            if seqTupL[0] == "ALA":
                count += 1
            else:
                if count > 9:
                    has_consecutive_ALA = True
                #
                count = 0
            #
        #
        if count > 9:
            has_consecutive_ALA = True
        #
        if count == len(seqAuth):
            return 1, authFObj.getEntitySourceMethod(), authFObj.getPolymerLinkingType(), len(seqAuth)
        elif has_consecutive_ALA:
            return 2, authFObj.getEntitySourceMethod(), authFObj.getPolymerLinkingType(), len(seqAuth)
        #
        return 0, authFObj.getEntitySourceMethod(), authFObj.getPolymerLinkingType(), len(seqAuth)

    def __getSummaryPageInfo(self, seqId="1", partIdList=None, altId=1):
        """
        """
        if partIdList is None:
            partIdList = []

        infoKeys = ("ENTITY_DESCRIPTION", "ENTITY_DESCRIPTION_ORIG", "SOURCE_ORGANISM", "SOURCE_ORGANISM_ORIG", "SOURCE_TAXID", "SOURCE_TAXID_ORIG")
        #
        retValD = {}
        for partId in partIdList:
            verList = self.__sds.getVersionIds(seqId, partId=partId, altId=altId, dataType="sequence", seqType="auth")
            valD = {}
            if len(verList) > 0:
                authFObj = self.__sds.getFeatureObj(seqId, "auth", partId=partId, altId=altId, version=verList[0])
                mapD = authFObj.get()
                for key in infoKeys:
                    if (key in mapD) and mapD[key]:
                        valD[key] = mapD[key]
                        #
                        additionlKey = ""
                        if key == "SOURCE_ORGANISM":
                            additionlKey = "SOURCE_STRAIN"
                        elif key == "SOURCE_ORGANISM_ORIG":
                            additionlKey = "SOURCE_STRAIN_ORIG"
                        #
                        if additionlKey and (additionlKey in mapD) and mapD[additionlKey]:
                            valD[key] = mapD[key] + " " + mapD[additionlKey]
                        #
                    #
                #
            #
            if ("SOURCE_TAXID" in valD) and ("SOURCE_TAXID_ORIG" in valD) and (valD["SOURCE_TAXID"] != valD["SOURCE_TAXID_ORIG"]):
                for key in ("SOURCE_ORGANISM_ORIG", "SOURCE_TAXID_ORIG"):
                    if (key in valD) and valD[key]:
                        valD[key] = '<span style="color:red">' + valD[key] + '</span>'
                    #
                #
            #
            for key in infoKeys:
                val = ""
                if (key in valD) and valD[key]:
                    val = valD[key]
                #
                if key in retValD:
                    retValD[key].append(val)
                else:
                    retValD[key] = [val]
                #
            #
        #
        return retValD

    def __getAuthSection(self, seqId0, partId0, altId, verList, detailsD, partIdListLength):
        """Assemble the data content for the author sequence summary view.

        Returns: summaryDataObject)

        ** Always show a single full sequence for each entity and store multi-part information as details

        """
        seqLabel = SequenceLabel()
        #
        rowIdList = []
        rowStatusList = []
        rowDataDictList = []
        #
        for ver in verList:
            seqAuth = self.__sds.getSequence(seqId=seqId0, seqType="auth", partId=partId0, altId=altId, version=ver)
            seqLabel.set(seqType="auth", seqInstId=seqId0, seqPartId=partId0, seqAltId=altId, seqVersion=ver)
            seqAuthId = seqLabel.pack()
            rowIdList.append(seqAuthId)
            #
            isSelected = seqAuthId in self.__summarySeqSelectList
            isAligned = seqAuthId in self.__summarySeqAlignList
            #
            # if op == 'reload':
            #     isSelected = (ver == max(verList))
            #     isAligned = (ver == max(verList))
            # #
            #
            rowStatusList.append((isSelected, isAligned))
            #
            rowDataDict = {}
            rowDataDict["ROW_VERSION"] = ver
            rowDataDict["ROW_SEQ_LENGTH"] = len(seqAuth)
            rowDataDict["ROW_DETAIL_STRING"] = ""
            rowDataDict.update(detailsD[ver])
            rowDataDictList.append(rowDataDict)
        #
        dT = {}
        dT["ROW_IDS"] = rowIdList
        dT["ROW_STATUS"] = rowStatusList
        dT["ROW_DATA_DICT"] = rowDataDictList
        dT["ENTITY_NUM_PARTS"] = partIdListLength
        #
        return dT

    def __markupAuthorAssignments(self, entityId):
        """Markup the author provided reference assignments -"""
        spStr = "<br />"
        authRefAssignText = ""
        dbAccessionList = []
        dbAccessionMap = {}
        depSeqAssignD = self.__sds.getDepositorReferenceAssignments()
        sADep = SequenceAssignDepositor(verbose=self.__verbose, log=self.__lfh)
        sADep.set(depSeqAssignD)
        refSeqAssignL = sADep.getReferenceList(entityId)

        if self.__debug:
            self.__lfh.write("+SummaryView.__markupAuthorAssignments() for entityId %r author reference count is %d\n" % (entityId, sADep.getReferenceCount(entityId)))
            sADep.printIt(log=self.__lfh)
            for ii, rsa in enumerate(refSeqAssignL):
                self.__lfh.write("+SummaryView.__markupAuthorAssignments() depositor reference  %d\n" % (ii + 1))
                rsa.printIt(self.__lfh)
        #

        for rsa in refSeqAssignL:
            dbName, dbCode, dbAccession = rsa.getDbReference()
            if (len(dbAccession) > 0) and (dbAccession not in [".", "?"]):
                if dbAccession not in dbAccessionList:
                    dbAccessionList.append(dbAccession)
                    if (dbAccession not in dbAccessionMap) and (len(dbName) > 0) and (dbName not in [".", "?"]):
                        dbAccessionMap[dbAccession] = dbName.upper()
                    #
                #
            #
            if ((len(dbAccession) > 0) and (dbAccession not in [".", "?"])) or ((len(dbCode) > 0) and (dbCode not in [".", "?"])):
                tRef = "<b>Depositor reference:&nbsp;</b> %s %s %s" % (dbName, dbAccession, dbCode)
                sp = spStr if len(authRefAssignText) > 0 else ""
                authRefAssignText += sp + tRef
            #
            refDetails = rsa.getDetails()
            if len(refDetails) > 0 and refDetails not in [".", "?"]:
                tDetails = "<b>Depositor details:&nbsp;</b> %s" % refDetails
                sp = spStr if len(authRefAssignText) > 0 else ""
                authRefAssignText += sp + tDetails
            #
        #
        if self.__verbose:
            self.__lfh.write("+SummaryView.__markupAuthorAssignments() depositor reference assignments for entity %s : %s\n" % (entityId, authRefAssignText))
        #
        if len(dbAccessionList) > 0:
            dbAccTaxMap = {}
            picklePath = os.path.join(self.__sessionPath, "dbAccTaxMap.pic")
            if os.access(picklePath, os.F_OK):
                fb = open(picklePath, "rb")
                dbAccTaxMap = pickle.load(fb)
                fb.close()
            #
            fetchSeqUtil = FetchSeqInfoUtils(siteId=str(self.__reqObj.getValue("WWPDB_SITE_ID")), verbose=self.__verbose, log=self.__lfh)
            tmpList = []
            for dbAccession in dbAccessionList:
                if dbAccession in dbAccTaxMap:
                    tmpList.append((dbAccession, dbAccTaxMap[dbAccession]))
                    continue
                #
                if dbAccession not in dbAccessionMap:
                    tmpList.append((dbAccession, ""))
                    continue
                #
                dbIsoform = ""
                dbAccessionS = dbAccession
                if dbAccessionMap[dbAccession] in ["UNP", "SP", "TR"]:
                    tL = dbAccession.split("-")
                    if len(tL) > 1:
                        dbIsoform = dbAccession
                        dbAccessionS = tL[0]
                    #
                #
                _accCode, refInfoD = fetchSeqUtil.getRefInfo(dbAccessionMap[dbAccession], dbAccessionS, dbIsoform, 0, 0, addMissingKeyFlag=False)
                if ("taxonomy_id" in refInfoD) and refInfoD["taxonomy_id"]:
                    tmpList.append((dbAccession, refInfoD["taxonomy_id"]))
                    dbAccTaxMap[dbAccession] = refInfoD["taxonomy_id"]
                else:
                    tmpList.append((dbAccession, ""))
                #
            #
            dbAccessionList = tmpList
            #
            fb = open(picklePath, "wb")
            pickle.dump(dbAccTaxMap, fb)
            fb.close()
        #
        return authRefAssignText, dbAccessionList

    def __checkPartRange(self, entityId, partInfoList):
        """Check part range definition and write out the warning message if there are errors"""
        seqLength = 0
        authVerList = self.__sds.getVersionIds(seqId=entityId, partId=1, altId=1, dataType="sequence", seqType="auth")
        if len(authVerList) > 0:
            seq = self.__sds.getSequence(seqId=entityId, seqType="auth", partId=1, altId=1, version=authVerList[0])
            seqLength = len(seq)
        #
        if seqLength == 0:
            return
        #
        ok, partText = self.__checkPartRangeError(seqLength, partInfoList)
        if ok:
            return
        #
        if self.__partRangeErrorMsg:
            self.__partRangeErrorMsg += "<br />\n"
        #
        self.__partRangeErrorMsg += "Entity '" + str(entityId) + "' (seq. length=" + str(seqLength) + ") : " + partText

    def __checkPartRangeError(self, seqLength, partList):
        """Check part range definition"""
        if len(partList) < 1:
            return False, "No part information defined."
        #
        try:
            status = True
            text = ""
            for i in range(0, len(partList)):
                if text:
                    text += ", "
                #
                text += "Part '" + str(partList[i][0]) + "' - ( " + str(partList[i][1]) + ", " + str(partList[i][2]) + " )"
                #
                seqNumBeg = int(partList[i][1])
                seqNumEnd = int(partList[i][2])
                if seqNumEnd < seqNumBeg:
                    status = False
                if (i == 0) and (seqNumBeg != 1):
                    status = False
                if (i == (len(partList) - 1)) and (seqNumEnd != seqLength):
                    status = False
                if i > 0:
                    prevNumEnd = int(partList[i - 1][2])
                    if seqNumBeg != (prevNumEnd + 1):
                        status = False
                    #
                #
            #
            return status, text
        except Exception as _e:  # noqa: F841
            text = ""
            for i in range(0, len(partList)):
                if text:
                    text += ", "
                #
                text += "Part '" + str(partList[i][0]) + "' - ( " + str(partList[i][1]) + ", " + str(partList[i][2]) + " )"
            #
            return False, text
        #


if __name__ == "__main__":
    pass
