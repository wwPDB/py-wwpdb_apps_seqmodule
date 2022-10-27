##
# File:  AlignmentDataStore.py
# Date:  25-Oct-2018
#
# Updates:
#   27-Oct-2022 zf   add getSummaryPageInfo() method
##
"""
Provide a storage interface for recording sequence alignments.
"""
__docformat__ = "restructuredtext en"
__author__ = "Zukang Feng"
__email__ = "zfeng@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.09"

try:
    import cPickle as pickle
except ImportError:
    import pickle as pickle
#
import copy
import os
import sys
import traceback
from operator import itemgetter

from wwpdb.apps.seqmodule.util.SequenceLabel import ResidueLabel, SequenceLabel, SequenceFeature
from wwpdb.apps.seqmodule.util.SequenceReferenceData import SequenceReferenceData
from wwpdb.apps.seqmodule.util.UpdateSequenceDataStoreUtils import UpdateSequenceDataStoreUtils
from wwpdb.io.locator.PathInfo import PathInfo


class AlignmentDataStore(UpdateSequenceDataStoreUtils):
    """Store and recover sequence alignment objects"""

    def __init__(self, reqObj=None, entityId=None, pathInfo=None, seqDataStore=None, deserializeFlag=True, verbose=False, log=sys.stderr):
        super(AlignmentDataStore, self).__init__(reqObj=reqObj, seqDataStore=seqDataStore, verbose=verbose, log=log)
        self._entityId = entityId
        self._pI = pathInfo
        self.__deserializeFlag = deserializeFlag
        #
        self._siteId = self._reqObj.getValue("WWPDB_SITE_ID")
        self._dataSetId = str(self._reqObj.getValue("identifier")).upper()
        self.__resLabel = ResidueLabel(verbose=self._verbose)
        self.__seqLabel = SequenceLabel(verbose=self._verbose)
        self._srd = SequenceReferenceData(verbose=self._verbose, log=self._lfh)
        self._gapSymbol = self._srd.getGapSymbol()
        #
        self.__errorMsg = ""
        self._authLabel = ""
        self._xyzAlignList = []
        self.__instIdMap = {}
        self.__partInfoDict = {}
        self.__xyzAlignIndexListMap = {}
        self.__xyzAlignIndexListIndices = {}
        self.__refAlignIndexListMap = {}
        self.__refAlignIndexListIndices = {}
        self._hasAlignmentFlag = False
        self._seqAlignLabelIndices = {}
        self._reverseSeqAlignLabelIndices = {}
        self._seqAlignList = []
        self._xyzAlignLabelIndices = {}
        self._missingSeqMap = {}
        self._authDefinedMutationList = []
        self._reverseXyzAlignLabelIndices = {}
        #
        self._clearAllPersistVariables()
        #
        self.__filePath = ""
        self.__setup()

    def _clearAllPersistVariables(self):
        """Clear all persist variables"""
        self._authLabel = ""
        self.__instIdMap = {}
        self.__partInfoDict = {}
        self.__xyzAlignIndexListMap = {}
        self.__xyzAlignIndexListIndices = {}
        self.__refAlignIndexListMap = {}
        self.__refAlignIndexListIndices = {}
        #
        self._resetAlignmentInfo()
        self._resetLogInfo()

    def _resetAlignmentInfo(self):
        """Clear current alignment info."""
        self._hasAlignmentFlag = False
        self._seqAlignLabelIndices = {}
        self._reverseSeqAlignLabelIndices = {}
        self._seqAlignList = []
        self._xyzAlignLabelIndices = {}
        self._reverseXyzAlignLabelIndices = {}
        self._xyzAlignList = []
        self._missingSeqMap = {}
        self._authDefinedMutationList = []

    def _resetLogInfo(self):
        """Clear current log info."""
        self.__errorMsg = ""

    def __setup(self):
        """ """
        try:
            self._sessionObj = self._reqObj.getSessionObj()
            self._sessionPath = self._sessionObj.getPath()
            #
            if not self._pI:
                self._pI = PathInfo(siteId=self._siteId, sessionPath=self._sessionPath, verbose=self._verbose, log=self._lfh)
            #
            self.__filePath = self._pI.getSequenceAlignFilePath(dataSetId=self._dataSetId, entityId=self._entityId, fileSource="session", versionId="latest")
            if self.__deserializeFlag and self.__filePath and os.access(self.__filePath, os.R_OK):
                self.deserialize()
            #
        except:  # noqa: E722 pylint: disable=bare-except
            if self._verbose:
                traceback.print_exc(file=self._lfh)
            #
            self._addErrorMessage(traceback.format_exc())
        #

    def getErrorMessage(self):
        """Return error message"""
        return self.__errorMsg

    def getSummaryPageInfo(self):
        """Return summary page info map"""
        if ("summary_page" in self._missingSeqMap) and self._missingSeqMap["summary_page"]:
            return self._missingSeqMap["summary_page"]
        #
        return {}

    def serialize(self):
        """Write alignment pickle file"""
        try:
            fb = open(self.__filePath, "wb")
            pickle.dump(self._authLabel, fb)
            pickle.dump(self._hasAlignmentFlag, fb)
            pickle.dump(self._seqAlignLabelIndices, fb)
            pickle.dump(self._reverseSeqAlignLabelIndices, fb)
            pickle.dump(self._seqAlignList, fb)
            pickle.dump(self._xyzAlignLabelIndices, fb)
            pickle.dump(self._reverseXyzAlignLabelIndices, fb)
            pickle.dump(self._xyzAlignList, fb)
            pickle.dump(self._missingSeqMap, fb)
            pickle.dump(self._authDefinedMutationList, fb)
            pickle.dump(self.__instIdMap, fb)
            pickle.dump(self.__partInfoDict, fb)
            pickle.dump(self.__xyzAlignIndexListMap, fb)
            pickle.dump(self.__xyzAlignIndexListIndices, fb)
            pickle.dump(self.__refAlignIndexListMap, fb)
            pickle.dump(self.__refAlignIndexListIndices, fb)
            pickle.dump(self.__errorMsg, fb)
            fb.close()
        except:  # noqa: E722 pylint: disable=bare-except
            if self._verbose:
                traceback.print_exc(file=self._lfh)
            #
            self._addErrorMessage(traceback.format_exc())
        #

    def updateEntityId(self, oldEntityId):
        """Update entityId for all seqIds"""
        self._authLabel = self.__replaceEntityId(self._authLabel, oldEntityId)
        #
        instIdMap = {}
        for k, v in self.__instIdMap.items():
            instIdMap[self.__replaceEntityId(k, oldEntityId)] = self.__replaceEntityId(v, oldEntityId)
        #
        self.__instIdMap = instIdMap
        #
        partInfoDict = {}
        for k, v in self.__partInfoDict.items():
            partInfoDict[self.__replaceEntityId(k, oldEntityId)] = v
        #
        self.__partInfoDict = partInfoDict
        #
        self.__xyzAlignIndexListIndices = self.__updateAlignIndexListIndices(self.__xyzAlignIndexListIndices, oldEntityId)
        self.__refAlignIndexListIndices = self.__updateAlignIndexListIndices(self.__refAlignIndexListIndices, oldEntityId)
        #
        self._seqAlignLabelIndices = self.__updateAlignLabelIndices(self._seqAlignLabelIndices, oldEntityId)
        self._reverseSeqAlignLabelIndices = self.__updateReverseAlignLabelIndices(self._reverseSeqAlignLabelIndices, oldEntityId)
        self._xyzAlignLabelIndices = self.__updateAlignLabelIndices(self._xyzAlignLabelIndices, oldEntityId)
        self._reverseXyzAlignLabelIndices = self.__updateReverseAlignLabelIndices(self._reverseXyzAlignLabelIndices, oldEntityId)

    def deserialize(self):
        """Read alignment pickle file"""
        try:
            fb = open(self.__filePath, "rb")
            self._authLabel = pickle.load(fb)
            self._hasAlignmentFlag = pickle.load(fb)
            self._seqAlignLabelIndices = pickle.load(fb)
            self._reverseSeqAlignLabelIndices = pickle.load(fb)
            self._seqAlignList = pickle.load(fb)
            self._xyzAlignLabelIndices = pickle.load(fb)
            self._reverseXyzAlignLabelIndices = pickle.load(fb)
            self._xyzAlignList = pickle.load(fb)
            self._missingSeqMap = pickle.load(fb)
            self._authDefinedMutationList = pickle.load(fb)
            self.__instIdMap = pickle.load(fb)
            self.__partInfoDict = pickle.load(fb)
            self.__xyzAlignIndexListMap = pickle.load(fb)
            self.__xyzAlignIndexListIndices = pickle.load(fb)
            self.__refAlignIndexListMap = pickle.load(fb)
            self.__refAlignIndexListIndices = pickle.load(fb)
            self.__errorMsg = pickle.load(fb)
            fb.close()
        except:  # noqa: E722 pylint: disable=bare-except
            if self._verbose:
                traceback.print_exc(file=self._lfh)
            #
            self._clearAllPersistVariables()
            # self._addErrorMessage(traceback.format_exc())
        #

    def _addErrorMessage(self, errMsg):
        """Concatenate error message to self._errorMsg"""
        if not errMsg:
            return
        #
        if self.__errorMsg:
            self.__errorMsg += "\n"
        #
        self.__errorMsg += errMsg

    def _insertAssociatedInstIds(self, newId, oldId):
        """Insert new/old instId relationship"""
        self.__instIdMap[newId] = oldId

    def _insertPartInfoDict(self, authId, partInfo):
        """Insert auth sequence parts information"""
        self.__partInfoDict[authId] = partInfo

    def _getPartInfoDict(self, authId):
        """Get auth sequence parts information"""
        if authId in self.__partInfoDict:
            return self.__partInfoDict[authId]
        else:
            orgId = self.__searchOriginalId(authId)
            if orgId and (orgId in self.__partInfoDict):
                return self.__partInfoDict[orgId]
            #
        #
        return {}

    def _insertXyzAlignIndexList(self, alignIdList, alignIndexList):
        """Insert xyz align index list info."""
        key = "|".join(sorted(alignIdList))
        if key in self.__xyzAlignIndexListIndices:
            self.__xyzAlignIndexListMap[self.__xyzAlignIndexListIndices[key][0]] = alignIndexList
            return
        #
        pos = len(self.__xyzAlignIndexListMap)
        self.__xyzAlignIndexListMap[pos] = alignIndexList
        self.__xyzAlignIndexListIndices[key] = (pos, alignIdList)

    def _getXyzAlignIndexList(self, alignIdList):
        """Get xyz align index list info."""
        _oldIdList, alignIndexList = self.__getAlignIndexList(alignIdList, self.__xyzAlignIndexListMap, self.__xyzAlignIndexListIndices)
        return alignIndexList

    def _insertRefAlignIndexList(self, alignIdList, alignIndexList):
        """Insert ref align index list info."""
        key = "|".join(sorted(alignIdList))
        if key in self.__refAlignIndexListIndices:
            self.__refAlignIndexListMap[self.__refAlignIndexListIndices[key][0]] = alignIndexList
            return
        #
        pos = len(self.__refAlignIndexListMap)
        self.__refAlignIndexListMap[pos] = alignIndexList
        self.__refAlignIndexListIndices[key] = (pos, alignIdList)

    def _getRefAlignIndexList(self, alignIdList):
        """Get ref align index list info."""
        _oldIdList, alignIndexList = self.__getAlignIndexList(alignIdList, self.__refAlignIndexListMap, self.__refAlignIndexListIndices)
        return alignIndexList

    def _copyXyzAlignmentToSeqAlignment(self):
        """Copy auth/coordinate alignment(s) to overall alignment(s)"""
        self._seqAlignLabelIndices = copy.deepcopy(self._xyzAlignLabelIndices)
        self._reverseSeqAlignLabelIndices = copy.deepcopy(self._reverseXyzAlignLabelIndices)
        self._seqAlignList = copy.deepcopy(self._xyzAlignList)

    def _getResLabelFromResLabelId(self, resLabelId):
        """Get ResidueLabel object"""
        if self.__resLabel.unpack(resLabelId):
            return self.__resLabel
        #
        return None

    def _getResLabelId(self, seqType="ref", seqInstId="", seqAltId=1, seqVersion=1, residueCode3="", residueLabelIndex=0, alignIndex=0, seqIndex=0, residueType="AA", seqPartId=1):
        """Get ResidueLabel Id"""
        self.__resLabel.set(
            seqType=seqType,
            seqInstId=seqInstId,
            seqAltId=seqAltId,
            seqVersion=seqVersion,
            residueCode3=residueCode3,
            residueLabelIndex=residueLabelIndex,
            alignIndex=alignIndex,
            seqIndex=seqIndex,
            residueType=residueType,
            seqPartId=seqPartId,
        )
        #
        return self.__resLabel.pack()

    def _getUnpackSeqLabel(self, seqLabelId):
        """Gte unpack sequence label"""
        if self.__seqLabel.unpack(seqLabelId):
            return self.__seqLabel.get()
        else:
            return ("", "", 1, 1, 1)
        #

    def _getSeqLabelId(self, tup):
        """SequenceLabel Id"""
        seqType, seqInstId, seqPartId, seqAltId, seqVersion = tup
        self.__seqLabel.set(seqType=seqType, seqInstId=seqInstId, seqPartId=seqPartId, seqAltId=seqAltId, seqVersion=seqVersion)
        return self.__seqLabel.pack()

    def _getSequenceByPackLabelFromDataStore(self, seqLabelId):
        """Get sequence from SequenceDataStore object"""
        (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = self._getUnpackSeqLabel(seqLabelId)
        if seqType:
            return self._getSequenceByUnpackLabelFromDataStore(seqType, seqInstId, seqPartId, seqAltId, seqVersion)
        else:
            return []
        #

    def _getSequenceByUnpackLabelFromDataStore(self, seqType, seqInstId, seqPartId, seqAltId, seqVersion):
        """Get sequence from SequenceDataStore object"""
        seqs = self.getSequence(seqType, seqInstId, seqPartId, seqAltId, seqVersion)
        if not seqs:
            label = "entity " + seqInstId
            if seqType == "xyz":
                label = "chain " + seqInstId
            elif seqType == "ref":
                label = "reference " + seqInstId
            #
            self._addErrorMessage("Failed to get sequence for '" + label + "' from SequencDataStore object.")
        #
        return seqs

    def _getFeatureByPackLabelFromDataStore(self, seqLabelId):
        """Get feature from SequenceDataStore object"""
        (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = self._getUnpackSeqLabel(seqLabelId)
        if seqType:
            return self.getFeature(seqType, seqInstId, seqPartId, seqAltId, seqVersion)
        else:
            return {}
        #

    def _getFeatureObjByPackLabelFromDataStore(self, seqLabelId):
        """Get feature object from SequenceDataStore object"""
        sfObj = SequenceFeature()
        sfDic = self._getFeatureByPackLabelFromDataStore(seqLabelId)
        if sfDic:
            sfObj.set(sfDic)
        #
        return sfObj

    def _getFeatureObjByUnpackLabelFromDataStore(self, seqType, seqInstId, seqPartId, seqAltId, seqVersion):
        """Get feature object from SequenceDataStore object"""
        sfObj = SequenceFeature()
        sfDic = self.getFeature(seqType, seqInstId, seqPartId, seqAltId, seqVersion)
        if sfDic:
            sfObj.set(sfDic)
        #
        return sfObj

    def _getProperAlignIdList(self, inputIdList, allSelectedIdList):
        """Check the input alignIds and return proper proper alignIds for doing alignment"""
        chainIdList = self.getGroup(self._entityId)
        #
        inputAuthList, inputXyzMap, inputRefMap, _inputSelfMap = self.__getSubSeqList(inputIdList, chainIdList)
        #
        selectAuthList, selectXyzMap, selectRefMap, selectSelfMap = self.__getSubSeqList(allSelectedIdList, chainIdList)
        #
        currXyzMap, currRefMap = self.__getCurrentSubSeqList()
        #
        authSeqId, extraAuthIdList = self.__getSeqId("auth", self._entityId, 1, self._authLabel, inputAuthList, selectAuthList)
        (seqType, seqInstId, _seqPartId, _seqAltId, _seqVersion) = self._getUnpackSeqLabel(authSeqId)
        #
        partIdList = self.getPartIdList(seqType, seqInstId)
        #
        retInputIdList = [authSeqId]
        selectIdList = [authSeqId]
        allAlignIdList = [authSeqId]
        xyzAlignList = [authSeqId]
        if len(extraAuthIdList) > 0:
            retInputIdList.extend(extraAuthIdList)
            allAlignIdList.extend(extraAuthIdList)
        #
        for chainId in chainIdList:
            iIdList = []
            if chainId in inputXyzMap:
                iIdList = inputXyzMap[chainId]
            #
            sIdList = []
            if chainId in selectXyzMap:
                sIdList = selectXyzMap[chainId]
            #
            currId = ""
            if chainId in currXyzMap:
                currId = currXyzMap[chainId]
            #
            xyzId, extraXyzList = self.__getSeqId("xyz", chainId, 1, currId, iIdList, sIdList)
            if len(xyzId) == 0:
                continue
            #
            if xyzId in iIdList:
                retInputIdList.append(xyzId)
            #
            selectIdList.append(xyzId)
            allAlignIdList.append(xyzId)
            xyzAlignList.append(xyzId)
            if len(extraXyzList) > 0:
                allAlignIdList.extend(extraXyzList)
                retInputIdList.extend(extraXyzList)
                xyzAlignList.extend(extraXyzList)
            #
        #
        refAlignList = []
        selfReflist = []
        for partId in partIdList:
            if str(partId) in selectSelfMap:
                selfReflist.append(selectSelfMap[str(partId)][0])
                continue
            #
            iIdList = []
            if str(partId) in inputRefMap:
                iIdList = inputRefMap[str(partId)]
            #
            sIdList = []
            if str(partId) in selectRefMap:
                sIdList = selectRefMap[str(partId)]
            #
            currId = ""
            if str(partId) in currRefMap:
                currId = currRefMap[str(partId)]
            #
            refId, extraRefList = self.__getSeqId("ref", self._entityId, partId, currId, iIdList, sIdList)
            if len(refId) == 0:
                continue
            #
            if refId in iIdList:
                retInputIdList.append(refId)
            #
            selectIdList.append(refId)
            allAlignIdList.append(refId)
            refAlignList.append(refId)
            if len(extraRefList) > 0:
                allAlignIdList.extend(extraRefList)
                retInputIdList.extend(extraRefList)
                refAlignList.extend(extraRefList)
            #
        #
        return retInputIdList, selectIdList, allAlignIdList, xyzAlignList, refAlignList, selfReflist, extraAuthIdList

    def _updateDefaultSelections(self, selectedIdList):
        """Update default selections in SequenceDataStore object for this entity"""
        chainIdList = self.getGroup(self._entityId)
        #
        defaultSelectIdList = self.getSelectedIds()
        #
        updatedSelectIdList = []
        for seqId in defaultSelectIdList:
            tL = str(seqId).strip().split("_")
            if len(tL) < 3:
                continue
            #
            elif (tL[0] in ("selfref", "auth", "ref")) and (tL[1] == self._entityId):
                continue
            elif (tL[0] == "xyz") and (tL[1] in chainIdList):
                continue
            #
            updatedSelectIdList.append(seqId)
        #
        updatedSelectIdList.extend(selectedIdList)
        self.setSelectedIds(updatedSelectIdList)
        self.saveSequenceDataStore()

    def __replaceEntityId(self, seqId, oldEntityId):
        """ """
        tL = str(seqId).split("_")
        if (len(tL) > 2) and (tL[0] in ("selfref", "auth", "ref")) and (tL[1] == oldEntityId):
            tL[1] = self._entityId
            return "_".join(tL)
        #
        return seqId

    def __updateAlignIndexListIndices(self, oldAlignIndexListIndices, oldEntityId):
        """ """
        newAlignIndexListIndices = {}
        for _k, tupL in oldAlignIndexListIndices.items():
            alignIdList = []
            for alignId in tupL[1]:
                alignIdList.append(self.__replaceEntityId(alignId, oldEntityId))
            #
            key = "|".join(sorted(alignIdList))
            newAlignIndexListIndices[key] = (tupL[0], alignIdList)
        #
        return newAlignIndexListIndices

    def __updateAlignLabelIndices(self, oldIndices, oldEntityId):
        """ """
        newIndices = {}
        for k, v in oldIndices.items():
            newIndices[self.__replaceEntityId(k, oldEntityId)] = v
        #
        return newIndices

    def __updateReverseAlignLabelIndices(self, oldIndices, oldEntityId):
        """ """
        newIndices = {}
        for k, v in oldIndices.items():
            newIndices[k] = self.__replaceEntityId(v, oldEntityId)
        #
        return newIndices

    def __searchOriginalId(self, newId):
        """Search associated previoud instId"""
        if newId not in self.__instIdMap:
            return newId
        else:
            return self.__searchOriginalId(self.__instIdMap[newId])
        #

    def __getAlignIndexList(self, alignIdList, alignIndexListMap, alignIndexListIndices):
        """Get align index list info."""
        if not alignIdList:
            return [], []
        #
        # Search by key
        #
        key = "|".join(sorted(alignIdList))
        if key in alignIndexListIndices:
            return alignIndexListIndices[key][1], alignIndexListMap[alignIndexListIndices[key][0]]
        #
        # Search with current instIds
        #
        for key, tupL in alignIndexListIndices.items():
            found = True
            for instId in alignIdList:
                if key.find(instId) == -1:
                    found = False
                    break
                #
            #
            if found:
                return tupL[1], alignIndexListMap[tupL[0]]
            #
        #
        # Search with previous instIds
        #
        for key, tupL in alignIndexListIndices.items():
            found = True
            idListMap = {}
            for instId in alignIdList:
                orgId = self.__searchOriginalId(instId)
                if key.find(orgId) == -1:
                    found = False
                    break
                #
                if orgId != instId:
                    idListMap[orgId] = instId
                #
            #
            if found:
                orgIdList = []
                for instId in tupL[1]:
                    if instId in idListMap:
                        orgIdList.append(idListMap[instId])
                    else:
                        orgIdList.append(instId)
                    #
                #
                return orgIdList, alignIndexListMap[tupL[0]]
            #
        #
        return [], []

    def __getSubSeqList(self, inputIdList, chainIdList):
        """ """
        sortAuthList = []
        xyzMap = {}
        refMap = {}
        selfMap = {}
        for seqId in inputIdList:
            tL = str(seqId).strip().split("_")
            if len(tL) < 3:
                continue
            #
            if (tL[0] == "selfref") and (tL[1] == self._entityId):
                if tL[2] in selfMap:
                    selfMap[tL[2]].append(seqId)
                else:
                    selfMap[tL[2]] = [seqId]
                #
                continue
            #
            seq = self._getSequenceByPackLabelFromDataStore(seqId)
            if len(seq) == 0:
                continue
            #
            if (tL[0] == "auth") and (tL[1] == self._entityId):
                sortAuthList.append((int(tL[4]), seqId))
            elif (tL[0] == "ref") and (tL[1] == self._entityId):
                if tL[2] in refMap:
                    refMap[tL[2]].append(seqId)
                else:
                    refMap[tL[2]] = [seqId]
                #
            elif (tL[0] == "xyz") and (tL[1] in chainIdList):
                if tL[1] in xyzMap:
                    xyzMap[tL[1]].append(seqId)
                else:
                    xyzMap[tL[1]] = [seqId]
                #
            #
        #
        if len(sortAuthList) > 1:
            sortAuthList.sort(key=itemgetter(0), reverse=True)
        #
        authList = []
        for sortTup in sortAuthList:
            authList.append(sortTup[1])
        #
        return authList, xyzMap, refMap, selfMap

    def __getCurrentSubSeqList(self):
        """ """
        xyzMap = {}
        refMap = {}
        for seqId in self._seqAlignLabelIndices.keys():
            tL = str(seqId).strip().split("_")
            if tL[0] == "ref":
                if tL[2] not in refMap:
                    refMap[tL[2]] = seqId
                #
            elif tL[0] == "xyz":
                if tL[1] not in xyzMap:
                    xyzMap[tL[1]] = seqId
                #
            #
        #
        return xyzMap, refMap

    def __getSeqId(self, seqType, seqInstId, seqPartId, currId, inputAuthList, selectAuthList):
        """ """
        if len(selectAuthList) == 1:
            return selectAuthList[0], self._substractedIdList(selectAuthList[0], inputAuthList)
        elif len(selectAuthList) > 0:
            if (len(inputAuthList) == 1) and (inputAuthList[0] in selectAuthList):
                return inputAuthList[0], inputAuthList[1:]
            else:
                return selectAuthList[0], self._substractedIdList(selectAuthList[0], inputAuthList)
            #
        elif len(inputAuthList) > 0:
            return inputAuthList[0], inputAuthList[1:]
        #
        if currId:
            return currId, self._substractedIdList(currId, inputAuthList)
        #
        altIdList = self.getAlternativeIdList(seqType, seqInstId, seqPartId)
        if len(altIdList) == 0:
            return "", inputAuthList
        #
        verList = self.getVersionIdList(seqType, seqInstId, seqPartId, altIdList[0])
        if len(verList) == 0:
            return "", inputAuthList
        #
        seqId = self._getSeqLabelId((seqType, seqInstId, seqPartId, altIdList[0], verList[0]))
        return seqId, self._substractedIdList(seqId, inputAuthList)

    def _substractedIdList(self, slectId, inputIdList):
        """ """
        if (len(inputIdList) == 0) or (slectId not in inputIdList):
            return inputIdList
        #
        sIdList = []
        for seqId in inputIdList:
            if seqId == slectId:
                continue
            #
            sIdList.append(seqId)
        #
        return sIdList
