##
# File:  AlignmentBackEndEditingTools.py
# Date:  10-Nov-2018
#
# Updates:
#  06-Jan-2023  zf  keep the reference sequence range information
#
"""
Methods to manage sequence alignment editing operations.
"""
__docformat__ = "restructuredtext en"
__author__ = "Zukang Feng"
__email__ = "zfeng@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.09"

import copy
import sys

from wwpdb.apps.seqmodule.align.AlignmentTools import AlignmentTools
from wwpdb.apps.seqmodule.align.AlignmentToolUtils import codeIndex, decodeIndex


class AlignmentBackEndEditingTools(AlignmentTools):
    """Manage sequence alignment editing"""

    def __init__(self, reqObj=None, entityId=None, pathInfo=None, seqDataStore=None, deserializeFlag=True, verbose=False, log=sys.stderr):
        super(AlignmentBackEndEditingTools, self).__init__(
            reqObj=reqObj, entityId=entityId, pathInfo=pathInfo, seqDataStore=seqDataStore, deserializeFlag=deserializeFlag, verbose=verbose, log=log
        )
        #
        self.__local_authLabel = self._authLabel
        self.__newPartInfoDict = {}
        self.__xyzAlignIndexList = []
        self.__editedType = []
        self.__editedSeqLabelMap = {}
        self.__newSeqLabelIdMap = {}

    def _update_details(self, resLabel, editInfoTuple):
        """Replace comment details: editInfoTuple ( 0: seqType, 1: seqLabelId, 2: alignTupleIndex, 3: editObj )"""
        self._seqAlignList[int(resLabel.getAlignmentIndex())][int(editInfoTuple[2])][5] = "ANNOTATED:" + str(editInfoTuple[3].getValueNew())
        #
        self.__updateEditingHistory(editInfoTuple[3].getEditType(), editInfoTuple[0], editInfoTuple[1])

    def _update_replace(self, resLabel, editInfoTuple):
        """Replace residue name: editInfoTuple ( 0: seqType, 1: seqLabelId, 2: alignTupleIndex, 3: editObj )"""
        newResidueName = editInfoTuple[3].getValueNew()[0]
        newOneLetterCode = self._srd.cnv3To1(newResidueName)
        if editInfoTuple[0] == "auth":
            for labelId, tupleIndex in self._xyzAlignLabelIndices.items():
                self._seqAlignList[int(resLabel.getAlignmentIndex())][int(tupleIndex)][6] = 0
                self._seqAlignList[int(resLabel.getAlignmentIndex())][int(tupleIndex)][7] = ""
                if labelId.startswith("xyz"):
                    if self._seqAlignList[int(resLabel.getAlignmentIndex())][int(tupleIndex)][1] == self._gapSymbol:
                        continue
                    #
                    if self._seqAlignList[int(resLabel.getAlignmentIndex())][int(tupleIndex)][1] != newResidueName:
                        self.__updateEditingHistory(editInfoTuple[3].getEditType(), "xyz", labelId)
                    #
                #
                self._seqAlignList[int(resLabel.getAlignmentIndex())][int(tupleIndex)][0] = newOneLetterCode
                self._seqAlignList[int(resLabel.getAlignmentIndex())][int(tupleIndex)][1] = newResidueName
            #
        else:
            self._seqAlignList[int(resLabel.getAlignmentIndex())][int(editInfoTuple[2])][0] = newOneLetterCode
            self._seqAlignList[int(resLabel.getAlignmentIndex())][int(editInfoTuple[2])][1] = newResidueName
        #
        self.__updateEditingHistory(editInfoTuple[3].getEditType(), editInfoTuple[0], editInfoTuple[1])

    def _update_replaceid(self, resLabel, editInfoTuple):
        """Exchange residue position: editInfoTuple ( 0: seqType, 1: seqLabelId, 2: alignTupleIndex, 3: editObj )"""
        newId = editInfoTuple[3].getNewElementId()
        newLabelObj = self._getResLabelFromResLabelId(newId)
        srcSeqPos = newLabelObj.getSequenceIndex()
        srcSeqNum = newLabelObj.getResidueLabelIndex()
        srcResNam = editInfoTuple[3].getValueNew()[0]
        oneLetterCode = self._srd.cnv3To1(srcResNam)
        #
        if editInfoTuple[0] == "xyz":
            if srcResNam != self._gapSymbol:
                for labelId, tupleIndex in self._xyzAlignLabelIndices.items():
                    if labelId.startswith("auth"):
                        if srcResNam != self._seqAlignList[int(resLabel.getAlignmentIndex())][int(tupleIndex)][1]:
                            oneLetterCode = self._seqAlignList[int(resLabel.getAlignmentIndex())][int(tupleIndex)][0]
                            srcResNam = self._seqAlignList[int(resLabel.getAlignmentIndex())][int(tupleIndex)][1]
                        #
                        break
                    #
                #
                try:
                    seq = self._getSequenceByPackLabelFromDataStore(editInfoTuple[1])
                    self._seqAlignList[int(resLabel.getAlignmentIndex())][int(editInfoTuple[2])][5] = seq[int(srcSeqPos)][2]
                except:  # noqa: E722 pylint: disable=bare-except
                    pass
                #
            else:
                self._seqAlignList[int(resLabel.getAlignmentIndex())][int(editInfoTuple[2])][5] = ""
            #
        #
        self._seqAlignList[int(resLabel.getAlignmentIndex())][int(editInfoTuple[2])][0] = oneLetterCode
        self._seqAlignList[int(resLabel.getAlignmentIndex())][int(editInfoTuple[2])][1] = srcResNam
        self._seqAlignList[int(resLabel.getAlignmentIndex())][int(editInfoTuple[2])][2] = srcSeqNum
        self._seqAlignList[int(resLabel.getAlignmentIndex())][int(editInfoTuple[2])][3] = srcSeqPos
        #
        self.__updateEditingHistory(editInfoTuple[3].getEditType(), editInfoTuple[0], editInfoTuple[1])
        #
        if (srcResNam != self._gapSymbol) and (editInfoTuple[0] == "auth"):
            for labelId, tupleIndex in self._xyzAlignLabelIndices.items():
                self._seqAlignList[int(resLabel.getAlignmentIndex())][int(tupleIndex)][6] = 0
                self._seqAlignList[int(resLabel.getAlignmentIndex())][int(tupleIndex)][7] = ""
                if labelId.startswith("xyz"):
                    if self._seqAlignList[int(resLabel.getAlignmentIndex())][int(tupleIndex)][1] == self._gapSymbol:
                        continue
                    #
                    if self._seqAlignList[int(resLabel.getAlignmentIndex())][int(tupleIndex)][1] != srcResNam:
                        self.__updateEditingHistory(editInfoTuple[3].getEditType(), "xyz", labelId)
                    #
                #
                self._seqAlignList[int(resLabel.getAlignmentIndex())][int(tupleIndex)][0] = oneLetterCode
                self._seqAlignList[int(resLabel.getAlignmentIndex())][int(tupleIndex)][1] = srcResNam
            #
        #

    def _update_insert(self, resLabel, editInfoTuple):
        """Insert residue(s): editInfoTuple ( 0: seqType, 1: seqLabelId, 2: alignTupleIndex, 3: editObj )"""
        alignList = []
        for alignTup in self._seqAlignList:
            if int(alignTup[int(editInfoTuple[2])][4]) == int(resLabel.getAlignmentIndex()):
                residueList = editInfoTuple[3].getValueNew()
                newResidueName = residueList[0]
                newOneLetterCode = self._srd.cnv3To1(newResidueName)
                alignTup[int(editInfoTuple[2])][0] = newOneLetterCode
                alignTup[int(editInfoTuple[2])][1] = newResidueName
                alignList.append(alignTup)
                for residueName in residueList[1:]:
                    oneLetterCode = self._srd.cnv3To1(residueName)
                    newAlignTup = []
                    for i in range(0, len(alignTup)):
                        if i == int(editInfoTuple[2]):
                            newAlignTup.append([oneLetterCode, residueName, "", "", "", "", 0, ""])
                        else:
                            newAlignTup.append([self._gapSymbol, self._gapSymbol, "", "", "", "", 0, ""])
                        #
                    #
                    alignList.append(newAlignTup)
                #
            else:
                alignList.append(alignTup)
            #
        #
        self._seqAlignList = alignList
        #
        self.__updateEditingHistory(editInfoTuple[3].getEditType(), editInfoTuple[0], editInfoTuple[1])

    def _update_delete(self, resLabel, editInfoTuple):
        """Delete residue: editInfoTuple ( 0: seqType, 1: seqLabelId, 2: alignTupleIndex, 3: editObj )"""
        self._seqAlignList[int(resLabel.getAlignmentIndex())][int(editInfoTuple[2])][0] = self._gapSymbol
        self._seqAlignList[int(resLabel.getAlignmentIndex())][int(editInfoTuple[2])][1] = self._gapSymbol
        self._seqAlignList[int(resLabel.getAlignmentIndex())][int(editInfoTuple[2])][2] = ""
        self._seqAlignList[int(resLabel.getAlignmentIndex())][int(editInfoTuple[2])][3] = ""
        self._seqAlignList[int(resLabel.getAlignmentIndex())][int(editInfoTuple[2])][5] = ""
        self._seqAlignList[int(resLabel.getAlignmentIndex())][int(editInfoTuple[2])][6] = 0
        self._seqAlignList[int(resLabel.getAlignmentIndex())][int(editInfoTuple[2])][7] = ""
        #
        self.__updateEditingHistory(editInfoTuple[3].getEditType(), editInfoTuple[0], editInfoTuple[1])

    def _updateAlignmentAndSequences(self, inputAlignIdList, inputSelectedIdList):
        """Consolidate alignment and update any sequence changes after alignment editing operation."""
        if len(self.__editedType) < 1:
            return
        #
        if "delete" in self.__editedType:
            reIndexingFlag = False
            if "insert" in self.__editedType:
                reIndexingFlag = True
            #
            self.__removeEmptyAlignRecord(reIndexingFlag)
        elif "insert" in self.__editedType:
            self.__reIndexingAlignment()
        #
        authIdx = self._seqAlignLabelIndices[self.__local_authLabel]
        self._checkPartStartEndPosMap(authIdx, self.__local_authLabel)
        #
        updatedSeqLabelIdList = []
        for seqType in ("auth", "xyz", "ref"):
            if seqType not in self.__editedSeqLabelMap:
                continue
            #
            for seqLabelId, editTypeList in self.__editedSeqLabelMap[seqType].items():
                self.__updateSequenceInfo(authIdx, seqLabelId, editTypeList, nextVersionFlag=True)
                updatedSeqLabelIdList.append(seqLabelId)
            #
        #
        for seqLabelId in inputAlignIdList:
            if seqLabelId in updatedSeqLabelIdList:
                continue
            #
            self.__updateSequenceInfo(authIdx, seqLabelId, "")
        #
        if ("details" in self.__editedType) or ("replace" in self.__editedType) or ("delete" in self.__editedType) or ("insert" in self.__editedType):
            self.saveSequenceDataStore()
        #
        outputAlignIdList, outputSelectedIdList = self.__updateSeqLabelIdWithNewVersion(inputAlignIdList, inputSelectedIdList)
        #
        if ("replaceid" in self.__editedType) or ("delete" in self.__editedType) or ("insert" in self.__editedType):
            self.__updateXyzAlignList()
            self.__updateRefAlignList()
        #
        # Insert new sequence parts information
        if len(self.__newPartInfoDict) > 0:
            self._insertPartInfoDict(self._authLabel, self.__newPartInfoDict)
        #
        # Insert new self.__xyzAlignIndexList
        if len(self.__xyzAlignIndexList) > 0:
            xyzLabelList = []
            for idx in sorted(self._reverseXyzAlignLabelIndices.keys()):
                xyzLabelList.append(self._reverseXyzAlignLabelIndices[idx])
            #
            self._insertXyzAlignIndexList(xyzLabelList, self.__xyzAlignIndexList)
        #
        return outputAlignIdList, outputSelectedIdList

    def __updateEditingHistory(self, edType, seqType, seqLabelId):
        """Record editing history"""
        if edType not in self.__editedType:
            self.__editedType.append(edType)
        #
        if seqType in self.__editedSeqLabelMap:
            if seqLabelId in self.__editedSeqLabelMap[seqType]:
                if edType not in self.__editedSeqLabelMap[seqType][seqLabelId]:
                    self.__editedSeqLabelMap[seqType][seqLabelId].append(edType)
                #
            else:
                self.__editedSeqLabelMap[seqType][seqLabelId] = [edType]
            #
        else:
            tmpD = {}
            tmpD[seqLabelId] = [edType]
            self.__editedSeqLabelMap[seqType] = tmpD
        #

    def __removeEmptyAlignRecord(self, reIndexingFlag):
        """Remove all empty align records after "delete" operation"""
        alignList = []
        for alignTup in self._seqAlignList:
            hasValue = False
            for alignRecord in alignTup:
                if alignRecord[0] and (alignRecord[0] != self._gapSymbol):
                    hasValue = True
                    break
                #
            #
            if hasValue:
                alignList.append(alignTup)
            #
        #
        hasChanged = False
        if len(alignList) != len(self._seqAlignList):
            self._seqAlignList = alignList
            hasChanged = True
        #
        if hasChanged or reIndexingFlag:
            self.__reIndexingAlignment()
        #

    def __reIndexingAlignment(self):
        """Re-indexing align position index after "delete" or "insert" operations"""
        for idx, alignTup in enumerate(self._seqAlignList):
            for alignRecord in alignTup:
                alignRecord[4] = idx
            #
        #

    def __updateSequenceInfo(self, authIdx, seqLabelId, editTypeList, nextVersionFlag=False):
        """Propogate the alignment edits to new versions of the sequences in the sequence data store."""
        if nextVersionFlag and ((len(editTypeList) < 1) or ((len(editTypeList) == 1) and (editTypeList[0] == "replaceid"))):
            return
        #
        (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = self._getUnpackSeqLabel(seqLabelId)
        orgSeqList = self._getSequenceByUnpackLabelFromDataStore(seqType, seqInstId, seqPartId, seqAltId, seqVersion)
        nextVersion = self.getNextVersionNumber(seqType, seqInstId, seqPartId, seqAltId)
        #
        begPos = 0
        endPos = len(self._seqAlignList) - 1
        if (seqType == "ref") and (int(seqPartId) in self._partPosDict):
            begPos = self._partPosDict[int(seqPartId)][0]
            endPos = self._partPosDict[int(seqPartId)][1]
        #
        alignIdx = self._seqAlignLabelIndices[seqLabelId]
        posPartIdDict = {}
        if (alignIdx == authIdx) and (("delete" in editTypeList) or ("insert" in editTypeList) or ("replace" in editTypeList)):
            for partId, posRange in self._partPosDict.items():
                for pos in range(posRange[0], posRange[1] + 1):
                    posPartIdDict[pos] = partId
                #
            #
        #
        alignLength = 0
        numMatch = 0
        numMatchGaps = 0
        #
        newSeqList = []
        for idx, alignTup in enumerate(self._seqAlignList):
            if (idx < begPos) or (idx > endPos):
                continue
            #
            if (alignTup[authIdx][1] == self._gapSymbol) and (alignTup[alignIdx][1] == self._gapSymbol):
                continue
            #
            alignLength += 1
            if alignTup[authIdx][1] == alignTup[alignIdx][1]:
                numMatch += 1
                numMatchGaps += 1
            if alignTup[alignIdx][1] == self._gapSymbol:
                numMatchGaps += 1
                continue
            #
            seqIdx = len(newSeqList)
            comment = ""
            if seqType == "auth":
                alignTup[alignIdx][2] = str(seqIdx + 1)
                intIdx, comment = decodeIndex(alignTup[alignIdx][3])
                if intIdx >= 0:
                    newSeqList.append((alignTup[alignIdx][1], alignTup[alignIdx][2], orgSeqList[intIdx][2], seqIdx + 1, alignTup[alignIdx][0], orgSeqList[intIdx][5]))
                else:
                    newSeqList.append((alignTup[alignIdx][1], alignTup[alignIdx][2], "", seqIdx + 1, alignTup[alignIdx][0], alignTup[alignIdx][1]))
                #
                if (alignIdx == authIdx) and (("delete" in editTypeList) or ("insert" in editTypeList) or ("replace" in editTypeList)):
                    if idx in posPartIdDict:
                        partId = posPartIdDict[idx]
                        if partId in self.__newPartInfoDict:
                            self.__newPartInfoDict[partId][1] = seqIdx + 1
                        else:
                            self.__newPartInfoDict[partId] = [seqIdx + 1, seqIdx + 1]
                        #
                    #
                #
            elif seqType == "xyz":
                intIdx, comment = decodeIndex(alignTup[alignIdx][3])
                if intIdx >= 0:
                    newSeqList.append(
                        (alignTup[alignIdx][1], alignTup[alignIdx][2], orgSeqList[intIdx][2], seqIdx + 1, alignTup[alignIdx][0], orgSeqList[intIdx][5], orgSeqList[intIdx][6])
                    )
                else:
                    newSeqList.append((alignTup[alignIdx][1], alignTup[alignIdx][2], "", seqIdx + 1, alignTup[alignIdx][0], alignTup[alignIdx][1], "0"))
                #
            else:
                newSeqList.append((alignTup[alignIdx][1], alignTup[alignIdx][2], "", seqIdx + 1, alignTup[alignIdx][0]))
            #
            alignTup[alignIdx][3] = codeIndex(seqIdx, comment)
        #
        if nextVersionFlag:
            self.setSequence(newSeqList, seqType, seqInstId, seqPartId, seqAltId, nextVersion)
        #
        if seqType == "auth":
            if alignIdx == authIdx:
                partIdList = self.getPartIdList(seqType, seqInstId)
                for partId in partIdList:
                    sFeature = self._getFeatureObjByUnpackLabelFromDataStore(seqType, seqInstId, partId, seqAltId, seqVersion)
                    sFeature.clearAlignDetails()
                    _jPartId, jSeqPartType = sFeature.getPartInfo()
                    if partId in self.__newPartInfoDict:
                        sFeature.setAuthPartDetails(partId, self.__newPartInfoDict[partId][0], self.__newPartInfoDict[partId][1], jSeqPartType)
                        #
                        altIdList = self.getAlternativeIdList("ref", seqInstId, partId)
                        for altId in altIdList:
                            fObj = self._getFeatureObjByUnpackLabelFromDataStore("ref", seqInstId, partId, altId, seqVersion)
                            fObj.setAuthPartDetails(partId, self.__newPartInfoDict[partId][0], self.__newPartInfoDict[partId][1], jSeqPartType)
                            self.setFeature(fObj.get(), "ref", seqInstId, partId, altId, seqVersion)
                        #
                    #
                    self.setFeature(sFeature.get(), seqType, seqInstId, partId, seqAltId, nextVersion)
                #
            #
        else:
            sFeature = self._getFeatureObjByUnpackLabelFromDataStore(seqType, seqInstId, seqPartId, seqAltId, seqVersion)
            matchBeg = sFeature.getItem("REF_MATCH_BEGIN")
            matchEnd = sFeature.getItem("REF_MATCH_END")
            sFeature.clearAlignDetails()
            if seqType == "xyz":
                sFeature.setAuthXyzAlignDetails(
                    seqLen=len(newSeqList), alignLen=int(alignLength), seqSim=float(numMatch) / float(alignLength), seqSimWithGaps=float(numMatchGaps) / float(alignLength)
                )
            else:
                sFeature.setAuthRefAlignDetails(
                    seqLen=len(newSeqList), alignLen=int(alignLength), seqSim=float(numMatch) / float(alignLength), seqSimWithGaps=float(numMatchGaps) / float(alignLength)
                )
                sFeature.setAuthRefAlignRange(refMatchBegin=matchBeg, refMatchEnd=matchEnd)
            #
            if nextVersionFlag:
                self.setFeature(sFeature.get(), seqType, seqInstId, seqPartId, seqAltId, nextVersion)
            else:
                self.setFeature(sFeature.get(), seqType, seqInstId, seqPartId, seqAltId, seqVersion)
            #
        #
        if nextVersionFlag:
            self.__newSeqLabelIdMap[seqLabelId] = self._getSeqLabelId((seqType, seqInstId, seqPartId, seqAltId, nextVersion))
        #

    def __updateXyzAlignList(self):
        """Re-create self._xyzAlignLabelIndices, self._reverseXyzAlignLabelIndices & self._xyzAlignList from
        self._seqAlignLabelIndices, self._reverseSeqAlignLabelIndices & self._seqAlignList
        """
        if self._authLabel not in self._seqAlignLabelIndices:
            return
        #
        indexList = []
        indexList.append((self._authLabel, self._seqAlignLabelIndices[self._authLabel]))
        for i in range(0, len(self._reverseSeqAlignLabelIndices)):
            seqLabelId = self._reverseSeqAlignLabelIndices[i]
            if not seqLabelId.startswith("xyz"):
                continue
            #
            indexList.append((seqLabelId, i))
        #
        self._xyzAlignLabelIndices = {}
        self._reverseXyzAlignLabelIndices = {}
        for idx, indexTuple in enumerate(indexList):
            self._xyzAlignLabelIndices[indexTuple[0]] = idx
            self._reverseXyzAlignLabelIndices[idx] = indexTuple[0]
        #
        self._xyzAlignList, self.__xyzAlignIndexList = self.__getAlignAndIndexList(inputSeqIndexList=indexList)

    def __getAlignAndIndexList(self, inputSeqIndexList=None, includeAlignmentFlag=True):
        """ """
        if inputSeqIndexList is None:
            inputSeqIndexList = []

        alignList = []
        alignIndexList = []
        for seqAlignTuple in self._seqAlignList:
            hasRealValue = False
            for indexTuple in inputSeqIndexList:
                intIdx, _comment = decodeIndex(seqAlignTuple[indexTuple[1]][3])
                if intIdx >= 0:
                    hasRealValue = True
                    break
                #
            #
            if not hasRealValue:
                continue
            #
            alignTuple = []
            alignIndexTuple = []
            for indexTuple in inputSeqIndexList:
                if includeAlignmentFlag:
                    alignTuple.append(copy.deepcopy(seqAlignTuple[indexTuple[1]]))
                #
                alignIndexTuple.append(seqAlignTuple[indexTuple[1]][3])
            #
            if includeAlignmentFlag:
                alignList.append(alignTuple)
            #
            alignIndexList.append(alignIndexTuple)
        #
        return alignList, alignIndexList

    def __updateRefAlignList(self):
        """ """
        if self._authLabel not in self._seqAlignLabelIndices:
            return
        #
        for i in range(0, len(self._reverseSeqAlignLabelIndices)):
            seqLabelId = self._reverseSeqAlignLabelIndices[i]
            if not seqLabelId.startswith("ref"):
                continue
            #
            indexList = [(self._authLabel, self._seqAlignLabelIndices[self._authLabel]), (seqLabelId, i)]
            #
            _alignList, refAlignIndexList = self.__getAlignAndIndexList(inputSeqIndexList=indexList, includeAlignmentFlag=False)
            if refAlignIndexList:
                self._insertRefAlignIndexList([self._authLabel, seqLabelId], refAlignIndexList)
            #
        #

    def __updateSeqLabelIdWithNewVersion(self, inputAlignIdList, inputSelectedIdList):
        """Update seqLabelId with next version if exists"""
        if not self.__newSeqLabelIdMap:
            return inputAlignIdList, inputSelectedIdList
        #
        # Insert new/old instId relationship
        #
        for oldId, newId in self.__newSeqLabelIdMap.items():
            self._insertAssociatedInstIds(newId, oldId)
        #
        # update seqIds with new version
        #
        if self.__local_authLabel in self.__newSeqLabelIdMap:
            self._authLabel = self.__newSeqLabelIdMap[self.__local_authLabel]
        #
        self._seqAlignLabelIndices, self._reverseSeqAlignLabelIndices = self.__reAssignAlignLabelIndices(self._seqAlignLabelIndices)
        self._xyzAlignLabelIndices, self._reverseXyzAlignLabelIndices = self.__reAssignAlignLabelIndices(self._xyzAlignLabelIndices)
        #
        # update input align ID list
        outputAlignIdList = self.__reAssignAlignIdList(inputAlignIdList)
        outputSelectedIdList = self.__reAssignAlignIdList(inputSelectedIdList)
        return outputAlignIdList, outputSelectedIdList

    def __reAssignAlignLabelIndices(self, inAlignLabelIndices):
        """Re-assign alignment seqLabelId indices with next version"""
        outAlignLabelIndices = {}
        outReverseAlignLabelIndices = {}
        for seqLabelId, index in inAlignLabelIndices.items():
            if seqLabelId in self.__newSeqLabelIdMap:
                outAlignLabelIndices[self.__newSeqLabelIdMap[seqLabelId]] = index
                outReverseAlignLabelIndices[index] = self.__newSeqLabelIdMap[seqLabelId]
            else:
                outAlignLabelIndices[seqLabelId] = index
                outReverseAlignLabelIndices[index] = seqLabelId
            #
        #
        return outAlignLabelIndices, outReverseAlignLabelIndices

    def __reAssignAlignIdList(self, inputSeqLabelIdList):
        """Update seqLabelId list with next version if exists"""
        outputSeqLabelIdList = []
        for seqLabelId in inputSeqLabelIdList:
            if seqLabelId in self.__newSeqLabelIdMap:
                outputSeqLabelIdList.append(self.__newSeqLabelIdMap[seqLabelId])
            else:
                outputSeqLabelIdList.append(seqLabelId)
            #
        #
        return outputSeqLabelIdList
