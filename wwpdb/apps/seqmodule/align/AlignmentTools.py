##
# File:  AlignmentTools.py
# Date:  25-Oct-2018
#
# Updates:
# 29-Sep-2022: zf added 'auth_numbering' as 4th value in seqTuple for coodinate sequence in __generateInputSeqInfoForPseudoMultiAlignFromSeqList() method.
#              added self.__conflictMap for new entity summary page.
#
"""
Methods to manage sequence alignments.
"""
__docformat__ = "restructuredtext en"
__author__ = "Zukang Feng"
__email__ = "zfeng@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.09"

import copy
import os
import sys
import traceback

from wwpdb.apps.seqmodule.align.AlignmentToolUtils import getSeqAlignment, mergeSeqAlignment, codeSeqIndex, decodeIndex, assignConflict
from wwpdb.apps.seqmodule.io.AlignmentDataStore import AlignmentDataStore
from wwpdb.utils.align.alignlib import PseudoMultiAlign  # pylint: disable=no-name-in-module
from wwpdb.io.file.mmCIFUtil import mmCIFUtil
from wwpdb.io.locator.ChemRefPathInfo import ChemRefPathInfo


class AlignmentTools(AlignmentDataStore):
    """Manage sequence alignments"""

    def __init__(self, reqObj=None, entityId=None, pathInfo=None, seqDataStore=None, deserializeFlag=True, verbose=False, log=sys.stderr):
        super(AlignmentTools, self).__init__(
            reqObj=reqObj, entityId=entityId, pathInfo=pathInfo, seqDataStore=seqDataStore, deserializeFlag=deserializeFlag, verbose=verbose, log=log
        )
        #
        self._longExpressionTagCountCutoff = 20
        self.__clearLocalVariables()
        self.__local_authLabel = self._authLabel
        #
        self.__standardList = (
            "ALA",
            "ARG",
            "ASN",
            "ASP",
            "CYS",
            "GLN",
            "GLU",
            "GLY",
            "HIS",
            "ILE",
            "LEU",
            "LYS",
            "MET",
            "PHE",
            "PRO",
            "PYL",
            "SEC",
            "SER",
            "THR",
            "TRP",
            "TYR",
            "VAL",
            "UNK",
            "A",
            "C",
            "G",
            "I",
            "T",
            "U",
            "DA",
            "DC",
            "DG",
            "DI",
            "DT",
            "DU",
        )
        #
        self.__parentCompIdMap = {}
        self.__refSeqVariantList = []
        self.__alignFlag = True
        self.__refLabelDict = {}
        self.__partInfoDict = {}
        self.__refAlignIndexDict = {}
        self.__selectedIdList = []
        self.__xyzLabel = []
        self.__newXyzAlignIndexListFlag = False
        self.__xyzAlignIndexList = []
        self.__extraAuthLabel = []
        # Should be safe to initialize here - not in base class
        self._partPosDict = {}
        self._selfRefPartIdList = []

    def __clearLocalVariables(self):
        """Clear all local variables"""
        self.__local_authLabel = ""
        self.__alignFlag = True
        self.__selectedIdList = []
        self.__xyzLabel = []
        self.__newXyzAlignIndexListFlag = False
        self.__xyzAlignIndexList = []
        self.__partInfoDict = {}
        self.__refLabelDict = {}
        self.__refAlignIndexDict = {}
        self.__extraAuthLabel = []
        self._selfRefPartIdList = []
        self._partPosDict = {}

    def setInputAlignData(self, alignD=None):
        """Input entity ID, auth, xyz & ref labels, and align indices data among auth seq vs coordinate seqs and auth seq vs ref seqs. etc

        alignD["alignids"] = [ "auth_label", "xyz_label_1", "xyz_label_2", ..., "ref_label_1", "ref_label_2", ... ]
        alignD["xyz_label"] = [ "xyz_label_1", "xyz_label_2", ...]
        alignD["auth_coord_align_index"] = [ [...], [...], ...]
        alignD["part_info"] = { partId_1 : [ beg_num_1, end_num_1], partId_2 : [ beg_num_2, end_num_2], ... }
        alignD["ref_label"] = { partId_1 : ["ref_label_1", "ref_label_2", ...], partId_2 : ["ref_label_3", "ref_label_4", ...], ... }
        alignD["auth_ref_align_index"] = { "ref_label_1" : [ [...], [...], ...], "ref_label_2" : [ [...], [...], ...], ... }
        """
        try:
            self._clearAllPersistVariables()
            self.__clearLocalVariables()
            #
            self.__local_authLabel = alignD["auth_label"]
            self.__selectedIdList = alignD["alignids"]
            self.__xyzLabel = alignD["xyz_label"]
            self.__xyzAlignIndexList = alignD["auth_coord_align_index"]
            self.__partInfoDict = alignD["part_info"]
            self.__refLabelDict = alignD["ref_label"]
            self.__refAlignIndexDict = alignD["auth_ref_align_index"]
            #
            self.__doAllAlignments()
            if self.__alignFlag:
                self._assignAllConflicts(self.__local_authLabel, self.__selectedIdList)
                #
                if (len(self.__xyzLabel) > 0) and (len(self.__xyzAlignIndexList) > 0):
                    self._insertXyzAlignIndexList(self.__xyzLabel, self.__xyzAlignIndexList)
                #
                if len(self.__partInfoDict) > 0:
                    self._insertPartInfoDict(self.__local_authLabel, self.__partInfoDict)
                #
                if len(self.__refAlignIndexDict) > 0:
                    for redId, alignIndexList in self.__refAlignIndexDict.items():
                        self._insertRefAlignIndexList([self.__local_authLabel, redId], alignIndexList)
                    #
                #
            #
            self._updateAlignmentDataStore(self.__local_authLabel)
        except:  # noqa: E722 pylint: disable=bare-except
            if self._verbose:
                traceback.print_exc(file=self._lfh)
            #
            self._addErrorMessage(traceback.format_exc())
        #

    def resetInputAlignData(self, alignD=None):
        """Input entity ID, auth, xyz & ref labels, and align indices data among auth seq vs coordinate seqs and auth seq vs ref seqs. etc

        alignD["alignids"] = [ "auth_label", "xyz_label_1", "xyz_label_2", ..., "ref_label_1", "ref_label_2", ... ]
        alignD["xyz_label"] = [ "xyz_label_1", "xyz_label_2", ...]
        alignD["part_info"] = { partId_1 : [ beg_num_1, end_num_1], partId_2 : [ beg_num_2, end_num_2], ... }
        alignD["ref_label"] = { partId_1 : ["ref_label_1", "ref_label_2", ...], partId_2 : ["ref_label_3", "ref_label_4", ...], ... }
        alignD["auth_ref_align_index"] = { "ref_label_1" : [ [...], [...], ...], "ref_label_2" : [ [...], [...], ...], ... }
        """
        try:
            self.__clearLocalVariables()
            #
            self.__local_authLabel = alignD["auth_label"]
            self.__selectedIdList = alignD["alignids"]
            self.__xyzLabel = alignD["xyz_label"]
            self.__partInfoDict = alignD["part_info"]
            self.__refLabelDict = alignD["ref_label"]
            self.__refAlignIndexDict = alignD["auth_ref_align_index"]
            #
            self.__resetAllAlignments()
            #
            if self.__alignFlag:
                self._assignAllConflicts(self.__local_authLabel, self.__selectedIdList)
                #
                if (len(self.__xyzLabel) > 0) and (len(self.__xyzAlignIndexList) > 0):
                    self._insertXyzAlignIndexList(self.__xyzLabel, self.__xyzAlignIndexList)
                #
                if len(self.__partInfoDict) > 0:
                    self._insertPartInfoDict(self.__local_authLabel, self.__partInfoDict)
                #
                if len(self.__refAlignIndexDict) > 0:
                    for redId, alignIndexList in self.__refAlignIndexDict.items():
                        self._insertRefAlignIndexList([self.__local_authLabel, redId], alignIndexList)
                    #
                #
                self._updateAlignmentDataStore(self.__local_authLabel)
            #
        except:  # noqa: E722 pylint: disable=bare-except
            if self._verbose:
                traceback.print_exc(file=self._lfh)
            #
            self._addErrorMessage(traceback.format_exc())
        #

    def getAlignRefSeqRange(self, authSeqId="", partId=0, refSeqs=None):
        """ """
        if refSeqs is None:
            refSeqs = []
        if len(refSeqs) == 0:
            return None, None
        #
        authSeqs = self.__getSequenceFromDataStore(authSeqId)
        if not authSeqs:
            return None, None
        #
        self.__checkPartInfoDict(authSeqId)
        if partId not in self.__partInfoDict:
            return None, None
        #
        alignIndexList, _alignLength, _numMatch, _numMatchGaps = self.__getAuthRefAlignIndexList(
            authSeqs, refSeqs, self.__partInfoDict[partId][0] - 1, self.__partInfoDict[partId][1] - 1
        )
        #
        alignStart = -1
        alignEnd = -1
        for rowIdx, alignIdx in enumerate(alignIndexList):
            if alignIdx[0] < 0:
                continue
            #
            if alignStart < 0:
                alignStart = rowIdx
            #
            alignEnd = rowIdx
        #
        if (alignStart >= 0) and (alignEnd >= 0):
            refSeqBeg = None
            refSeqEnd = None
            for idx in range(alignStart, alignEnd + 1):
                if alignIndexList[idx][1] < 0:
                    continue
                #
                if refSeqBeg is None:
                    refSeqBeg = alignIndexList[idx][1] + 1
                #
                refSeqEnd = alignIndexList[idx][1] + 1
            #
            return refSeqBeg, refSeqEnd
        #
        return None, None

    def addRefAlignIndices(self, authSeqId="", refSeqId="", partId=0):
        """Add new ref sequence alignmet index"""
        if (authSeqId == "") or (refSeqId == "") or (partId == 0):
            return
        #
        indexList = self._getRefAlignIndexList([authSeqId, refSeqId])
        if indexList:
            return
        #
        authSeqs = self.__getSequenceFromDataStore(authSeqId)
        if not authSeqs:
            return
        #
        refSeqs = self.__getSequenceFromDataStore(refSeqId)
        if not refSeqs:
            return
        #
        self.__checkPartInfoDict(authSeqId)
        if partId not in self.__partInfoDict:
            return
        #
        alignIndexList, alignLength, numMatch, numMatchGaps = self.__getAuthRefAlignIndexList(
            authSeqs, refSeqs, self.__partInfoDict[partId][0] - 1, self.__partInfoDict[partId][1] - 1
        )
        if alignIndexList:
            self._insertRefAlignIndexList([authSeqId, refSeqId], alignIndexList)
            #
            (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = self._getUnpackSeqLabel(refSeqId)
            sFeature = self._getFeatureObjByUnpackLabelFromDataStore(seqType, seqInstId, seqPartId, seqAltId, seqVersion)
            sFeature.setAuthRefAlignDetails(
                seqLen=len(refSeqId), alignLen=int(alignLength), seqSim=float(numMatch) / float(alignLength), seqSimWithGaps=float(numMatchGaps) / float(alignLength)
            )
            self.setFeature(sFeature.get(), seqType, seqInstId, seqPartId, seqAltId, seqVersion)
        #
        self.serialize()

    def checkResidueNameReplaceCompatibility(self, alignIndex, origResName, resLabelIndex, newResName):
        """Check if residue name replacement is allowed based on the aligned coordiante residue(s)"""
        if (alignIndex < 0) or (alignIndex >= len(self._seqAlignList)):
            return "Compatibility checking for residue name replacement failed.<br />\n"
        #
        error = self.__isAuthResidueNameCompatible(alignIndex, newResName)
        errorMsg = ""
        if error:
            errorMsg = "'" + origResName + " " + resLabelIndex + "' can not be replaced by '" + newResName + "':<br /><br />\n\n" + error
        #
        return errorMsg

    def checkResidueMovingCompatibility(self, srcResidueId, dstResidueId):
        """Check if a residue moving is allowed based on the aligned coordiante residue(s)"""
        srcResidueLabel = copy.deepcopy(self._getResLabelFromResLabelId(srcResidueId))
        dstResidueLabel = copy.deepcopy(self._getResLabelFromResLabelId(dstResidueId))
        #
        srcAlignPos = srcResidueLabel.getAlignmentIndex()
        dstAlignPos = dstResidueLabel.getAlignmentIndex()
        if (srcAlignPos < 0) or (srcAlignPos >= len(self._seqAlignList)) or (dstAlignPos < 0) or (dstAlignPos >= len(self._seqAlignList)):
            return "Compatibility checking for residue moving failed.<br />\n"
        #
        seqType = srcResidueLabel.getSequenceType()
        instId = srcResidueLabel.getSequenceInstId()
        ResName = srcResidueLabel.getResidueCode3()
        ResNum = srcResidueLabel.getResidueLabelIndex()
        #
        error = ""
        if seqType == "auth":
            error = self.__isAuthResidueNameCompatible(dstAlignPos, ResName)
        elif seqType == "xyz":
            error = self.__isXyzResidueNameCompatible(instId, srcAlignPos, dstAlignPos)
        #
        errorMsg = ""
        if error:
            errorMsg = "'" + ResName + " " + ResNum + "' can not be moved to position '" + str(dstAlignPos) + "':<br /><br />\n\n" + error
        #
        return errorMsg

    def getGlobalEditList(self):
        """Generate editlist for global input form"""
        editList = []
        unCompatibleList = []
        chainIdList = str(self._reqObj.getValue("chainids")).split(",")
        for chainId in chainIdList:
            chainIdx = -1
            chainSeqId = ""
            for seqId, idx in self._seqAlignLabelIndices.items():
                tL = str(seqId).strip().split("_")
                if len(tL) < 3:
                    continue
                #
                if (tL[0] == "xyz") and (tL[1] == chainId):
                    chainIdx = idx
                    chainSeqId = str(seqId).strip()
                    break
                #
            #
            if chainIdx < 0:
                continue
            #
            start_position = str(self._reqObj.getValue("start_position_" + chainId))
            if not start_position:
                continue
            #
            end_position = str(self._reqObj.getValue("end_position_" + chainId))
            if not end_position:
                continue
            #
            move_to = str(self._reqObj.getValue("move_to_" + chainId))
            if not move_to:
                continue
            #
            start = int(start_position) - 1
            end = int(end_position) - 1
            target = int(move_to) - 1
            pairList, gapList = self.__getPairAndGapList(start, end, target)
            edPairList, unPairList = self.__getPairEditList(self._seqAlignLabelIndices[self._authLabel], chainIdx, chainSeqId, pairList)
            if len(edPairList) > 0:
                editList.extend(edPairList)
            #
            if len(unPairList) > 0:
                unCompatibleList.extend(unPairList)
                continue
            #
            edGapList = self.__getGapEditList(chainIdx, chainSeqId, gapList)
            if len(edGapList) > 0:
                editList.extend(edGapList)
            #
        #
        return editList, unCompatibleList

    def getGlobalEditAuthSeqList(self):
        """Generate editlist for insertion missing author sequence part(s) based on reference sequence."""
        authIdx = -1
        authSeqId = ""
        refIdx = -1
        for seqId, idx in self._seqAlignLabelIndices.items():
            tL = str(seqId).strip().split("_")
            if len(tL) < 3:
                continue
            #
            if tL[0] == "auth":
                authIdx = idx
                authSeqId = str(seqId).strip()
            elif tL[0] == "ref":
                refIdx = idx
            #
        #
        if (authIdx < 0) or (refIdx < 0):
            return "Can not find correct sequence alignment.", []
        #
        blockList = []
        selected_blockid = str(self._reqObj.getValue("blockid"))
        if selected_blockid == "all":
            repblocknum = int(self._reqObj.getValue("repblocknum"))
            for i in range(0, repblocknum):
                blockList.append(i)
            #
        else:
            blockList.append(int(selected_blockid))
        #
        (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = self._getUnpackSeqLabel(authSeqId)
        featureD = self.getFeature(seqType, seqInstId, seqPartId, seqAltId, seqVersion)
        #
        polymerTypeCode = "AA"
        if "POLYMER_TYPE" in featureD:
            polymerTypeCode = featureD["POLYMER_TYPE"]
        #
        errMsg = ""
        editList = []
        for blockid in blockList:
            start_position = str(self._reqObj.getValue("repstartposition_" + str(blockid)))
            if not start_position:
                continue
            #
            end_position = str(self._reqObj.getValue("rependposition_" + str(blockid)))
            if not end_position:
                continue
            #
            for i in range(int(start_position) - 1, int(end_position)):
                if (i < 0) or (i >= len(self._seqAlignList)):
                    continue
                #
                authRecord = self._seqAlignList[i][authIdx]
                refRecord = self._seqAlignList[i][refIdx]
                if (authRecord[1] != "") and (authRecord[1] != "."):
                    if errMsg:
                        errMsg += "\n"
                    errMsg += "Author sequence already has residue '" + authRecord[1] + "' at position '" + str(i + 1) + "'."
                    continue
                elif (refRecord[1] == "") or (refRecord[1] == "."):
                    if errMsg:
                        errMsg += "\n"
                    errMsg += "There is no residue in reference sequence at position '" + str(i + 1) + "'."
                    continue
                #
                currId = self._getResLabelId(
                    seqType=seqType,
                    seqInstId=seqInstId,
                    seqAltId=seqAltId,
                    seqVersion=seqVersion,
                    residueCode3=authRecord[1],
                    residueLabelIndex=authRecord[2],
                    alignIndex=i,
                    seqIndex=codeSeqIndex(authRecord[3]),
                    residueType=polymerTypeCode,
                    seqPartId=seqPartId,
                )
                newId = self._getResLabelId(
                    seqType=seqType,
                    seqInstId=seqInstId,
                    seqAltId=seqAltId,
                    seqVersion=seqVersion,
                    residueCode3=refRecord[1],
                    residueLabelIndex=authRecord[2],
                    alignIndex=i,
                    seqIndex=codeSeqIndex(authRecord[3]),
                    residueType=polymerTypeCode,
                    seqPartId=seqPartId,
                )
                editList.append((currId, refRecord[0], refRecord[1], authRecord[0], newId))
            #
        #
        if not editList:
            if errMsg:
                return errMsg, []
            return "Can not find missing residue(s) in author sequence.", []
        #
        return errMsg, editList

    def _checkAndUpdateAlignment(self, inputIdList, selectedIdList):
        """Check if input alignIds match the default selected alignIds saved in alignment pickle file and
        update alignemnt if necessary
        """
        retInputIdList, retSelectedIdList, allAlignIdList, self.__xyzLabel, refIdList, self._selfRefPartIdList, self.__extraAuthLabel = self._getProperAlignIdList(
            inputIdList, selectedIdList
        )
        #
        if self.__checkCurrentAlignIdMatch(allAlignIdList, self._seqAlignLabelIndices):
            return retInputIdList, retSelectedIdList
        #
        self.__local_authLabel = retSelectedIdList[0]
        #
        if not self.__checkCurrentAlignIdMatch(self.__xyzLabel, self._xyzAlignLabelIndices):
            self.__getXyzAlignIndexList()
            if not self.__getCurrentAuthCoordAlignment():
                self.__alignFlag = False
                return retInputIdList, retSelectedIdList
            #
        #
        self._copyXyzAlignmentToSeqAlignment()
        #
        if (len(refIdList) > 0) and self.__getRefLabelDictFromRefIdList(refIdList):
            self.__checkPartInfoDict(self.__local_authLabel)
            self.__getCurrentAuthRefAlignment()
        #
        if self.__alignFlag and (len(self.__extraAuthLabel) > 0):
            self.__alignExtraAuthSequence()
        #
        if self.__alignFlag and self.__newXyzAlignIndexListFlag and (len(self.__xyzLabel) > 0) and (len(self.__xyzAlignIndexList) > 0):
            self._insertXyzAlignIndexList(self.__xyzLabel, self.__xyzAlignIndexList)
        #
        self._updateAlignmentDataStore(self.__local_authLabel)
        return retInputIdList, retSelectedIdList

    def _isCompatible(self, comment, newResName):
        """Verify the coordinate residue side-chain pattern defined in comment is compatible with new residue name"""
        compatibilityMap = {
            "ala-gly-like": [
                "ALA",
                "ARG",
                "ASN",
                "ASP",
                "CYS",
                "GLN",
                "GLU",
                "GLY",
                "HIS",
                "ILE",
                "LEU",
                "LYS",
                "MET",
                "MSE",
                "PHE",
                "PRO",
                "PYL",
                "SEC",
                "SER",
                "THR",
                "TRP",
                "TYR",
                "UNK",
                "VAL",
            ],
            "c-gamma-like": ["ARG", "ASN", "ASP", "GLN", "GLU", "HIS", "ILE", "LEU", "LYS", "MET", "MSE", "PHE", "PYL", "TRP", "THR", "TYR", "VAL"],
            "c-delta-like": ["ARG", "GLN", "GLU", "HIS", "ILE", "LYS", "PHE", "PYL", "TRP", "TYR"],
            "ile-val-like": ["ILE", "VAL"],
            "ser-thr-like": ["SER", "THR"],
            "phe-tyr-like": ["PHE", "TYR"],
            "lys-pyl-like": ["LYS", "PYL"],
        }

        #
        if not comment:
            return False
        #
        for patternType, allowList in compatibilityMap.items():
            if (comment.find(patternType) != -1) and (newResName in allowList):
                return True
            #
        #
        return False

    def _getNaturalSourceWarningMessage(self, sourceType, eelCommentL):
        """Create warning for natural source entity with [ "engineered mutation", "expression tag", "linker" ] sequence annotation"""
        if (sourceType.strip().upper() != "NAT") or (len(eelCommentL) < 1):
            return ""
        #
        warningMsg = "Entity '" + self._entityId + "' with natural source has '" + "', '".join(sorted(set(eelCommentL))) + "' sequence annotation information.<br />"
        return warningMsg

    def _clearAllConflicts(self, authIdx):
        """Clear all numeric conflicts and conflicting comments"""
        selfRefPosPartIdDict = {}
        if len(self._selfRefPartIdList) > 0:
            for selfRefId in self._selfRefPartIdList:
                tup = selfRefId.strip().split("_")
                if int(tup[2]) in self._partPosDict:
                    for pos in range(self._partPosDict[int(tup[2])][0], self._partPosDict[int(tup[2])][1] + 1):
                        selfRefPosPartIdDict[pos] = int(tup[2])
                    #
                #
            #
        #
        for rowIdx, alignTuple in enumerate(self._seqAlignList):
            for columnIdx, alignRecord in enumerate(alignTuple):
                alignRecord[6] = 0
                alignRecord[7] = ""
                if columnIdx == authIdx:
                    if len(alignRecord[5]) < 1:
                        continue
                    #
                    # Remove default conflicting comments or all conflicting comments for self reference part
                    #
                    if (not alignRecord[5].startswith("ANNOTATED:")) or (rowIdx in selfRefPosPartIdDict):
                        alignRecord[5] = ""
                    #
                #
            #
        #

    def _assignAllConflicts(self, authLabel, selectedIdList, writeConflictFlag=True):
        """Assign numeric conflicts and conflicting comments"""
        self._missingSeqMap = {}
        #
        authIdx = self._seqAlignLabelIndices[authLabel]
        #
        self.__conflictMap = {}  # pylint: disable=attribute-defined-outside-init
        totalSeqCoodConflict = 0
        for otherIdx in sorted(self._reverseSeqAlignLabelIndices.keys()):
            otherLabel = self._reverseSeqAlignLabelIndices[otherIdx]
            if (not otherLabel.startswith("xyz")) or (otherLabel not in selectedIdList):
                continue
            #
            numConflict, missingAuthSeqList = self.__assignNumericConflicts(authIdx, otherIdx, "xyz", 0, len(self._seqAlignList) - 1)
            totalSeqCoodConflict += numConflict
            if missingAuthSeqList:
                (_seqType, seqInstId, seqPartId, _seqAltId, _seqVersion) = self._getUnpackSeqLabel(otherLabel)
                self._missingSeqMap[seqInstId] = missingAuthSeqList
            #
        #
        self._checkPartStartEndPosMap(authIdx, authLabel)
        #
        for otherIdx in sorted(self._reverseSeqAlignLabelIndices.keys()):
            otherLabel = self._reverseSeqAlignLabelIndices[otherIdx]
            if (not otherLabel.startswith("ref")) or (otherLabel not in selectedIdList):
                continue
            #
            (_seqType, seqInstId, seqPartId, _seqAltId, _seqVersion) = self._getUnpackSeqLabel(otherLabel)
            #
            start_position = 0
            end_position = len(self._seqAlignList) - 1
            if int(seqPartId) in self._partPosDict:
                start_position = self._partPosDict[int(seqPartId)][0]
                end_position = self._partPosDict[int(seqPartId)][1]
            #
            self.__getRefSeqVariants(otherLabel)
            numConflict, missingAuthSeqList = self.__assignNumericConflicts(authIdx, otherIdx, "ref", start_position, end_position)
            self.__annotateConflictingComments(authIdx, otherIdx, start_position, end_position, writeConflictFlag)
        #
        if writeConflictFlag:
            if totalSeqCoodConflict > 0:
                self.__conflictMap["mismatch"] = totalSeqCoodConflict
            #
            self._missingSeqMap["summary_page"] = self.__conflictMap
        #
        return totalSeqCoodConflict

    def _checkPartStartEndPosMap(self, authIdx, authLabel):
        """Get auth sequence part start & end pos range in alignment if it does not exist"""
        if self._partPosDict:
            return
        #
        self.__checkPartInfoDict(authLabel)
        #
        if not self.__partInfoDict:
            self._partPosDict = {1: [0, len(self._seqAlignList) - 1]}
            return
        #
        SeqIdxPartIdMap = {}
        for partId, rangePair in self.__partInfoDict.items():
            for i in range(rangePair[0] - 1, rangePair[1]):
                SeqIdxPartIdMap[i] = int(partId)
            #
        #
        partId = 1
        self._partPosDict = {}
        self._partPosDict[partId] = [0, 0]
        #
        for idx, alignTup in enumerate(self._seqAlignList):
            seqIdx, _comment = decodeIndex(alignTup[authIdx][3])
            if seqIdx not in SeqIdxPartIdMap:
                continue
            #
            oldPartId = SeqIdxPartIdMap[seqIdx]
            if oldPartId != partId:
                self._partPosDict[partId][1] = idx - 1
                partId = oldPartId
                self._partPosDict[partId] = [idx, idx]
            #
        #
        self._partPosDict[partId][1] = len(self._seqAlignList) - 1

    def _updateAlignmentDataStore(self, newAuthLabel):
        """Update persist alignment storage"""
        if self.__alignFlag:
            self._authLabel = newAuthLabel
            self._hasAlignmentFlag = True
        else:
            self._resetAlignmentInfo()
            self._addErrorMessage("Failed to generate alignment for entity '" + self._entityId + "'.")
        #
        self.serialize()

    def _decodeComment(self, comment):
        """Remove prefix "DEFAULT:" or ""ANNOTATED:" from comment string"""
        if str(comment).strip().startswith("ANNOTATED:"):
            return "ANNOTATED", str(comment).strip()[10:].lower()
        elif str(comment).strip().startswith("SELECTED:"):
            return "SELECTED", str(comment).strip()[9:].lower()
        elif str(comment).strip().startswith("DEFAULT:"):
            return "DEFAULT", str(comment).strip()[8:].lower()
        else:
            return "", str(comment).strip().lower()
        #

    def __checkCurrentAlignIdMatch(self, alignIdList, alignIndices):
        """ """
        if len(alignIdList) != len(alignIndices):
            return False
        #
        for seqId in alignIdList:
            if seqId not in alignIndices:
                return False
            #
        #
        return True

    def __doAllAlignments(self):
        """Generate pseduo multiple alignment among auth, coordiates and reference sequences"""
        self._resetAlignmentInfo()
        self._resetLogInfo()
        #
        if not self.__alignFlag:
            return
        #
        authSeqs = self.__getSequenceFromDataStore(self.__local_authLabel)
        if (not authSeqs) or (not self.__getCurrentAuthCoordAlignment(existAuthSeqs=authSeqs)):
            self.__alignFlag = False
            return
        #
        self._copyXyzAlignmentToSeqAlignment()
        #
        if not self.__getCurrentAuthRefAlignment(authSeqs):
            self.__alignFlag = False
            return
        #
        self.__alignExtraAuthSequence(existAuthSeqs=authSeqs)

    def __resetAllAlignments(self):
        """Generate pseduo multiple alignment among auth, coordiates and reference sequences"""
        authSeqs = self.__getSequenceFromDataStore(self.__local_authLabel)
        if not authSeqs:
            self.__alignFlag = False
            return
        #
        #       foundFalg = True
        #       for xyzLabel in  self.__xyzLabel:
        #           if xyzLabel not in self._xyzAlignLabelIndices:
        #               foundFalg = False
        #               break
        #           #
        #       #
        foundFalg = False
        if not foundFalg:
            self._resetAlignmentInfo()
            self._resetLogInfo()
            self.__getXyzAlignIndexList()
            if not self.__alignFlag:
                return
            #
            if not self.__getCurrentAuthCoordAlignment(existAuthSeqs=authSeqs):
                self.__alignFlag = False
                return
            #
        #
        self._copyXyzAlignmentToSeqAlignment()
        #
        if not self.__getCurrentAuthRefAlignment(authSeqs):
            self.__alignFlag = False
        #

    def __getSequenceFromDataStore(self, idLabel):
        """Get sequence from sequence data storage using packed id label"""
        seq = self._getSequenceByPackLabelFromDataStore(idLabel)
        if self._verbose:
            self._lfh.write("Instance '%r' sequence length=%d\n" % (idLabel, len(seq)))
        #
        if (not seq) and idLabel.startswith("auth_"):
            self.__alignFlag = False
        #
        return seq

    def __getCurrentAuthCoordAlignment(self, existAuthSeqs=None):
        """Generate auth/coordiates sequence alignment"""
        if existAuthSeqs is None:
            existAuthSeqs = []

        authSeqs = existAuthSeqs
        if not authSeqs:
            authSeqs = self.__getSequenceFromDataStore(self.__xyzLabel[0])
        #
        if not self.__alignFlag:
            return False
        #
        self._xyzAlignLabelIndices = {}
        self._reverseXyzAlignLabelIndices = {}
        seqList = []
        self._xyzAlignLabelIndices[self.__xyzLabel[0]] = len(seqList)
        self._reverseXyzAlignLabelIndices[len(seqList)] = self.__xyzLabel[0]
        seqList.append(authSeqs)
        #
        for label in self.__xyzLabel[1:]:
            seqs = self.__getSequenceFromDataStore(label)
            self._xyzAlignLabelIndices[label] = len(seqList)
            self._reverseXyzAlignLabelIndices[len(seqList)] = label
            seqList.append(seqs)
        #
        if len(seqList) < 2:
            if self._verbose:
                self._lfh.write("AlignmentTools.__getCurrentAuthCoordAlignment seqList=%d\n" % len(seqList))
            #
            return False
        #
        errMsg, self._xyzAlignList = getSeqAlignment(seqList, self.__xyzAlignIndexList, self._reverseXyzAlignLabelIndices, "xyz", self._gapSymbol)
        if errMsg:
            if self._verbose:
                self._lfh.write("AlignmentTools.__getCurrentAuthCoordAlignment errMsg=%r\n" % errMsg)
            #
            self._addErrorMessage(errMsg)
            return False
        #
        self.__getAuthDefinedMutations()
        return True

    def __getCurrentAuthRefAlignment(self, existAuthSeqs=None):
        """Generate auth/ref sequence alignment(s)"""
        if existAuthSeqs is None:
            existAuthSeqs = []

        if not self.__refLabelDict:
            return True
        #
        authSeqs = existAuthSeqs
        if not authSeqs:
            authSeqs = self.__getSequenceFromDataStore(self.__local_authLabel)
        #
        if not self.__alignFlag:
            return False
        #
        multiAlignList = []
        multiAlignList.append(self._seqAlignList)
        #
        pA = PseudoMultiAlign()
        pA.setAuthScore()
        pA.setAuthSequence(self.__generateInputSeqInfoForPseudoMultiAlignFromAlignList(self._seqAlignList))
        #
        for partId in sorted(self.__refLabelDict.keys()):
            for label in self.__refLabelDict[partId]:
                seqs = self.__getSequenceFromDataStore(label)
                if not seqs:
                    return False
                #
                if label not in self.__refAlignIndexDict:
                    if self._verbose:
                        self._lfh.write("No reference sequence align index for %s\n" % label)
                    #
                    self._addErrorMessage("No reference sequence align index for '" + label + "'.")
                    return False
                #
                reverseIndices = {}
                reverseIndices[0] = self.__local_authLabel
                reverseIndices[1] = label
                seqList = [authSeqs, seqs]
                #
                errMsg, seqAlignList = getSeqAlignment(seqList, self.__refAlignIndexDict[label], reverseIndices, "ref", self._gapSymbol)
                if errMsg:
                    if self._verbose:
                        self._lfh.write("AlignmentTools.__getCurrentAuthRefAlignment errMsg=%r\n" % errMsg)
                    #
                    self._addErrorMessage(errMsg)
                    return False
                #
                multiAlignList.append(seqAlignList)
                pA.addAlignSequence(self.__generateInputSeqInfoForPseudoMultiAlignFromAlignList(seqAlignList))
                #
                idxPos = len(self._seqAlignLabelIndices)
                self._seqAlignLabelIndices[label] = idxPos
                self._reverseSeqAlignLabelIndices[idxPos] = label
            #
        #
        errMsg, self._seqAlignList = mergeSeqAlignment(multiAlignList, pA.getAlignIndices(), self._gapSymbol)
        if errMsg:
            self._addErrorMessage(errMsg)
            return False
        #
        return True

    def __alignExtraAuthSequence(self, existAuthSeqs=None):
        """ """
        if existAuthSeqs is None:
            existAuthSeqs = []

        if not self.__extraAuthLabel:
            return
        #
        authSeqs = existAuthSeqs
        if not authSeqs:
            authSeqs = self.__getSequenceFromDataStore(self.__local_authLabel)
        #
        if not authSeqs:
            return
        #
        seqList = []
        reverseIndices = {}
        reverseIndices[len(seqList)] = self.__local_authLabel
        seqList.append(authSeqs)
        #
        pA = PseudoMultiAlign()
        pA.setAuthScore()
        pA.setAuthSequence(self.__generateInputSeqInfoForPseudoMultiAlignFromSeqList(inputSeqList=authSeqs))
        #
        for label in self.__extraAuthLabel:
            seqs = self.__getSequenceFromDataStore(label)
            if not seqs:
                return
            #
            reverseIndices[len(seqList)] = label
            seqList.append(seqs)
            pA.addAlignSequence(self.__generateInputSeqInfoForPseudoMultiAlignFromSeqList(inputSeqList=seqs))
            #
            idxPos = len(self._seqAlignLabelIndices)
            self._seqAlignLabelIndices[label] = idxPos
            self._reverseSeqAlignLabelIndices[idxPos] = label
        #
        indexList = pA.getAlignIndices()
        errMsg, seqAlignList = getSeqAlignment(seqList, indexList, reverseIndices, "auth", self._gapSymbol)
        if errMsg:
            if self._verbose:
                self._lfh.write("AlignmentTools.__alignExtraAuthSequence errMsg=%r\n" % errMsg)
            #
            self._addErrorMessage(errMsg)
            self.__alignFlag = False
            return
        #
        multiAlignList = []
        multiAlignList.append(self._seqAlignList)
        #
        pA.clear()
        pA.setAuthScore()
        pA.setAuthSequence(self.__generateInputSeqInfoForPseudoMultiAlignFromAlignList(self._seqAlignList))
        #
        multiAlignList.append(seqAlignList)
        pA.addAlignSequence(self.__generateInputSeqInfoForPseudoMultiAlignFromAlignList(seqAlignList))
        #
        errMsg, self._seqAlignList = mergeSeqAlignment(multiAlignList, pA.getAlignIndices(), self._gapSymbol)
        if errMsg:
            if self._verbose:
                self._lfh.write("AlignmentTools.__alignExtraAuthSequence errMsg=%r\n" % errMsg)
            #
            self._addErrorMessage(errMsg)
            self.__alignFlag = False
        #

    def __getXyzAlignIndexList(self):
        """Re-generate self.__xyzAlignIndexList based on selected self.__xyzLabel list"""
        if len(self.__xyzLabel) < 2:
            self.__alignFlag = False
            return
        #
        self.__xyzAlignIndexList = self._getXyzAlignIndexList(self.__xyzLabel)
        if self.__xyzAlignIndexList:
            if self._verbose:
                self._lfh.write("AlignmentTools.__getXyzAlignIndexList found existing alignment indices for id list=%r\n" % self.__xyzLabel)
            #
            return
        #
        authSeqs = self.__getSequenceFromDataStore(self.__xyzLabel[0])
        if not authSeqs:
            self.__alignFlag = False
            return
        #
        pA = PseudoMultiAlign()
        pA.setAuthSequence(self.__generateInputSeqInfoForPseudoMultiAlignFromSeqList(inputSeqList=authSeqs))
        #
        for label in self.__xyzLabel[1:]:
            seqs = self.__getSequenceFromDataStore(label)
            if not seqs:
                self.__alignFlag = False
                return
            #
            pA.addAlignSequence(self.__generateInputSeqInfoForPseudoMultiAlignFromSeqList(seqType="xyz", inputSeqList=seqs))
        #
        self.__xyzAlignIndexList = pA.getAlignIndices()
        self.__newXyzAlignIndexListFlag = True

    def __getAuthRefAlignIndexList(self, authSeqs, refSeqs, startIdx, endIdx):
        """ """
        alignIndexList = []
        for i in range(0, startIdx):
            alignIndexList.append([i, -1])
        #
        pA = PseudoMultiAlign()
        pA.setRefScore()
        pA.setAuthSequence(self.__generateInputSeqInfoForPseudoMultiAlignFromSeqList(inputSeqList=authSeqs, start=startIdx, end=endIdx))
        pA.addAlignSequence(self.__generateInputSeqInfoForPseudoMultiAlignFromSeqList(inputSeqList=refSeqs))
        indexList = pA.getAlignIndices()
        #
        alignLength = 0
        numMatch = 0
        numMatchGaps = 0
        for indexTuple in indexList:
            alignIndexList.append(indexTuple)
            if (indexTuple[0] < 0) and (indexTuple[1] < 0):
                continue
            #
            alignLength += 1
            if indexTuple[0] >= 0:
                if indexTuple[1] < 0:
                    numMatchGaps += 1
                elif str(authSeqs[indexTuple[0]][0]).strip().upper() == str(refSeqs[indexTuple[1]][0]).strip().upper():
                    numMatch += 1
                    numMatchGaps += 1
                #
            #
        #
        for i in range(endIdx + 1, len(authSeqs)):
            alignIndexList.append([i, -1])
        #
        return alignIndexList, alignLength, numMatch, numMatchGaps

    def __getRefLabelDictFromRefIdList(self, refIdList):
        """Get self.__refLabelDict from input refIdList"""
        self.__refLabelDict = {}
        self.__refAlignIndexDict = {}
        #
        if (not refIdList) or (not self.__alignFlag):
            return False
        #
        for refId in refIdList:
            indexList = self._getRefAlignIndexList([self.__local_authLabel, refId])
            #
            (_seqType, _seqInstId, seqPartId, _seqAltId, _seqVersion) = self._getUnpackSeqLabel(refId)
            #
            if len(indexList) == 0:
                authSeqs = self.__getSequenceFromDataStore(self.__local_authLabel)
                refSeqs = self.__getSequenceFromDataStore(refId)
                self.__checkPartInfoDict(self.__local_authLabel)
                if authSeqs and refSeqs and (seqPartId in self.__partInfoDict):
                    indexList, _alignLength, _numMatch, _numMatchGaps = self.__getAuthRefAlignIndexList(
                        authSeqs, refSeqs, self.__partInfoDict[seqPartId][0] - 1, self.__partInfoDict[seqPartId][1] - 1
                    )
                    if len(indexList) > 0:
                        self._insertRefAlignIndexList([self.__local_authLabel, refId], indexList)
                    #
                #
            #
            if len(indexList) == 0:
                self.__refLabelDict = {}
                self.__refAlignIndexDict = {}
                self.__alignFlag = False
                return False
            #
            if seqPartId in self.__refLabelDict:
                self.__refLabelDict[seqPartId].append(refId)
            else:
                self.__refLabelDict[seqPartId] = [refId]
            #
            self.__refAlignIndexDict[refId] = indexList
        #
        return True

    def __getAuthDefinedMutations(self):
        """Extract all author input mutation information (e.g. V178A ) from _entity.pdbx_mutation"""
        featureD = self._getFeatureByPackLabelFromDataStore(self.__local_authLabel)
        if "ENTITY_MUTATION_DETAILS" not in featureD:
            return
        #
        self._authDefinedMutationList = []
        mutMap = {}
        #
        mutation_details = featureD["ENTITY_MUTATION_DETAILS"].strip().upper().replace(",", " ")
        for val in mutation_details.split(" "):
            if len(val) < 3:
                continue
            #
            hasDigit = False
            allDigit = True
            for i in range(1, len(val) - 1):
                if val[i].isdigit():
                    hasDigit = True
                else:
                    allDigit = False
                #
            #
            if hasDigit and allDigit and val[0].isalpha() and val[len(val) - 1].isalpha():
                mutMap[val[1:]] = val[:1]
                self._authDefinedMutationList.append(val)
            #
            if mutMap:
                for alignTup in self._xyzAlignList:
                    for i in range(1, len(alignTup)):
                        if (alignTup[0][0] == alignTup[i][0]) and (alignTup[0][2] != alignTup[i][2]) and ((alignTup[i][2] + alignTup[i][0]) in mutMap):
                            mut = mutMap[alignTup[i][2] + alignTup[i][0]] + alignTup[0][2] + alignTup[0][0]
                            if mut not in self._authDefinedMutationList:
                                self._authDefinedMutationList.append(mut)
                            #
                        #
                    #
                #
            #
        #

    def __getRefSeqVariants(self, refLabel):
        """ """
        self.__refSeqVariantList = []
        #
        featureD = self._getFeatureByPackLabelFromDataStore(refLabel)
        if ("REF_SEQ_VARIANT" not in featureD) or (not featureD["REF_SEQ_VARIANT"]):
            return
        #
        self.__refSeqVariantList = featureD["REF_SEQ_VARIANT"].split(",")

    def __generateInputSeqInfoForPseudoMultiAlignFromSeqList(self, seqType="auth", inputSeqList=None, start=0, end=0):
        """Generate input sequence information for PseudoMultiAlign program"""
        if inputSeqList is None:
            inputSeqList = []

        if end == 0:
            end = len(inputSeqList) - 1
        #
        returnSeqList = []
        for (i, seqTup) in enumerate(inputSeqList):
            if (i < start) or (i > end):
                continue
            #
            if seqType == "xyz":
                # added 'auth_numbering' as 4th value in seqTuple. See setCoordinateInstanceInfo() method in
                # util/UpdateSequenceDataStoreUtils.py for the definition of all values in coordinate seqTup
                returnSeqList.append((str(seqTup[0]).upper(), str(i), str(seqTup[6]), str(seqTup[1])))
            else:
                returnSeqList.append((str(seqTup[0]).upper(), str(i)))
            #
        #
        return returnSeqList

    def __generateInputSeqInfoForPseudoMultiAlignFromAlignList(self, alignList):
        """Generate input sequence information for PseudoMultiAlign program"""
        seqList = []
        for (i, alignTup) in enumerate(alignList):
            seqList.append((str(alignTup[0][1]).upper(), str(i)))
        #
        return seqList

    def __checkPartInfoDict(self, authLabel):
        """Get Parts information if it does not exist"""
        if self.__partInfoDict:
            return
        #
        self.__partInfoDict = self._getPartInfoDict(authLabel)
        #
        if self.__partInfoDict:
            return
        #
        (seqType, seqInstId, _seqPartId, seqAltId, seqVersion) = self._getUnpackSeqLabel(authLabel)
        parIdList = self.getPartIdList(seqType, seqInstId)
        for partId in parIdList:
            seqFeatureD = self.getFeature(seqType, seqInstId, partId, seqAltId, seqVersion)
            self.__partInfoDict[int(seqFeatureD["AUTH_SEQ_PART_ID"])] = (int(seqFeatureD["AUTH_SEQ_NUM_BEGIN"]), int(seqFeatureD["AUTH_SEQ_NUM_END"]))
        #
        self._insertPartInfoDict(authLabel, self.__partInfoDict)

    def __assignNumericConflicts(self, authIdx, otherIdx, otherType, start_position, end_position):
        """Assign numeric conflicts"""
        numConflict = 0
        missingAuthSeqList = []
        for intIdx, alignTup in enumerate(self._seqAlignList):
            if (intIdx < start_position) or (intIdx > end_position):
                continue
            #
            if (alignTup[authIdx][1] == self._gapSymbol) and (alignTup[otherIdx][1] == self._gapSymbol):
                continue
            #
            isConflict, authConflict, otherConflict = assignConflict("auth", alignTup[authIdx][1], otherType, alignTup[otherIdx][1], self._gapSymbol)
            if isConflict:
                if alignTup[otherIdx][1] != self._gapSymbol:
                    numConflict += 1
                #
                if otherType == "xyz" and (alignTup[authIdx][1] == self._gapSymbol):
                    missingAuthSeqList.append((alignTup[otherIdx][1], alignTup[otherIdx][2]))
                #
                alignTup[authIdx][6] = authConflict[0]
                alignTup[authIdx][7] = authConflict[1]
                alignTup[otherIdx][6] = otherConflict[0]
                alignTup[otherIdx][7] = otherConflict[1]
            #
        #
        return numConflict, missingAuthSeqList

    def __annotateConflictingComments(self, authIdx, refIdx, beg_pos, end_pos, writeConflictFlag):
        """Add default 'expression tag annotation comments' for leading or trailing gap residues, and
        comments for any conflicts which can be attributed to a residue modification.
        """
        first_fragment = False
        if beg_pos < 10:
            first_fragment = True
        #
        last_fragment = False
        if (len(self._seqAlignList) - end_pos) < 10:
            last_fragment = True
        #
        first_pair_position = -1
        for i in range(beg_pos, end_pos):
            if (self._seqAlignList[i][authIdx][1] != self._gapSymbol) and (self._seqAlignList[i][refIdx][1] != self._gapSymbol):
                first_pair_position = i
                break
            #
        #
        last_pair_position = beg_pos + 1
        for i in range(end_pos, beg_pos, -1):
            if (self._seqAlignList[i][authIdx][1] != self._gapSymbol) and (self._seqAlignList[i][refIdx][1] != self._gapSymbol):
                last_pair_position = i
                break
            #
        #
        for i in range(beg_pos, end_pos + 1):
            if self._seqAlignList[i][authIdx][1] == self._seqAlignList[i][refIdx][1]:
                self._seqAlignList[i][authIdx][5] = ""
                self._seqAlignList[i][refIdx][5] = ""
                continue
            #
            comment = ""
            if self._seqAlignList[i][refIdx][1] == self._gapSymbol:
                if (i == 0) and (self._seqAlignList[i][authIdx][1] in ("MET", "MSE")):
                    comment = "initiating methionine"
                elif i < first_pair_position:
                    if first_fragment:
                        comment = "expression tag"
                    else:
                        comment = "linker"
                    #
                elif i > last_pair_position:
                    if last_fragment:
                        comment = "expression tag"
                    else:
                        comment = "linker"
                    #
                else:
                    comment = "insertion"
                #
            elif self._seqAlignList[i][authIdx][1] == self._gapSymbol:
                comment = "deletion"
            else:
                comment = self.__annotateRealConflict(self._seqAlignList[i][authIdx], self._seqAlignList[i][refIdx])
            #
            finalComment = self.__consolidateConflict(comment, self._seqAlignList[i][authIdx][5], self._seqAlignList[i][refIdx][5])
            self._seqAlignList[i][authIdx][5] = finalComment
            if writeConflictFlag:
                _commentType, commentValue = self._decodeComment(finalComment)
                if commentValue in self.__conflictMap:
                    self.__conflictMap[commentValue] += 1
                else:
                    self.__conflictMap[commentValue] = 1
                #
            #
        #

    def __annotateRealConflict(self, authAlignTuple, refAlignTuple):
        """Annotate "modified residue", "engineered mutation", and "conflict" """
        comment = ""
        if authAlignTuple[1] == refAlignTuple[1]:
            comment = ""
        elif self.__isModificationOf(refAlignTuple[1], authAlignTuple[1]):
            # check if "microheterogeneity/modified residue" is applied
            comment = "modified residue"
        else:
            mut1 = refAlignTuple[0] + authAlignTuple[2] + authAlignTuple[0]
            mut2 = refAlignTuple[0] + refAlignTuple[2] + authAlignTuple[0]
            if mut2 in self.__refSeqVariantList:
                comment = "variant"
            elif (mut1 in self._authDefinedMutationList) or (mut2 in self._authDefinedMutationList):
                comment = "engineered mutation"
            else:
                comment = "conflict"
            #
        #
        return comment

    def __consolidateConflict(self, defaultComment, authComment, refComment):
        """Consolidate default & annotated conflicts"""
        _aType, aComment = self._decodeComment(authComment)
        if aComment == "chromophore":
            return authComment
        #
        _rType, rComment = self._decodeComment(refComment)
        if rComment == "chromophore":
            return refComment
        #
        if len(defaultComment) < 1:
            return ""
        #
        if (defaultComment == "insertion") or (defaultComment == "modified residue") or ((len(authComment) < 1) and (len(refComment) < 1)):
            return "DEFAULT:" + defaultComment
        #
        if defaultComment in ("engineered mutation", "conflict", "variant", "microheterogeneity", "linker"):
            if aComment in ("engineered mutation", "conflict", "variant", "microheterogeneity", "linker"):
                return authComment
            elif rComment in ("engineered mutation", "conflict", "variant", "microheterogeneity", "linker"):
                return refComment
            #
        #
        if defaultComment in ("expression tag", "initiating methionine", "insertion", "linker", "cloning artifact", "acetylation", "amidation"):
            if aComment in ("expression tag", "initiating methionine", "insertion", "linker", "cloning artifact", "acetylation", "amidation"):
                return authComment
            elif rComment in ("expression tag", "initiating methionine", "insertion", "linker", "cloning artifact", "acetylation", "amidation"):
                return refComment
            #
        #
        return "DEFAULT:" + defaultComment

    def __isModificationOf(self, refCompId, authCompId):
        """Check if refCompId is the parent residue of authCompId"""
        rCompId = refCompId.upper().strip()
        aCompId = authCompId.upper().strip()
        if (aCompId in self.__parentCompIdMap) and (self.__parentCompIdMap[aCompId] == rCompId):
            return True
        #
        if aCompId in self.__standardList:
            return False
        #
        crpi = ChemRefPathInfo(siteId=self._siteId, verbose=self._verbose, log=self._lfh)
        ccPath = crpi.getFilePath(aCompId, "CC")
        if ccPath and os.access(ccPath, os.F_OK):
            cifObj = mmCIFUtil(filePath=ccPath)
            parentCompId = cifObj.GetSingleValue("chem_comp", "mon_nstd_parent_comp_id")
            if parentCompId:
                pCompId = parentCompId.upper().strip()
                self.__parentCompIdMap[aCompId] = pCompId
                if pCompId == rCompId:
                    return True
                #
            #
        #
        return False

    def __isAuthResidueNameCompatible(self, alignIndex, ResName):
        """ """
        alignTup = self._seqAlignList[alignIndex]
        #
        error = ""
        for idx in range(0, len(alignTup)):
            if idx not in self._reverseSeqAlignLabelIndices:
                continue
            #
            (seqType, seqInstId, _seqPartId, _seqAltId, _seqVersion) = self._getUnpackSeqLabel(self._reverseSeqAlignLabelIndices[idx])
            if seqType != "xyz":
                continue
            #
            if self.__isCompatible(alignTup[idx], ResName):
                continue
            #
            error += "Residue '" + seqInstId + " " + alignTup[idx][1] + " " + alignTup[idx][2] + "' can not be changed to '" + ResName + "'.<br />\n"
        #
        return error

    def __isXyzResidueNameCompatible(self, chainId, srcAlignPos, dstAlignPos):
        """ """
        chainIdx = -1
        for seqId, idx in self._seqAlignLabelIndices.items():
            tL = str(seqId).strip().split("_")
            if len(tL) < 3:
                continue
            #
            if (tL[0] == "xyz") and (tL[1] == chainId):
                chainIdx = idx
                break
            #
        #
        if chainIdx < 0:
            return ""
        #
        srcRecord = self._seqAlignList[srcAlignPos][chainIdx]
        authRecord = self._seqAlignList[dstAlignPos][self._seqAlignLabelIndices[self._authLabel]]
        if not self.__isCompatible(srcRecord, authRecord[1]):
            return (
                "Residue '"
                + chainId
                + " "
                + srcRecord[1]
                + " "
                + srcRecord[2]
                + "' can not be moved to position '"
                + str(dstAlignPos + 1)
                + "' where it is '"
                + authRecord[1]
                + "' in the auth sequence.<br/>\n"
            )
        #
        return ""

    def __isCompatible(self, alignRecord, authResName):
        """Verify the coordinate residue side-chain pattern defined in comment is compatible with new residue name"""
        if ((authResName == "") or (authResName == ".")) and (alignRecord[1] != "") and (alignRecord[1] != "."):
            return False
        #
        if (alignRecord[1] == ".") and (alignRecord[2] == ""):
            return True
        elif authResName == alignRecord[1]:
            return True
        elif (len(alignRecord[5]) > 0) and self._isCompatible(alignRecord[5], authResName):
            return True
        elif (authResName == "ASP") and (alignRecord[1] == "ASN"):
            return True
        elif (authResName == "GLU") and (alignRecord[1] == "GLN"):
            return True
        elif (authResName == "MSE") and ((alignRecord[1] == "MSE") or (alignRecord[1] == "MET")):
            return True
        #
        return False

    def __getPairAndGapList(self, start, end, target):
        """Get moving pair list and leaving gap list"""
        pairList = []
        gapList = []
        if target > end:
            i = end
            j = target
            while i >= start:
                pairList.append((i, j))
                i -= 1
                j -= 1
            #
            gapEnd = end
            if j < gapEnd:
                gapEnd = j
            #
            i = gapEnd
            while i >= start:
                gapList.append(i)
                i -= 1
            #
        else:
            i = start
            j = target
            while i <= end:
                pairList.append((i, j))
                i += 1
                j += 1
            #
            gapStart = start
            if j > gapStart:
                gapStart = j
            #
            i = gapStart
            while i <= end:
                gapList.append(i)
                i += 1
            #
        #
        return pairList, gapList

    def __getPairEditList(self, authIdx, chainIdx, chainSeqId, pairList):
        """Get paired block editlist or un-compatible list"""
        if not pairList:
            return [], []
        #
        (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = self._getUnpackSeqLabel(chainSeqId)
        featureD = self.getFeature(seqType, seqInstId, seqPartId, seqAltId, seqVersion)
        #
        polymerTypeCode = "AA"
        if "POLYMER_TYPE" in featureD:
            polymerTypeCode = featureD["POLYMER_TYPE"]
        #
        editList = []
        unCompatibleList = []
        for pairTup in pairList:
            srcRecord = self._seqAlignList[pairTup[0]][chainIdx]
            tgtRecord = self._seqAlignList[pairTup[1]][chainIdx]
            authRecord = self._seqAlignList[pairTup[1]][authIdx]
            #
            currId = self._getResLabelId(
                seqType=seqType,
                seqInstId=seqInstId,
                seqAltId=seqAltId,
                seqVersion=seqVersion,
                residueCode3=tgtRecord[1],
                residueLabelIndex=tgtRecord[2],
                alignIndex=pairTup[1],
                seqIndex=codeSeqIndex(tgtRecord[3]),
                residueType=polymerTypeCode,
                seqPartId=seqPartId,
            )
            #
            if self.__isCompatible(srcRecord, authRecord[1]):
                newId = self._getResLabelId(
                    seqType=seqType,
                    seqInstId=seqInstId,
                    seqAltId=seqAltId,
                    seqVersion=seqVersion,
                    residueCode3=srcRecord[1],
                    residueLabelIndex=srcRecord[2],
                    alignIndex=pairTup[1],
                    seqIndex=codeSeqIndex(srcRecord[3]),
                    residueType=polymerTypeCode,
                    seqPartId=seqPartId,
                )
                #
                editList.append((currId, srcRecord[0], srcRecord[1], tgtRecord[0], newId))
            else:
                unCompatibleList.append((currId, srcRecord[0], seqInstId, srcRecord[1], srcRecord[2], authRecord[1], str(pairTup[1] + 1)))
            #
        #
        return editList, unCompatibleList

    def __getGapEditList(self, chainIdx, chainSeqId, gapList):
        """Get gap block editlist"""
        if not gapList:
            return []
        #
        (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = self._getUnpackSeqLabel(chainSeqId)
        featureD = self.getFeature(seqType, seqInstId, seqPartId, seqAltId, seqVersion)
        #
        polymerTypeCode = "AA"
        if "POLYMER_TYPE" in featureD:
            polymerTypeCode = featureD["POLYMER_TYPE"]
        #
        editList = []
        for alignPos in gapList:
            xyzRecord = self._seqAlignList[alignPos][chainIdx]
            #
            currId = self._getResLabelId(
                seqType=seqType,
                seqInstId=seqInstId,
                seqAltId=seqAltId,
                seqVersion=seqVersion,
                residueCode3=xyzRecord[1],
                residueLabelIndex=xyzRecord[2],
                alignIndex=alignPos,
                seqIndex=codeSeqIndex(xyzRecord[3]),
                residueType=polymerTypeCode,
                seqPartId=seqPartId,
            )
            #
            newId = self._getResLabelId(
                seqType=seqType,
                seqInstId=seqInstId,
                seqAltId=seqAltId,
                seqVersion=seqVersion,
                residueCode3=self._gapSymbol,
                residueLabelIndex="",
                alignIndex=alignPos,
                seqIndex="",
                residueType=polymerTypeCode,
                seqPartId=seqPartId,
            )
            #
            editList.append((currId, self._gapSymbol, self._gapSymbol, xyzRecord[0], newId))
        #
        return editList
