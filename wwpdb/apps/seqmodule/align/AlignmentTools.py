##
# File:  AlignmentTools.py
# Date:  25-Oct-2018
#
"""
Methods to manage sequence alignments.
"""
__docformat__ = "restructuredtext en"
__author__ = "Zukang Feng"
__email__ = "zfeng@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.09"

import os, sys, traceback

from wwpdb.utils.config.ConfigInfo import ConfigInfo
from wwpdb.apps.seqmodule.align.AlignmentToolUtils import getSeqAlignment,mergeSeqAlignment,decodeIndex,assignConflict
from wwpdb.apps.seqmodule.io.AlignmentDataStore import AlignmentDataStore
from wwpdb.utils.align.alignlib import PseudoMultiAlign
from wwpdb.io.file.mmCIFUtil  import mmCIFUtil

class AlignmentTools(AlignmentDataStore):
    """ Manage sequence alignments
    """
    def __init__(self, reqObj=None, entityId=None, pathInfo=None, seqDataStore=None, deserializeFlag=True, verbose=False, log=sys.stderr):
        super(AlignmentTools, self).__init__(reqObj=reqObj, entityId=entityId, pathInfo=pathInfo, seqDataStore=seqDataStore, \
                                             deserializeFlag=deserializeFlag, verbose=verbose, log=log)
        #
        self.__clearLocalVariables()
        self.__local_authLabel = self._authLabel
        #
        self.__cI = ConfigInfo(siteId=self._siteId, verbose=self._verbose, log=self._lfh)
        self.__ccTopPath = os.path.join(self.__cI.get("SITE_REFDATA_TOP_CVS_SB_PATH"), self.__cI.get("SITE_REFDATA_PROJ_NAME_CC"))
        self.__standardList = ( "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", \
                                "LYS", "MET", "PHE", "PRO", "PYL", "SEC", "SER", "THR", "TRP", "TYR", "VAL", \
                                "UNK", "A", "C", "G", "I", "T", "U", "DA", "DC", "DG", "DI", "DT", "DU" )
        #
        self.__parentCompIdMap = {}

    def __clearLocalVariables(self):
        """ Clear all local variables
        """
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
        """ Input entity ID, auth, xyz & ref labels, and align indices data among auth seq vs coordinate seqs and auth seq vs ref seqs. etc

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
            self.__local_authLabel   = alignD["auth_label"]
            self.__selectedIdList    = alignD["alignids"]
            self.__xyzLabel          = alignD["xyz_label"]
            self.__xyzAlignIndexList = alignD["auth_coord_align_index"]
            self.__partInfoDict      = alignD["part_info"]
            self.__refLabelDict      = alignD["ref_label"]
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
                    for redId,alignIndexList in self.__refAlignIndexDict.items():
                        self._insertRefAlignIndexList([ self.__local_authLabel, redId ], alignIndexList)
                    #
                #
            #
            self._updateAlignmentDataStore(self.__local_authLabel)
        except:
            if (self._verbose):
                traceback.print_exc(file=self._lfh)
            #
            self._addErrorMessage(traceback.format_exc())
        #

    def resetInputAlignData(self, alignD=None):
        """ Input entity ID, auth, xyz & ref labels, and align indices data among auth seq vs coordinate seqs and auth seq vs ref seqs. etc

            alignD["alignids"] = [ "auth_label", "xyz_label_1", "xyz_label_2", ..., "ref_label_1", "ref_label_2", ... ]
            alignD["xyz_label"] = [ "xyz_label_1", "xyz_label_2", ...]
            alignD["part_info"] = { partId_1 : [ beg_num_1, end_num_1], partId_2 : [ beg_num_2, end_num_2], ... }
            alignD["ref_label"] = { partId_1 : ["ref_label_1", "ref_label_2", ...], partId_2 : ["ref_label_3", "ref_label_4", ...], ... }
            alignD["auth_ref_align_index"] = { "ref_label_1" : [ [...], [...], ...], "ref_label_2" : [ [...], [...], ...], ... }
        """
        try:
            self.__clearLocalVariables()
            #
            self.__local_authLabel   = alignD["auth_label"]
            self.__selectedIdList    = alignD["alignids"]
            self.__xyzLabel          = alignD["xyz_label"]
            self.__partInfoDict      = alignD["part_info"]
            self.__refLabelDict      = alignD["ref_label"]
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
                    for redId,alignIndexList in self.__refAlignIndexDict.items():
                        self._insertRefAlignIndexList([ self.__local_authLabel, redId ], alignIndexList)
                    #
                #
                self._updateAlignmentDataStore(self.__local_authLabel)
            #
        except:
            if (self._verbose):
                traceback.print_exc(file=self._lfh)
            #
            self._addErrorMessage(traceback.format_exc())
        #

    def addRefAlignIndices(self, authSeqId="", refSeqId="", partId=0):
        """ Add new ref sequence alignmet index
        """
        if (authSeqId == "") or (refSeqId =="") or (partId == 0):
            return
        #
        indexList = self._getRefAlignIndexList([ authSeqId, refSeqId ])
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
        if not partId in self.__partInfoDict:
            return
        #
        alignIndexList,alignLength,numMatch,numMatchGaps = self.__getAuthRefAlignIndexList(authSeqs, refSeqs, self.__partInfoDict[partId][0] - 1, \
                                                                                           self.__partInfoDict[partId][1] - 1)
        if alignIndexList:
            self._insertRefAlignIndexList([ authSeqId, refSeqId ], alignIndexList)
            #
            (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = self._getUnpackSeqLabel(refSeqId)
            sFeature = self._getFeatureObjByUnpackLabelFromDataStore(seqType, seqInstId, seqPartId, seqAltId, seqVersion)
            sFeature.setAuthRefAlignDetails(seqLen=len(refSeqId), alignLen=int(alignLength), seqSim=float(numMatch) / float(alignLength), \
                                            seqSimWithGaps=float(numMatchGaps) / float(alignLength))
            self.setFeature(sFeature.get(), seqType, seqInstId, seqPartId, seqAltId, seqVersion)
        #
        self.serialize()

    def checkResidueNameReplaceCompatibility(self, alignIndex, origResName, resLabelIndex, newResName):
        """ Check if residue name replacement is allowed based on the aligned coordiante residue(s) 
        """
        if alignIndex >= len(self._seqAlignList):
            return "Compatibility checking for residue name replacement failed.<br />\n"
        #
        alignTup = self._seqAlignList[alignIndex]
        #
        error = ""
        for idx in range(0, len(alignTup)):
            if idx not in self._reverseSeqAlignLabelIndices:
                continue
            #
            (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = self._getUnpackSeqLabel(self._reverseSeqAlignLabelIndices[idx])
            if seqType != "xyz":
                continue
            #
            if (alignTup[idx][1] == ".") and (alignTup[idx][2] == ""):
                continue
            elif newResName == alignTup[idx][1]:
                continue
            elif (len(alignTup[idx][5]) > 0) and self._isCompatible(alignTup[idx][5], newResName):
                continue
            elif (newResName == "ASP") and (alignTup[idx][1] == "ASN"):
                continue
            elif (newResName == "GLU") and (alignTup[idx][1] == "GLN"):
                continue
            elif (newResName == "MSE") and ((alignTup[idx][1] == "MSE") or (alignTup[idx][1] == "MET")):
                continue
            #
            error += "Residue '" + seqInstId + " " + alignTup[idx][1] + " " + alignTup[idx][2] + "' can not be changed to '" + newResName + "'.<br />\n"
        #
        errorMsg = ""
        if error:
            errorMsg = "'" + origResName + " " + resLabelIndex + "' can not be replaced by '" + newResName + "':<br /><br />\n\n" + error
        #
        return errorMsg

    def _checkAndUpdateAlignment(self, inputIdList, selectedIdList):
        """ Check if input alignIds match the default selected alignIds saved in alignment pickle file and
            update alignemnt if necessary
        """
        retInputIdList,retSelectedIdList,allAlignIdList,self.__xyzLabel,refIdList,self._selfRefPartIdList,self.__extraAuthLabel = \
                   self._getProperAlignIdList(inputIdList, selectedIdList)
        #
        if self.__checkCurrentAlignIdMatch(allAlignIdList, self._seqAlignLabelIndices):
            return retInputIdList,retSelectedIdList
        #
        self.__local_authLabel = retSelectedIdList[0]
        #
        if not self.__checkCurrentAlignIdMatch(self.__xyzLabel, self._xyzAlignLabelIndices):
            self.__getXyzAlignIndexList()
            if not self.__getCurrentAuthCoordAlignment():
                self.__alignFlag = False
                return retInputIdList,retSelectedIdList
            #
            self._copyXyzAlignmentToSeqAlignment()
        #
        if (len(refIdList) > 0) and self.__getRefLabelDictFromRefIdList(refIdList):
            self._copyXyzAlignmentToSeqAlignment()
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
        return retInputIdList,retSelectedIdList

    def _isCompatible(self, comment, newResName):
        """ Verify the coordinate residue side-chain pattern defined in comment is compatible with new residue name
        """
        compatibilityMap = { "ala-gly-like" : [ "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", \
                                                "MET", "MSE", "PHE", "PRO", "PYL", "SEC", "SER", "THR", "TRP", "TYR", "UNK", "VAL" ], \
                             "c-gamma-like" : [ "ARG", "ASN", "ASP", "GLN", "GLU", "HIS", "ILE", "LEU", "LYS", "MET", "MSE", "PHE", \
                                                "PYL", "TRP", "THR", "TYR", "VAL" ], \
                             "c-delta-like" : [ "ARG", "GLN", "GLU", "HIS", "ILE", "LYS", "PHE", "PYL", "TRP", "TYR" ], \
                             "ile-val-like" : [ "ILE", "VAL" ], \
                             "ser-thr-like" : [ "SER", "THR" ], \
                             "phe-tyr-like" : [ "PHE", "TYR" ], \
                             "lys-pyl-like" : [ "LYS", "PYL" ] }
    
        #
        if not comment:
            return False
        #
        for patternType,allowList in compatibilityMap.items():
            if (comment.find(patternType) != -1) and (newResName in allowList):
                 return True
            #
        #
        return False

    def _getNaturalSourceWarningMessage(self, sourceType, eelCommentL):
        """ Create warning for natural source entity with [ "engineered mutation", "expression tag", "linker" ] sequence annotation
        """
        if (sourceType.strip().upper() != "NAT") or (len(eelCommentL) < 1):
            return ""
        #
        warningMsg = "Entity '" + self._entityId + "' with natural source has '" + "', '".join(sorted(set(eelCommentL))) \
                   + "' sequence annotation information.<br />"
        return warningMsg

    def _clearAllConflicts(self, authIdx):
        """ Clear all numeric conflicts and conflicting comments
        """
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
        for rowIdx,alignTuple in enumerate(self._seqAlignList):
            for columnIdx,alignRecord in enumerate(alignTuple):
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

    def _assignAllConflicts(self, authLabel, selectedIdList):
        """ Assign numeric conflicts and conflicting comments
        """
        self._missingSeqMap = {}
        #
        authIdx = self._seqAlignLabelIndices[authLabel]
        #
        totalSeqCoodConflict = 0
        for otherIdx in sorted(self._reverseSeqAlignLabelIndices.keys()):
            otherLabel = self._reverseSeqAlignLabelIndices[otherIdx]
            if (not otherLabel.startswith("xyz")) or (otherLabel not in selectedIdList):
                continue
            #
            numConflict,missingAuthSeqList = self.__assignNumericConflicts(authIdx, otherIdx, "xyz", 0, len(self._seqAlignList) - 1)
            totalSeqCoodConflict += numConflict
            if missingAuthSeqList:
                (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = self._getUnpackSeqLabel(otherLabel)
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
            (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = self._getUnpackSeqLabel(otherLabel)
            #
            start_position = 0
            end_position = len(self._seqAlignList) - 1
            if int(seqPartId) in self._partPosDict:
                start_position = self._partPosDict[int(seqPartId)][0]
                end_position = self._partPosDict[int(seqPartId)][1]
            #
            numConflict,missingAuthSeqList = self.__assignNumericConflicts(authIdx, otherIdx, "ref", start_position, end_position)
            self.__annotateConflictingComments(authIdx, otherIdx, start_position, end_position)
        #
        return totalSeqCoodConflict

    def _checkPartStartEndPosMap(self, authIdx, authLabel):
        """ Get auth sequence part start & end pos range in alignment if it does not exist
        """
        if self._partPosDict:
            return
        #
        self.__checkPartInfoDict(authLabel)
        #
        if not self.__partInfoDict:
            self._partPosDict = { 1: [ 0, len(self._seqAlignList) - 1 ] }
            return
        #
        SeqIdxPartIdMap = {}
        for partId,rangePair in self.__partInfoDict.items():
            for i in range(rangePair[0] - 1, rangePair[1]):
                SeqIdxPartIdMap[i] = int(partId)
            #
        #
        partId = 1
        self._partPosDict = {}
        self._partPosDict[partId] = [ 0, 0 ]
        #
        for idx,alignTup in enumerate(self._seqAlignList):
            seqIdx,comment = decodeIndex(alignTup[authIdx][3])
            if not seqIdx in SeqIdxPartIdMap:
                continue
            #
            oldPartId = SeqIdxPartIdMap[seqIdx]
            if oldPartId != partId:
                self._partPosDict[partId][1] = idx - 1
                partId = oldPartId
                self._partPosDict[partId] = [ idx, idx ]
            #
        #
        self._partPosDict[partId][1] = len(self._seqAlignList) - 1

    def _updateAlignmentDataStore(self, newAuthLabel):
        """ Update persist alignment storage
        """
        if self.__alignFlag:
            self._authLabel = newAuthLabel
            self._hasAlignmentFlag = True
        else:
            self._resetAlignmentInfo()
            self._addErrorMessage("Failed to generate alignment for entity '" + self._entityId + "'.")
        #
        self.serialize()

    def _decodeComment(self, comment):
        """ Remove prefix "DEFAULT:" or ""ANNOTATED:" from comment string
        """
        if str(comment).strip().startswith("ANNOTATED:"):
            return "ANNOTATED",str(comment).strip()[10:].lower()
        elif str(comment).strip().startswith("SELECTED:"):
            return "SELECTED",str(comment).strip()[9:].lower()
        elif str(comment).strip().startswith("DEFAULT:"):
            return "DEFAULT",str(comment).strip()[8:].lower()
        else:
            return "",str(comment).strip().lower()
        #

    def __checkCurrentAlignIdMatch(self, alignIdList, alignIndices):
        """
        """
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
        """ Generate pseduo multiple alignment among auth, coordiates and reference sequences
        """
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
        """ Generate pseduo multiple alignment among auth, coordiates and reference sequences
        """
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
        """ Get sequence from sequence data storage using packed id label
        """
        seq = self._getSequenceByPackLabelFromDataStore(idLabel)
        if self._verbose:
            self._lfh.write("Instance '%r' sequence length=%d\n" % (idLabel, len(seq)))
        #
        if not seq:
            self.__alignFlag = False
        #
        return seq

    def __getCurrentAuthCoordAlignment(self, existAuthSeqs=[]):
        """ Generate auth/coordiates sequence alignment
        """
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
            if not seqs:
                return False
            #
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
        errMsg,self._xyzAlignList = getSeqAlignment(seqList, self.__xyzAlignIndexList, self._reverseXyzAlignLabelIndices, "xyz", self._gapSymbol)
        if errMsg:
            if self._verbose:
                self._lfh.write("AlignmentTools.__getCurrentAuthCoordAlignment errMsg=%r\n" % errMsg)
            #
            self._addErrorMessage(errMsg)
            return False
        #
        self.__getAuthDefinedMutations()
        return True

    def __getCurrentAuthRefAlignment(self, existAuthSeqs=[]):
        """ Generate auth/ref sequence alignment(s)
        """
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
                if not label in self.__refAlignIndexDict:
                    if self._verbose:
                        self._lfh.write("No reference sequence align index for %s\n" % label)
                    #
                    self._addErrorMessage("No reference sequence align index for '" + label + "'.")
                    return False;
                #
                reverseIndices = {}
                reverseIndices[0] = self.__local_authLabel
                reverseIndices[1] = label
                seqList = [ authSeqs, seqs ]
                #
                errMsg,seqAlignList = getSeqAlignment(seqList, self.__refAlignIndexDict[label], reverseIndices, "ref", self._gapSymbol)
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
        errMsg,self._seqAlignList = mergeSeqAlignment(multiAlignList, pA.getAlignIndices(), self._gapSymbol)
        if errMsg:
            self._addErrorMessage(errMsg)
            return False
        #
        return True

    def __alignExtraAuthSequence(self, existAuthSeqs=[]):
        """
        """
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
        errMsg,seqAlignList = getSeqAlignment(seqList, indexList, reverseIndices, [], "auth", self._gapSymbol)
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
        errMsg,self._seqAlignList = mergeSeqAlignment(multiAlignList, pA.getAlignIndices(), self._gapSymbol)
        if errMsg:
            if self._verbose:
                self._lfh.write("AlignmentTools.__alignExtraAuthSequence errMsg=%r\n" % errMsg)
            #
            self._addErrorMessage(errMsg)
            self.__alignFlag = False
        #

    def __getXyzAlignIndexList(self):
        """ Re-generate self.__xyzAlignIndexList based on selected self.__xyzLabel list
        """
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
        """
        """
        alignIndexList = []
        for i in range(0, startIdx):
            alignIndexList.append([ i, -1 ])
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
            alignIndexList.append([ i, -1 ])
        #
        return alignIndexList,alignLength,numMatch,numMatchGaps

    def __getRefLabelDictFromRefIdList(self, refIdList):
        """ Get self.__refLabelDict from input refIdList
        """
        self.__refLabelDict = {}
        self.__refAlignIndexDict = {}
        #
        if (not refIdList) or (not self.__alignFlag):
            return False
        #
        for refId in refIdList:
            indexList = self._getRefAlignIndexList([self.__local_authLabel, refId])
            #
            (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = self._getUnpackSeqLabel(refId)
            #
            if len(indexList) == 0:
                authSeqs = self.__getSequenceFromDataStore(self.__local_authLabel)
                refSeqs = self.__getSequenceFromDataStore(refId)
                self.__checkPartInfoDict(self.__local_authLabel)
                if authSeqs and refSeqs and (seqPartId in self.__partInfoDict):
                    indexList,alignLength,numMatch,numMatchGaps = self.__getAuthRefAlignIndexList(authSeqs, refSeqs, \
                                   self.__partInfoDict[seqPartId][0] - 1, self.__partInfoDict[seqPartId][1] - 1)
                    if len(indexList) > 0:
                        self._insertRefAlignIndexList([ self.__local_authLabel, refId ], indexList)
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
                self.__refLabelDict[seqPartId] = [ refId ]
            #
            self.__refAlignIndexDict[refId] = indexList
        #
        return True

    def __getAuthDefinedMutations(self):
        """ Extract all author input mutation information (e.g. V178A ) from _entity.pdbx_mutation
        """
        featureD = self._getFeatureByPackLabelFromDataStore(self.__local_authLabel)
        if not "ENTITY_MUTATION_DETAILS" in featureD:
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

    def __generateInputSeqInfoForPseudoMultiAlignFromSeqList(self, seqType="auth", inputSeqList=[], start=0, end=0):
        """ Generate input sequence information for PseudoMultiAlign program
        """
        if end == 0:
            end = len(inputSeqList) - 1
        #
        returnSeqList = []
        for (i, seqTup) in enumerate(inputSeqList):
            if (i < start) or (i > end):
                continue
            #
            if seqType == "xyz":
                returnSeqList.append( ( str(seqTup[0]).upper(), str(i), str(str(seqTup[6])) ) )
            else:
                returnSeqList.append( ( str(seqTup[0]).upper(), str(i) ) )
            #
        #
        return returnSeqList

    def __generateInputSeqInfoForPseudoMultiAlignFromAlignList(self, alignList):
        """ Generate input sequence information for PseudoMultiAlign program
        """
        seqList = []
        for (i, alignTup) in enumerate(alignList):
            seqList.append( ( str(alignTup[0][1]).upper(), str(i) ) )
        #
        return seqList

    def __checkPartInfoDict(self, authLabel):
        """ Get Parts information if it does not exist
        """
        if self.__partInfoDict:
            return
        #
        self.__partInfoDict = self._getPartInfoDict(authLabel)
        #
        if self.__partInfoDict:
            return
        #
        (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = self._getUnpackSeqLabel(authLabel)
        parIdList = self.getPartIdList(seqType, seqInstId)
        for partId in parIdList:
            seqFeatureD = self.getFeature(seqType, seqInstId, partId, seqAltId, seqVersion)
            self.__partInfoDict[int(seqFeatureD['AUTH_SEQ_PART_ID'])] = ( int(seqFeatureD['AUTH_SEQ_NUM_BEGIN']), int(seqFeatureD['AUTH_SEQ_NUM_END']) )
        #
        self._insertPartInfoDict(authLabel, self.__partInfoDict)
        
    def __assignNumericConflicts(self, authIdx, otherIdx, otherType, start_position, end_position):
        """ Assign numeric conflicts
        """
        numConflict = 0
        missingAuthSeqList = []
        for intIdx,alignTup in enumerate(self._seqAlignList):
            if (intIdx < start_position) or (intIdx > end_position):
                continue
            #
            if (alignTup[authIdx][1] == self._gapSymbol) and (alignTup[otherIdx][1] == self._gapSymbol):
                continue
            #
            isConflict,authConflict,otherConflict = assignConflict("auth", alignTup[authIdx][1], otherType, alignTup[otherIdx][1], self._gapSymbol)
            if isConflict:
                if alignTup[otherIdx][1] != self._gapSymbol:
                    numConflict += 1
                #
                if otherType == "xyz" and (alignTup[authIdx][1] == self._gapSymbol):
                    missingAuthSeqList.append(( alignTup[otherIdx][1], alignTup[otherIdx][2] ))
                #
                alignTup[authIdx][6] = authConflict[0]
                alignTup[authIdx][7] = authConflict[1]
                alignTup[otherIdx][6] = otherConflict[0]
                alignTup[otherIdx][7] = otherConflict[1]
            #
        #
        return numConflict,missingAuthSeqList

    def __annotateConflictingComments(self, authIdx, refIdx, beg_pos, end_pos):
        """ Add default 'expression tag annotation comments' for leading or trailing gap residues, and
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
        #

    def __annotateRealConflict(self, authAlignTuple, refAlignTuple):
        """ Annotate "modified residue", "engineered mutation", and "conflict"
        """
        comment = ""
        if authAlignTuple[1] == refAlignTuple[1]:
            comment = ""
        elif self.__isModificationOf(refAlignTuple[1], authAlignTuple[1]):
            # check if "microheterogeneity/modified residue" is applied
            comment = "modified residue"
        else:
            mut1 = refAlignTuple[0] + authAlignTuple[2] + authAlignTuple[0]
            mut2 = refAlignTuple[0] + refAlignTuple[2] + authAlignTuple[0]
            if (mut1 in self._authDefinedMutationList) or (mut2 in self._authDefinedMutationList):
                comment = "engineered mutation"
            else:
                comment = "conflict"
            #
        #
        return comment

    def __consolidateConflict(self, defaultComment, authComment, refComment):
        """ Consolidate default & annotated conflicts
        """
        if len(defaultComment) < 1:
            return ""
        #
        if (defaultComment == "insertion") or (defaultComment == "modified residue") or ((len(authComment) < 1) and (len(refComment) < 1)):
            return "DEFAULT:" + defaultComment
        #
        aType,aComment = self._decodeComment(authComment)
        rType,rComment = self._decodeComment(refComment)
        #
        if defaultComment in ( "engineered mutation", "conflict", "variant" ):
            if aComment in ( "engineered mutation", "conflict", "variant" ):
                return authComment
            elif rComment in ( "engineered mutation", "conflict", "variant" ):
                return refComment
            #
        #
        if defaultComment in ( "expression tag", "initiating methionine", "insertion", "linker", "cloning artifact", "acetylation", "amidation" ):
            if aComment in ( "expression tag", "initiating methionine", "insertion", "linker", "cloning artifact", "acetylation", "amidation" ):
                return authComment
            elif rComment in ( "expression tag", "initiating methionine", "insertion", "linker", "cloning artifact", "acetylation", "amidation" ):
                return refComment
            #
        #
        return "DEFAULT:" + defaultComment

    def __isModificationOf(self, refCompId, authCompId):
        """ Check if refCompId is the parent residue of authCompId
        """
        rCompId = refCompId.upper().strip()
        aCompId = authCompId.upper().strip()
        if (aCompId in self.__parentCompIdMap) and (self.__parentCompIdMap[aCompId] == rCompId):
            return True
        #
        if aCompId in self.__standardList:
            return False
        #
        ccPath = os.path.join(self.__ccTopPath, aCompId[0], aCompId, aCompId + ".cif")
        if os.access(ccPath, os.F_OK):
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
