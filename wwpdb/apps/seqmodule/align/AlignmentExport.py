##
# File:  AlignmentExport_v2.py
# Date:  31-Oct-2018
#
##
"""
Export alignment details to production RCSB pipeline.

"""
__docformat__ = "restructuredtext en"
__author__ = "Zukang Feng"
__email__ = "zfeng@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.07"

import sys

from wwpdb.apps.seqmodule.align.AlignmentTools import AlignmentTools
from wwpdb.apps.seqmodule.align.AlignmentToolUtils import decodeIndex


class AlignmentExport(AlignmentTools):
    """This class encapsulates all of the data export operations of alignment data to the RCSB/WF data pipeline."""

    def __init__(self, reqObj=None, entityId=None, pathInfo=None, seqDataStore=None, verbose=False, log=sys.stderr):
        super(AlignmentExport, self).__init__(reqObj=reqObj, entityId=entityId, pathInfo=pathInfo, seqDataStore=seqDataStore, verbose=verbose, log=log)

    def doExport(self, idAuthSeq, idListRef, idListXyz, selfRefList, sourceType):
        """Export alignment data"""
        alignIdList = []
        selectedIdList = []
        alignIdList.append(idAuthSeq)
        selectedIdList.append(idAuthSeq)
        if len(idListRef) > 0:
            alignIdList.extend(idListRef)
            selectedIdList.extend(idListRef)
        #
        if len(idListXyz) > 0:
            alignIdList.extend(idListXyz)
            selectedIdList.extend(idListXyz)
        #
        if len(selfRefList) > 0:
            selectedIdList.extend(selfRefList)
        #
        self._checkAndUpdateAlignment(alignIdList, selectedIdList)
        authIdx = self._seqAlignLabelIndices[idAuthSeq]
        self._checkPartStartEndPosMap(authIdx, idAuthSeq)
        self._clearAllConflicts(authIdx)
        self._assignAllConflicts(idAuthSeq, alignIdList)
        #
        errorFlag = False
        for instId in alignIdList:
            if instId not in self._seqAlignLabelIndices:
                errorFlag = True
                self._lfh.write("Instance '%s' does not exist in alignment for entity '%s'\n" % (instId, self._entityId))
            #
        #
        if errorFlag:
            return [], [], [], [], [], "", 0
        #
        rptRefL, rptCommentL, rptCommentModL, rptDeleteL, warningMsg = self.__alignRefReport(self._seqAlignLabelIndices[idAuthSeq], idListRef, sourceType)
        rptXyzL, numConflicts = self.__alignXyzReport(self._seqAlignLabelIndices[idAuthSeq], idListXyz)
        #
        return rptRefL, rptCommentL, rptCommentModL, rptDeleteL, rptXyzL, warningMsg, numConflicts

    def __alignRefReport(self, authIdx, idListRef, sourceType):
        """Export all reference alignments"""
        rptRefL = []
        rptCommentL = []
        rptCommentModL = []
        rptDeleteL = []
        if len(idListRef) < 1:
            return rptRefL, rptCommentL, rptCommentModL, rptDeleteL, ""
        #
        eelCommentL = []
        refMap = {}
        expressionTagCount = 0
        foundLongExpressionTag = False
        #
        for idRef in idListRef:
            refMap[idRef] = self._getUnpackSeqLabel(idRef)
        #
        for idx, alignTup in enumerate(self._seqAlignList):
            entity_mon_id = str(alignTup[authIdx][1])
            entity_seq_num = str(alignTup[authIdx][2])
            #
            tstCompId = self._gapSymbol
            tstSeqNum = "."
            tstPartId = "."
            featureD = {}
            for idRef in idListRef:
                if (idRef not in self._seqAlignLabelIndices) or (alignTup[self._seqAlignLabelIndices[idRef]][1] == self._gapSymbol):
                    continue
                #
                (seqTypeT, seqInstIdT, seqPartIdT, seqAltIdT, seqVersionT) = refMap[idRef]
                tstCompId = alignTup[self._seqAlignLabelIndices[idRef]][1]
                tstSeqNum = alignTup[self._seqAlignLabelIndices[idRef]][2]
                tstPartId = seqPartIdT
                featureD = self.getFeature(seqTypeT, seqInstIdT, seqPartIdT, seqAltIdT, seqVersionT)
                break
            #
            if len(alignTup[authIdx][5]) > 0:
                _cType, comment = self._decodeComment(alignTup[authIdx][5])
                #
                if comment == "expression tag":
                    expressionTagCount += 1
                else:
                    if expressionTagCount >= self._longExpressionTagCountCutoff:
                        foundLongExpressionTag = True
                    #
                    expressionTagCount = 0
                #
                if sourceType.strip().upper() == "NAT":
                    for ccType in ("engineered mutation", "expression tag", "linker"):
                        if comment.find(ccType) != -1:
                            eelCommentL.append(ccType)
                        #
                    #
                #
                if comment in ["modified residue", "microheterogeneity/modified residue"]:
                    irow = len(rptCommentModL) + 1
                    rptCommentModL.append([str(irow), str(self._entityId), entity_mon_id, entity_seq_num, "modified residue"])
                #
                if (comment == "chromophore") and ((not entity_mon_id) or entity_mon_id == ".") and ((not entity_seq_num) or entity_seq_num == "."):
                    start = idx - 2
                    if start < 0:
                        start = 0
                    #
                    end = idx + 3
                    if end > len(self._seqAlignList):
                        end = len(self._seqAlignList)
                    #
                    for idx1 in range(start, end):
                        _cType1, comment1 = self._decodeComment(self._seqAlignList[idx1][authIdx][5])
                        if (
                            (comment1 == "chromophore")
                            and str(self._seqAlignList[idx1][authIdx][1])
                            and (str(self._seqAlignList[idx1][authIdx][1]) != ".")
                            and str(self._seqAlignList[idx1][authIdx][2])
                            and (str(self._seqAlignList[idx1][authIdx][2]) != ".")
                        ):
                            entity_mon_id = str(self._seqAlignList[idx1][authIdx][1])
                            entity_seq_num = str(self._seqAlignList[idx1][authIdx][2])
                            break
                        #
                    #
                elif comment not in ("modified residue", "deletion"):
                    if comment == "microheterogeneity/modified residue":
                        comment = "microheterogeneity"
                    #
                    irow = len(rptCommentL) + 1
                    rptCommentL.append([str(irow), str(self._entityId), entity_mon_id, entity_seq_num, comment])
                elif (comment in ("Deletion", "deletion")) and featureD and (tstCompId != self._gapSymbol) and (tstSeqNum != ".") and (tstPartId != "."):
                    irow = len(rptDeleteL) + 1
                    rptDeleteL.append(
                        [
                            str(irow),
                            str(self._entityId),
                            str(tstPartId),
                            self._srd.convertDbNameToResource(featureD["DB_NAME"]),
                            featureD["DB_ACCESSION"],
                            featureD["DB_ISOFORM"],
                            tstCompId,
                            tstSeqNum,
                        ]
                    )
                #
            else:
                if expressionTagCount >= self._longExpressionTagCountCutoff:
                    foundLongExpressionTag = True
                #
                expressionTagCount = 0
            #
            irow = len(rptRefL) + 1
            rptRefL.append([str(irow), str(self._entityId), entity_mon_id, entity_seq_num, str(self._entityId), str(tstCompId), str(tstSeqNum), str(tstPartId)])
            #
        #
        if expressionTagCount >= self._longExpressionTagCountCutoff:
            foundLongExpressionTag = True
        #
        warningMsg = ""
        if foundLongExpressionTag:
            warningMsg += "Entity '" + self._entityId + "' has a long expression tag.<br />"
        #
        warningMsg += self._getNaturalSourceWarningMessage(sourceType, eelCommentL)
        #
        return rptRefL, rptCommentL, rptCommentModL, rptDeleteL, warningMsg

    def __alignXyzReport(self, authIdx, idListXyz):
        """Export all coordinate alignments"""
        rptXyzL = []
        numConflicts = 0
        for idXyz in idListXyz:
            if idXyz not in self._seqAlignLabelIndices:
                continue
            #
            xyzSeqList = self._getSequenceByPackLabelFromDataStore(idXyz)
            xyzIdx = self._seqAlignLabelIndices[idXyz]
            (_seqTypeT, seqInstIdT, _seqPartIdT, _seqAltIdT, _seqVersionT) = self._getUnpackSeqLabel(idXyz)
            #
            for alignTup in self._seqAlignList:
                if (str(alignTup[authIdx][1]) in (".", "?", "")) and (str(alignTup[xyzIdx][1]) in (".", "?", "")):
                    continue
                #
                seqIdx, _comment = decodeIndex(alignTup[xyzIdx][3])
                orgName = "."
                if (seqIdx >= 0) and (seqIdx < len(xyzSeqList)):
                    orgName = xyzSeqList[seqIdx][5]
                #
                if str(alignTup[xyzIdx][1]) in (".", "?", ""):
                    tstResId = self._gapSymbol
                    resIdx = "."
                    insCode = "."
                else:
                    tstResId = str(alignTup[xyzIdx][1])
                    # tstResIdx = alignTup[xyzIdx][2]
                    if str(alignTup[xyzIdx][2])[-1].isdigit():
                        insCode = "."
                        resIdx = str(alignTup[xyzIdx][2])
                    else:
                        insCode = str(alignTup[xyzIdx][2])[-1]
                        resIdx = str(alignTup[xyzIdx][2])[:-1]
                    #
                #
                if (len(alignTup[xyzIdx][5]) > 0) and (alignTup[xyzIdx][5].find("hetero") != -1):
                    ii = alignTup[xyzIdx][5].find("hetero-")
                    tS = alignTup[xyzIdx][5][ii + 7 :]
                    hL = tS.split(":")
                    for h in hL:
                        irow = len(rptXyzL) + 1
                        rptXyzL.append([str(irow), str(self._entityId), str(alignTup[authIdx][1]), str(alignTup[authIdx][2]), str(seqInstIdT), h, resIdx, insCode, "Y", orgName])
                    #
                else:
                    if (tstResId is not self._gapSymbol) and (str(alignTup[authIdx][1]) != tstResId):
                        numConflicts += 1
                    #
                    irow = len(rptXyzL) + 1
                    rptXyzL.append([str(irow), str(self._entityId), str(alignTup[authIdx][1]), str(alignTup[authIdx][2]), str(seqInstIdT), tstResId, resIdx, insCode, "N", orgName])
                #
            #
        #
        return rptXyzL, numConflicts
