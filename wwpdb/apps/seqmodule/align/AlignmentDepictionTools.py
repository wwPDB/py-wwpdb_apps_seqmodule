##
# File:  AlignmentDepictionTools.py
# Date:  27-Oct-2018
#
"""
Controlling class for the production of sequence alignment views.
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

from wwpdb.apps.seqmodule.align.AlignmentBackEndEditingTools import AlignmentBackEndEditingTools
from wwpdb.apps.seqmodule.align.AlignmentToolUtils import codeSeqIndex
from wwpdb.apps.seqmodule.io.SequenceEditStore import getEditStoreFilename, SequenceEditStore
from wwpdb.utils.session.UtilDataStore import UtilDataStore


class AlignmentDepictionTools(AlignmentBackEndEditingTools):
    """Render the input alignment view and conflict table."""

    def __init__(self, reqObj=None, entityId=None, pathInfo=None, seqDataStore=None, verbose=False, log=sys.stderr):
        super(AlignmentDepictionTools, self).__init__(reqObj=reqObj, entityId=entityId, pathInfo=pathInfo, seqDataStore=seqDataStore, verbose=verbose, log=log)
        #
        self.__sessionId = self._reqObj.getSessionId()
        self.__alignViewOrder = self._reqObj.getAlignmentOrdering()
        self.__alignIdList = self._reqObj.getAlignList()
        self.__allSelectedIdList = str(self._reqObj.getValue("selectids")).split(",")
        self.__operation = str(self._reqObj.getValue("operation"))
        self.__identifier = str(self._reqObj.getValue("identifier"))
        #
        self.__seqTypeAuth = ""
        self.__seqInstIdAuth = ""
        self.__seqPartIdAuth = 1
        self.__seqAltIdAuth = 1
        self.__seqVersionAuth = 1
        self.__polymerTypeCodeAuth = "AA"
        #
        self.__orderAlignIndexList = []
        self.__orderAlignRefIndexList = []
        self.__misMatchTypes = []
        self.__alignmentBlock = []
        self.__pdbChainIdList = []
        self.__missingAuthSeqMap = {}
        self.__warningMsg = ""
        #
        self.__refAlignIdx = -1
        self.__authAlignIdx = -1
        #
        self.__selectedIdList = []
        self.__annotationDict = {}

    def doRender(self):
        """Render alignment, conflict table, annotations & warning"""
        self.__setup()
        authIdx = self._seqAlignLabelIndices[self._authLabel]
        self._checkPartStartEndPosMap(authIdx, self._authLabel)
        self.__buildAlignIndexOrder()
        self._clearAllConflicts(authIdx)
        totalSeqCoodConflict = self._assignAllConflicts(self._authLabel, self.__selectedIdList)
        self.serialize()
        if len(self.__alignIdList) > len(self.__selectedIdList):
            self._clearAllConflicts(authIdx)
            _numSeqCoodConflict = self._assignAllConflicts(self._authLabel, self.__alignIdList, writeConflictFlag=False)  # noqa: F841
        #
        self.__renderAlignment()
        self.__renderConflictTable(authIdx)
        self.__renderAnnotations()
        self.__renderWarning(authIdx)
        # self.__writeMisMatchTypeText()
        self.__updateSeqCoodConflict(totalSeqCoodConflict)
        self.__updateMiscDataStore()
        self.__updateUtilDataStore()
        if len(self._selfRefPartIdList) > 0:
            selectedIdList = copy.deepcopy(self.__selectedIdList)
            selectedIdList.extend(self._selfRefPartIdList)
            self._updateDefaultSelections(selectedIdList)
            if self._verbose:
                self._lfh.write("self.__selectedIdList=%r\n" % self.__selectedIdList)
                self._lfh.write("self._selfRefPartIdList=%r\n" % self._selfRefPartIdList)
                self._lfh.write("selectedIdList=%r\n" % selectedIdList)
            #
        else:
            self._updateDefaultSelections(self.__selectedIdList)
            if self._verbose:
                self._lfh.write("self.__selectedIdList=%r\n" % self.__selectedIdList)
            #
        #
        return True

    def setLogHandle(self, log=sys.stderr):
        """Reset the stream for logging output."""
        try:
            self._lfh = log
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            return False
        #

    def __setup(self):
        """Check if input alignIds match the default selected alignIds saved in alignment pickle file and
        build proper display order index based input alignment view order
        """
        chainIdList = self.getGroup(self._entityId)
        if self._verbose:
            self._lfh.write("AlignmentDepictionTools.__setup identifier=%r entityId=%r operation=%r\n" % (self.__identifier, self._entityId, self.__operation))
            self._lfh.write("alignIdList=%r\n" % self.__alignIdList)
            self._lfh.write("allSelectedIdList=%r\n" % self.__allSelectedIdList)
            self._lfh.write("chainIdList=%r\n" % chainIdList)
        #
        if self.__operation == "re-load":
            self.__selectedIdList = []
            self._selfRefPartIdList = []
            for seqId in self.__allSelectedIdList:
                tL = str(seqId).strip().split("_")
                if len(tL) < 3:
                    continue
                #
                if ((tL[0] in ("auth", "ref")) and (tL[1] == self._entityId)) or ((tL[0] == "xyz") and (tL[1] in chainIdList)):
                    self.__selectedIdList.append(seqId)
                elif (tL[0] == "selfref") and (tL[1] == self._entityId):
                    self._selfRefPartIdList.append(seqId)
                #
            #
            self.__editAndUpdateAlignment()
        else:
            self.__alignIdList, self.__selectedIdList = self._checkAndUpdateAlignment(self.__alignIdList, self.__allSelectedIdList)
        #
        if self._verbose:
            self._lfh.write("self.__selectedIdList=%r\n" % self.__selectedIdList)
        #

    def __editAndUpdateAlignment(self):
        """Apply any stored edits to the current alignment in the input sequence alignment list."""
        esfn = getEditStoreFilename(self._entityId)
        ses = SequenceEditStore(sessionObj=self._reqObj.getSessionObj(), fileName=esfn, verbose=self._verbose, log=self._lfh)
        edObjList = ses.getList()
        if self._verbose:
            self._lfh.write("+AlignmentDepictionTools.__editAndUpdateAlignment() edit store %s edit list length %d\n" % (esfn, len(edObjList)))
        #
        edObjMap = {}
        for edObj in edObjList:
            edType = edObj.getEditType()
            resLabelId = edObj.getTargetElementId()
            resLabelObj = self._getResLabelFromResLabelId(resLabelId)
            if not resLabelObj:
                if self._verbose:
                    self._lfh.write("+AlignmentDepictionTools.__editAndUpdateAlignment() unpack resLabelId failed: edType=%r getTargetElementId()=%r\n" % (edType, resLabelId))
                #
                continue
            #
            if edType == "replaceid":
                newId = edObj.getNewElementId()
                newLabelObj = self._getResLabelFromResLabelId(newId)
                if not newLabelObj:
                    if self._verbose:
                        self._lfh.write("+AlignmentDepictionTools.__editAndUpdateAlignment() unpack resLabelId failed: edType=%r getNewElementId()=%r\n" % (edType, newId))
                    #
                    continue
                #
            #
            seqLabelId = self._getSeqLabelId(resLabelObj.getSeq())
            if seqLabelId not in self._seqAlignLabelIndices:
                if self._verbose:
                    self._lfh.write("+AlignmentDepictionTools.__editAndUpdateAlignment() seqLabelId %r does not exist in alignment.\n" % seqLabelId)
                #
                continue
            #
            seqType = resLabelObj.getSequenceType()
            if edType in edObjMap:
                if seqType in edObjMap[edType]:
                    edObjMap[edType][seqType].append((resLabelId, (seqType, seqLabelId, self._seqAlignLabelIndices[seqLabelId], edObj)))
                else:
                    edObjMap[edType][seqType] = [(resLabelId, (seqType, seqLabelId, self._seqAlignLabelIndices[seqLabelId], edObj))]
                #
            else:
                myD = {}
                myD[seqType] = [(resLabelId, (seqType, seqLabelId, self._seqAlignLabelIndices[seqLabelId], edObj))]
                edObjMap[edType] = myD
            #
        #
        if not edObjMap:
            return
        #
        for edType in ("details", "replace", "replaceid", "delete", "insert"):
            if edType not in edObjMap:
                continue
            #
            for seqType in ("auth", "xyz", "ref"):
                if seqType not in edObjMap[edType]:
                    continue
                #
                for (resLabelId, editInfoTuple) in edObjMap[edType][seqType]:
                    try:
                        getattr(self, "_update_%s" % edType)(self._getResLabelFromResLabelId(resLabelId), editInfoTuple)
                    except:  # noqa: E722 pylint: disable=bare-except
                        traceback.print_exc(file=self._lfh)
                    #
                #
            #
        #
        self.__alignIdList, self.__selectedIdList = self._updateAlignmentAndSequences(self.__alignIdList, self.__selectedIdList)
        self._updateAlignmentDataStore(self._authLabel)
        ses.deleteEditList(edObjList)

    def __buildAlignIndexOrder(self):
        """Build display order index based input alignment view order"""
        orderedAlignIdList = []
        for viewOrder in self.__alignViewOrder.split("-"):
            foundAlignIdList = []
            for alignId in self.__alignIdList:
                if not alignId.startswith(viewOrder):
                    continue
                #
                if alignId not in self._seqAlignLabelIndices:
                    self._addErrorMessage("Missing '" + alignId + "' in alignment.")
                    continue
                #
                foundAlignIdList.append([self._seqAlignLabelIndices[alignId], alignId])
            #
            if not foundAlignIdList:
                continue
            #
            if len(foundAlignIdList) > 1:
                foundAlignIdList.sort(key=itemgetter(0))
            #
            orderedAlignIdList.extend(foundAlignIdList)
        #
        if not orderedAlignIdList:
            return
        #
        num_refs = 0
        num_auths = 0
        for orderedAlignIdTup in orderedAlignIdList:
            (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = self._getUnpackSeqLabel(orderedAlignIdTup[1])
            featureD = self.getFeature(seqType, seqInstId, seqPartId, seqAltId, seqVersion)
            #
            polymerTypeCode = "AA"
            if "POLYMER_TYPE" in featureD:
                polymerTypeCode = featureD["POLYMER_TYPE"]
            #
            if seqType == "ref":
                num_refs += 1
                self.__refAlignIdx = orderedAlignIdTup[0]
                #
                dbName = ""
                if "DB_NAME" in featureD:
                    dbName = featureD["DB_NAME"]
                #
                displayCode = ""
                uniprotCode = ""
                if "DB_ACCESSION" in featureD:
                    displayCode = featureD["DB_ACCESSION"]
                    if dbName in ("UNP", "TR", "SP"):
                        uniprotCode = featureD["DB_ACCESSION"]
                    #
                #
                if (dbName in ("UNP", "TR", "SP")) and ("DB_ISOFORM" in featureD) and (len(featureD["DB_ISOFORM"]) > 0):
                    displayCode = featureD["DB_ISOFORM"]
                #
                label = dbName + ":" + displayCode
                if uniprotCode != "":
                    label = '<a href="http://www.uniprot.org/uniprot/%s" target="blank"><span class="nowrap">%s</span></a>' % (uniprotCode, dbName + ":" + displayCode)
                #
            elif seqType == "auth":
                num_auths += 1
                self.__authAlignIdx = orderedAlignIdTup[0]
                #
                label = seqType.upper() + " Entity: " + self._entityId + " V:" + str(seqVersion)
                #
                self.__seqTypeAuth = seqType
                self.__seqInstIdAuth = seqInstId
                self.__seqPartIdAuth = seqPartId
                self.__seqAltIdAuth = seqAltId
                self.__seqVersionAuth = seqVersion
                self.__polymerTypeCodeAuth = polymerTypeCode
            else:
                label = seqType.upper() + " Chain: " + seqInstId + " V:" + str(seqVersion)
                self.__pdbChainIdList.append(seqInstId)
            #
            self.__orderAlignIndexList.append((orderedAlignIdTup[0], seqType, seqInstId, seqPartId, seqAltId, seqVersion, label, polymerTypeCode))
            if seqType == "ref":
                begPos = 0
                endPos = len(self._seqAlignList) - 1
                if int(seqPartId) in self._partPosDict:
                    begPos = self._partPosDict[int(seqPartId)][0]
                    endPos = self._partPosDict[int(seqPartId)][1]
                #
                self.__orderAlignRefIndexList.append((orderedAlignIdTup[0], seqType, seqInstId, seqPartId, seqAltId, seqVersion, label, polymerTypeCode, begPos, endPos))
            #
        #
        if (num_refs > 1) or (num_auths > 1):
            self.__refAlignIdx = -1
            self.__authAlignIdx = -1
        #

    def __renderAlignment(self):
        """Render the alignment list according the selected option type"""
        htmlFilePath = os.path.join(self._sessionPath, "current-alignment-" + self._entityId + ".html")
        self.__removeFile(htmlFilePath)
        #
        alignLength = len(self._seqAlignList)
        if alignLength < 1:
            return
        #
        found_non_ala_residue = False
        found_ala_residue = False
        for alignTupList in self._seqAlignList:
            for alignTup in alignTupList:
                if (alignTup[1] == self._gapSymbol) or (alignTup[1].upper() == "ALA"):
                    if alignTup[1].upper() == "ALA":
                        found_ala_residue = True
                    #
                    continue
                #
                found_non_ala_residue = True
                break
            #
            if found_non_ala_residue:
                break
            #
        #
        all_ala_flag = False
        if (not found_non_ala_residue) and found_ala_residue:
            all_ala_flag = True
        #
        fp = open(htmlFilePath, "w")
        #
        alignmentLine = 0
        ibeg = 0
        resPerLine = 60
        self.__alignmentBlock = []
        #
        cssClassPickable = "pickable"
        cssClassTerminalResidue = "trmnlrsdue"
        #
        # contructing reusable legend markup fragment
        legendL = []
        legendL.append('<ul class="legend">\n<li>|</li>')
        for idx in range(2, resPerLine + 1):
            if idx % 10 == 0:
                legendL.append("<li>|</li>")
            elif idx % 5 == 0:
                legendL.append("<li>+</li>")
            else:
                legendL.append("<li>-</li>")
        legendL.append("</ul>\n")
        legend = "".join(legendL)
        #
        while ibeg < alignLength:
            iend = min(ibeg + resPerLine, alignLength)
            #
            # Alternate background style (odd lines)
            if alignmentLine % 2:
                cssClassBg = "greybg"
            else:
                cssClassBg = "whitebg"
            #
            fp.write('<div id="AL%s" class="%s">\n' % (str(alignmentLine), cssClassBg))
            self.__alignmentBlock.append(str(alignmentLine))
            #
            for (alignIdx, seqType, seqInstId, seqPartId, seqAltId, seqVersion, label, polymerTypeCode) in self.__orderAlignIndexList:
                idL = seqType + "_" + seqInstId + "_" + str(alignmentLine)
                fp.write('<div id="%s" class="%s">%s &nbsp;</div>\n' % (idL, seqType, label))
                #
                cssT = cssClassPickable + " " + cssClassBg
                fp.write('<ul id="%s" class="%s">\n' % (idL, cssT))
                #
                for sPos in range(ibeg, iend):
                    # id contains type(ref,auth,coordinate) + chain_id + 3-letter-code +
                    # orginal residue label index + position in alignment  + position in sequence
                    # type + '_' + chainId + '_' + compId + '_'+ resNum + '_' + str(sPos) + '_' + seqIdx
                    #
                    idS = self._getResLabelId(
                        seqType=seqType,
                        seqInstId=seqInstId,
                        seqAltId=seqAltId,
                        seqVersion=seqVersion,
                        residueCode3=self._seqAlignList[sPos][alignIdx][1],
                        residueLabelIndex=self._seqAlignList[sPos][alignIdx][2],
                        alignIndex=sPos,
                        seqIndex=codeSeqIndex(self._seqAlignList[sPos][alignIdx][3]),
                        residueType=polymerTypeCode,
                        seqPartId=seqPartId,
                    )
                    #
                    cssClassType = ""
                    if polymerTypeCode == "RNA":
                        cssClassType = " bgcolrna "
                    elif polymerTypeCode == "DNA":
                        cssClassType = " bgcoldna "
                    elif polymerTypeCode == "XNA":
                        cssClassType = " bgcolrna "
                        if self._srd.isDNA(self._seqAlignList[sPos][alignIdx][1]):
                            cssClassType = " bgcoldna "
                        #
                    #
                    if (seqType == "xyz") or (seqType == "ref"):
                        cssClassEditable = "draggable"
                    elif (seqType == "ref") and (self._seqAlignList[sPos][alignIdx][1] != self._gapSymbol):
                        cssClassEditable = ""
                    else:
                        cssClassEditable = "dblclick draggable"
                    #
                    if seqType == "ref" and (sPos == 0 or sPos == (alignLength - 1)):
                        cssClassEditable += " " + cssClassTerminalResidue
                    #
                    if self._seqAlignList[sPos][alignIdx][6] != 0:
                        cssPosClassBg = (
                            cssClassEditable
                            + " "
                            + self.__assignConflictCssStyle(
                                (self._seqAlignList[sPos][alignIdx][6], self._seqAlignList[sPos][alignIdx][7]), seqType, self._seqAlignList[sPos][alignIdx][5]
                            )
                        )
                    else:
                        cssPosClassBg = cssClassEditable + cssClassType
                    #
                    # Add highlighting of exceptional features with the coordinate sequence.
                    #
                    if (
                        (seqType == "xyz")
                        and (len(self._seqAlignList[sPos][alignIdx][5]) > 0)
                        and (not (cssPosClassBg.find("cf-misc-") != -1))
                        and (not (cssPosClassBg.find("cf-ala-gly") != -1))
                        and (not (cssPosClassBg.find("cf-glu-gln") != -1))
                        and (not (cssPosClassBg.find("cf-asp-asn") != -1))
                    ):
                        if self._seqAlignList[sPos][alignIdx][5].find("link") != -1:
                            cssPosClassBg += " bgxyzlink "
                        #
                        if self._seqAlignList[sPos][alignIdx][5].find("disorder") != -1:
                            cssPosClassBg += " bgxyzdisorder "
                        #
                        if self._seqAlignList[sPos][alignIdx][5].find("occ") != -1:
                            cssPosClassBg += " bgxyzocc  "
                        #
                    #
                    if (seqType == "xyz") and (len(self._seqAlignList[sPos][alignIdx][5]) > 0) and (self._seqAlignList[sPos][alignIdx][5].find("hetero") != -1):
                        cssPosClassBg += " bgxyzhetero  "
                        ii = self._seqAlignList[sPos][alignIdx][5].find("hetero")
                        idS += "_" + self._seqAlignList[sPos][alignIdx][5][ii + 7 :]
                    #
                    if all_ala_flag and (seqType == "auth"):
                        cssPosClassBg += " cf-ala-unk "
                    #
                    # Add highlighting for linker regions in the sample sequence --
                    #
                    if (seqType == "auth") and (len(self._seqAlignList[sPos][alignIdx][5]) > 0):
                        if self._seqAlignList[sPos][alignIdx][5].find("linker") != -1:
                            cssPosClassBg += " cf-misc-ref "
                            self._seqAlignList[sPos][alignIdx][6] = 1
                        #
                    #
                    if self._seqAlignList[sPos][alignIdx][1] == self._gapSymbol:
                        cssPosClassBg += " cf-gap-test "
                    #
                    # for tool tips -
                    cssPosClassBg += " bt "
                    #
                    fp.write('<li id="%s" class="%s">%s</li>' % (idS, cssPosClassBg, self._seqAlignList[sPos][alignIdx][0]))
                #
                fp.write('</ul>\n<div class="clearfloat"></div>\n')
            #
            # output legend line -
            fp.write('<div class="legendcount">%s</div>\n' % str(ibeg + 1))
            fp.write(legend)
            fp.write('<div class="clearfloat"></div>\n')
            fp.write("<br />\n")
            #
            alignmentLine += 1
            ibeg = alignmentLine * resPerLine
            #
            fp.write("</div>\n")
        #
        fp.close()

    def __renderConflictTable(self, authIdx):
        """Render the conflict table from the alignment"""
        htmlFilePath = os.path.join(self._sessionPath, "conflict-report-" + self._entityId + ".html")
        self.__removeFile(htmlFilePath)
        #
        conflictList = self.__getConflictTable(authIdx)
        if not conflictList:
            return
        #
        fp = open(htmlFilePath, "w")
        fp.write('<table id="conflicttable">\n')
        fp.write("<thead>\n")
        fp.write("<tr>\n")
        fp.write("<th>%s</th>" % ("Alignment<br />Position"))
        fp.write("<th>%s</th>" % ("Auth<br />Entity:" + self._entityId))
        fp.write("<th>%s</th>" % ("Aligned<br />Sequence"))
        fp.write("<th>%s</th>" % ("Residue"))
        fp.write("<th>%s</th>" % ("Annotation<br />Details"))
        fp.write("</tr>\n")
        fp.write("</thead>\n")
        fp.write("<tbody>\n")
        fp.write("%s" % "".join(conflictList))
        fp.write("</tbody>\n</table>\n")
        fp.close()

    def __renderAnnotations(self):
        """Get selected annotations for display in alignment view display"""
        htmlFilePath = os.path.join(self._sessionPath, "current-alignment-annotation-" + self._entityId + ".html")
        self.__removeFile(htmlFilePath)
        #
        if len(self._seqAlignList) < 1:
            return
        #
        featureObj = self._getFeatureObjByPackLabelFromDataStore(self._authLabel)
        #
        self.__annotationDict = {}
        self.__annotationDict["mutation"] = featureObj.getEntityMutationDetails()
        self.__annotationDict["mutationOrig"] = featureObj.getEntityMutationDetailsOrig()
        self.__annotationDict["description"] = featureObj.getEntityDescription()
        self.__annotationDict["descriptionOrig"] = featureObj.getEntityDescriptionOrig()
        self.__annotationDict["details"] = featureObj.getEntityDetails()
        self.__annotationDict["detailsOrig"] = featureObj.getEntityDetailsOrig()
        self.__annotationDict["source_method"] = featureObj.getEntityDescriptionOrig()
        #
        fp = open(htmlFilePath, "w")
        fp.write("<table>")
        fp.write("<tr><th>&nbsp;</th><th>Current</th><th>Author Provided</th></tr>")
        fp.write("<tr>")
        fp.write("<td><b>Name</b></td>")
        fp.write("<td> %s </td> " % self.__annotationDict["description"])
        fp.write("<td> %s </td> " % self.__annotationDict["descriptionOrig"])
        fp.write("</tr><tr>")
        fp.write("<td><b>Mutation</b></td>")
        fp.write("<td>%s</td>" % self.__annotationDict["mutation"])
        fp.write("<td>%s</td>" % self.__annotationDict["mutationOrig"])
        fp.write("</tr><tr>")
        fp.write("<td><b>Entity details</b></td>")
        fp.write("<td>%s</td>" % self.__annotationDict["details"])
        fp.write("<td>%s</td>" % self.__annotationDict["detailsOrig"])
        fp.write("</tr>")
        fp.write("</table>")
        fp.close()

    def __renderWarning(self, authIdx):
        """Render engineered mutation warning for nature source sequence"""
        #       htmlFilePath = os.path.join(self._sessionPath, "warning-alignment-" + self._entityId + ".html")
        #       self.__removeFile(htmlFilePath)
        #
        if len(self._seqAlignList) < 1:
            return
        #
        if self.__annotationDict["source_method"].upper() != "NAT":
            return
        #
        eelCommentL = []
        for alignTup in self._seqAlignList:
            if len(alignTup[authIdx][5]) < 1:
                continue
            #
            _aType, comment = self._decodeComment(alignTup[authIdx][5])
            for cType in ("engineered mutation", "expression tag", "linker"):
                if comment.find(cType) != -1:
                    eelCommentL.append(cType)
                #
            #
        #
        self.__warningMsg += self._getNaturalSourceWarningMessage(self.__annotationDict["source_method"].upper(), eelCommentL)

    #       if not self.__warningMsg:
    #           return
    #       #
    #       fp = open(htmlFilePath, "w")
    #       fp.write("%s\n" % self.__warningMsg)
    #       fp.close()

    # def __writeMisMatchTypeText(self):
    #     """Write over all mismatch type into the text file"""
    #     textFilePath = os.path.join(self._sessionPath, "mismatch-type.txt")
    #     fp = open(textFilePath, "w")
    #     if not self.__misMatchTypes:
    #         fp.write("%s" % "no-mismatch")
    #     elif len(self.__misMatchTypes) == 1:
    #         fp.write("%s" % self.__misMatchTypes[0])
    #     else:
    #         fp.write("%s" % "cf-all-mismatch")
    #     #
    #     fp.close()

    def __updateSeqCoodConflict(self, numConflict):
        """ """
        warningD = {}
        picklePath = self._pI.getFilePath(self.__identifier, contentType="mismatch-warning", formatType="pic", fileSource="session")
        if os.access(picklePath, os.F_OK):
            fb = open(picklePath, "rb")
            warningD = pickle.load(fb)
            fb.close()
        #
        misMatchList = []
        if "mismatch" in warningD:
            misMatchList = warningD["mismatch"]
        #
        if numConflict > 0:
            misMatchList.append(self._entityId)
        else:
            tmpList = []
            for entityId in misMatchList:
                if entityId == self._entityId:
                    continue
                #
                tmpList.append(entityId)
            #
            misMatchList = tmpList
        #
        if len(misMatchList) > 0:
            sortedEntityIdList = []
            try:
                intEntityIdList = []
                for entityId in misMatchList:
                    if int(entityId) in intEntityIdList:
                        continue
                    #
                    intEntityIdList.append(int(entityId))
                #
                intEntityIdList.sort()
                for iEntityId in intEntityIdList:
                    sortedEntityIdList.append(str(iEntityId))
                #
            except:  # noqa: E722 pylint: disable=bare-except
                sortedEntityIdList = sorted(set(misMatchList))
            #
            misMatchList = sortedEntityIdList
        #
        warningD["mismatch"] = misMatchList
        fb = open(picklePath, "wb")
        pickle.dump(warningD, fb)
        fb.close()

    def __updateMiscDataStore(self):
        """Write self.__alignIdList, self.__selectedIdList, self.__warningMsg & self.__misMatchTypes values to pickle file"""
        pickleFilePath = os.path.join(self._sessionPath, "alignment-" + self._entityId + "-misc.pic")
        self.__removeFile(pickleFilePath)
        #
        myD = {}
        myD["alignids"] = ",".join(self.__alignIdList)
        if len(self._selfRefPartIdList) > 0:
            selectedIdList = copy.deepcopy(self.__selectedIdList)
            selectedIdList.extend(self._selfRefPartIdList)
            myD["selectids"] = ",".join(selectedIdList)
        else:
            myD["selectids"] = ",".join(self.__selectedIdList)
        #
        myD["alignmentblocklist"] = ",".join(self.__alignmentBlock)
        myD["missingauthseqmap"] = self.__missingAuthSeqMap
        myD["blockedithtml"] = self.__getBlockEditFormHtml()
        repdelhtml = self.__getRepopulateDeletionFormHtml()
        if repdelhtml:
            myD["repdelhtml"] = repdelhtml
        #
        if not self.__misMatchTypes:
            myD["gedittype"] = "no-mismatch"
        elif len(self.__misMatchTypes) == 1:
            myD["gedittype"] = self.__misMatchTypes[0]
        else:
            myD["gedittype"] = "cf-all-mismatch"
        #
        if self.__warningMsg:
            myD["warning"] = self.__warningMsg
        #
        fb = open(pickleFilePath, "wb")
        pickle.dump(myD, fb)
        fb.close()

    def __updateUtilDataStore(self):
        """Update UtilDataStore object"""
        eD = {}
        for key in ("STRUCT_TITLE", "PDB_ID"):
            eD[key] = self.getEntryDetail(key)
        #
        uds = UtilDataStore(reqObj=self._reqObj, verbose=self._verbose, log=self._lfh)
        uds.set("title", eD["STRUCT_TITLE"])
        uds.set("pdbid", eD["PDB_ID"])
        uds.set("rev-allalignids", self.__alignIdList)
        if len(self._selfRefPartIdList) > 0:
            selectedIdList = copy.deepcopy(self.__selectedIdList)
            selectedIdList.extend(self._selfRefPartIdList)
            uds.set("rev-selectids", selectedIdList)
        else:
            uds.set("rev-selectids", self.__selectedIdList)
        #
        uds.serialize()

    def __removeFile(self, filePath):
        """remove existing file"""
        if os.access(filePath, os.F_OK):
            os.remove(filePath)
        #

    def __assignConflictCssStyle(self, conflictTup, seqType, comment):
        """Assign a CSS Sytle for the input conflict -"""
        styleS = ""
        if conflictTup[0] == 0:
            styleS = ""
        elif conflictTup[0] == 1:
            styleS = " cf-misc-ref "
        elif conflictTup[0] == 2:
            styleS = " cf-misc-test "
            if (conflictTup[1] == "UNK") or ((len(comment) > 0) and self._isCompatible(comment, conflictTup[1])):
                styleS = " cf-misc-test cf-other-mismatch cf-all-mismatch "
                if "cf-other-mismatch" not in self.__misMatchTypes:
                    self.__misMatchTypes.append("cf-other-mismatch")
                #
            #
        elif conflictTup[0] == 3:
            styleS = " cf-gap-test "
        elif conflictTup[0] == 4:
            styleS = " cf-gap-ref "
        elif conflictTup[0] == 5:
            styleS = " cf-glu-gln cf-all-mismatch "
            if "cf-glu-gln" not in self.__misMatchTypes:
                self.__misMatchTypes.append("cf-glu-gln")
            #
        elif conflictTup[0] == 6:
            styleS = " cf-asp-asn cf-all-mismatch "
            if "cf-asp-asn" not in self.__misMatchTypes:
                self.__misMatchTypes.append("cf-asp-asn")
            #
        elif conflictTup[0] == 7 and seqType == "xyz":
            styleS = " cf-ala-gly cf-all-mismatch "
            if "cf-ala-gly" not in self.__misMatchTypes:
                self.__misMatchTypes.append("cf-ala-gly")
            #
        elif conflictTup[0] == 8:
            styleS = " cf-met-mse cf-misc-ref cf-all-mismatch "
            if "cf-met-mse" not in self.__misMatchTypes:
                self.__misMatchTypes.append("cf-met-mse")
            #
        elif conflictTup[0] == 9:
            if len(comment) > 0:
                styleS = " cf-met-mse cf-misc-ref cf-all-mismatch "
            else:
                styleS = " cf-met-mse "
            #
        elif conflictTup[0] == 10:
            styleS = " cf-misc-test "
        else:
            styleS = ""
        #
        # JDW JDW - generalize editing
        if (len(styleS) > 0) and (seqType != "ref"):
            styleS += " cf-rep-" + conflictTup[1]
        #
        return styleS

    def __getConflictTable(self, authIdx):
        """Get conflict table content from the alignment

        alignTup[0: 1_L_code, 1: 3_L_code, 2: numbering, 3: sequence index, 4: align position, 5: comment, 6: conflict number, 7: conflict comment ]
        """
        conflictList = []
        # seqPosAuth = 0
        expressionTagCount = 0
        foundLongExpressionTag = False
        #
        for alPos, alignTup in enumerate(self._seqAlignList):
            if len(alignTup[authIdx][5]) < 1:
                if expressionTagCount >= self._longExpressionTagCountCutoff:
                    foundLongExpressionTag = True
                #
                expressionTagCount = 0
                continue
            #
            _cType, comment = self._decodeComment(alignTup[authIdx][5])
            if comment == "expression tag":
                expressionTagCount += 1
            else:
                if expressionTagCount >= self._longExpressionTagCountCutoff:
                    foundLongExpressionTag = True
                #
                expressionTagCount = 0
            #
            idAuth = self._getResLabelId(
                seqType=self.__seqTypeAuth,
                seqInstId=self.__seqInstIdAuth,
                seqAltId=self.__seqAltIdAuth,
                seqVersion=self.__seqVersionAuth,
                residueCode3=alignTup[authIdx][1],
                residueLabelIndex=alignTup[authIdx][2],
                alignIndex=alPos,
                seqIndex=codeSeqIndex(alignTup[authIdx][3]),
                residueType=self.__polymerTypeCodeAuth,
                seqPartId=self.__seqPartIdAuth,
            )
            #
            inPart = False
            for (alignIdx, seqType, seqInstId, seqPartId, seqAltId, seqVersion, label, polymerTypeCode, begPos, endPos) in self.__orderAlignRefIndexList:
                if (alPos < begPos) or (alPos > endPos):
                    continue
                #
                inPart = True
                #
                if alignTup[authIdx][1] != alignTup[alignIdx][1]:
                    if alignTup[authIdx][1] == self._gapSymbol:
                        if alignTup[alignIdx][0] == "X":
                            self.__missingAuthSeqMap[str(alPos)] = "(" + alignTup[alignIdx][1] + ")"
                        else:
                            self.__missingAuthSeqMap[str(alPos)] = alignTup[alignIdx][0]
                        #
                    #
                    idRef = self._getResLabelId(
                        seqType=seqType,
                        seqInstId=seqInstId,
                        seqAltId=seqAltId,
                        seqVersion=seqVersion,
                        residueCode3=alignTup[alignIdx][1],
                        residueLabelIndex=alignTup[alignIdx][2],
                        alignIndex=alPos,
                        seqIndex=codeSeqIndex(alignTup[alignIdx][3]),
                        residueType=polymerTypeCode,
                        seqPartId=seqPartId,
                    )
                    conflictList.extend(
                        self.__wrConflictTableRow(
                            alPos, idAuth, alignTup[authIdx][1], idRef, label, seqType, alignTup[alignIdx][1], (alignTup[alignIdx][6], alignTup[alignIdx][7]), comment
                        )
                    )
                #
            #
            if not inPart:
                conflictList.extend(self.__wrConflictTableRowOther(alPos, idAuth, alignTup[authIdx][1], comment))
            #
        #
        if expressionTagCount >= self._longExpressionTagCountCutoff:
            foundLongExpressionTag = True
        #
        if foundLongExpressionTag:
            self.__warningMsg += "Entity '" + self._entityId + "' has a long expression tag.<br />"
        #
        return conflictList

    def __wrConflictTableRow(self, sPos, idRef, compIdRef, idTest, labelTest, seqTypeTest, compIdTest, cLTest, details):
        """ """
        cssEditDetails = "dblclickselect"
        cL = []
        cL.append("<tr><td>%d</td>" % (int(sPos) + 1))
        cL.append('<td><span id="%s_R" class="%s">%s</span></td>' % (idRef, "", compIdRef))
        #
        cL.append('<td><span class="nowrap">%s</span></td>' % labelTest)
        cL.append("<td>")
        #
        cssConflict = self.__assignConflictCssStyle(cLTest, seqTypeTest, details)
        cL.append('<span id="%s_R" class="%s nowrap">%s</span>' % (idTest, cssConflict, compIdTest))
        cL.append("</td>")
        cL.append('<td><span id="%s_D" class="%s">%s</span></td>' % (idRef, cssEditDetails, details))
        cL.append("</tr>\n")
        return cL

    def __wrConflictTableRowOther(self, sPos, idRef, compIdRef, details):
        """ """
        cssEditDetails = "dblclickselect"
        cL = []
        cL.append("<tr><td>%d</td>" % (int(sPos) + 1))
        cL.append('<td><span id="%s_R" class="%s">%s</span></td>' % (idRef, " ", compIdRef))
        #
        cL.append("<td>%s</td>" % " ")
        cL.append("<td>")
        #
        cssConflict = " "
        cL.append('<span id="%s_R" class="%s">%s</span>' % (idRef, cssConflict, " "))
        cL.append("</td>")
        cL.append('<td><span id="%s_D" class="%s">%s</span></td>' % (idRef, cssEditDetails, details))
        cL.append("</tr>\n")
        return cL

    def __getBlockEditFormHtml(self):
        """ """
        form_template = """
            <input type="hidden" id="chainids" name="chainids" value="%(chainids)s" />
            <input type="hidden" name="identifier" value="%(identifier)s" />
            <input type="hidden" name="aligntag" value="%(entityid)s" />
            <input type="hidden" name="sessionid" value="%(sessionid)s" />
            <table>
            <tr><th>XYZ Chain ID</th><th>Start Position</th><th>End Posotion</th><th>Move to Position</th></tr>
            %(table_row)s
            </table>
        """
        #
        myD = {}
        myD["chainids"] = ",".join(self.__pdbChainIdList)
        myD["identifier"] = self._dataSetId
        myD["entityid"] = self._entityId
        myD["sessionid"] = self.__sessionId
        myD["table_row"] = ""
        for chainID in self.__pdbChainIdList:
            myD["table_row"] += "<tr><td>" + chainID + "</td>"
            myD["table_row"] += '<td><input type="text" id="start_position_' + chainID + '" name="start_position_' + chainID + '" value="" + size="20" /></td>'
            myD["table_row"] += '<td><input type="text" id="end_position_' + chainID + '" name="end_position_' + chainID + '" value="" + size="20" /></td>'
            myD["table_row"] += '<td><input type="text" id="move_to_' + chainID + '" name="move_to_' + chainID + '" value="" + size="20" /></td></tr>\n'
        #
        return form_template % myD

    def __getRepopulateDeletionFormHtml(self):
        """ """
        if (self.__refAlignIdx < 0) or (self.__authAlignIdx < 0):
            return ""
        #
        deleteList = []
        start = -1
        end = -1
        for i in range(0, len(self._seqAlignList)):
            if self._seqAlignList[i][self.__authAlignIdx][1] != self._gapSymbol:
                if (start >= 0) and (end >= 0):
                    deleteList.append((start, end))
                #
                start = -1
                end = -1
            elif self._seqAlignList[i][self.__refAlignIdx][1] != self._gapSymbol:
                if start < 0:
                    start = i + 1
                #
                end = i + 1
            #
        #
        if (start >= 0) and (end >= 0):
            deleteList.append((start, end))
        #
        if not deleteList:
            return ""
        #
        form_template = """
            <input type="hidden" id="repblocknum" name="repblocknum" value="%(repblocknum)s" />
            <input type="hidden" name="identifier" value="%(identifier)s" />
            <input type="hidden" name="aligntag" value="%(entityid)s" />
            <input type="hidden" name="sessionid" value="%(sessionid)s" />
            <table>
            <tr><th>Block #</th><th>Start Position</th><th>End Posotion</th><th>Action</th></tr>
            %(table_row)s
            </table>
        """
        #
        myD = {}
        myD["repblocknum"] = str(len(deleteList))
        myD["identifier"] = self._dataSetId
        myD["entityid"] = self._entityId
        myD["sessionid"] = self.__sessionId
        myD["table_row"] = ""
        #
        row = 0
        for delTup in deleteList:
            myD["table_row"] += '<input type="hidden" id="actual_repstartposition_' + str(row) + '" value="' + str(delTup[0]) + '" />'
            myD["table_row"] += '<input type="hidden" id="actual_rependposition_' + str(row) + '" value="' + str(delTup[1]) + '" />'
            #
            myD["table_row"] += "<tr><td>" + str(row + 1) + "</td>"
            myD["table_row"] += (
                '<td><input type="text" id="repstartposition_' + str(row) + '" name="repstartposition_' + str(row) + '" value="' + str(delTup[0]) + '" size="20" /></td>'
            )
            myD["table_row"] += (
                '<td><input type="text" id="rependposition_' + str(row) + '" name="rependposition_' + str(row) + '" value="' + str(delTup[1]) + '" size="20" /></td>'
            )
            myD["table_row"] += (
                '<td><input type="button" id="repsubmit_' + str(row) + '" name="repsubmit_' + str(row) + '" class="submit_populatefrm" value="Update" /></td></tr>\n'
            )
            #
            row += 1
        #
        myD["table_row"] += (
            '<tr><td style="border-style:none"></td><td style="border-style:none"></td><td style="border-style:none"></td>'
            + '<td style="border-style:none"><input type="button" id="repsubmit_all" name="repsubmit_all" value="Update All"'
            + ' class="submit_populatefrm" /></td></tr>\n'
        )
        return form_template % myD
