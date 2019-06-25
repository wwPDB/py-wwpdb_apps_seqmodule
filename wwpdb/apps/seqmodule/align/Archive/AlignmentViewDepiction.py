##
# File:    AlignmentViewDepiction.py
# Date:    20-Apr-2010
#
# Updates:
# 2010-04-20 jdw Ported to module seqmodule.
# 2010-04-20 jdw Separated the rendering modules from AlignView.
# 2010-07-25 RPS Modified so that DB reference sequence types are immutable except for deletions of terminal ranges
# 2010-07-27 RPS Added support for accommodating different ordering of sequence types as per user preferences
# 2010-08-10 RPS Predefined global edits to be effective only for non-DBref sequences
# 2010-08-11 RPS Updated __renderConflictTableOrg so that no longer itemizing gap conflicts
# 2010-10-19 RPS Updated __renderConflictTableOrg so that user is allowed to input annotation details for auth vs. ref discrepancies
# 2013-03-10 jdw Added new rendering of coordinate residues with notable features
# 2013-03-30 jdw Adust css for drag & drop operations around gap regions
# 2013-05-21 jdw fix sequence indexing in conflict table
# 2013-06-12 jdw change zero based to one-based index for filtering sequence parts in render conflict table method
# 2013-06-13 jdw add new display class "bt" to enable tool tips for all sequence types in alignment view.
# 2014-02-14 jdw fix treatment of existing conflicts in default conflict annotation.
#                Add conflict style to initiating methionine  -
# 2014-02-14 jdw adjust rendering of methionine.
# 2014-05-14 jdw adjust conflict table to provide more details for 'deletion'
# 2017-09-08 zf  treat 'cf-ala-gly', 'cf-glu-gln', 'cf-asp-asn' coloring same as 'cf-misc-ref' & 'cf-misc-test'
##
"""
Controlling class for the production of sequence alignment views.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.14"

import sys
import copy


from wwpdb.apps.seqmodule.util.SequenceLabel import ResidueLabel, SequenceFeature
from wwpdb.apps.seqmodule.util.SequenceReferenceData import SequenceReferenceData

# for part boundaries
from wwpdb.apps.seqmodule.align.AlignmentUtils import AlignmentUtils,IsCompatible


class AlignmentViewDepiction(object):

    """ Render the input alignment view and conflict table.

    """

    def __init__(self, reqObj, verbose=False, log=sys.stderr):
        self.__reqObj = reqObj
        self.__verbose = verbose
        self.__debug = True
        self.__lfh = log
        self.__alignGroupId = None

        self.__srd = SequenceReferenceData(verbose=self.__verbose, log=self.__lfh)
        self.__gapSymbol = self.__srd.getGapSymbol()
        self.__sds = None
        self.__misMatchTypes = []

    def renderAnnotations(self, aD):
        """ Get selected annotations for display in alignment view display --
        """
        description = aD['description']
        descriptionOrig = aD['descriptionOrig']
        oL = []
        oL.append('<table>')
        oL.append('<tr><th>&nbsp;</th><th>Current</th><th>Author Provided</th></tr>')
        oL.append('<tr>')
        oL.append('<td><b>Name</b></td>')
        oL.append('<td> %s </td> ' % description)
        oL.append('<td> %s </td> ' % descriptionOrig)
        oL.append('</tr><tr>')
        oL.append('<td><b>Mutation</b></td>')
        mutation = aD['mutation']
        mutationOrig = aD['mutationOrig']
        oL.append('<td>%s</td>' % mutation)
        oL.append('<td>%s</td>' % mutationOrig)
        oL.append('</tr>')
        oL.append('</table>')
        return oL

    def renderAlignment(self, alignSeqList=[], alignGroupId='1', type='original', viewOrderCode='auth-xyz-ref'):
        """Render the input alignment list in HTML according the selected option type.
        """
        self.__alignGroupId = alignGroupId
        if self.__verbose:
            self.__lfh.write("+AlignViewDepiction.renderAlignment() starting for align group %s\n" % self.__alignGroupId)

        if (type == "original"):
            return self.__renderAlignment(alignSeqList, viewOrderCode)
        else:
            return self.__renderAlignment(alignSeqList, viewOrderCode)

    def renderConflictTable(self, alignSeqList=[], alignGroupId='1', type='original'):
        """Render the conflict table from the input alignment list in HTML according the selected option type.
        """
        self.__alignGroupId = alignGroupId
        if (type == "original"):
            return self.__renderConflictTable(alignSeqList)
        else:
            return self.__renderConflictTable(alignSeqList)

    def getMisMatchType(self):
        if not self.__misMatchTypes:
            return 'no-mismatch'
        elif len(self.__misMatchTypes) == 1:
            return self.__misMatchTypes[0]
        else:
            return 'cf-all-mismatch'
        #

    def __generateSeqListForDisplay(self, alignSeqList, viewOrderCode):
        """ RPS, 2010-07-27: From input alignSeqList generate a local re-ordered copy
            that satisfies the user preferences for display order
        """

        returnList = []
        viewOrderL = viewOrderCode.split("-")

        for sqncType in viewOrderL:

            for aTup in alignSeqList:

                seqLabel = aTup[1]
                (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = seqLabel.get()

                if (seqType == sqncType):
                    returnList.append(aTup)

        return returnList

    def __renderAlignment(self, alignSeqList, viewOrderCode):
        """ Render alignment list in HTML.  (orignal version)

            Returns:

            List of lines of HTML markup.

        """

        # verify length consistency -
        #
        alignLength = len(alignSeqList[0][2])
        #
        #
        # obtain locally reordered version of alignSeqlist to meet
        # user preferences for display order
        alignSeqDisplayList = self.__generateSeqListForDisplay(alignSeqList, viewOrderCode)
        #
        resLabel = ResidueLabel()
        alignmentLine = 0
        ibeg = 0
        resPerLine = 60
        cssClassPickable = "pickable"
        cssClassSeqLabel = "seqLabel"
        cssClassTerminalResidue = "trmnlrsdue"

        # contructing reusable legend markup fragment
        legendL = []
        legendL.append('<ul class="legend">\n<li>|</li>')
        for idx in range(2, resPerLine + 1):
            if (idx % 10 == 0):
                legendL.append('<li>|</li>')
            elif (idx % 5 == 0):
                legendL.append('<li>+</li>')
            else:
                legendL.append('<li>-</li>')
        legendL.append('</ul>\n')
        # ---
        legend = "".join(legendL)
        #
        outL = []
        #
        while (ibeg < alignLength):
            iend = min(ibeg + resPerLine, alignLength)

            # Alternate background style (odd lines)
            if (alignmentLine % 2):
                cssClassBg = "greybg"
            else:
                cssClassBg = "whitebg"

            outL.append('<div id="AL%s" class="%s">\n' % (str(alignmentLine), cssClassBg))

            for aTup in alignSeqDisplayList:
                seqLabel = aTup[1]
                # seqLabel.printIt(self.__lfh)
                (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = seqLabel.get()

                fD = aTup[4]
                if (seqType == "ref"):  # aD['LABEL']
                    #label = fD['DB_NAME'] + ":" + fD['DB_ACCESSION'] + " (" + "R" + str(seqAltId) +  ",V" + str(seqVersion) + ")"
                    displayCode = fD['DB_ACCESSION']
                    if len(fD['DB_ISOFORM']) > 0:
                        displayCode = fD['DB_ISOFORM']

                    label = fD['DB_NAME'] + ":" + displayCode
                elif (seqType == "auth"):
                    label = seqType.upper() + " Entity: " + self.__alignGroupId + " V:" + str(seqVersion)
                else:
                    label = seqType.upper() + " Chain: " + seqInstId + " V:" + str(seqVersion)
                type = seqType        # aD['TYPE']
                chainId = seqInstId      # aD['CHAIN_ID']
                aL = aTup[2]  # ['ALIGN_LIST']
                cL = aTup[3]
                polymerTypeCode = fD['POLYMER_TYPE']
                #
                cssClass = type
                #
                idL = type + '_' + chainId + '_' + str(alignmentLine)
                outL.append('<div id="%s" class="%s">%s &nbsp;</div>\n' % (idL, cssClass, label))
                #
                cssT = cssClassPickable + " " + cssClassBg
                outL.append('<ul id="%s" class="%s">\n' % (idL, cssT))

                for sPos in range(ibeg, iend):
                    rT = aL[sPos]

                    # id contains type(ref,auth,coordinate) + chain_id + 3-letter-code +
                    # orginal residue label index + position in alignment  + position in sequence
                    # type + '_' + chainId + '_' + rT[1] + '_'+ tT[2] + '_' + str(sPos) + '_' + tT[3]
                    #
                    resLabel.set(seqType=seqType, seqInstId=seqInstId, seqAltId=seqAltId, seqVersion=seqVersion,
                                 residueCode3=rT[1], residueLabelIndex=rT[2],
                                 alignIndex=sPos, seqIndex=rT[3], residueType=polymerTypeCode, seqPartId=seqPartId)
                    #
                    idS = resLabel.pack()
                    cssClassType = ""
                    if (polymerTypeCode == "RNA"):
                        cssClassType = " bgcolrna "
                    elif (polymerTypeCode == "DNA"):
                        cssClassType = " bgcoldna "
                    elif (polymerTypeCode == "XNA"):
                        cssClassType = " bgcolrna "
                        if self.__srd.isDNA(rT[1]):
                            cssClassType = " bgcoldna "

                    if type == "xyz":
                        cssClassEditable = "draggable"
                    elif((type == "ref") and (rT[0] != self.__gapSymbol)):
                        cssClassEditable = ""
                    else:
                        cssClassEditable = "dblclick draggable"

                    if(type == "ref" and (sPos == 0 or sPos == (alignLength - 1))):
                        cssClassEditable += " " + cssClassTerminalResidue

                    if cL[sPos][0] != 0:
                        cssPosClassBg = cssClassEditable + " " + self.__assignConflictCssStyle(cL[sPos], type, rT[5])
                    else:
                        cssPosClassBg = cssClassEditable + cssClassType

                    #
                    # Add highlighting of exceptional features with the coordinate sequence.
                    #
                    if ((type == 'xyz') and (len(rT[5]) > 0) and (not (cssPosClassBg.find("cf-misc-") != -1)) and \
                       (not (cssPosClassBg.find("cf-ala-gly") != -1)) and (not (cssPosClassBg.find("cf-glu-gln") != -1)) and \
                       (not (cssPosClassBg.find("cf-asp-asn") != -1))):
                        if rT[5].find("link") != -1:
                            cssPosClassBg += " bgxyzlink "
                        if rT[5].find("disorder") != -1:
                            cssPosClassBg += " bgxyzdisorder "
                        if rT[5].find("occ") != -1:
                            cssPosClassBg += " bgxyzocc  "

                    if ((type == 'xyz') and (len(rT[5]) > 0) and (rT[5].find("hetero") != -1)):
                        cssPosClassBg += " bgxyzhetero  "
                        ii = rT[5].find("hetero")
                        idS += '_' + rT[5][ii + 7:]

                    # Add highlighting for linker regions in the sample sequence --
                    #
                    if ((type == 'auth') and (len(rT[5]) > 0)):
                        if rT[5].find("linker") != -1:
                            cssPosClassBg += " cf-misc-ref "
                            tt = list(cL[sPos])
                            tt[0] = 1
                            cL[sPos] = tuple(tt)

                    # if (( rT[0] == self.__gapSymbol) and ((type == 'xyz') or (type == 'ref'))):
                    if (rT[0] == self.__gapSymbol):
                        cssPosClassBg += " cf-gap-test "

                    # for tool tips -
                    cssPosClassBg += " bt "

                    outL.append('<li id="%s" class="%s">%s</li>' % (idS, cssPosClassBg, rT[0]))

                outL.append('</ul>\n<div class="clearfloat"></div>\n')

            # output legend line -
            outL.append('<div class="legendcount">%s</div>\n' % str(ibeg + 1))
            outL.append(legend)
            outL.append('<div class="clearfloat"></div>\n')
            outL.append('<br />\n')

            alignmentLine += 1
            ibeg = alignmentLine * resPerLine
            #
            outL.append('</div>\n')
        #
        return outL

    def __assignConflictCssStyle(self, conflictTup, seqType, comment):
        """  Assign a CSS Sytle for the input conflict -

        """
        styleS = ''
        if (conflictTup[0] == 0):
            styleS = ''
        elif (conflictTup[0] == 1):
            styleS = " cf-misc-ref "
        elif (conflictTup[0] == 2):
            styleS = " cf-misc-test "
            if (conflictTup[1] == 'UNK') or ((len(comment) > 0) and IsCompatible(comment, conflictTup[1])):
                styleS = " cf-misc-test cf-other-mismatch cf-all-mismatch "
                if 'cf-other-mismatch' not in self.__misMatchTypes:
                    self.__misMatchTypes.append('cf-other-mismatch')
                #
            #
        elif (conflictTup[0] == 3):
            styleS = " cf-gap-test "
        elif (conflictTup[0] == 4):
            styleS = " cf-gap-ref "
        elif (conflictTup[0] == 5):
            styleS = " cf-glu-gln cf-all-mismatch "
            if 'cf-glu-gln' not in self.__misMatchTypes:
                self.__misMatchTypes.append('cf-glu-gln')
            #
        elif (conflictTup[0] == 6):
            styleS = " cf-asp-asn cf-all-mismatch "
            if 'cf-asp-asn' not in self.__misMatchTypes:
                self.__misMatchTypes.append('cf-asp-asn')
            #
        elif (conflictTup[0] == 7 and seqType == 'xyz'):
            styleS = " cf-ala-gly cf-all-mismatch "
            if 'cf-ala-gly' not in self.__misMatchTypes:
                self.__misMatchTypes.append('cf-ala-gly')
            #
        elif (conflictTup[0] == 8):
            styleS = " cf-met-mse cf-misc-ref cf-all-mismatch "
            if 'cf-met-mse' not in self.__misMatchTypes:
                self.__misMatchTypes.append('cf-met-mse')
            #
        elif (conflictTup[0] == 9):
            if len(comment) > 0:
                styleS = " cf-met-mse cf-misc-ref cf-all-mismatch "
            else:
                styleS = " cf-met-mse "
            #
        elif (conflictTup[0] == 10):
            styleS = " cf-misc-test "
        else:
            styleS = ''
        #
        # JDW JDW - generalize editing
        if (len(styleS) > 0) and (seqType != "ref"):
            styleS += " cf-rep-" + conflictTup[1]
        #
        return styleS

    def __wrConflictTableRow(self, sPos, idRef, compIdRef, idTest, labelTest, seqTypeTest, compIdTest, cLTest, details):
        """
        """
        cssEditDetails = "dblclickselect"
        cL = []
        cL.append('<tr><td>%d</td>' % (int(sPos) + 1))
        cL.append('<td><span id="%s_R" class="%s">%s</span></td>' % (idRef, "", compIdRef))
        #
        cL.append('<td><span class="nowrap">%s</span></td>' % labelTest)
        cL.append('<td>')

        if(seqTypeTest == "ref"):
            cssEdit = ""
        else:
            cssEdit = "dblclick"

        cssConflict = cssEdit + " " + self.__assignConflictCssStyle(cLTest[sPos], seqTypeTest, details)
        cL.append('<span id="%s_R" class="%s nowrap">%s</span>' % (idTest, cssConflict, compIdTest))
        cL.append('</td>')
        cL.append('<td><span id="%s_D" class="%s">%s</span></td>' % (idRef, cssEditDetails, details))
        cL.append('</tr>\n')
        return cL

    def __wrConflictTableRowOther(self, sPos, idRef, compIdRef, details):
        """
        """
        cssEditDetails = "dblclickselect"
        cL = []
        cL.append('<tr><td>%d</td>' % (int(sPos) + 1))
        cL.append('<td><span id="%s_R" class="%s">%s</span></td>' % (idRef, " ", compIdRef))
        #
        cL.append('<td>%s</td>' % " ")
        cL.append('<td>')

        cssEdit = ""
        # cssEdit="dblclick"

        #cssConflict =cssEdit +  " " + self.__assignConflictCssStyle(cLTest[sPos],seqTypeTest)
        #cssConflict =cssEdit +  " " + self.__assignConflictCssStyle(cLTest[sPos],seqTypeTest)
        cssConflict = " "
        cL.append('<span id="%s_R" class="%s">%s</span>' % (idRef, cssConflict, " "))
        cL.append('</td>')
        cL.append('<td><span id="%s_D" class="%s">%s</span></td>' % (idRef, cssEditDetails, details))
        cL.append('</tr>\n')
        if self.__verbose:
            self.__lfh.write("+AlignViewDepiction.wrConflictTableRowOther() sPos %s idRef %s compIdRef %s details %s\n" % (sPos, idRef, compIdRef, details))
        return cL

    def __renderConflictTable(self, alignSeqList):
        """ Render a summary of all conflicts using the input alignment list.

            Storage model - input data as organized in a tuple format used in the
                            internal alignment list.

            Aligned sequences have storage model -

            (one-letter-code, 3-letter-code, original label residue index, position in sequence, position in alignment, comment )

            In building the conflict table, conflict row identifier and conflict details (e.g. the comment above) are taken from
            reference sequence.

            Returns:

            List of lines of HTML markup for a skeleton table of conflict data.

            Return an empty list if no conflicts are found

            # Note - The identifiers associated with residues and details in the conflict
                     table are simply derivable from the residue identifiers used in the
                     alignment depiction.

        """
        seqFeature = SequenceFeature()
        aU = AlignmentUtils(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        #
        # Get the reference for the alignment and all associated details -
        #
        alTupRef = alignSeqList[0]
        #
        seqLabelRef = alTupRef[1]
        fDRef = alTupRef[4]
        aLRef = alTupRef[2]
        cLRef = alTupRef[3]

        (seqTypeRef, seqInstIdRef, seqPartIdRef, seqAltIdRef, seqVersionRef) = seqLabelRef.get()
        #
        if seqTypeRef == 'auth':
            authPartD = aU.getAuthPartDetails(seqLabelRef.pack())

        if (seqTypeRef == "ref"):
            labelRef = fDRef['DB_NAME'] + "<br />" + fDRef['DB_ACCESSION']
        elif (seqTypeRef in ["xyz"]):
            labelRef = seqTypeRef.upper() + "<br />Chain:" + seqInstIdRef
        elif (seqTypeRef in ["auth"]):
            labelRef = seqTypeRef.upper() + "<br />Entity:" + self.__alignGroupId
        else:
            labelRef = "SEQ ID: " + seqInstId

        #

        resLabel = ResidueLabel()
        # --------------------
        iCount = 0
        cL = []
        cL.append('<table id="conflicttable">\n')

        cL.append('<thead>\n')
        cL.append('<tr>\n')
        cL.append('<th>%s</th>' % ("Alignment<br />Position"))
        cL.append('<th>%s</th>' % (labelRef))
        cL.append('<th>%s</th>' % ("Aligned<br />Sequence"))
        cL.append('<th>%s</th>' % ("Residue"))
        cL.append('<th>%s</th>' % ("Annotation<br />Details"))
        cL.append('</tr>\n')
        cL.append('</thead>\n')
        cL.append('<tbody>\n')

        seqPosRef = 0
        for alPos in range(0, len(aLRef)):

            if cLRef[alPos][0] != 0:
                # There is some conflict with the current reference sequence -
                rTRef = aLRef[alPos]
                compIdRef = rTRef[1]
                labelIdxRef = rTRef[2]
                seqIdxRef = rTRef[3]

                if compIdRef != self.__gapSymbol:
                    seqPosRef = seqIdxRef + 1
                detailsRef = rTRef[5]
                #
                #  Add to show deletions
                if ((compIdRef == self.__gapSymbol) and (detailsRef in ['deletion', 'DELETION'])):
                    # position in alignment -- JDW --
                    seqPosRef = rTRef[4] + 1
                #
                resLabel.set(seqType=seqTypeRef, seqInstId=seqInstIdRef, seqAltId=seqAltIdRef, seqVersion=seqVersionRef,
                             residueCode3=compIdRef, residueLabelIndex=labelIdxRef,
                             alignIndex=alPos, seqIndex=seqIdxRef, seqPartId=seqPartIdRef)
                idRef = resLabel.pack()

                inPart = False

                for alTupTest in alignSeqList[1:]:
                    seqLabelTest = alTupTest[1]
                    (seqTypeTest, seqInstIdTest, seqPartIdTest, seqAltIdTest, seqVersionTest) = seqLabelTest.get()
                    #
                    # we do not annotate conflicts with coordinates here.
                    if seqTypeTest in ['XYZ', 'xyz']:
                        continue
                    #
                    fDTest = alTupTest[4]
                    seqFeature.set(fDTest)
                    partId, seqPartBeg, seqPartEnd, SeqPartType = seqFeature.getAuthPartDetails()
                    #
                    seqPartBeg = authPartD[seqPartIdTest][1]
                    seqPartEnd = authPartD[seqPartIdTest][2]

                    if (self.__verbose):
                        self.__lfh.write("+AlignViewDepiction.renderConflictTable() idTest %s seqPosRef %r seqPartBeg %r seqPartEnd %r\n" %
                                         (alTupTest[0], seqPosRef, seqPartBeg, seqPartEnd))

                    if seqPosRef < seqPartBeg or seqPosRef > seqPartEnd:
                        continue

                    inPart = True

                    #
                    aLTest = alTupTest[2]
                    rTTest = aLTest[alPos]
                    compIdTest = rTTest[1]
                    labelIdxTest = rTTest[2]
                    seqIdxTest = rTTest[3]
                    detailsTest = rTTest[5]
                    #
                    cLTest = alTupTest[3]

                    if (self.__verbose):
                        self.__lfh.write("+AlignViewDepiction.renderConflictTable() compIdRef %s compIdTest %s seqTypeRef %s seqTypeTest %s\n" %
                                         (compIdRef, compIdTest, seqTypeRef, seqTypeTest))

                    if ((compIdRef != compIdTest) and not ((compIdTest == self.__gapSymbol) and (seqTypeTest in ['XYZ', 'xyz']))):
                        #
                        if (seqTypeTest == "ref"):
                            displayCode = fDTest['DB_ACCESSION']
                            if ((fDTest['DB_NAME'] in ['UNP', 'TR', 'SP']) and (len(fDTest['DB_ISOFORM']) > 0)):
                                displayCode = fDTest['DB_ISOFORM']
                            labelTest = fDTest['DB_NAME'] + ":" + displayCode
                        else:
                            labelTest = seqTypeTest.upper() + " Chain:" + seqInstIdTest + " V:" + str(seqVersionTest)

                        resLabel.set(seqType=seqTypeTest, seqInstId=seqInstIdTest, seqAltId=seqAltIdTest, seqVersion=seqVersionTest,
                                     residueCode3=compIdTest, residueLabelIndex=labelIdxTest,
                                     alignIndex=alPos, seqIndex=seqIdxTest, seqPartId=seqPartIdTest)
                        idTest = resLabel.pack()
                        #
                        cL.extend(self.__wrConflictTableRow(alPos, idRef, compIdRef, idTest, labelTest, seqTypeTest, compIdTest, cLTest, detailsRef))
                        iCount += 1
                    elif (compIdRef != compIdTest):
                        if (self.__verbose):
                            self.__lfh.write("+AlignViewDepiction.renderConflictTable() skipped compIdRef %s compIdTest %s\n" % (compIdRef, compIfTest))

                if not inPart and (seqTypeTest not in ['XYZ', 'xyz']):
                    cL.extend(self.__wrConflictTableRowOther(alPos, idRef, compIdRef, detailsRef))
                    iCount += 1
                if ((not inPart) and (cLRef[alPos][1] == 'out-of-part') and (seqTypeRef == 'auth')):
                    cL.extend(self.__wrConflictTableRowOther(alPos, idRef, compIdRef, detailsRef))
                    iCount += 1
        #
        cL.append('</tbody>\n</table>\n')
        #
        if iCount > 0:
            return cL
        else:
            return []
