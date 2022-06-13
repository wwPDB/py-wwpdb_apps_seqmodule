##
# File:  AlignmentToolUtils.py
# Date:  27-Aug-2018
#
##
__docformat__ = "restructuredtext en"
__author__ = "Zukang Feng"
__email__ = "zfeng@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.09"

import sys


def getSeqAlignment(seqList, alignIndicesList, labelIndexDict, testSeqType, gapSymbol):  # pylint: disable=unused-argument
    """Materialize sequences alignment based sequences alignment index"""
    for alignIndex in alignIndicesList:
        if len(alignIndex) != len(seqList):
            return "getSeqAlignment: Align indices number is not same as sequence number.", []
        #
        for i in range(0, len(alignIndex)):
            intIdx, annotatedComment = decodeIndex(alignIndex[i])
            #
            if intIdx < 0:
                continue
            #
            if intIdx >= len(seqList[i]):
                return "getSeqAlignment: Align index out of range for '" + labelIndexDict[i] + "': index=" + str(intIdx) + ", sequence length=" + str(len(seqList[i])) + ".", []
            #
        #
    #
    alignList = []
    for _idx, alignIndex in enumerate(alignIndicesList):
        alignPos = len(alignList)
        alignTup = []
        for i in range(0, len(alignIndex)):
            intIdx, annotatedComment = decodeIndex(alignIndex[i])
            #
            if intIdx < 0:
                alignTup.append([gapSymbol, gapSymbol, "", "", alignPos, annotatedComment, 0, ""])
            else:
                # [ 0: 1_L_code, 1: 3_L_code, 2: numbering, 3: sequence index, 4: align position, 5: comment, 6: conflict number, 7: conflict comment ]
                comment = seqList[i][intIdx][2]
                if (annotatedComment != "") and ((comment == "") or (comment.startswith("DEFAULT:") and annotatedComment.startswith("ANNOTATED:"))):
                    comment = annotatedComment
                #
                alignTup.append([seqList[i][intIdx][4], seqList[i][intIdx][0], seqList[i][intIdx][1], intIdx, alignPos, comment, 0, ""])
            #
        #
        alignList.append(alignTup)
    #
    return "", alignList


def mergeSeqAlignment(seqList, indicesList, gapSymbol):
    """Copy all aligned sequences from seqList[0][idx], only copy seqList[i][idx][1:] sequence where i >= 1"""
    for alignIndex in indicesList:
        if len(alignIndex) != len(seqList):
            return "mergeSeqAlignment: Align indices number is not same as sequence alignment number.", []
        #
        for i in range(0, len(alignIndex)):
            intIdx, _annotatedComment = decodeIndex(alignIndex[i])
            #
            if intIdx < 0:
                continue
            #
            if intIdx >= len(seqList[i]):
                return "mergeSeqAlignment: Align index out of range: index=" + str(intIdx) + ", sequence length=" + str(len(seqList[i])) + ".", []
            #
        #
    #
    alignList = []
    for alignIndex in indicesList:
        alignPos = len(alignList)
        alignTup = []
        for i in range(0, len(alignIndex)):
            intIdx, _annotatedComment = decodeIndex(alignIndex[i])
            #
            if intIdx < 0:
                if i == 0:
                    for seqTup in seqList[i][0]:
                        alignTup.append([gapSymbol, gapSymbol, "", "", alignPos, "", 0, ""])
                    #
                else:
                    for seqTup in seqList[i][0][1:]:
                        alignTup.append([gapSymbol, gapSymbol, "", "", alignPos, "", 0, ""])
                    #
                #
            else:
                if i == 0:
                    for seqTup in seqList[i][intIdx]:
                        alignTup.append([seqTup[0], seqTup[1], seqTup[2], seqTup[3], alignPos, seqTup[5], seqTup[6], seqTup[7]])
                    #
                else:
                    for seqTup in seqList[i][intIdx][1:]:
                        alignTup.append([seqTup[0], seqTup[1], seqTup[2], seqTup[3], alignPos, seqTup[5], seqTup[6], seqTup[7]])
                    #
                #
            #
        #
        alignList.append(alignTup)
    #
    return "", alignList


def codeIndex(idx, comment):
    """'idx' is an integer index related an alignment or sequence.
    'comment' is any string comment(s) related to the property of the element at 'idx' position.
    """
    if comment != "":
        return (idx, comment)
    #
    return idx


def codeSeqIndex(index):
    """index could be a tuple like ( idx, comment ) or simply an integer or string
    return an integer if there is an positive integer or empty string for anything else
    """
    if index == "":
        return index
    #
    intIdx, _comment = decodeIndex(index)
    if intIdx < 0:
        return ""
    #
    return intIdx


def decodeIndex(index):
    """index could be a tuple like ( idx, comment ) or simply an integer or string"""
    intIdx = -1
    comment = ""
    inputType = str(type(index)).lower()
    if inputType.find("int") > 0:
        intIdx = index
    elif (inputType.find("list") > 0) or (inputType.find("tuple") > 0):
        intIdx = index[0]
        comment = index[1]
    else:
        try:
            intIdx = int(index)
        except:  # noqa: E722 pylint: disable=bare-except
            intIdx = -1
        #
    #
    return intIdx, comment


def assignConflict(refSeqType, refCompId, testSeqType, testCompId, gapSymbol):
    """On input:
              refSeqType,testSeqType = ["auth"|"ref"|"xyz"]
              refCompId,testCompId   =  component 3-letter-codes

    Returns:  isConflict, refConflictTup (code,correction), testConflictTup (code,correction)

    Methionine --  Any reference methione is candidate for met/mse -
    """
    refConflict = (0, "")
    testConflict = (0, "")
    isConflict = False

    if (testSeqType.lower() in ["xyz"]) and (testCompId == gapSymbol) and (refCompId == "MET"):
        return True, (9, "MSE"), testConflict
    #
    if (testSeqType.lower() in ["xyz"]) and (testCompId == gapSymbol):
        return isConflict, refConflict, testConflict
    #
    if (refSeqType in ["auth"]) and (testSeqType in ["ref"]) and (refCompId != testCompId) and (testCompId == gapSymbol):
        refConflict, testConflict = assignConflictType(refCompId, refSeqType, testCompId, testSeqType, gapSymbol)
        isConflict = True
    elif (refCompId != testCompId) and (refCompId == gapSymbol):
        refConflict, testConflict = assignConflictType(refCompId, refSeqType, testCompId, testSeqType, gapSymbol)
        isConflict = True
    elif (refCompId != testCompId) and not ((refCompId == gapSymbol) or (testCompId == gapSymbol)):
        refConflict, testConflict = assignConflictType(refCompId, refSeqType, testCompId, testSeqType, gapSymbol)
        isConflict = True
    #
    return isConflict, refConflict, testConflict


def assignConflictType(r3Ref, refSeqType, r3Test, testSeqType, gapSymbol):  # pylint: disable=unused-argument
    """Assign conflict type as  -

     0 - None
     1 - non-specific in ref sequence
     2 - non-specific in test sequence

     3 - r3Ref = not gap r3Test = gap
     4 - r3Ref = gap     r3Test = not gap

     5 - r3Ref = GLU     r3Test = GLN (type == "xyz")
    10 - r3Ref = GLU     r3Test = GLN (type != "xyz")

     6 - r3Ref - ASP     r3Test = ASN (type == "xyz")
    10 - r3Ref - ASP     r3Test = ASN (type != "xyz")

     7 - r3Ref = Any     r3Test = ALA/GLY (type == "xyz")
    10 - r3Ref = Any     r3Test = ALA/GLY (type != "xyz")

     8 - r3Ref = MET     r3Test = MSE
     9 - r3Ref = MET     r3Test = gap
    10 - r3Ref = MSE     r3Test = MET

     and return a tuple of conflict integer type and
     the likely correction.  ref,test   (code,r3),(code,r3)
    """
    if r3Ref == r3Test:
        return ((0, r3Ref), (0, r3Ref))
    elif r3Ref == "MET" and r3Test == gapSymbol:
        return ((9, "MSE"), (0, r3Test))
    elif r3Ref != gapSymbol and r3Test == gapSymbol:
        return ((1, r3Ref), (3, r3Ref))
    elif r3Ref == gapSymbol and r3Test != gapSymbol:
        return ((1, r3Ref), (4, r3Test))
    elif r3Ref == "GLU" and r3Test == "GLN":
        if testSeqType == "xyz":
            return ((1, r3Ref), (5, r3Ref))
        else:
            return ((1, r3Ref), (10, r3Test))
        #
    elif r3Ref == "ASP" and r3Test == "ASN":
        if testSeqType == "xyz":
            return ((1, r3Ref), (6, r3Ref))
        else:
            return ((1, r3Ref), (10, r3Test))
        #
    elif r3Test == "ALA" or r3Test == "GLY":
        if testSeqType == "xyz":
            return ((1, r3Ref), (7, r3Ref))
        else:
            return ((1, r3Ref), (10, r3Test))
        #
    elif r3Ref == "MET" and r3Test == "MSE":
        return ((8, r3Test), (1, r3Test))
    elif r3Ref == "MSE" and r3Test == "MET":
        return ((1, r3Ref), (10, r3Test))
    elif testSeqType == "xyz":
        return ((1, r3Ref), (2, r3Ref))
    #
    return ((1, r3Ref), (10, r3Test))


def printAlignment(lfh=sys.stderr, alignList=None, shortFlag=True):
    """Print alignment"""
    if alignList is None:
        alignList = []
    lfh.write("Alignment length = %d\n" % len(alignList))
    for alignTup in alignList:
        first = True
        for sTup in alignTup:
            if not first:
                lfh.write(" + ")
            #
            if shortFlag:
                lfh.write(
                    "(%s %5s %6s %6s %6s %4s %s)"
                    % (
                        "'" + sTup[0] + "'",
                        "'" + sTup[1] + "'",
                        "'" + sTup[2] + "'",
                        "'" + str(sTup[3]) + "'",
                        "'" + str(sTup[4]) + "'",
                        "'" + str(sTup[6]) + "'",
                        "'" + sTup[7] + "'",
                    )
                )
            else:
                lfh.write(
                    "(%s %5s %6s %6s %6s %35s %4s %s)"
                    % (
                        "'" + sTup[0] + "'",
                        "'" + sTup[1] + "'",
                        "'" + sTup[2] + "'",
                        "'" + str(sTup[3]) + "'",
                        "'" + str(sTup[4]) + "'",
                        "'" + str(sTup[5]) + "'",
                        "'" + str(sTup[6]) + "'",
                        "'" + sTup[7] + "'",
                    )
                )
            #
            first = False
        #
        lfh.write("\n")
    #


def printCondenseAlignment(lfh, instIdList, alignList):
    """Print alignment"""
    lfh.write("Alignment length = %d\n" % len(alignList))
    for InstId in instIdList:
        lfh.write("%-21s" % InstId)
    #
    lfh.write("\n")
    #
    for alignTup in alignList:
        first = True
        for sTup in alignTup:
            if not first:
                lfh.write(" + ")
            #
            lfh.write("(%s %5s %6s)" % ("'" + sTup[0] + "'", "'" + sTup[1] + "'", "'" + sTup[2] + "'"))
            first = False
        #
        lfh.write("\n")
    #
