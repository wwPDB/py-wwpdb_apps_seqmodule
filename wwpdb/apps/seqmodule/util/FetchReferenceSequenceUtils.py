##
# File:  FetchReferenceSequenceUtils.py
# Date:  21-Nov-2018
#
# Updates:
#  3-Oct-2022  zf  update fetchReferenceSequenceWithSeqMatch() methods for better handling author provided reference sequence cases.
#                  add runSeqAlignment() & __toList() methods
# 11-Nov-2022  zf  move __getRefInfo() method to FetchSeqInfoUtils class.
#  2-Jun-2024  zf  add persistent object self.__fetchSeqUtil.
##
"""
Methods to get reference sequence data from reference database based database name and identifier.

"""
__docformat__ = "restructuredtext en"
__author__ = "Zukang Feng"
__email__ = "zfeng@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.09"

import os
import string
import sys

from wwpdb.apps.seqmodule.io.FetchSeqInfoUtils import FetchSeqInfoUtils
from wwpdb.apps.seqmodule.util.SequenceReferenceData import SequenceReferenceData
from wwpdb.utils.align.alignlib import PseudoMultiAlign  # pylint: disable=no-name-in-module


class FetchReferenceSequenceUtils(object):
    """Fetch reference sequence data."""

    def __init__(self, siteId="WWPDB_DEPLOY_TEST", seqReferenceData=None, verbose=False, log=sys.stderr):
        self.__verbose = verbose
        self.__lfh = log
        self.__siteId = siteId
        self.__srd = seqReferenceData
        if not self.__srd:
            self.__srd = SequenceReferenceData(self.__verbose, self.__lfh)
        self.__accCode = ""
        self.__refInfoD = {}
        #
        self.__fetchSeqUtil = FetchSeqInfoUtils(siteId=self.__siteId, seqReferenceData=self.__srd, verbose=self.__verbose, log=self.__lfh)

    def fetchReferenceSequence(self, dbName, dbAccession, dbIsoform, polyTypeCode="AA", refSeqBeg=None, refSeqEnd=None):
        """return error message, sequence feature dictionary, sequence list
        """
        try:
            start = int(str(refSeqBeg))
        except:  # noqa: E722 pylint: disable=bare-except
            start = 0
        #
        try:
            end = int(str(refSeqEnd))
        except:  # noqa: E722 pylint: disable=bare-except
            end = 0
        #
        self.__accCode, self.__refInfoD = self.__fetchSeqUtil.getRefInfo(dbName, dbAccession, dbIsoform, start, end)
        #
        if (not self.__refInfoD) or ("sequence" not in self.__refInfoD):
            error = "Fetch reference sequence [ dbName=" + dbName + ", Accession=" + self.__accCode + "] failed."
            self.__accCode = ""
            self.__refInfoD = {}
            #
            return error, {}, []
        #
        # for DAOTHER-6304
        if polyTypeCode == "DNA":
            self.__refInfoD["sequence"] = self.__refInfoD["sequence"].replace("U", "T")
        #
        seqLength = len(self.__refInfoD["sequence"])
        #
        try:
            refSeqBeg = int(str(refSeqBeg))
        except:  # noqa: E722 pylint: disable=bare-except
            refSeqBeg = 1
        #
        try:
            refSeqEnd = int(str(refSeqEnd))
        except:  # noqa: E722 pylint: disable=bare-except
            refSeqEnd = seqLength
        #
        # Added reverse order seqeunce case: DAOTHER-4455
        #
        reverseOrder = False
        if refSeqEnd < refSeqBeg:
            reverseOrder = True
            refSeqBeg, refSeqEnd = refSeqEnd, refSeqBeg
        #
        errorMessage = ""
        if (refSeqBeg is not None) and (refSeqBeg < 1):
            if errorMessage:
                errorMessage += "\n"
            #
            errorMessage += "Invalid SEQ BEGIN number: " + str(refSeqBeg) + "."
        #
        if (refSeqEnd is not None) and (refSeqEnd > seqLength):
            if errorMessage:
                errorMessage += "\n"
            #
            errorMessage += "Invalid SEQ END number: " + str(refSeqEnd) + " > " + str(seqLength) + " ( sequence length of " + self.__accCode + " )."
        #
        if errorMessage:
            return errorMessage, {}, []
        #
        if refSeqBeg is None:
            refSeqBeg = 1
        #
        if refSeqEnd is None:
            refSeqEnd = seqLength  # len(sequence)
        #
        self.__refInfoD["db_length"] = seqLength
        if reverseOrder:
            self.__refInfoD["hitFrom"] = refSeqEnd
            self.__refInfoD["hitTo"] = refSeqBeg
        else:
            self.__refInfoD["hitFrom"] = refSeqBeg
            self.__refInfoD["hitTo"] = refSeqEnd
        #
        return "", self.__refInfoD, self.__getReferenceList(self.__refInfoD["sequence"], polyTypeCode, refSeqBeg, refSeqEnd, reverseOrder)

    def fetchReferenceSequenceWithSeqMatch(self, dbName, dbAccession, dbCode, taxId, authSeq, seqNumBeg, mutationList):  # pylint: disable=unused-argument
        """ """
        self.__accCode = ""
        self.__refInfoD = {}
        #
        dbIsoform = ""
        if dbName in ["UNP", "SP", "TR"]:
            tL = dbAccession.split("-")
            if len(tL) > 1:
                dbIsoform = dbAccession
                dbAccession = tL[0]
            #
        #
        self.__accCode, self.__refInfoD = self.__fetchSeqUtil.getRefInfo(dbName, dbAccession, dbIsoform, 0, 0)
        #
        if (not self.__refInfoD) or ("sequence" not in self.__refInfoD) or (not self.__refInfoD["sequence"]):
            if dbCode:
                self.__accCode, self.__refInfoD = self.__fetchSeqUtil.getRefInfo(dbName, dbCode, "", 0, 0)
                if (not self.__refInfoD) or ("sequence" not in self.__refInfoD) or (not self.__refInfoD["sequence"]):
                    return False, False, {}
                #
            else:
                return False, False, {}
            #
        #
# Removed per ticket DAOTHER-7903
#       # Skip reference sequence with different Taxonomy ID
#       if ("taxonomy_id" in self.__refInfoD) and (self.__refInfoD["taxonomy_id"]) and taxId and (self.__refInfoD["taxonomy_id"] != taxId):
#           return False, {}
#       #
        mutationMap = {}
        for val in mutationList:
            mut = val[:1] + "_" + val[-1]
            if mut in mutationMap:
                mutationMap[mut] += 1
            else:
                mutationMap[mut] = 1
            #
        #
        if ("variant" in self.__refInfoD) and self.__refInfoD["variant"]:
            for var in self.__refInfoD["variant"].split(","):
                if var not in mutationList:
                    mutationList.append(var)
                #
            #
        #
        autoMatchStatus, skipBlastSearch, alignInfoD = self.runSeqAlignment(self.__refInfoD["sequence"], authSeq, seqNumBeg, mutationList, mutationMap)
        if alignInfoD:
            # Check Taxonomy ID information
            if ("taxonomy_id" in self.__refInfoD) and (self.__refInfoD["taxonomy_id"]) and taxId and (self.__refInfoD["taxonomy_id"] != taxId):
                autoMatchStatus = False
            #
            self.__refInfoD.update(alignInfoD)
            return autoMatchStatus, skipBlastSearch, self.__refInfoD
        else:
            return False, False, {}
        #

    def runSeqAlignment(self, refSeq, authSeq, seqNumBeg, mutationList, mutationMap):
        """ Run Smith–Waterman sequence alignment algorithm with wwpdb.utils.align.alignlib.PseudoMultiAlign utility
        """
        refSeqList = self.__toList(refSeq)
        authSeqList = self.__toList(authSeq)
        #
        pA = PseudoMultiAlign()
        pA.setRefScore()
        pA.setAuthSequence(refSeqList)
        pA.addAlignSequence(authSeqList)
        alignIndexList = pA.getAlignIndices()
        #
        blockList = []
        start = -1
        end = -1
        for idx, alignIdx in enumerate(alignIndexList):
            ref = "-"  # pylint: disable=unused-variable
            if alignIdx[0] >= 0:
                ref = refSeqList[alignIdx[0]][0]  # noqa: F841
            #
            auth = "-"  # pylint: disable=unused-variable
            if alignIdx[1] >= 0:
                auth = authSeqList[alignIdx[1]][0]  # noqa: F841
            #
            if (alignIdx[0] >= 0) and (alignIdx[1] >= 0):
                if start < 0:
                    start = idx
                #
                end = idx
            else:
                # Only insert aligned block with more than one residue
                if (start >= 0) and (end > start):
                    blockList.append((start, end))
                #
                start = -1
                end = -1
            #
        #
        if (start >= 0) and (end > start):
            blockList.append((start, end))
        #
        if not blockList:
            return False, False, {}
        #
        start = blockList[0][0]
        end = blockList[-1][1]
        if (start < 0) or ((2 * (end - start + 1)) < len(authSeqList)):
            return False, False, {}
        #
        identity = 0
        mutation = 0
        gaps = 0
        midline = ""
        query = ""
        subject = ""
        for idx in range(start, end + 1):
            if (alignIndexList[idx][0] >= 0) and (alignIndexList[idx][1] >= 0):
                if authSeqList[alignIndexList[idx][1]][0] == refSeqList[alignIndexList[idx][0]][0]:
                    identity += 1
                    midline += authSeqList[alignIndexList[idx][1]][0]
                else:
                    isMutation = False
                    mut = refSeqList[alignIndexList[idx][0]][0] + str(alignIndexList[idx][0] + 1) + authSeqList[alignIndexList[idx][1]][0]
                    if mut in mutationList:
                        isMutation = True
                    #
                    mut_ = refSeqList[alignIndexList[idx][0]][0] + "_" + authSeqList[alignIndexList[idx][1]][0]
                    if mut_ in mutationMap:
                        isMutation = True
                        mutationMap[mut_] -= 1
                        if mutationMap[mut_] == 0:
                            del mutationMap[mut_]
                        #
                    #
                    if isMutation:
                        mutation += 1
                    #
                    midline += " "
                #
                query += authSeqList[alignIndexList[idx][1]][0]
                subject += refSeqList[alignIndexList[idx][0]][0]
            elif alignIndexList[idx][0] >= 0:
                gaps += 1
                midline += " "
                query += "-"
                subject += refSeqList[alignIndexList[idx][0]][0]
            elif alignIndexList[idx][1] >= 0:
                gaps += 1
                midline += " "
                query += authSeqList[alignIndexList[idx][1]][0]
                subject += "-"
            else:
                gaps += 1
                midline += " "
                query += "-"
                subject += "-"
            #
        #
        seq_sim = float(identity) / float(end - start + 1)
# Remove the similarity cutoff per ticket DAOTHER-7903
#       if seq_sim < 0.7:
#           return False, {}
#       #
        retD = {}
        retD["db_length"] = str(len(refSeqList))
        retD["query_length"] = str(len(authSeqList))
        retD["identity"] = str(identity)
        retD["positive"] = str(identity)
        retD["gaps"] = str(gaps)
        retD["midline"] = midline
        retD["query"] = query
        retD["queryFrom"] = str(alignIndexList[start][1] + int(seqNumBeg))
        retD["queryTo"] = str(alignIndexList[end][1] + int(seqNumBeg))
        retD["subject"] = subject
        retD["hitFrom"] = str(alignIndexList[start][0] + 1)
        retD["hitTo"] = str(alignIndexList[end][0] + 1)
        retD["alignLen"] = str(end - start + 1)
        retD["match_length"] = str(end - start + 1)
        retD["seq_sim"] = seq_sim
        retD["sort_metric"] = "404"
        retD["sort_order"] = "1"
        #
        if (identity + mutation) == len(authSeqList):
            # Add for DAOTHER-9536
            if gaps == 0:
                return True, True, retD
            else:
                return False, True, retD
            #
        else:
            return False, False, retD
        #

    def __getReferenceList(self, sequence, polyTypeCode, refSeqBeg, refSeqEnd, reverseOrder):
        """Convert the one-letter code sequence from the reference resource to internal indexed list
        format seqIdx=[(3-letter-code, ref-db-index, comment, position in sequence (1-length), 1-letter code]
        """
        if reverseOrder:
            if (polyTypeCode == "RNA") or (polyTypeCode == "XNA") or (polyTypeCode == "DNA"):
                complimentSeq = self.__srd.compliment1NA(sequence[refSeqBeg - 1 : refSeqEnd][::-1], polyTypeCode)
                return self.__srd.cnv1To3ListIdx(complimentSeq, refSeqEnd, polyTypeCode, indexStep=-1)
            else:
                return self.__srd.cnv1To3ListIdx(sequence[refSeqBeg - 1 : refSeqEnd][::-1], refSeqEnd, polyTypeCode, indexStep=-1)
            #
        else:
            return self.__srd.cnv1To3ListIdx(sequence[refSeqBeg - 1 : refSeqEnd], refSeqBeg, polyTypeCode)
        #

    def __toList(self, strIn):
        """ Convert one letter sequence string into list
        """
        sL = []
        count = 0
        for ss in strIn:
            if ss in string.whitespace:
                continue
            #
            sL.append((ss, str(count)))
            count += 1
        #
        return sL


def testmain():
    fetchUtil = FetchReferenceSequenceUtils(siteId=os.getenv("WWPDB_SITE_ID"), verbose=True)
    # myD=fetchUtil.fetchReferenceSequenceWithSeqMatch("UNP", "SYUA_HUMAN", "VVHGVATVAEKTK")
    err, myD, myList = fetchUtil.fetchReferenceSequence("UNP", "A0A2X2RSX5", None)
    if err:
        print(err)
    #
    for k, v in myD.items():
        print(k + " = " + str(v))
    #
    for myTuple in myList:
        print(myTuple)
    #


if __name__ == "__main__":
    testmain()
