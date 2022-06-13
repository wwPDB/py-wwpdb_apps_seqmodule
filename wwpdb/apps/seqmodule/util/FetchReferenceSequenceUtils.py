##
# File:  FetchReferenceSequenceUtils.py
# Date:  21-Nov-2018
#
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
import sys

from wwpdb.apps.seqmodule.io.FetchSeqInfoUtils import fetchUniProt, fetchNcbiGi
from wwpdb.apps.seqmodule.util.SequenceReferenceData import SequenceReferenceData


class FetchReferenceSequenceUtils(object):
    """Fetch reference sequence data."""

    def __init__(self, siteId="WWPDB_DEPLOY_TEST", seqReferenceData=None, verbose=False, log=sys.stderr):
        self.__verbose = verbose
        self.__lfh = log
        self.__siteId = siteId
        self.__srd = seqReferenceData
        if not self.__srd:
            self.__srd = SequenceReferenceData(self.__verbose, self.__lfh)
        #
        self.__accCode = ""
        self.__refInfoD = {}

    def fetchReferenceSequence(self, dbName, dbAccession, dbIsoform, polyTypeCode="AA", refSeqBeg=None, refSeqEnd=None):
        """return error message, sequence feature dictionary, sequence list"""
        if not self.__refInfoD:
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
            self.__accCode, self.__refInfoD = self.__getRefInfo(dbName, dbAccession, dbIsoform, start, end)
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

    def fetchReferenceSequenceWithSeqMatch(self, dbName, dbAccession, authSeqs):
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
        self.__accCode, self.__refInfoD = self.__getRefInfo(dbName, dbAccession, dbIsoform, 0, 0)
        #
        if (not self.__refInfoD) or ("sequence" not in self.__refInfoD):
            return {}
        #
        idx = self.__refInfoD["sequence"].find(authSeqs)
        if idx == -1:
            return {}
        #
        refSeqBeg = idx + 1
        refSeqEnd = refSeqBeg + len(authSeqs) - 1
        if dbName in ["UNP", "SP", "TR"]:
            self.__accCode, self.__refInfoD = self.__getRefInfo(dbName, dbAccession, dbIsoform, refSeqBeg, refSeqEnd)
        #
        self.__refInfoD["db_length"] = len(self.__refInfoD["sequence"])
        self.__refInfoD["identity"] = len(authSeqs)
        self.__refInfoD["positive"] = len(authSeqs)
        self.__refInfoD["gaps"] = 0
        self.__refInfoD["midline"] = authSeqs
        self.__refInfoD["query"] = authSeqs
        self.__refInfoD["queryFrom"] = 1
        self.__refInfoD["queryTo"] = len(authSeqs)
        self.__refInfoD["subject"] = authSeqs
        self.__refInfoD["hitFrom"] = refSeqBeg
        self.__refInfoD["hitTo"] = refSeqEnd
        self.__refInfoD["alignLen"] = len(authSeqs)
        self.__refInfoD["match_length"] = len(authSeqs)
        self.__refInfoD["sort_metric"] = 404
        self.__refInfoD["sort_order"] = 1
        #
        return self.__refInfoD

    def __getRefInfo(self, dbName, dbAccession, dbIsoform, start, end):
        """Fetch sequence data from Uniprot or GeneBank database"""
        dbResource = self.__srd.convertDbNameToResource(dbName)
        #
        if dbResource in ["UNP"]:
            idCode = str(dbAccession)
            if dbIsoform is not None and len(dbIsoform) > 0:
                idCode = str(dbIsoform)
            #
            dt = fetchUniProt(idTupleList=[(idCode, start, end)], verbose=self.__verbose, log=self.__lfh)
            #
            infoD = {}
            if (idCode, start, end) in dt:
                infoD = dt[(idCode, start, end)]
            elif (dbAccession, start, end) in dt:
                infoD = dt[(dbAccession, start, end)]
            elif len(dt.values()) == 1:
                if ("db_code" in dt.values()[0]) and (dt.values()[0]["db_code"] == dbAccession):
                    infoD = dt.values()[0]
                #
            #
            if infoD and ("db_name" not in infoD) and ("db_accession" in infoD):
                # guess --
                if infoD["db_accession"][0] in ["P", "Q", "O"]:
                    infoD["db_name"] = "SP"
                else:
                    infoD["db_name"] = "TR"
                #
            #
            return idCode, self.__addingMissingKey(infoD)
        elif dbResource in ["GB", "DBJ", "EMB", "EMBL", "REF"]:
            infoD = fetchNcbiGi(dbAccession, xmlPath=None, siteId=self.__siteId)
            if infoD:
                infoD["db_accession"] = dbAccession
                infoD["db_name"] = dbName
            #
            return dbAccession, self.__addingMissingKey(infoD)
        #
        return dbAccession, {}

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

    def __addingMissingKey(self, myD):
        """ """
        defaultKeys = (
            "db_name",
            "db_accession",
            "db_code",
            "db_isoform",
            "db_description",
            "db_isoform_description",
            "name",
            "keyword",
            "sequence",
            "comments",
            "synonyms",
            "source_scientific",
            "source_strain",
            "taxonomy_id",
            "gene",
            "source_common",
            "ec",
        )
        for key in defaultKeys:
            if key not in myD:
                myD[key] = ""
            #
        #
        return myD


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
