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

import os, sys, traceback

from wwpdb.apps.seqmodule.io.FetchSeqInfoUtils import fetchUniProt,fetchNcbiGi
from wwpdb.apps.seqmodule.util.SequenceReferenceData import SequenceReferenceData

class FetchReferenceSequenceUtils(object):
    """ Fetch reference sequence data.
    """
    def __init__(self, siteId="WWPDB_DEPLOY_TEST", seqReferenceData=None, verbose=False, log=sys.stderr):
        self.__verbose = verbose
        self.__lfh = log
        self.__siteId = siteId
        self.__srd = seqReferenceData
        if not self.__srd:
            self.__srd = SequenceReferenceData(self.__verbose, self.__lfh)
        #

    def fetchReferenceSequence(self, dbName, dbAccession, dbIsoform, polyTypeCode="AA", refSeqBeg=None, refSeqEnd=None):
        """ return error message, sequence feature dictionary, sequence list
        """
        accCode,refInfoD = self.__getRefInfo(dbName, dbAccession, dbIsoform)
        if (not refInfoD) or ("sequence" not in refInfoD):
            return "Fetch reference sequence [ dbName=" + dbName + ", Accession=" + accCode + "] failed.",{},[]
        #
        seqLength = len(refInfoD["sequence"])
        #
        try:
            refSeqBeg = int(str(refSeqBeg))
        except:
            refSeqBeg = 1
        #
        try:
            refSeqEnd = int(str(refSeqEnd))
        except:
            refSeqEnd = seqLength
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
            errorMessage += "Invalid SEQ END number: " + str(refSeqEnd) + " > " + str(seqLength) + " ( sequence length of " + accCode + " )."
        #
        if errorMessage:
            return errorMessage,{},[]
        #
        if refSeqBeg is None:
            refSeqBeg = 1
        #
        if refSeqEnd is None:
            refSeqEnd = len(sequence)
        #
        refInfoD["db_length"] = seqLength
        refInfoD["hitFrom"] = refSeqBeg
        refInfoD["hitTo"] = refSeqEnd
        return "",refInfoD,self.__getReferenceList(refInfoD["sequence"], polyTypeCode, refSeqBeg, refSeqEnd)

    def __getRefInfo(self, dbName, dbAccession, dbIsoform):
        """ Fetch sequence data from Uniprot or GeneBank database
        """
        dbResource = self.__srd.convertDbNameToResource(dbName)
        #
        if dbResource in ["UNP"]:
            idCode = str(dbAccession)
            if dbIsoform is not None and len(dbIsoform) > 0:
                idCode = str(dbIsoform)
            #
            dt = fetchUniProt(siteId=self.__siteId, verbose=self.__verbose, log=self.__lfh, idCodeList=[idCode])
            #
            infoD = {}
            if idCode in dt:
                infoD = dt[idCode]
            elif dbAccession in dt:
                infoD = dt[dbAccession]
            #
            if infoD and ("db_name" not in infoD):
                # guess --
                if dbAccession[0] in ["P", "Q", "O"]:
                    infoD["db_name"] = "SP"
                else:
                    infoD["db_name"] = "TR"
                #
            #
            return idCode,self.__addingMissingKey(infoD)
        elif dbResource in ["GB", "DBJ", "EMB"]:
            infoD = fetchNcbiGi(dbAccession, xmlPath=None)
            if infoD:
                infoD["db_accession"] = dbAccession
                infoD["db_name"] = dbName
            #
            return dbAccession,self.__addingMissingKey(infoD)
        #
        return dbAccession,{}

    def __getReferenceList(self, sequence, polyTypeCode, refSeqBeg, refSeqEnd):
        """  Convert the one-letter code sequence from the reference resource to internal indexed list 
             format seqIdx=[(3-letter-code, ref-db-index, comment, position in sequence (1-length), 1-letter code]
        """
        return self.__srd.cnv1To3ListIdx(sequence[refSeqBeg-1:refSeqEnd], refSeqBeg, polyTypeCode)

    def __addingMissingKey(self, myD):
        """
        """
        defaultKeys = ( "db_name", "db_accession", "db_code", "db_isoform", "db_description", "db_isoform_description", "name", "keyword", "sequence", \
                        "comments", "synonyms", "source_scientific", "source_strain", "taxonomy_id", "gene", "source_common", "ec" )
        for key in defaultKeys:
            if key not in myD:
                myD[key] = ""
            #
        #
        return myD

if __name__ == "__main__":
    siteId = os.getenv("WWPDB_SITE_ID")
    fetchUtil = FetchReferenceSequenceUtils(siteId=os.getenv("WWPDB_SITE_ID"), verbose=True)
    err,myD,myList=fetchUtil.fetchReferenceSequence("UNP", "A0A2X2RSX5", None)
    if err:
        print err
    #
    for k,v in myD.items():
        print k + ' = ' + str(v)
    #
    for myTuple in myList:
        print myTuple
    #
