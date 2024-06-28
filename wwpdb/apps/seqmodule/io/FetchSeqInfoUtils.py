##
# File:  FetchSeqInfoUtils.py
# Date:  21-Feb-2013
#
# Updates:
# 21-Feb-2013  jdw include methods to execute Blast search services.
# 24-Feb-2013  jdw added tests for fetch methods for UniProt and  NCBI sequence/taxonomy entries.
# 25-Feb-2013  jdw add siteId to constructor to support TaxonmyUtils class
# 17-Apr-2013  jdw Add methods for local sequence search
# 19-Apr-2013  jdw Local protein sequence search now turned on
# 22-Apr-2013  jdw Revert to entry summary for RNA sequences.
# 22-Apr-2013  jdw Use ncbi taxonomy database to lookup missing source organism names after
#                  nucleotide search
# 04-Nov-2013  jdw update sort score
# 15-Dec-2013  jdw trap cases missing taxonomy in either input or in returned reference list
# 27-Jun-2014  jdw overhaul the runBlastLocal() and optimize for children of taxId 562.
#                  reduce the number of id lookups and try to find distant matching SP entries.
#  1-Aug-2014  jdw update handling of annotation fetch failures
#  8-Dec-2015  jdw change isoform annoation processing -
# 30-Dec-2020   zf Use FetchUnpXml instead of FetchUniProtEntry class. Allow the location specific feature names if exist.
# 29-Sep-2022   zf Added input 'begin' & 'end' values checking from idTupleList in fetchUniProt() method.
# 11-Nov-2022   zf Create FetchSeqInfoUtils class with fetchUniProt(), fetchNcbiGi(), fetchNcbiSummary(), and getRefInfo() methods.
#                  getRefInfo() method is moved from FetchReferenceSequenceUtils.__getRefInfo() method.
#  2-Jun-2024   zf Added persistent objects self.__resultDicts & self.__multiResultDicts. Fixed the bug related to missing start & end sequence numbers.
##
"""
Utility methods to retrieve information from NCBI & Uniprot databases
"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.09"

import sys

from wwpdb.apps.seqmodule.util.SequenceReferenceData import SequenceReferenceData
from wwpdb.utils.config.ConfigInfo import ConfigInfo
from wwpdb.utils.seqdb_v2.FetchNcbiXml import FetchFullNcbiXml, FetchNcbiXml
from wwpdb.utils.seqdb_v2.FetchUnpXml import FetchUnpXml


class FetchSeqInfoUtils(object):
    """Fetch reference sequence data.
    """
    def __init__(self, siteId="WWPDB_DEPLOY_TEST", seqReferenceData=None, verbose=False, log=sys.stderr):
        """
        """
        self.__verbose = verbose
        self.__lfh = log
        self.__siteId = siteId
        self.__cI = ConfigInfo(self.__siteId)
        self.__apikey = self.__cI.get("NCBI_API_KEY", None)
        self.__srd = seqReferenceData
        if not self.__srd:
            self.__srd = SequenceReferenceData(self.__verbose, self.__lfh)
        #
        self.__resultDicts = {}
        self.__multiResultDicts = {}

    def getRefInfo(self, dbName, dbAccession, dbIsoform, start, end, addMissingKeyFlag=True):
        """Fetch sequence data from Uniprot or GeneBank database.
        """
        dbResource = self.__srd.convertDbNameToResource(dbName)
        #
        if dbResource in ["UNP"]:
            idCode = str(dbAccession)
            if dbIsoform is not None and len(dbIsoform) > 0:
                idCode = str(dbIsoform)
            #
            dt = self.fetchUniProt(idTupleList=[(idCode, start, end)])
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
            if not infoD:
                return idCode, {}
            #
            if infoD and ("db_name" not in infoD) and ("db_accession" in infoD):
                # guess --
                if infoD["db_accession"][0] in ["P", "Q", "O"]:
                    infoD["db_name"] = "SP"
                else:
                    infoD["db_name"] = "TR"
                #
            #
            if addMissingKeyFlag:
                return idCode, self.__addingMissingKey(infoD)
            else:
                return idCode, infoD
            #
        elif dbResource in ["GB", "DBJ", "EMB", "EMBL", "REF"]:
            infoD = self.fetchNcbiGi(dbAccession)
            if not infoD:
                return dbAccession, {}
            #
            if infoD:
                infoD["db_accession"] = dbAccession
                infoD["db_name"] = dbName
            #
            if addMissingKeyFlag:
                return dbAccession, self.__addingMissingKey(infoD)
            else:
                return dbAccession, infoD
            #
        #
        return dbAccession, {}

    def fetchUniProt(self, idTupleList=None, filePath=None):
        """ Fetch Uniprot reference sequence data.
        """
        if idTupleList is None:
            idTupleList = []
        #
        d = {}
        if not idTupleList:
            return d
        #
        idCodeList = []
        foundCache = False
        for idTuple in idTupleList:
            if idTuple[0] in self.__resultDicts:
                foundCache = True
                continue
            #
            if idTuple[0] not in idCodeList:
                idCodeList.append(idTuple[0])
            #
        #
        if (not foundCache) and (not idCodeList):
            return d
        #
        self.__fetchUniProtXmlFile(idCodeList, filePath)
        #
        # filter any redundant annotations --
        #
        for idTuple in idTupleList:
            if idTuple[0] not in self.__resultDicts:
                continue
            #
            vd = self.__resultDicts[idTuple[0]]
            seqLength = 0
            if ("sequence" in vd) and vd["sequence"] and (len(vd["sequence"]) > 1):
                seqLength = len(vd["sequence"])
            #
            if (idTuple[0] in self.__multiResultDicts) and self.__multiResultDicts[idTuple[0]]:
                try:
                    idTuple1 = int(idTuple[1])
                except:  # noqa: E722 pylint: disable=bare-except
                    idTuple1 = 0
                #
                try:
                    idTuple2 = int(idTuple[2])
                except:  # noqa: E722 pylint: disable=bare-except
                    idTuple2 = 0
                #
                if ((idTuple1 == 0) or (idTuple2 == 0)) and (seqLength > 1):
                    idTuple1 = 1
                    idTuple2 = seqLength
                #
                found = False
                if (idTuple1 > 0) and (idTuple2 > 0) and (idTuple1 <= idTuple2):
                    diff = -1
                    for retD in self.__multiResultDicts[idTuple[0]]:
                        if ("begin" not in retD) or ("end" not in retD) or (idTuple1 and (idTuple1 < retD["begin"])) or \
                           (idTuple2 and (idTuple2 > retD["end"])):
                            continue
                        #
                        if diff == -1:
                            found = True
                            vd = retD
                            try:
                                diff = (int(idTuple1) - int(retD["begin"])) + (int(retD["end"]) - int(idTuple2))
                            except:  # noqa: E722 pylint: disable=bare-except
                                pass
                            #
                        else:
                            try:
                                diff1 = (int(idTuple1) - int(retD["begin"])) + (int(retD["end"]) - int(idTuple2))
                                if diff1 < diff:
                                    diff = diff1
                                    vd = retD
                                #
                            except:  # noqa: E722 pylint: disable=bare-except
                                pass
                            #
                        #
                    #
                #
                if (not found) and ("all" in self.__multiResultDicts[idTuple[0]][-1]) and (self.__multiResultDicts[idTuple[0]][-1]["all"] == "yes"):
                    vd = self.__multiResultDicts[idTuple[0]][-1]
                #
            #
            for k in vd.keys():
                if k in ["ec", "comments", "synonyms"]:
                    v = vd[k]
                    oL = []
                    tL = v.split(",")
                    for it in tL:
                        sV = str(it).strip()
                        if k == "ec":
                            sVL = sV.split(".")
                            oL1 = []
                            for s1 in sVL:
                                if str(s1).startswith("n"):
                                    oL1.append("-")
                                else:
                                    oL1.append(s1)
                                #
                            #
                            sV = ".".join(oL1)
                        #
                        if sV in oL:
                            continue
                        else:
                            oL.append(sV)
                        #
                    #
                    vd[k] = ",".join(oL)
                #
            #
            d[idTuple] = vd
        #
        return d

    def fetchNcbiGi(self, giIdCode, xmlPath=None):
        """ Fetch GeneBank reference data.
        """
        fetchobj = FetchFullNcbiXml(giIdCode, "Nucleotide", apikey=self.__apikey)
        if xmlPath is not None:
            fetchobj.WriteNcbiXml(filename=xmlPath)
        #
        return fetchobj.ParseNcbiXmlData()

    def fetchNcbiSummary(self, giIdCode, xmlPath=None):
        """ Fetch GeneBank reference summary.
        """
        fetchobj = FetchNcbiXml(giIdCode, "Nucleotide", apikey=self.__apikey)
        if xmlPath is not None:
            fetchobj.WriteNcbiXml(filename=xmlPath)
        #
        return fetchobj.ParseNcbiXmlData()

    def __addingMissingKey(self, myD):
        """ Add missing key items
        """
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

    def __fetchUniProtXmlFile(self, idCodeList=None, filePath=None):
        """ Fetch Uniprot Xml File
        """
        # Handles empty list as well
        if not idCodeList:
            return
        #
        fobj = FetchUnpXml(verbose=self.__verbose, log=self.__lfh)
        ok = fobj.fetchList(idCodeList)
        if filePath is not None:
            fobj.writeUnpXml(filePath)
        #
        if ok:
            self.__resultDicts.update(fobj.getResult())
            self.__multiResultDicts.update(fobj.getMultipleResultDict())
        #
