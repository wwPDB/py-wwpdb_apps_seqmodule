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
from wwpdb.utils.seqdb_v2.FetchNcbiXml import FetchFullNcbiXml, FetchNcbiXml
from wwpdb.utils.seqdb_v2.FetchUnpXml import FetchUnpXml
from wwpdb.utils.config.ConfigInfo import ConfigInfo


def fetchUniProt(idTupleList=None, filePath=None, verbose=False, log=sys.stderr):
    """ """
    if idTupleList is None:
        idTupleList = []
    d = {}
    if not idTupleList:
        return d
    #
    idCodeList = []
    for idTuple in idTupleList:
        if idTuple[0] not in idCodeList:
            idCodeList.append(idTuple[0])
        #
    #
    if not idCodeList:
        return d
    #
    fobj = FetchUnpXml(verbose=verbose, log=log)
    ok = fobj.fetchList(idCodeList)
    if filePath is not None:
        fobj.writeUnpXml(filePath)
    #
    if ok:
        resultDicts = fobj.getResult()
        multiResultDicts = fobj.getMultipleResultDict()
        # filter any redundant annotations --
        #
        for idTuple in idTupleList:
            if idTuple[0] not in resultDicts:
                continue
            #
            vd = resultDicts[idTuple[0]]
            #
            if (idTuple[0] in multiResultDicts) and multiResultDicts[idTuple[0]]:
                found = False
                diff = -1
                for retD in multiResultDicts[idTuple[0]]:
                    if ("begin" not in retD) or ("end" not in retD) or (idTuple[1] and (idTuple[1] < retD["begin"])) or \
                       (idTuple[2] and (idTuple[2] > retD["end"])):
                        continue
                    #
                    if diff == -1:
                        found = True
                        vd = retD
                        try:
                            diff = (int(idTuple[1]) - int(retD["begin"])) + (int(retD["end"]) - int(idTuple[2]))
                        except:  # noqa: E722 pylint: disable=bare-except
                            pass
                        #
                    else:
                        try:
                            diff1 = (int(idTuple[1]) - int(retD["begin"])) + (int(retD["end"]) - int(idTuple[2]))
                            if diff1 < diff:
                                diff = diff1
                                vd = retD
                            #
                        except:  # noqa: E722 pylint: disable=bare-except
                            pass
                        #
                    #
                #
                if (not found) and ("all" in multiResultDicts[idTuple[0]][-1]) and (multiResultDicts[idTuple[0]][-1]["all"] == "yes"):
                    vd = multiResultDicts[idTuple[0]][-1]
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
    #
    return d


def fetchNcbiGi(giIdCode, xmlPath=None, siteId=None):
    """ """
    cI = ConfigInfo(siteId)
    apikey = cI.get("NCBI_API_KEY", None)
    fetchobj = FetchFullNcbiXml(giIdCode, "Nucleotide", apikey=apikey)
    if xmlPath is not None:
        fetchobj.WriteNcbiXml(filename=xmlPath)
    #
    return fetchobj.ParseNcbiXmlData()


def fetchNcbiSummary(giIdCode, xmlPath=None, siteId=None):
    """ """
    cI = ConfigInfo(siteId)
    apikey = cI.get("NCBI_API_KEY", None)
    fetchobj = FetchNcbiXml(giIdCode, "Nucleotide", apikey=apikey)
    if xmlPath is not None:
        fetchobj.WriteNcbiXml(filename=xmlPath)
    #
    return fetchobj.ParseNcbiXmlData()
