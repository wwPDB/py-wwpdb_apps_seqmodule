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
##
"""
Utility methods to retrieve information from NCBI & Uniprot databases
"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.09"

import os, sys
from wwpdb.utils.seqdb_v2.FetchNcbiXml import FetchFullNcbiXml, FetchNcbiXml
from wwpdb.utils.seqdb_v2.FetchUniProtEntry import FetchUniProtEntry
from wwpdb.utils.config.ConfigInfo import ConfigInfo

def fetchUniProt(siteId=None, verbose=False, log=sys.stderr, idCodeList=None, filePath=None):
    """
    """
    d = {}
    fobj = FetchUniProtEntry(siteId=siteId, verbose=verbose, log=log)
    ok = fobj.fetchList(idCodeList)
    if filePath is not None:
        fobj.writeUnpXml(filePath)
    #
    if ok:
        d = fobj.getResult()
        # filter any redundant annotations --
        #
        for (acId, vd) in d.items():
            for k in vd.keys():
                if k in ['ec', 'comments', 'synonyms']:
                    v = vd[k]
                    oL = []
                    tL = v.split(',')
                    for it in tL:
                        sV = str(it).strip()
                        if k == 'ec':
                            sVL = sV.split('.')
                            oL1 = []
                            for s1 in sVL:
                                if str(s1).startswith('n'):
                                    oL1.append('-')
                                else:
                                    oL1.append(s1)
                                #
                            #
                            sV = '.'.join(oL1)
                        #
                        if sV in oL:
                            continue
                        else:
                            oL.append(sV)
                        #
                    #
                    vd[k] = ','.join(oL)
                #
            #
        #
    #
    return d

def fetchNcbiGi(giIdCode, xmlPath=None, siteId=None):
    """
    """
    cI = ConfigInfo(siteId)
    apikey = cI.get('NCBI_API_KEY', None)  
    fetchobj = FetchFullNcbiXml(giIdCode, 'Nucleotide', apikey=apikey)
    if xmlPath is not None:
        fetchobj.WriteNcbiXml(filename=xmlPath)
    #
    return fetchobj.ParseNcbiXmlData()

def fetchNcbiSummary(giIdCode, xmlPath=None, siteId=None):
    """
    """
    cI = ConfigInfo(siteId)
    apikey = cI.get('NCBI_API_KEY', None)  
    fetchobj = FetchNcbiXml(giIdCode, 'Nucleotide', apikey=apikey)
    if (xmlPath is not None):
        fetchobj.WriteNcbiXml(filename=xmlPath)
    #
    return fetchobj.ParseNcbiXmlData()
