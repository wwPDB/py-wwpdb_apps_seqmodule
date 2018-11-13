##
# File:  ReferenceSequenceUtils.py
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
Methods to manage local and external sequence search services and package matching reference sequence data.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.09"

import sys
import os
import time
import traceback
from operator import itemgetter
from wwpdb.apps.seqmodule.io.PdbxIoUtils import ReferenceSequenceIo, PdbxFileIo
from wwpdb.apps.seqmodule.io.TaxonomyUtils import TaxonomyUtils
from wwpdb.apps.seqmodule.util.SequenceReferenceData import SequenceReferenceData

from wwpdb.utils.seqdb_v2.UnpBlastService import UnpBlastService
from wwpdb.utils.seqdb_v2.ReadUnpBlastXml import ReadUnpBlastXmlString
from wwpdb.utils.seqdb_v2.NcbiBlastService import NcbiBlastService
from wwpdb.utils.seqdb_v2.ReadNcbiBlastXml import ReadNcbiBlastXmlString

from wwpdb.utils.seqdb_v2.FetchUnpXml import FetchUnpXml
from wwpdb.utils.seqdb_v2.FetchNcbiXml import FetchFullNcbiXml, FetchNcbiXml
from wwpdb.utils.seqdb_v2.FetchUniProtEntry import FetchUniProtEntry
from wwpdb.utils.dp.RcsbDpUtility import RcsbDpUtility
from wwpdb.apps.seqmodule.io.BlastPlusReader import BlastPlusReader


def compareElement0(a, b):
    if a[0] > b[0]:
        return True
    return False


class ReferenceSequenceUtils(object):
    """
    Execute search service and assemble reference sequence data.

    """

    def __init__(self, siteId="WWPDB_DEPLOY_TEST", verbose=False, log=sys.stderr):
        self.__verbose = verbose
        self.__lfh = log
        self.__siteId = siteId
        self.__taxUtils = TaxonomyUtils(siteId=self.__siteId, verbose=self.__verbose, log=self.__lfh)
        #
        self.__minRnaSequenceSearchLength = 50
        # JDW
        self.__cleanUp = True
        self.__debug = False
        #

    def getMatchingSequenceFilePath(self, entryId, entityId):
        try:
            return str(entryId) + "." + str(entityId) + '-info.cif'
        except:
            return ''

    def searchEntities(self, entityD=None, saveBlast=True, filePath='.', filePrefix="refseq", addSortMetric=True, localSearch=True):
        """  Process the a search for the entity described in the input dictionary.   Sequence searches
             are performed for each entity part in the input feature list.

             LocalSearch = boolean toggles between local and remove sequence search services -

            Return a list of hitLists for matching reference sequences extracted from the BLAST search results
        """
        rL = []
        seqType = entityD['POLYMER_LINKING_TYPE']
        oneLetterCodeSeq = entityD['SEQ_ENTITY_1_CAN']
        entityId = entityD['ENTITY_ID']
        #
        maxHits = 50
        #
        if self.__verbose:
            self.__lfh.write("\n+ReferenceSequenceUtils.searchEntities() STARTING with polymer type  %s sequence = %s\n" % (seqType, oneLetterCodeSeq))

        for partNum, fD in enumerate(entityD['PART_LIST'], start=1):
            maxHitsSearch = 100
            seqNumBeg = fD['SEQ_NUM_BEG']
            seqNumEnd = fD['SEQ_NUM_END']
            taxId = fD['SOURCE_TAXID']
            #
            #  Use a deeper search if we are
            extendSearch = False
            anD = self.__taxUtils.getAncestors(taxId)
            if (('p_id' in anD and anD['p_id'] in ['562']) or ('gp_id' in anD and anD['gp_id'] in ['562'])):
                extendSearch = True

            if (seqType in ['polypeptide(L)', 'polypeptide(D)', 'polypeptide']):
                maxHitsSearch = 5000
            else:
                maxHitsSearch = 5000

            seqPartType = fD['SEQ_PART_TYPE']
            seqPartId = fD['SEQ_PART_ID']
            if not str(seqPartType).lower().startswith('biol'):
                if self.__verbose:
                    self.__lfh.write("+ReferenceSequenceUtils.searchEntities() skipping sequence part %r type %r\n" % (seqPartId, seqPartType))
                continue

            iBeg = int(seqNumBeg) - 1
            iEnd = int(seqNumEnd)
            sequence = oneLetterCodeSeq[iBeg:iEnd]
            if self.__verbose:
                self.__lfh.write("\n+ReferenceSequenceUtils.searchEntities() starting database search for entityId %r partId %r range begin %r end %r taxId %r sequence length = %d maxHits %d maxHitsSearch %d\n" %
                                 (entityId, partNum, iBeg, iEnd, taxId, len(sequence), maxHits, maxHitsSearch))

            if saveBlast:
                xmlFilePath = os.path.join(filePath, filePrefix + "_blast-match_E" + entityId + "_P" + str(partNum) + ".xml")
            else:
                xmlFilePath = None

            if (localSearch):
                hitList = self.__runBlastLocal(sequence, polyType=seqType, blastResultPath=xmlFilePath, tmpPath=filePath, maxHitsSearch=maxHitsSearch, maxHitsSave=maxHits)
            else:
                hitList = self.__runBlastService(sequence, polyType=seqType, xmlFilePath=xmlFilePath)

            for hit in hitList:
                hit['beg_seq_num'] = seqNumBeg
                hit['end_seq_num'] = seqNumEnd
                hit['fragment_id'] = str(seqPartId)

            #
            if (addSortMetric):
                self.__sortHitList(hitList, authTaxId=taxId, authAccessionCode=None)
            #
            if self.__verbose:
                self.__lfh.write(
                    "+ReferenceSequenceUtils.searchEntities() completed for entityId %r partId %r taxId %r (hitList length = %d save limit is %d)  \n" %
                    (entityId, partNum, taxId, len(hitList), maxHits))
            if len(hitList) > maxHits:
                rL.append(hitList[len(hitList) - maxHits:])
            else:
                rL.append(hitList)

        if self.__verbose:
            self.__lfh.write("+ReferenceSequenceUtils.searchEntities() COMPLETED search for entity %r returning results for %d part(s)\n" % (entityId, len(rL)))

        return rL

    def __runBlastLocal(self, oneLetterCodeSeq, polyType, blastResultPath=None, tmpPath='.', maxHitsSearch=200, maxHitsSave=50, maxIdLookup=100):
        """ Internal method to execute the sequence search according to input polymer type.

            List of match results --

        """
        if self.__verbose:
            self.__lfh.write("+ReferenceSequenceUtils.__runBlastLocal() STARTING %s Launch search for %s sequence = %s\n" % (self.__siteId, polyType, oneLetterCodeSeq))
            self.__lfh.flush()
        #
        timeBegin = time.time()
        hitList = []

        resultPath = blastResultPath
        if resultPath is None:
            resultPath = os.path.join(tmpPath, "blast-result.xml")

        resultPath = os.path.abspath(resultPath)
        tmpPathAbs = os.path.abspath(tmpPath)

        if polyType in ['polypeptide(L)', 'polypeptide(D)', 'polypeptide']:
            dp = RcsbDpUtility(tmpPath=tmpPathAbs, siteId=self.__siteId, verbose=True)
            dp.addInput(name="one_letter_code_sequence", value=oneLetterCodeSeq)
            dp.addInput(name="db_name", value="my_uniprot_all")
            dp.addInput(name="num_threads", value="4")
            dp.addInput(name="max_hits", value=maxHitsSearch)
            dp.op("seq-blastp")
            dp.exp(resultPath)
            #
            if (self.__verbose):
                timeEnd = time.time()
                self.__lfh.write("+ReferenceSequenceUtils.__runBlastLocal() Blast search completed after %.2f seconds\n" % (timeEnd - timeBegin))

            bpr = BlastPlusReader(verbose=self.__verbose, log=self.__lfh)
            searchHitList = bpr.readFile(filePath=resultPath)
            if (self.__verbose):
                timeEnd = time.time()
                self.__lfh.write("+ReferenceSequenceUtils.__runBlastLocal() Read %d hits from blast output file after %.2f seconds\n" % (len(searchHitList), timeEnd - timeBegin))
            #
            if (self.__cleanUp):
                dp.cleanup()

            hitList = []
            idCodeList = []
            iCount = 0
            iTotal = 0
            for hit in searchHitList:
                if 'db_accession' in hit and 'db_name' in hit:
                    iCount += 1
                    if ((iCount > maxHitsSave) and (hit['db_name'] not in ['SP', 'sp'])):
                        continue
                    if iTotal > maxIdLookup:
                        break
                    if 'db_isoform' in hit and (len(hit['db_isoform']) > 0) and (hit['db_isoform'] not in ['.', '?']):
                        idCodeList.append(hit['db_isoform'])
                    else:
                        idCodeList.append(hit['db_accession'])
                    iTotal += 1
            idCodeList = list(set(idCodeList))
            if (self.__verbose):
                timeEnd = time.time()
                self.__lfh.write("+ReferenceSequenceUtils.__runBlastLocal() Blast search completed after %.2f seconds\n" % (timeEnd - timeBegin))
                self.__lfh.write("+ReferenceSequenceUtils.__runBlastLocal() Starting fetch of sequence database hits for %d reference idcodes\n" % len(idCodeList))

            #
            # Incorporate non-overlapping additional details about each matching sequence entry.
            #   ( Still using the REST API to get the entry information )
            #
            if (self.__verbose):
                self.__lfh.write("+ReferenceSequenceUtils.__runBlastLocal() Fetch sequence database entries for %d filtered reference idcodes\n" % len(idCodeList))
            if len(idCodeList) > 0:
                unpD = self.fetchUniProt(idCodeList)
                for hit in searchHitList:
                    acId = None
                    isIso = False
                    if 'db_isoform' in hit and (len(hit['db_isoform']) > 0) and (hit['db_isoform'] not in ['.', '?']):
                        acId = hit['db_isoform']
                        isIso = True
                    elif 'db_accession' in hit:
                        acId = hit['db_accession']
                    if self.__debug:
                        self.__lfh.write("+ReferenceSequenceUtils.__runBlastLocal() Updating blast hit[accession] %s with uniprot code %s \n" % (hit['db_accession'], acId))
                    if acId is not None:
                        dd = {}
                        if acId in unpD:
                            dd = unpD[acId]
                        elif isIso and hit['db_accession'] in unpD:
                            dd = unpD[hit['db_accession']]
                        for (k, v) in dd.items():
                            if k not in hit:
                                if (self.__debug):
                                    self.__lfh.write("+ReferenceSequenceUtils.__runBlastLocal() accCode %s adding key %r and value %r\n" % (acId, k, v))
                                hit[k] = v
                        hitList.append(hit)
            if (self.__verbose):
                self.__lfh.write("+ReferenceSequenceUtils.__runBlastLocal() Final hit list length is %d\n" % len(hitList))

        elif polyType == 'polyribonucleotide' and (len(oneLetterCodeSeq) > self.__minRnaSequenceSearchLength):
            dp = RcsbDpUtility(tmpPath=tmpPathAbs, siteId=self.__siteId, verbose=True)
            dp.addInput(name="one_letter_code_sequence", value=oneLetterCodeSeq)
            dp.addInput(name="db_name", value="my_ncbi_nt")
            dp.addInput(name="num_threads", value="4")
            dp.addInput(name="max_hits", value=maxHitsSave)
            dp.op("seq-blastn")
            dp.exp(resultPath)
            if (self.__cleanUp and not self.__debug):
                dp.cleanup()
            #
            bpr = BlastPlusReader(verbose=self.__verbose, log=self.__lfh)
            bpr.setSequenceType(polyType)
            hitList = bpr.readFile(filePath=resultPath)
            idCodeList = []
            for hit in hitList:
                if 'db_accession' in hit:
                    idCodeList.append(hit['db_accession'])
            #
            idCodeList = list(set(idCodeList))
            if (self.__verbose):
                timeEnd = time.time()
                self.__lfh.write("+ReferenceSequenceUtils.__runBlastLocal() Blast search completed in %.2f seconds\n" % (timeEnd - timeBegin))
                self.__lfh.write("+ReferenceSequenceUtils.__runBlastLocal() Starting fetch of sequence database entries for %d reference idcodes\n" % len(idCodeList))
            #
            # Incorporate non-overlapping additional details about each matching sequence entry.
            #
            if len(idCodeList) > 0:
                giD = {}
                for idCode in idCodeList:
                    summaryPath = None
                    if self.__debug:
                        summaryPath = os.path.join(tmpPathAbs, idCode + '-summary.xml')
                    giD[idCode] = self.fetchNcbiSummary(idCode, xmlPath=summaryPath)
                for hit in hitList:
                    if 'db_accession' in hit:
                        acId = hit['db_accession']
                        if acId in giD:
                            dd = giD[acId]
                            for (k, v) in dd.items():
                                if k not in hit:
                                    if (self.__debug):
                                        self.__lfh.write("+ReferenceSequenceUtils.__runBlastLocal() idCode %s adding key %r and value %r\n" % (idCode, k, v))
                                    hit[k] = v
                    if (('source_scientific' not in hit or len(hit['source_scientific']) < 2) and 'taxonomy_id' in hit):
                        # try to lookup the scientific name - Take the first hit --
                        snL = self.__taxUtils.lookUpSource(taxId=hit['taxonomy_id'])
                        hit['source_scientific'] = ''
                        if len(snL) > 0:
                            hit['source_scientific'] = snL[0]
                        if (self.__debug):
                            self.__lfh.write("+ReferenceSequenceUtils.__runBlastLocal() idCode %s assigning source from taxonomy db %r\n" % (idCode, hit['source_scientific']))

        else:
            self.__lfh.write("+ReferenceSequenceUtils.__runBlastLocal() Search failed for unknown type    = %s\n" % polyType)

        timeEnd = time.time()
        if self.__verbose:
            self.__lfh.write("+ReferenceSequenceUtils.__runBlastLocal() Search processing completed after %.2f seconds\n" % (timeEnd - timeBegin))
            self.__lfh.write("+ReferenceSequenceUtils.__runBlastLocal() COMPLETED returning hit list length       = %d\n" % len(hitList))
            self.__lfh.flush()

        return hitList

    def __runBlastService(self, oneLetterCodeSeq, polyType, xmlFilePath=None):
        """ Internal method to execute the sequence search using remote webservices according to input polymer type.

            List of match
        """
        if self.__verbose:
            self.__lfh.write("+ReferenceSequenceUtils.__runBlastService() Launch search for %s sequence = %s\n" % (polyType, oneLetterCodeSeq))
            self.__lfh.flush()
        #
        timeBegin = time.clock()
        blast_match_result = []

        # run Uniprot BLAST service for protein sequences
        if polyType in ['polypeptide(L)', 'polypeptide(D)', 'polypeptide']:
            service = UnpBlastService(oneLetterCodeSeq, verbose=self.__verbose, log=self.__lfh)
            service.RunService()
            # fetch the raw XML result from the Blast search
            xmlResult = service.GetResult()
            if xmlFilePath is not None:
                service.WriteResultFile(xmlFilePath)
            #
            if xmlResult:
                blastresult = ReadUnpBlastXmlString(xmlResult, verbose=self.__verbose, log=self.__lfh)
                blast_match_result = blastresult.GetResult()

        # Run NCBI BLAST service for RNA sequences
        elif polyType == 'polyribonucleotide' and (len(oneLetterCodeSeq) > self.__minRnaSequenceSearchLength):
            service = NcbiBlastService(oneLetterCodeSeq)
            service.RunService()
            # fetch the raw XML result from the Blast search
            xmlResult = service.GetResult()
            if xmlFilePath is not None:
                service.WriteResultFile(xmlFilePath)
            timeBlast = time.time()
            if xmlResult:
                blastresult = ReadNcbiBlastXmlString(xmlResult)
                blast_match_result = blastresult.GetResult()
        else:
            self.__lfh.write("+ReferenceSequenceUtils.__runBlastService() Search failed for unknown type    = %s\n" % polyType)

        timeEnd = time.clock()
        if self.__verbose:
            self.__lfh.write("+ReferenceSequenceUtils.__runBlastService() Search processing completed in %d seconds\n" % (timeEnd - timeBegin))
            self.__lfh.write("+ReferenceSequenceUtils.__runBlastService() Hit list length       = %d\n" % len(blast_match_result))
            self.__lfh.flush()

        return blast_match_result

    def __runBlastLocalNoNuc(self, oneLetterCodeSeq, polyType, blastResultPath=None, tmpPath='.'):
        """ Internal method to execute the sequence search according to input polymer type.

            List of match results --
        """
        if self.__verbose:
            self.__lfh.write("+ReferenceSequenceUtils.__runBlastLocal() Launch search for %s sequence = %s\n" % (polyType, oneLetterCodeSeq))
            self.__lfh.flush()
        #
        timeBegin = time.clock()
        hitList = []

        resultPath = blastResultPath
        if resultPath is None:
            resultPath = os.path.join(tmpPath, "blast-result.xml")

        resultPath = os.path.abspath(resultPath)

        if self.__verbose:
            self.__lfh.write("+ReferenceSequenceUtils.__runBlastLocal() blast result path set to = %s\n" % resultPath)

        if polyType in ['polypeptide(L)', 'polypeptide(D)', 'polypeptide']:
            dp = RcsbDpUtility(tmpPath=tmpPath, siteId=self.__siteId, verbose=True)
            dp.addInput(name="one_letter_code_sequence", value=oneLetterCodeSeq)
            dp.addInput(name="db_name", value="my_uniprot_all")
            dp.addInput(name="num_threads", value="4")
            dp.addInput(name="evalue", value="0.001")
            dp.op("seq-blastp")
            dp.exp(resultPath)
            if (self.__cleanUp):
                dp.cleanup()
            #
            bpr = BlastPlusReader(verbose=self.__verbose, log=self.__lfh)
            hitList = bpr.readFile(filePath=resultPath)
            #
            idCodeList = []
            for hit in hitList:
                if 'db_accession' in hit:
                    idCodeList.append(hit['db_accession'])
            #
            # Incorporate non-overlapping additional details about each matching sequence entry.
            #   ( Still using the REST API to get the entry information )
            #
            if len(idCodeList) > 0:
                unpD = self.fetchUniProt(idCodeList)
                for hit in hitList:
                    if 'db_accession' in hit:
                        acId = hit['db_accession']
                        if acId in unpD:
                            dd = unpD[acId]
                            for (k, v) in dd.items():
                                if k not in hit:
                                    hit[k] = v

        elif polyType == 'polyribonucleotide' and (len(oneLetterCodeSeq) > self.__minRnaSequenceSearchLength):
            #
            service = NcbiBlastService(oneLetterCodeSeq)
            service.RunService()
            # fetch the raw XML result from the Blast search
            xmlResult = service.GetResult()
            if resultPath is not None:
                service.WriteResultFile(resultPath)
            timeBlast = time.time()
            if xmlResult:
                blastresult = ReadNcbiBlastXmlString(xmlResult)
                hitList = blastresult.GetResult()
        else:
            self.__lfh.write("+ReferenceSequenceUtils.__runBlastLocal() Search failed for unknown type = %s\n" % polyType)

        timeEnd = time.clock()
        if self.__verbose:
            self.__lfh.write("+ReferenceSequenceUtils.__runBlastLocal() Search processing completed in %d seconds\n" % (timeEnd - timeBegin))
            self.__lfh.write("+ReferenceSequenceUtils.__runBlastLocal() Hit list length       = %d\n" % len(hitList))
            self.__lfh.flush()

        return hitList

    def fetchUniProt(self, idCodeList, filePath=None):
        """
        """
        self.__debug = False
        if (self.__debug):
            self.__lfh.write("+ReferenceSequenceUtils.fetchUniProt() starting with idCodeList %r\n" % idCodeList)
        d = {}
        # -- JDW updating --
        #fobj = FetchUnpXml(verbose=self.__verbose,log=self.__lfh)
        fobj = FetchUniProtEntry(siteId=self.__siteId, verbose=self.__verbose, log=self.__lfh)
        ok = fobj.fetchList(idCodeList)
        if filePath is not None:
            fobj.writeUnpXml(filePath)
        if (self.__debug):
            self.__lfh.write("+ReferenceSequenceUtils.fetchUniProt() return status for idCodeList is %r\n" % ok)
        if ok:
            d = fobj.getResult()
            if (self.__debug):
                for (acId, vd) in d.items():
                    for k, v in vd.items():
                        self.__lfh.write("+ReferenceSequenceUtils.fetchUniProt() id %s content key %r value %r\n" % (acId, k, v))

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
                                sV = '.'.join(oL1)
                            if sV in oL:
                                continue
                            else:
                                oL.append(sV)
                        vd[k] = ','.join(oL)

            if (self.__debug):
                for (acId, vd) in d.items():
                    for k, v in vd.items():
                        self.__lfh.write("+ReferenceSequenceUtils.__fetchUniProt() id %s content key %r value %r\n" % (acId, k, v))
        return d

    def fetchNcbiGi(self, giIdCode, xmlPath=None):
        """
        """
        d = {}
        fetchobj = FetchFullNcbiXml(giIdCode, 'Nucleotide')
        d = fetchobj.ParseNcbiXmlData()
        if xmlPath is not None:
            fetchobj.WriteNcbiXml(filename=xmlPath)
        if (self.__debug):
            for k, v in d.items():
                self.__lfh.write("+ReferenceSequenceUtils.__fetchGiNcbi id content key %r value %r\n" % (k, v))
        return d

    def fetchNcbiSummary(self, giIdCode, xmlPath=None):
        """
        """
        d = {}
        fetchobj = FetchNcbiXml(giIdCode, 'Nucleotide')
        if (xmlPath is not None):
            fetchobj.WriteNcbiXml(filename=xmlPath)

        d = fetchobj.ParseNcbiXmlData()
        if (self.__debug):
            for k, v in d.items():
                self.__lfh.write("+ReferenceSequenceUtils.__fetchNcbiSummary id content key %r value %r\n" % (k, v))
        return d

    def fetchNcbiTaxId(self, taxId):
        fetchobj = FetchNcbiXml(taxId, 'taxonomy')
        d = fetchobj.ParseNcbiXmlData()
        if (self.__debug):
            for (k, v) in d.items():
                self.__lfh.write("+ReferenceSequenceUtils.__fetchTaxIdNcbi key %s value %r\n" % (k, v))
        return d

    def __sortHitList(self, hitList, authTaxId=None, authAccessionCode=None):
        """ Add a sorting index to the dictionary of sequence correspondence hitLists obtained from
            the BLAST search.   The sort index is based on heurestic which includes sequence matching
            metrics and taxonomy data.

            Return hitList[]=d{}  -> with additional keys 'sort_order' & 'sort_metric'
                                     ordered by increasing (better matchint) score.

        """
        if not hitList:
            return

        authAncestorD = self.__taxUtils.getAncestors(authTaxId)
        if (self.__debug):
            self.__lfh.write("+ReferenceSequenceUtils.__sortHitList() length %d authTaxId %s authAncestorD %r\n" % (len(hitList), authTaxId, authAncestorD.items()))

        sList = []

        for i, hit in enumerate(hitList):
            score = (float(hit['identity']) - float(hit['gaps'])) * 100.0 / float(hit['query_length'])
            if self.__debug:
                self.__lfh.write("+ReferencSequenceUtils.__sortHitList() i %r  identity %r gaps %r query_length %r  score %r\n" %
                                 (i, hit['identity'], hit['gaps'], hit['query_length'], score))
            if hit['db_name'] == 'SP':
                sp_score = 1
            else:
                sp_score = 0
            if authAccessionCode and (hit['db_code'] == authAccessionCode or hit['db_accession'] == authAccessionCode):
                score += 1
            #
            targetAncestorD = {}
            if 'taxonomy_id' in hit:
                targetAncestorD = self.__taxUtils.getAncestors(hit['taxonomy_id'])
            #
            taxidMatch = 0
            if authAncestorD and targetAncestorD and 'id' in authAncestorD and 'id' in targetAncestorD:
                if authAncestorD['id'] == targetAncestorD['id']:
                    taxidMatch = 3
                elif 'p_id' in authAncestorD and authAncestorD['p_id'] == targetAncestorD['id']:
                    taxidMatch = 2
                elif 'p_id' in targetAncestorD and targetAncestorD['p_id'] == authAncestorD['id']:
                    taxidMatch = 2
                elif 'gp_id' in authAncestorD and authAncestorD['gp_id'] == targetAncestorD['id']:
                    taxidMatch = 1
                elif 'gp_id' in targetAncestorD and targetAncestorD['gp_id'] == authAncestorD['id']:
                    taxidMatch = 1
            #
            score = score * 4 + taxidMatch + sp_score

            if self.__debug:
                if authAncestorD and targetAncestorD and 'id' in authAncestorD and 'id' in targetAncestorD:
                    self.__lfh.write("+ReferencSequenceUtils.__sortHitList() taxidMatch %d  authAncestorD[id] %s targetAncestorD[id] %s final score %r\n" %
                                     (taxidMatch, authAncestorD['id'], targetAncestorD['id'], score))
            hit['sort_metric'] = score
            #
        #
        # hitList.sort(key=itemgetter('sort_metric'),reverse=True)
        hitList.sort(key=itemgetter('sort_metric'))
        #
        for ii, hit in enumerate(hitList, start=1):
            hit['sort_order'] = str(ii)

    def writeFasta(self, filePath, sequence, comment="myquery"):
        num_per_line = 60
        l = len(sequence) / num_per_line
        x = len(sequence) % num_per_line
        m = l
        if x:
            m = l + 1

        seq = '>' + str(comment).strip() + '\n'
        for i in range(m):
            n = num_per_line
            if i == l:
                n = x
            seq += sequence[i * num_per_line:i * num_per_line + n]
            if i != (m - 1):
                seq += '\n'
        try:
            ofh = open(filePath, 'w')
            ofh.write(seq)
            ofh.close()
            return True
        except:
            pass

        return False
