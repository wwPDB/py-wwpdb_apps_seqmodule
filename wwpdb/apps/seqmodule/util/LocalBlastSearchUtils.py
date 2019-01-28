##
# File:  LocalBlastSearchUtils.py
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
# 24-Aug-2018  zf
##
"""
Methods to manage local and sequence search services and package matching reference sequence data.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.09"

try:
    import cPickle as pickle
except ImportError:
    import pickle as pickle
#
import multiprocessing, os, sys, time, traceback
from operator import itemgetter

from wwpdb.apps.seqmodule.io.BlastPlusReader import BlastPlusReader
from wwpdb.apps.seqmodule.io.FetchSeqInfoUtils import fetchNcbiSummary,fetchUniProt
from wwpdb.apps.seqmodule.io.TaxonomyUtils import TaxonomyUtils
from wwpdb.apps.seqmodule.util.SequenceReferenceData import SequenceReferenceData
from rcsb.utils.multiproc.MultiProcUtil import MultiProcUtil
from wwpdb.io.locator.PathInfo import PathInfo
from wwpdb.utils.dp.RcsbDpUtility import RcsbDpUtility

class LocalBlastSearchUtils(object):
    """ Execute search service and assemble reference sequence data.
    """

    def __init__(self, siteId="WWPDB_DEPLOY_TEST", sessionPath=".", pathInfo=None, doRefSearchFlag=False, verbose=False, log=sys.stderr):
        self.__verbose = verbose
        self.__lfh = log
        self.__siteId = siteId
        self.__sessionPath = sessionPath
        self.__pI = pathInfo
        self.__doRefSearchFlag = doRefSearchFlag
        #
        self.__minRnaSequenceSearchLength = 50
        self.__maxHitsSearch = 100
        self.__maxHitsSave = 50
        self.__cleanUp = True
        self.__debug = False
        #
        self.__dataSetId = ''
        self.__entityD = {}
        self.__fullOneLetterSeq = []
        self.__seqType = ''
        self.__entityId = ''
        self.__entitySeq = ''
        self.__blastResultFile = ''
        #
        if not self.__pI:
            self.__pI = PathInfo(siteId=self.__siteId, sessionPath=self.__sessionPath, verbose=self.__verbose, log=self.__lfh)
        #
        self.__taxUtils = TaxonomyUtils(siteId=self.__siteId, verbose=self.__verbose, log=self.__lfh)
        self.__srd = SequenceReferenceData(self.__verbose, self.__lfh)
        #
        self.__matchEntityAttribNameList = ['id', 'fragment_id', 'beg_seq_num', 'end_seq_num', 'sort_order', 'sort_metric', 'db_name', 'db_code', \
                                            'db_accession', 'db_isoform', 'match_length', 'queryFrom', 'queryTo', 'hitFrom', 'hitTo', 'identity', \
                                            'positive', 'gaps', 'alignLen', 'query', 'subject', 'midline', 'query_length', 'name', \
                                            'source_scientific', 'source_common', 'taxonomy_id', 'gene', 'synonyms', 'comments', 'keyword', 'ec', \
                                            'db_length', 'db_description', 'db_isoform_description']
        #

    def searchSeqReference(self, dataSetId=None, entityD=None):
        """ Reference sequence search process for the entity described in the input dictionary.
            Return a list of hitLists for matching reference sequences extracted from the BLAST search results
        """
        self.__dataSetId = dataSetId
        self.__entityD = entityD
        if (not self.__dataSetId) or (not self.__entityD):
            return {}
        #
        rD = self.__getExistingSearchResult()
        if rD:
            return rD
        #
        rD = self.__runBlastSearch()
        if rD:
            self.__writeSearchResult(rD)
            #self.__writeBlastSearchResultCifFile(rD)
        #
        return rD

    def __getExistingSearchResult(self):
        """ Get prvious blast search result
        """
        self.__entityId = self.__entityD['ENTITY_ID']
        currSeq = self.__entityD['SEQ_ENTITY_1']
        self.__entitySeq = currSeq.replace(' ', '').replace('\n', '').replace('\t', '')
        self.__blastResultFile = self.__pI.getFilePath(self.__dataSetId, contentType='seqdb-match', formatType='pic', fileSource='session', \
                                                       partNumber=self.__entityId)
        #
        if (not self.__doRefSearchFlag) and os.access(self.__blastResultFile, os.R_OK):
            if (self.__verbose):
                self.__lfh.write("+LocalBlastSearchUtils.__getExistingSearchResult() found previous ReferenceSequenceFile %s for entity id %s\n" \
                              % (self.__blastResultFile, self.__entityId))
            #
            pickleObj = {}
            try:
                fb = open(self.__blastResultFile, 'rb')
                pickleObj = pickle.load(fb)
                fb.close()
            except:
                pickleObj = {}
            #
            if ('sequence' not in pickleObj) or ('part_info' not in pickleObj) or ('ref_data' not in pickleObj) or \
               (self.__entitySeq != pickleObj['sequence']) or (len(self.__entityD['PART_LIST']) != len(pickleObj['part_info'])):
                return {}
            #
            for partNum, fD in enumerate(self.__entityD['PART_LIST'], start=1):
                if (partNum not in pickleObj['part_info']) or (fD['SEQ_NUM_BEG'] != pickleObj['part_info'][partNum][0]) or \
                   (fD['SEQ_NUM_END'] != pickleObj['part_info'][partNum][1]):
                    return {}
                #
            #
            if (self.__verbose):
                self.__lfh.write("+LocalBlastSearchUtils.__getExistingSearchResult() using previous search result for entity id %s\n" % self.__entityId)
            #
            return pickleObj['ref_data']
        #
        return {}

    def __runBlastSearch(self):
        """ Perform sequence search for each entity part in the input feature list described in the input dictionary.
        """
        oneLetterCodeSeq = str(self.__entityD['SEQ_ENTITY_1_CAN']).strip().upper()
        self.__fullOneLetterSeq = self.__srd.toList(oneLetterCodeSeq)
        self.__seqType = self.__entityD['POLYMER_LINKING_TYPE']
        #
        if self.__verbose:
            self.__lfh.write("\n+LocalBlastSearchUtils.__runBlastSearch() STARTING with polymer type  %s sequence = %s\n" % (self.__seqType, oneLetterCodeSeq))
        #
        if not self.__checkPartRange(len(self.__fullOneLetterSeq), self.__entityD['PART_LIST']):
            return {}
        #
        partList = []
        for partNum, fD in enumerate(self.__entityD['PART_LIST'], start=1):
            seqPartType = fD['SEQ_PART_TYPE']
            seqPartId = fD['SEQ_PART_ID']
            if not str(seqPartType).lower().startswith('biol'):
                if self.__verbose:
                    self.__lfh.write("+LocalBlastSearchUtils.__runBlastSearch() skipping sequence entity %r part %r type %r\n" %
                                    (self.__entityId, seqPartId, seqPartType))
                #
                continue
            #
            seqNumBeg = fD['SEQ_NUM_BEG']
            seqNumEnd = fD['SEQ_NUM_END']
            if (int(seqNumBeg) < 1) or (int(seqNumEnd) > len(self.__fullOneLetterSeq)):
                if self.__verbose:
                    self.__lfh.write("+LocalBlastSearchUtils.__runBlastSearch() skipping sequence entity %r part %r type %r BegNum %r EndNum %r\n" % \
                                    (self.__entityId, seqPartId, seqPartType, seqNumBeg, seqNumEnd))
                #
                continue
            #
            # Skip blast search for UNK sequence DAOTHER-3639
            countUNK = 0
            countTotal = 0
            for oneLetterCode in self.__fullOneLetterSeq[(int(seqNumBeg) - 1):int(seqNumEnd)]:
                if oneLetterCode == "X":
                    countUNK += 1
                #
                countTotal += 1
            #
            if (10 * countUNK) > (9 * countTotal):
                 if self.__verbose:
                     self.__lfh.write("+LocalBlastSearchUtils.__runBlastSearch() skipping unknown sequence part %r type %r BegNum %r EndNum %r\n" % \
                                     (seqPartId, seqPartType, seqNumBeg, seqNumEnd))
                 #
                 continue
            #
            #
            taxId = fD['SOURCE_TAXID']
            sequence = oneLetterCodeSeq[(int(seqNumBeg) - 1):int(seqNumEnd)]
            #
            xmlFilePath = os.path.join(self.__sessionPath, self.__dataSetId + "_blast-match_E" + self.__entityId + "_P" + str(partNum) + ".xml")
            if os.access(xmlFilePath, os.F_OK):
                os.remove(xmlFilePath)
            #
            partList.append(( str(partNum), sequence, seqNumBeg, seqNumEnd, taxId, xmlFilePath ))
        #
        rD = {}
        if not partList:
            return rD
        #
        try:
            numProc = multiprocessing.cpu_count() / 2
            mpu = MultiProcUtil(verbose = True)
            mpu.set(workerObj = self, workerMethod = "runMultiLocalBlasts")
            mpu.setWorkingDir(self.__sessionPath)
            ok,failList,retLists,diagList = mpu.runMulti(dataList = partList, numProc = numProc, numResults = 1)

            id_count = 1
            for partTup in partList:
                if not os.access(partTup[5], os.F_OK):
                    continue
                #
                hitList = self.__readBlastResultXmlFile(seqNumBeg=partTup[2], seqNumEnd=partTup[3], seqPartId=partTup[0], authTaxId=partTup[4], \
                                                        blastResultXmlPath=partTup[5])
                #
                if hitList:
                    SeqAlignmentMap = self.__getRefSeqAlignments(hitList)
                    for hit in hitList:
                        if not hit['sort_order'] in SeqAlignmentMap:
                            continue
                        #
                        hit['alignment'] = SeqAlignmentMap[hit['sort_order']][0]
                        hit['seq_tup_list'] = SeqAlignmentMap[hit['sort_order']][1]
                        hit['statistics'] = ( SeqAlignmentMap[hit['sort_order']][2], SeqAlignmentMap[hit['sort_order']][3], SeqAlignmentMap[hit['sort_order']][4] )
                        #
                        hit['id'] = id_count
                        id_count += 1
                        for key in self.__matchEntityAttribNameList:
                            if not key in hit:
                                hit[key] = ''
                            #
                        #
                        hit['seq_sim'] = 0.0
                        if ('alignLen' in hit) and hit['alignLen']:
                            hit['seq_sim'] = float(hit['identity']) / float(hit['alignLen'])
                        #
                        if partTup[0] in rD:
                            rD[partTup[0]].append(hit)
                        else:
                            rD[partTup[0]] = [ hit ]
                        #
                    #
                #
            #
            if self.__verbose:
                self.__lfh.write("+LocalBlastSearchUtils.__runBlastSearch() COMPLETED search for entity %r\n" % self.__entityId)
        #
        except:
            if (self.__verbose):
                traceback.print_exc(file=self.__lfh)
                self.__lfh.write("+LocalBlastSearchUtils.__runBlastSearch() search for entity %r failed\n" % self.__entityId)
            #
        #
        return rD

    def __writeSearchResult(self, rD):
        """ Write out the current search result pickle file
        """
        partInfoD = {}
        for partNum, fD in enumerate(self.__entityD['PART_LIST'], start=1):
            partInfoD[partNum] = ( fD['SEQ_NUM_BEG'], fD['SEQ_NUM_END'] )
        #
        pickleObj = {}
        pickleObj['sequence'] = self.__entitySeq
        pickleObj['ref_data'] = rD
        pickleObj['part_info'] = partInfoD
        #
        fb = open(self.__blastResultFile, 'wb')
        pickle.dump(pickleObj, fb)
        fb.close()

    def __writeBlastSearchResultCifFile(self, rD):
        """ Write out the current search result mmcif file for debuging purpose
        """
        #
        mR = []
        for partId,rList in rD.items():
            mR.append(rList)
        #
        fn = self.__pI.getReferenceSequenceFilePath(self.__dataSetId, entityId=self.__entityId, fileSource='session')
        #
        from wwpdb.apps.seqmodule.io.PdbxIoUtils import ReferenceSequenceIo
        rsio = ReferenceSequenceIo(verbose=self.__verbose, log=self.__lfh)
        rsio.writeMatchResults(self.__entityD, outFilePath=fn, matchResults=mR)

    def __checkPartRange(self, seqLength, partList):
        """ Check part range definition
        """
        if len(partList) < 1:
            return False
        #
        try:
            for i in range(0, len(partList)):
                seqNumBeg = int(partList[i]['SEQ_NUM_BEG'])
                seqNumEnd = int(partList[i]['SEQ_NUM_END'])
                if seqNumEnd < seqNumBeg:
                    return False
                if (i == 0) and (seqNumBeg != 1):
                    return False
                if (i == (len(partList) - 1)) and (seqNumEnd != seqLength):
                    return False
                if i > 0: 
                    prevNumEnd = int(partList[i-1]['SEQ_NUM_END'])
                    if seqNumBeg != (prevNumEnd + 1):
                        return False
                    #
                #
            #
            return True
        except:
            return False
        #

    def runMultiLocalBlasts(self, dataList, procName, optionsD, workingDir):
        """ Miltiple blast search processing API
        """
        rList = []
        for tupL in dataList:
            if self.__verbose:
                self.__lfh.write("\n+LocalBlastSearchUtils.runMultiLocalBlasts() starting database search for entityId %r partId %r range begin %r end %r taxId %r sequence length = %d\n" \
                              % (self.__entityId, tupL[0], tupL[2], tupL[3], tupL[4], len(tupL[1])))
            #
            self.__runSingleLocalBlast(oneLetterCodeSeq=tupL[1], blastResultXmlPath=tupL[5], partNum=tupL[0], taxId=tupL[4])
            rList.append(tupL[0])
        #
        return rList,rList,[]

    def __runSingleLocalBlast(self, oneLetterCodeSeq, blastResultXmlPath, partNum, taxId):
        """ Internal method to execute the sequence search according to input polymer type.
        """
        if self.__verbose:
            self.__lfh.write("+LocalBlastSearchUtils.__runSingleLocalBlast() STARTING %s Launch search for %s sequence = %s\n" % (self.__siteId, self.__seqType, oneLetterCodeSeq))
            self.__lfh.flush()
        #
        timeBegin = time.time()
        resultPath = os.path.abspath(blastResultXmlPath)
        tmpPathAbs = os.path.abspath(self.__sessionPath)
        #
        if self.__seqType in ['polypeptide(L)', 'polypeptide(D)', 'polypeptide']:
            dp = RcsbDpUtility(tmpPath=tmpPathAbs, siteId=self.__siteId, verbose=True)
            dp.addInput(name="one_letter_code_sequence", value=oneLetterCodeSeq)
            dp.addInput(name="db_name", value="my_uniprot_all")
            dp.addInput(name="num_threads", value="4")
            dp.addInput(name="max_hits", value=self.__maxHitsSearch)
            dp.op("seq-blastp")
            dp.exp(resultPath)
            if (self.__cleanUp and not self.__debug):
                dp.cleanup()
            #
        elif self.__seqType == 'polyribonucleotide':
            if len(oneLetterCodeSeq) > self.__minRnaSequenceSearchLength:
                dp = RcsbDpUtility(tmpPath=tmpPathAbs, siteId=self.__siteId, verbose=True)
                dp.addInput(name="one_letter_code_sequence", value=oneLetterCodeSeq)
                dp.addInput(name="db_name", value="my_ncbi_nt")
                dp.addInput(name="num_threads", value="4")
                dp.addInput(name="max_hits", value=self.__maxHitsSearch)
                dp.op("seq-blastn")
                dp.exp(resultPath)
                if (self.__cleanUp and not self.__debug):
                    dp.cleanup()
                #
            #
        else:
            self.__lfh.write("+LocalBlastSearchUtils.__runSingleLocalBlast() Search failed for unknown type = %s\n" % self.__seqType)
        #
        timeEnd = time.time()
        if self.__verbose:
            self.__lfh.write("+LocalBlastSearchUtils.__runSingleLocalBlast() completed for entityId %r partId %r taxId %r\n" % (self.__entityId, partNum, taxId))
            self.__lfh.write("+LocalBlastSearchUtils.__runSingleLocalBlast() Search processing completed after %.2f seconds\n" % (timeEnd - timeBegin))
            self.__lfh.flush()
        #

    def __readBlastResultXmlFile(self, seqNumBeg, seqNumEnd, seqPartId, authTaxId, blastResultXmlPath):
        """ Read blast result
        """
        hitList = []
        bpr = BlastPlusReader(verbose=self.__verbose, log=self.__lfh)
        searchHitList = bpr.readFile(filePath=blastResultXmlPath)
        #
        if self.__seqType in ['polypeptide(L)', 'polypeptide(D)', 'polypeptide']:
            idCodeList = []
            for hit in searchHitList:
                if 'db_accession' in hit and 'db_name' in hit:
                    if 'db_isoform' in hit and (len(hit['db_isoform']) > 0) and (hit['db_isoform'] not in ['.', '?']):
                        idCodeList.append(hit['db_isoform'])
                    else:
                        idCodeList.append(hit['db_accession'])
                    #
                #
            #
            idCodeList = list(set(idCodeList))
            #
            # Incorporate non-overlapping additional details about each matching sequence entry.
            #   ( Still using the REST API to get the entry information )
            #
            if (self.__verbose):
                self.__lfh.write("+LocalBlastSearchUtils.__readBlastResultXmlFile() Fetch sequence database entries for %d filtered reference idcodes\n" % len(idCodeList))
            if len(idCodeList) > 0:
                unpD = fetchUniProt(siteId=self.__siteId, verbose=self.__verbose, log=self.__lfh, idCodeList=idCodeList)
                for hit in searchHitList:
                    acId = None
                    isIso = False
                    if 'db_isoform' in hit and (len(hit['db_isoform']) > 0) and (hit['db_isoform'] not in ['.', '?']):
                        acId = hit['db_isoform']
                        isIso = True
                    elif 'db_accession' in hit:
                        acId = hit['db_accession']
                    #
                    if acId is not None:
                        dd = {}
                        if acId in unpD:
                            dd = unpD[acId]
                        elif isIso and hit['db_accession'] in unpD:
                            dd = unpD[hit['db_accession']]
                        #
                        for (k, v) in dd.items():
                            if (k != "sequence") and (k not in hit):
                                hit[k] = v
                            #
                        #
                        hitList.append(hit)
                    #
                #
            #
        elif self.__seqType == 'polyribonucleotide':
            hitList = searchHitList
            idCodeList = []
            for hit in hitList:
                if 'db_accession' in hit:
                    idCodeList.append(hit['db_accession'])
                #
            #
            idCodeList = list(set(idCodeList))
            #
            # Incorporate non-overlapping additional details about each matching sequence entry.
            #
            if len(idCodeList) > 0:
                giD = {}
                for idCode in idCodeList:
                    giD[idCode] = fetchNcbiSummary(idCode)
                #
                for hit in hitList:
                    if 'db_accession' in hit:
                        acId = hit['db_accession']
                        if acId in giD:
                            dd = giD[acId]
                            for (k, v) in dd.items():
                                if k not in hit:
                                    hit[k] = v
                                #
                            #
                        #
                    #
                    if (('source_scientific' not in hit or len(hit['source_scientific']) < 2) and 'taxonomy_id' in hit):
                        # try to lookup the scientific name - Take the first hit --
                        snL = self.__taxUtils.lookUpSource(taxId=hit['taxonomy_id'])
                        hit['source_scientific'] = ''
                        if len(snL) > 0:
                            hit['source_scientific'] = snL[0]
                        #
                    #
                #
            #
        #
        if (self.__verbose):
            self.__lfh.write("+LocalBlastSearchUtils.__readBlastResultXmlFile() hit list length is %d\n" % len(hitList))
        #
        for hit in hitList:
            hit['beg_seq_num'] = seqNumBeg
            hit['end_seq_num'] = seqNumEnd
            hit['fragment_id'] = seqPartId
        #
        return self.__sortHitList(hitList, authTaxId=authTaxId)

    def __sortHitList(self, hitList, authTaxId=None, authAccessionCode=None):
        """ Add a sorting index to the dictionary of sequence correspondence hitLists obtained from
            the BLAST search.   The sort index is based on heurestic which includes sequence matching
            metrics and taxonomy data.

            Return hitList[]=d{}  -> with additional keys 'sort_order' & 'sort_metric'
                                     ordered by increasing (better matchint) score.

        """
        if not hitList:
            return []
        #
        authAncestorD = self.__taxUtils.getAncestors(authTaxId)
        if (self.__debug):
            self.__lfh.write("+LocalBlastSearchUtils.__sortHitList() length %d authTaxId %s authAncestorD %r\n" % (len(hitList), authTaxId, authAncestorD.items()))
        #
        scoreList = []
        highest_identity_score = 0
        for i, hit in enumerate(hitList):
            identity_score = (float(hit['identity']) - float(hit['gaps'])) * 100.0 / float(hit['query_length'])
            if self.__debug:
                self.__lfh.write("+ReferencSequenceUtils.__sortHitList() i %r  identity %r gaps %r query_length %r  score %r\n" % \
                                 (i, hit['identity'], hit['gaps'], hit['query_length'], identity_score))
                #
            if hit['db_name'] == 'SP':
                sp_score = 1
            else:
                sp_score = 0
            #
            if authAccessionCode and (hit['db_code'] == authAccessionCode or hit['db_accession'] == authAccessionCode):
                db_code_score = 1
            else:
                db_code_score = 0
            #
            targetAncestorD = {}
            taxonomy_id = ''
            if 'taxonomy_id' in hit:
                taxonomy_id = hit['taxonomy_id']
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
            #
            if (i >= self.__maxHitsSave) and taxidMatch == 0:
                continue
            #
            if taxidMatch == 3:
                if identity_score > highest_identity_score:
                    highest_identity_score = identity_score
                #
            #
            final_score = (identity_score + db_code_score) * 4 + taxidMatch + sp_score
            scoreList.append( { "idx" : i, "identity_score" : identity_score, "sort_metric" : final_score, "taxid_match" : taxidMatch, \
                                "code_score" : db_code_score, "db_score" : sp_score } )
            #
            if self.__debug:
                self.__lfh.write("+ReferencSequenceUtils.__sortHitList() hit i=%d db_code=%s authTaxId=%s taxonomy_id=%s taxidMatch=%d sp_score=%d score=%.2f\n" % \
                        (i, hit['db_code'], authTaxId, taxonomy_id, taxidMatch, sp_score, final_score))
                if authAncestorD and targetAncestorD and 'id' in authAncestorD and 'id' in targetAncestorD:
                    self.__lfh.write("+ReferencSequenceUtils.__sortHitList() taxidMatch %d  authAncestorD[id] %s targetAncestorD[id] %s final score %r\n" % \
                                     (taxidMatch, authAncestorD['id'], targetAncestorD['id'], final_score))
                #
            #
        #
        sortScoreList = []
        for scoreD in scoreList:
            if highest_identity_score > 0:
                if (scoreD["taxid_match"] == 0) or (scoreD["identity_score"] < (0.85 * highest_identity_score)):
                    continue
                #
            #
            sortScoreList.append(scoreD)
        #
        #sortScoreList.sort(key=itemgetter('taxid_match', 'code_score', 'identity_score', 'code_score', 'db_score'),reverse=True)
        sortScoreList.sort(key=itemgetter('taxid_match', 'code_score', 'identity_score', 'code_score', 'db_score'))
        #
        finalList = []
        for i,scoreD in enumerate(sortScoreList, start=1):
            hitList[scoreD["idx"]]["sort_metric"] = scoreD["sort_metric"]
            hitList[scoreD["idx"]]["sort_order"] = str(i)
            finalList.append(hitList[scoreD["idx"]])
        #
        if len(finalList) > self.__maxHitsSave:
            return finalList[:self.__maxHitsSave]
        else:
            return finalList
        #

    def __getRefSeqAlignments(self, hitList):
        """ Get Auth/Ref sequence alignments
        """
        alignMap = {}
        try:
            numProc = multiprocessing.cpu_count() / 2
            mpu = MultiProcUtil(verbose = True)
            mpu.set(workerObj = self, workerMethod = 'getMultiSeqAlignmentProcess')
            mpu.setWorkingDir(self.__sessionPath)
            ok,failList,retLists,diagList = mpu.runMulti(dataList = hitList, numProc = numProc, numResults = 1)
            #
            for tupList in retLists:
                for tup in tupList:
                    alignMap[tup[0]] = ( tup[1], tup[2], tup[3], tup[4], tup[5] )
                #
            #
        except:
            if (self.__verbose):
                traceback.print_exc(file=self.__lfh)
            #
        #
        return alignMap

    def getMultiSeqAlignmentProcess(self, dataList, procName, optionsD, workingDir):
        """ Get Auth/Ref sequence alignments MultiProcUtil API
        """
        rList = []
        eList = []
        for hitD in dataList:
            alignIndex,sTup3L,alignLength,seqSim,seqSimWithGaps,error = self.__getSeqAlignIndex(hitD)
            if error:
                eList.append((hitD['sort_order'], error))
            else:
                rList.append(( hitD['sort_order'], alignIndex, sTup3L, alignLength, seqSim, seqSimWithGaps ))
            #
        #
        return rList,rList,eList

    def __getSeqAlignIndex(self, hitD):
        """ Get Auth/Ref sequence alignment
        """
        querySeq = self.__srd.toList(hitD['query'])
        refSeq = self.__srd.toList(hitD['subject'])
        if len(querySeq) != len(refSeq):
            return [],[],0,0.0,0.0,'Alignment length error'
        #
        hBegin = int(hitD['hitFrom'])
        hEnd = int(hitD['hitTo'])
        if hBegin < hEnd:
            sTup3L = self.__srd.cnv1To3ListIdx(hitD['subject'], hBegin, self.__entityD['POLYMER_TYPE_CODE'])
        else:
            sTup3L = self.__srd.cnv1To3ListIdx(hitD['subject'], hBegin, self.__entityD['POLYMER_TYPE_CODE'], indexStep=-1)
        #
        align_start_position = int(hitD['beg_seq_num']) - 1
        align_end_position = align_start_position + (int(hitD['queryFrom']) - int(hitD['beg_seq_num']))
        #
        qIndexBegin = int(hitD['beg_seq_num']) + int(hitD['queryFrom']) - 2
        alignIndex = []
        for i in range(0, qIndexBegin):
            alignIndex.append([ i, -1 ])
        #
        usedCount = 0
        hIndexBegin = 0
        for i in range(0, len(querySeq)):
            if (querySeq[i] == '-') and (refSeq[i] == '-'):
                continue
            #
            if (querySeq[i] != '-') and (querySeq[i] != self.__fullOneLetterSeq[qIndexBegin]):
                return [],[],0,0.0,0.0,'query sequence mismatch'
            #
            if (refSeq[i] != '-') and (refSeq[i] != sTup3L[hIndexBegin][4]):
                return [],[],0,0.0,0.0,'hit sequence mismatch'
            #
            align_end_position += 1
            if querySeq[i] == '-':
                alignIndex.append([ -1, hIndexBegin ])
                hIndexBegin += 1
            elif refSeq[i] == '-':
                alignIndex.append([ qIndexBegin, -1 ])
                qIndexBegin += 1
                usedCount += 1
            else:
                alignIndex.append([ qIndexBegin, hIndexBegin ])
                qIndexBegin += 1
                hIndexBegin += 1
                usedCount += 1
            #
        #
        align_end_position += (int(hitD['end_seq_num']) - int(hitD['queryFrom']) - usedCount) + 1
        #
        for i in range(qIndexBegin, len(self.__fullOneLetterSeq)):
            alignIndex.append([ i, -1 ])
        #
        beg_pos,end_pos = self.__getAlignmentPositionIndex(alignIndex, int(hitD['beg_seq_num']) - 1, int(hitD['end_seq_num']) - 1)
        if (beg_pos < 0) or (end_pos < 0):
            return [],[],0,0.0,0.0,'Alignment length error'
        #
        refSeq = []
        for sTup in sTup3L:
            refSeq.append(sTup[4])
        #
        self.__processTerminalMismatch(alignIndex, beg_pos, end_pos, refSeq)
        alignLength,seqSim,seqSimWithGaps = self.__getAlignmentStatistics(alignIndex, beg_pos, end_pos, refSeq)
        #
        return alignIndex,sTup3L,alignLength,seqSim,seqSimWithGaps,''

    def __getAlignmentPositionIndex(self, alignIndex, beg_seq_num, end_seq_num):
        """ Get the starting & ending aligned positions
        """
        beg_pos = -1
        end_pos = -1
        for i in range(0, len(alignIndex)):
            if (alignIndex[i][0] != '') and (int(alignIndex[i][0]) == beg_seq_num):
                beg_pos = i
            #
            if (alignIndex[i][0] != '') and (int(alignIndex[i][0]) == end_seq_num):
                end_pos = i
            #
        #
        return beg_pos,end_pos

    def __processTerminalMismatch(self, alignIndex, beg_pos, end_pos, refSeq):
        """ Re-align with previous residue (for N-terminal) or next residue (for C-terminal) if they are the same
        """
        first = -1
        for i in range(beg_pos, end_pos):
            if (alignIndex[i][0] != -1) and (alignIndex[i][1] != -1):
                if self.__fullOneLetterSeq[alignIndex[i][0]] != refSeq[alignIndex[i][1]]:
                    first = i
                #
                break
            #
        #
        if (first > beg_pos) and (alignIndex[first - 1][0] >= 0) and (alignIndex[first - 1][1] < 0) and \
           (self.__fullOneLetterSeq[alignIndex[first - 1][0]] == refSeq[alignIndex[first][1]]):
            alignIndex[first - 1][1] = alignIndex[first][1]
            alignIndex[first][1] = -1
        #
        last = end_pos + 1
        for i in range(end_pos, beg_pos, -1):
            if (alignIndex[i][0] != -1) and (alignIndex[i][1] != -1):
                if self.__fullOneLetterSeq[alignIndex[i][0]] != refSeq[alignIndex[i][1]]:
                    last = i
                #
                break
            #
        #
        if (last < end_pos) and (alignIndex[last + 1][0] >= 0) and (alignIndex[last + 1][1] < 0) and \
           (self.__fullOneLetterSeq[alignIndex[last + 1][0]] == refSeq[alignIndex[last][1]]):
            alignIndex[last + 1][1] = alignIndex[last][1];
            alignIndex[last][1] = -1;
        #

    def __getAlignmentStatistics(self, alignIndex, beg_pos, end_pos, refSeq):
        """ Get alignment length and similarities
        """
        alignLength = 0
        numMatch = 0
        numMatchGaps = 0
        for i in range(beg_pos, end_pos + 1):
            if (alignIndex[i][0] < 0) and (alignIndex[i][1] < 0):
                continue
            #
            alignLength += 1
            if alignIndex[i][0] < 0:
                continue
            elif alignIndex[i][1] < 0:
                numMatchGaps += 1
            elif self.__fullOneLetterSeq[alignIndex[i][0]] == refSeq[alignIndex[i][1]]:
                numMatch += 1
                numMatchGaps += 1
            #
        #
        if alignLength == 0:
            return 0,0.0,0.0
        else:
            return alignLength,float(numMatch) / float(alignLength),float(numMatchGaps) / float(alignLength)
        #
