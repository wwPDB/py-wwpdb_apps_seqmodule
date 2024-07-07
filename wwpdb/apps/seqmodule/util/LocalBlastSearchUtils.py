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
#  3-Oct-2022  zf  update __runBlastSearch() & __fetchAuthProvidedRefSequence() methods for better handling author provided
#                  reference sequence cases.
#  6-Jan-2024  zf  add getPredefinedRefSequence(), __readPredefinedRefSequenceInfo() and __processPredefinedRefSequence() methods
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
import multiprocessing
import copy
import os
import sys
import time
import traceback
from operator import itemgetter

from wwpdb.apps.seqmodule.io.BlastPlusReader import BlastPlusReader
from wwpdb.apps.seqmodule.io.FetchSeqInfoUtils import FetchSeqInfoUtils
from wwpdb.apps.seqmodule.io.TaxonomyUtils import TaxonomyUtils
from wwpdb.apps.seqmodule.util.FetchReferenceSequenceUtils import FetchReferenceSequenceUtils
from wwpdb.apps.seqmodule.util.SequenceReferenceData import SequenceReferenceData
from rcsb.utils.multiproc.MultiProcUtil import MultiProcUtil
from wwpdb.io.locator.PathInfo import PathInfo
from wwpdb.utils.config.ConfigInfoApp import ConfigInfoAppCommon
from wwpdb.utils.dp.RcsbDpUtility import RcsbDpUtility


class LocalBlastSearchUtils(object):
    """Execute search service and assemble reference sequence data."""

    def __init__(self, siteId="WWPDB_DEPLOY_TEST", sessionPath=".", pathInfo=None, doRefSearchFlag=False, verbose=False, log=sys.stderr, ncbilock=None):
        self.__verbose = verbose
        self.__lfh = log
        self.__siteId = siteId
        self.__sessionPath = sessionPath
        self.__pI = pathInfo
        self.__doRefSearchFlag = doRefSearchFlag
        self.__ncbilock = ncbilock
        #
        self.__minRnaSequenceSearchLength = 50
        self.__maxHitsSearch = 250
        self.__maxHitsSave = 50
        self.__shortSequenceLengthLimit = 20
        self.__cleanUp = True
        self.__debug = False
        #
        self.__dataSetId = ""
        self.__entityD = {}
        self.__authRefIdList = []
        self.__authRefDataList = []
        self.__fullOneLetterSeq = []
        self.__seqType = ""
        self.__entityId = ""
        self.__entitySeq = ""
        self.__blastResultFile = ""
        #
        if not self.__pI:
            self.__pI = PathInfo(siteId=self.__siteId, sessionPath=self.__sessionPath, verbose=self.__verbose, log=self.__lfh)
        #
        self.__taxUtils = TaxonomyUtils(siteId=self.__siteId, verbose=self.__verbose, log=self.__lfh)
        self.__srd = SequenceReferenceData(self.__verbose, self.__lfh)
        #
        self.__cICommon = ConfigInfoAppCommon(self.__siteId)
        self.__seqDbPath = self.__cICommon.get_site_refdata_sequence_db_path()
        #
        self.__matchEntityAttribNameList = [
            "id",
            "fragment_id",
            "beg_seq_num",
            "end_seq_num",
            "sort_order",
            "sort_metric",
            "db_name",
            "db_code",
            "db_accession",
            "db_isoform",
            "match_length",
            "queryFrom",
            "queryTo",
            "hitFrom",
            "hitTo",
            "identity",
            "positive",
            "gaps",
            "alignLen",
            "query",
            "subject",
            "midline",
            "query_length",
            "name",
            "source_scientific",
            "source_common",
            "taxonomy_id",
            "gene",
            "synonyms",
            "comments",
            "keyword",
            "ec",
            "db_length",
            "db_description",
            "db_isoform_description",
        ]
        #

    def searchSeqReference(self, dataSetId=None, entityD=None, authRefList=()):
        """Reference sequence search process for the entity described in the input dictionary.
        Return a list of hitLists for matching reference sequences extracted from the BLAST search results
        """
        self.__dataSetId = dataSetId
        self.__entityD = entityD
        if authRefList:
            self.__authRefIdList = authRefList[0]
            self.__authRefDataList = authRefList[1]
        #
        if (not self.__dataSetId) or (not self.__entityD):
            return {}
        #
        autoMatchStatus, rD = self.__getExistingSearchResult()
        if rD:
            return autoMatchStatus, rD
        #
        autoMatchStatus, rD = self.__runBlastSearch()
        if rD:
            self.__writeSearchResult(rD)
            # self.__writeBlastSearchResultCifFile(rD)
        #
        return autoMatchStatus, rD

    def getPredefinedRefSequence(self, dataSetId=None, entityD=None):
        """ Get reference sequence information for the sequence builder application
        """
        self.__dataSetId = dataSetId
        self.__entityD = entityD
        if not self.__entityD:
            return {}
        #
        oneLetterCodeSeq = str(self.__entityD["SEQ_ENTITY_1_CAN"]).strip().upper()
        self.__fullOneLetterSeq = self.__srd.toList(oneLetterCodeSeq)
        if not self.__checkPartRange(len(self.__fullOneLetterSeq), self.__entityD["PART_LIST"]):
            return {}
        #
        refSeqInfoD, matchedRefList = self.__readPredefinedRefSequenceInfo()
        if (not refSeqInfoD) or (not matchedRefList):
            return {}
        #
        return self.__processPredefinedRefSequence(oneLetterCodeSeq, refSeqInfoD, matchedRefList)

    def __getExistingSearchResult(self):
        """Get prvious blast search result"""
        self.__entityId = self.__entityD["ENTITY_ID"]
        currSeq = self.__entityD["SEQ_ENTITY_1"]
        self.__entitySeq = currSeq.replace(" ", "").replace("\n", "").replace("\t", "")
        self.__blastResultFile = self.__pI.getFilePath(self.__dataSetId, contentType="seqdb-match", formatType="pic", fileSource="session", partNumber=self.__entityId)
        #
        if (not self.__doRefSearchFlag) and os.access(self.__blastResultFile, os.R_OK):
            if self.__verbose:
                self.__lfh.write(
                    "+LocalBlastSearchUtils.__getExistingSearchResult() found previous ReferenceSequenceFile %s for entity id %s\n" % (self.__blastResultFile, self.__entityId)
                )
            #
            pickleObj = {}
            try:
                fb = open(self.__blastResultFile, "rb")
                pickleObj = pickle.load(fb)
                fb.close()
            except:  # noqa: E722 pylint: disable=bare-except
                pickleObj = {}
            #
            if (
                ("sequence" not in pickleObj)
                or ("part_info" not in pickleObj)
                or ("ref_data" not in pickleObj)
                or (self.__entitySeq != pickleObj["sequence"])
                or (len(self.__entityD["PART_LIST"]) != len(pickleObj["part_info"]))
            ):
                return False, {}
            #
            for partNum, fD in enumerate(self.__entityD["PART_LIST"], start=1):
                if (
                    (partNum not in pickleObj["part_info"])
                    or (fD["SEQ_NUM_BEG"] != pickleObj["part_info"][partNum][0])
                    or (fD["SEQ_NUM_END"] != pickleObj["part_info"][partNum][1])
                ):
                    return False, {}
                #
            #
            if self.__verbose:
                self.__lfh.write("+LocalBlastSearchUtils.__getExistingSearchResult() using previous search result for entity id %s\n" % self.__entityId)
            #
            return self.__getReSortHitInfo(pickleObj["ref_data"])
        #
        return False, {}

    def __getReSortHitInfo(self, inD):
        """ """
        try:
            autoMatchList = []
            outD = {}
            for partNum, fD in enumerate(self.__entityD["PART_LIST"], start=1):
                partId = str(partNum)
                if partId not in inD:
                    continue
                #
                _start_index, hitList = self.__sortHitList(inD[partId], fD["SOURCE_TAXID"], False)
                outD[partId] = hitList
                if (len(hitList) == 1) and ("auto_match_status" in hitList[0]) and hitList[0]["auto_match_status"]:
                    autoMatchList.append(True)
                #
            #
            autoMatchStatus = False
            if (len(autoMatchList) > 0) and (len(autoMatchList) == len(self.__entityD["PART_LIST"])):
                autoMatchStatus = True
            #
            return autoMatchStatus, outD
        except:  # noqa: E722 pylint: disable=bare-except
            return False, {}
        #

    def __runBlastSearch(self):
        """Perform sequence search for each entity part in the input feature list described in the input dictionary."""
        oneLetterCodeSeq = str(self.__entityD["SEQ_ENTITY_1_CAN"]).strip().upper()
        #
        if "SEQ_ENTITY_1_SEARCH" in self.__entityD:
            oneLetterSearchSeq = str(self.__entityD["SEQ_ENTITY_1_SEARCH"]).strip().upper()
            # Use one letter code 'X' instead of multiple one letter codes defined in CCD for chromophore residue (DAOTHER-9059)
            if len(oneLetterSearchSeq) != len(oneLetterCodeSeq):
                oneLetterSearchSeq = oneLetterCodeSeq
            #
        else:
            # for backward compatible
            oneLetterSearchSeq = oneLetterCodeSeq
        #
        self.__fullOneLetterSeq = self.__srd.toList(oneLetterCodeSeq)
        self.__seqType = self.__entityD["POLYMER_LINKING_TYPE"]
        #
        if self.__verbose:
            self.__lfh.write("\n+LocalBlastSearchUtils.__runBlastSearch() STARTING with polymer type  %s sequence = %s\n" % (self.__seqType, oneLetterCodeSeq))
        #
        authRefD = {}
        if not self.__checkPartRange(len(self.__fullOneLetterSeq), self.__entityD["PART_LIST"]):
            return False, authRefD
        #
        foundPerfectMatch = True
        partDataList = []
        blastDataList = []
        for partNum, fD in enumerate(self.__entityD["PART_LIST"], start=1):
            seqPartType = fD["SEQ_PART_TYPE"]
            seqPartId = fD["SEQ_PART_ID"]
            if not str(seqPartType).lower().startswith("biol"):
                if self.__verbose:
                    self.__lfh.write("+LocalBlastSearchUtils.__runBlastSearch() skipping sequence entity %r part %r type %r\n" % (self.__entityId, seqPartId, seqPartType))
                #
                continue
            #
            seqNumBeg = fD["SEQ_NUM_BEG"]
            seqNumEnd = fD["SEQ_NUM_END"]
            if (int(seqNumBeg) < 1) or (int(seqNumEnd) > len(self.__fullOneLetterSeq)):
                if self.__verbose:
                    self.__lfh.write(
                        "+LocalBlastSearchUtils.__runBlastSearch() skipping sequence entity %r part %r type %r BegNum %r EndNum %r\n"
                        % (self.__entityId, seqPartId, seqPartType, seqNumBeg, seqNumEnd)
                    )
                #
                continue
            #
            # Skip blast search for poly-UNK/ALA sequence (DAOTHER-3639 & DAOTHER-8135)
            countUNK = 0
            countALA = 0
            countTotal = 0
            for oneLetterCode in self.__fullOneLetterSeq[(int(seqNumBeg) - 1) : int(seqNumEnd)]:
                if oneLetterCode == "X":
                    countUNK += 1
                elif oneLetterCode == "A":
                    countALA += 1
                #
                countTotal += 1
            #
            if ((10 * countUNK) > (9 * countTotal)) or (self.__seqType.startswith("polypeptide") and ((10 * countALA) > (9 * countTotal))):
                if self.__verbose:
                    self.__lfh.write("+LocalBlastSearchUtils.__runBlastSearch() skipping unknown/poly-ALA sequence part %r type %r BegNum %r EndNum %r\n" %
                                     (seqPartId, seqPartType, seqNumBeg, seqNumEnd))
                    self.__lfh.flush()
                #
                if self.__seqType.startswith("polypeptide") and ((10 * countALA) > (9 * countTotal)):
                    foundPerfectMatch = False
                #
                continue
            #
            taxId = fD["SOURCE_TAXID"]
            searchSequence = oneLetterSearchSeq[(int(seqNumBeg) - 1) : int(seqNumEnd)]
            #
            searchResultFileList = []
            xmlFilePath = os.path.join(self.__sessionPath, self.__dataSetId + "_blast-match_E" + self.__entityId + "_P" + str(partNum) + ".xml")
            if os.access(xmlFilePath, os.F_OK):
                os.remove(xmlFilePath)
            #
            searchResultFileList.append(xmlFilePath)
            if self.__seqType in ["polypeptide(L)", "polypeptide(D)", "polypeptide"]:
                xmlFilePath = os.path.join(self.__sessionPath, self.__dataSetId + "_blast-match_E" + self.__entityId + "_P" + str(partNum) + ".1.xml")
                #
                if os.access(xmlFilePath, os.F_OK):
                    os.remove(xmlFilePath)
                #
                searchResultFileList.append(xmlFilePath)
            #
            partDataList.append((str(partNum), searchSequence, seqNumBeg, seqNumEnd, taxId, searchResultFileList))
            #
            mutationList = []
            if ("MUTATION_DETAILS" in self.__entityD) and self.__entityD["MUTATION_DETAILS"]:
                for val in self.__entityD["MUTATION_DETAILS"].strip().upper().replace(",", " ").split(" "):
                    if (not val) or (len(val) < 3) or (not val[0].isalpha()) or (not val[len(val) - 1].isalpha()):
                        continue
                    #
                    hasDigit = False
                    allDigit = True
                    for i in range(1, len(val) - 1):
                        if val[i].isdigit():
                            hasDigit = True
                        else:
                            allDigit = False
                        #
                    #
                    if hasDigit and allDigit:
                        mutationList.append(val)
                    #
                #
            #
            autoMatchStatus, skipBlastSearch, authHitList = self.__fetchAuthProvidedRefSequence(searchSequence, seqNumBeg, seqNumEnd, str(partNum),
                                                                                                taxId, mutationList, len(self.__entityD["PART_LIST"]))
            if not autoMatchStatus:
                foundPerfectMatch = False
            #
            if authHitList:
                authRefD[str(partNum)] = authHitList
                # Skip blast search if found perfect author provided reference sequence match
                if autoMatchStatus or skipBlastSearch:
                    continue
                #
                # Skip blast search for short sequence
#               if len(searchSequence) <= self.__shortSequenceLengthLimit:
#                   continue
#               #
            else:
                foundPerfectMatch = False
            #
            blastDataList.append((str(partNum), searchSequence, seqNumBeg, seqNumEnd, taxId, searchResultFileList))
        #
        if not partDataList:
            return foundPerfectMatch, authRefD
        #
        rD = {}
        try:
            if blastDataList:
                numProc = int(multiprocessing.cpu_count() / 2)
                mpu = MultiProcUtil(verbose=True)
                mpu.set(workerObj=self, workerMethod="runMultiLocalBlasts")
                mpu.setWorkingDir(self.__sessionPath)
                mpu.setOptions({"ncbilock": self.__ncbilock})
                _ok, _failList, _retLists, _diagList = mpu.runMulti(dataList=blastDataList, numProc=numProc, numResults=1)
            #
            id_count = 1
            for partTup in partDataList:
                mergedhitList = []
                uniqueIdList = []
                for resultFile in partTup[5]:
                    if not os.access(resultFile, os.F_OK):
                        continue
                    #
                    returnHitList = self.__readBlastResultXmlFile(seqNumBeg=partTup[2], seqNumEnd=partTup[3], seqPartId=partTup[0],
                                                                  blastResultXmlPath=resultFile)
                    if len(returnHitList) == 0:
                        continue
                    #
                    for hitD in returnHitList:
                        refID = ""
                        if ("db_isoform" in hitD) and hitD["db_isoform"]:
                            refID = hitD["db_isoform"]
                        elif ("db_accession" in hitD) and hitD["db_accession"]:
                            refID = hitD["db_accession"]
                        #
                        if (refID == "") or (refID in uniqueIdList):
                            continue
                        #
                        uniqueIdList.append(refID)
                        mergedhitList.append(hitD)
                    #
                #
                hitList = self.__getHitList(mergedhitList, partTup[4])
                if hitList:
                    SeqAlignmentMap = self.__getRefSeqAlignments(hitList)
                    for hit in hitList:
                        if not hit["sort_order"] in SeqAlignmentMap:
                            continue
                        #
                        hit["alignment"] = SeqAlignmentMap[hit["sort_order"]][0]
                        hit["seq_tup_list"] = SeqAlignmentMap[hit["sort_order"]][1]
                        hit["statistics"] = (SeqAlignmentMap[hit["sort_order"]][2], SeqAlignmentMap[hit["sort_order"]][3], SeqAlignmentMap[hit["sort_order"]][4])
                        #
                        hit["id"] = id_count
                        id_count += 1
                        for key in self.__matchEntityAttribNameList:
                            if key not in hit:
                                hit[key] = ""
                            #
                        #
                        hit["seq_sim"] = 0.0
                        if ("alignLen" in hit) and hit["alignLen"]:
                            hit["seq_sim"] = float(hit["identity"]) / float(hit["alignLen"])
                        #
                        if partTup[0] in rD:
                            rD[partTup[0]].append(hit)
                        else:
                            rD[partTup[0]] = [hit]
                        #
                    #
                #
                # Merge author provided reference sequence match result
                #
                if (partTup[0] in authRefD) and authRefD[partTup[0]]:
                    if partTup[0] not in rD:
                        rD[partTup[0]] = authRefD[partTup[0]]
                    else:
                        rD[partTup[0]] = self.__mergeAuthorProvidedRefSeqMath(rD[partTup[0]], authRefD[partTup[0]], partTup[4])
                    #
                #
            #
            if self.__verbose:
                self.__lfh.write("+LocalBlastSearchUtils.__runBlastSearch() COMPLETED search for entity %r\n" % self.__entityId)
            #
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
                self.__lfh.write("+LocalBlastSearchUtils.__runBlastSearch() search for entity %r failed\n" % self.__entityId)
            #
        #
        return foundPerfectMatch, rD

    def __fetchAuthProvidedRefSequence(self, sequence, seqNumBeg, seqNumEnd, seqPartId, taxId, mutationList, partNumbers):
        """Fetch reference sequence based on author provided ref IDs"""
        if (not self.__authRefDataList) or (len(self.__authRefDataList) > partNumbers):
            return False, False, []
        #
        refData = []
        if partNumbers > 1:
            for aRefData in self.__authRefDataList:
                if aRefData[3] and aRefData[4] and (int(seqNumBeg) <= int(aRefData[3])) and (int(seqNumEnd) >= int(aRefData[4])):
                    refData = aRefData
                    break
                #
            #
        else:
            refData = self.__authRefDataList[0]
        #
        if not refData:
            return False, False, []
        #
        fetchUtil = FetchReferenceSequenceUtils(siteId=self.__siteId, seqReferenceData=self.__srd, verbose=self.__verbose, log=self.__lfh)
        autoMatchStatus, skipBlastSearch, hitD = fetchUtil.fetchReferenceSequenceWithSeqMatch(refData[0], refData[1], refData[2], taxId, sequence, seqNumBeg, mutationList)
        if hitD:
            hitD["beg_seq_num"] = seqNumBeg
            hitD["end_seq_num"] = seqNumEnd
            hitD["fragment_id"] = seqPartId
            #
            alignIndex, sTup3L, alignLength, seqSim, seqSimWithGaps, error = self.__getSeqAlignIndex(hitD)
            if not error:
                hitD["alignment"] = alignIndex
                hitD["seq_tup_list"] = sTup3L
                hitD["statistics"] = (alignLength, seqSim, seqSimWithGaps)
                hitD["author_provided_id"] = True
                if autoMatchStatus:
                    hitD["auto_match_status"] = True
                #
                return autoMatchStatus, skipBlastSearch, [hitD]
            else:
                self.__lfh.write("+LocalBlastSearchUtils.__fetchAuthProvidedRefSequence() get sequence alignment index error: %s\n" % error)
            #
        #
        return False, False, []

    def __writeSearchResult(self, rD):
        """Write out the current search result pickle file"""
        partInfoD = {}
        for partNum, fD in enumerate(self.__entityD["PART_LIST"], start=1):
            partInfoD[partNum] = (fD["SEQ_NUM_BEG"], fD["SEQ_NUM_END"])
        #
        pickleObj = {}
        pickleObj["sequence"] = self.__entitySeq
        pickleObj["ref_data"] = rD
        pickleObj["part_info"] = partInfoD
        #
        fb = open(self.__blastResultFile, "wb")
        pickle.dump(pickleObj, fb)
        fb.close()

    def __readPredefinedRefSequenceInfo(self):
        """ Read reference sequence information generated in UpdatePolymerEntityPartitions.seqBuilderResponder()
        """
        #
        # Read reference sequence information generated in UpdatePolymerEntityPartitions.seqBuilderResponder()
        #
        refSeqInfoD = {}
        refSeqPickleFile = os.path.join(self.__sessionPath, "fetchedRefSeqs.pic")
        if os.access(refSeqPickleFile, os.F_OK):
            try:
                fb = open(refSeqPickleFile, "rb")
                refSeqInfoD = pickle.load(fb)
                fb.close()
            except:  # noqa: E722 pylint: disable=bare-except
                refSeqInfoD = {}
            #
        #
        if not refSeqInfoD:
            return {}, []
        #
        # Read matched sequence information generated in UpdatePolymerEntityPartitions.seqBuilderResponder()
        #
        matchedRefList = []
        self.__entityId = self.__entityD["ENTITY_ID"]
        matchRefSeqPickleFile = os.path.join(self.__sessionPath, "Entity-" + self.__entityId + "-MatchedRefSeqs.pic")
        if os.access(matchRefSeqPickleFile, os.F_OK):
            try:
                fb = open(matchRefSeqPickleFile, "rb")
                matchedRefList = pickle.load(fb)
                fb.close()
            except:  # noqa: E722 pylint: disable=bare-except
                matchedRefList = []
            #
        #
        return refSeqInfoD, matchedRefList

    def __processPredefinedRefSequence(self, oneLetterCodeSeq, refSeqInfoD, matchedRefList):
        """ Create Auth/Seq sequences alignment
        """
        predefinedRefD = {}
        for partNum, fD in enumerate(self.__entityD["PART_LIST"], start=1):
            try:
                seqNumBeg = int(fD["SEQ_NUM_BEG"])
                seqNumEnd = int(fD["SEQ_NUM_END"])
                foundRefFlag = False
                for matchedRefD in matchedRefList:
                    if (matchedRefD["queryFrom"] < seqNumBeg) or (matchedRefD["queryTo"] > seqNumEnd):
                        continue
                    #
                    if matchedRefD["db_accession"] not in refSeqInfoD:
                        continue
                    #
                    authSeq = oneLetterCodeSeq[(seqNumBeg - 1) : seqNumEnd]
                    matchedAuthSeq = authSeq[(matchedRefD["queryFrom"] - seqNumBeg) : (matchedRefD["queryTo"] - seqNumBeg + 1)]
                    # matchedAuthSeq = oneLetterCodeSeq[(queryFrom - 1) : queryTo]
                    matchedRefSeq = refSeqInfoD[matchedRefD["db_accession"]]["sequence"][(matchedRefD["hitFrom"] - 1) : matchedRefD["hitTo"]]
                    if matchedAuthSeq != matchedRefSeq:
                        return {}
                    #
                    hitD = copy.deepcopy(refSeqInfoD[matchedRefD["db_accession"]])
                    hitD["query_length"] = str(len(authSeq))
                    hitD["identity"] = str(len(matchedAuthSeq))
                    hitD["positive"] = str(len(matchedAuthSeq))
                    hitD["gaps"] = "0"
                    hitD["query"] = matchedAuthSeq
                    hitD["queryFrom"] = str(matchedRefD["queryFrom"] - seqNumBeg + 1)
                    hitD["queryTo"] = str(matchedRefD["queryTo"] - seqNumBeg + 1)
                    hitD["subject"] = matchedRefSeq
                    hitD["hitFrom"] = str(matchedRefD["hitFrom"])
                    hitD["hitTo"] = str(matchedRefD["hitTo"])
                    hitD["alignLen"] = str(len(matchedAuthSeq))
                    hitD["match_length"] = str(len(matchedAuthSeq))
                    hitD["seq_sim"] = 1.0
                    hitD["sort_metric"] = "404"
                    hitD["sort_order"] = "1"
                    hitD["beg_seq_num"] = fD["SEQ_NUM_BEG"]
                    hitD["end_seq_num"] = fD["SEQ_NUM_END"]
                    hitD["fragment_id"] = str(partNum)
                    #
                    alignIndex, sTup3L, alignLength, seqSim, seqSimWithGaps, error = self.__getSeqAlignIndex(hitD)
                    if not error:
                        hitD["alignment"] = alignIndex
                        hitD["seq_tup_list"] = sTup3L
                        hitD["statistics"] = (alignLength, seqSim, seqSimWithGaps)
                        predefinedRefD[str(partNum)] = [hitD]
                    else:
                        self.__lfh.write("+LocalBlastSearchUtils.__processPredefinedRefSequence() get sequence alignment index error: %s\n" % error)
                        return {}
                    #
                    foundRefFlag = True
                    break
                #
                if not foundRefFlag:
                    return {}
                #
            except:  # noqa: E722 pylint: disable=bare-except
                if self.__verbose:
                    traceback.print_exc(file=self.__lfh)
                    self.__lfh.write("+LocalBlastSearchUtils.__processPredefinedRefSequence() for entity %r failed\n" % self.__entityId)
                #
                return {}
            #
        #
        if predefinedRefD:
            retD = {}
            retD[self.__entityId] = predefinedRefD
            #
            currSeq = self.__entityD["SEQ_ENTITY_1"]
            self.__entitySeq = currSeq.replace(" ", "").replace("\n", "").replace("\t", "")
            self.__blastResultFile = self.__pI.getFilePath(self.__dataSetId, contentType="seqdb-match", formatType="pic",
                                                           fileSource="session", partNumber=self.__entityId)
            self.__writeSearchResult(retD)
            #
            return retD
        #
        return {}

    # def __writeBlastSearchResultCifFile(self, rD):
    #     """Write out the current search result mmcif file for debuging purpose"""
    #     #
    #     mR = []
    #     for _partId, rList in rD.items():
    #         mR.append(rList)
    #     #
    #     fn = self.__pI.getReferenceSequenceFilePath(self.__dataSetId, entityId=self.__entityId, fileSource="session")
    #     #
    #     from wwpdb.apps.seqmodule.io.PdbxIoUtils import ReferenceSequenceIo

    #     rsio = ReferenceSequenceIo(verbose=self.__verbose, log=self.__lfh)
    #     rsio.writeMatchResults(self.__entityD, outFilePath=fn, matchResults=mR)

    def __checkPartRange(self, seqLength, partList):
        """Check part range definition"""
        if len(partList) < 1:
            return False
        #
        try:
            for i in range(0, len(partList)):
                seqNumBeg = int(partList[i]["SEQ_NUM_BEG"])
                seqNumEnd = int(partList[i]["SEQ_NUM_END"])
                if seqNumEnd < seqNumBeg:
                    return False
                if (i == 0) and (seqNumBeg != 1):
                    return False
                if (i == (len(partList) - 1)) and (seqNumEnd != seqLength):
                    return False
                if i > 0:
                    prevNumEnd = int(partList[i - 1]["SEQ_NUM_END"])
                    if seqNumBeg != (prevNumEnd + 1):
                        return False
                    #
                #
            #
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            return False
        #

    def runMultiLocalBlasts(self, dataList, procName, optionsD, workingDir):  # pylint: disable=unused-argument
        """Multiple blast search processing API"""
        rList = []
        ncbilock = optionsD.get("ncbilock", None)
        for tupL in dataList:
            if self.__verbose:
                self.__lfh.write(
                    "\n+LocalBlastSearchUtils.runMultiLocalBlasts() starting database search for entityId %r partId %r range begin %r end %r taxId %r sequence length = %d\n"
                    % (self.__entityId, tupL[0], tupL[2], tupL[3], tupL[4], len(tupL[1]))
                )
                self.__lfh.write("+LocalBlastSearchUtils.runMultiLocalBlasts() ncbilock %s\n" % ncbilock)
            #
            if self.__seqType in ["polypeptide(L)", "polypeptide(D)", "polypeptide"]:
                swissprotDBPath = os.path.join(self.__seqDbPath, "my_swissprot_all.pjs")
                if os.access(swissprotDBPath, os.F_OK):
                    self.__runSingleLocalBlast(oneLetterCodeSeq=tupL[1], searchResultFile=tupL[5][0], partNum=tupL[0], taxId=tupL[4], dbName="my_swissprot_all")
                    if (not self.__findPerfectMath(begSeqNum=tupL[2], endSeqNum=tupL[3], partId=tupL[0], taxId=tupL[4], searchResultFile=tupL[5][0])) and \
                       (len(tupL[1]) > self.__shortSequenceLengthLimit):
                        # only search whole uniprot for long sequence based on DAOTHER-7190
                        self.__runSingleLocalBlast(oneLetterCodeSeq=tupL[1], searchResultFile=tupL[5][1], partNum=tupL[0], taxId=tupL[4], dbName="my_uniprot_all")
                    #
                else:
                    self.__runSingleLocalBlast(oneLetterCodeSeq=tupL[1], searchResultFile=tupL[5][0], partNum=tupL[0], taxId=tupL[4], dbName="my_uniprot_all")
                #
            else:
                self.__runSingleLocalBlast(oneLetterCodeSeq=tupL[1], searchResultFile=tupL[5][0], partNum=tupL[0], taxId=tupL[4], dbName="my_ncbi_nt")
            #
            rList.append(tupL[0])
        #
        return rList, rList, []

    def __runSingleLocalBlast(self, oneLetterCodeSeq, searchResultFile, partNum, taxId, dbName):
        """Internal method to execute the sequence search according to input polymer type."""
        if self.__verbose:
            self.__lfh.write("+LocalBlastSearchUtils.__runSingleLocalBlast() STARTING %s Launch search for %s sequence = %s\n" % (self.__siteId, self.__seqType, oneLetterCodeSeq))
            self.__lfh.flush()
        #
        timeBegin = time.time()
        resultPath = os.path.abspath(searchResultFile)
        tmpPathAbs = os.path.abspath(self.__sessionPath)
        #
        if self.__seqType in ["polypeptide(L)", "polypeptide(D)", "polypeptide"]:
            dp = RcsbDpUtility(tmpPath=tmpPathAbs, siteId=self.__siteId, verbose=True)
            dp.addInput(name="one_letter_code_sequence", value=oneLetterCodeSeq)
            dp.addInput(name="db_name", value=dbName)
            dp.addInput(name="num_threads", value="4")
            dp.addInput(name="max_hits", value=self.__maxHitsSearch)
            if len(oneLetterCodeSeq) <= self.__shortSequenceLengthLimit:
                dp.addInput(name="short_seq", value=True)
                dp.addInput(name="evalue", value="0.02")
            #
            if taxId:
                dp.addInput(name="tax_id", value=taxId)
            #
            dp.op("seq-blastp")
            dp.exp(resultPath)
            if self.__cleanUp and not self.__debug:
                dp.cleanup()
            #
        # elif self.__seqType == 'polyribonucleotide':
        # for DAOTHER-6304
        elif self.__seqType in ["polydeoxyribonucleotide", "polyribonucleotide"]:
            if len(oneLetterCodeSeq) > self.__minRnaSequenceSearchLength:
                dp = RcsbDpUtility(tmpPath=tmpPathAbs, siteId=self.__siteId, verbose=True)
                dp.addInput(name="one_letter_code_sequence", value=oneLetterCodeSeq)
                dp.addInput(name="db_name", value="my_ncbi_nt")
                dp.addInput(name="num_threads", value="4")
                dp.addInput(name="max_hits", value=self.__maxHitsSearch)
                if taxId:
                    dp.addInput(name="tax_id", value=taxId)
                #
                dp.op("seq-blastn")
                dp.exp(resultPath)
                if self.__cleanUp and not self.__debug:
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

    def __findPerfectMath(self, begSeqNum="", endSeqNum="", partId="", taxId="", searchResultFile=None):
        """ Check if there is a perfect match
        """
        if (not begSeqNum) or (not endSeqNum) or (not partId) or (not taxId) or (not searchResultFile) or (not os.access(searchResultFile, os.F_OK)):
            return False
        #
        hitList = self.__readBlastResultXmlFile(seqNumBeg=begSeqNum, seqNumEnd=endSeqNum, seqPartId=partId, blastResultXmlPath=searchResultFile)
        for hitD in hitList:
            if ("query_length" not in hitD) or (not hitD["query_length"]) or ("identity" not in hitD) or (not hitD["identity"]) or \
               ("taxonomy_id" not in hitD) or (not hitD["taxonomy_id"]):
                continue
            #
            if (hitD["query_length"] == hitD["identity"]) and (hitD["taxonomy_id"] == taxId):
                return True
            #
        #
        return False

    def __readBlastResultXmlFile(self, seqNumBeg, seqNumEnd, seqPartId, blastResultXmlPath):
        """Read blast result"""
        hitList = []
        bpr = BlastPlusReader(verbose=self.__verbose, log=self.__lfh)
        if self.__seqType == "polyribonucleotide":
            bpr.setSequenceType(self.__seqType)
        #
        searchHitList = bpr.readFile(filePath=blastResultXmlPath)
        #
        fetchSeqUtil = FetchSeqInfoUtils(siteId=self.__siteId, seqReferenceData=self.__srd, verbose=self.__verbose, log=self.__lfh)
        if self.__seqType in ["polypeptide(L)", "polypeptide(D)", "polypeptide"]:
            idCodeList = []
            for hit in searchHitList:
                hit["non_isoform_score"] = 1
                if "db_accession" in hit and "db_name" in hit:
                    try:
                        start = int(str(hit["hitFrom"]))
                    except:  # noqa: E722 pylint: disable=bare-except
                        start = 0
                    #
                    try:
                        end = int(str(hit["hitTo"]))
                    except:  # noqa: E722 pylint: disable=bare-except
                        end = 0
                    #
                    if "db_isoform" in hit and (len(hit["db_isoform"]) > 0) and (hit["db_isoform"] not in [".", "?"]):
                        idCodeList.append((hit["db_isoform"], start, end))
                    else:
                        idCodeList.append((hit["db_accession"], start, end))
                    #
                #
            #
            # Incorporate non-overlapping additional details about each matching sequence entry.
            #   ( Still using the REST API to get the entry information )
            #
            if self.__verbose:
                self.__lfh.write("+LocalBlastSearchUtils.__readBlastResultXmlFile() Fetch sequence database entries for %d filtered reference idcodes\n" % len(idCodeList))
            #
            if len(idCodeList) > 0:
                try:
                    unpD = fetchSeqUtil.fetchUniProt(idTupleList=idCodeList)
                except:  # noqa: E722 pylint: disable=bare-except
                    traceback.print_exc(file=self.__lfh)
                    self.__lfh.flush()
                #
                for hit in searchHitList:
                    acId = None
                    isIso = False
                    if "db_isoform" in hit and (len(hit["db_isoform"]) > 0) and (hit["db_isoform"] not in [".", "?"]):
                        acId = hit["db_isoform"]
                        hit["non_isoform_score"] = 0
                        isIso = True
                    elif "db_accession" in hit:
                        acId = hit["db_accession"]
                    #
                    if acId is not None:
                        try:
                            start = int(str(hit["hitFrom"]))
                        except:  # noqa: E722 pylint: disable=bare-except
                            start = 0
                        #
                        try:
                            end = int(str(hit["hitTo"]))
                        except:  # noqa: E722 pylint: disable=bare-except
                            end = 0
                        #
                        dd = {}
                        if (acId, start, end) in unpD:
                            dd = unpD[(acId, start, end)]
                        elif isIso and ((hit["db_accession"], start, end) in unpD):
                            dd = unpD[(hit["db_accession"], start, end)]
                        #
                        if dd:
                            for (k, v) in dd.items():
                                if (k != "sequence") and ((k not in hit) or (v and (k in ("name", "source_scientific", "strain", "taxonomy_id", "gene")))):
                                    hit[k] = v
                                #
                            #
                        #
                        hitList.append(hit)
                    #
                #
            #
        # elif self.__seqType == 'polyribonucleotide':
        # for DAOTHER-6304
        elif self.__seqType in ["polydeoxyribonucleotide", "polyribonucleotide"]:
            hitList = searchHitList
            idCodeList = []
            for hit in hitList:
                hit["non_isoform_score"] = 1
                if "db_accession" in hit:
                    idCodeList.append(hit["db_accession"])
                #
            #
            idCodeList = list(set(idCodeList))
            #
            # Incorporate non-overlapping additional details about each matching sequence entry.
            #
            if len(idCodeList) > 0:
                giD = {}
                for idCode in idCodeList:
                    # Ensure do not max out request rate
                    if self.__ncbilock:
                        self.__ncbilock.waitnext()
                    #
                    giD[idCode] = fetchSeqUtil.fetchNcbiSummary(idCode)
                #
                for hit in hitList:
                    if "db_accession" in hit:
                        acId = hit["db_accession"]
                        if acId in giD:
                            dd = giD[acId]
                            for (k, v) in dd.items():
                                if k not in hit:
                                    hit[k] = v
                                #
                            #
                        #
                    #
                    if ("source_scientific" not in hit or len(hit["source_scientific"]) < 2) and "taxonomy_id" in hit:
                        # try to lookup the scientific name - Take the first hit --
                        snL = self.__taxUtils.lookUpSource(taxId=hit["taxonomy_id"])
                        hit["source_scientific"] = ""
                        if len(snL) > 0:
                            hit["source_scientific"] = snL[0]
                        #
                    #
                #
            #
        #
        if self.__verbose:
            self.__lfh.write("+LocalBlastSearchUtils.__readBlastResultXmlFile() hit list length is %d\n" % len(hitList))
        #
        for hit in hitList:
            hit["beg_seq_num"] = seqNumBeg
            hit["end_seq_num"] = seqNumEnd
            hit["fragment_id"] = seqPartId
        #
        return hitList

    def __getHitList(self, hitList, authTaxId):
        """ """
        if not hitList:
            return []
        #
        start_index, finalList = self.__sortHitList(hitList, authTaxId, True)
        cutoff_length = len(finalList) - start_index
        if cutoff_length > self.__maxHitsSave:
            cutoff_length = self.__maxHitsSave
        #
        if len(finalList) > cutoff_length:
            return finalList[(len(finalList) - cutoff_length) :]
        else:
            return finalList
        #

    def __sortHitList(self, hitList, authTaxId, fullScoreFlag):
        """Add a sorting index to the dictionary of sequence correspondence hitLists obtained from
        the BLAST search.   The sort index is based on heurestic which includes sequence matching
        metrics and taxonomy data.

        Return hitList[]=d{}  -> with additional keys 'sort_order' & 'sort_metric'
                                 ordered by increasing (better matchint) score.

        """
        authAncestorD = self.__taxUtils.getAncestors(authTaxId)
        if self.__debug:
            self.__lfh.write("+LocalBlastSearchUtils.__sortHitList() length %d authTaxId %s authAncestorD %r\n" % (len(hitList), authTaxId, authAncestorD.items()))
        #
        highest_identity_score_value = 0.0
        highest_identity_score_index = 0
        same_taxid_indice_list = []
        highest_identity_score_with_taxid_match = 0
        lowest_identity_score_with_taxid_match = 101
        lowest_identity_score_with_refid_match = 101
        has_author_provided_id = False
        for i, hit in enumerate(hitList):
            if "non_isoform_score" not in hit:
                hit["non_isoform_score"] = 1
            #
            hit["db_score"] = 0
            if fullScoreFlag or ("identity_score" not in hit):
                hit["identity_score"] = (float(hit["identity"]) - float(hit["gaps"])) * 100.0 / float(hit["query_length"])
                hit["sorting_score"] = hit["identity_score"]
                if hit["identity_score"] > highest_identity_score_value:
                    highest_identity_score_value = hit["identity_score"]
                    highest_identity_score_index = i
                #
                if self.__debug:
                    self.__lfh.write(
                        "+ReferencSequenceUtils.__sortHitList() i %r  identity %r gaps %r query_length %r  score %r\n"
                        % (i, hit["identity"], hit["gaps"], hit["query_length"], hit["identity_score"])
                    )
                    #
                if hit["db_name"] == "SP":
                    hit["db_score"] = 1
                #
            #
            if self.__authRefIdList and ((hit["db_code"].strip().upper() in self.__authRefIdList) or (hit["db_accession"].strip().upper() in
               self.__authRefIdList) or (hit["db_isoform"].strip().upper() in self.__authRefIdList)):
                if hit["identity_score"] < lowest_identity_score_with_refid_match:
                    lowest_identity_score_with_refid_match = hit["identity_score"]
                #
                hit["author_provided_id"] = True
                has_author_provided_id = True
                hit["code_score"] = 1
            else:
                hit["code_score"] = 0
            #
            targetAncestorD = {}
            taxonomy_id = ""
            if ("taxonomy_id" in hit) and hit["taxonomy_id"]:
                taxonomy_id = hit["taxonomy_id"]
                targetAncestorD = self.__taxUtils.getAncestors(hit["taxonomy_id"])
            #
            taxidMatch = 0
            if authAncestorD and targetAncestorD and "id" in authAncestorD and "id" in targetAncestorD:
                if authAncestorD["id"] == targetAncestorD["id"]:
                    taxidMatch = 3
                    same_taxid_indice_list.append(i)
                elif "p_id" in authAncestorD and authAncestorD["p_id"] == targetAncestorD["id"]:
                    taxidMatch = 2
                elif "p_id" in targetAncestorD and targetAncestorD["p_id"] == authAncestorD["id"]:
                    taxidMatch = 2
                elif "gp_id" in authAncestorD and authAncestorD["gp_id"] == targetAncestorD["id"]:
                    taxidMatch = 1
                elif "gp_id" in targetAncestorD and targetAncestorD["gp_id"] == authAncestorD["id"]:
                    taxidMatch = 1
                #
            #
            # Add the following modification per request by DAOTHER-8863
            #
            if (hit["db_score"] == 1) and (taxidMatch > 0):
                taxidMatch = 3
            #
            hit["taxid_match"] = taxidMatch
            if fullScoreFlag and (taxidMatch == 3):
                if hit["identity_score"] > highest_identity_score_with_taxid_match:
                    highest_identity_score_with_taxid_match = hit["identity_score"]
                #
                if hit["identity_score"] < lowest_identity_score_with_taxid_match:
                    lowest_identity_score_with_taxid_match = hit["identity_score"]
                #
            #
            hit["sort_metric"] = (hit["identity_score"] + hit["code_score"]) * 4 + taxidMatch + hit["db_score"]
            #
            if self.__debug:
                self.__lfh.write(
                    "+ReferencSequenceUtils.__sortHitList() hit i=%d db_code=%s authTaxId=%s taxonomy_id=%s taxidMatch=%d sp_score=%d score=%.2f\n"
                    % (i, hit["db_code"], authTaxId, taxonomy_id, taxidMatch, hit["db_score"], hit["sort_metric"])
                )
                if authAncestorD and targetAncestorD and "id" in authAncestorD and "id" in targetAncestorD:
                    self.__lfh.write(
                        "+ReferencSequenceUtils.__sortHitList() taxidMatch %d  authAncestorD[id] %s targetAncestorD[id] %s final score %r\n"
                        % (taxidMatch, authAncestorD["id"], targetAncestorD["id"], hit["sort_metric"])
                    )
                #
            #
        #
        cutoff_identity_score = 0
        if highest_identity_score_with_taxid_match > 0:
            cutoff_identity_score = 0.85 * highest_identity_score_with_taxid_match
            if lowest_identity_score_with_refid_match < cutoff_identity_score:
                cutoff_identity_score = lowest_identity_score_with_refid_match - 0.1
            #
            if lowest_identity_score_with_taxid_match < cutoff_identity_score:
                cutoff_identity_score = lowest_identity_score_with_taxid_match - 0.1
            #
        #
        # Improvement based on entries D_1000212253/5CVM, D_1000236156/6EEB with highest hit has conflicts
        if (len(same_taxid_indice_list) > 0) and (highest_identity_score_index not in same_taxid_indice_list) and \
           ("match_length" in hitList[highest_identity_score_index]) and (int(hitList[highest_identity_score_index]["identity"])
                                                                          < int(hitList[highest_identity_score_index]["match_length"])):
            for idx in same_taxid_indice_list:
                if "match_length" not in hitList[idx]:
                    continue
                #
                allowed_diffs = int(int(hitList[idx]["match_length"]) / 100) + 1
                if ((int(hitList[idx]["match_length"]) - int(hitList[idx]["identity"])) < allowed_diffs) and \
                   ((highest_identity_score_value - hitList[idx]["identity_score"]) < 1.0):
                    hitList[idx]["sorting_score"] = highest_identity_score_value
                #
            #
        #
        # The highest score hit is at the bottom of the list.
        #
        hitList.sort(key=itemgetter("sorting_score", "taxid_match", "code_score", "db_score", "non_isoform_score"))
        #
        # reorder list if found author provided reference id
        #
        if has_author_provided_id and (len(hitList) > 2):
            idx = -1
            for i, hit in enumerate(hitList):
                if ("author_provided_id" in hit) and hit["author_provided_id"]:
                    idx = i
                #
            #
            if (idx >= 0) and (idx < (len(hitList) - 2)):
                hitList.insert(-1, hitList.pop(idx))
            #
        #
        start_index = 0
        first = True
        for i, hit in enumerate(hitList):
            hit["sort_order"] = str(i + 1)
            if ("author_provided_id" in hit) and hit["author_provided_id"]:
                continue
            #
            if first and (hit["identity_score"] > cutoff_identity_score):
                start_index = i
                first = False
            #
        #
        if self.__debug:
            self.__lfh.write("ReferencSequenceUtils.__sortHitList() hitList=%d maxHitsSave=%d\n" % (len(hitList), self.__maxHitsSave))
            for hit in hitList:
                taxonomy_id = ""
                if "taxonomy_id" in hit:
                    taxonomy_id = hit["taxonomy_id"]
                #
                db_code = ""
                if "db_code" in hit:
                    db_code = hit["db_code"]
                #
                self.__lfh.write(
                    "+ReferencSequenceUtils.__sortHitList() final hit sort_order=%r sort_metric=%r db_code=%r authTaxId=%r taxonomy_id=%r\n"
                    % (hit["sort_order"], hit["sort_metric"], db_code, authTaxId, taxonomy_id)
                )
            self.__lfh.flush()
            #
        #
        return start_index, hitList

    def __getRefSeqAlignments(self, hitList):
        """Get Auth/Ref sequence alignments"""
        alignMap = {}
        try:
            numProc = int(multiprocessing.cpu_count() / 2)
            mpu = MultiProcUtil(verbose=True)
            mpu.set(workerObj=self, workerMethod="getMultiSeqAlignmentProcess")
            mpu.setWorkingDir(self.__sessionPath)
            _ok, _failList, retLists, _diagList = mpu.runMulti(dataList=hitList, numProc=numProc, numResults=1)
            #
            for tupList in retLists:
                for tup in tupList:
                    alignMap[tup[0]] = (tup[1], tup[2], tup[3], tup[4], tup[5])
                #
            #
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
            #
        #
        return alignMap

    def __mergeAuthorProvidedRefSeqMath(self, blastHitList, authHitList, authTaxId):
        """ Merge author provided reference sequence match with blast search match(es)
            See ticket DAOTHER-7903 the rule for displaying row(s)
        """
        mergedHitList = []
        for authD in authHitList:
            found = False
            for blastD in blastHitList:
                if ("db_accession" in authD) and ("db_isoform" in authD) and ("db_accession" in blastD) and ("db_isoform" in blastD) and \
                   (authD["db_accession"] == blastD["db_accession"]) and (authD["db_isoform"] == blastD["db_isoform"]):
                    if ("auto_match_status" in authD) and authD["auto_match_status"]:
                        blastD["auto_match_status"] = True
                    #
                    found = True
                    break
                #
            #
            if not found:
                mergedHitList.append(authD)
            #
        #
        if not mergedHitList:
            return blastHitList
        #
        mergedHitList.extend(blastHitList)
        _start_index, finalList = self.__sortHitList(mergedHitList, authTaxId, True)
        return finalList

    def getMultiSeqAlignmentProcess(self, dataList, procName, optionsD, workingDir):  # pylint: disable=unused-argument
        """Get Auth/Ref sequence alignments MultiProcUtil API"""
        rList = []
        eList = []
        for hitD in dataList:
            alignIndex, sTup3L, alignLength, seqSim, seqSimWithGaps, error = self.__getSeqAlignIndex(hitD)
            if error:
                eList.append((hitD["sort_order"], error))
            else:
                rList.append((hitD["sort_order"], alignIndex, sTup3L, alignLength, seqSim, seqSimWithGaps))
            #
        #
        return rList, rList, eList

    def __getSeqAlignIndex(self, hitD):
        """Get Auth/Ref sequence alignment"""
        querySeq = self.__srd.toList(hitD["query"])
        refSeq = self.__srd.toList(hitD["subject"])
        if len(querySeq) != len(refSeq):
            return [], [], 0, 0.0, 0.0, "Alignment length error"
        #
        hBegin = int(hitD["hitFrom"])
        hEnd = int(hitD["hitTo"])
        if hBegin < hEnd:
            sTup3L = self.__srd.cnv1To3ListIdx(hitD["subject"], hBegin, self.__entityD["POLYMER_TYPE_CODE"])
        else:
            sTup3L = self.__srd.cnv1To3ListIdx(hitD["subject"], hBegin, self.__entityD["POLYMER_TYPE_CODE"], indexStep=-1)
        #
        align_start_position = int(hitD["beg_seq_num"]) - 1
        align_end_position = align_start_position + (int(hitD["queryFrom"]) - int(hitD["beg_seq_num"]))
        #
        qIndexBegin = int(hitD["beg_seq_num"]) + int(hitD["queryFrom"]) - 2
        alignIndex = []
        for i in range(0, qIndexBegin):
            alignIndex.append([i, -1])
        #
        usedCount = 0
        hIndexBegin = 0
        for i in range(0, len(querySeq)):
            if (querySeq[i] == "-") and (refSeq[i] == "-"):
                continue
            #
            if (querySeq[i] != "-") and (querySeq[i] != self.__fullOneLetterSeq[qIndexBegin]) and (self.__fullOneLetterSeq[qIndexBegin] != "X"):
                return [], [], 0, 0.0, 0.0, "query sequence mismatch"
            #
            if (refSeq[i] != "-") and (refSeq[i] != sTup3L[hIndexBegin][4]):
                return [], [], 0, 0.0, 0.0, "hit sequence mismatch"
            #
            align_end_position += 1
            if querySeq[i] == "-":
                alignIndex.append([-1, hIndexBegin])
                hIndexBegin += 1
            elif refSeq[i] == "-":
                alignIndex.append([qIndexBegin, -1])
                qIndexBegin += 1
                usedCount += 1
            else:
                alignIndex.append([qIndexBegin, hIndexBegin])
                qIndexBegin += 1
                hIndexBegin += 1
                usedCount += 1
            #
        #
        align_end_position += (int(hitD["end_seq_num"]) - int(hitD["queryFrom"]) - usedCount) + 1
        #
        for i in range(qIndexBegin, len(self.__fullOneLetterSeq)):
            alignIndex.append([i, -1])
        #
        beg_pos, end_pos = self.__getAlignmentPositionIndex(alignIndex, int(hitD["beg_seq_num"]) - 1, int(hitD["end_seq_num"]) - 1)
        if (beg_pos < 0) or (end_pos < 0):
            return [], [], 0, 0.0, 0.0, "Alignment length error"
        #
        refSeq = []
        for sTup in sTup3L:
            refSeq.append(sTup[4])
        #
        self.__processTerminalMismatch(alignIndex, beg_pos, end_pos, refSeq)
        alignLength, seqSim, seqSimWithGaps = self.__getAlignmentStatistics(alignIndex, beg_pos, end_pos, refSeq)
        #
        return alignIndex, sTup3L, alignLength, seqSim, seqSimWithGaps, ""

    def __getAlignmentPositionIndex(self, alignIndex, beg_seq_num, end_seq_num):
        """Get the starting & ending aligned positions"""
        beg_pos = -1
        end_pos = -1
        for i in range(0, len(alignIndex)):
            if (alignIndex[i][0] != "") and (int(alignIndex[i][0]) == beg_seq_num):
                beg_pos = i
            #
            if (alignIndex[i][0] != "") and (int(alignIndex[i][0]) == end_seq_num):
                end_pos = i
            #
        #
        return beg_pos, end_pos

    def __processTerminalMismatch(self, alignIndex, beg_pos, end_pos, refSeq):
        """Re-align with previous residue (for N-terminal) or next residue (for C-terminal) if they are the same"""
        first = -1
        for i in range(beg_pos, end_pos):
            if (alignIndex[i][0] != -1) and (alignIndex[i][1] != -1):
                if self.__fullOneLetterSeq[alignIndex[i][0]] != refSeq[alignIndex[i][1]]:
                    first = i
                #
                break
            #
        #
        if (
            (first > beg_pos)
            and (alignIndex[first - 1][0] >= 0)
            and (alignIndex[first - 1][1] < 0)
            and (self.__fullOneLetterSeq[alignIndex[first - 1][0]] == refSeq[alignIndex[first][1]])
        ):
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
        if (
            (last < end_pos)
            and (alignIndex[last + 1][0] >= 0)
            and (alignIndex[last + 1][1] < 0)
            and (self.__fullOneLetterSeq[alignIndex[last + 1][0]] == refSeq[alignIndex[last][1]])
        ):
            alignIndex[last + 1][1] = alignIndex[last][1]
            alignIndex[last][1] = -1
        #

    def __getAlignmentStatistics(self, alignIndex, beg_pos, end_pos, refSeq):
        """Get alignment length and similarities"""
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
            return 0, 0.0, 0.0
        else:
            return alignLength, float(numMatch) / float(alignLength), float(numMatchGaps) / float(alignLength)
        #
