##
# File:  GetSameSeqAnnotation.py
# Date:  23-Mar-2018
# Updates:
#
# 19-Oct-2022 zf add author provided reference sequence information
#
##
"""

This software was developed as part of the World Wide Protein Data Bank
Common Deposition and Annotation System Project

Copyright (c) 2018 wwPDB

This software is provided under a Creative Commons Attribution 3.0 Unported
License described at http://creativecommons.org/licenses/by/3.0/.

"""
__docformat__ = "restructuredtext en"
__author__ = "Zukang Feng"
__email__ = "zfeng@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.07"

import os
import sys
import traceback
from operator import itemgetter

from wwpdb.io.locator.PathInfo import PathInfo
from wwpdb.io.file.mmCIFUtil import mmCIFUtil
from wwpdb.apps.seqmodule.align.AlignmentToolUtils import codeIndex
from wwpdb.apps.seqmodule.io.FetchSeqInfoUtils import FetchSeqInfoUtils
from wwpdb.apps.seqmodule.util.SequenceReferenceData import SequenceReferenceData


class GetSameSeqAnnotation(object):
    def __init__(self, siteId=None, sessionPath=None, pathInfo=None, verbose=False, log=sys.stderr):
        """ """
        self.__siteId = siteId
        self.__sessionPath = sessionPath
        self.__pI = pathInfo
        self.__verbose = verbose
        self.__lfh = log
        self.__oneLetterCodeSeqList = []
        self.__threeLetterCodeSeqList = []
        #
        if not self.__pI:
            self.__pI = PathInfo(siteId=self.__siteId, sessionPath=self.__sessionPath, verbose=False, log=self.__lfh)
        #
        self.__srd = SequenceReferenceData(verbose=self.__verbose, log=self.__lfh)
        self.__gapSymbol = self.__srd.getGapSymbol()

    def setEntitySeq(self, seq="", polyTypeCode="AA"):
        """Set three letter sequence from input one letter seq"""
        if (not seq) or (len(seq) < 1):
            return
        #
        (self.__oneLetterCodeSeqList, self.__threeLetterCodeSeqList) = self.__srd.cnv1ListPlus3List(seq, polyTypeCode)

    def getSeqAnnotationFromAssignFile(self, retList=None, authPartsTaxIdInfoList=None, authRefList=None, includeSelfRef=False):
        """retList[][0]: depID, retList[][1]: entityID, retList[][2]: pdbID, retList[][3]: AnnInitial, retList[][4]: statusCode,
        retList[][5]: date_begin_processing
        """
        if (not retList) or (not authPartsTaxIdInfoList):
            return {}
        #
        if len(authPartsTaxIdInfoList) == 1:
            authPartsTaxIdInfoList[0][1] = "1"
            authPartsTaxIdInfoList[0][2] = len(self.__threeLetterCodeSeqList)
        #
        assignFileList = []
        for entityTup in retList:
            assignFile = self.__pI.getFilePath(
                dataSetId=entityTup[0], wfInstanceId=None, contentType="seq-assign", formatType="pdbx", fileSource="archive", versionId="latest", partNumber="1"
            )
            if (not assignFile) or (not os.access(assignFile, os.F_OK)):
                continue
            #
            statinfo = os.stat(assignFile)
            #
            entityTup.append(assignFile)
            entityTup.append(statinfo.st_mtime)
            assignFileList.append(entityTup)
        #
        if not assignFileList:
            return {}
        #
        if len(assignFileList) > 1:
            assignFileList.sort(key=itemgetter(7), reverse=True)
        #
        autoMatchList = []
        otherMatchList = []
        for entityTup in assignFileList:
            annInfoMap = self.__getSeqAnnotationFromAssignFile(entityTup, authPartsTaxIdInfoList, authRefList, includeSelfRef)
            if not annInfoMap:
                continue
            #
            if ("auto_match_status" in annInfoMap) and annInfoMap["auto_match_status"]:
                autoMatchList.append(annInfoMap)
            else:
                otherMatchList.append(annInfoMap)
            #
        #
        if len(autoMatchList) > 0:
            return autoMatchList[0]
        elif len(otherMatchList) > 0:
            return otherMatchList[0]
        else:
            return {}
        #

    def testGetSeqAnnotationFromAssignFile(self, assignFile, entityId, authPartsTaxIdInfoList):
        """ """
        try:
            cifObj = mmCIFUtil(filePath=assignFile)
            seqAnntationInfoMap, partsTaxIdInfo = self.__getEntityDetailsMap(cifObj, entityId, authPartsTaxIdInfoList)
            self.__lfh.write("seqAnntationInfoMap=%d\n" % len(seqAnntationInfoMap))
            self.__lfh.write("partsTaxIdInfo=%d\n" % len(partsTaxIdInfo))
            for partId, infoD in seqAnntationInfoMap.items():
                self.__lfh.write("seqAnntationInfoMap.partId=%r\n" % partId)
                self.__lfh.write("seqAnntationInfoMap.infoD=%r\n" % infoD)
            #
            alignmentMap = self.__getAlignmentMap(cifObj, entityId, partsTaxIdInfo)
            for partId, alignmentTup in alignmentMap.items():
                self.__lfh.write("alignmentMap.partId=%r\n" % partId)
                self.__lfh.write("alignmentMap.alignment:\n")
                for alignTup in alignmentTup[0]:
                    print(alignTup)
                #
                print(alignmentTup[1])
                self.__lfh.write("alignmentMap.seq_tup_list:\n")
                for seqTup in alignmentTup[2]:
                    print(seqTup)
                #
            #
            dbRefMap = self.__getDbRefMap(cifObj, entityId)
            for partId, infoD in dbRefMap.items():
                self.__lfh.write("dbRefMap.partId=%r\n" % partId)
                self.__lfh.write("dbRefMap.partIdinfoD=%r\n" % infoD)
            #
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
            #
        #

    def __getSeqAnnotationFromAssignFile(self, entityInfo, authPartsTaxIdInfoList, authRefList, includeSelfRef):
        """ """
        try:
            cifObj = mmCIFUtil(filePath=entityInfo[6])
            seqAnntationInfoMap, partsTaxIdInfo = self.__getEntityDetailsMap(cifObj, entityInfo[1], authPartsTaxIdInfoList)
            if not seqAnntationInfoMap:
                return seqAnntationInfoMap
            #
            alignmentMap = self.__getAlignmentMap(cifObj, entityInfo[1], partsTaxIdInfo)
            #           if not alignmentMap:
            #               return {}
            #           #
            dbRefMap = self.__getDbRefMap(cifObj, entityInfo[1], authRefList)
            selfRefmap = self.__getSelfRefmap(cifObj, entityInfo[1], includeSelfRef)
            auto_match_status = True
            for partId, infoD in seqAnntationInfoMap.items():
                if partId in selfRefmap:
                    infoD["selfref"] = True
                    continue
                #
                if partId not in dbRefMap:
                    return {}
                #
                if "auto_match_status" not in dbRefMap[partId]:
                    auto_match_status = False
                #
                for k, v in dbRefMap[partId].items():
                    infoD[k] = v
                #
                if partId not in alignmentMap:
                    return {}
                #
                infoD["alignment"] = alignmentMap[partId][0]
                infoD["statistics"] = alignmentMap[partId][1]
                infoD["seq_tup_list"] = alignmentMap[partId][2]
                infoD["seq_sim"] = alignmentMap[partId][1][1]
                infoD["sort_metric"] = (infoD["seq_sim"] * 100.0 + 1) * 4 + 4
                #
                infoD["REF_ENTRY_ID"] = entityInfo[0]
                infoD["REF_ENTRY_ENTITY_ID"] = entityInfo[1]
                infoD["REF_PDB_ID"] = entityInfo[2]
                infoD["REF_ENTRY_STATUS"] = entityInfo[4]
                infoD["REF_ENTRY_ANN"] = entityInfo[3]
            #
            if auto_match_status:
                seqAnntationInfoMap["auto_match_status"] = True
            #
            return seqAnntationInfoMap
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
            #
        #
        return {}

    def __getEntityDetailsMap(self, cifObj, entityId, authPartsTaxIdInfoList):
        """ """
        itemList = (
            ("seq_part_id", "SEQ_PART_ID"),
            ("seq_part_beg", "SEQ_NUM_BEG"),
            ("seq_part_end", "SEQ_NUM_END"),
            ("seq_part_type", "SEQ_PART_TYPE"),
            ("entity_description", "name"),
            ("entity_synonyms", "synonyms"),
            ("gene_name", "gene"),
            ("taxonomy_id", "taxonomy_id"),
            ("source_scientific_name", "source_scientific"),
            ("source_strain", "strain"),
            ("source_common_name", "source_common"),
            ("entity_enzyme_class", "ec"),
            ("entity_fragment_details", "entity_fragment_details"),
        )  # variant
        #
        partItemList = ("SEQ_PART_ID", "SEQ_NUM_BEG", "SEQ_NUM_END", "SEQ_PART_TYPE", "taxonomy_id")
        #
        annInfoMap = {}
        partsTaxIdInfo = []
        #
        retList = cifObj.GetValue("pdbx_seqtool_entity_details")
        for retDic in retList:
            if ("entity_id" not in retDic) or (retDic["entity_id"] != entityId):
                continue
            #
            infoDic = {}
            for itemTup in itemList:
                if (itemTup[0] in retDic) and retDic[itemTup[0]]:
                    infoDic[itemTup[1]] = retDic[itemTup[0]]
                else:
                    infoDic[itemTup[1]] = ""
                #
            #
            if ("SEQ_PART_ID" in infoDic) and infoDic["SEQ_PART_ID"]:
                annInfoMap[infoDic["SEQ_PART_ID"]] = infoDic
            #
            partTaxIdTup = []
            for partItem in partItemList:
                if (partItem not in infoDic) or (not infoDic[partItem]):
                    continue
                #
                partTaxIdTup.append(infoDic[partItem])
            #
            if len(partTaxIdTup) == len(partItemList):
                partsTaxIdInfo.append(partTaxIdTup)
            #
        #
        if (len(annInfoMap) != len(partsTaxIdInfo)) or (len(partsTaxIdInfo) != len(authPartsTaxIdInfoList)):
            return {}, []
        #
        for i in range(0, len(partsTaxIdInfo)):
            if (
                (str(partsTaxIdInfo[i][1]) != str(authPartsTaxIdInfoList[i][1]))
                or (str(partsTaxIdInfo[i][2]) != str(authPartsTaxIdInfoList[i][2]))
                or (str(partsTaxIdInfo[i][4]) != str(authPartsTaxIdInfoList[i][4]))
            ):
                return {}, []
            #
        #
        return annInfoMap, partsTaxIdInfo

    def __getAlignmentMap(self, cifObj, entityId, partsTaxIdInfo):
        """ """
        alignList = []
        threeLetterSeqList = []
        numMap = {}
        authMap = {}
        refMap = {}
        seq_count = 0
        #
        # Read pdbx_seqtool_mapping_ref category
        #
        retList = cifObj.GetValue("pdbx_seqtool_mapping_ref")
        for retDic in retList:
            if ("entity_id" not in retDic) or (retDic["entity_id"] != entityId):
                continue
            #
            # alignTup[0]: entity_mon_id - three letter code
            # alignTup[1]: entity_seq_num
            # alignTup[2]: ref_mon_id - three letter code
            # alignTup[3]: ref_mon_num
            # alignTup[4]: entity_part_id
            # alignTup[5]: comment
            # alignTup[6]: entity_mon_id - one letter code
            # alignTup[7]: ref_mon_id - one letter code
            #
            alignTup = []
            for item in ("entity_mon_id", "entity_seq_num", "ref_mon_id", "ref_mon_num", "entity_part_id"):
                if (item in retDic) and retDic[item]:
                    alignTup.append(retDic[item])
                else:
                    alignTup.append(self.__gapSymbol)
                #
            #
            alignTup.append("")
            alignTup.append(self.__srd.cnv3To1(alignTup[0]))
            alignTup.append(self.__srd.cnv3To1(alignTup[2]))
            #
            if alignTup[1] == self.__gapSymbol:
                alignTup[1] = ""
            #
            if alignTup[3] == self.__gapSymbol:
                alignTup[3] = ""
            #
            if alignTup[1]:
                numMap[alignTup[1]] = len(alignList)
            #
            if (alignTup[0] != self.__gapSymbol) and (alignTup[1] != self.__gapSymbol):
                if self.__threeLetterCodeSeqList[seq_count] != alignTup[0]:
                    return {}
                #
                seq_count += 1
                authMap[alignTup[0] + "_" + alignTup[1]] = len(alignList)
            #
            if (alignTup[2] != self.__gapSymbol) and alignTup[3] and (alignTup[4] != self.__gapSymbol):
                refMap[alignTup[2] + "_" + alignTup[3] + "_" + alignTup[4]] = len(alignList)
            #
            alignList.append(alignTup)
            threeLetterSeqList.append((alignTup[0], alignTup[6]))
        #
        if (not alignList) or (seq_count != len(self.__threeLetterCodeSeqList)):
            return {}
        #
        for partTup in partsTaxIdInfo:
            if (not partTup[1] in numMap) or (not partTup[2] in numMap):
                return {}
            #
            for i in range(numMap[partTup[1]], numMap[partTup[2]] + 1):
                alignList[i][4] = partTup[0]
            #
        #
        for alignTup in alignList:
            if alignTup[4] == self.__gapSymbol:
                return {}
            #
        #
        # Read pdbx_seqtool_mapping_comment category
        #
        retList = cifObj.GetValue("pdbx_seqtool_mapping_comment")
        for retDic in retList:
            if ("entity_id" not in retDic) or (retDic["entity_id"] != entityId):
                continue
            #
            if (
                ("entity_mon_id" not in retDic)
                or (not retDic["entity_mon_id"])
                or ("entity_seq_num" not in retDic)
                or (not retDic["entity_seq_num"])
                or ("comment" not in retDic)
                or (not retDic["comment"])
            ):
                continue
            #
            if (retDic["entity_mon_id"] + "_" + retDic["entity_seq_num"]) in authMap:
                alignList[authMap[retDic["entity_mon_id"] + "_" + retDic["entity_seq_num"]]][5] = "SELECTED:" + str(retDic["comment"])
            #
        #
        # Read pdbx_seqtool_ref_deletions category
        #
        retList = cifObj.GetValue("pdbx_seqtool_ref_deletions")
        for retDic in retList:
            if ("entity_id" not in retDic) or (retDic["entity_id"] != entityId):
                continue
            #
            if (
                ("ref_mon_id" not in retDic)
                or (not retDic["ref_mon_id"])
                or ("ref_mon_num" not in retDic)
                or (not retDic["ref_mon_num"])
                or ("part_id" not in retDic)
                or (not retDic["part_id"])
            ):
                continue
            #
            if (retDic["ref_mon_id"] + "_" + retDic["ref_mon_num"] + "_" + retDic["part_id"]) in refMap:
                if not alignList[refMap[retDic["ref_mon_id"] + "_" + retDic["ref_mon_num"] + "_" + retDic["part_id"]]][5]:
                    alignList[refMap[retDic["ref_mon_id"] + "_" + retDic["ref_mon_num"] + "_" + retDic["part_id"]]][5] = "deletion"
                #
            #
        #
        return self.__getAuthRefAlignmentMap(alignList, partsTaxIdInfo)

    def __getDbRefMap(self, cifObj, entityId, authRefList=None):
        """
        """
        itemList = (
            ("db_name", "db_name"),
            ("db_code", "db_code"),
            ("db_accession", "db_accession"),
            ("db_isoform", "db_isoform"),
            ("match_begin", "hitFrom"),
            ("match_end", "hitTo"),
        )
        #
        __authRefList = []
        if authRefList:
            __authRefList = authRefList[0]
        #
        unpIdList = []
        gbIdList = []
        dbRefMap = {}
        retList = cifObj.GetValue("pdbx_seqtool_db_ref")
        for retDic in retList:
            if ("entity_id" not in retDic) or (retDic["entity_id"] != entityId) or ("entity_part_id" not in retDic) or (not retDic["entity_part_id"]):
                continue
            #
            infoDic = {}
            for itemTup in itemList:
                if (itemTup[0] in retDic) and retDic[itemTup[0]]:
                    infoDic[itemTup[1]] = retDic[itemTup[0]]
                else:
                    infoDic[itemTup[1]] = ""
                #
            #
            auto_match_status = False
            if __authRefList:
                if (("db_accession" in infoDic) and (infoDic["db_accession"] in __authRefList)) or \
                   (("db_isoform" in infoDic) and (infoDic["db_isoform"] in __authRefList)):
                    auto_match_status = True
                    infoDic["author_provided_id"] = True
                #
            else:
                auto_match_status = True
            #
            if auto_match_status:
                infoDic["auto_match_status"] = True
            #

            if "db_name" in infoDic:
                if infoDic["db_name"] == "UNP":
                    try:
                        start = int(str(infoDic["hitFrom"]))
                    except:  # noqa: E722 pylint: disable=bare-except
                        start = 0
                    #
                    try:
                        end = int(str(infoDic["hitTo"]))
                    except:  # noqa: E722 pylint: disable=bare-except
                        end = 0
                    #
                    if ("db_isoform" in infoDic) and infoDic["db_isoform"]:
                        unpIdList.append((str(infoDic["db_isoform"]), start, end))
                    elif ("db_accession" in infoDic) and infoDic["db_accession"]:
                        unpIdList.append((str(infoDic["db_accession"]), start, end))
                    #
                elif infoDic["db_name"] == "GB":
                    if ("db_accession" in infoDic) and infoDic["db_accession"]:
                        gbIdList.append(infoDic["db_accession"])
                    #
                #
            #
            infoDic["id"] = 101
            infoDic["sort_order"] = 101
            #
            dbRefMap[retDic["entity_part_id"]] = infoDic
        #
        fetchSeqUtil = FetchSeqInfoUtils(siteId=self.__siteId, seqReferenceData=self.__srd, verbose=self.__verbose, log=self.__lfh)
        unpD = {}
        gbD = {}
        if len(unpIdList) > 0:
            unpD = fetchSeqUtil.fetchUniProt(idTupleList=unpIdList)
        elif len(gbIdList) > 0:
            gbIdList = list(set(gbIdList))
            for idCode in gbIdList:
                # No need to rate limit here as this fetches the whole entry - is likely a genome and takes seconds.
                gbD[idCode] = fetchSeqUtil.fetchNcbiGi(idCode)
            #
        #
        for _partId, infoDic in dbRefMap.items():
            dbDic = {}
            if "db_name" in infoDic:
                if infoDic["db_name"] == "UNP":
                    try:
                        start = int(str(infoDic["hitFrom"]))
                    except:  # noqa: E722 pylint: disable=bare-except
                        start = 0
                    #
                    try:
                        end = int(str(infoDic["hitTo"]))
                    except:  # noqa: E722 pylint: disable=bare-except
                        end = 0
                    #
                    if ("db_isoform" in infoDic) and infoDic["db_isoform"]:
                        if (infoDic["db_isoform"], start, end) in unpD:
                            dbDic = unpD[(infoDic["db_isoform"], start, end)]
                        #
                    elif ("db_accession" in infoDic) and infoDic["db_accession"]:
                        if (infoDic["db_accession"], start, end) in unpD:
                            dbDic = unpD[(infoDic["db_accession"], start, end)]
                        #
                    #
                elif infoDic["db_name"] == "GB":
                    if ("db_accession" in infoDic) and infoDic["db_accession"]:
                        if infoDic["db_accession"] in gbD:
                            dbDic = gbD[infoDic["db_accession"]]
                        #
                    #
                #
            #
            if not dbDic:
                return {}
            #
            # Remove merging "ec" number
            # for item in ( 'db_isoform_description', 'db_description', 'ec', 'db_description', 'comments', 'keyword' ):
            for item in ("db_isoform_description", "db_description", "comments", "keyword"):
                if item in dbDic:
                    infoDic[item] = dbDic[item]
                else:
                    infoDic[item] = ""
                #
            #
            if "sequence" in dbDic:
                infoDic["db_length"] = len(dbDic["sequence"])
            else:
                infoDic["db_length"] = 0
            #
        #
        return dbRefMap

    def __getSelfRefmap(self, cifObj, entityId, includeSelfRef):
        """ """
        selfRefmap = {}
        #
        if not includeSelfRef:
            return selfRefmap
        #
        retList = cifObj.GetValue("pdbx_seqtool_self_ref")
        for retDic in retList:
            if ("entity_id" not in retDic) or (retDic["entity_id"] != entityId) or ("entity_part_id" not in retDic) or (not retDic["entity_part_id"]):
                continue
            #
            selfRefmap[retDic["entity_part_id"]] = True
        #
        if selfRefmap and (not self.__checkSelfRefSeq(cifObj, entityId)):
            return {}
        #
        return selfRefmap

    def __getAuthRefAlignmentMap(self, matchList, partsTaxIdInfo):
        """ """
        alignmentMap = {}
        for partTup in partsTaxIdInfo:
            authIdx = -1
            refIdx = -1
            alignLength = 0
            seqLength = 0
            numMatch = 0
            numMatchGaps = 0
            alignmentList = []
            sTup3L = []
            for matchTup in matchList:
                currAuthIdx = -1
                if matchTup[0] != self.__gapSymbol:
                    authIdx += 1
                    currAuthIdx = authIdx
                    seqLength += 1
                #
                alignTup = []
                if matchTup[4] == partTup[0]:
                    alignLength += 1
                    if matchTup[0] == matchTup[2]:
                        numMatch += 1
                        numMatchGaps += 1
                    elif matchTup[2] == self.__gapSymbol:
                        numMatchGaps += 1
                    #
                    currRefIdx = -1
                    if matchTup[2] != self.__gapSymbol:
                        refIdx += 1
                        currRefIdx = refIdx
                        sTup3L.append((matchTup[2], matchTup[3], matchTup[5], currRefIdx + 1, matchTup[7]))
                    #
                    alignTup = [codeIndex(currAuthIdx, matchTup[5]), codeIndex(currRefIdx, matchTup[5])]
                else:
                    alignTup = [currAuthIdx, -1]
                #
                alignmentList.append(alignTup)
            #
            alignmentMap[partTup[0]] = (alignmentList, (alignLength, float(numMatch) / float(alignLength), float(numMatchGaps) / float(seqLength)), sTup3L)
        #
        return alignmentMap

    def __checkSelfRefSeq(self, cifObj, entityId):
        """ """
        threeLetterSeqList = []
        auth_asym_id = ""
        #
        # Read pdbx_seqtool_mapping_xyz category
        #
        retList = cifObj.GetValue("pdbx_seqtool_mapping_xyz")
        for retDic in retList:
            if ("entity_id" not in retDic) or (retDic["entity_id"] != entityId) or ("entity_mon_id" not in retDic) or ("auth_asym_id" not in retDic):
                continue
            #
            if not auth_asym_id:
                auth_asym_id = retDic["auth_asym_id"]
            #
            if auth_asym_id != retDic["auth_asym_id"]:
                break
            #
            threeLetterSeqList.append((retDic["entity_mon_id"], self.__srd.cnv3To1(retDic["entity_mon_id"])))
        #
        if (not threeLetterSeqList) or (len(threeLetterSeqList) != len(self.__threeLetterCodeSeqList)):
            return False
        #
        for i in range(0, len(threeLetterSeqList)):
            if str(threeLetterSeqList[i][0]) != str(self.__threeLetterCodeSeqList[i]):
                if str(threeLetterSeqList[i][1]) != str(self.__oneLetterCodeSeqList[i]):
                    return False
                #
            #
        #
        return True


def __mymain():
    # from wwpdb.utils.config.ConfigInfo import ConfigInfo
    # from wwpdb.apps.seqmodule.webapp.SeqModWebRequest import SeqModInputRequest

    #
    # siteId = os.getenv("WWPDB_SITE_ID")
    # cI = ConfigInfo(siteId)
    #
    # myReqObj = SeqModInputRequest({}, verbose=True, log=sys.stderr)
    # myReqObj.setValue("TopSessionPath", cI.get("SITE_WEB_APPS_TOP_SESSIONS_PATH"))
    # myReqObj.setValue("TopPath", cI.get("SITE_WEB_APPS_TOP_PATH"))
    # myReqObj.setValue("WWPDB_SITE_ID", siteId)
    # myReqObj.setValue("sessionid", "447a9d00e63720c7df7ca0bf755f004dec2acb44")
    annObj = GetSameSeqAnnotation(verbose=True, log=sys.stderr)
    annObj.testGetSeqAnnotationFromAssignFile(sys.argv[1], sys.argv[2], [sys.argv[3]])


if __name__ == "__main__":
    __mymain()
