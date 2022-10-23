##
# File:  SeqReferenceSearchUtils.py
# Date:  22-Nov-2018
#
# Updates:
#
# 19-Oct-2022 zf add author provided reference sequence information
#
##
"""
Wrap class for blast sequence search and 100% matched previous processed sequence search
"""
__docformat__ = "restructuredtext en"
__author__ = "Zukang Feng"
__email__ = "zfeng@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.09"

import sys
import traceback

from wwpdb.apps.seqmodule.util.LocalBlastSearchUtils import LocalBlastSearchUtils
from wwpdb.apps.seqmodule.util.GetSameSeqAnnotationWithTaxIdOnly import GetSameSeqAnnotation
from wwpdb.apps.seqmodule.util.SearchEntityPolySeqs import SearchEntityPolySeqs


class SeqReferenceSearchUtils(object):
    """ """

    def __init__(self, siteId=None, sessionPath=None, searchUtil=None, pathInfo=None, verbose=False, log=sys.stderr, ncbilock=None):
        """ """
        self.__siteId = siteId
        self.__sessionPath = sessionPath
        self.__sepsUtil = searchUtil
        self.__pI = pathInfo
        self.__verbose = verbose
        self.__lfh = log
        self.__ncbilock = ncbilock

    def run(self, dataSetId, eD, authDefinedRefList, refRefSearchFlag, forceBlastSearchFlag):
        """ """
        sasUtil = SeqAnnotationSearchUtils(
            siteId=self.__siteId, sessionPath=self.__sessionPath, searchUtil=self.__sepsUtil, pathInfo=self.__pI, verbose=self.__verbose, log=self.__lfh
        )
        selfRefD, sameSeqRefD = sasUtil.getSameSeqRefInfo(dataSetId, eD, authDefinedRefList)
        #
        eRefD = {}
        if (len(selfRefD) == 0) or forceBlastSearchFlag:
            lbsUtil = LocalBlastSearchUtils(
                siteId=self.__siteId,
                sessionPath=self.__sessionPath,
                pathInfo=self.__pI,
                doRefSearchFlag=refRefSearchFlag,
                verbose=self.__verbose,
                log=self.__lfh,
                ncbilock=self.__ncbilock,
            )
            _autoMatchStatus, eRefD = lbsUtil.searchSeqReference(dataSetId=dataSetId, entityD=eD, authRefList=authDefinedRefList)
        #
        return eRefD, selfRefD, sameSeqRefD


class SeqAnnotationSearchUtils(object):
    """ """

    def __init__(self, siteId=None, sessionPath=None, searchUtil=None, pathInfo=None, verbose=False, log=sys.stderr):
        """ """
        self.__siteId = siteId
        self.__sessionPath = sessionPath
        self.__sepsUtil = searchUtil
        self.__pI = pathInfo
        self.__verbose = verbose
        self.__lfh = log
        #
        if not self.__sepsUtil:
            self.__sepsUtil = SearchEntityPolySeqs(siteId=self.__siteId, sessionPath=self.__sessionPath, verbose=self.__verbose, log=self.__lfh)
        #

    def getSameSeqRefInfo(self, dataSetId, eD, authRefList):
        """ """
        if not eD:
            return {}, {}
        #
        try:
            selfList, sameSeqEntryList = self.__sepsUtil.searchSameSequenceEntity(entryId=dataSetId, seq=eD["SEQ_ENTITY_1"])
            selfInfoMap = {}
            sameSeqInfoMap = {}
            if (len(selfList) > 0) or (len(sameSeqEntryList) > 0):
                partsTaxIdInfoList = []
                for (partNo, pD) in enumerate(eD["PART_LIST"], start=1):
                    if ("SOURCE_TAXID" in pD) and pD["SOURCE_TAXID"] and ("SEQ_NUM_BEG" in pD) and pD["SEQ_NUM_BEG"] and \
                       ("SEQ_NUM_END" in pD) and pD["SEQ_NUM_END"]:
                        partsTaxIdInfoList.append([str(partNo), pD["SEQ_NUM_BEG"], pD["SEQ_NUM_END"], "", pD["SOURCE_TAXID"]])
                    elif (len(eD["PART_LIST"]) == 1) and ("SOURCE_TAXID" in pD) and pD["SOURCE_TAXID"]:
                        partsTaxIdInfoList.append([str(partNo), "", "", "", pD["SOURCE_TAXID"]])
                    #
                #
                annObj = GetSameSeqAnnotation(siteId=self.__siteId, sessionPath=self.__sessionPath, pathInfo=self.__pI, verbose=self.__verbose, log=self.__lfh)
                annObj.setEntitySeq(seq=eD["SEQ_ENTITY_1"], polyTypeCode=eD["POLYMER_TYPE_CODE"])
                if len(selfList) > 0:
                    selfInfoMap = annObj.getSeqAnnotationFromAssignFile(retList=selfList, authPartsTaxIdInfoList=partsTaxIdInfoList,
                                                                        authRefList=authRefList, includeSelfRef=True)
                #
                if len(sameSeqEntryList) > 0:
                    sameSeqInfoMap = annObj.getSeqAnnotationFromAssignFile(retList=sameSeqEntryList, authPartsTaxIdInfoList=partsTaxIdInfoList,
                                                                           authRefList=authRefList)
                #
            #
            return selfInfoMap, sameSeqInfoMap
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            return {}, {}
        #
