##
# File:  SeqReferenceSearchUtils.py
# Date:  22-Nov-2018
#
# Updates:
##
"""
Wrap class for blast sequence search and 100% matched previous processed sequence search 
"""
__docformat__ = "restructuredtext en"
__author__ = "Zukang Feng"
__email__ = "zfeng@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.09"

import os, sys, traceback

from wwpdb.apps.seqmodule.util.LocalBlastSearchUtils import LocalBlastSearchUtils
from wwpdb.apps.seqmodule.util.GetSameSeqAnnotationWithTaxIdOnly import GetSameSeqAnnotation
from wwpdb.apps.seqmodule.util.SearchEntityPolySeqs import SearchEntityPolySeqs

class SeqReferenceSearchUtils(object):
    """
    """
    def __init__(self, siteId=None, sessionPath=None, searchUtil=None, pathInfo=None, verbose=False, log=sys.stderr, ncbilock=None):
        """
        """
        self.__siteId = siteId
        self.__sessionPath = sessionPath
        self.__sepsUtil = searchUtil
        self.__pI = pathInfo
        self.__verbose = verbose
        self.__lfh = log
        self.__ncbilock = ncbilock

    def run(self, dataSetId, eD, refRefSearchFlag, forceBlastSearchFlag):
        """
        """
        sasUtil = SeqAnnotationSearchUtils(siteId=self.__siteId, sessionPath=self.__sessionPath, searchUtil=self.__sepsUtil, pathInfo=self.__pI, \
                                           verbose=self.__verbose, log=self.__lfh)
        selfRefD,sameSeqRefD = sasUtil.getSameSeqRefInfo(dataSetId, eD)
        #
        eRefD = {}
        if (len(selfRefD) == 0) or forceBlastSearchFlag:
            lbsUtil = LocalBlastSearchUtils(siteId=self.__siteId, 
                                            sessionPath=self.__sessionPath, 
                                            pathInfo=self.__pI, 
                                            doRefSearchFlag=refRefSearchFlag, 
                                            verbose=self.__verbose, log=self.__lfh,
                                            ncbilock=self.__ncbilock)
            eRefD = lbsUtil.searchSeqReference(dataSetId=dataSetId, entityD=eD)
        #
        return eRefD,selfRefD,sameSeqRefD

class SeqAnnotationSearchUtils(object):
    """
    """
    def __init__(self, siteId=None, sessionPath=None, searchUtil=None, pathInfo=None, verbose=False, log=sys.stderr):
        """
        """
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

    def getSameSeqRefInfo(self, dataSetId, eD):
        """
        """
        if not eD:
            return {},{}
        #
        try:
            selfList,sameSeqEntryList = self.__sepsUtil.searchSameSequenceEntity(entryId=dataSetId, seq=eD["SEQ_ENTITY_1"])
            selfInfoMap = {}
            sameSeqInfoMap = {}
            if (len(selfList) > 0) or (len(sameSeqEntryList) > 0):
                taxIdList = []
                for (partNo, pD) in enumerate(eD["PART_LIST"], start=1):
                    if ("SOURCE_TAXID" in pD) and pD["SOURCE_TAXID"] and (pD["SOURCE_TAXID"] not in taxIdList):
                        taxIdList.append(pD["SOURCE_TAXID"])
                    #
                #
                annObj = GetSameSeqAnnotation(siteId=self.__siteId, sessionPath=self.__sessionPath, pathInfo=self.__pI, verbose=self.__verbose, log=self.__lfh)
                annObj.setEntitySeq(seq=eD["SEQ_ENTITY_1"], polyTypeCode=eD["POLYMER_TYPE_CODE"])
                if len(selfList) > 0:
                    selfInfoMap = annObj.getSeqAnnotationFromAssignFile(retList=selfList, TaxIdList=taxIdList)
                #
                if len(sameSeqEntryList) > 0:
                    sameSeqInfoMap = annObj.getSeqAnnotationFromAssignFile(retList=sameSeqEntryList, TaxIdList=taxIdList)
                #
            #
            return selfInfoMap,sameSeqInfoMap
        except:
            traceback.print_exc(file=self.__lfh)
            return {},{}
        #
