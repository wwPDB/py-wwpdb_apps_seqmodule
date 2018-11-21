##
# File:  ChemCompUtils.py
# Date:  23-Apr-2013
#
# Updates:
#   24-Apr-2013 jdw Add method to create all required object and index files -
##
"""
Accessors for chemical component index of residue modifications.

Python style serialization and deserialization of chemical component index structures is provided
and coordinated within the project reference data storage model.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.07"


import os
import sys
import traceback
import time
import cPickle
from wwpdb.utils.config.ConfigInfo import ConfigInfo
from wwpdb.utils.cc_dict_util.persist.PdbxChemCompDictUtil import PdbxChemCompDictUtil
from wwpdb.utils.cc_dict_util.persist.PdbxChemCompDictIndex import PdbxChemCompDictIndex

class ChemCompUtils(object):
    """ Accessors for chemical component index of residue modifications.

       Python style serialization and deserialization of chemical component index structures is provided
       and coordinated within the project reference data storage model.
    """

    def __init__(self, siteId='WWPDB_DEPLOY_TEST', verbose=True, log=sys.stderr):
        self.__siteId = siteId
        self.__verbose = verbose
        self.__debug = False
        self.__lfh = log
        #
        self.__cI = ConfigInfo(self.__siteId)
        self.__pathChemCompParentIndexFile = os.path.join(self.__cI.get('SITE_REFDATA_CHEM_COMP_INDEX_PATH'), 'chem-comp-parent-index.pic')
        self.__pathChemCompIndexFile = os.path.join(self.__cI.get('SITE_REFDATA_CHEM_COMP_INDEX_PATH'), 'chem-comp-index.pic')
        #
        self.__pathChemCompDictStoreFile = os.path.join(self.__cI.get('SITE_REFDATA_CHEM_COMP_INDEX_PATH'), 'chemcomp-all.db')
        self.__pathChemCompDictFile = os.path.join(self.__cI.get('SITE_CC_DICT_PATH'), 'Components-all-v3.cif')
        #
        self.__pickleProtocol = cPickle.HIGHEST_PROTOCOL
        #
        self.__parentD = None
        self.__childD = None
        self.__setup()
        #

    def __setup(self):
        self.__readChemCompParentIndex()

    def __createChemCompStore(self):
        startTime = time.time()
        try:
            dUtil = PdbxChemCompDictUtil(verbose=self.__verbose, log=self.__lfh)
            dUtil.makeStoreFromFile(dictPath=self.__pathChemCompDictFile, storePath=self.__pathChemCompDictStoreFile)
            ok = True
        except:
            ok = False
            traceback.print_exc(file=self.__lfh)

        endTime = time.time()
        if (self.__verbose):
            self.__lfh.write("+ChemCompUtils(__createChemCompStore) completed with status %r at %s (%d seconds)\n" %
                             (ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime))
        return ok

    def __createChemCompIndices(self):
        self.__createChemCompIndex()
        self.__createChemCompParentIndex()

    def __createChemCompIndex(self):
        startTime = time.time()
        try:
            if (not os.access(self.__pathChemCompDictStoreFile, os.R_OK)):
                self.__createChemCompStore()
            dIndx = PdbxChemCompDictIndex(verbose=self.__verbose, log=self.__lfh)
            dIndx.makeIndex(storePath=self.__pathChemCompDictStoreFile, indexPath=self.__pathChemCompIndexFile)
            ok = True
        except:
            ok = False
            traceback.print_exc(file=self.__lfh)

        endTime = time.time()
        if (self.__verbose):
            self.__lfh.write("+ChemCompUtils(__createChemCompIndex) completed with status %r at %s (%d seconds)\n" %
                             (ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime))
        return ok

    def __createChemCompParentIndex(self):
        startTime = time.time()
        try:
            if (not os.access(self.__pathChemCompDictStoreFile, os.R_OK)):
                self.__createChemCompStore()
            dIndx = PdbxChemCompDictIndex(verbose=self.__verbose, log=self.__lfh)
            pD, cD = dIndx.makeParentComponentIndex(storePath=self.__pathChemCompDictStoreFile, indexPath=self.__pathChemCompParentIndexFile)
            ok = True
        except:
            ok = False
            traceback.print_exc(file=self.__lfh)

        endTime = time.time()
        if (self.__verbose):
            self.__lfh.write("+ChemCompUtils(__createChemCompParentIndex) completed with status %r at %s (%d seconds)\n" %
                             (ok, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime))
        return ok

    def __readChemCompParentIndex(self):
        try:
            if (not os.access(self.__pathChemCompParentIndexFile, os.R_OK)):
                self.__createChemCompIndices()
            dIndx = PdbxChemCompDictIndex(verbose=self.__verbose, log=self.__lfh)
            self.__parentD, self.__childD = dIndx.readParentComponentIndex(indexPath=self.__pathChemCompParentIndexFile)
            ok = True
            if (self.__debug):
                self.__lfh.write("+ChemCompUtils(__readChemCompParentIndex) recovered parent dictionary %r\n" % self.__parentD.items())
                self.__lfh.write("+ChemCompUtils(__readChemCompParentIndex) recovered child dictionary %r\n" % self.__childD.items())
        except:
            self.__lfh.write("+ChemCompUtils(__readChemCompParentIndex) reading failed %s\n" % self.__pathChemCompParentIndexFile)
            ok = False
            traceback.print_exc(file=self.__lfh)
        return ok

    def getParentComponentList(self, modCompId):
        try:
            if self.__childD is None:
                self.__readChemCompParentIndex()
            return self.__childD[modCompId]
        except:
            return None

    def getModificationList(self, parentCompId):
        try:
            if self.__parentD is None:
                self.__readChemCompParentIndex()
            return self.__parentD[parentCompId]
        except:
            return []

    def isModificationOf(self, parentCompId, modCompId):
        try:
            if self.__parentD is None:
                self.__readChemCompParentIndex()
            return (modCompId in self.__parentD[parentCompId])
        except:
            if (self.__debug):
                self.__lfh.write("+ChemCompUtils(isModificationOf) Failed for parentCompId %s modCompId %s\n" % (parentCompId, modCompId))
                traceback.print_exc(file=self.__lfh)
        return False
