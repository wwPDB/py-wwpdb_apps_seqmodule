##
# File:    ReferenceSequenceUtilsTests.py
# Date:    24-Feb-2013
#
# Updates:
# 24-Feb-2013  jdw added tests for fetch methods for UniProt and  NCBI sequence/taxonomy entries.
#
##
"""
Test cases for performing reference sequence database searches and processing search results.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.07"


import sys
import unittest
import traceback
import time
import os
import os.path

from wwpdb.apps.seqmodule.io.ReferenceSequenceUtils import ReferenceSequenceUtils
from wwpdb.apps.seqmodule.io.ModelSequenceUtils import ModelSequenceUtils
from wwpdb.apps.seqmodule.io.PdbxIoUtils import PdbxFileIo, ReferenceSequenceIo
from wwpdb.io.misc.FormatOut import FormatOut
from wwpdb.io.locator.PathInfo import PathInfo
from wwpdb.utils.config.ConfigInfo import ConfigInfo, getSiteId


class ReferenceSequenceUtilsTests(unittest.TestCase):

    def setUp(self):
        #
        self.__verbose = True
        self.__debug = True
        self.__lfh = sys.stdout
        self.__siteId = getSiteId(defaultSiteId="WWPDB_DEPLOY_TEST")
        self.__pathExamplesRel = "../tests"
        self.__pathExamples = os.path.abspath(self.__pathExamplesRel)
        #
        self.__examFileList = ['1cbs.cif', '3rer.cif', 'rcsb056751.cif']
        #
        self.__unpIdList = ['P20937', 'P21877', 'P22868', 'P23832', 'P25665', 'P26562', 'P27614', 'E3PJ86', 'P42284-1', 'Q9H8M2-1']
        self.__giIdList = ['28948918', '9845486', '116077928']
        self.__taxIdList = ['9606', '270523']

    def testSearchEntitySequences(self):
        """ Test reference sequence database search for each entity sequence extracted from the input
            model coordinate file.  The appropriate reference sequence database service is selected.
            Matching sequences are extracted from service output files and stored along with key
            reference sequence features in PDBx data files.
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            pI = PathInfo(siteId=self.__siteId, sessionPath='.', verbose=self.__verbose, log=self.__lfh)
            for f in self.__examFileList:
                fN = os.path.join(self.__pathExamples, f)
                c0 = PdbxFileIo(verbose=self.__verbose, log=self.__lfh).getContainer(fN)
                msu = ModelSequenceUtils(dataContainer=c0, verbose=self.__verbose, log=self.__lfh)
                entityD = msu.getEntitySequenceDetails()
                if self.__debug:
                    fOut = FormatOut()
                    fOut.autoFormat("Entity dictionary", entityD, 3, 3)
                    fOut.writeStream(self.__lfh)
                #
                rsu = ReferenceSequenceUtils(verbose=self.__verbose, log=self.__lfh)
                rsio = ReferenceSequenceIo(verbose=self.__verbose, log=self.__lfh)
                for eId, eD in entityD.items():
                    entryIdL = str(eD['ENTRY_ID']).lower()
                    if len(eD['SEQ_ENTITY_1']) > 50:
                        mR = rsu.searchEntities(entityD=eD, saveBlast=True, filePrefix=entryIdL)
                        fn = pI.getReferenceSequenceFilePath(entryIdL, entityId=eId, fileSource='session')
                        rsio.writeMatchResults(eD, outFilePath=fn, matchResults=mR)
                        rc0 = PdbxFileIo(verbose=self.__verbose, log=self.__lfh).getContainer(fn)
                        #
                        rsin = ReferenceSequenceIo(dataContainer=rc0, verbose=self.__verbose, log=self.__lfh)
                        rL = rsin.readMatchResults()
                        self.__lfh.write("Entity %s Results %r\n" % (eId, rL))

        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testFetchUniProtEntry(self):
        """ Test fetch UniProt entry - -
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            rsu = ReferenceSequenceUtils(verbose=self.__verbose, log=self.__lfh)
            d = rsu.fetchUniProt(self.__unpIdList, filePath="Test-unp.xml")
            for k, v in d.items():
                self.__lfh.write("\nKey %s value %sr\n" % (k, v))
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testFetchGiEntry(self):
        """ Test fetch UniProt entry - -
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            rsu = ReferenceSequenceUtils(verbose=self.__verbose, log=self.__lfh)
            for gi in self.__giIdList:
                d = rsu.fetchNcbiGi(gi)
                for k, v in d.items():
                    self.__lfh.write("gi %s Key %s value %r\n" % (gi, k, v))
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testFetchTaxIdEntry(self):
        """ Test fetch UniProt entry - -
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            rsu = ReferenceSequenceUtils(verbose=self.__verbose, log=self.__lfh)
            for taxId in self.__taxIdList:
                d = rsu.fetchNcbiTaxId(taxId)
                for k, v in d.items():
                    self.__lfh.write("Key %s value %s\n" % (k, v))
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()


def suiteSeqSearchTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ReferenceSequenceUtilsTests("testSearchEntitySequences"))
    return suiteSelect


def suiteSeqFetchTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ReferenceSequenceUtilsTests("testFetchUniProtEntry"))
    # suiteSelect.addTest(ReferenceSequenceUtilsTests("testFetchTaxIdEntry"))
    # suiteSelect.addTest(ReferenceSequenceUtilsTests("testFetchGiEntry"))
    return suiteSelect


if __name__ == '__main__':
    if (False):
        mySuite = suiteSeqSearchTests()
        unittest.TextTestRunner(verbosity=2).run(mySuite)
    if (True):
        mySuite = suiteSeqFetchTests()
        unittest.TextTestRunner(verbosity=2).run(mySuite)
