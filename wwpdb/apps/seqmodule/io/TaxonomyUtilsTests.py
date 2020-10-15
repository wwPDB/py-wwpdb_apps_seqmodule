##
# File:    TaxonomyUtilsTests.py
# Date:    21-Feb-2013
#
# Updates:
#  24-Feb-2013 jdw  remove the dependency on the reqObj.
##
"""
Test cases for taxonomy reference data management.

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
import time

from wwpdb.apps.seqmodule.io.TaxonomyUtils import TaxonomyUtils
from wwpdb.io.misc.FormatOut import FormatOut
from wwpdb.utils.config.ConfigInfo import ConfigInfo, getSiteId


class TaxonomyUtilsTests(unittest.TestCase):

    def setUp(self):
        self.__verbose = True
        self.__lfh = sys.stdout
        self.__siteId = getSiteId(defaultSiteId="WWPDB_DEPLOY_TEST")

    def tearDown(self):
        pass

    def testReadTaxonomyNames(self):
        startTime = time.time()
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__,
                                                       sys._getframe().f_code.co_name,
                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime())))
        try:
            tu = TaxonomyUtils(siteId=self.__siteId, verbose=self.__verbose, log=self.__lfh)
            nameD = tu.readTaxonomyNameData(serialize=False)
            fOut = FormatOut()
            fOut.autoFormat("Taxonomy names", nameD, 3, 3)
            fOut.writeStream(self.__lfh)
            fOut.clear()

        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write("\nCompleted %s %s at %s (%.2f seconds)\n" % (self.__class__.__name__,
                                                                       sys._getframe().f_code.co_name,
                                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                                       endTime - startTime))

    def testReadAndSerializeTaxonomyNames(self):
        startTime = time.time()
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__,
                                                       sys._getframe().f_code.co_name,
                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime())))
        try:
            tu = TaxonomyUtils(siteId=self.__siteId, verbose=self.__verbose, log=self.__lfh)
            nameD = tu.readTaxonomyNameData(serialize=True)
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write("\nCompleted %s %s at %s (%.2f seconds)\n" % (self.__class__.__name__,
                                                                       sys._getframe().f_code.co_name,
                                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                                       endTime - startTime))

    def testDeserializeTaxonomyNames(self):
        startTime = time.time()
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__,
                                                       sys._getframe().f_code.co_name,
                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime())))
        try:
            tu = TaxonomyUtils(siteId=self.__siteId, verbose=self.__verbose, log=self.__lfh)
            d = tu.getNames()
            self.__lfh.write("testDeserializeTaxonomyNodes() name dictionary length %d\n" % len(d))
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()
        endTime = time.time()
        self.__lfh.write("\nCompleted %s %s at %s (%.2f seconds)\n" % (self.__class__.__name__,
                                                                       sys._getframe().f_code.co_name,
                                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                                       endTime - startTime))

    def testLookupNames(self):
        startTime = time.time()
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__,
                                                       sys._getframe().f_code.co_name,
                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime())))
        try:
            tu = TaxonomyUtils(siteId=self.__siteId, verbose=self.__verbose, log=self.__lfh)
            for taxId in ['83333', '9606', 9606]:
                self.__lfh.write("testLookupNames() taxid %r names %r\n" % (taxId, tu.lookUpSource(taxId)))

        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()
        endTime = time.time()
        self.__lfh.write("\nCompleted %s %s at %s (%.2f seconds)\n" % (self.__class__.__name__,
                                                                       sys._getframe().f_code.co_name,
                                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                                       endTime - startTime))

    def testLookupParents(self):
        startTime = time.time()
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__,
                                                       sys._getframe().f_code.co_name,
                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime())))
        try:
            tu = TaxonomyUtils(siteId=self.__siteId, verbose=self.__verbose, log=self.__lfh)
            for taxId in ['562', '83333', '316407', '679895', '1318715', '595496', '531853', '511145', '527799', '1211845', '1110693', '1403831']:
                anD = tu.getAncestors(taxId)
                self.__lfh.write("testLookupNames() taxid %r names %r parent dict %r\n" % (taxId, tu.lookUpSource(taxId), anD))

        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()
        endTime = time.time()
        self.__lfh.write("\nCompleted %s %s at %s (%.2f seconds)\n" % (self.__class__.__name__,
                                                                       sys._getframe().f_code.co_name,
                                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                                       endTime - startTime))

    def testReadTaxonomyNodes(self):
        startTime = time.time()
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__,
                                                       sys._getframe().f_code.co_name,
                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime())))
        try:
            tu = TaxonomyUtils(siteId=self.__siteId, verbose=self.__verbose, log=self.__lfh)
            nodeD = tu.readTaxonomyNodeData(serialize=False)
            fOut = FormatOut()
            fOut.autoFormat("Taxonomy nodes", nodeD, 3, 3)
            fOut.writeStream(self.__lfh)
            fOut.clear()

        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write("\nCompleted %s %s at %s (%.2f seconds)\n" % (self.__class__.__name__,
                                                                       sys._getframe().f_code.co_name,
                                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                                       endTime - startTime))

    def testReadAndSerializeTaxonomyNodes(self):
        startTime = time.time()
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__,
                                                       sys._getframe().f_code.co_name,
                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime())))
        try:
            tu = TaxonomyUtils(siteId=self.__siteId, verbose=self.__verbose, log=self.__lfh)
            nodeD = tu.readTaxonomyNodeData(serialize=True)
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write("\nCompleted %s %s at %s (%.2f seconds)\n" % (self.__class__.__name__,
                                                                       sys._getframe().f_code.co_name,
                                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                                       endTime - startTime))

    def testDeserializeTaxonomyNodes(self):
        startTime = time.time()
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__,
                                                       sys._getframe().f_code.co_name,
                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime())))
        try:
            tu = TaxonomyUtils(siteId=self.__siteId, verbose=self.__verbose, log=self.__lfh)
            d = tu.getParentDict()
            self.__lfh.write("testDeserializeTaxonomyNodes() node dictionary length %d\n" % len(d))
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write("\nCompleted %s %s at %s (%.2f seconds)\n" % (self.__class__.__name__,
                                                                       sys._getframe().f_code.co_name,
                                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                                       endTime - startTime))


def suiteTaxonomyNodeTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(TaxonomyUtilsTests("testReadTaxonomyNodes"))
    suiteSelect.addTest(TaxonomyUtilsTests("testReadAndSerializeTaxonomyNodes"))
    suiteSelect.addTest(TaxonomyUtilsTests("testDeserializeTaxonomyNodes"))
    return suiteSelect


def suiteTaxonomyNameTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(TaxonomyUtilsTests("testReadTaxonomyNames"))
    suiteSelect.addTest(TaxonomyUtilsTests("testReadAndSerializeTaxonomyNames"))
    suiteSelect.addTest(TaxonomyUtilsTests("testDeserializeTaxonomyNames"))
    return suiteSelect


def suiteLookupNameTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(TaxonomyUtilsTests("testLookupNames"))
    suiteSelect.addTest(TaxonomyUtilsTests("testLookupParents"))
    return suiteSelect

if __name__ == '__main__':
    if (True):
        mySuite = suiteTaxonomyNodeTests()
        unittest.TextTestRunner(verbosity=2).run(mySuite)

        mySuite = suiteTaxonomyNameTests()
        unittest.TextTestRunner(verbosity=2).run(mySuite)

    mySuite = suiteLookupNameTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
