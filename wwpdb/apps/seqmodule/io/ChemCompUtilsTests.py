##
# File:    ChemCompUtilsTests.py
# Date:    23-Apr-2013
#
# Updates:
##
"""
Test cases for chemical component index construction and access -

Here the parent index is created/exercised.

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

from wwpdb.apps.seqmodule.io.ChemCompUtils import ChemCompUtils
from wwpdb.utils.rcsb.FormatOut import FormatOut
from wwpdb.api.facade.ConfigInfo import ConfigInfo, getSiteId


class ChemCompUtilsTests(unittest.TestCase):

    def setUp(self):
        self.__verbose = True
        self.__lfh = sys.stdout
        self.__siteId = getSiteId(defaultSiteId="WWPDB_DEPLOY_TEST")

    def tearDown(self):
        pass

    def testChemCompParentIndex(self):
        startTime = time.clock()
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__,
                                                       sys._getframe().f_code.co_name,
                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime())))
        try:
            ccU = ChemCompUtils(siteId=self.__siteId, verbose=self.__verbose, log=self.__lfh)
            for compId in ['ALA', 'TRP', 'MET', 'LYS']:
                self.__lfh.write("+testChemCompParentIndex- Modification list for component %s : %r\n" %
                                 (compId, ccU.getModificationList(compId)))

            mL = ccU.getModificationList('MET')
            for compId in mL:
                pCompId = ccU.getParentComponentList(compId)
                self.__lfh.write("+testChemCompParentIndex-  modified component %s parent list %s\n" % (compId, pCompId))
                self.__lfh.write("+testChemCompParentIndex-  parentCompId  %s modified compId  %s isModfication() %r\n" %
                                 (pCompId[0], compId, ccU.isModificationOf(pCompId[0], compId)))
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.clock()
        self.__lfh.write("\nCompleted %s %s at %s (%.2f seconds)\n" %
                         (self.__class__.__name__, sys._getframe().f_code.co_name,
                          time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime))


def suiteParentIndexTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ChemCompUtilsTests("testChemCompParentIndex"))
    return suiteSelect

if __name__ == '__main__':

    mySuite = suiteParentIndexTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
