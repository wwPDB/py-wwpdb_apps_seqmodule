##
# File:    AlignmentStatisticsTests.py
# Date:    26-Feb-2013
#
# Updates:
# 04-Mar-2013  jdw added test case for multipart entity.
# 07-Mar-2013  jdw archive and upload examples
##
"""
Test cases for computing alignment statistics for an existing sequence data store.

The sequence data store is created by the  SequenceDataAssembly() class.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.07"


import sys
import unittest
import shutil
import traceback
import time
import os
import os.path
import inspect

from wwpdb.apps.seqmodule.align.AlignmentStatistics import AlignmentStatistics
from wwpdb.apps.seqmodule.control.SequenceDataAssemble_v2 import SequenceDataAssemble
from wwpdb.io.locator.PathInfo import PathInfo
from wwpdb.utils.config.ConfigInfo import ConfigInfo, getSiteId
from wwpdb.apps.seqmodule.webapp.SeqModWebRequest import SeqModInputRequest


class AlignmentStatisticsTests(unittest.TestCase):
    def setUp(self):
        #
        self.__verbose = True
        self.__lfh = sys.stderr
        self.__pathExamplesRel = "../tests"
        self.__pathExamples = os.path.abspath(self.__pathExamplesRel)
        #
        # self.__exampleFileList    = ['1cbs.cif','3rer.cif','rcsb056751.cif']
        self.__exampleFileList = ["rcsb056751.cif"]
        self.__dsList = ["D_000000"]
        #
        self.__setup()

    def __setup(self, sessionId="seq-alignstats-tests"):
        """Simulate the web application environment for managing session storage of
        temporaty data files.
        """
        self.__siteId = getSiteId(defaultSiteId="WWPDB_DEPOLY_TEST")
        self.__cI = ConfigInfo(self.__siteId)
        self.__topPath = self.__cI.get("SITE_WEB_APPS_TOP_PATH")
        self.__topSessionPath = self.__cI.get("SITE_WEB_APPS_TOP_SESSIONS_PATH")
        #
        self.__reqObj = SeqModInputRequest({}, verbose=self.__verbose, log=self.__lfh)
        self.__templatePath = os.path.join(self.__topPath, "htdocs", "seqmodule")
        self.__reqObj.setValue("TopSessionPath", self.__topSessionPath)
        self.__reqObj.setValue("TemplatePath", self.__templatePath)
        self.__reqObj.setValue("TopPath", self.__topPath)
        self.__reqObj.setValue("WWPDB_SITE_ID", self.__siteId)
        os.environ["WWPDB_SITE_ID"] = self.__siteId

        if sessionId is not None:
            self.__reqObj.setValue("sessionid", sessionId)

        # self.__sessionId = self.__reqObj.getSessionId()
        self.__sessionObj = self.__reqObj.newSessionObj()
        self.__sessionPath = self.__sessionObj.getPath()

    def testSearchAndAssembleFromUpload(self):
        """Test search each entity sequence against appropriate reference sequence database
        service storing the matching sequences.

        Using upload file source.
        """
        startTime = time.time()
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__, inspect.currentframe().f_code.co_name, time.strftime("%Y %m %d %H:%M:%S", time.localtime())))

        try:
            pI = PathInfo(siteId=self.__siteId, sessionPath=self.__sessionPath, verbose=self.__verbose, log=self.__lfh)
            for f in self.__exampleFileList:
                (idCode, _fExt) = os.path.splitext(f)
                pdbxFilePath = pI.getModelPdbxFilePath(idCode, fileSource="session")
                inpFilePath = os.path.join(self.__pathExamples, f)
                shutil.copyfile(inpFilePath, pdbxFilePath)
                self.__reqObj.setValue("identifier", idCode)
                #
                sda = SequenceDataAssemble(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
                sda.doAssemble()
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write(
            "\nCompleted %s %s at %s (%.2f seconds)\n"
            % (self.__class__.__name__, inspect.currentframe().f_code.co_name, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        )

    def testSearchAndAssembleFromArchive(self):
        """Test search each entity sequence against appropriate reference sequence database
        service storing the matching sequences.

        Using archive file source.
        """
        startTime = time.time()
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__, inspect.currentframe().f_code.co_name, time.strftime("%Y %m %d %H:%M:%S", time.localtime())))

        try:
            for dsId in self.__dsList:
                self.__reqObj.setValue("identifier", dsId)
                #
                sda = SequenceDataAssemble(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
                sda.doAssemble()
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write(
            "\nCompleted %s %s at %s (%.2f seconds)\n"
            % (self.__class__.__name__, inspect.currentframe().f_code.co_name, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        )

    def testAlignStatsFromUpload(self):
        """Test compute the sequence alignment statistics between author sequence and coordinate sequences
        and reference database sequences.

        Uses the sequence data prepared by the testSearchAndAssemblefromUpload() test above.
        """
        startTime = time.time()
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__, inspect.currentframe().f_code.co_name, time.strftime("%Y %m %d %H:%M:%S", time.localtime())))

        try:
            for f in self.__exampleFileList:
                (idCode, _fExt) = os.path.splitext(f)
                self.__reqObj.setValue("identifier", idCode)
                alstat = AlignmentStatistics(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
                alstat.doUpdate()
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write(
            "\nCompleted %s %s at %s (%.2f seconds)\n"
            % (self.__class__.__name__, inspect.currentframe().f_code.co_name, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        )

    def testAlignStatsFromArchive(self):
        """Test compute the sequence alignment statistics between author sequence and coordinate sequences
        and reference database sequences.

        Uses the sequence data prepared by the testSearchAndAssembleFromArchive() test above.
        """
        startTime = time.time()
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__, inspect.currentframe().f_code.co_name, time.strftime("%Y %m %d %H:%M:%S", time.localtime())))

        try:
            for dsId in self.__dsList:
                self.__reqObj.setValue("identifier", dsId)
                alstat = AlignmentStatistics(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
                alstat.doUpdate()
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write(
            "\nCompleted %s %s at %s (%.2f seconds)\n"
            % (self.__class__.__name__, inspect.currentframe().f_code.co_name, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        )


def suiteSearchAndAssembleTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(AlignmentStatisticsTests("testSearchAndAssembleFromUpload"))
    suiteSelect.addTest(AlignmentStatisticsTests("testSearchAndAssembleFromArchive"))
    return suiteSelect


def suiteAlignmentTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(AlignmentStatisticsTests("testAlignStatsFromUpload"))
    suiteSelect.addTest(AlignmentStatisticsTests("testAlignStatsFromArchive"))
    return suiteSelect


if __name__ == "__main__":
    if True:  # pylint: disable=using-constant-test
        mySuite = suiteSearchAndAssembleTests()
        unittest.TextTestRunner(verbosity=2).run(mySuite)
    if True:  # pylint: disable=using-constant-test
        mySuite = suiteAlignmentTests()
        unittest.TextTestRunner(verbosity=2).run(mySuite)
