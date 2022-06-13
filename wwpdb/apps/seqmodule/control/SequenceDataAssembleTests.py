##
# File:    SequenceDataAssembleTests.py
# Date:    26-Feb-2013
#
# Updates:
#
##
"""
Test cases for search and assembly of sequence data for input to sequence editor tool.

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
import inspect

from wwpdb.apps.seqmodule.control.SequenceDataAssemble_v2 import SequenceDataAssemble
from wwpdb.io.locator.PathInfo import PathInfo
from wwpdb.utils.config.ConfigInfo import ConfigInfo, getSiteId
from wwpdb.apps.seqmodule.webapp.SeqModWebRequest import SeqModInputRequest
from wwpdb.utils.dp.DataFileAdapter import DataFileAdapter
from wwpdb.apps.seqmodule.align.AlignmentStatistics import AlignmentStatistics

from wwpdb.apps.seqmodule.io.ModelSequenceUtils import ModelSequenceUtils
from wwpdb.apps.seqmodule.io.PdbxIoUtils import PdbxFileIo
from wwpdb.io.misc.FormatOut import FormatOut


class SequenceDataAssembleTests(unittest.TestCase):
    def setUp(self):
        #
        self.__verbose = True
        self.__lfh = sys.stdout
        self.__pathExamplesRel = "../tests"
        self.__pathExamples = os.path.abspath(self.__pathExamplesRel)
        self.__exampleFileType = ""
        #
        # self.__examFileList    = ['1cbs.cif','3rer.cif','rcsb056751.cif']
        #
        self.__examFileList = ["rcsb056751.cif"]
        # self.__examFileList    = ['rcsb076237.cif']
        self.__dsList = ["D_000000"]
        #

        self.mySetup()

    def mySetup(self, sessionId=None):
        """Simulate the web application environment for managing session storage of
        temporaty data files.
        """
        self.__maxRefAlign = 100
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

        #

    def testSearchAndAssembleFromUpload(self):
        """Test search each entity sequence against appropriate reference sequence database
        service storing the matching sequences.

        Using upload file source. rcsb-mmcif
        """
        startTime = time.time()
        self.__lfh.write("\n\n========================================================================================================\n")
        self.__lfh.write("Starting %s %s at %s\n" % (self.__class__.__name__, inspect.currentframe().f_code.co_name, time.strftime("%Y %m %d %H:%M:%S", time.localtime())))

        try:
            self.__exampleFileType = "rcsb-mmcif"
            pI = PathInfo(siteId=self.__siteId, sessionPath=self.__sessionPath, verbose=self.__verbose, log=self.__lfh)
            for f in self.__examFileList:
                (idCode, _fExt) = os.path.splitext(f)
                pdbxFilePath = pI.getModelPdbxFilePath(idCode, fileSource="session")
                inpFilePath = os.path.join(self.__pathExamples, f)
                self.__lfh.write("+testSearchAndAssembleFromUpload() Starting with id %s \n + input file %s \n   + pdbxfile %s\n" % (idCode, inpFilePath, pdbxFilePath))
                dfa = DataFileAdapter(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
                _ok = dfa.modelConvertToPdbx(filePath=inpFilePath, fileType=self.__exampleFileType, pdbxFilePath=pdbxFilePath)  # noqa: F841
                self.__reqObj.setValue("identifier", idCode)
                #
                if not os.access(inpFilePath, os.F_OK):
                    self.__lfh.write("+testSearchAndAssembleFromUpload() input file missing %s\n" % inpFilePath)
                    self.fail()
                    break
                if not os.access(pdbxFilePath, os.F_OK):
                    self.__lfh.write("+testSearchAndAssembleFromUpload() format conversion failed for %s\n" % idCode)
                    self.fail()
                    break
                #
                sda = SequenceDataAssemble(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
                # sda.doAssemble(fileSource="local-upload")
                sda.doAssemble()
                alstat = AlignmentStatistics(reqObj=self.__reqObj, maxRefAlign=self.__maxRefAlign, verbose=self.__verbose, log=self.__lfh)
                alstat.doUpdate()
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write(
            "\nCompleted %s %s at %s (%.2f seconds)\n"
            % (self.__class__.__name__, inspect.currentframe().f_code.co_name, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        )

    def testExtractDataFromUpload(self):
        """Test entry data extraction in preparation for sequence search

        Using upload file source. rcsb-mmcif
        """
        startTime = time.time()
        self.__lfh.write("\n\n========================================================================================================\n")
        self.__lfh.write("Starting %s %s at %s\n" % (self.__class__.__name__, inspect.currentframe().f_code.co_name, time.strftime("%Y %m %d %H:%M:%S", time.localtime())))

        try:
            self.__exampleFileType = "rcsb-mmcif"
            pI = PathInfo(siteId=self.__siteId, sessionPath=self.__sessionPath, verbose=self.__verbose, log=self.__lfh)
            for f in self.__examFileList:
                (idCode, _fExt) = os.path.splitext(f)
                pdbxFilePath = pI.getModelPdbxFilePath(idCode, fileSource="session")
                inpFilePath = os.path.join(self.__pathExamples, f)
                self.__lfh.write("+testSearchAndAssembleFromUpload() Starting with id %s \n + input file %s \n   + pdbxfile %s\n" % (idCode, inpFilePath, pdbxFilePath))
                dfa = DataFileAdapter(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
                _ok = dfa.modelConvertToPdbx(filePath=inpFilePath, fileType=self.__exampleFileType, pdbxFilePath=pdbxFilePath)  # noqa: F841
                self.__reqObj.setValue("identifier", idCode)
                #
                if not os.access(inpFilePath, os.F_OK):
                    self.__lfh.write("+testSearchAndAssembleFromUpload() input file missing %s\n" % inpFilePath)
                    self.fail()
                    break
                if not os.access(pdbxFilePath, os.F_OK):
                    self.__lfh.write("+testSearchAndAssembleFromUpload() format conversion failed for %s\n" % idCode)
                    self.fail()
                    break
                #
                c0 = PdbxFileIo(verbose=self.__verbose, log=self.__lfh).getContainer(pdbxFilePath)
                msu = ModelSequenceUtils(dataContainer=c0, verbose=self.__verbose, log=self.__lfh)
                entityD = msu.getEntitySequenceDetails()
                instD = msu.getCoordinateSequenceDetails()
                fOut = FormatOut()
                fOut.autoFormat("Entity dictionary", entityD, 3, 3)
                fOut.autoFormat("Instance  dictionary", instD, 3, 3)
                fOut.writeStream(self.__lfh)
                fOut.clear()
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
                # sda.doAssemble(fileSource="archive")
                sda.doAssemble()
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write(
            "\nCompleted %s %s at %s (%.2f seconds)\n"
            % (self.__class__.__name__, inspect.currentframe().f_code.co_name, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        )


def suiteSearchAndAssembleFromUploadTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(SequenceDataAssembleTests("testSearchAndAssembleFromUpload"))
    return suiteSelect


def suiteSearchAndAssembleFromArchiveTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(SequenceDataAssembleTests("testSearchAndAssembleFromArchive"))
    return suiteSelect


def suiteExtractDataFromUploadTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(SequenceDataAssembleTests("testExtractDataFromUpload"))
    return suiteSelect


if __name__ == "__main__":
    if True:  # pylint: disable=using-constant-test
        mySuite = suiteSearchAndAssembleFromUploadTests()
        unittest.TextTestRunner(verbosity=2).run(mySuite)

    if False:  # pylint: disable=using-constant-test
        mySuite = suiteSearchAndAssembleFromArchiveTests()
        unittest.TextTestRunner(verbosity=2).run(mySuite)

    if False:  # pylint: disable=using-constant-test
        mySuite = suiteExtractDataFromUploadTests()
        unittest.TextTestRunner(verbosity=2).run(mySuite)
