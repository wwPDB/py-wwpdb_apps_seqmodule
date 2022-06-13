##
# File:    ModelSequenceUtilsTests.py
# Date:    21-Feb-2013
#
# Updates:
#  23-Feb-2013 jdw generalize to multi-part entities added test case rcsb056751.
#  28-Feb-2013 jdw update
##
"""
Test cases for extracting entity, instance, and sequence details from deposited entry files.

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

from wwpdb.apps.seqmodule.io.ModelSequenceUtils import ModelSequenceUtils
from wwpdb.apps.seqmodule.io.PdbxIoUtils import PdbxFileIo
from wwpdb.io.misc.FormatOut import FormatOut
from wwpdb.io.locator.PathInfo import PathInfo
from wwpdb.utils.config.ConfigInfo import ConfigInfo, getSiteId
from wwpdb.apps.seqmodule.webapp.SeqModWebRequest import SeqModInputRequest
from wwpdb.utils.dp.DataFileAdapter import DataFileAdapter


class ModelSequenceUtilsTests(unittest.TestCase):
    def setUp(self):
        #
        self.__verbose = True
        self.__lfh = sys.stdout
        self.__pathExamplesRel = "../tests"
        self.__pathExamples = os.path.abspath(self.__pathExamplesRel)
        #
        # self.__examFileList    = ['1cbs.cif','3rer.cif','rcsb056751.cif','rcsb076237.cif', 'rcsb057171.cif']

        self.__exampleFileType = "rcsb-mmcif"
        self.__exampleFileList = [
            "rcsb052365.cif",
            "rcsb053659.cif",
            "rcsb055453.cif",
            "rcsb056215.cif",
            "rcsb056751.cif",
            "rcsb057171.cif",
            "rcsb057584.cif",
            "rcsb057750.cif",
            "rcsb076237.cif",
            "rcsb095269.cif",
        ]
        self.__exampleFileList = ["rcsb056751.cif", "rcsb076237.cif"]

        self.mySetup()

    def mySetup(self, sessionId=None):
        """Simulate the web application environment for managing session storage of
        temporaty data files.
        """
        self.__siteId = getSiteId(defaultSiteId="WWPDB_DEPLOY_TEST")
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

        self.__sessionObj = self.__reqObj.newSessionObj()
        self.__sessionPath = self.__sessionObj.getPath()

    def testGetModelSequence(self):
        """Test extraction of author and coordinate sequence data."""
        startTime = time.time()
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__, inspect.currentframe().f_code.co_name, time.strftime("%Y %m %d %H:%M:%S", time.localtime())))
        try:
            pI = PathInfo(siteId=self.__siteId, sessionPath=self.__sessionPath, verbose=self.__verbose, log=self.__lfh)
            for f in self.__exampleFileList:
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

                c0 = PdbxFileIo(verbose=self.__verbose, log=self.__lfh).getContainer(pdbxFilePath)
                msu = ModelSequenceUtils(dataContainer=c0, verbose=self.__verbose, log=self.__lfh)
                entityD = msu.getEntitySequenceDetails()
                instD = msu.getCoordinateSequenceDetails()
                #
                depSeqAssign, seqAssign = msu.getSequenceAssignmentDetails()
                #
                fOut = FormatOut()
                fOut.autoFormat("Entity dictionary", entityD, 3, 3)
                fOut.autoFormat("Instance  dictionary", instD, 3, 3)
                fOut.autoFormat("Depositor sequence assignment dictionary", depSeqAssign, 3, 3)
                fOut.autoFormat("Archive sequence assignment dictionary", seqAssign, 3, 3)

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


def suiteGetModelSeqTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ModelSequenceUtilsTests("testGetModelSequence"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = suiteGetModelSeqTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
