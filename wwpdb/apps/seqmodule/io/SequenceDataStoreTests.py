##
# File:    SequenceDataStoreTests.py
# Date:    16-Dec-2009
#
# Updates:
#   20-Apr-2010 jdw Ported to module seqmodule.
#   27-Feb-2013 jdw updated test cases for extensions in sequence identification.
##
"""
Test cases for sequence data management.

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

from wwpdb.utils.config.ConfigInfo import ConfigInfo, getSiteId
from wwpdb.apps.seqmodule.webapp.SeqModWebRequest import SeqModInputRequest
from wwpdb.apps.seqmodule.io.SequenceDataStore import SequenceDataStore
from wwpdb.apps.seqmodule.util.SequenceExamples import SequenceExamples


class SequenceDataStoreTests(unittest.TestCase):
    def setUp(self):
        self.__lfh = sys.stdout
        self.__verbose = True
        # Load up some test data -
        self.__sE = SequenceExamples()
        self.refSL = {}
        self.xyzSL = {}
        self.authSL = {}
        self.refFD = {}
        for sid in ["A", "B"]:
            self.refSL[sid] = self.__sE.getRefSequenceWithIndexList(sid)
            self.authSL[sid] = self.__sE.getAuthSequenceWithIndexList(sid)
            self.xyzSL[sid] = self.__sE.getXyzSequenceWithIndexList(sid)
            self.refFD[sid] = self.__sE.getRefFeatureDict(sid)

        #
        # Create additional auth test sequences with random insertions and deletions
        #
        for sid in ["C", "D", "E", "F", "G"]:
            self.refSL[sid] = self.__sE.getRefSequenceWithIndexList("A")
            self.authSL[sid] = self.__sE.getAuthSequenceTestWithIndexList("A")
            self.xyzSL[sid] = self.__sE.getXyzSequenceWithIndexList("A")
            self.refFD[sid] = self.__sE.getRefFeatureDict("A")

        self.idList = ["A", "B", "C", "D", "E", "F", "G"]
        self.entryDict = {}
        self.entryDict["PDB_ID"] = self.__sE.getIdCode()
        #
        self.groupDict = {}
        self.groupDict[1] = ["A", "C", "D", "E", "F", "G"]
        self.groupDict[2] = ["B"]
        #
        # Setup the session environment
        self.__setup()
        #
        # The following identifer must be set to support the common tool naming functions in PathInfo()
        #
        self.__reqObj.setValue("identifier", "testcase")

    def __setup(self, sessionId="seqdatastore-testing"):
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

        # self.__sessionId = self.__reqObj.getSessionId()
        self.__sObj = self.__reqObj.newSessionObj()  # pylint: disable=unused-private-member

    def tearDown(self):
        pass

    def testStoreExampleData(self):
        """Test creating store, serializing and deserializing data."""
        startTime = time.time()
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__, inspect.currentframe().f_code.co_name, time.strftime("%Y %m %d %H:%M:%S", time.localtime())))

        try:
            sda = SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            sda.reset()
            for mid in self.idList:
                sda.setSequence(self.xyzSL[mid], mid, "xyz", partId=1, altId=1, version=1)
                sda.setSequence(self.authSL[mid], mid, "auth", partId=1, altId=1, version=1)
                for altId in range(1, 10):
                    sda.setFeature(self.refFD[mid], mid, "ref", partId=1, altId=altId, version=1)
                    sda.setSequence(self.refSL[mid], mid, "ref", partId=1, altId=altId, version=1)

            #
            for k, v in self.groupDict.items():
                sda.setGroup(k, v)
            for k, v in self.entryDict.items():
                sda.setEntryDetail(k, v)

            sda.serialize()
            #
            sda.dump(self.__lfh)
            sda.reset()
            sda.deserialize()
            sda.dump(self.__lfh)

        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write(
            "\nCompleted %s %s at %s (%.2f seconds)\n"
            % (self.__class__.__name__, inspect.currentframe().f_code.co_name, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        )

    def testReadExampleData(self):
        """Test read existing store -"""
        startTime = time.time()
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__, inspect.currentframe().f_code.co_name, time.strftime("%Y %m %d %H:%M:%S", time.localtime())))

        try:
            sda = SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            sda.dump(self.__lfh)

            sTypes = sda.getSequenceTypes()
            for st in sTypes:
                self.__lfh.write("Data for sequence type: %5s\n" % st)
                ids = sda.getIds(dataType="sequence", seqType=st)
                for sid in ids:
                    sL = sda.getSequence(sid, st, partId=1, altId=1, version=1)
                    self.__lfh.write("Sequence type %5s id %5s length %5d\n" % (st, sid, len(sL)))

        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write(
            "\nCompleted %s %s at %s (%.2f seconds)\n"
            % (self.__class__.__name__, inspect.currentframe().f_code.co_name, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        )

    def testAccessMethodsExampleData(self):
        """SequenceDataStore() accessor tests -"""
        startTime = time.time()
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__, inspect.currentframe().f_code.co_name, time.strftime("%Y %m %d %H:%M:%S", time.localtime())))

        try:
            sda = SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)

            sTypes = sda.getSequenceTypes()
            for st in sTypes:
                self.__lfh.write("Data for sequence type: %5s\n" % st)
                ids = sda.getIds(dataType="sequence", seqType=st)
                for mid in ids:
                    verList = sda.getVersionIds(seqId=mid, partId=1, altId=1, dataType="sequence", seqType=st)
                    self.__lfh.write(" Sequence id : %5s version length %d\n" % (mid, len(verList)))
                    for ver in verList:
                        sL = sda.getSequence(mid, st, partId=1, altId=1, version=ver)
                        self.__lfh.write("   Sequence type %5s id %5s version %d length %5d\n" % (st, mid, ver, len(sL)))
                    altList = sda.getAlternativeIds(seqId=mid, dataType="sequence", seqType=st)
                    self.__lfh.write(" Sequence id : %5s alternative list  length %d\n" % (mid, len(altList)))
                    for altId in altList:
                        sL = sda.getSequence(seqId=mid, seqType=st, partId=1, altId=altId, version=1)
                        self.__lfh.write("   Sequence type %5s id %5s altId %s length %5d\n" % (st, mid, altId, len(sL)))

            sda.dump(self.__lfh)
            sda.filterIndex(seqId="F", dataType="sequence", seqType="ref")
            sda.dump(self.__lfh)
            sda.filterIndex(seqId="A", dataType="sequence", seqType="ref")
            sda.dump(self.__lfh)

        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write(
            "\nCompleted %s %s at %s (%.2f seconds)\n"
            % (self.__class__.__name__, inspect.currentframe().f_code.co_name, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        )


def suiteWriteAndReadTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(SequenceDataStoreTests("testStoreExampleData"))
    suiteSelect.addTest(SequenceDataStoreTests("testReadExampleData"))
    suiteSelect.addTest(SequenceDataStoreTests("testAccessMethodsExampleData"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = suiteWriteAndReadTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
