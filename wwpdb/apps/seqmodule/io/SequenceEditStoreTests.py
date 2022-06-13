##
# File:    SequenceEditStoreTests.py
# Date:    13-Feb-2010
#
# Updates:
# 20-Apr-2010 jdw  Ported to module seqmodule.
# 23-Feb-2013 jdw  Updated session module and SequenceDataStore() calling interface.
##
"""
Test cases for sequence editing operations.

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

from wwpdb.apps.seqmodule.io.SequenceEditStore import SequenceEditStore, SequenceEdit
from wwpdb.apps.seqmodule.io.SequenceDataStore import SequenceDataStore
from wwpdb.apps.seqmodule.util.SequenceExamples import SequenceExamples
from wwpdb.apps.seqmodule.util.SequenceLabel import ResidueLabel


class SequenceEditStoreTests(unittest.TestCase):
    def setUp(self):
        self.__verbose = True
        self.__lfh = sys.stderr
        #
        # Load up some test data -
        self.__sE = SequenceExamples()

        self.refSL = {}
        self.xyzSL = {}
        self.authSL = {}
        self.refFD = {}
        for cid in ["A", "B"]:
            self.refSL[cid] = self.__sE.getRefSequenceWithIndexList(cid)
            self.authSL[cid] = self.__sE.getAuthSequenceWithIndexList(cid)
            self.xyzSL[cid] = self.__sE.getXyzSequenceWithIndexList(cid)
            self.refFD[cid] = self.__sE.getRefFeatureDict(cid)

        #
        # Create additional auth test sequences with random insertions and deletions
        #
        for cid in ["C", "D", "E", "F", "G"]:
            self.refSL[cid] = self.__sE.getRefSequenceWithIndexList("A")
            self.authSL[cid] = self.__sE.getAuthSequenceTestWithIndexList("A")
            self.xyzSL[cid] = self.__sE.getXyzSequenceWithIndexList("A")
            self.refFD[cid] = self.__sE.getRefFeatureDict("A")

        self.idList = ["A", "B", "C", "D", "E", "F", "G"]

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
        self.__sObj = self.__reqObj.newSessionObj()

    def tearDown(self):
        pass

    def testStoreAndEdit(self):
        """ """
        startTime = time.time()
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__, inspect.currentframe().f_code.co_name, time.strftime("%Y %m %d %H:%M:%S", time.localtime())))

        try:
            sda = SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            for sid in self.idList:
                sda.setSequence(self.refSL[sid], sid, "ref", version=1)
                sda.setSequence(self.xyzSL[sid], sid, "xyz", version=1)
                sda.setSequence(self.authSL[sid], sid, "auth", version=1)
                sda.setFeature(self.refFD[sid], sid, "ref", version=1)
            #
            sda.serialize()
            #
            sda.dump(self.__lfh)
            sda.reset()
            sda.deserialize()
            sda.dump(self.__lfh)

            ses = SequenceEditStore(sessionObj=self.__sObj, verbose=self.__verbose)
            resLabel = ResidueLabel()
            for cid in self.idList:
                seqType = "xyz"
                seqAltId = 1
                seqVersion = 1
                seqInstId = cid
                aL = self.xyzSL[cid]
                ibeg = 10
                iend = 20
                for sPos in range(ibeg, iend):
                    rT = aL[sPos]

                    # idS contains: type(ref,auth,coordinate) + chain_id + 3-letter-code +
                    #               orginal residue label index + position in alignment  + position in sequence
                    #
                    resLabel.set(
                        seqType=seqType, seqInstId=seqInstId, seqAltId=seqAltId, seqVersion=seqVersion, residueCode3=rT[0], residueLabelIndex=rT[1], alignIndex=sPos, seqIndex=sPos
                    )
                    idS = resLabel.pack()
                    if seqType == "auth" or seqType == "xyz":
                        editOpLast = ses.getLastEditOp()
                        editOpNext = int(editOpLast) + 1
                        newVal = ["TRP"]
                        priorVal = rT[0]
                        sE = SequenceEdit(self.__verbose)
                        sE.setValueNew(newVal)
                        sE.setValuePrevious(priorVal)
                        sE.setEditType("replace")
                        sE.setTargetElementId(idS)
                        sE.setEditOpId(editOpNext)
                        ses.storeEdit(sE)
                        editOpNext += 1
                        sE = SequenceEdit(self.__verbose)
                        sE.setValueNew("Preliminary comment about " + idS)
                        sE.setValuePrevious("none")
                        sE.setEditType("details")
                        sE.setTargetElementId(idS)
                        sE.setEditOpId(editOpNext)
                        ses.storeEdit(sE)
                        #
                        editOpNext += 1
                        sE = SequenceEdit(self.__verbose)
                        sE.setValueNew("Final comment about " + idS)
                        sE.setValuePrevious("none")
                        sE.setEditType("details")
                        sE.setTargetElementId(idS)
                        sE.setEditOpId(editOpNext)
                        ses.storeEdit(sE)

            ses.printIt(self.__lfh)
            dD = ses.getDetailsByTarget()
            for k, v in dD.items():
                self.__lfh.write("Residue id %10s: comment  %s\n" % (k, v))

        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write(
            "\nCompleted %s %s at %s (%.2f seconds)\n"
            % (self.__class__.__name__, inspect.currentframe().f_code.co_name, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        )

    def testRestore(self):
        """ """
        startTime = time.time()
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__, inspect.currentframe().f_code.co_name, time.strftime("%Y %m %d %H:%M:%S", time.localtime())))

        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % inspect.currentframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            sda = SequenceDataStore(reqObj=self.__reqObj, verbose=True, log=self.__lfh)
            sda.reset()
            sda.deserialize()
            sda.dump(self.__lfh)

            sTypes = sda.getSequenceTypes()
            for st in sTypes:
                ids = sda.getIds(seqType=st, dataType="sequence")
                for cid in ids:
                    sL = sda.getSequence(cid, st, version=1)
                    self.__lfh.write("Sequence type %5s id %5s length %5d\n" % (st, cid, len(sL)))

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
    suiteSelect.addTest(SequenceEditStoreTests("testStoreAndEdit"))
    suiteSelect.addTest(SequenceEditStoreTests("testRestore"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = suiteWriteAndReadTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
