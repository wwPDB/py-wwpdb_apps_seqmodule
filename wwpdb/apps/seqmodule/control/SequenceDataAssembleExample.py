##
# File:  SequenceDataAssembleExample.py
# Date:  05-Mar-2013
#
# Updates:
#
##
"""
Assemble sequence and other required data to start/restart the sequence editor tool.
"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.09"

import sys
import os
import shutil

from wwpdb.apps.seqmodule.io.SequenceDataStore import SequenceDataStore
from wwpdb.apps.seqmodule.util.SequenceExamples import SequenceExamples


class SequenceDataAssembleExample(object):
    """
    This class loads the sequence data store with example data
    from the SequenceExamples() class.

    This class is provided primarily for testing and development.

    Storage model - imported data is loaded into the sequence data store
                    where it is managed by the SequenceDataStore() class.

    """

    def __init__(self, reqObj=None, verbose=False, log=sys.stderr):
        self.__verbose = verbose
        self.__reqObj = reqObj
        self.__sessionObj = None
        self.__lfh = log
        self.__localTest = False
        #
        self.__sessionPath = "."
        #
        self.__setup()

    def __setup(self):
        try:
            self.__sessionObj = self.__reqObj.getSessionObj()
            self.__sessionPath = self.__sessionObj.getPath()
            if self.__verbose:
                self.__lfh.write("+SequenceDataAssembleExample.__setup() - session id %s\n" % self.__sessionObj.getId())
                self.__lfh.write("+SequenceDataAssembleExample.__setup() - session path %s\n" % self.__sessionObj.getPath())
        except:  # noqa: E722 pylint: disable=bare-except
            pass

    def __copyExampleStructureFile(self, did, fileType="pdb"):

        fn = str(did).lower() + "." + fileType

        if self.__verbose:
            self.__lfh.write("+SequenceDataAssembleExample.__copyExampleStructureFile () - copying file name  %s\n" % fn)
        dst = os.path.join(self.__sessionPath, fn)
        topPath = self.__reqObj.getValue("TopPath")
        src = os.path.join(topPath, "xyz", fn)
        if self.__verbose:
            self.__lfh.write("+SequenceDataAssembleExample.__copyExampleStructureFile() - src  %s dst %s\n" % (src, dst))
        try:
            shutil.copyfile(src, dst)
        except:  # noqa: E722 pylint: disable=bare-except
            pass
        return dst

    def doImport(self):
        return self.__loadExampleData()

    def __loadExampleData(self):
        """Store example sequence data in a SequenceDataStore repository
        in the current session directory.
        """
        #
        if self.__verbose:
            self.__lfh.write("+SequenceDataAssembleExample.__loadExample()  sessionId %s\n" % (self.__sessionObj.getId()))

        # Load up some test data -
        sE = SequenceExamples()

        refSL = {}
        xyzSL = {}
        authSL = {}
        refFD = {}
        xyzFD = {}
        authFD = {}
        for did in ["A", "B"]:
            refSL[did] = sE.getRefSequenceWithIndexList(did)
            authSL[did] = sE.getAuthSequenceWithIndexList(did)
            xyzSL[did] = sE.getXyzSequenceWithIndexList(did)
            refFD[did] = sE.getRefFeatureDict(did)
            authFD[did] = sE.getAuthFeatureDict(did)
            xyzFD[did] = sE.getXyzFeatureDict(did)

        if self.__localTest:
            #
            # Create additional auth test sequences with random insertions and deletions
            #
            for did in ["C", "D", "E", "F", "G"]:
                refSL[did] = sE.getRefSequenceTestWithIndexList("A")
                authSL[did] = sE.getAuthSequenceTestWithIndexList("A")
                xyzSL[did] = sE.getXyzSequenceTestWithIndexList("A")
                refFD[did] = sE.getRefFeatureDict("A")
                authFD[did] = sE.getAuthFeatureDict("A")
                xyzFD[did] = sE.getXyzFeatureDict("A")

            idList = ["A", "B", "C", "D", "E", "F", "G"]
            groupDict = {}
            groupDict[1] = ["A", "C", "D", "E", "F", "G"]
            groupDict[2] = ["B"]

        else:
            idList = ["A", "B"]
            groupDict = {}
            groupDict[1] = ["A"]
            groupDict[2] = ["B"]

        #
        sds = SequenceDataStore(reqObj=self.__sessionObj, verbose=self.__verbose, log=self.__lfh)
        sds.reset()
        for did in idList:
            sds.setSequence(authSL[did], did, "auth", altId=1, version=1)
            sds.setFeature(authFD[did], did, "auth", altId=1, version=1)
            sds.setSequence(xyzSL[did], did, "xyz", altId=1, version=1)
            sds.setFeature(xyzFD[did], did, "xyz", altId=1, version=1)
            for altId in range(1, 10):
                sds.setSequence(refSL[did], did, "ref", altId=altId, version=1)
                sds.setFeature(refFD[did], did, "ref", altId=altId, version=1)
        #
        #
        entryDict = {}
        entryDict["PDB_ID"] = sE.getIdCode()
        entryDict["RCSB_ID"] = "RCSB054485"
        #
        for k, v in groupDict.items():
            sds.setGroup(k, v)
        for k, v in entryDict.items():
            sds.setEntryDetail(k, v)

        sds.serialize()

        self.__copyExampleStructureFile(did=sE.getIdCode(), fileType="pdb")

        return {}
