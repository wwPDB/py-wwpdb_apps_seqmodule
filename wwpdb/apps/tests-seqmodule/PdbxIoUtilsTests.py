##
# File:    PdbxIoUtilsTests.py
# Date:    21-Feb-2013
#
# Updates:
#
#  05-Mar-2013 jdw  tests add getSequenceFeaturesFromAtomSite()
##
"""
Test cases for extracting entity, chain, and sequence details from deposited entry files.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.07"

import os
import os.path
import sys
import inspect
import traceback
import unittest

from wwpdb.apps.seqmodule.io.PdbxIoUtils import ModelFileIo, PdbxFileIo


class PdbxIoUtilsTests(unittest.TestCase):
    def setUp(self):
        #
        self.__verbose = True
        self.__lfh = sys.stdout

        # Old examples -
        HERE = os.path.abspath(os.path.dirname(__file__))
        self.__pathExamples = os.path.join(HERE, "data")
        #
        self.__examFileList = ["4ec0.cif", "3rer.cif"]

    def tearDown(self):
        pass

    def testSequenceDictionary(self):
        """Test extraction of author and coordinate sequence details."""
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % inspect.currentframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            for f in self.__examFileList:
                fN = os.path.join(self.__pathExamples, f)
                c0 = PdbxFileIo(verbose=self.__verbose, log=self.__lfh).getContainer(fN)
                sdf = ModelFileIo(dataContainer=c0, verbose=self.__verbose, log=self.__lfh)
                self.__lfh.write("File %s\n" % fN)
                for eId in sdf.getPolymerEntityList("polymer"):
                    self.__lfh.write("  entity id= %s\n" % eId)
                    cList = sdf.getPdbChainIdList(eId)
                    kd = {}
                    for cId in cList:
                        kd["chainId"] = cId
                        sqEntity = sdf.getEntityCanSequence1Auth(kd)
                        sqXyz = sdf.getSequence1Xyz(cId)
                        self.__lfh.write("  Chain Id %s entity sequence: %s\n" % (cId, str(sqEntity).strip()))
                        self.__lfh.write("  Chain Id %s   xyz  sequence: %s\n" % (cId, str(sqXyz).strip()))

        except Exception as _e:  # noqa: F841
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testGetEntityCounts(self):
        """Get entity counts  -"""
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % inspect.currentframe().f_code.co_name)
        self.__lfh.write(" ------------------------\n")
        try:
            for f in self.__examFileList:
                fN = os.path.join(self.__pathExamples, f)
                c0 = PdbxFileIo(verbose=self.__verbose, log=self.__lfh).getContainer(fN)
                sdf = ModelFileIo(dataContainer=c0, verbose=self.__verbose, log=self.__lfh)
                self.__lfh.write("------------------------------------------------------\n")
                self.__lfh.write("  polymers         = %d\n" % sdf.getEntityCount("polymer"))
                self.__lfh.write("  non-polymers     = %d\n" % sdf.getEntityCount("non-polymer"))
                self.__lfh.write("  (L) polypeptide polymers = %d\n" % sdf.getPolyPeptideLEntityCount())
                self.__lfh.write("  (D) polypeptide polymers = %d\n" % sdf.getPolyPeptideDEntityCount())
        except Exception as _e:  # noqa: F841
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testGetEntityLists(self):
        """Get entity lists  -"""
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % inspect.currentframe().f_code.co_name)
        self.__lfh.write(" ------------------------\n")
        try:
            for f in self.__examFileList:
                fN = os.path.join(self.__pathExamples, f)
                c0 = PdbxFileIo(verbose=self.__verbose, log=self.__lfh).getContainer(fN)
                sdf = ModelFileIo(dataContainer=c0, verbose=self.__verbose, log=self.__lfh)
                self.__lfh.write("------------------------------------------------------\n")
                self.__lfh.write("File %s id  name : %s\n" % (fN, sdf.getDbCode("PDB")))
                self.__lfh.write("  polymers         = %d\n" % sdf.getEntityCount("polymer"))
                #
                for eId in sdf.getPolymerEntityList("polymer"):
                    self.__lfh.write("  entity id= %s\n" % eId)
                self.__lfh.write("  non-polymers     = %d\n" % sdf.getEntityCount("non-polymer"))
                for eId in sdf.getPolymerEntityList("non-polymer"):
                    self.__lfh.write("  entity id= %s\n" % eId)

        except Exception as _e:  # noqa: F841
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testGetInstIdLists(self):
        """Get entity and chain correspondences -"""
        self.__lfh.write("\n-- -------------------- ")
        self.__lfh.write("Starting test function  %s" % inspect.currentframe().f_code.co_name)
        self.__lfh.write(" -----------------------\n")
        try:
            for f in self.__examFileList:
                fN = os.path.join(self.__pathExamples, f)
                c0 = PdbxFileIo(verbose=self.__verbose, log=self.__lfh).getContainer(fN)
                sdf = ModelFileIo(dataContainer=c0, verbose=self.__verbose, log=self.__lfh)
                self.__lfh.write("------------------------------------------------------\n")
                self.__lfh.write("File %s id  name : %s\n" % (fN, sdf.getDbCode("PDB")))
                self.__lfh.write("  polymers         = %d\n" % sdf.getEntityCount("polymer"))
                for eId in sdf.getPolymerEntityList("polymer"):
                    self.__lfh.write("  entity id= %s\n" % eId)
                    cList = sdf.getPdbChainIdList(eId)
                    for cId in cList:
                        self.__lfh.write("  chain id = %s\n" % cId)
                self.__lfh.write("------------------------------------------------------\n")
                ed = sdf.getPolymerEntityChainDict()
                for eId, iList in ed.items():
                    for instId in iList:
                        self.__lfh.write("   Entity %s  instance id = %s\n" % (eId, instId))
                #
                cd = sdf.getChainPolymerEntityDict()
                for k, v in cd.items():
                    self.__lfh.write("   Instance id = %s entity %s  \n" % (k, v))

        except Exception as _e:  # noqa: F841
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testGetSeq3Lists(self):
        """Get 3-letter-code sequences"""
        self.__lfh.write("\n-------------------- ")
        self.__lfh.write("Starting test function  %s" % inspect.currentframe().f_code.co_name)
        self.__lfh.write(" --------------------\n")
        try:
            for f in self.__examFileList:
                fN = os.path.join(self.__pathExamples, f)
                c0 = PdbxFileIo(verbose=self.__verbose, log=self.__lfh).getContainer(fN)
                sdf = ModelFileIo(dataContainer=c0, verbose=self.__verbose, log=self.__lfh)
                self.__lfh.write("------------------------------------------------------\n")
                self.__lfh.write("File %s id  name : %s\n" % (fN, sdf.getDbCode("PDB")))

                self.__lfh.write("  polymers         = %d\n" % sdf.getEntityCount("polymer"))
                for eId in sdf.getPolymerEntityList("polymer"):
                    self.__lfh.write("  entity id= %s\n" % eId)
                    iList = sdf.getPdbChainIdList(eId)
                    for instId in iList:
                        self.__lfh.write("------------------------------------------------------\n")
                        self.__lfh.write("  Sequence for chain id = %s\n" % instId)
                        seq3List = sdf.getSequence3AlignList(instId)
                        for monTup in seq3List:
                            self.__lfh.write("Index %s Entity monomer  %s  XYZ monomer = %s Pdb Index %s \n" % (monTup[2], monTup[1], monTup[0], monTup[3]))
        except Exception as _e:  # noqa: F841
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testEntitySourceDetails(self):
        """Get 3-letter-code sequences and print any alignment conflicts -"""
        self.__lfh.write("\n------------------ ")
        self.__lfh.write("Starting test function  %s" % inspect.currentframe().f_code.co_name)
        self.__lfh.write(" ------------------\n")
        try:
            for f in self.__examFileList:

                fN = os.path.join(self.__pathExamples, f)
                c0 = PdbxFileIo(verbose=self.__verbose, log=self.__lfh).getContainer(fN)
                sdf = ModelFileIo(dataContainer=c0, verbose=self.__verbose, log=self.__lfh)
                self.__lfh.write("------------------------------------------------------\n")
                self.__lfh.write("File %s id  name : %s\n" % (fN, sdf.getDbCode("PDB")))
                self.__lfh.write("  polymers         = %d\n" % sdf.getEntityCount("polymer"))
                for eId in sdf.getPolymerEntityList("polymer"):
                    sL = []
                    self.__lfh.write("  Entity id= %s\n" % eId)
                    seqOneLetterCode = sdf.getSequence(eId)
                    #
                    seqLength = len(seqOneLetterCode)
                    catList = ["entity_src_gen", "entity_src_nat", "pdbx_entity_src_syn"]
                    #
                    for cat in catList:
                        sL.extend(sdf.getSourceDetailsList(eId, sourceCategoryName=cat, seqLength=seqLength))
                    #
                    self.__lfh.write("Source details entity %s:\n" % eId)
                    for sD in sL:
                        self.__lfh.write("  -------------------------------  Part -------------------\n")
                        for k, v in sD.items():
                            self.__lfh.write("  key = %30s  value=%s\n" % (k, v))
        except Exception as _e:  # noqa: F841
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testGetSeq3AlignConficts(self):
        """Get 3-letter-code sequences and print any alignment conflicts -"""
        self.__lfh.write("\n------------------ ")
        self.__lfh.write("Starting test function  %s" % inspect.currentframe().f_code.co_name)
        self.__lfh.write(" ------------------\n")
        try:
            for f in self.__examFileList:
                fN = os.path.join(self.__pathExamples, f)
                c0 = PdbxFileIo(verbose=self.__verbose, log=self.__lfh).getContainer(fN)
                sdf = ModelFileIo(dataContainer=c0, verbose=self.__verbose, log=self.__lfh)
                self.__lfh.write("------------------------------------------------------\n")
                self.__lfh.write("File %s id  name : %s\n" % (fN, sdf.getDbCode("PDB")))
                self.__lfh.write("  polymers         = %d\n" % sdf.getEntityCount("polymer"))
                for eId in sdf.getPolymerEntityList("polymer"):
                    self.__lfh.write("  entity id= %s\n" % eId)
                    iList = sdf.getPdbChainIdList(eId)
                    for instId in iList:
                        self.__lfh.write("------------------------------------------------------\n")
                        self.__lfh.write("  Sequence conflicts for chain id = %s\n" % instId)
                        seq3List = sdf.getSequence3AlignList(instId)
                        numC = 0
                        for monTup in seq3List:
                            if monTup[1] != monTup[0]:
                                numC += 1
                                self.__lfh.write("Conflict at index %s Entity monomer  %s  XYZ monomer = %s Pdb Index %s \n" % (monTup[2], monTup[1], monTup[0], monTup[3]))
                        if numC == 0:
                            self.__lfh.write("No conflicts\n")
        except Exception as _e:  # noqa: F841
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testAAAGetSeqFromAtomSite(self):
        """Get sequences, occupancy and disorder from atom site -"""
        self.__lfh.write("\n-------------------- ")
        self.__lfh.write("Starting test function  %s" % inspect.currentframe().f_code.co_name)
        self.__lfh.write(" --------------------\n")
        try:
            for f in self.__examFileList:
                fN = os.path.join(self.__pathExamples, f)
                c0 = PdbxFileIo(verbose=self.__verbose, log=self.__lfh).getContainer(fN)
                sdf = ModelFileIo(dataContainer=c0, verbose=self.__verbose, log=self.__lfh)
                self.__lfh.write("------------------------------------------------------\n")
                self.__lfh.write("File %s id  name : %s\n" % (fN, sdf.getDbCode("PDB")))
                self.__lfh.write("  polymers         = %d\n" % sdf.getEntityCount("polymer"))
                instD = sdf.getSequenceFeaturesFromAtomSite()
                for instId, rList in instD.items():
                    self.__lfh.write(" ------ Chain id= %s\n" % instId)
                    for rTup in rList:
                        self.__lfh.write(" %s %s %s %s %d\n" % (instId, rTup[0], rTup[1], rTup[2], rTup[3]))
        except Exception as _e:  # noqa: F841
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testGetEntityDetails(self):
        """Get entity lists  -"""
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % inspect.currentframe().f_code.co_name)
        self.__lfh.write(" ------------------------\n")
        try:
            for f in self.__examFileList:
                fN = os.path.join(self.__pathExamples, f)
                c0 = PdbxFileIo(verbose=self.__verbose, log=self.__lfh).getContainer(fN)
                sdf = ModelFileIo(dataContainer=c0, verbose=self.__verbose, log=self.__lfh)
                self.__lfh.write("------------------------------------------------------\n")
                self.__lfh.write("File %s id  idcode : %s\n" % (fN, sdf.getDbCode("PDB")))
                self.__lfh.write("  polymers         = %d\n" % sdf.getEntityCount("polymer"))
                #
                for eId in sdf.getPolymerEntityList("polymer"):
                    desc = sdf.getEntityDescription(entityId=eId)
                    name = sdf.getEntityName(entityId=eId)
                    self.__lfh.write("  entity id= %s  \n     +description=%s \n     +name=%s\n" % (eId, desc, name))

                self.__lfh.write("  non-polymers     = %d\n" % sdf.getEntityCount("non-polymer"))

        except Exception as _e:  # noqa: F841
            traceback.print_exc(file=self.__lfh)
            self.fail()


def suiteReadModelTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(PdbxIoUtilsTests("testAAAGetSeqFromAtomSite"))
    return suiteSelect


def suiteReadModelSourceTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(PdbxIoUtilsTests("testEntitySourceDetails"))
    return suiteSelect


def suiteReadEntityDetailsTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(PdbxIoUtilsTests("testGetEntityDetails"))
    return suiteSelect


if __name__ == "__main__":

    if True:  # pylint: disable=using-constant-test
        mySuite = suiteReadModelTests()
        unittest.TextTestRunner(verbosity=2).run(mySuite)

        mySuite = suiteReadModelSourceTests()
        unittest.TextTestRunner(verbosity=2).run(mySuite)

    if True:  # pylint: disable=using-constant-test
        mySuite = suiteReadEntityDetailsTests()
        unittest.TextTestRunner(verbosity=2).run(mySuite)
