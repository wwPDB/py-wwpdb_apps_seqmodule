##
# File:    PairwiseAlignTests.py
# Date:    11-Jan-2009
#
# Updates:
# 20-Apr-2010 jdw Ported to module seqmodule.
##
"""
Test cases for pairwise sequence alignment wrapper class.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.06"

import sys
import traceback
import unittest

from wwpdb.utils.align.alignlib import PairwiseAlign

from wwpdb.apps.seqmodule.util.SequenceExamples import SequenceExamples


class PairwiseAlignTests(unittest.TestCase):
    def setUp(self):
        self.__verbose = True
        sE = SequenceExamples()
        self.seqRef3L = sE.getRefSequence3List('A')
        self.seqAuth3L = sE.getAuthSequenceList('A')
        #
        # Test sequence with random insertions and deletions
        #
        self.sTests = {}
        for tt in ['T1', 'T2', 'T3', 'T4', 'T5']:
            self.sTests[tt] = sE.getAuthSequenceListTest('A')
        #

    def tearDown(self):
        pass

    def testAlign(self):
        """ Run internal alignment test embedded in the class -  This is a basic santity check.
        """
        sys.stdout.write("\n------------------------ ")
        sys.stdout.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        sys.stdout.write(" -------------------------\n")
        try:
            pA = PairwiseAlign()
            pA.testExample()
        except:
            traceback.print_exc(file=sys.stdout)
            self.fail()

    def testAlign2(self):
        """ Run author vs reference sequence alignment returning a copy of the alignment
            via getAlignment() -         
        """
        sys.stdout.write("\n------------------------ ")
        sys.stdout.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        sys.stdout.write(" -------------------------\n")
        try:
            pA = PairwiseAlign()
            pA.setVerbose(self.__verbose)
            pA.setReferenceSequence(self.seqRef3L, "REFA")
            pA.addTestSequence(self.seqAuth3L, "A")
            pA.doAlign()
            pA.prAlignmentConflicts("A")
            sys.stdout.write("Length of reference sequence = %d\n" % len(self.seqRef3L))
            sys.stdout.write("Length of    author sequence = %d\n" % len(self.seqAuth3L))
            myAlign = pA.getAlignment("A")
            sys.stdout.write("Length   of alignment     = %d\n" % len(myAlign))
            ii = 0
            for myPr in myAlign:
                if myPr[0] != myPr[1]:
                    sys.stdout.write("Py - conflict at alignment position %d  -  %s - %s\n" % (ii, myPr[0], myPr[1]))
                ii += 1
        except:
            traceback.print_exc(file=sys.stdout)
            self.fail()

    def testAlign3(self):
        """  Consensus alignment for author and reference sequences.
             Returning the name of any sequence that is not part of the consensus.
        """
        sys.stdout.write("\n------------------------ ")
        sys.stdout.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        sys.stdout.write(" -------------------------\n")
        try:
            pA = PairwiseAlign()
            pA.setVerbose(self.__verbose)
            pA.setReferenceSequence(self.seqRef3L, "REFA")
            pA.addTestSequence(self.seqAuth3L, "A")
            #
            myFails = pA.doAlignConsensus()
            sys.stdout.write("Failed sequences = %d\n" % len(myFails))
            pA.prAlignmentConflicts("A")
            sys.stdout.write("Length reference sequence = %d\n" % len(self.seqRef3L))
            sys.stdout.write("Length    author sequence = %d\n" % len(self.seqAuth3L))
            myAlign = pA.getAlignment("A")
            sys.stdout.write("Length   of alignment     = %d\n" % len(myAlign))
            ii = 0
            for myPr in myAlign:
                if myPr[0] != myPr[1]:
                    sys.stdout.write("Py - conflict position %8d  -  %3s - %3s\n" % (ii, myPr[0], myPr[1]))
                ii += 1
        except:
            traceback.print_exc(file=sys.stdout)
            self.fail()

    def testAlign4(self):
        """  Consensus alignment for reference and 5 test sequences

        """
        sys.stdout.write("\n------------------------ ")
        sys.stdout.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        sys.stdout.write(" -------------------------\n")
        try:
            pA = PairwiseAlign()
            pA.setVerbose(self.__verbose)
            pA.setReferenceSequence(self.seqRef3L, "REFA")
            for k, v in self.sTests.items():
                sys.stdout.write("Added test sequence %10s len = %8d\n" % (k, len(v)))
                pA.addTestSequence(v, k)
            #
            sys.stdout.flush()
            myFails = pA.doAlignConsensus()
            sys.stdout.write("Failed sequences = %d\n" % len(myFails))
            sys.stdout.write("Length reference sequence = %d\n" % len(self.seqRef3L))
            sys.stdout.flush()
            pA.prAlignmentFull()
        except:
            traceback.print_exc(file=sys.stdout)
            self.fail()

    def testAlign5(self):
        """  Consensus alignment for reference and 5 test sequences

        """
        sys.stdout.write("\n------------------------ ")
        sys.stdout.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        sys.stdout.write(" -------------------------\n")
        try:
            pA = PairwiseAlign()
            pA.setVerbose(self.__verbose)
            pA.setReferenceSequence(self.seqRef3L, "REFA")
            sys.stdout.write("Length reference sequence = %d\n" % len(self.seqRef3L))
            #
            for k, v in self.sTests.items():
                sys.stdout.write("Added test sequence %10s len = %8d\n" % (k, len(v)))
                pA.addTestSequence(v, k)
            #
            sys.stdout.flush()
            myFails = pA.doAlignConsensus()
            sys.stdout.write("Failed sequences = %d\n" % len(myFails))

            sys.stdout.flush()
            ofh = sys.stdout
            ofh.write("This example includes %d test sequences.\n" % len(self.sTests))
            for k, v in self.sTests.items():
                ofh.write("Conflict list for case %s\n" % k)
                pA.prAlignmentConflicts(k)
            ofh.write("Alignment Details \n")
            pA.prAlignmentFull()
            ofh.write("Done\n")
    
        except:
            traceback.print_exc(file=sys.stdout)
            self.fail()


def suite():
    return unittest.makeSuite(PairwiseAlignTests, 'test')


if __name__ == '__main__':
    unittest.main()
