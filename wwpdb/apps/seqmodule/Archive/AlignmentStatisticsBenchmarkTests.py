##
# File:    AlignmentStatisticsBenchmarkTests.py
# Author:  jdw
# Date:    1-Feb-2010
#
# Updates:
# 12-Mar-2010  Revised example handling.
# 20-Apr-2010 jdw - ported to module seqmodule.
##

"""
Test cases for the update of alignment statistics.  These tests first
load data to the sequence data store and then perform an update of the
alignment statistics in the feature dictionaries for coordinate and
reference sequences.

Note - these test cases require the RCSB environment -

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"

import sys, unittest, traceback
import time, os, os.path

from wwpdb.apps.seqmodule.webapp.WebRequest         import SequenceInputRequest
from wwpdb.apps.seqmodule.io.SequenceDataImport     import SequenceDataImportRcsb
from wwpdb.apps.seqmodule.align.AlignmentStatistics import AlignmentStatistics

class AlignmentStatisticsTests(unittest.TestCase):
    def setUp(self):
        # Create a session object
        #
        self.__verbose=True
        self.__lfh = sys.stderr
        #
        # Create a session object and session directories for test cases
        self.__reqObj=SequenceInputRequest(paramDict={},verbose=self.__verbose,log=self.__lfh)
        self.__reqObj.setValue("SessionPath","../sessions")
        self.__sobj=self.__reqObj.newSessionObj()
        self.__sessionPath=self.__sobj.getPath()
        #
        self.__rcsbIdExample="RCSB101544"
        self.__idList=['rcsb005868','rcsb030617','rcsb037577','rcsb042832','rcsb043474','rcsb050503','rcsb052365','rcsb053095','rcsb053659',
                       'rcsb055430','rcsb055453','rcsb056215','rcsb057171','rcsb057525','rcsb057584','rcsb057620','rcsb057630','rcsb057776',
                       'rcsb100925','rcsb101544']
    def tearDown(self):
        pass
    
    def testStatsRcsbExample(self): 
        """ 
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            self.__reqObj.setValue("RcsbDataPath","../examples/rcsb-data")
            self.__reqObj.setValue("RcsbReferenceSequencePath","../examples/rcsb-sequence")
            
            sdi=SequenceDataImportRcsb(reqObj=self.__reqObj,rcsbId=self.__rcsbIdExample,fileSource='rcsb-repository',
                                       verbose=self.__verbose,log=self.__lfh,fetchPdbFile=False)
            sdi.doImport()
            #sdi.printIt()
            alstat=AlignmentStatistics(sessionObj=self.__sobj,verbose=self.__verbose,log=self.__lfh)
            alstat.doUpdate()
            
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testStatsExampleList(self): 
        """ 
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            for id in self.__idList:
                self.__reqObj.setValue("RcsbDataPath","../examples/rcsb-data")
                self.__reqObj.setValue("RcsbReferenceSequencePath","../examples/rcsb-sequence")                
                sdi=SequenceDataImportRcsb(reqObj=self.__reqObj,rcsbId=id,fileSource='rcsb-repository',
                                           verbose=self.__verbose,log=self.__lfh,fetchPdbFile=False)
                sdi.doImport()
                sdi.printIt()
                alstat=AlignmentStatistics(sessionObj=self.__sobj,verbose=self.__verbose,log=self.__lfh)
                alstat.doUpdate()
            
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()


def suite():
    return unittest.makeSuite(AlignmentStatisticsTests,'test')

if __name__ == '__main__':
    unittest.main()

