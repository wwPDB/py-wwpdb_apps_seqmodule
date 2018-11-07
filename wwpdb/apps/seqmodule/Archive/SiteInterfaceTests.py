##
# File:    SiteInterfaceTests.py
# Date:    5-Feb-2010
#
# Updates:
# 20-Apr-2010 jdw Ported to module seqmodule.
##
"""
Test cases for site integration utilities which find and copy files
from the current data processing pipelines.

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"


import sys, unittest, traceback
import time, os, os.path

from wwpdb.apps.seqmodule.io.SiteInterface     import InterfaceSequenceRcsb, InterfaceArchiveRcsb
from wwpdb.apps.seqmodule.webapp.WebRequest    import SequenceInputRequest

class SiteInterfaceTests(unittest.TestCase):
    def setUp(self):
        self.__verbose=True
        self.__lfh=sys.stderr
        #
        # Create a session object and session directories for test cases
        self.__reqObj=SequenceInputRequest(paramDict={},verbose=self.__verbose,log=self.__lfh)
        self.__reqObj.setValue("SessionPath","../sessions")
        self.__sobj=self.__reqObj.newSessionObj()
        #
        # Load up some test data - 
        #
        self.__idList=['RCSB101544','rcsb056608', 'RCSB054845']
        

    def tearDown(self):
        pass
    

    def testArchiveRcsbList(self): 
        """ 
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            sI=InterfaceArchiveRcsb(self.__sobj,verbose=self.__verbose)
            for id in self.__idList:
                sI.setId(id)
                pth=sI.copyStructureFile("cif")
                self.__lfh.write("Archive file copied to path: %s\n" % pth)
                pth=sI.copyStructureFile("eps")
                self.__lfh.write("Encapsulated file copied to path: %s\n" % pth)                
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()


    def testSequenceRcsbList(self): 
        """ 
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            sI=InterfaceSequenceRcsb(self.__sobj,verbose=self.__verbose)
            for id in self.__idList:
                sI.setId(id)
                pth=sI.copyEntityFile(entityId='1')
                self.__lfh.write("Sequence file copied to path: %s\n" % pth)                
            
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()


def suite():
    return unittest.makeSuite(SiteInterfaceTests,'test')

if __name__ == '__main__':
    unittest.main()

