##
# File:    RcsbPath.py
# Date:    4-Jun-2010
#
##
"""
Test cases for site integration utilities which find files from the current data processing pipelines.

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"


import sys, unittest, traceback
import time, os, os.path

from wwpdb.apps.seqmodule.io.RcsbPath          import RcsbPath

class RcsbPathTests(unittest.TestCase):
    def setUp(self):
        self.__verbose=True
        self.__lfh=sys.stderr
        #
        #
        self.__idList=['RCSB101544','rcsb056608', 'RCSB054845']
        

    def tearDown(self):
        pass
    

    def testGetPath(self): 
        """ 
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            sI=RcsbPath(verbose=self.__verbose)
            for id in self.__idList:
                sI.setId(id)
                for format in ['cif','pdb','eps','sf']:
                    pth=sI.getFilePath(format)
                    self.__lfh.write("Rcsb %s file path: %s\n" % (format,pth))

        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()


def suite():
    return unittest.makeSuite(RcsbPathTests,'test')

if __name__ == '__main__':
    unittest.main()

