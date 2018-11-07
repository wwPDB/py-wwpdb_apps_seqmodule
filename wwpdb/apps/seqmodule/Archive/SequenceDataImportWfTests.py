##
# File:    SequenceDataImportWfTests.py
# Date:    25-Apr-2010
#
# Updates:
##
"""
Test cases for sequence data import. 

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"

import sys, unittest, traceback,shutil
import time, os, os.path
#
from wwpdb.apps.seqmodule.webapp.WebRequest       import SequenceInputRequest
from wwpdb.apps.seqmodule.io.SessionManager       import SessionManager
from wwpdb.apps.seqmodule.io.SequenceDataImport   import SequenceDataImportWf
from wwpdb.apps.seqmodule.io.PolymerLinkDepiction import PolymerLinkDepiction

class SequenceDataImportWfTests(unittest.TestCase):
    def setUp(self):
        #
        self.__verbose=True
        self.__lfh=sys.stderr
        #
        # Create a session object and session directories for test cases
        self.__reqObj=SequenceInputRequest(paramDict={},verbose=self.__verbose,log=self.__lfh)
        self.__reqObj.setValue("SessionPath","../sessions")
        self.__sobj=self.__reqObj.newSessionObj()
        self.__sessionPath=self.__sobj.getPath()
        #
        self.__depDataSetIdList=['D_055430','D_055453','D_056215','D_057171','D_057525','D_057584','D_057620','D_057630',
                                 'D_057750','D_057776','D_058195','D_058198','D_058417','D_1009416','D_101544','D_101653',
                                 'D_1040975','D_1043050','D_1043325','D_1043518']
        self.__wfInstanceId='W_000009'

        self.__depDataSeqIdExample='D_056215'

    def tearDown(self):
        pass
    

    def testImportArchiveOne(self):
        """ 
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            # 
            depDataSetId=self.__depDataSeqIdExample
            sdi=SequenceDataImportWf(reqObj=self.__reqObj,depDataSetId=depDataSetId,verbose=self.__verbose,log=self.__lfh)
            sdi.doImport()
            pL=sdi.getAtypicalPolymerLinkages()
            pld=PolymerLinkDepiction(upperBound=-100., lowerBound=100. ,verbose=self.__verbose,log=self.__lfh)
            oL=pld.buildPolymerLinkageTable(pL)
            self.__lfh.write("%s" % "\n".join(oL))            
        except:
            traceback.print_exc(file=sys.stderr)
            self.fail()


    def xtestImportWfInstanceOne(self):
        """ 
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            # 
            depDataSetId=self.__depDataSeqIdExample
            wfInstanceId=self.__wfInstanceId
            sdi=SequenceDataImportWf(reqObj=self.__reqObj,depDataSetId=depDataSetId,wfInstanceId=wfInstanceId,verbose=self.__verbose,log=self.__lfh)
            sdi.doImportWithCheck()
            pL=sdi.getAtypicalPolymerLinkages()
            pld=PolymerLinkDepiction(upperBound=-100., lowerBound=100. ,verbose=self.__verbose,log=self.__lfh)
            oL=pld.buildPolymerLinkageTable(pL)
            self.__lfh.write("%s" % oL.join("\n"))
            
            
        except:
            traceback.print_exc(file=sys.stderr)
            self.fail()

    def xtestImportWfList(self): 
        """ 
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            for depId in self.__depDataSetIdList:
                sdi=SequenceDataImportWf(reqObj=self.__reqObj,depDataSetId=depId,verbose=self.__verbose,log=self.__lfh)
                sdi.doImportWithCheck()
                sdi.printIt()
                self.__sobj.remakeSessionPath()
        except:
            traceback.print_exc(file=sys.stderr)
            self.fail()



def suite():
    suite = unittest.TestSuite()
    suite.addTest(SequenceDataImportWfTests("testImportArchiveOne"))
    return suite
    
def suiteAll():
    return unittest.makeSuite(SequenceDataImportWfTests,'test')

if __name__ == '__main__':
    unittest.main()

