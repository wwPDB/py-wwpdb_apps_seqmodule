##
# File:    SummaryViewWfTests.py
# Date:    05-May-2010
#
# Updates:
##
"""
Test cases for summary view data generation.

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"

import sys, unittest, traceback
import time, os, os.path

from wwpdb.apps.seqmodule.control.SummaryView          import SummaryView
from wwpdb.apps.seqmodule.control.SummaryViewDepiction import SummaryViewDepiction
#from wwpdb.apps.seqmodule.control.DataImportView      import DataImportView
from wwpdb.apps.seqmodule.control.DataImporter         import DataImporter
from wwpdb.apps.seqmodule.webapp.WebRequest            import SequenceInputRequest

class SummaryViewWfTests(unittest.TestCase):
    def setUp(self):
        self.__verbose=True
        self.__lfh=sys.stderr
        #
        #
        self.__depDataSetIdList=['D_055430','D_055453','D_056215','D_057171','D_057525','D_057584','D_057620','D_057630',
                                 'D_057750','D_057776','D_058195','D_058198','D_058417','D_1009416','D_101544','D_101653',
                                 'D_1040975','D_1043050','D_1043325','D_1043518']
        self.__wfInstanceId='W_000009'

        self.__depDataSeqIdExample='D_056215'
        self.__fileSource='archive'
        
        # Create a session object and session directories for test cases
        self.__reqObj=SequenceInputRequest(paramDict={},verbose=self.__verbose,log=self.__lfh)
        self.__reqObj.setValue("SessionPath","../sessions")
        self.__sobj=self.__reqObj.newSessionObj()
        self.__sessionPath=self.__sobj.getPath()

    def tearDown(self):
        pass
    
    def testViewSummaryWfArchiveOne(self): 
        """ Verify data in summary data object for example test cases - 
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            self.__sobj.remakeSessionPath()
            self.__reqObj.setValue("identifier", self.__depDataSeqIdExample)
            self.__reqObj.setValue("instance",   self.__wfInstanceId)
            self.__reqObj.setValue("filesource", self.__fileSource)                        
            self.__reqObj.printIt(self.__lfh)
            #
            dI=DataImporter(reqObj=self.__reqObj, fileSource='archive', verbose=self.__verbose, log=self.__lfh)
            ok=dI.loadData()
            #
            sV=SummaryView(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            sumObj=sV.loadSummary(operation='load')
            sV.formatTable(sumObj,self.__lfh)
            sVD=SummaryViewDepiction(verbose=self.__verbose,log=self.__lfh)
            oL=sVD.buildSummaryView(sumObj)
            self.__lfh.write("%s" % '\n'.join(oL))            
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()


    def xtestViewSummaryRcsbExampleOne(self): 
        """ Verify data in summary data object for example test cases - 
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            self.__sobj.remakeSessionPath()
            self.__reqObj.setValue("RcsbDataPath","../examples/rcsb-data")
            self.__reqObj.setValue("RcsbReferenceSequencePath","../examples/rcsb-sequence")
            self.__rcsbIdExample='RCSB101544'
            self.__reqObj.setValue("identifier",self.__rcsbIdExample)                        
            #self.__reqObj.setValue("operation","loadrcsb")
            self.__reqObj.printIt(self.__lfh)
            #sV=DataImportView(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            #dObj=sV.loadData()                        
            #
            dI=DataImporter(reqObj=self.__reqObj, fileSource='rcsb-repository', verbose=self.__verbose, log=self.__lfh)
            ok=dI.loadData()
            #
            sV=SummaryView(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            sumObj=sV.loadSummary(operation='load')
            sV.formatTable(sumObj,self.__lfh)
            sVD=SummaryViewDepiction(verbose=self.__verbose,log=self.__lfh)
            oL=sVD.buildSummaryView(sumObj)
            self.__lfh.write("%s" % '\n'.join(oL))                        
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()


def suite():
    return unittest.makeSuite(SummaryViewWfTests,'test')

if __name__ == '__main__':
    unittest.main()

