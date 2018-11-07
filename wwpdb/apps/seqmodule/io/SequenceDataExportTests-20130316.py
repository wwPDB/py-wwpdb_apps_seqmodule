##
# File:    SequenceDataExportTests.py
# Date:    5-Feb-2010
# Updates:
#   12-Mar-2010 jdw Revised example handling.
#   20-Apr-2010 jdw Ported to module seqmodule.
#    2-May-2010 jdw Add SequenceSelection()
#                   Add export file path.
##
"""
Test cases for sequence data import. 

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"

import sys, unittest, traceback
import time, os, os.path
#
from wwpdb.apps.seqmodule.webapp.WebRequest      import SequenceInputRequest
from wwpdb.apps.seqmodule.io.SequenceDataImport  import SequenceDataImportRcsb
from wwpdb.apps.seqmodule.io.SequenceSelection   import SequenceSelection
from wwpdb.apps.seqmodule.io.SequenceDataExport  import SequenceDataExport
from wwpdb.apps.seqmodule.align.AlignmentStatistics import AlignmentStatistics

class SequenceDataExportTests(unittest.TestCase):
    def setUp(self):
        # Create a session object
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
        # Local example data -
        #
        #self.__rcsbIdExample='RCSB101544'
        # Load up some test data - 
        #
        # Multi-domain match test example
        #self.__rcsbId='RCSB101544'
        # ribosome example
        self.__rcsbIdExample= 'RCSB056215'
        #self.__rcsbId= 'RCSB042832'        
        self.__rcsbIdRepositoryList=['RCSB160027','RCSB054845','RCSB056183','RCSB056213','RCSB056218',
                                     'RCSB056042', 'RCSB056211', 'RCSB056214',  'RCSB056294', 'RCSB029053', 'RCSB029052']
        #
        self.__rcsbIdExampleList=['rcsb005868','rcsb030617','rcsb037577','rcsb042832','rcsb043474','rcsb050503','rcsb052365','rcsb053095','rcsb053659',
                                  'rcsb055430','rcsb055453','rcsb056215','rcsb057171','rcsb057525','rcsb057584','rcsb057620','rcsb057630','rcsb057776',
                                  'rcsb100925','rcsb101544']

    def tearDown(self):
        pass
    

    def testExportExampleOne(self): 
        """ 
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            self.__reqObj.setValue("RcsbDataPath","../examples/rcsb-data")
            self.__reqObj.setValue("RcsbReferenceSequencePath","../examples/rcsb-sequence")
            sdi=SequenceDataImportRcsb(reqObj=self.__reqObj,rcsbId=self.__rcsbIdExample,verbose=self.__verbose,log=self.__lfh,fetchPdbFile=False)
            sdi.doImport()
            sdi.printIt()
            alstat=AlignmentStatistics(sessionObj=self.__sobj,verbose=self.__verbose,log=self.__lfh)
            alstat.doUpdate()
            #
            sSel=SequenceSelection(reqObj=self.__reqObj,verbose=self.__verbose,log=self.__lfh)            
            selectList=sSel.makeDefaultSelection()            
            ePath=os.path.join(self.__sessionPath,'export-file.cif')                
            sdx=SequenceDataExport(reqObj=self.__reqObj,exportList=selectList, exportFilePath=ePath, verbose=self.__verbose,log=self.__lfh)
            sdx.doExportTrial()
            
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def xtestExportExampleList(self): 
        """ 
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            self.__reqObj.setValue("RcsbDataPath","../examples/rcsb-data")
            self.__reqObj.setValue("RcsbReferenceSequencePath","../examples/rcsb-sequence")
            for rcsbId in self.__rcsbIdExampleList:
                self.__lfh.write("++++Starting Export Test for Example ID %s\n" % rcsbId)
                sdi=SequenceDataImportRcsb(reqObj=self.__reqObj,rcsbId=rcsbId,verbose=self.__verbose,log=self.__lfh,fetchPdbFile=False)
                sdi.doImport()
                #sdi.printIt()
                #
                sSel=SequenceSelection(reqObj=self.__reqObj,verbose=self.__verbose,log=self.__lfh)            
                selectList=sSel.makeDefaultSelection()
                ePath=os.path.join(self.__sessionPath,'export-file.cif')                
                sdx=SequenceDataExport(reqObj=self.__reqObj,exportList=selectList, exportFilePath=ePath,verbose=self.__verbose,log=self.__lfh)
                sdx.doExportTrial()
                self.__sobj.remakeSessionPath()                
            
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()


    def xtestExportRepositoryOne(self): 
        """ 
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            self.__reqObj.setValue("RcsbDataPath","/annotation")
            self.__reqObj.setValue("RcsbReferenceSequencePath","/www-rcsb/supertool/blast/rcsb")                
            sdi=SequenceDataImportRcsb(reqObj=self.__reqObj,rcsbId=self.__rcsbId,verbose=self.__verbose,log=self.__lfh,fetchPdbFile=False)
            sdi.doImport()
            sdi.printIt()

            sSel=SequenceSelection(reqObj=self.__reqObj,verbose=self.__verbose,log=self.__lfh)            
            selectList=sSel.makeDefaultSelection()
            
            ePath=os.path.join(self.__sessionPath,'export-file.cif')
            sdx=SequenceDataExport(reqObj=self.__reqObj,exportList=selectList, exportFilePath=ePath,verbose=self.__verbose,log=self.__lfh)
            sdx.doExportTrial()
            
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()


def suite():
    return unittest.makeSuite(SequenceDataExportTests,'test')

if __name__ == '__main__':
    unittest.main()

