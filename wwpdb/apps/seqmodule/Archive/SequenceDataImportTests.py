##
# File:    SequenceDataImportTests.py
# Date:    5-Feb-2010
#
# Updates:
#   12-Mar-2010 jdw Revised example handling.
#   20-Apr-2010 jdw Ported to module seqmodule.
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
from wwpdb.apps.seqmodule.webapp.WebRequest      import SequenceInputRequest
from wwpdb.apps.seqmodule.io.SessionManager      import SessionManager
from wwpdb.apps.seqmodule.io.SequenceDataImport  import SequenceDataImportRcsb, SequenceDataImportExample

class SequenceDataImportTests(unittest.TestCase):
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
        
        ## ID details - 
        #
        # These directories mimic the rcsb repository and sequence data file system organization.
        #self.__reqObj.setValue("RcsbDataPath","../examples/rcsb-data")
        #self.__reqObj.setValue("RcsbReferenceSequencePath","../examples/rcsb-sequence")
        #self.__rcsbIdExample='RCSB101544'
        self.__rcsbIdExample='RCSB058417'
        #
        # 04-13 test examples pdb_extract
        # rcsb101653 rcsb058198 rcsb058195
        #
        # real path to rcsb repository
        #self.__reqObj.setValue("RcsbDataPath","/annotation")
        #self.__reqObj.setValue("RcsbReferenceSequencePath","/www-rcsb/supertool/blast/rcsb")                
        self.__rcsbIdRepositoryList=['RCSB160027','RCSB054845','RCSB056183','RCSB056213','RCSB056218',
                                     'RCSB056042', 'RCSB056211', 'RCSB056214',  'RCSB056294', 'RCSB029053', 'RCSB029052'] 
        self.__rcsbId='RCSB042832'

        # Some other test data ids.
        #
        # Multi-domain match
        #
        # ribosome
        #self.__rcsbId= 'RCSB056215'
        # 
        #self.__rcsbId='rcsb056901'
        #self.__rcsbId='rcsb057630'
        # 


    def tearDown(self):
        pass
    
    def testImportExampleOne(self):
        """ 
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            self.__reqObj.setValue("RcsbDataPath","../examples/rcsb-data")
            self.__reqObj.setValue("RcsbReferenceSequencePath","../examples/rcsb-sequence")
            rcsbIdExample=self.__rcsbIdExample
            sdi=SequenceDataImportRcsb(reqObj=self.__reqObj,rcsbId=rcsbIdExample,verbose=self.__verbose,log=self.__lfh)
            sdi.doImport()
        except:
            traceback.print_exc(file=sys.stdout)
            self.fail()

    def testImportRepositoryOne(self):
        """ 
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            # real path
            self.__reqObj.setValue("RcsbDataPath","/annotation")
            self.__reqObj.setValue("RcsbReferenceSequencePath","/www-rcsb/supertool/blast/rcsb")                
            rcsbId='RCSB101544'
            sdi=SequenceDataImportRcsb(reqObj=self.__reqObj,rcsbId=rcsbId,verbose=self.__verbose,log=self.__lfh)
            sdi.doImport()
        except:
            traceback.print_exc(file=sys.stdout)
            self.fail()

    def testImportRepositoryList(self): 
        """ 
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            self.__reqObj.setValue("RcsbDataPath","/annotation")
            self.__reqObj.setValue("RcsbReferenceSequencePath","/www-rcsb/supertool/blast/rcsb")                
            for rcsbId in self.__rcsbIdRepositoryList:
                sdi=SequenceDataImportRcsb(reqObj=self.__reqObj,rcsbId=rcsbId,verbose=self.__verbose,log=self.__lfh,fetchPdbFile=True)
                sdi.doImport()
                sdi.printIt()
                self.__sobj.remakeSessionPath()
        except:
            traceback.print_exc(file=sys.stdout)
            self.fail()

    def testUploadOne(self): 
        """ 
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            rcsbId=self.__rcsbIdExample
            # emulate the web file upload by copying data file to the session area.
            self.__reqObj.setValue("RcsbDataPath","../examples/rcsb-data")
            self.__reqObj.setValue("RcsbReferenceSequencePath","../examples/rcsb-sequence")
            #
            fType='cif'
            idLc=str(rcsbId).lower()
            fName=idLc+'.cif'
            src=os.path.join(self.__reqObj.getValue("RcsbDataPath"),'test',idLc,fName)
            dst=os.path.join(self.__sessionPath,fName)
            shutil.copyfile(src,dst)
            #
            self.__reqObj.setValue("UploadFileType",fType)
            self.__reqObj.setValue("UploadFileName",fName)
            sdi=SequenceDataImportRcsb(reqObj=self.__reqObj,rcsbId=rcsbId,fileSource='upload',verbose=self.__verbose,log=self.__lfh,fetchPdbFile=True)            
            sdi.doImport()
            sdi.printIt()
            
        except:
            traceback.print_exc(file=sys.stdout)
            self.fail()


    def testImportEmbeddedExample1(self): 
        """Test older embedded example data.
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            sdi=SequenceDataImportExample(reqObj=self.__reqObj,verbose=self.__verbose,log=self.__lfh)
            sdi.doImport()
            
        except:
            traceback.print_exc(file=sys.stdout)
            self.fail()



def suite():
    suite = unittest.TestSuite()
    suite.addTest(SequenceDataImportTests("testImportExampleOne"))
    return suite
    
def suiteAll():
    return unittest.makeSuite(SequenceDataImportTests,'test')

if __name__ == '__main__':
    unittest.main()

