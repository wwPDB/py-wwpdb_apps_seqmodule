##
# File:    Tom.py
# Date:    09-Apr-2013
#
# Updates:
#
##
"""
Examples for Tom

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"


import sys, unittest, traceback
import time, os, os.path

from wwpdb.apps.seqmodule.io.PdbxIoUtils import ModelFileIo,PdbxFileIo

class PdbxIoUtilsTests(unittest.TestCase):
    def setUp(self):
        #
        self.__verbose=True
        self.__lfh=sys.stdout
        ## TOM -- Put your path and test files here -- 
        self.__pathExamplesRel = "../tests"
        self.__pathExamples    = os.path.abspath(self.__pathExamplesRel)                
        self.__examFileList    = ['1cbs.cif']

    def tearDown(self):
        pass
    

    def testGetEntityCounts(self): 
        """Get entity counts  - 
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" ------------------------\n")        
        try:
            for f in self.__examFileList:
                fN=os.path.join(self.__pathExamples,f)
                c0 = PdbxFileIo(verbose=self.__verbose,log=self.__lfh).getContainer(fN)
                sdf = ModelFileIo(dataContainer=c0,verbose=self.__verbose,log=self.__lfh)
                self.__lfh.write("------------------------------------------------------\n")                
                self.__lfh.write("  polymers         = %d\n" % sdf.getEntityCount('polymer'))
                self.__lfh.write("  non-polymers     = %d\n" % sdf.getEntityCount('non-polymer'))
                self.__lfh.write("  (L) polypeptide polymers = %d\n" % sdf.getPolyPeptideLEntityCount())
                self.__lfh.write("  (D) polypeptide polymers = %d\n" % sdf.getPolyPeptideDEntityCount())
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testGetEntityLists(self): 
        """Get entity lists  - 
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" ------------------------\n")        
        try:
            for f in self.__examFileList:
                fN=os.path.join(self.__pathExamples,f)
                c0 = PdbxFileIo(verbose=self.__verbose,log=self.__lfh).getContainer(fN)
                sdf = ModelFileIo(dataContainer=c0,verbose=self.__verbose,log=self.__lfh)
                self.__lfh.write("------------------------------------------------------\n")
                self.__lfh.write("File %s id  name : %s\n" % (fN,sdf.getDbCode('PDB') ))                
                self.__lfh.write("  polymers         = %d\n" % sdf.getEntityCount('polymer'))
                #
                for eId in sdf.getPolymerEntityList('polymer'):
                    self.__lfh.write("  entity id= %s\n" % eId)
                self.__lfh.write("  non-polymers     = %d\n" % sdf.getEntityCount('non-polymer'))
                for eId in sdf.getPolymerEntityList('non-polymer'):
                    self.__lfh.write("  entity id= %s\n" % eId)                

        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()


    def testCategoryExists(self): 
        """Get category list and test for category existence -- 
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" ------------------------\n")        
        try:
            for f in self.__examFileList:
                fN=os.path.join(self.__pathExamples,f)
                c0 = PdbxFileIo(verbose=self.__verbose,log=self.__lfh).getContainer(fN)
                catList = c0.getObjNameList()
                self.__lfh.write(" Category list  = %r\n" % catList)                                
                #
                catNameList=['atom_site','pdbx_chem_comp_instance_depositor_info','symmetry']
                for catName in catNameList:
                    if (c0.exists(catName)):
                        self.__lfh.write("  Found category = %s\n" % catName)                
        
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()



def suiteReadEntityDetailsTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(PdbxIoUtilsTests("testGetEntityCounts"))
    suiteSelect.addTest(PdbxIoUtilsTests("testGetEntityLists"))
    suiteSelect.addTest(PdbxIoUtilsTests("testCategoryExists"))
    return suiteSelect


if __name__ == '__main__':

    if (True):
        mySuite=suiteReadEntityDetailsTests()
        unittest.TextTestRunner(verbosity=2).run(mySuite)
