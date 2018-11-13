##
# File:    BlastProcessTests.py
# Author:  j. westbrook
# Date:    2-Jan-2013
# Version: 0.001
#
# Update:

##
"""
Test cases from funning reference sequence database blast searches and processing the search results.

"""

import sys, unittest, os, os.path, traceback

from wwpdb.utils.config.ConfigInfo          import ConfigInfo,getSiteId
from wwpdb.utils.seqdb_v2.BlastProcess    import BlastProcess


class BlastProcessTests(unittest.TestCase):
    def setUp(self):
        self.__verbose=True
        self.__lfh=sys.stderr        
        # Pick up site information from the environment or failover to the development site id.
        self.__siteId=getSiteId(defaultSiteId='WWPDB_DEPLOY_TEST')
        self.__lfh.write("\nTesting with site environment for:  %s\n" % self.__siteId)
        #
        cI=ConfigInfo(self.__siteId)
        self.__testFilePath='./data'
        self.__testFileCif='4ec0.cif'
        self.__testFileFragmentsCif='rcsb056751.cif'        
        self.__taxonomyDataFile='nodes.dmp'
        
            
    def tearDown(self):
        pass

    def testGetPolymerEntityDetails(self): 
        """
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__, sys._getframe().f_code.co_name))
        try:
            for fn in [self.__testFileCif,self.__testFileFragmentsCif]:
                cifFilePath=os.path.join(self.__testFilePath,fn)
                taxonomyFilePath=os.path.join(self.__testFilePath,self.__taxonomyDataFile)
                bp= BlastProcess(cifFilePath=cifFilePath, taxonomyFilePath=taxonomyFilePath,verbose=self.__verbose,log=self.__lfh)
                ped=bp.getPolymerEntityDetails()
                for (entityId,polyDetails) in ped.items():
                    self.__lfh.write("\n----Entry %s  Entity Id = %s -------  \n" % (fn,entityId))
                    self.__lfh.write("one-letter-code = %s\n" % polyDetails['seq'])
                    self.__lfh.write("length          = %d\n" % len(polyDetails['seq'])) 
                    self.__lfh.write("type            = %s\n" % polyDetails['type'])
                    self.__lfh.write("fragments       = %r\n" % polyDetails['fragments'])
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testGetPolymerEntityDetailsFragments(self): 
        """
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__, sys._getframe().f_code.co_name))
        try:
            for fn in [self.__testFileCif,self.__testFileFragmentsCif]:            
                cifFilePath=os.path.join(self.__testFilePath,fn)
                taxonomyFilePath=os.path.join(self.__testFilePath,self.__taxonomyDataFile)
                bp= BlastProcess(cifFilePath=cifFilePath, taxonomyFilePath=taxonomyFilePath,verbose=self.__verbose,log=self.__lfh)
                ped=bp.getPolymerEntityDetails()
                for (entityId,polyDetails) in ped.items():
                    self.__lfh.write("\n----Entry %s  Entity Id = %s -------  \n" % (fn,entityId))                
                    self.__lfh.write("one-letter-code = %s\n" % polyDetails['seq'])
                    self.__lfh.write("length          = %d\n" % len(polyDetails['seq'])) 
                    self.__lfh.write("type            = %s\n" % polyDetails['type'])
                    self.__lfh.write("fragments       = %r\n" % polyDetails['fragments'])
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()


    def testPolymerSearch1(self): 
        """
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__, sys._getframe().f_code.co_name))
        try:
            for fn in [self.__testFileCif,self.__testFileFragmentsCif]:            
                cifFilePath=os.path.join(self.__testFilePath,fn)            
                taxonomyFilePath=os.path.join(self.__testFilePath,self.__taxonomyDataFile)
                bp= BlastProcess(cifFilePath=cifFilePath, taxonomyFilePath=taxonomyFilePath,verbose=self.__verbose,log=self.__lfh)
                ped=bp.getPolymerEntityDetails()
                for (entityId,polyDetails) in ped.items():            
                    resultlist = bp.Run(entityId=entityId)
                    self.__lfh.write("\n----Entry %s  Entity Id = %s -------  \n" % (fn,entityId))                                    
                    for ii,rL in enumerate(resultlist):
                        for r in rL:
                            for k in r.keys():
                                self.__lfh.write("%d: %s:  %s\n" % (ii,k,r[k]))
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testPolymerSearchAndStore(self): 
        """
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__, sys._getframe().f_code.co_name))
        try:
            for fn in [self.__testFileCif,self.__testFileFragmentsCif]:
                entryId,fExt=os.path.splitext(fn)
                cifFilePath=os.path.join(self.__testFilePath,fn)            
                taxonomyFilePath=os.path.join(self.__testFilePath,self.__taxonomyDataFile)
                bp= BlastProcess(cifFilePath=cifFilePath, taxonomyFilePath=taxonomyFilePath,verbose=self.__verbose,log=self.__lfh)
                bp.saveBlastResults(blastPath='.',blastFileNamePrefix=entryId)
                ped=bp.getPolymerEntityDetails()
                for (entityId,polyDetails) in ped.items():
                    ofn=entryId+"_seqdb-match_P"+str(entityId)+".cif"
                    ok = bp.RunAndSave(entityId=entityId,fName=ofn)
                    self.__lfh.write("\n----Entry %s  Entity Id = %s ofn = %s status = %r -------  \n" % (fn,entityId,ofn,ok))                                        


        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()


def suiteSearchTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(BlastProcessTests("testGetPolymerEntityDetails"))
    #suiteSelect.addTest(BlastProcessTests("testPolymerSearch1"))
    suiteSelect.addTest(BlastProcessTests("testPolymerSearchAndStore"))
    #
    return suiteSelect

    
if __name__ == '__main__':
    # Run all tests -- 
    # unittest.main()
    #
    mySuite=suiteSearchTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
    #
