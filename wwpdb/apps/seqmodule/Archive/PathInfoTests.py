##
# File:    PathInfoTests.py
# Date:    26-Feb-2013
#
# Updates:
##
"""
Test cases for creating standard file names for sequence resources and data files.

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"


import sys, unittest, traceback
import time, os, os.path


from wwpdb.apps.seqmodule.util.PathInfo import PathInfo

class PathInfoTests(unittest.TestCase):
    def setUp(self):
        #
        self.__verbose=True
        self.__lfh=sys.stdout
        
    def tearDown(self):
        pass
    
    def testGetStandardPaths(self): 
        """ Test getting standadard file names within session paths.
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")        
        try:
            # fileSource, id, eId, wfInst
            tests=[('session', "1abc", '1', None, 'latest'), 
                   ('session', 'e1', '1', None, 'latest'),
                   ('archive', 'D_000111', '1', None, '1'), 
                   ('archive', 'D_000111', '1', 'W_010','1'),
                   ('wf-instance', 'D_000111', '1', 'W_010','1')]
            dataSetId="1abc"
            eId='1'
            fileSourceList=['session']
            wfInst=None

            for (fs, dataSetId, eId, wfInst,vId) in tests: 
                self.__lfh.write("\n\n----------------------------------\n")
                self.__lfh.write("File source %s dataSetId %s entityId %s wfInst %s version %s\n" % (fs,dataSetId,eId,wfInst,vId))

                pI=PathInfo(siteId="WWPDB_DEPLOY_TEST",sessionPath=".",verbose=self.__verbose,log=self.__lfh)
                #
                fp=pI.getModelPdbxFilePath(dataSetId,wfInstanceId=wfInst, fileSource=fs, versionId=vId)
                self.__lfh.write("Model path (PDBx):   %s\n" % fp)

                fp=pI.getModelPdbFilePath(dataSetId,wfInstanceId=wfInst,fileSource=fs, versionId=vId)
                self.__lfh.write("Model path (PDB):    %s\n" % fp)

                fp=pI.getPolyLinkFilePath(dataSetId,wfInstanceId=wfInst,fileSource=fs, versionId=vId)
                self.__lfh.write("Link dist  (PDBx):   %s\n" % fp)

                fp=pI.getSequenceStatsFilePath(dataSetId,wfInstanceId=wfInst,fileSource=fs, versionId=vId)
                self.__lfh.write("Sequence stats (PIC):   %s\n" % fp)
                
                fp=pI.getReferenceSequenceFilePath(dataSetId,entityId=eId,wfInstanceId=wfInst,fileSource=fs, versionId=vId)
                self.__lfh.write("Reference match entity %s (PDBx):   %s\n" % (eId,fp))
                
                fp=pI.getSequenceAssignmentFilePath(dataSetId,wfInstanceId=wfInst,fileSource=fs, versionId=vId)
                self.__lfh.write("Sequence assignment (PDBx):   %s\n" % fp)

                if fs in ['session','wf-session']:
                    for ePart in ['1','2','3']:
                        fp=pI.getBlastMatchFilePath(entityId=eId,entityPart=ePart,wfInstanceId=wfInst,fileSource=fs, versionId=vId)
                        self.__lfh.write("Blast match entity %s  (xml):   %s\n" % (eId,fp))

                        
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()


def suiteStandardPathTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(PathInfoTests("testGetStandardPaths"))
    return suiteSelect


if __name__ == '__main__':
    if (True):
        mySuite=suiteStandardPathTests()
        unittest.TextTestRunner(verbosity=2).run(mySuite)

