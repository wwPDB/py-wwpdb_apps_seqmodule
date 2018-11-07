##
# File:    ReferenceSequenceUtilsLocalTests.py
# Date:    17-Apr-2013
#
# Updates:
#   20-Apr-2013  jdw add ribosome example -- 
#
##
"""
Test cases for performing local reference sequence database searches and processing search results.

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"


import sys, unittest, traceback
import time, os, os.path

from wwpdb.apps.seqmodule.io.ReferenceSequenceUtils import ReferenceSequenceUtils
from wwpdb.apps.seqmodule.io.ModelSequenceUtils     import ModelSequenceUtils
from wwpdb.apps.seqmodule.io.PdbxIoUtils            import PdbxFileIo,ReferenceSequenceIo
from wwpdb.utils.rcsb.FormatOut                     import FormatOut
from wwpdb.utils.rcsb.PathInfo                      import PathInfo
from wwpdb.api.facade.ConfigInfo                    import ConfigInfo,getSiteId

class ReferenceSequenceUtilsLocalTests(unittest.TestCase):
    def setUp(self):
        #
        self.__verbose=True
        self.__debug=True
        self.__lfh=sys.stdout
        self.__siteId=getSiteId(defaultSiteId="WWPDB_DEPLOY_TEST")
        self.__pathExamplesRel = "../tests"
        self.__pathExamples    = os.path.abspath(self.__pathExamplesRel)                
        #
        # isoform example
        self.__examFileList    = ['rcsb052365.cif']
        #  ribosome test case
        #self.__examFileList    = ['rcsb056215.cif']
        #
    
    def testSearchEntitySequences(self): 
        """ Test reference sequence database search for each entity sequence extracted from the input 
            model coordinate file.  The appropriate reference sequence database service is selected.
            Matching sequences are extracted from service output files and stored along with key
            reference sequence features in PDBx data files.
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")        
        try:
            pI=PathInfo(siteId=self.__siteId,sessionPath='.',verbose=self.__verbose,log=self.__lfh)
            for f in self.__examFileList:
                fN=os.path.join(self.__pathExamples,f)
                c0 = PdbxFileIo(verbose=self.__verbose,log=self.__lfh).getContainer(fN)
                msu=ModelSequenceUtils(dataContainer=c0,verbose=self.__verbose,log=self.__lfh)
                entityD=msu.getEntitySequenceDetails()
                if self.__debug:
                    fOut=FormatOut()
                    fOut.autoFormat("Entity dictionary", entityD,3,3)
                    fOut.writeStream(self.__lfh)
                #
                rsu=ReferenceSequenceUtils(verbose=self.__verbose,log=self.__lfh)
                rsio=ReferenceSequenceIo(verbose=self.__verbose,log=self.__lfh)
                for eId,eD in entityD.items():
                    entryIdL=str(eD['ENTRY_ID']).lower()
                    if len(eD['SEQ_ENTITY_1']) > 10:
                        mR=rsu.searchEntities(entityD=eD,saveBlast=True,filePrefix=entryIdL,localSearch=True)
                        if len(mR)>0:
                            fn=pI.getReferenceSequenceFilePath(entryIdL,entityId=eId,fileSource='session')
                            rsio.writeMatchResults(eD,outFilePath=fn,matchResults=mR)
                            rc0 = PdbxFileIo(verbose=self.__verbose,log=self.__lfh).getContainer(fn)
                            #
                            rsin=ReferenceSequenceIo(dataContainer=rc0,verbose=self.__verbose,log=self.__lfh)
                            rL=rsin.readMatchResults()
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()


def suiteSeqSearchTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ReferenceSequenceUtilsLocalTests("testSearchEntitySequences"))
    return suiteSelect



if __name__ == '__main__':
    if (True):
        mySuite=suiteSeqSearchTests()
        unittest.TextTestRunner(verbosity=2).run(mySuite)


