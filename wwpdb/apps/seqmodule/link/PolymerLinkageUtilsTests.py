##
# File:    PolymerLinkageUtilsTests.py
# Date:    21-Feb-2013
#
# Updates:
# 27-Feb-2013  jdw update site configuration discovery
##
"""
Test cases for calculating polymer linkage distances and reading calculated linkage results.

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"


import sys, unittest, traceback
import time, os, os.path


from wwpdb.apps.seqmodule.link.PolymerLinkageUtils  import PolymerLinkageUtils
from wwpdb.apps.seqmodule.link.PolymerLinkageDepict import PolymerLinkageDepict
from wwpdb.io.misc.FormatOut  import FormatOut

from wwpdb.utils.config.ConfigInfo                        import ConfigInfo,getSiteId
from wwpdb.apps.seqmodule.webapp.SeqModWebRequest       import SeqModInputRequest
from wwpdb.io.locator.PathInfo                          import PathInfo


class PolymerLinkageUtilsTests(unittest.TestCase):
    def setUp(self):
        #
        self.__verbose=True
        self.__lfh=sys.stdout
        self.__pathExamplesRel = "../tests"
        self.__pathExamples    = os.path.abspath(self.__pathExamplesRel)                
        #
        self.__examFileList    = ['1cbs.cif','3rer.cif']
        self.__setup()

    def __setup(self,sessionId="testing"):
        """  Simulate the web application environment for managing session storage of 
             temporaty data files.
        """
        self.__siteId=getSiteId(defaultSiteId="WWPDB_DEPOLY_TEST")
        self.__cI=ConfigInfo(self.__siteId)
        self.__topPath=self.__cI.get('SITE_WEB_APPS_TOP_PATH')
        self.__topSessionPath=self.__cI.get('SITE_WEB_APPS_TOP_SESSIONS_PATH')
        #
        self.__reqObj=SeqModInputRequest({},verbose=self.__verbose,log=self.__lfh)
        self.__templatePath = os.path.join(self.__topPath,"htdocs","seqmodule")
        self.__reqObj.setValue("TopSessionPath", self.__topSessionPath)
        self.__reqObj.setValue("TemplatePath",self.__templatePath)
        self.__reqObj.setValue("TopPath", self.__topPath)
        self.__reqObj.setValue("WWPDB_SITE_ID", self.__siteId)
        os.environ["WWPDB_SITE_ID"]=self.__siteId
        
        if sessionId is not None:
            self.__reqObj.setValue("sessionid", sessionId)

        self.__sessionId = self.__reqObj.getSessionId()   
        self.__sObj=self.__reqObj.newSessionObj()

    def tearDown(self):
        pass
    
    
    def testGetPolymerLinkages(self): 
        """ Calculate and read polymer linkage data files. 
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")        
        try:
            for f in self.__examFileList:
                fN=os.path.join(self.__pathExamples,f)
                (fId,fExt)=os.path.splitext(f)

                pI=PathInfo(siteId="WWPDB_DEPLOY_TEST",sessionPath=".",verbose=self.__verbose,log=self.__lfh)
                plFile=pI.getPolyLinkFilePath(fId,fileSource='session')
                self.__lfh.write("Link dist  (PDBx):   %s\n" % plFile)

                plu=PolymerLinkageUtils(reqObj=self.__reqObj,verbose=self.__verbose, log=self.__lfh)
                linkList=plu.calc(fN,plFile)
                
                fOut=FormatOut()
                fOut.autoFormat("Link list", linkList,3,3)
                fOut.writeStream(self.__lfh)
                fOut.clear()
                #
                # conver linkList to sequence list - 
                #
                linkInstD={}
                for lD in linkList:
                    seqId=lD['auth_seq_id_1']
                    compId=lD['auth_comp_id_1']
                    authAsymId=lD['auth_asym_id_1']
                    dist=lD['dist']
                    if float(dist) > 1.38:
                        linkInstD[(authAsymId,compId,seqId)]= (compId,seqId,'link=%s'%dist,seqId)
                #
                self.__lfh.write('link dictionary %r\n' % linkInstD.items())

                
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()




def suiteGetPolyLinkTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(PolymerLinkageUtilsTests("testGetPolymerLinkages"))
    return suiteSelect


if __name__ == '__main__':
    if (True):
        mySuite=suiteGetPolyLinkTests()
        unittest.TextTestRunner(verbosity=2).run(mySuite)

