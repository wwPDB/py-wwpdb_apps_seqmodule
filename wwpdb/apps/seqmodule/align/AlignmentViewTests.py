##
# File:    AlignmentViewTests.py
# Date:    26-Feb-2013
#
# Updates:
# 04-Mar-2013  jdw added test case for multipart entity.
#
##
"""
Test cases for computing alignment statistics for an existing sequence data store.  

The sequence data store is created by the  SequenceDataAssembly() class.  

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"

import time, os, os.path
import sys, unittest, shutil, traceback
#
from wwpdb.apps.seqmodule.align.AlignmentView          import AlignmentView
from wwpdb.apps.seqmodule.align.AlignmentViewDepiction import AlignmentViewDepiction
from wwpdb.apps.seqmodule.align.AlignmentStatistics    import AlignmentStatistics
from wwpdb.apps.seqmodule.control.SequenceDataAssemble_v2 import SequenceDataAssemble
#
from wwpdb.io.misc.FormatOut                     import FormatOut
from wwpdb.io.locator.PathInfo                      import PathInfo
from wwpdb.utils.config.ConfigInfo                    import ConfigInfo,getSiteId
from wwpdb.apps.seqmodule.webapp.SeqModWebRequest   import SeqModInputRequest

class AlignmentViewTests(unittest.TestCase):
    def setUp(self):
        #
        self.__verbose=True
        self.__lfh=sys.stderr
        self.__pathExamplesRel = "../tests"
        self.__pathExamples    = os.path.abspath(self.__pathExamplesRel)                
        #
        #self.__exampleFileList    = ['1cbs.cif','3rer.cif','rcsb056751.cif']
        self.__exampleFileList    = ['rcsb056751.cif']
        self.__dsList=['D_000000']
        self.__exampleIdList      = ['rcsb056751']
        #
        self.__maxRefAlign=100
        #
        self.__setup()

    def __setup(self,sessionId='seq-alignview-tests'):
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

        self.__sessionId   = self.__reqObj.getSessionId()   
        self.__sessionObj  = self.__reqObj.newSessionObj()
        self.__sessionPath = self.__sessionObj.getPath()
        self.__seqDataCachePath=os.path.join(self.__reqObj.getSessionPath(),'SEQUENCE-DATA-CACHE')        

    def testSearchAndAssembleFromUpload(self): 
        """ Test search each entity sequence against appropriate reference sequence database
            service storing the matching sequences.  

            Using upload file source.
        """
        startTime=time.clock()        
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__,
                                                       sys._getframe().f_code.co_name,
                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime())))

        try:
            pI=PathInfo(siteId=self.__siteId,sessionPath=self.__sessionPath,verbose=self.__verbose,log=self.__lfh)
            for f in self.__exampleFileList:
                (idCode,fExt)=os.path.splitext(f)
                pdbxFilePath=pI.getModelPdbxFilePath(idCode,fileSource='session')
                inpFilePath=os.path.join(self.__pathExamples,f)
                shutil.copyfile(inpFilePath,pdbxFilePath)
                self.__reqObj.setValue("identifier", idCode)
                #
                sda=SequenceDataAssemble(reqObj=self.__reqObj,verbose=self.__verbose,log=self.__lfh)
                sda.doAssemble(fileSource='local-upload')
                alstat=AlignmentStatistics(reqObj=self.__reqObj,maxRefAlign=self.__maxRefAlign,verbose=self.__verbose,log=self.__lfh)
                alstat.doUpdate()
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime=time.clock()
        self.__lfh.write("\nCompleted %s %s at %s (%.2f seconds)\n" % (self.__class__.__name__,
                                                                     sys._getframe().f_code.co_name,
                                                                     time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                                     endTime-startTime))

    def testSearchAndAssembleFromArchive(self): 
        """ Test search each entity sequence against appropriate reference sequence database
            service storing the matching sequences.  

            Using archive file source.
        """
        startTime=time.clock()        
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__,
                                                       sys._getframe().f_code.co_name,
                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime())))

        try:
            for dsId in self.__dsList:
                self.__reqObj.setValue("identifier", dsId)
                sda=SequenceDataAssemble(reqObj=self.__reqObj,verbose=self.__verbose,log=self.__lfh)
                sda.doAssemble(fileSource='archive')
                alstat=AlignmentStatistics(reqObj=self.__reqObj,maxRefAlign=self.__maxRefAlign,verbose=self.__verbose,log=self.__lfh)
                alstat.doUpdate()
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime=time.clock()
        self.__lfh.write("\nCompleted %s %s at %s (%.2f seconds)\n" % (self.__class__.__name__,
                                                                     sys._getframe().f_code.co_name,
                                                                     time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                                     endTime-startTime))


    def testAssembleFromRespositoryCache(self): 
        """ Test search each entity sequence against appropriate reference sequence database
            service storing the matching sequences.  

            Using repository cache with 'session' file source.
        """
        startTime=time.clock()        
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__,
                                                       sys._getframe().f_code.co_name,
                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime())))

        try:
            for dsId in self.__exampleIdList:
                self.__reqObj.setValue("identifier", dsId)
                sda=SequenceDataAssemble(reqObj=self.__reqObj,verbose=self.__verbose,log=self.__lfh)
                sda.doAssemble(fileSource='local-repository',cachePath=self.__seqDataCachePath)
                alstat=AlignmentStatistics(reqObj=self.__reqObj,maxRefAlign=self.__maxRefAlign,verbose=self.__verbose,log=self.__lfh)
                alstat.doUpdate()
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime=time.clock()
        self.__lfh.write("\nCompleted %s %s at %s (%.2f seconds)\n" % (self.__class__.__name__,
                                                                     sys._getframe().f_code.co_name,
                                                                     time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                                     endTime-startTime))



    def testSimpleAlignView(self): 
        """ Create simple alignment view for a single triple or sequence identifiers
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            for f in self.__exampleFileList:
                (idCode,fExt)=os.path.splitext(f)
                idString="['auth_A_1_1_1,xyz_A_1_1_1,ref_A_1_1_1']"
                self.__reqObj.setValue("operation","load")                        
                self.__reqObj.setValue("alignids",idString)         
                self.__reqObj.setValue("identifier",idCode)
                alignViewOrder = self.__reqObj.getAlignmentOrdering()       
                #
                aV=AlignmentView(reqObj=self.__reqObj, verbose=self.__verbose,log=self.__lfh)
                alignSeqList=aV.loadAlign()
            
                #
                avd=AlignmentViewDepiction(self.__verbose,self.__lfh)
                oL=avd.renderAlignment(alignSeqList=alignSeqList,type='original',viewOrderCode=alignViewOrder)
                cL=avd.renderConflictTable(alignSeqList=alignSeqList,type='original')                    

                self.__lfh.write("\n\nRendered alignment length %d:\n" % len(oL))
                self.__lfh.write("%s\n" % "\n".join(oL))
                #
                self.__lfh.write("\n\nRendered conflict report length %d:\n" % len(cL))
                self.__lfh.write("%s\n" % "\n".join(cL))            

        except:
            traceback.print_exc(file=sys.stderr)
            self.fail()

    def testAlignViewFromRespositoryCache(self): 
        """ Test construction of a summary view from a sequence data store containing pairwise alignment stats.

            Use example/test case setup from repository cache - 
        """
        startTime=time.clock()        
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__,
                                                       sys._getframe().f_code.co_name,
                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime())))
        try:
            idString="['auth_A_1_1_1,ref_A_1_1_1,ref_A_2_1_1']"
            #idString="['auth_A_1_1_1,xyz_A_1_1_1,ref_A_1_1_1']"
            self.__reqObj.setValue("identifier", self.__exampleIdList[0])
            #
            self.__reqObj.setValue("operation","load")                        
            self.__reqObj.setValue("alignids",idString)  

            alignViewOrder = self.__reqObj.getAlignmentOrdering()       
            #
            aV=AlignmentView(reqObj=self.__reqObj, verbose=self.__verbose,log=self.__lfh)
            alignSeqList=aV.loadAlign()
            
            #
            avd=AlignmentViewDepiction(self.__verbose,self.__lfh)
            oL=avd.renderAlignment(alignSeqList=alignSeqList,type='original',viewOrderCode=alignViewOrder)
            cL=avd.renderConflictTable(alignSeqList=alignSeqList,type='original')                    

            self.__lfh.write("\n\nRendered alignment length %d:\n" % len(oL))
            self.__lfh.write("%s\n" % "\n".join(oL))
            #
            self.__lfh.write("\n\nRendered conflict report length %d:\n" % len(cL))
            self.__lfh.write("%s\n" % "\n".join(cL))            

        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime=time.clock()
        self.__lfh.write("\nCompleted %s %s at %s (%.2f seconds)\n" % (self.__class__.__name__,
                                                                     sys._getframe().f_code.co_name,
                                                                     time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                                     endTime-startTime))



def suiteSearchAndAssembleTests():
    suiteSelect = unittest.TestSuite()
    #suiteSelect.addTest(AlignmentViewTests("testSearchAndAssembleFromUpload"))
    #suiteSelect.addTest(AlignmentViewTests("testSearchAndAssembleFromArchive"))
    suiteSelect.addTest(AlignmentViewTests("testAssembleFromRespositoryCache"))
    return suiteSelect


def suiteAlignViewSimpleTests():
    suiteSelect = unittest.TestSuite()
    #suiteSelect.addTest(AlignmentViewTests("testSimpleAlignView"))
    suiteSelect.addTest(AlignmentViewTests("testAlignViewFromRespositoryCache"))

    return suiteSelect

if __name__ == '__main__':
    if (True):
        mySuite=suiteSearchAndAssembleTests()
        unittest.TextTestRunner(verbosity=2).run(mySuite)
    if (True):    
        mySuite=suiteAlignViewSimpleTests()
        unittest.TextTestRunner(verbosity=2).run(mySuite)
