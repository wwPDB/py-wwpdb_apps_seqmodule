##
# File:    ReferenceSequenceDataUpdateTests.py
# Date:    18-Mar-2013
#
# Updates:
#
##
"""
Test cases for out of band search and updates to reference sequence data in the sequence data store.

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"


import sys, unittest, shutil, traceback
import time, os, os.path

from wwpdb.apps.seqmodule.control.SequenceDataAssemble           import SequenceDataAssemble
from wwpdb.apps.seqmodule.align.AlignmentStatistics              import AlignmentStatistics
from wwpdb.apps.seqmodule.control.ReferenceSequenceDataUpdate    import ReferenceSequenceDataUpdate

from wwpdb.utils.rcsb.FormatOut                         import FormatOut
from wwpdb.api.facade.ConfigInfo                        import ConfigInfo,getSiteId
from wwpdb.apps.seqmodule.webapp.SeqModWebRequest       import SeqModInputRequest
from wwpdb.utils.rcsb.PathInfo                          import PathInfo
from wwpdb.utils.rcsb.DataFileAdapter                   import DataFileAdapter

class ReferenceSequenceDataUpdateTests(unittest.TestCase):
    def setUp(self):
        #
        self.__verbose=True
        self.__lfh=sys.stderr
        self.__pathExamplesRel = "../tests"
        self.__pathExamples    = os.path.abspath(self.__pathExamplesRel)                
        #
        ##self.__exampleIdList    = ['1cbs','3rer','rcsb056751']
        self.__dsList             = ['D_000000']
        self.__exampleFileType    ='rcsb-mmcif'
        self.__exampleIdList      = ['rcsb056751']
        #self.__exampleIdList =['rcsb056215']
        #
        self.__myConfig(sessionId='seq-ref-update-tests')

    def __myConfig(self,sessionId='seq-ref-update-tests'):
        """  Simulate the web application environment for managing session storage of 
             temporaty data files.
        """
        self.__maxRefAlign=100
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
        #
        self.__seqDataCachePath=os.path.join(self.__reqObj.getSessionPath(),'SEQUENCE-DATA-CACHE')        

    def testSearchAndAssembleFromUpload(self): 
        """ Test search each entity sequence against appropriate reference sequence database
            service storing the matching sequences.  

            Using upload file source. rcsb-mmcif
        """
        startTime=time.clock()        
        self.__lfh.write("\n\n========================================================================================================\n")
        self.__lfh.write("Starting %s %s at %s\n" % (self.__class__.__name__,
                                                       sys._getframe().f_code.co_name,
                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime())))

        try:
            pI=PathInfo(siteId=self.__siteId,sessionPath=self.__sessionPath,verbose=self.__verbose,log=self.__lfh)
            for idCode in self.__exampleIdList:
                pdbxFilePath=pI.getModelPdbxFilePath(idCode,fileSource='session')
                inpFilePath=os.path.join(self.__pathExamples,idCode+".cif")
                self.__lfh.write("+testSearchAndAssembleFromUpload() Starting with id %s \n + input file %s \n   + pdbxfile %s\n" % 
                                 (idCode,inpFilePath,pdbxFilePath))
                dfa=DataFileAdapter(reqObj=self.__reqObj,verbose=self.__verbose,log=self.__lfh)
                ok=dfa.modelConvertToPdbx(filePath=inpFilePath,fileType=self.__exampleFileType,pdbxFilePath=pdbxFilePath)
                self.__reqObj.setValue("identifier", idCode)
                #
                if not os.access(inpFilePath,os.F_OK):
                    self.__lfh.write("+testSearchAndAssembleFromUpload() input file missing %s\n" % inpFilePath)
                    self.fail()
                    break
                if not os.access(pdbxFilePath,os.F_OK):
                    self.__lfh.write("+testSearchAndAssembleFromUpload() format conversion failed for %s\n" % idCode)
                    self.fail()
                    break
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

    def testSearchAndAssembleFromUpload(self): 
        """ Test search each entity sequence against appropriate reference sequence database
            service storing the matching sequences.  

            Using examples files in to simulate upload file source.
        """
        startTime=time.clock()        
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__,
                                                       sys._getframe().f_code.co_name,
                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime())))

        try:
            pI=PathInfo(siteId=self.__siteId,sessionPath=self.__sessionPath,verbose=self.__verbose,log=self.__lfh)
            for idCode in self.__exampleIdList:
                pdbxFilePath=pI.getModelPdbxFilePath(idCode,fileSource='session')
                inpFilePath=os.path.join(self.__pathExamples,idCode+'.cif')
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

    def testUpdateOperations(self): 
        """ Test update of reference sequence data store 

            Use example/test case setup from repository cache - 
        """
        startTime=time.clock()        
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__,
                                                       sys._getframe().f_code.co_name,
                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime())))
        try:
            identifier = self.__exampleIdList[0]
            self.__reqObj.setValue("identifier", identifier)
            eId=str(1)
            rsdU=ReferenceSequenceDataUpdate(reqObj=self.__reqObj,verbose=self.__verbose,log=self.__lfh)
            rsdU.dump()
            rsdU.doUpdate(entityId=eId)
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
    #suiteSelect.addTest(ReferenceSequenceDataUpdateTests("testSearchAndAssembleFromArchive"))
    #suiteSelect.addTest(ReferenceSequenceDataUpdateTests("testSearchAndAssembleFromUpload"))
    suiteSelect.addTest(ReferenceSequenceDataUpdateTests("testAssembleFromRespositoryCache"))
    return suiteSelect


def suiteReferenceSequenceDataUpdateTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(ReferenceSequenceDataUpdateTests("testUpdateOperations"))
    return suiteSelect



if __name__ == '__main__':
    if (True):
        mySuite=suiteSearchAndAssembleTests()
        unittest.TextTestRunner(verbosity=2).run(mySuite)

    if (True):
        mySuite=suiteReferenceSequenceDataUpdateTests()
        unittest.TextTestRunner(verbosity=2).run(mySuite)

