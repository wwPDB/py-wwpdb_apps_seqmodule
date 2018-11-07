##
# File:    SequenceDataReviewTests.py
# Date:    16-Sep-2013  J. Westbrook
#
# Updates:
#
##
"""
Test cases and examples for access and update of sequence annotation data.

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"


import sys, unittest, traceback, shutil
import time, os, os.path

from wwpdb.apps.seqmodule.control.SequenceDataAssemble  import SequenceDataAssemble
from wwpdb.apps.seqmodule.io.SequenceDataStore          import SequenceDataStore
from wwpdb.apps.seqmodule.util.SequenceLabel            import SequenceLabel
from wwpdb.api.facade.ConfigInfo                        import ConfigInfo,getSiteId
from wwpdb.apps.seqmodule.webapp.SeqModWebRequest       import SeqModInputRequest
from wwpdb.apps.seqmodule.align.AlignmentStatistics     import AlignmentStatistics
from wwpdb.utils.rcsb.PathInfo               import PathInfo

class SequenceDataReviewTests(unittest.TestCase):
    def setUp(self):
        #
        self.__verbose=True
        self.__debug=False
        self.__lfh=sys.stdout

        #  nat src
        #self.__dsList= ['D_058993']
        self.__dsList= ['D_1100200023']
        #
        # input parameters -- 
        #
        # ref_id = auth_entity_1
        # selectids = [ auth_A_1_1_2,xyz_A_1_1_1,ref_A_1_100_1 ]
        self.mySetup()

    def mySetup(self,sessionId=None):
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

        self.__sessionId = self.__reqObj.getSessionId()   
        self.__sessionObj=self.__reqObj.newSessionObj()
        self.__sessionPath = self.__sessionObj.getPath()
        self.__seqDataCachePath=os.path.join(self.__reqObj.getSessionPath(),'SEQUENCE-DATA-CACHE')
        #
        #self.__pI=PathInfo(siteId=self.__siteId,sessionPath=self.__sessionPath,verbose=self.__verbose,log=self.__lfh)
        #self.__filePath=seqDataStatsPath=self.__pI.getSequenceStatsFilePath(self.__identifier,fileSource='session')            

        #

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
            for dsId in self.__dsList:
                self.__reqObj.setValue("identifier", dsId)
                sda=SequenceDataAssemble(reqObj=self.__reqObj,verbose=self.__verbose,log=self.__lfh)
                sda.doAssemble(fileSource='local-repository',cachePath=self.__seqDataCachePath)
                alstat=AlignmentStatistics(reqObj=self.__reqObj,maxRefAlign=self.__maxRefAlign,verbose=self.__verbose,log=self.__lfh)
                alstat.doUpdate()
                if self.__debug:
                    self.__testAccessMethods()
                self.__testGetAuthDetails(entityId='1')
                self.__testGetXyzDetails(entityId='1')
                self.__testGetRefDetails(refId='ref_A_1_100_1')
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime=time.clock()
        self.__lfh.write("\nCompleted %s %s at %s (%.2f seconds)\n" % (self.__class__.__name__,
                                                                     sys._getframe().f_code.co_name,
                                                                     time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                                     endTime-startTime))

    def __testGetRefDetails(self,refId): 
        """  SequenceDataStore() accessor tests -  seqId = <type>_<inst>_<part>_<altId>_<versionId>
        """
        startTime=time.clock()        
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__,
                                                       sys._getframe().f_code.co_name,
                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime())))

        try:
            self.__sds=SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            refAssignD=self.__sds.getReferenceAssignments()
            refAssignDepositorD=self.__sds.getDepositorReferenceAssignments()

            if (self.__verbose):
                self.__lfh.write("+SequenceDataReview.__getRefDetails() refId %s instance refAssignD          %r\n" % (refId,refAssignD.items()))
                self.__lfh.write("+SequenceDataReview.__getRefDetails() refId %s instance refAssignDepositorD %r\n" % (refId,refAssignDepositorD.items()))

            sLab=SequenceLabel()
            sLab.unpack(refId)
            seqType,seqId0,partId,altId,versionId=sLab.get()
            if (self.__verbose):
                self.__lfh.write("+SequenceDataReview.__getAuthDetails() refId %s instance id %s\n" % (refId,seqId0))
                self.__lfh.write("   Sequence id %5s partId %s altId %s version %d \n" % (seqId0,partId,altId,versionId))
            pfD=self.__sds.getFeature(seqId0,seqType=seqType,partId=partId,altId=altId,version=versionId)
            fObj=self.__sds.getFeatureObj(seqId0,seqType=seqType,partId=partId,altId=altId,version=versionId)
            fObj.printIt(log=self.__lfh)
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime=time.clock()
        self.__lfh.write("\nCompleted %s %s at %s (%.2f seconds)\n" % (self.__class__.__name__,
                                                                     sys._getframe().f_code.co_name,
                                                                     time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                                       endTime-startTime))

    def __testGetAuthDetails(self,entityId): 
        """  SequenceDataStore() accessor tests for 'auth' seqType and input entityId.
        """
        startTime=time.clock()        
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__,
                                                       sys._getframe().f_code.co_name,
                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime())))

        try:
            self.__sds=SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            seqIdList=self.__sds.getGroup(entityId)
            if (self.__verbose):
                self.__lfh.write("+SequenceDataReview.__getAuthDetails() entityId %s instance list %s\n" % (entityId,seqIdList))

            if len(seqIdList) < 1:
                #  fail here
                pass

            seqId0=seqIdList[0]
            partIdList=self.__sds.getPartIds(seqId0, dataType="sequence", seqType='auth')
            for partId in partIdList:
                vL=self.__sds.getVersionIds(seqId0,partId=partId, altId=1, dataType="sequence", seqType="auth")
                if len(vL)>0:
                    # grab features for latest version -
                    self.__lfh.write("   Sequence id %5s partId %s version %d \n" % (seqId0,partId,vL[0]))
                    pfD=self.__sds.getFeature(seqId0,seqType="auth",partId=partId,altId=1,version=vL[0])
                    fObj=self.__sds.getFeatureObj(seqId0,seqType='auth',partId=partId, altId=1,version=vL[0])
                    fObj.printIt(log=self.__lfh)
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime=time.clock()
        self.__lfh.write("\nCompleted %s %s at %s (%.2f seconds)\n" % (self.__class__.__name__,
                                                                     sys._getframe().f_code.co_name,
                                                                     time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                                       endTime-startTime))


    def __testGetXyzDetails(self,entityId): 
        """  SequenceDataStore() accessor tests for 'xyz' seqType and input entityId.
        """
        startTime=time.clock()        
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__,
                                                       sys._getframe().f_code.co_name,
                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime())))

        try:
            self.__sds=SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            seqIdList=self.__sds.getGroup(entityId)
            if (self.__verbose):
                self.__lfh.write("+SequenceDataReview.__getXyzDetails() entityId %s instance list %s\n" % (entityId,seqIdList))

            if len(seqIdList) < 1:
                #  fail here
                pass

            for seqId0 in seqIdList:
                partId=1
                altId=1
                vL=self.__sds.getVersionIds(seqId0,partId=partId, altId=altId, dataType="sequence", seqType="xyz")
                if len(vL)>0:
                    self.__lfh.write("   Sequence id %5s partId %s version %d \n" % (seqId0,partId,vL[0]))
                    pfD=self.__sds.getFeature(seqId0,seqType="xyz",partId=partId,altId=altId,version=vL[0])
                    fObj=self.__sds.getFeatureObj(seqId0,seqType='xyz',partId=partId, altId=altId,version=vL[0])
                    fObj.printIt(log=self.__lfh)
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime=time.clock()
        self.__lfh.write("\nCompleted %s %s at %s (%.2f seconds)\n" % (self.__class__.__name__,
                                                                     sys._getframe().f_code.co_name,
                                                                     time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                                       endTime-startTime))



    def __testAccessDetails(self): 
        """  SequenceDataStore() accessor tests - 
        """
        startTime=time.clock()        
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__,
                                                       sys._getframe().f_code.co_name,
                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime())))

        try:
            sda=SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)

            sTypes = sda.getSequenceTypes()
            for st in sTypes:
                self.__lfh.write("Data for sequence type: %5s\n" % st)                
                ids = sda.getIds(dataType="sequence", seqType=st)
                for id in ids:
                    verList=sda.getVersionIds(seqId=id, partId=1, altId=1, dataType="sequence", seqType=st)
                    self.__lfh.write(" Sequence id : %5s version length %d\n" % (id,len(verList)))                                    
                    for ver in verList:
                        sL=sda.getSequence(id,st,partId=1, altId=1,version=ver)
                        self.__lfh.write("   Sequence type %5s id %5s version %d length %5d\n" % (st,id,ver,len(sL)))
                        fObj=sda.getFeatureObj(id,st,partId=1, altId=1,version=ver)
                        fObj.printIt(log=self.__lfh)

                    altList=sda.getAlternativeIds(seqId=id, dataType="sequence", seqType=st)
                    self.__lfh.write(" Sequence id : %5s alternative list  length %d\n" % (id,len(altList)))                                       
                    for altId in altList[:5]:
                        sL=sda.getSequence(seqId=id,seqType=st,partId=1, altId=altId,version=1)
                        self.__lfh.write("   Sequence type %5s id %5s altId %s version %d length %5d\n" % (st,id,altId,ver,len(sL)))
                        fObj=sda.getFeatureObj(seqId=id,seqType=st,partId=1, altId=altId,version=1)
                        fObj.printIt(log=self.__lfh)

        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime=time.clock()
        self.__lfh.write("\nCompleted %s %s at %s (%.2f seconds)\n" % (self.__class__.__name__,
                                                                       sys._getframe().f_code.co_name,
                                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                                                                       endTime-startTime))


def suiteSearchAndAssembleFromArchiveTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(SequenceDataReviewTests("testAssembleFromRespositoryCache"))
    #suiteSelect.addTest(SequenceDataReviewTests("testAccessMethods"))
    return suiteSelect


if __name__ == '__main__':

    if (True):
        mySuite=suiteSearchAndAssembleFromArchiveTests()
        unittest.TextTestRunner(verbosity=2).run(mySuite)


