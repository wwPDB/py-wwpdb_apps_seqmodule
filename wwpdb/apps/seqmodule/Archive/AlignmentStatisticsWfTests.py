##
# File:    SequenceDataImportWfTests.py
# Date:    25-Apr-2010
#
# Updates:
##
"""
Test cases for workflow alignment assessment and sequence selection export .

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
from wwpdb.apps.seqmodule.io.SequenceDataImport  import SequenceDataImportWf
from wwpdb.apps.seqmodule.align.AlignmentStatistics import AlignmentStatistics
#
from wwpdb.apps.seqmodule.io.SequenceDataExport  import SequenceDataExportRcsb

from wwpdb.wwpdb.utils.wf.DataReference  import DataFileReference
from wwpdb.utils.config.ConfigInfo import ConfigInfo
from wwpdb.utils.rcsb.RcsbDpUtil import RcsbDpUtil

class SequenceDataImportWfTests(unittest.TestCase):
    def setUp(self):
        #
        self.__verbose=True
        self.__lfh=sys.stderr
        self.__cleanUp=False
        #
        self.__depDataSetIdList=['D_055453','D_056215','D_057171','D_057525','D_057584','D_057620','D_057630',
                                 'D_057750','D_057776','D_058195','D_058198','D_058417','D_1009416','D_101544','D_101653',
                                 'D_1040975','D_1043050','D_1043325','D_1043518']        
        #
        self.__wfInstanceId='W_000001'
        self.__depDataSeqIdExample='D_057171'
        self.__fileSource='wf-archive'
        #
        # Create a session object and session directories for test cases
        
        self.__reqObj=SequenceInputRequest(paramDict={},verbose=self.__verbose,log=self.__lfh)

        #
        # We are now using this request object and session path and this may change to be id specific.
        #
        self.__reqObj.setValue("SessionPath","../sessions")
        self.__sobj=self.__reqObj.newSessionObj()
        self.__sessionPath=self.__sobj.getPath()
        #

    def tearDown(self):
        pass
    

    def testFromArchiveOne(self):
        """ 
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:

            ##  This provide all information to define source data -- 
            wfInstanceId=self.__wfInstanceId
            depDataSetId=self.__depDataSeqIdExample
            fileSource='archive'
            ##
            ##
            dfRef=DataFileReference()
            dfRef.setDepositionDataSetId(depDataSetId)
            dfRef.setStorageType(fileSource)
            dfRef.setContentTypeAndFormat('model','pdbx')
            dfRef.setVersionId('latest')
            ##
            ##
            if (dfRef.isReferenceValid()):                  
                dP=dfRef.getDirPathReference()
                fP=dfRef.getFilePathReference()
                if (self.__verbose):                
                    self.__lfh.write("+AlignmentStatisticsWfTests() PDBx model directory path: %s\n" % dP)
                    self.__lfh.write("+AlignmentStatisticsWftests() PDBx model file      path: %s\n" % fP)
            
            modelPath=fP
            sdi=SequenceDataImportWf(reqObj=self.__reqObj,depDataSetId=depDataSetId,
                                     wfInstanceId=wfInstanceId,fileSource=fileSource,
                                     verbose=self.__verbose,log=self.__lfh)
            sdi.doImport()
            #sdi.printIt()
            alstat=AlignmentStatistics(sessionObj=self.__sobj,verbose=self.__verbose,log=self.__lfh)
            alstat.doUpdate()
            sdsPath=sdi.getSequenceDataStorePath()
            #
            self.__lfh.write("+AlignmentStatisticsWfTests() Sequence data store file path: %s\n" % sdsPath)            
            #
            dfRef=DataFileReference()
            dfRef.setDepositionDataSetId(depDataSetId)
            dfRef.setStorageType(fileSource)
            dfRef.setContentTypeAndFormat('seqed-summary-stats','pic')
            dfRef.setVersionId('latest')
            ##
            ##
            if (dfRef.isReferenceValid()):                  
                dP=dfRef.getDirPathReference()
                fP=dfRef.getFilePathReference()
                if (self.__verbose):                
                    self.__lfh.write("+AlignmentStatisticsWfTests() Summary stats directory path: %s\n" % dP)
                    self.__lfh.write("+AlignmentStatisticsWftests() Summary status file     path: %s\n" % fP)
            

            sdx0=SequenceDataExportRcsb(reqObj=self.__reqObj,verbose=self.__verbose,log=self.__lfh)
            selectList=sdx0.makeDefaultSelection()
            sdx=SequenceDataExportRcsb(reqObj=self.__reqObj,exportList=selectList, verbose=self.__verbose,log=self.__lfh)
            sdx.doExportTrial()
            expPath=sdx.getFilePath()
            ##
            self.__lfh.write("+AlignmentStatisticsWfTests() Export mapping file  path: %s\n" % expPath)
            #
            # Prepare file for maxit
            #
            oModelPath=os.path.join(self.__sessionPath,'model-prep.cif')
            ofh=open(oModelPath,'w')
            ifh=open(modelPath,'r')
            ofh.writelines(ifh.readlines())
            ifh.close()            
            #
            ifh=open(expPath,'r')
            lines=ifh.readlines()
            ofh.writelines(lines[2:])
            #
            ofh.close()

            #
            updModelPath=os.path.join(self.__sessionPath,'model-upd.cif')            
            #
            dirPath =self.__sessionPath
            #
            cI=ConfigInfo()
            rcsbPath=cI.get("SITE_RCSB_APPS_PATH")
            dp=RcsbDpUtil(rcsbPath=rcsbPath,tmpPath=dirPath)
            dp.imp(oModelPath)
            dp.op("cif-seqed2cif-pdbx")
            dp.exp(updModelPath)
            if (self.__cleanUp): dp.cleanup()                        
            if (self.__verbose):
                self.__lfh.write("+AlignmentStatisticsWfTests() - Model prep    file path: %s\n" % oModelPath)
                self.__lfh.write("+AlignmentStatisticsWfTests() - Updated model file path: %s\n" % updModelPath)                                

        except:
            traceback.print_exc(file=sys.stderr)
            self.fail()


    def testFromArchiveAll(self):
        """ 
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            for depId in self.__depDataSetIdList:
                self.__reqObj=SequenceInputRequest(paramDict={},verbose=self.__verbose,log=self.__lfh)

                self.__reqObj.setValue("SessionPath","../sessions")
                self.__sobj=self.__reqObj.newSessionObj()
                self.__sessionPath=self.__sobj.getPath()
                #
                
                wfInstanceId=self.__wfInstanceId
                depDataSetId=depId
                fileSource='archive'
                ##
                ##
                dfRef=DataFileReference()
                dfRef.setDepositionDataSetId(depDataSetId)
                dfRef.setStorageType(fileSource)
                dfRef.setContentTypeAndFormat('model','pdbx')
                dfRef.setVersionId('latest')
                ##
                ##
                if (dfRef.isReferenceValid()):                  
                    dP=dfRef.getDirPathReference()
                    fP=dfRef.getFilePathReference()
                    if (self.__verbose):                
                        self.__lfh.write("+AlignmentStatisticsWfTests() PDBx model directory path: %s\n" % dP)
                        self.__lfh.write("+AlignmentStatisticsWftests() PDBx model file      path: %s\n" % fP)

                modelPath=fP
                sdi=SequenceDataImportWf(reqObj=self.__reqObj,depDataSetId=depDataSetId,
                                         wfInstanceId=wfInstanceId,fileSource=fileSource,
                                         verbose=self.__verbose,log=self.__lfh)
                sdi.doImport()
                #sdi.printIt()
                alstat=AlignmentStatistics(sessionObj=self.__sobj,verbose=self.__verbose,log=self.__lfh)
                alstat.doUpdate()
                sdsPath=sdi.getSequenceDataStorePath()
                #
                self.__lfh.write("+AlignmentStatisticsWfTests() Sequence data store file path: %s\n" % sdsPath)            
                #
                dfRef=DataFileReference()
                dfRef.setDepositionDataSetId(depDataSetId)
                dfRef.setStorageType(fileSource)
                dfRef.setContentTypeAndFormat('seqed-summary-stats','pic')
                dfRef.setVersionId('latest')
                ##
                ##
                if (dfRef.isReferenceValid()):                  
                    dP=dfRef.getDirPathReference()
                    fP=dfRef.getFilePathReference()
                    if (self.__verbose):                
                        self.__lfh.write("+AlignmentStatisticsWfTests() Summary stats directory path: %s\n" % dP)
                        self.__lfh.write("+AlignmentStatisticsWftests() Summary status file     path: %s\n" % fP)


                sdx0=SequenceDataExportRcsb(reqObj=self.__reqObj,verbose=self.__verbose,log=self.__lfh)
                selectList=sdx0.makeDefaultSelection()
                sdx=SequenceDataExportRcsb(reqObj=self.__reqObj,exportList=selectList, verbose=self.__verbose,log=self.__lfh)
                sdx.doExportTrial()
                expPath=sdx.getFilePath()
                ##
                self.__lfh.write("+AlignmentStatisticsWfTests() Export mapping file  path: %s\n" % expPath)
                #
                # Prepare file for maxit
                #
                oModelPath=os.path.join(self.__sessionPath,'model-prep.cif')
                ofh=open(oModelPath,'w')
                ifh=open(modelPath,'r')
                ofh.writelines(ifh.readlines())
                ifh.close()            
                #
                ifh=open(expPath,'r')
                lines=ifh.readlines()
                ofh.writelines(lines[2:])
                #
                ofh.close()

                #
                updModelPath=os.path.join(self.__sessionPath,'model-upd.cif')            
                #
                dirPath =self.__sessionPath
                #
                cI=ConfigInfo()
                rcsbPath=cI.get("SITE_RCSB_APPS_PATH")
                dp=RcsbDpUtil(rcsbPath=rcsbPath,tmpPath=dirPath)
                dp.imp(oModelPath)
                dp.op("cif-seqed2cif-pdbx")
                dp.exp(updModelPath)
                if (self.__cleanUp): dp.cleanup()                        
                if (self.__verbose):
                    self.__lfh.write("+AlignmentStatisticsWfTests() - Model prep    file path: %s\n" % oModelPath)
                    self.__lfh.write("+AlignmentStatisticsWfTests() - Updated model file path: %s\n" % updModelPath)                                

        except:
            traceback.print_exc(file=sys.stderr)
            self.fail()


    def xtestImportWfList(self): 
        """ 
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function  %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            for depId in self.__depDataSetIdList:
                sdi=SequenceDataImportWf(reqObj=self.__reqObj,depDataSetId=depId,verbose=self.__verbose,log=self.__lfh)
                sdi.doImport()
                sdi.printIt()
                self.__sobj.remakeSessionPath()
        except:
            traceback.print_exc(file=sys.stderr)
            self.fail()



def suite():
    suite = unittest.TestSuite()
    suite.addTest(AlignmentStatisticsWfTests("testFromArchiveOne"))
    return suite
    
def suiteAll():
    return unittest.makeSuite(AlignmentStatisticsWfTests,'test')

if __name__ == '__main__':
    unittest.main()

