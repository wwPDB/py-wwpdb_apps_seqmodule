##
# File:  PathInfo.py
# Date:  25-Feb-2013
#
# Updated:
#  26-Feb-2013  jdw   implement new session storage options in DataFileReference()
## 
"""
Common methods for finding path information for sequence module resource and data files.

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"

import sys, os, os.path
from wwpdb.api.facade.DataReference  import DataFileReference

class PathInfo(object):
    """ Common methods for finding path information for sequence resources and data files.
    """    
    def __init__(self, siteId="WWPDB_DEPLOY_TEST", sessionPath='.', verbose=False, log=sys.stderr):
        """ 
        """
        self.__verbose=verbose
        self.__lfh=log
        self.__debug=False
        self.__siteId=siteId
        self.__sessionPath = sessionPath
        #
        
    def getModelPdbxFilePath(self,dataSetId,wfInstanceId=None,fileSource="archive",versionId="latest"):
        return self.__getStandardPath(dataSetId=dataSetId,
                                      wfInstanceId=wfInstanceId,
                                      fileSource=fileSource,
                                      versionId=versionId,
                                      contentType='model',
                                      formatType='pdbx')
    
    def getModelPdbFilePath(self,dataSetId,wfInstanceId=None,fileSource="archive",versionId="latest"):
        return self.__getStandardPath(dataSetId=dataSetId,
                                      wfInstanceId=wfInstanceId,
                                      fileSource=fileSource,
                                      versionId=versionId,
                                      contentType='model',
                                      formatType='pdb')
    
    def getPolyLinkFilePath(self,dataSetId,wfInstanceId=None,fileSource="archive",versionId="latest"):
        return self.__getStandardPath(dataSetId=dataSetId,
                                      wfInstanceId=wfInstanceId,
                                      fileSource=fileSource,
                                      versionId=versionId,
                                      contentType='polymer-linkage-distances',
                                      formatType='pdbx')

    def getSequenceStatsFilePath(self,dataSetId,wfInstanceId=None,fileSource="archive",versionId="latest"):
        return self.__getStandardPath(dataSetId=dataSetId,
                                      wfInstanceId=wfInstanceId,
                                      fileSource=fileSource,
                                      versionId=versionId,
                                      contentType='seq-data-stats',
                                      formatType='pic')


    def getReferenceSequenceFilePath(self,dataSetId,entityId='1',wfInstanceId=None,fileSource="archive",versionId="latest"):
        return self.__getStandardPath(dataSetId=dataSetId,
                                      wfInstanceId=wfInstanceId,
                                      fileSource=fileSource,
                                      versionId=versionId,
                                      partNumber=entityId,
                                      contentType='seqdb-match',
                                      formatType='pdbx')

    def getSequenceAssignmentFilePath(self,dataSetId,wfInstanceId=None,fileSource="archive",versionId="latest"):
        return self.__getStandardPath(dataSetId=dataSetId,
                                      wfInstanceId=wfInstanceId,
                                      fileSource=fileSource,
                                      versionId=versionId,
                                      contentType='seq-assign',
                                      formatType='pdbx')

    def getBlastMatchFilePath(self,entityId,entityPart='1',wfInstanceId=None,fileSource="archive",versionId="latest"):
        return self.__getStandardPath(dataSetId=entityId,
                                      wfInstanceId=wfInstanceId,
                                      fileSource=fileSource,
                                      versionId=versionId,
                                      partNumber=entityPart,
                                      contentType='blast-match',
                                      formatType='xml')


    def __getStandardPath(self,dataSetId,wfInstanceId=None,contentType=None,formatType=None,fileSource="archive",versionId="latest",partNumber='1'):
        """   Get WF conforming path/file names.
        """                
        #
        retPath=None
        if (fileSource in ['archive','wf-archive']):
            dfRef=DataFileReference(siteId=self.__siteId,verbose=self.__verbose,log=self.__lfh)
            dfRef.setDepositionDataSetId(dataSetId)
            dfRef.setStorageType('archive')
            dfRef.setContentTypeAndFormat(contentType,formatType)
            dfRef.setPartitionNumber(partNumber)
            dfRef.setVersionId(versionId)
        elif (fileSource =='wf-instance'):
            dfRef=DataFileReference(siteId=self.__siteId,verbose=self.__verbose,log=self.__lfh)
            dfRef.setDepositionDataSetId(dataSetId)
            dfRef.setWorkflowInstanceId(wfInstanceId)
            dfRef.setStorageType('wf-instance')
            dfRef.setContentTypeAndFormat(contentType,formatType)
            dfRef.setPartitionNumber(partNumber)
            dfRef.setVersionId(versionId) 
        elif (fileSource in ['session','wf-session']):
            dfRef=DataFileReference(siteId=self.__siteId,verbose=self.__verbose,log=self.__lfh)
            dfRef.setSessionPath(self.__sessionPath)
            dfRef.setSessionDataSetId(dataSetId)
            dfRef.setStorageType('session')
            dfRef.setContentTypeAndFormat(contentType,formatType)
            dfRef.setPartitionNumber(partNumber)
            dfRef.setVersionId(versionId)

        else:
            self.__lfh.write("+PathInfo.__getStandardPath() bad file source %s for id %s wf id %s\n" % (fileSource,dataSetId,wfInstanceId))
            return retPath

        if (dfRef.isReferenceValid()):                  
            dP=dfRef.getDirPathReference()
            fP=dfRef.getFilePathReference()
            retPath=fP
            if (self.__debug):                
                self.__lfh.write("+PathInfo.__getStandardPath() directory path: %s\n" % dP)
                self.__lfh.write("+PathInfo.__getStandardPath() file path:      %s\n" % fP)
        else:
            self.__lfh.write("+PathInfo.__getStandardPath() invalid file path for %s for id %s wf id %s\n" % (fileSource,dataSetId,wfInstanceId))

        #
        return retPath
