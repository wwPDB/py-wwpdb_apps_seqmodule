##
# File:  SequenceDataImportWf.py
# Date:  21-Feb-2013
#
# Updates:
#   20-Apr-2010 jdw Ported to module seqmodule.
#   24-Apr-2010 jdw Add wwpdb-repository options
#   04-May-2010 jdw Add checks for existing stats file and polymer linkage data
#               jdw add doImportWithCheck()
#   08-May-2010 jdw add polymer linkage calculation to RCSB data import version
#   17-Jun-2010 Fix residue case and filter waters.
#   10-Aug-2010 RPS: SequenceDataImportWf.__loadSequenceDataStore() now logging 
#                    via self.__lfh.write() instead of sys.stderr.write()
#   18-Aug-2010 RPS: Updated log messages output by 
#                    SequenceDataImportWf.__getReferenceSequenceFilePath() and 
#                    SequenceDataImportWf__readReferenceSequences()
#                    regarding use of data files containing reference database
#                    sequence matches.
#   27-Jun-2011 jdw  switch RcsbDpUtility for WF version
#   21-Feb-2013 jdw  refactoring
##
"""
Import sequence and other required data for the sequence editor tool.
     
"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.09"

import sys, os, string, shutil

from wwpdb.apps.seqmodule.io.SiteInterface           import InterfaceSequenceRcsb, InterfaceArchiveRcsb
from wwpdb.apps.seqmodule.io.DepositionDataFile      import DepositionDataFile, CifFile, ReferenceSequenceDataFile
# 
from wwpdb.apps.seqmodule.io.DepositionDataFilePdbx  import DepositionDataFilePdbx, PdbxFile, ReferenceSequenceDataFilePdbx,PolymerLinkageDistanceFilePdbx
from wwpdb.apps.seqmodule.io.SequenceDataStore       import SequenceDataStore
from wwpdb.apps.seqmodule.util.SequenceReferenceData import SequenceReferenceData
from wwpdb.apps.seqmodule.util.SequenceLabel         import SequenceFeature
from wwpdb.apps.seqmodule.util.SequenceExamples      import SequenceExamples

from wwpdb.utils.rcsb.FormatOut               import FormatOut
from wwpdb.utils.rcsb.RcsbDpUtil              import RcsbDpUtil
# for WF version
from wwpdb.utils.rcsb.RcsbDpUtility           import RcsbDpUtility
from wwpdb.api.facade.ConfigInfo              import ConfigInfo

from wwpdb.utils.pair_align.wrapper.libPairwiseAlignPackage import PairwiseAlign, ostream

from wwpdb.api.facade.DataReference  import DataFileReference

def editDistance(s1, s2):
    """
    Compute the Damerau-Levenshtein distance between two given
    strings (s1 and s2)
    """    
    d = {}
    lenstr1 = len(s1)
    lenstr2 = len(s2)
    for i in xrange(-1,lenstr1+1):
        d[(i,-1)] = i+1
    for j in xrange(-1,lenstr2+1):
        d[(-1,j)] = j+1
 
    for i in xrange(0,lenstr1):
        for j in xrange(0,lenstr2):
            if s1[i] == s2[j]:
                cost = 0
            else:
                cost = 1
            d[(i,j)] = min(
                           d[(i-1,j)] + 1, # deletion
                           d[(i,j-1)] + 1, # insertion
                           d[(i-1,j-1)] + cost, # substitution
                          )
            if i>1 and j>1 and s1[i]==s2[j-1] and s1[i-1] == s2[j]:
                d[(i,j)] = min (d[(i,j)], d[i-2,j-2] + cost) # transposition
 
    return d[lenstr1-1,lenstr2-1]

def multikeysort(items, columns):
    """
    Sort list of dictionaries on multiple keys -

    Example 
    uu = [{"a": 1, "b": 2, "c": "tiger" },
    {"a": 2, "b": 1, "c": "tiger" },
    {"a": 3, "b": 5, "c": "bear" },
    {"a": 4, "b": 4, "c": "tiger" },
    {"a": 5, "b": 1, "c": "bear" }
    ]
    result = multikeysort(uu, ['c', 'b', 'a'])
    """    
    from operator import itemgetter
    comparers = [ ((itemgetter(col[1:].strip()), -1) if col.startswith('-') else (itemgetter(col.strip()), 1)) for col in columns]
    def sign(a, b):
        if   a < b:  return -1
        elif a > b:  return 1
        else:        return 0
    def comparer(left, right):
        for fn, mult in comparers:
            result = sign(fn(left), fn(right))
            if result:
                return mult * result
        else:
            return 0
    return sorted(items, cmp=comparer)


class SequenceDataImportWf(object):
    """
     This class encapsulates all of the data import operations 
     of sequence and other data from the wwPDB Workflow data pipeline required
     for the sequence editor tool.

     Storage model - imported data is loaded into the sequence data store
                     where it is managed by the SequenceDataStore() class. 

     Methods in this class extract the author and coordinate sequence data
     and source information from PDBx CIF.

     Reference sequence data is extracted from the processed BLAST results
     (ie. wwPDB seqdb-match data files) for each polymer entity.

     A PDB format file is produced from the current PDBx CIF file
     via software translation (Maxit).  This file is required for the 3D
     visualization.
     
    """
    def __init__(self,reqObj=None, depDataSetId=None, wfInstanceId=None, fileSource='archive', verbose=False,log=sys.stderr,fetchPdbFile=True):

        self.__reqObj=reqObj
        self.__depDataSetId=depDataSetId
        self.__wfInstanceId=wfInstanceId        
        self.__fileSource=fileSource
        self.__verbose=verbose        
        self.__lfh=log
        self.__fetchPdbFile=fetchPdbFile        
        #
        self.__sequenceDataStorePath=None
        #
        self.__debug=False        
        self.__pdbId=None

        #
        self.__sessionObj=None        
        self.__sessionPath='.'
        self.__cifPath=None
        self.__cifEpsPath=None
        self.__pdbPath=None
        self.__polyLinkPath=None
        self.__polymerLinkDistList=[]
        #
        self.__seqDataStatsPath=None
        #
        self.__polyEntityList=[]
        self.__entityD={}
        self.__eRefD={}
        self.__chD={}
        self.__statsLoadedFromRestart=False
        #
        #
        self.__setup()
        
    def __setup(self):
        #
        self.__sessionObj  = self.__reqObj.getSessionObj()
        self.__sessionPath = self.__sessionObj.getPath()

        #
        if (self.__verbose):
            self.__lfh.write("+SequenceDataImportWf.__setup() - session id %s \n" % (self.__sessionObj.getId()))
            self.__lfh.write("+SequenceDataImportWf.__setup() - session path %s\n" % (self.__sessionObj.getPath()))
            self.__lfh.write("+SequenceDataImportWf.__setup() - Dep Data set ID  %s file source %s\n" % (self.__depDataSetId, self.__fileSource))
            self.__lfh.write("+SequenceDataImportWf.__setup() - operation  %s\n" % self.__reqObj.getValue("operation"))
            self.__lfh.write("+SequenceDataImportWf.__setup() - identifier %s\n" % self.__reqObj.getValue("identifier"))
            #
            #self.__lfh.write("+SequenceDataImportWf.__setup() - UploadFileType  %s\n" % self.__reqObj.getValue("UploadFileType"))
            #self.__lfh.write("+SequenceDataImportWf.__setup() - UploadFileName  %s\n" % self.__reqObj.getValue("UploadFileName"))
            #self.__lfh.write("+SequenceDataImportWf.__setup() - RcsbDataPath  %s\n" % self.__rcsbDataPath)
            #self.__lfh.write("+SequenceDataImportWf.__setup() - RcsbReferenceSequencePath  %s\n" % self.__rcsbReferenceSequencePath)

            self.__lfh.flush()
            self.__getModelPdbxFilePath()
            self.__getModelPdbFilePath()
            #
            self.__getSequenceStatsFilePath()
            self.__getPolyLinkFilePath()
            if self.__fetchPdbFile:
                self.__fetchModelPdbFile()
                
            self.__lfh.flush()            


    def statsLoaded(self):
        return (self.__statsLoadedFromRestart)
    
    def getAtypicalPolymerLinkages(self):
        return self.__polymerLinkDistList
        
    def getSequenceDataStorePath(self):
        return self.__sequenceDataStorePath
    
    def doImport(self):
        """ Perform all data import operations  - 

            Extract the author and coordinate sequence data
            and source information from the current RCSB internal CIF file.

            Coordinate sequence data is extracted from the encapsulated coordinate data.

            Reference sequence data is extracted from the processed BLAST results
            (ie. rcsb pre-blast data files) for each polymer entity.

           
        """

        self.__readPolymerLinkageDistances()
        if (self.__verbose):
            self.__lfh.write("+SequenceDataImportWf.__doImport() polymer linkage atypical distance count is %d\n" % len(self.__polymerLinkDistList))

            
        # first extract data from CIF file if this provide -
        if (len(self.__cifPath) > 0):
            # read entity details from the CIF -
            ft="CIF"
            if (self.__verbose):
                self.__lfh.write("+SequenceDataImportWf.__doImport() (%s) %s coordinate extracted from data from file %s\n"
                                 % (ft,self.__depDataSetId, self.__cifPath ))                                    
            self.__readCIF(self.__cifPath)
            self.__getEntityDetails()
            #
            #  First try to get the residue sequence from the CIF -
            #
            if (self.__sdf.getResidueTableLength() > 20):
                self.__getCoordinateSequences()
                if (self.__verbose):
                    self.__lfh.write("+SequenceDataImportWf.__doImport() (%s) %s coordinate sequences extracted from residue table for %d chains\n"
                                                  % (ft,self.__depDataSetId,len(self.__chD)))
            else:
                #
                nAtoms=self.__sdf.getAtomSiteTableLength()
                if (nAtoms > 5):
                    if (self.__verbose):
                        self.__lfh.write("+SequenceDataImportWf.__doImport() (%s) %s atom site table length = %d\n" %
                                     (ft,self.__depDataSetId,nAtoms))                    
                    self.__chD=self.__sdf.getSequenceFromAtomSite()
                    if (self.__verbose):
                        for k,v in self.__chD.items():
                            self.__lfh.write("+SequenceDataImportWf.__doImport() chain %s sequence length %d\n" % (k,len(v)))
                else:
                    #
                    # Failing this try for the encapsulated coordinates in the CIF
                    #
                    xyzS=self.__sdf.getEncapsulatedCoordinates()
                    if (self.__verbose):
                        self.__lfh.write("+SequenceDataImportWf.__doImport() (%s) %s coordinate data string length = %d\n" %
                                     (ft,self.__depDataSetId,len(xyzS)))
                    self.__getCoordinateSequencesFromString(xyzS)
                    if (self.__verbose):
                        self.__lfh.write("+SequenceDataImportWf.__doImport() (%s) %s coordinate sequences extracted from string for %d chains\n"
                                                  % (ft,self.__depDataSetId,len(self.__chD)))                        
        else:
            # no coordinate file found
            if (self.__verbose):
                self.__lfh.write("+SequenceDataImportWf.__doImport() NO coordinate files for data import.\n")
            return False

        #
        self.__readReferenceSequences()


        self.__lfh.flush()                
        self.__loadSequenceDataStore()
        self.__lfh.flush()                
        return True

    def doImportWithCheck(self):
        """ Perform all data import operations  - 

            Extract the author and coordinate sequence data
            and source information from the current RCSB internal CIF file.

            Coordinate sequence data is extracted from the encapsulated coordinate data.

            Reference sequence data is extracted from the processed BLAST results
            (ie. rcsb pre-blast data files) for each polymer entity.

           *NOTE*  this module is for use with the live editor version  jdw
            
        """

        self.__readPolymerLinkageDistances()
        if (self.__verbose):
            self.__lfh.write("+SequenceDataImportWf.__doImport() polymer linkage atypical distance count is %d\n" % len(self.__polymerLinkDistList))


        if (os.access( self.__seqDataStatsPath,os.R_OK)):
            fileName='sequenceDataStore.pic'            
            filePath = os.path.join(self.__sessionPath,fileName)
            shutil.copyfile(self.__seqDataStatsPath,filePath)
            if (self.__verbose):
                self.__lfh.write("+SequenceDataImportWf.__doImport() recovered existing sequence data store %s\n" % filePath)
            # jump out here using the the copy of the sequence data store file
            self.__statsLoadedFromRestart=True
            self.__sequenceDataStorePath=filePath
            return True
            

            
        # first extract data from CIF file if this provide -
        if (len(self.__cifPath) > 0):
            # read entity details from the CIF -
            ft="CIF"
            if (self.__verbose):
                self.__lfh.write("+SequenceDataImportWf.__doImport() (%s) %s coordinate extracted from data from file %s\n"
                                 % (ft,self.__depDataSetId, self.__cifPath ))                                    
            self.__readCIF(self.__cifPath)
            self.__getEntityDetails()
            #
            #  First try to get the residue sequence from the CIF -
            #
            if (self.__sdf.getResidueTableLength() > 20):
                self.__getCoordinateSequences()
                if (self.__verbose):
                    self.__lfh.write("+SequenceDataImportWf.__doImport() (%s) %s coordinate sequences extracted from residue table for %d chains\n"
                                                  % (ft,self.__depDataSetId,len(self.__chD)))
            else:
                #
                nAtoms=self.__sdf.getAtomSiteTableLength()
                if (nAtoms > 5):
                    if (self.__verbose):
                        self.__lfh.write("+SequenceDataImportWf.__doImport() (%s) %s atom site table length = %d\n" %
                                     (ft,self.__depDataSetId,nAtoms))                    
                    self.__chD=self.__sdf.getSequenceFromAtomSite()
                    if (self.__verbose):
                        for k,v in self.__chD.items():
                            self.__lfh.write("+SequenceDataImportWf.__doImport() chain %s sequence length %d\n" % (k,len(v)))
                else:
                    #
                    # Failing this try for the encapsulated coordinates in the CIF
                    #
                    xyzS=self.__sdf.getEncapsulatedCoordinates()
                    if (self.__verbose):
                        self.__lfh.write("+SequenceDataImportWf.__doImport() (%s) %s coordinate data string length = %d\n" %
                                     (ft,self.__depDataSetId,len(xyzS)))
                    self.__getCoordinateSequencesFromString(xyzS)
                    if (self.__verbose):
                        self.__lfh.write("+SequenceDataImportWf.__doImport() (%s) %s coordinate sequences extracted from string for %d chains\n"
                                                  % (ft,self.__depDataSetId,len(self.__chD)))                        
        else:
            # no coordinate file found
            if (self.__verbose):
                self.__lfh.write("+SequenceDataImportWf.__doImport() NO coordinate files for data import.\n")
            return False

        #
        self.__readReferenceSequences()
        self.__lfh.flush()                
        self.__loadSequenceDataStore()
        self.__lfh.flush()                
        return True


    def __getModelPdbxFilePath(self):
        """Get model file path (PDBx format).
        """                
        #
        # Get PDBx model file -
        #
        if (self.__fileSource in ['archive','wf-archive']):
            dfRef=DataFileReference()
            dfRef.setDepositionDataSetId(self.__depDataSetId)
            dfRef.setStorageType('archive')
            dfRef.setContentTypeAndFormat('model','pdbx')
            dfRef.setVersionId('latest')
              
            if (dfRef.isReferenceValid()):                  
                dP=dfRef.getDirPathReference()
                fP=dfRef.getFilePathReference()
                if (self.__verbose):                
                    self.__lfh.write("+SequenceDataImportWf.__setup() PDBx model directory path: %s\n" % dP)
                    self.__lfh.write("+SequenceDataImportWf.__setup() PDBx model file      path: %s\n" % fP)
            else:
                self.__lfh.write("+SequenceDataImportWf.__setup() Bad archival reference for id %s \n" % self.__depDataSetId)

        elif (self.__fileSource =='wf-instance'):
            dfRef=DataFileReference()
            dfRef.setDepositionDataSetId(self.__depDataSetId)
            dfRef.setWorkflowInstanceId(self.__wfInstanceId)            
            dfRef.setStorageType('wf-instance')
            dfRef.setContentTypeAndFormat('model','pdbx')
            dfRef.setVersionId('latest')            

            if (dfRef.isReferenceValid()):
                dP=dfRef.getDirPathReference()
                fP=dfRef.getFilePathReference()
                if (self.__verbose):                                
                    self.__lfh.write("+SequenceDataImportWf.__setup() PDBx model directory path: %s\n" % dP)
                    self.__lfh.write("+SequenceDataImportWf.__setup() PDBx nodel file      path: %s\n" % fP)                
            else:
                self.__lfh.write("+SequenceDataImportWf.__setup() Bad wf-instance reference for id %s wf id %s\n" %
                                 (self.__depDataSetId,self.__wfInstanceId))                

        else:
            self.__lfh.write("+SequenceDataImportWf.__setup() Bad file source for id %s wf id %s\n" %
                             (self.__depDataSetId,self.__wfInstanceId))
        #
        self.__lfh.flush()        
        #
        self.__cifPath=fP

    def __getModelPdbFilePath(self):
        """Get model file path (PDB format).
        """        
        #
        # Get PDB model file -
        #
        if (self.__fileSource in ['archive', 'wf-archive']):
            dfRef=DataFileReference()
            dfRef.setDepositionDataSetId(self.__depDataSetId)
            dfRef.setStorageType('archive')
            dfRef.setContentTypeAndFormat('model','pdb')
            dfRef.setVersionId('latest')
              
            if (dfRef.isReferenceValid()):                  
                dP=dfRef.getDirPathReference()
                fP=dfRef.getFilePathReference()
                sP=dfRef.getSitePrefix()
                if (self.__verbose):
                    self.__lfh.write("+SequenceDataImportWf.__setup() site prefix             : %s\n" % sP)
                    self.__lfh.write("+SequenceDataImportWf.__setup() PDB model directory path: %s\n" % dP)
                    self.__lfh.write("+SequenceDataImportWf.__setup() PDB model file      path: %s\n" % fP)
            else:
                self.__lfh.write("+SequenceDataImportWf.__setup() Bad archival reference for id %s \n" % self.__depDataSetId)

        elif (self.__fileSource =='wf-instance'):
            dfRef=DataFileReference()
            dfRef.setDepositionDataSetId(self.__depDataSetId)
            dfRef.setWorkflowInstanceId(self.__wfInstanceId)            
            dfRef.setStorageType('wf-instance')
            dfRef.setContentTypeAndFormat('model','pdb')
            dfRef.setVersionId('latest')
            
            if (dfRef.isReferenceValid()):
                dP=dfRef.getDirPathReference()
                fP=dfRef.getFilePathReference()
                sP=dfRef.getSitePrefix()                
                if (self.__verbose):
                    self.__lfh.write("+SequenceDataImportWf.__setup() site prefix             : %s\n" % sP)                    
                    self.__lfh.write("+SequenceDataImportWf.__setup() PDB model directory path: %s\n" % dP)
                    self.__lfh.write("+SequenceDataImportWf.__setup() PDB model file      path: %s\n" % fP)                
            else:
                self.__lfh.write("+SequenceDataImportWf.__setup() Bad wf-instance reference for id %s wf id %s\n" %
                                 (self.__depDataSetId,self.__wfInstanceId))                

        else:
            self.__lfh.write("+SequenceDataImportWf.__setup() Bad file soure for id %s wf id %s\n" %
                             (self.__depDataSetId,self.__wfInstanceId))
        self.__lfh.flush()
        #
        #
        self.__pdbPath=fP

        self.__lfh.flush()

    def __fetchModelPdbFile(self):
        """  Make a local copy of the PDB model file for 3D visualization -- 
        """
        #
        # Copy the latest version of the model file from the workflow storage if this exists
        #
        doCopy=True
        if (self.__verbose):
            self.__lfh.write("+SequenceDataImportWf.__fetchModelPdbFile() - fetching/creating PDB model file in session %s\n" % self.__sessionPath)
        
        if (doCopy and (self.__pdbPath is not None) and (len(self.__pdbPath) > 0) and (os.access(self.__pdbPath,os.R_OK))):
            try:
                dst=os.path.join(self.__sessionPath,str(self.__depDataSetId).lower()+".pdb")                            
                shutil.copyfile(self.__pdbPath,dst)
                if (self.__verbose):
                    self.__lfh.write("+SequenceDataImportWf.__fetchModelPdbFile() - copied PDB model file in: %s\n" % dst)                
            except:
                self.__lfh.write("+SequenceDataImportWf.__fetchModelPdbFile() PDB file copy failed: %s\n" % fP)                                
                pass
        elif ((self.__cifPath is not None) and (len(self.__cifPath)>0) and os.access(self.__cifPath,os.R_OK) ):
            pdbPath=os.path.join(self.__sessionPath,str(self.__depDataSetId).lower()+".pdb")
            logPath=os.path.join(self.__sessionPath,str(self.__depDataSetId).lower()+"-cif2pdb.log")            
            siteId=self.__reqObj.getValue("WWPDB_SITE_ID")           
            dp=RcsbDpUtility(tmpPath=self.__sessionPath,siteId=siteId,verbose=self.__verbose, log=self.__lfh)
            dp.imp(self.__cifPath)
            dp.op("cif2pdb")
            dp.exp(pdbPath)
            dp.expLog(logPath)            
            dp.cleanup()
            if (self.__verbose):
                self.__lfh.write("+SequenceDataImportWf.__fetchModelPdbFile() - site %s created PDB model file in: %s\n" % (siteId,pdbPath))
        #
        #
        
        
    def __getPolyLinkFilePath(self):
        """Get polymer linkage file path.
        """
        if (self.__fileSource in ['archive','wf-archive']):
            dfRef=DataFileReference()
            dfRef.setDepositionDataSetId(self.__depDataSetId)
            dfRef.setStorageType('archive')
            dfRef.setContentTypeAndFormat('polymer-linkage-distances','pdbx')
            dfRef.setVersionId('latest')
              
            if (dfRef.isReferenceValid()):                  
                dP=dfRef.getDirPathReference()
                fP=dfRef.getFilePathReference()
                if (self.__verbose):                
                    self.__lfh.write("+SequenceDataImportWf.__setup() Polymer linkage file directory path: %s\n" % dP)
                    self.__lfh.write("+SequenceDataImportWf.__setup() Polymer linkage file           path: %s\n" % fP)
            else:
                self.__lfh.write("+SequenceDataImportWf.__setup() Bad archival reference for id %s \n" % self.__depDataSetId)

        elif (self.__fileSource =='wf-instance'):
            dfRef=DataFileReference()
            dfRef.setDepositionDataSetId(self.__depDataSetId)
            dfRef.setWorkflowInstanceId(self.__wfInstanceId)            
            dfRef.setStorageType('wf-instance')
            dfRef.setContentTypeAndFormat('polymer-linkage-distances','pdbx')
            dfRef.setVersionId('latest')            

            if (dfRef.isReferenceValid()):
                dP=dfRef.getDirPathReference()
                fP=dfRef.getFilePathReference()
                if (self.__verbose):                                
                    self.__lfh.write("+SequenceDataImportWf.__setup() Polymer linkage directory path: %s\n" % dP)
                    self.__lfh.write("+SequenceDataImportWf.__setup() Polymer linkage file      path: %s\n" % fP)                
            else:
                self.__lfh.write("+SequenceDataImportWf.__setup() Bad wf-instance reference for id %s wf id %s\n" %
                                 (self.__depDataSetId,self.__wfInstanceId))                

        else:
            self.__lfh.write("+SequenceDataImportWf.__setup() Bad file source for id %s wf id %s\n" %
                             (self.__depDataSetId,self.__wfInstanceId))
        #
        self.__lfh.flush()        
        #
        self.__polyLinkPath=fP


    def __getSequenceStatsFilePath(self):
        """Get sequence stats file path.
        """
        if (self.__fileSource in ['archive','wf-archive']):
            dfRef=DataFileReference()
            dfRef.setDepositionDataSetId(self.__depDataSetId)
            dfRef.setStorageType('archive')
            dfRef.setContentTypeAndFormat('seq-data-stats','pic')
            dfRef.setVersionId('latest')
              
            if (dfRef.isReferenceValid()):                  
                dP=dfRef.getDirPathReference()
                fP=dfRef.getFilePathReference()
                if (self.__verbose):                
                    self.__lfh.write("+SequenceDataImportWf.__setup() Sequence data stats file directory path: %s\n" % dP)
                    self.__lfh.write("+SequenceDataImportWf.__setup() Sequence data file           path: %s\n" % fP)
            else:
                self.__lfh.write("+SequenceDataImportWf.__setup() Bad archival reference for id %s \n" % self.__depDataSetId)

        elif (self.__fileSource =='wf-instance'):
            dfRef=DataFileReference()
            dfRef.setDepositionDataSetId(self.__depDataSetId)
            dfRef.setWorkflowInstanceId(self.__wfInstanceId)            
            dfRef.setStorageType('wf-instance')
            dfRef.setContentTypeAndFormat('seq-data-stats','pic')
            dfRef.setVersionId('latest')            

            if (dfRef.isReferenceValid()):
                dP=dfRef.getDirPathReference()
                fP=dfRef.getFilePathReference()
                if (self.__verbose):                                
                    self.__lfh.write("+SequenceDataImportWf.__setup() Sequence data stats directory path: %s\n" % dP)
                    self.__lfh.write("+SequenceDataImportWf.__setup() Sequence data stats file      path: %s\n" % fP)                
            else:
                self.__lfh.write("+SequenceDataImportWf.__setup() Bad wf-instance reference for id %s wf id %s\n" %
                                 (self.__depDataSetId,self.__wfInstanceId))                

        else:
            self.__lfh.write("+SequenceDataImportWf.__setup() Bad file source for id %s wf id %s\n" %
                             (self.__depDataSetId,self.__wfInstanceId))
        #
        self.__lfh.flush()        
        #
        self.__seqDataStatsPath=fP

    def __readCIF(self,cifPath):
        pdbxFile = PdbxFile(cifPath)
        self.__sdf = DepositionDataFilePdbx(pdbxFile)


    def __getCoordinateSequences(self):
        """
        """
        self.__chD={}
        for id,eD in self.__entityD.items():
            chainIdList=eD['CHAIN_LIST']
            for chainId in chainIdList:
                sL=self.__sdf.getCoordinateSequenceList(chainId)
                self.__chD[chainId]=sL
            
            
    def __getCoordinateSequencesFromString(self, xyzS):
        """ Extract the polymer sequence from the PDB coordinate records in the input
            string.

            self.__chD{} =[(a3,orglabel index, comment/details,  align index placeholder),(),..] is a dictionary with chain id key
                          containing sequences stored as a list of tuples.

            Each tuple has the value (3-letter-code, original label residue index, comment/details, alignment index )            
        """

        #if (self.__verbose):
        #    self.__lfh.write("%s\n" % xyzS)        
        #        ATOM   1151  ND1 HIS A  76A
        #        0123456789012345678901234567890
        #                  1         2
        tupPrev     = ('','','')
        chainIdPrev = ''
        firstChain = True
        rList = []
        self.__chD={}
        idx=1
        for line in xyzS.splitlines():
            if line[0:5] == "ATOM ":
                chainId = str(line[21:22]).strip()
                if firstChain:
                    firstChain=False
                    chainIdPrev=chainId
                if chainId != chainIdPrev:
                    # new chain
                    self.__chD[chainIdPrev]=rList
                    chainIdPrev = chainId
                    rList=[]
                else:
                    resId   = str(line[17:20]).strip()
                    resNo   = str(line[22:27]).strip()
                    tup=(resId,resNo,'')                    
                    if tup != tupPrev:
                        rList.append((resId,resNo,'',idx))
                        idx+=1
                        tupPrev = tup

            elif line[0:7] == "ENDMDL":
                break

        if not self.__chD.has_key(chainIdPrev):
            self.__chD[chainIdPrev]=rList


    def __getEntityDetails(self):
        self.__pdbId = self.__sdf.getDbCode("PDB")        
        #self.__polyEntityList=self.__sdf.getPolymerEntityList('polymer')
        self.__polyEntityList=self.__sdf.getEntityPolyList()        
        self.__entityD={}
        for entityId in self.__polyEntityList:
            if (self.__verbose): self.__lfh.write("_SequenceDataImportWf.__getEntityDetails() entity id= %s\n" % entityId)
            sD={}
            sD['CHAIN_LIST'] = self.__sdf.getPdbChainIdList(entityId)
            if (self.__verbose): self.__lfh.write("_SequenceDataImportWf.__getEntityDetails() entity id= %s chain list %r\n" %
                                                  (entityId,sD['CHAIN_LIST']))
            xS = str(self.__sdf.getSequence(entityId)).upper()
            sD['SEQ_ENTITY_1']=string.join(xS.split(),"").upper()
            if (self.__verbose): self.__lfh.write("_SequenceDataImportWf.__getEntityDetails() entity id= %s sequence %s\n"
                                                  % (entityId,sD['SEQ_ENTITY_1']) )
            sD['SOURCE_METHOD']=str(self.__sdf.getSourceMethod(entityId)).upper()
            if (self.__verbose): self.__lfh.write("_SequenceDataImportWf.__getEntityDetails() entity id= %s source method %s\n"
                                                  % (entityId,sD['SOURCE_METHOD']) )
                    
            sD['POLYMER_TYPE']=self.__sdf.getPolymerEntityType(entityId)
            name=''
            strain=''
            taxId=''
            if (len(sD['SOURCE_METHOD']) > 1): 
                if sD['SOURCE_METHOD'] == "NAT":
                    (name,strain,taxId)=self.__sdf.getSourceNatDetails(entityId)
                elif sD['SOURCE_METHOD'] == "MAN":
                    (name,strain,taxId)=self.__sdf.getSourceGenDetails(entityId)
                elif sD['SOURCE_METHOD'] == "SYN":
                    (name,strain,taxId)=self.__sdf.getSourceSynDetails(entityId)
            else:
                (name,strain,taxId)=self.__sdf.getSourceGenDetails(entityId)
                if len(name) < 2:
                    (name,strain,taxId)=self.__sdf.getSourceNatDetails(entityId)
                    if len(name) < 2:
                        (name,strain,taxId)=self.__sdf.getSourceSynDetails(entityId) 
                    
            sD['SOURCE_NAME']=name
            sD['SOURCE_STRAIN']=strain
            sD['SOURCE_TAXID']=taxId
            #
            # Check missing polymer type assignment 
            #
            if len(sD['POLYMER_TYPE']) < 2:
                sdr=SequenceReferenceData(self.__verbose,self.__lfh)                
                sD['POLYMER_TYPE']=sdr.guessPolymerType(sD['SEQ_ENTITY_1'])
            #
            if (self.__verbose): self.__lfh.write("_SequenceDataImportWf.__getEntityDetails() entity id= %s polymer type %s\n"
                                                  % (entityId,sD['POLYMER_TYPE']) )            
            
            self.__entityD[entityId]=sD


    def __getReferenceSequenceFilePath(self,eId):
        fP=None
        dP=None
        #
        # Get file path for matching reference sequence data.
        #
        if (self.__fileSource in ['archive','wf-archive']):
            dfRef=DataFileReference()
            dfRef.setDepositionDataSetId(self.__depDataSetId)
            dfRef.setStorageType('archive')
            dfRef.setContentTypeAndFormat('seqdb-match','pdbx')
            dfRef.setPartitionNumber(eId)
            dfRef.setVersionId('latest')
              
            if (dfRef.isReferenceValid()):                  
                dP=dfRef.getDirPathReference()
                fP=dfRef.getFilePathReference()
                if (self.__verbose):                
                    self.__lfh.write("+SequenceDataImportWf.__setup() Matching reference sequences--seeking directory path: %s\n" % dP)
                    self.__lfh.write("+SequenceDataImportWf.__setup() Matching reference sequences--seeking file      path: %s\n" % fP)
            else:
                self.__lfh.write("+SequenceDataImportWf.__setup() Bad archival seqdb-match reference for id %s \n" % self.__depDataSetId)

        elif (self.__fileSource =='wf-instance'):
            dfRef=DataFileReference()
            dfRef.setDepositionDataSetId(self.__depDataSetId)
            dfRef.setWorkflowInstanceId(self.__wfInstanceId)            
            dfRef.setStorageType('wf-instance')
            dfRef.setContentTypeAndFormat('seqdb-match','pdbx')
            dfRef.setVersionId('latest')            
            dfRef.setPartitionNumber(eId)
            if (dfRef.isReferenceValid()):
                dP=dfRef.getDirPathReference()
                fP=dfRef.getFilePathReference()
                if (self.__verbose):                                
                    self.__lfh.write("+SequenceDataImportWf.__setup() Matching reference sequences--seeking directory path: %s\n" % dP)
                    self.__lfh.write("+SequenceDataImportWf.__setup() Matching reference sequences--seeking file      path: %s\n" % fP)                
            else:
                self.__lfh.write("+SequenceDataImportWf.__setup() Bad wf-instance seqdb-match reference for id %s wf id %s\n" %
                                 (self.__depDataSetId,self.__wfInstanceId))                

        else:
            self.__lfh.write("+SequenceDataImportWf.__setup() Bad file source for id %s wf id %s\n" %
                             (self.__depDataSetId,self.__wfInstanceId))

        self.__lfh.flush()        
        return fP

    def __readReferenceSequences(self):
        """  Read all the matching reference sequences and order these by taxonomy, identity, database(sp,tr,gb)
        """
        #
        self.__eRefD={}
        #
        for eId in self.__entityD.keys():
            #fn=sI.copyEntityFile(eId)
            fn=self.__getReferenceSequenceFilePath(eId)
            if ((fn is not None) and (len(fn) > 0) and os.access(fn,os.R_OK)):
                pdbxFile = PdbxFile(fn)
                self.__rsf = ReferenceSequenceDataFilePdbx(pdbxFile, verbose=self.__verbose,log=self.__lfh)
                tL=self.__rsf.getReferenceDetails()
                #self.__eRefD[eId] = multikeysort(tL, ['-db_name',  '-seq_sim', 'taxonomy_id'])
                self.__eRefD[eId] = multikeysort(tL, ['taxonomy_id', '-seq_sim', '-db_name'])
                # for ii,tD in enumerate( self.__eRefD[eId]):
                #     self.__lfh.write("+SequenceDataImportWf.__readReferenceSequences() return entity %s list element %d : %r\n\n" % (eId,ii,tD))                
            else:
                self.__lfh.write("+SequenceDataImportWf.__readReferenceSequences() File of matching ref DB sequences for entity %s cannot be read/found at: %s\n" % (eId,fn))
                self.__eRefD[eId]=[]


    def __readPolymerLinkageDistances(self):
        """  Read the file of polymer linkage distaces for the current model.
        """
        #
        try:
            fn=self.__polyLinkPath
            if ((fn is not None) and (len(fn) > 0) and os.access(fn,os.R_OK)):
                pdbxFile = PdbxFile(fn)
                self.__rsf = PolymerLinkageDistanceFilePdbx(pdbxFile)
                self.__polymerLinkDistList=self.__rsf. getPolymerLinkDistances()
            else:
                self.__polymerLinkDistList=[]
        except:
            self.__polymerLinkDistList=[]

    
    def __OldreadReferenceSequences(self):
        #
        self.__eRefD={}
        #
        for eId in self.__entityD.keys():
            #fn=sI.copyEntityFile(eId)
            fn=self.__getReferenceSequenceFilePath(eId)
            if ((fn is not None) and (len(fn) > 0) and os.access(fn,os.R_OK)):
                pdbxFile = PdbxFile(fn)
                self.__rsf = ReferenceSequenceDataFilePdbx(pdbxFile)
                tL=self.__rsf.getReferenceDetails()
                self.__eRefD[eId] = multikeysort(tL, ['-db_name',  '-seq_sim', 'taxonomy_id'])
            else:
                self.__eRefD[eId]=[]


            #        myList  = ['id', 'db_name', 'db_code', 'db_accession', 'match_length', 'queryFrom', 'queryTo',
            #           'hitFrom', 'hitTo', 'identity', 'positive', 'gaps', 'alignLen', 'query', 'subject',
            #           'midline', 'query_length', 'name', 'source_scientific', 'source_common', 'taxonomy_id',
            #           'gene', 'synonyms', 'comments', 'keyword', 'ec']

            authTaxid=str(self.__entityD[eId]['SOURCE_TAXID'])
            if len(authTaxid) > 1:
                for rS in self.__eRefD[eId]:
                    # reorder list by proximity to author supplied taxId
                    sTaxid=str(rS['taxonomy_id'])
                    if ((sTaxid is not None) and (len(sTaxid) > 1)):
                        #d=editDistance(sTaxid,authTaxid)
                        d=abs(int(authTaxid) - int(sTaxid))
                        rS['taxid_string_distance']=d
                        if (self.__debug):
                            self.__lfh.write("+SequenceDataImport.__readReferenceSequences() auth taxid %s reference taxid %s distance %d\n" %(authTaxid,sTaxid,d))
                    else:
                        d=100
                        rS['taxid_string_distance']=d
                        #if self.__debug:
                        #    self.__lfh.write("+SequenceDataImport.__readReferenceSequences() auth taxid %s without reference distance %d\n" %(authTaxid,d))               
                tL = self.__eRefD[eId]
                self.__eRefD[eId] = multikeysort(tL, ['-db_name', 'taxid_string_distance', '-seq_sim'])
                continue

            authSrc=str(self.__entityD[eId]['SOURCE_NAME'])
            if len(authSrc) > 1:
                dmin=10000
                for rS in self.__eRefD[eId]:
                    sName=rS['source_scientific']
                    if ((sName is not None) and (len(sName) > 1)):
                        d=editDistance(sName,authSrc)
                    else:
                        d=100
                    if (d < dmin):
                        dmin=d

                if (self.__verbose):
                    self.__lfh.write("+SequenceDataImport.__readReferenceSequences() auth source %s minimum distance is %d\n" %(authSrc,dmin))
                if (float(dmin)/float(len(authSrc)) > 0.25):
                    for rS in self.__eRefD[eId]:
                        # reorder list by proximity to author supplied source name
                        sName=rS['source_scientific']
                        if ((sName is not None) and (len(sName) > 1)):
                            d=editDistance(sName,authSrc)
                            rS['source_string_distance']=d
                            if (self.__verbose):
                                self.__lfh.write("+SequenceDataImport.__readReferenceSequences() auth source %s reference source %s distance %d\n" %(authSrc,sName,d))
                        else:
                            d=100
                            rS['source_string_distance']=d
                    tL = self.__eRefD[eId]
                    self.__eRefD[eId] = multikeysort(tL, ['-db_name', 'source_string_distance', '-seq_sim'])
                else:
                    tL = self.__eRefD[eId]
                    self.__eRefD[eId] = multikeysort(tL, ['-db_name', '-seq_sim', 'taxonomy_id'])


            

        
    def __loadSequenceDataStore(self):
        """ Store sequence data in a SequenceDataStore repository
            in the current session directory.
        """
        #
        if (self.__verbose):
            self.__lfh.write("+SequenceDataImportWf.__loadSequenceDataStore()  sessionId %s\n" %
                             (self.__sessionObj.getId()))        

        #
        sdr=SequenceReferenceData(self.__verbose,self.__lfh)
        groupDict={}
        chainType={}
        for id,eD in self.__entityD.items():
            polyType=eD['POLYMER_TYPE']
            polyTypeCode=sdr.getPolymerTypeCode(polyType)
            groupDict[id]=eD['CHAIN_LIST']
            for ch in eD['CHAIN_LIST']:
                chainType[ch]=polyTypeCode
                if (self.__verbose):
                    self.__lfh.write("+SequenceDataImportWf.__loadSequenceDataStore()  chain %s  type %s\n" %
                                     (ch,chainType[ch] ))        
                
                
        #
        seqFeature=SequenceFeature(verbose=self.__verbose)
        sds=SequenceDataStore(sessionObj=self.__sessionObj,verbose=self.__verbose,log=self.__lfh)
        sds.reset()

        #
        # Load coordinate sequences & features -
        #
        for id,S3L in self.__chD.items():
            #
            # Only store sequences if we know the polymer type details -
            if not chainType.has_key(id):
                continue
            sds.setSequence(S3L, id, 'xyz', altId=1, version=1)
            seqFeature.clear()
            seqFeature.setId(dbName='PDB',dbCode=self.__pdbId, dbAccession=self.__pdbId)
            seqFeature.setPolymerType(chainType[id])                            
            #
            sds.setFeature( seqFeature.get(),  id, 'xyz',  altId=1, version=1)                        


        #
        # Load author/entity sequence and features  -
        #

        for id,eD in self.__entityD.items():
            if (len(eD['CHAIN_LIST']) < 1):
                continue
            chainId0=eD['CHAIN_LIST'][0]
            polyType=eD['POLYMER_TYPE']
            polyTypeCode=sdr.getPolymerTypeCode(polyType)
            seqS=eD['SEQ_ENTITY_1']
            (r1L,r3L)=sdr.parseSequence(seqS,polyTypeCode)
            sTupL=[]
            ir=1
            for r3 in r3L:
                sTupL.append( (r3,str(ir),'',ir))
                ir += 1
            sds.setSequence(sTupL, chainId0, 'auth', altId=1, version=1)
            seqFeature.clear()
            seqFeature.setPolymerType(polyTypeCode)
            seqFeature.setId(dbName='PDB',dbCode=self.__pdbId, dbAccession=self.__pdbId)            
            seqFeature.setSource(organism=eD['SOURCE_NAME'],strain=eD['SOURCE_STRAIN'], taxid=eD['SOURCE_TAXID'])
            sds.setFeature( seqFeature.get(), chainId0, 'auth', altId=1, version=1)            
        #
        # Load reference sequences and features -
        #
        for eId,eD in self.__entityD.items():
            if (len(eD['CHAIN_LIST']) < 1):
                continue            
            chainId0=eD['CHAIN_LIST'][0]
            polyType=eD['POLYMER_TYPE']            
            polyTypeCode=sdr.getPolymerTypeCode(polyType)            
            rList = self.__eRefD[eId]

                
            altId=1
            for rD  in rList:
                myOrderId   = int(rD['id'])
                iBegin      = int(rD['hitFrom'])
                iEnd        = int(rD['hitTo'])                
                seqS=rD['subject']
                sTup3L= sdr.cnv1To3ListIdx(seqS,iBegin,polyTypeCode)
                sds.setSequence(sTup3L,  chainId0, 'ref',  altId=altId, version=1)
                seqFeature.clear()
                seqFeature.setId(dbName=rD['db_name'],dbCode=rD['db_code'], dbAccession=rD['db_accession'])
                seqFeature.setSource(organism=rD['source_scientific'], taxid=rD['taxonomy_id'])
                seqFeature.setItem('REF_MATCH_BEGIN',iBegin)
                seqFeature.setItem('REF_MATCH_END',iEnd)
                seqFeature.setItem('ORG_ORDER_ID', myOrderId)
                seqFeature.setPolymerType(polyTypeCode)                
                #
                seqFeature.setItem('AUTH_REF_SEQ_SIM_BLAST',rD['seq_sim'])

                sds.setFeature( seqFeature.get(),  chainId0, 'ref',  altId=altId, version=1)
                altId+=1

        #
        #        myList  = ['id', 'db_name', 'db_code', 'db_accession', 'match_length', 'queryFrom', 'queryTo',
        #           'hitFrom', 'hitTo', 'identity', 'positive', 'gaps', 'alignLen', 'query', 'subject',
        #           'midline', 'query_length', 'name', 'source_scientific', 'source_common', 'taxonomy_id',
        #           'gene', 'synonyms', 'comments', 'keyword', 'ec']
        #


        
        entryDict={}
        entryDict['PDB_ID'] =self.__pdbId
        entryDict['RCSB_ID']=self.__depDataSetId
        entryDict['DEPOSITION_DATA_SET_ID']=self.__depDataSetId        
        #
        for k,v in groupDict.items():
            sds.setGroup(k,v)
        for k,v in entryDict.items():
            sds.setEntryDetail(k,v)

        sds.serialize()
        self.__sequenceDataStorePath=sds.getFilePath()
        #
        if (self.__debug):
            sds.dump(self.__lfh)
            
        return {}


     
    def printIt(self):
        fOut=FormatOut()
        fOut.autoFormat("Entity dictionary", self.__entityD,3,3)
        fOut.autoFormat("Chain  dictionary", self.__chD,3,3)
        fOut.autoFormat("Reference dictionary", self.__eRefD,3,3)                
        fOut.writeStream(self.__lfh)
        fOut.clear()

