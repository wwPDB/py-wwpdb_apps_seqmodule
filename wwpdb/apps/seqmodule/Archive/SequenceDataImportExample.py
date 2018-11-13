##
# File:  SequenceDataImportExample.py
# Date:  21-Feb-2013
#
# Updates:
#
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

from wwpdb.io.misc.FormatOut               import FormatOut
from wwpdb.utils.rcsb.RcsbDpUtil              import RcsbDpUtil
# for WF version
from wwpdb.utils.dp.RcsbDpUtility           import RcsbDpUtility
from wwpdb.utils.config.ConfigInfo              import ConfigInfo

from wwpdb.utils.pair_align.wrapper.libPairwiseAlignPackage import PairwiseAlign, ostream

from wwpdb.wwpdb.utils.wf.DataReference  import DataFileReference

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


class SequenceDataImportExample(object):
    """
     This class loads the sequence data store with example data
     from the SequenceExamples() class.

     This class is provided primarily for testing and development.

     Storage model - imported data is loaded into the sequence data store
                     where it is managed by the SequenceDataStore() class. 

    """
    def __init__(self,reqObj=None, verbose=False,log=sys.stderr):
        self.__verbose=verbose
        self.__reqObj=reqObj
        self.__sessionObj=None
        self.__lfh=log
        self.__localTest=False
        #
        self.__sessionPath='.'
        #
        self.__setup()
        
    def __setup(self):
        try:
            self.__sessionObj  = self.__reqObj.getSessionObj()            
            self.__sessionPath = self.__sessionObj.getPath()
            if (self.__verbose):
                self.__lfh.write("+SequenceDataImportExampleRcsb.__setup() - session id %s\n" % self.__sessionObj.getId())
                self.__lfh.write("+SequenceDataImportExampleRcsb.__setup() - session path %s\n" %  self.__sessionObj.getPath() )
        except:
            pass

    def __copyExampleStructureFile(self,id, fileType="pdb"):
        
        fn=str(id).lower() + "." + fileType
        
        if (self.__verbose):
            self.__lfh.write("+SequenceDataImportExampleRcsb.__copyExampleStructureFile () - copying file name  %s\n" %  fn)                
        dst=os.path.join(self.__sessionPath,fn)
        topPath = self.__reqObj.getValue("TopPath")
        src=os.path.join(topPath,'xyz',fn)
        if (self.__verbose):
            self.__lfh.write("+SequenceDataImportExampleRcsb.__copyExampleStructureFile() - src  %s dst %s\n" %  (src,dst))
        try:
            shutil.copyfile(src,dst)
        except:
            pass
        return dst

    def doImport(self):
        return self.__loadExampleData()
    
    def __loadExampleData(self):
        """ Store example sequence data in a SequenceDataStore repository
            in the current session directory.
        """
        #
        if (self.__verbose):
            self.__lfh.write("+SequenceDataImportExampleRcsb.__loadExample()  sessionId %s\n" %
                             (self.__sessionObj.getId()))        

        # Load up some test data - 
        sE=SequenceExamples()

        refSL ={}
        xyzSL ={}
        authSL={}
        refFD ={}
        xyzFD = {}
        authFD = {}        
        for id in ['A','B']:
            refSL[id]  = sE.getRefSequenceWithIndexList(id)
            authSL[id] = sE.getAuthSequenceWithIndexList(id)
            xyzSL[id]  = sE.getXyzSequenceWithIndexList(id)
            refFD[id]  = sE.getRefFeatureDict(id)
            authFD[id] = sE.getAuthFeatureDict(id)
            xyzFD[id]  = sE.getXyzFeatureDict(id)            

        if (self.__localTest):
            #
            # Create additional auth test sequences with random insertions and deletions
            #
            for id in ['C','D','E','F','G']:
                refSL[id]  = sE.getRefSequenceTestWithIndexList('A')           
                authSL[id] = sE.getAuthSequenceTestWithIndexList('A')
                xyzSL[id]  = sE.getXyzSequenceTestWithIndexList('A')
                refFD[id]  = sE.getRefFeatureDict('A')
                authFD[id] = sE.getAuthFeatureDict('A')
                xyzFD[id]  = sE.getXyzFeatureDict('A')                                    
                
            idList= ['A','B','C','D','E','F','G']
            groupDict={}
            groupDict[1]=['A','C','D','E','F','G']
            groupDict[2]=['B']
            
        else:
            idList= ['A','B']
            groupDict={}
            groupDict[1]=['A']
            groupDict[2]=['B']
            

        #
        sds=SequenceDataStore(sessionObj=self.__sessionObj,verbose=self.__verbose,log=self.__lfh)
        sds.reset()
        for id in idList:
            sds.setSequence(authSL[id], id, 'auth', altId=1, version=1)
            sds.setFeature( authFD[id], id, 'auth', altId=1, version=1)            
            sds.setSequence(xyzSL[id],  id, 'xyz',  altId=1, version=1)
            sds.setFeature( xyzFD[id],  id, 'xyz',  altId=1, version=1)            
            for altId in range(1,10):
                sds.setSequence(refSL[id],  id, 'ref',  altId=altId, version=1)
                sds.setFeature( refFD[id],  id, 'ref',  altId=altId, version=1)                                    
        #
        #
        entryDict={}
        entryDict['PDB_ID']=sE.getIdCode()
        entryDict['RCSB_ID']='RCSB054485'
        #
        for k,v in groupDict.items():
            sds.setGroup(k,v)
        for k,v in entryDict.items():
            sds.setEntryDetail(k,v)

        sds.serialize()

        self.__copyExampleStructureFile(id=sE.getIdCode(),fileType="pdb")
        
        return {}

