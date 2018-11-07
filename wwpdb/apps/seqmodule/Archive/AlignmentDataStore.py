##
# File:    AlignmentDataStore.py
# Date:    14-Jan-2010
#
# Updates:
#   15-Feb-2010 jdw Add list of aligned sequences as separate storeage element.
#   20-Apr-2010 jdw Ported to module seqmodule.
#   27-Feb-2013 jdw  obsoleted --  functionality transfered 
##
"""
Provides a storage utilities for sequence alignment for use by the sequence editing tool.

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"

import sys, cPickle, time, os.path,  pprint
from  wwpdb.apps.seqmodule.util.Autodict import Autodict

        
            
class AlignmentDataStore(object):
    """Storage utility class for sequence alignment data.

       All sequences are stored as lists::
    
       - These objects are opaque to this class which only manages their storage.

      Objects are identified by  -

        dataType:  Supported data types are 'sequence' and 'feature'
        seqType:   Sequence types are: auth (author provided), xyz (coordinate),
                   ref (from reference sequence database).
        seqId:     Identifies an instance of a sequence or a feature object within a 
                   sequence type.  Sequence and feature instances are versioned.  This 
                   identifier serves can be a PDB ChainId or the PDBx (_struct_asym.id) 
                   for the author and coordinate sequence types.
        partId:    Distinguishes continuous regions of a sequence with unique/separately
                   specified features. This identifier is used to maintain correpondences
                   between portions of an author provided sequence and related reference
                   sequences.    This is an integer value with default value 1. 
        altId:     Distinguishes alternatives associated with a particular sequence
                   instance.   This may be used to support micro-hetereogeneity and
                   alternative reference sequences (default=1)
        version:   integer revision id where 1 is the original version.

        self.__I[dataType][seqType][seqId][partId][altId][version]=((dataId,timeStamp)
        
        self.__D[dataId] = sequence list [...] or  feature dictionary {...} 

        where dataId is a concatenated identifier= seqType_dataType_seqId_altId_version

        self.__L[seqIds,...]  sequence ids from class SequenceLabel()
        
    """
    def __init__(self,sessionObj,fileName='alignmentDataStore.pic',verbose=False,log=sys.stderr):
        self.__verbose=verbose
        self.__lfh=log
        self.__sessionObj=sessionObj
        self.__fileName=fileName        
        self.__D={}
        self.__I = Autodict()
        self.__L=[]
        #
        self.__sessionPath='.'
        self.__filePath=self.__fileName
        #
        #self.__pickleProtocol = cPickle.HIGHEST_PROTOCOL
        self.__pickleProtocol=0
        #
        self.__setup()        

    def __setup(self):
        try:
            self.__sessionPath=self.__sessionObj.getPath()
            self.__filePath = os.path.join(self.__sessionPath,self.__fileName)
            if (self.__verbose):
                self.__lfh.write("+AlignmentDataStore.__setup() - session id %s path %s filePath %s\n" %
                                 (self.__sessionObj.getId(),self.__sessionObj.getPath(),self.__filePath))

            self.deserialize()
        except:
            self.__lfh.write("+AlignmentDataStore.__setup() - Failed to open alignment store for session id %s alignment store filePath %s\n" %
                             (self.__sessionObj.getId(), self.__filePath))            

        
    def reset(self):
        """Clear internal data store.
        """
        self.__D={}
        self.__I=Autodict()
        self.__L=[]
       
    def serialize(self):
        """Store aligned sequence data using cPickle serializer.
        """
        try:
            fb=open(self.__filePath,'wb')
            #cPickle.dump(self.__E,fb,self.__pickleProtocol)
            cPickle.dump(self.__L,fb,self.__pickleProtocol)                        
            cPickle.dump(self.__I,fb,self.__pickleProtocol)            
            cPickle.dump(self.__D,fb,self.__pickleProtocol)
            fb.close()
        except:
            pass
            
    def deserialize(self):
        """Recover aligned sequence data from persistent store.
        """
        try:
            fb=open(self.__filePath,'rb')
            #self.__E=cPickle.load(fb)
            self.__L=cPickle.load(fb)
            self.__I=cPickle.load(fb)            
            self.__D=cPickle.load(fb)
            fb.close()
        except:
            pass


    def setAlignIdList(self,alignIdList):
        """Set the list of identifiers for aligned sequences.
        """
        self.__L=alignIdList

    def getAlignIdList(self):
        """Return the list of identifiers for the stored aligned sequences.
        """
        return self.__L


    def __makeId(self,dataType,seqType,seqId,altId=1,version=1):
        id="%s_%s_%s_%d_%d" % (dataType,seqType,seqId,altId,version)
        return id

    def __updateIndex(self,id, dataType, seqType,seqId,altId=1,version=1):
        lt = time.strftime("%Y %m %d %H:%M:%S", time.localtime())
        self.__I[dataType][seqType][seqId][altId][version]=(id,lt)

        
    def setSequence(self,sL,seqId,seqType,altId=1,version=1):
        """Store a sequence alignment object corresponding to the input identifiers.
        """                                        
        try:
            id=self.__makeId(dataType="sequence", seqType=seqType,seqId=seqId,altId=altId,version=version)
            self.__updateIndex(id, dataType="sequence", seqType=seqType,seqId=seqId,altId=altId,version=version)
            self.__D[id]=sL
            return True
        except:
            return False
        
    def getSequence(self,seqId,seqType,altId=1,version=1):
        """Return a sequence alignment object corresponding to the input identifiers.
        """                                
        try:
            id=self.__makeId(dataType="sequence", seqType=seqType,seqId=seqId,altId=altId,version=version)
            return self.__D[id]
        except:
            return []

    def setFeature(self,fD,seqId,seqType,altId=1,version=1):
        """Store a feature object corresponding to the input identifiers.
        """                        
        try:
            id=self.__makeId(dataType="feature",seqType=seqType,seqId=seqId,altId=altId,version=version)            
            self.__updateIndex(id,dataType="feature",seqType=seqType,seqId=seqId,altId=altId,version=version)
            self.__D[id]=fD
            return True
        except: 
            return False 

    def getFeature(self,seqId,seqType,altId=1,version=1):
        """Return a feature object corresponding to the input identifiers.
        """                
        try:
            id=self.__makeId(dataType="feature",seqType=seqType,seqId=seqId,altId=altId,version=version)
            return self.__D[id]
        except:
            return {}
        
    def getSequenceTypes(self,dataType='sequence'):
        """Return the list of sequence types.
        """        
        try:
            return(self.__I[dataType].keys())
        except:
            return []

    def getDataTypes(self):
        """Return the list of data types.
        """
        try:
            return(self.__I.keys())
        except:
            return []

    def getIds(self, dataType="sequence", seqType="ref"):
        """Return the list of sequence or feature identifiers of the input type.
        """                
        try:
            return(self.__I[dataType][seqType].keys())
        except:
            return []

    def getAlternativeIds(self,seqId, dataType="sequence", seqType="ref"):
        """Return the list of alternate identifiers for the input sequence or feature.
        """        
        try:
            return(self.__I[dataType][seqType][seqId].keys())
        except:
            return []

    def getVersionIds(self,seqId, altId=1, dataType="sequence", seqType="ref"):
        """Return the list of version identifiers for the input sequence or feature.
        """
        try:
            return(self.__I[dataType][seqType][seqId][altId].keys())
        except:
            return []
        
        
    def dump(self,ofh):
        """Formatted output of internal indices used to store sequence and features data.

           The index storage model is::
           
             self.__I[dataType][seqType][seqId][altId][version]=((dataId,timeStamp)
             
        """
        ofh.write("\nAligned Sequence Id List:\n")
        for id in self.__L:
            ofh.write("   %s\n" %  id)
        ofh.write("\nAligned Sequence Index:\n")
        for type,v0 in self.__I["sequence"].items():  
            for id,v1 in v0.items():
                for altId, v2 in v1.items():
                    for ver, ival in v2.items():
                        ofh.write("Sequence type %5s id %4s version %2s alternative %2d updated %12s sequence length %10d\n" %
                                  (type,id,ver,altId,ival[1],len(self.__D[ival[0]]) ))

        ofh.write("\nAligned Sequence Feature Index:\n")
        for type,v0 in self.__I['feature'].items():
            for id,v1 in v0.items():
                for altId,v2 in v1.items():
                    for ver, ival in v2.items():                    
                        ofh.write("Sequence type %5s id %4s version %2s alternative %2d updated %12s feature length %10d\n" %
                                  (type,id,ver,altId,ival[1],len(self.__D[ival[0]]) ))
        
        
    def dumpData(self,ofh,seqId,seqType,dataType="sequence"):
        """Formatted output of sequence or feature data.
        """
        altIds = self.__I[dataType][seqType][seqId].keys()
        for altId, vOb in self.__I[dataType][seqType][seqId].items():
            for ver, ival in vOb.items():
                id=self.__makeId(dataType=dataType,seqType=seqType,seqId=seqId,altId=altId,version=ver)            
                ofh.write("Data contents for %10s type %5s id %5s alternative %2d version %2d\n" %
                          (dataType,seqType,seqId,altId,ver))
                if self.__D.has_key(id):
                    pprint.pprint(self.__D[id],stream=ofh)
                else:
                    ofh.write(" NO DATA FOUND\n")
