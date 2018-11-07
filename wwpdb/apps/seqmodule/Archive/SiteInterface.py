##
# File:  SiteInterface.py
# Date:  4-Feb-2010
#
# Updates:
# 20-Apr-2010 jdw Ported to module seqmodule.
##
"""
Utility methods for accessing data files in current data production systems.
     
"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"

    
import sys, shutil,  os.path


class InterfaceArchiveRcsb(object):
    """
        Utilities for accessing data files in the RCSB data production system.
        
    """

    def __init__(self,sessionObj,topPath="/annotation",verbose=False,log=sys.stderr):
        self.__verbose = verbose
        self.__sessionObj=sessionObj
        self.__topPath = topPath
        self.__lfh=log
        #
        self.__rcsbId=None
        self.__dirPath=None
        self.__sessionPath='.'
        self.__setup()
        
    def __setup(self):
        try:
            self.__sessionPath=self.__sessionObj.getPath()
            if (self.__verbose):
                self.__lfh.write("+InterfaceArchiveRcsb.__setup() - session id %s path %s\n" %
                                 (self.__sessionObj.getId(),self.__sessionObj.getPath()))
        except:
            pass


    def setId(self,rcsbId):
        self.__rcsbId=str(rcsbId).lower()
        self.__dirPath=self.__getPath()
        return (self.__dirPath != None)

    def existsPath(self,rcsbId):
        return (self.__getPath() != None)
    
    def __getPath(self):
        """ Enumerate the possible archive paths and return the first match or None
        """
        oPth=None
        for sDir in ['prot','nmr','ndb','test']:
            pth=os.path.join(self.__topPath,sDir,self.__rcsbId)
            if (self.__verbose):
                self.__lfh.write("+InterfaceArchiveRcsb.__getPath() - trying path %s\n" %  pth)        
            if os.access(pth, os.R_OK):
                oPth=pth

        if (self.__verbose):
            self.__lfh.write("+InterfaceArchiveRcsb.__getPath() - returning path %s\n" %  str(oPth))        
        return oPth

    def __exists(self,pth):
        if os.access(pth, os.R_OK):
            return True
        return False
            
    def __getStructureFileName(self,fileType):
        """ 
        """
        fn=""
        if (fileType == "cif"):
            fn=self.__rcsbId+'.cif'
        elif (fileType == "eps"):
            fn=self.__rcsbId+'.cifeps'            
        elif (fileType == "pdb"):
            fn=self.__rcsbId+'.pdb'            
        elif (fileType == "sf"):                        
            fn=self.__rcsbId+'-sf.cif'            
        else:
            pass
        
        return fn

    def copyStructureFile(self,fileType="cif"):
        fn=self.__getStructureFileName(fileType)
        if (self.__verbose):
            self.__lfh.write("+InterfaceArchiveRcsb.__copyStructureFile () - copying file name  %s\n" %  fn)                
        dst=os.path.join(self.__sessionPath,fn)
        src=os.path.join(self.__dirPath,fn)
        if (self.__exists(src)):
            if (self.__verbose):
                self.__lfh.write("+InterfaceArchiveRcsb.__copyStructureFile() - src  %s dst %s\n" %  (src,dst))                        
            shutil.copyfile(src,dst)
            return dst
        else:
            return ''


class InterfaceSequenceRcsb(object):
    """
        Utilities for accessing sequence data files in the RCSB data production system.
        
    """

    def __init__(self,sessionObj, topPath="/www-rcsb/supertool/blast/rcsb",verbose=False,log=sys.stderr):
        self.__verbose = verbose
        self.__sessionObj=sessionObj
        self.__lfh=log
        self.__topPath = topPath
        self.__dirPath = None
        self.__rcsbId=None
        self.__sessionPath='.'
        self.__setup()

    def __setup(self):
        try:
            self.__sessionPath=self.__sessionObj.getPath()
            if (self.__verbose):
                self.__lfh.write("+InterfaceSequenceRcsb.__setup() - session id %s path %s\n" %
                                 (self.__sessionObj.getId(),self.__sessionObj.getPath()))
        except:
            pass

    def setId(self,rcsbId):
        self.__rcsbId=str(rcsbId).upper()
        self.__dirPath=self.__getPath()
        return (self.__dirPath != None)

    def existsPath(self,rcsbId):
        return (self.__getPath() != None)
    
    def __getPath(self):
        """ Check for the existence of data directory and return path or None
        """
        pth=os.path.join(self.__topPath,self.__rcsbId)
        if (self.__verbose):
            self.__lfh.write("+InterfaceSequenceRcsb.__getPath() - trying path %s\n" %  pth)                
        if os.access(pth, os.R_OK):
            return pth
        return None

    def __getEntityFileName(self,entityId='1'):
        fn = self.__rcsbId + '.' + str(entityId).strip() + '.info.cif'
        return fn
        
    def existsEntityFile(self,entityId='1'):
        pth=os.path.join(self.__topPath,self.__rcsbId,self.__getEntityFileName(entityId))
        if os.access(pth, os.R_OK):        
            return True
        return False
        
    def copyEntityFile(self,entityId='1'):
        fn=self.__getEntityFileName(entityId)
        dst=os.path.join(self.__sessionPath,fn)
        src=os.path.join(self.__topPath,self.__rcsbId,fn)
        if (os.access(src,os.R_OK)):
            shutil.copyfile(src,dst)
            return dst
        else:
            return ''

    def putExportFile(self,fn):
        src=os.path.join(self.__sessionPath,fn)
        dst=os.path.join(self.__topPath,self.__rcsbId,fn)
        if (os.access(src,os.R_OK)):
            shutil.copyfile(src,dst)
            return dst
        else:
            return ''
    
    

    
        
        

    
