##
# File:  RcsbPath.py
# Date:  4-June-2010
#
# Updates:
# 3-May-2011  Add file type of cifeps as an alias for eps
##
"""
Utility methods for accessing data files in current data production systems.

Stripped down from SiteInterface in SeqModule.
     
"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"

    
import sys, os.path


class RcsbPath(object):
    """
        Utilities for accessing data files in the RCSB data production system.
        
    """

    def __init__(self,topPath="/annotation",verbose=False,log=sys.stderr):
        self.__verbose = verbose
        self.__topPath = topPath
        self.__lfh=log
        #
        self.__rcsbId=None
        self.__dirPath=None
        self.__sessionPath='.'

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
            #if (self.__verbose):
            #    self.__lfh.write("+RcsbPath.__getPath() - trying path %s\n" %  pth)        
            if os.access(pth, os.R_OK):
                oPth=pth

        if (self.__verbose):
            self.__lfh.write("+RcsbPath.__getPath() - returning path %s\n" %  str(oPth))        
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
        elif (fileType == "cifeps"):
            fn=self.__rcsbId+'.cifeps'                        
        elif (fileType == "pdb"):
            fn=self.__rcsbId+'.pdb'            
        elif (fileType == "sf"):                        
            fn=self.__rcsbId+'-sf.cif'            
        else:
            pass
        
        return fn

    def getFilePath(self,fileType="cif"):
        fn=self.__getStructureFileName(fileType)
        src=os.path.join(self.__dirPath,fn)
        if (self.__exists(src)):
            return src
        else:
            return None
