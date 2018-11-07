##
# File:    DataFileAdapter.py
# Date:    16-July-2012
#
# Updates:
#
#  4-Jan-2013 -jdw  add filters to remove derived categories after format conversion
# 12-Feb-2013  jdw  add option flags for stipping derived categories -- 
# 25-Feb-2013  jdw  add cif2pdb translation.
##
"""
Encapsulate data model format type conversions.

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"


import sys, os.path, traceback
from wwpdb.utils.rcsb.RcsbDpUtility       import RcsbDpUtility

class DataFileAdapter(object):
    """  Manage data model format type conversions.
        
    """
    def __init__(self,reqObj, verbose=False,log=sys.stderr):
        self.__reqObj=reqObj
        self.__verbose=verbose
        self.__lfh=log
        self.__debug=True
        self.__siteId=self.__reqObj.getValue("WWPDB_SITE_ID")
        self.__sObj  =self.__reqObj.getSessionObj()
        self.__sessionId=self.__sObj.getId()
        self.__sessionPath=self.__sObj.getPath()

    def rcsb2Pdbx(self,inpPath,outPath,stripFlag=False): 
        """  RCSB CIF -> PDBx conversion  (Using the smaller application in the annotation package)
        """
        try:
            dp = RcsbDpUtility(tmpPath=self.__sessionPath, siteId=self.__siteId, verbose=self.__verbose,log=self.__lfh)
            dp.imp(inpPath)
            if stripFlag:
                dp.op("annot-rcsb2pdbx-strip")
            else:
                dp.op("annot-rcsb2pdbx")
                
            dp.expLog("annot-rcsb2pdbx.log")
            dp.exp(outPath)
            if (not self.__debug):
                dp.cleanup()
        except:
            traceback.print_exc(file=self.__lfh)
            return False
        #
        return True        

    def rcsbEps2Pdbx(self,inpPath,outPath,stripFlag=False): 
        """  RCSB CIFEPS -> PDBx conversion (This still requires using the full maxit application)
        """
        try:
            dp = RcsbDpUtility(tmpPath=self.__sessionPath, siteId=self.__siteId, verbose=self.__verbose,log=self.__lfh)
            dp.imp(inpPath)
            if stripFlag:
                dp.op("annot-rcsbeps2pdbx-strip")
            else:
                dp.op("annot-cif2cif")                

            dp.expLog("annot-rcsbeps2pdbx.log")
            dp.exp(outPath)
            if (not self.__debug):            
                dp.cleanup()
        except:
            traceback.print_exc(file=self.__lfh)
            return False
        #
        return True                        

    def cif2Pdb(self,inpPath,outPath): 
        """   CIF -> PDB conversion  (Using the smaller application in the annotation package)
        """
        try:
            dp = RcsbDpUtility(tmpPath=self.__sessionPath, siteId=self.__siteId, verbose=self.__verbose,log=self.__lfh)
            dp.imp(inpPath)
            dp.op("annot-cif2pdb")
            dp.expLog("annot-cif2pdb.log")
            dp.exp(outPath)
            if (not self.__debug):
                dp.cleanup()
        except:
            traceback.print_exc(file=self.__lfh)
            return False
        #
        return True        
    

