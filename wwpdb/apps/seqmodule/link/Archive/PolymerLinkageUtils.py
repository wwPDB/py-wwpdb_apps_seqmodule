##
# File:  PolymerLinkageUtils.py
# Date:  21-Feb-2013
#
# Updates:
#   23-Feb-2013 jdw  Use ioadapter compatible IO
#    4-Apr-2013 jdw  Add method to perform more precise filtering including the cardinal sequence index.
#    1-Dec-2013 jdw  Update the linkage instance dictionary with the residue order (first|second) of residue distance outliers.
#                    Add insertion codes to link index
#    2-Dec-2013 jdw  Simplify the tuple returned by in the link instance dictionary -
##
"""
Calculate and assemble polymer linkage details from the model coordinate data.
     
"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.09"

import sys, os, string, traceback

from wwpdb.apps.seqmodule.io.PdbxIoUtils    import PolymerLinkageIo,PdbxFileIo

from wwpdb.utils.dp.RcsbDpUtility                     import RcsbDpUtility

class PolymerLinkageUtils(object):
    """
    Assemble polymer linkage details from the model coordinate data.
    """
    def __init__(self,reqObj=None,verbose=False,log=sys.stderr):
        self.__reqObj=reqObj
        self.__verbose=verbose        
        self.__lfh=log
        #
        self.__siteId=self.__reqObj.getValue("WWPDB_SITE_ID")
        self.__sObj=self.__reqObj.getSessionObj()
        self.__sessionPath=self.__sObj.getPath()
        #
        self.__cleanUp=True
        #

    def calc(self,modelFilePath,polyLinkPath="polymer-link-distances.cif"):
        if (self.__calcPolymerLinkages(modelFilePath,polyLinkPath)):
            return self.__readPolymerLinkageDistances(polyLinkPath)
        else:
            return []

    def __calcPolymerLinkages(self,modelFilePath,polyLinkPath):
        try:
            dp = RcsbDpUtility(tmpPath=self.__sessionPath, siteId=self.__siteId,verbose=self.__verbose,log=self.__lfh)
            dp.imp(modelFilePath)
            dp.op("annot-poly-link-dist")
            dp.exp(polyLinkPath)
            if (self.__cleanUp): dp.cleanup()                        
            if (self.__verbose):
                self.__lfh.write("+PolymerLinkageUtils.__calcPolymerLinkages() - Polymer link distance file copied to: %s\n" % polyLinkPath)
            return True
        except:
            if (self.__verbose):
                self.__lfh.write("+PolymerLinkageUtils.__calcPolymerLinkages() - Polymer link distance calculation failed for file %s\n" % modelFilePath)
                traceback.print_exc(file=self.__lfh)
            return False

    def makeLinkDict(self,polymerLinkDistList=None,minD=1.2,maxD=1.75):
        """  Convert linkList to a dictionary of link outliers with structure  - 

             d[authAsymId,compId,authSeqId+insCode] = (compId,authSeqId+insCode,float(dist),seqId,'begin|end',longFlag)
        """
        linkInstD={}
        for lD in polymerLinkDistList:
            dist=lD['dist']
            if float(dist) > maxD or float(dist) < minD :
                longFlag = float(dist) > maxD 
                seqId1=lD['auth_seq_id_1']
                compId1=lD['auth_comp_id_1']
                authAsymId1=lD['auth_asym_id_1']
                lblSeqId1=lD['label_seq_id_1']
                insCode1=lD['PDB_ins_code_1']

                if len(insCode1)>0 and insCode1 not in ['.','?']: 
                    seqId1 = seqId1+insCode1

                linkInstD[(authAsymId1,compId1,seqId1)]= (compId1,seqId1,float(dist),lblSeqId1,'begin',longFlag)

                seqId2=lD['auth_seq_id_2']
                compId2=lD['auth_comp_id_2']
                authAsymId2=lD['auth_asym_id_2']
                lblSeqId2=lD['label_seq_id_2']
                insCode2=lD['PDB_ins_code_2']
                if len(insCode2)>0 and insCode2 not in ['.','?']: 
                    seqId2 = seqId2+insCode2
                linkInstD[(authAsymId2,compId2,seqId2)]= (compId2,seqId2,float(dist),lblSeqId2,'end',longFlag)

        if self.__verbose:
            self.__lfh.write('+PolymerLinkageUtils.makeLinkDict() link dictionary %r\n' % linkInstD.items())
            
        return linkInstD



    def __readPolymerLinkageDistances(self,polyLinkPath):
        """  Read the file of polymer linkage distaces for the current model.
        """
        #
        polymerLinkDistList=[]
        try:
            if os.access(polyLinkPath,os.R_OK):
                c0 = PdbxFileIo(verbose=self.__verbose,log=self.__lfh).getContainer(polyLinkPath)
                pl = PolymerLinkageIo(dataContainer=c0,verbose=self.__verbose,log=self.__lfh)
                polymerLinkDistList=pl.getPolymerLinkDistances()
        except:
            traceback.print_exc(file=self.__lfh)

        #
        return polymerLinkDistList

    def makeLinkDictWithFilter(self,polymerLinkDistList=None,minD=1.2,maxD=1.75):
        """  Convert linkList to a link dictionary with keys - 

        Select cases   -    if residues are adjacent and dist < minD && dist > maxD 

                         or if residues are not adjacent and dist < maxD

        Adjacency is approximately determined by author residue numbering.
        
             d[authAsymId,compId,authSeqId+insCode] = (compId,authSeqId+insCode,float(dist),seqId,'begin|end',longFlag)

        """
        linkInstD={}
        for lD in polymerLinkDistList:
            dist=lD['dist']
            
            seqId1=lD['auth_seq_id_1']
            lblSeqId1=lD['label_seq_id_1']
            compId1=lD['auth_comp_id_1']
            authAsymId1=lD['auth_asym_id_1']
            insCode1=lD['PDB_ins_code_1']

            if len(insCode1)>0 and insCode1 not in ['.','?']: 
                seqId1 = seqId1+insCode1
            
            seqId2=lD['auth_seq_id_2']
            lblSeqId2=lD['label_seq_id_2']
            compId2=lD['auth_comp_id_2']
            authAsymId2=lD['auth_asym_id_2']
            insCode2=lD['PDB_ins_code_2']
            if len(insCode2)>0 and insCode2 not in ['.','?']: 
                seqId2 = seqId2+insCode2

            if self.__verbose:
                if ((float(dist) > maxD) or (float(dist) < minD)) :
                    self.__lfh.write('+PolymerLinkageUtils.makeLinkWithFilterDict() seqId1 %r seqId2 %r lblSeqId1 %r lblSeqId2 %r  dist %r\n' %
                                     (seqId1,seqId2,lblSeqId1,lblSeqId2,dist))
            #
            flag=False
            if (  abs(int(lblSeqId1) - int(lblSeqId2)) > 1 ):
                # non-adjacent residues
                if float(dist) < maxD:
                    flag=True
                    
            else:
                # adjacent residues
                if ((float(dist) > maxD) or (float(dist) < minD)) :
                    flag=True

            longFlag=(float(dist) > maxD)

            if (flag):
                linkInstD[(authAsymId1,compId1,seqId1)]= (compId1,seqId1,float(dist),lblSeqId1,'begin',longFlag)
                linkInstD[(authAsymId2,compId2,seqId2)]= (compId2,seqId2,float(dist),lblSeqId2,'end',longFlag)



        if self.__verbose:
            self.__lfh.write('+PolymerLinkageUtils.makeLinkWithFilterDict() filtered link dictionary %r\n' % linkInstD.items())


            
        return linkInstD
