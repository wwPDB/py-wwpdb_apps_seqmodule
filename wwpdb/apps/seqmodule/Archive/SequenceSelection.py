##
# File:  SequenceSelection.py
# Date:  02-May-2010
#
# Updates:
# 02-May-2010  jdw Moved from SequenceDataExport().
# 05-May-2010  jdw Select all coordinate sequences.
#
##
"""
Manage sequence selections for sequence summary and sequence export operations.
     
"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"

import sys, os.path, string, shutil

from wwpdb.apps.seqmodule.io.SequenceDataStore       import SequenceDataStore
from wwpdb.apps.seqmodule.util.SequenceLabel         import SequenceLabel
from wwpdb.apps.seqmodule.util.MiscUtils             import multikeysort

#
#
class SequenceSelection(object):
    """Manage sequence selections for sequence summary and sequence export operations.
       

    """
    def __init__(self,reqObj=None, exportList=[], verbose=False,log=sys.stderr):
        self.__verbose=verbose
        self.__reqObj=reqObj
        self.__lfh=log
        self.__debug=False
        self.__sessionObj=None
        self.__sessionPath='.'
        self.__sds = None
        #
        self.__groupStatusD={}
        #
        self.__setup()
        #

        
    def __setup(self):
        try:
            self.__sessionObj  = self.__reqObj.getSessionObj()
            self.__sessionPath = self.__sessionObj.getPath()
            #
            self.__sds=SequenceDataStore(reqObj=self.__reqObj,verbose=self.__verbose,log=self.__lfh)
            #
            if (self.__verbose):
                self.__lfh.write("+SequenceSelection.__setup() - session id %s path %s\n" %
                                 (self.__sessionObj.getId(),self.__sessionObj.getPath()))
        except:
            self.__lfh.write("+SequenceSelection.__setup() failed for session id %s\n" % (self.__sessionObj.getId()))


    def groupStatus(self):
        return self.__groupStatusD
    

    def makeDefaultSelection(self):
        """ Make a guess at the sequence selection and return the default selection of sequence identifiers.

            Returns:
        """

        if (self.__verbose):
            self.__lfh.write("+SequenceSelection.__makeDefaultSelection() with sessionId %s\n" % (self.__sessionObj.getId()))
        #
        defSelectList=[]
        self.__groupStatusD={}
        #
        # Get the groups list -
        #
        gIdList=self.__sds.getGroupIds()
        if len(gIdList) < 1:
            return defSelectList

        gIdList.sort(lambda x, y: int(x) - int(y))
        seqLabel=SequenceLabel()
        #
        for gId in gIdList:
            self.__groupStatusD[gId]=True
            #
            # Get seqIds in group -
            #
            seqIdList=self.__sds.getGroup(gId)
            if (self.__verbose):
                self.__lfh.write("+SequenceSelection.__makeDefaultSelection() group %s sequence list %r\n" % (gId,seqIdList))
            
            if len(seqIdList) < 1:
                if (self.__verbose):
                    self.__lfh.write("+SequenceSelection.__makeDefaultSelection() sessionId %s group %s is empty\n" %  (self.__sObj.getId(),gId))
                continue
            #
            seqId0=seqIdList[0]
            
            #
            # Author sequence - latest version 
            #
            verList=self.__sds.getVersionIds(seqId0, altId=1, dataType="sequence", seqType='auth')
            if len(verList) > 0:
                altId=1
                ver=verList[0]
                seqLabel.set( seqType='auth',seqInstId=seqId0, seqAltId=altId, seqVersion=ver)
                idAuthSeq=seqLabel.pack()
                defSelectList.append(idAuthSeq)

            #
            # Get properties of the author sequence --
            #
            seqAuthFD =self.__sds.getFeature(seqId=seqId0,seqType='auth',altId=altId,version=ver)
            seqAuthIdx=self.__sds.getSequence(seqId=seqId0,seqType='auth',altId=altId,version=ver)            
            authPolyTypeCode=seqAuthFD['POLYMER_TYPE']
            authTaxId=seqAuthFD['SOURCE_TAXID']
            authSeqLen=len(seqAuthIdx)

            if authTaxId is None:
                haveAuthTaxId=False                
            if ( (authTaxId is not None) and  ( len(authTaxId) > 0 ) and (authTaxId != '?') and (authTaxId != '.')):
                haveAuthTaxId=True
            else:
                haveAuthTaxId=False                

            if self.__verbose:
                self.__lfh.write("SequenceSelection.makeDefaultSelection() Author feature dictionary for sequence group %s ID  %s\n" % (gId,idAuthSeq))
                self.__lfh.write("SequenceSelection.makeDefaultSelection() Author polymer type code %s\n" % authPolyTypeCode)
                self.__lfh.write("SequenceSelection.makeDefaultSelection() Author TaxId %s have TaxId %r\n" % (authTaxId,haveAuthTaxId))
                self.__lfh.write("SequenceSelection.makeDefaultSelection() Author sequence length %d\n" % authSeqLen)
            #
            
            
            #
            # Coordinate sequences - latest version of each chain
            #
            tL=[]
            altId=1
            for seqId in seqIdList:
                verList=self.__sds.getVersionIds(seqId, altId=altId, dataType="sequence", seqType='xyz')
                for ver in verList:
                    seqLabel.set( seqType='xyz',seqInstId=seqId, seqAltId=altId, seqVersion=ver)
                    idXyzSeq=seqLabel.pack()
                    seqXyzFD=self.__sds.getFeature(seqId, 'xyz',altId=altId,version=ver)
                    #
                    lengthMatchXyz =seqXyzFD['MATCH_LENGTH']
                    alignLengthXyz =seqXyzFD['ALIGN_LENGTH'] 
                    identityXyz    =seqXyzFD['AUTH_XYZ_SEQ_SIM']
                    identityXyzGap =seqXyzFD['AUTH_XYZ_SEQ_SIM_WITH_GAPS']

                    dd={}
                    dd['ALIGN_LENGTH']               =seqXyzFD['ALIGN_LENGTH']
                    dd['MATCH_LENGTH']               =seqXyzFD['MATCH_LENGTH']                    
                    dd['AUTH_XYZ_SEQ_SIM']           =seqXyzFD['AUTH_XYZ_SEQ_SIM']
                    dd['AUTH_XYZ_SEQ_SIM_WITH_GAPS'] =seqXyzFD['AUTH_XYZ_SEQ_SIM_WITH_GAPS']
                    dd['ID'] = idXyzSeq
                    tL.append(dd)
                    defSelectList.append(idXyzSeq)                    


                xyzList  = multikeysort(tL, ['-AUTH_XYZ_SEQ_SIM_WITH_GAPS', '-MATCH_LENGTH', '-ALIGN_LENGTH'])

                if self.__verbose:
                    self.__lfh.write("SequenceSelection.makeDefaultSelection()  Ordered xyz list for GROUP ID %s\n" % gId)                                     
                    for xyz  in xyzList:
                        self.__lfh.write("SequenceSelection.makeDefaultSelection()  ID %s %r\n" %  (xyz['ID'],xyz.items()))

                # The top of xyzList is the best match  -
                #
                if len(xyzList) > 0:
                    xyzIdentity=xyzList[0]['AUTH_XYZ_SEQ_SIM_WITH_GAPS']
                    id=xyzList[0]['ID']
                    #defSelectList.append(id)
                    if self.__verbose:                
                        self.__lfh.write("SequenceSelection.makeDefaultSelection() xyz sequence assigned for chain %s group %s\n" % (id,gId))                    
                    
                    if ( abs(xyzIdentity - 1.000) < 0.001):
                        if self.__verbose:                
                            self.__lfh.write("SequenceSelection.makeDefaultSelection() Assigning matching xyz sequence (%8.4f) for %s sequence group %s\n" % (xyzIdentity,id,gId))
                    else:
                        self.__groupStatusD[gId]=False

                else:
                    self.__groupStatusD[gId]=False                    
                    if self.__verbose:                
                        self.__lfh.write("SequenceSelection.makeDefaultSelection() Missing xyz sequence data for chain %s group %s (NOT REFERENCED)\n" % (seqId,gId))

            ##
            ##




            # --------------------------------------------------------------------------------------------------------------------------
            # Special cases -- for self reference -- 
            #
            # No author TaxId
            if (len(authTaxId) < 1  or authTaxId == '0'):
                defSelectList.append('selfref_'+str(gId))
                if self.__verbose:                
                    self.__lfh.write("SequenceSelection.makeDefaultSelection() Self reference (no taxid) sequence group %s\n" % gId)
                #
                self.__groupStatusD[gId]=False
                continue

            #
            # Nucleic acid polymers 
            if ((authPolyTypeCode in ['RNA', 'DNA']) and (authSeqLen < 50)) :
                defSelectList.append('selfref_'+str(gId))
                if self.__verbose:                
                    self.__lfh.write("SequenceSelection.makeDefaultSelection() Self reference (NA polymer) sequence group %s\n" % gId)
                    continue                
            
            #
            # Reference sequences - 
            #
            
            # List of reference sequences for this group only for the leading sequence - 
            #
            tL=[]            
            altIdList=self.__sds.getAlternativeIds(seqId0, dataType="sequence", seqType="ref")
            for altId in altIdList:
                #
                verList=self.__sds.getVersionIds(seqId=seqId0, altId=altId, dataType="sequence", seqType='ref')
                for ver in verList:                
                    #
                    seqRef=self.__sds.getSequence(seqId=seqId0,seqType='ref',altId=altId,version=ver)
                    seqLabel.set( seqType='ref',seqInstId=seqId0, seqAltId=altId, seqVersion=ver)
                    idRefSeq=seqLabel.pack()
                    seqRefFD=self.__sds.getFeature(seqId0,'ref',altId=altId,version=ver)
                    lengthMatchRef =seqRefFD['MATCH_LENGTH']
                    identityRef    =seqRefFD['AUTH_REF_SEQ_SIM']
                    identityRefGap =seqRefFD['AUTH_REF_SEQ_SIM_WITH_GAPS']
                    refTaxId =seqRefFD['SOURCE_TAXID']
                    refDbName  =seqRefFD['DB_NAME']
                    refDbAcc   =seqRefFD['DB_ACCESSION']                                                            
    
                    if (haveAuthTaxId and (str(authTaxId) == str(refTaxId))):
                        if self.__verbose:                        
                            self.__lfh.write("SequenceSelection.makeDefaultSelection()  %s %s %s %s  match length %d identity %8.4f identity w/gap %8.4f \n" %
                                             (idRefSeq, refTaxId,refDbName, refDbAcc,lengthMatchRef,identityRef,identityRefGap))
                        dd={}
                        dd['MATCH_LENGTH']               =seqRefFD['MATCH_LENGTH']
                        dd['AUTH_REF_SEQ_SIM']           =seqRefFD['AUTH_REF_SEQ_SIM']
                        dd['AUTH_REF_SEQ_SIM_WITH_GAPS'] =seqRefFD['AUTH_REF_SEQ_SIM_WITH_GAPS']
                        dd['DB_NAME']                    =seqRefFD['DB_NAME']
                        dd['ID'] = idRefSeq
                        tL.append(dd)
                    else:
                        if self.__verbose:                        
                            self.__lfh.write("SequenceSelection.makeDefaultSelection() (without authTaxid)  %s %s %s %s  match length %d identity %8.4f identity w/gap %8.4f \n" %
                                             (idRefSeq, refTaxId,refDbName, refDbAcc,lengthMatchRef,identityRef,identityRefGap))
                        dd={}
                        dd['MATCH_LENGTH']               =seqRefFD['MATCH_LENGTH']
                        dd['AUTH_REF_SEQ_SIM']           =seqRefFD['AUTH_REF_SEQ_SIM']
                        dd['AUTH_REF_SEQ_SIM_WITH_GAPS'] =seqRefFD['AUTH_REF_SEQ_SIM_WITH_GAPS']
                        dd['DB_NAME']                    =seqRefFD['DB_NAME']
                        dd['ID'] = idRefSeq
                        tL.append(dd)
                        

            rfList  = multikeysort(tL, ['-AUTH_REF_SEQ_SIM_WITH_GAPS', '-MATCH_LENGTH', '-DB_NAME'])

            if self.__verbose:
                self.__lfh.write("SequenceSelection.makeDefaultSelection()  Ordered reference list for GROUP ID %s\n" % gId)                                     
                for rf in rfList:
                    self.__lfh.write("SequenceSelection.makeDefaultSelection()  ID %s %r\n" %  (rf['ID'],rf.items()))

            # The top of rfList is the best match  -
            #
            if len(rfList) > 0:
                rfIdentity=rfList[0]['AUTH_REF_SEQ_SIM_WITH_GAPS']
                id=rfList[0]['ID']
                if (rfIdentity < 0.8):
                    if self.__verbose:                
                        self.__lfh.write("SequenceSelection.makeDefaultSelection() Self reference ( %r match identity < .80) for sequence group %s\n"
                                         % (rfIdentity,gId))
                    defSelectList.append('selfref_'+str(gId))                    
                else:
                    defSelectList.append(id)

                if ( abs(rfIdentity - 1.000) > 0.001):
                    self.__groupStatusD[gId]=False                    
            else:
                defSelectList.append('selfref_'+str(gId))
                self.__groupStatusD[gId]=False                
                if self.__verbose:                
                    self.__lfh.write("SequenceSelection.makeDefaultSelection() Self reference (no matching ref seq) for sequence group %s\n" % gId)

        if (self.__verbose):
            self.__lfh.write("SequenceSelection.makeDefaultSelection() SUMMARY ---- SUMMARY ---- SUMMARY\n")                            
            for ss in defSelectList:
                self.__lfh.write("SequenceSelection.makeDefaultSelection() %s\n" % ss)                
            for gId in gIdList:
                if not self.__groupStatusD[gId]:
                    self.__lfh.write("SequenceSelection.makeDefaultSelection() Default assignment status is false for %s\n" % gId)
                

        return defSelectList


    def makeDefaultSelectionOrg(self):
        """ Make a guess at the sequence selection and return the default selection of sequence identifiers.
        """

        if (self.__verbose):
            self.__lfh.write("+SequenceSelection.__makeDefaultSelection() sessionId %s\n" %  (self.__sessionObj.getId()))
        #
        defSelectList=[]        
        # Get the groups list -
        #

        #
        gIdList=self.__sds.getGroupIds()
        if len(gIdList) < 1:
            return defSelectList

        gIdList.sort(lambda x, y: int(x) - int(y))
        seqLabel=SequenceLabel()
        #
        for gId in gIdList:
            #
            # Get seqIds in group -
            #
            seqIdList=self.__sds.getGroup(gId)
            if (self.__verbose):
                self.__lfh.write("+SequenceSelection.__makeDefaultSelection() group %s sequence list %r\n" % (gId,seqIdList))
            
            if len(seqIdList) < 1:
                if (self.__verbose):
                    self.__lfh.write("+SequenceSelection.__makeDefaultSelection() sessionId %s group %s is empty\n" %  (self.__sObj.getId(),gId))
                continue
            #
            seqId0=seqIdList[0]
            #
            # Author sequence - latest version 
            #
            verList=self.__sds.getVersionIds(seqId0, altId=1, dataType="sequence", seqType='auth')
            if len(verList) > 0:
                altId=1
                ver=verList[0]
                seqLabel.set( seqType='auth',seqInstId=seqId0, seqAltId=altId, seqVersion=ver)
                idAuthSeq=seqLabel.pack()
                defSelectList.append(idAuthSeq)

            #
            # Reference sequence - 
            #
            altId=1
            verList=self.__sds.getVersionIds(seqId0, altId=altId, dataType="sequence", seqType='ref')
            if len(verList) > 0:
                ver=verList[0]
                seqLabel.set( seqType='ref',seqInstId=seqId0, seqAltId=altId, seqVersion=ver)
                idRefSeq=seqLabel.pack()
                defSelectList.append(idRefSeq)
                
            #
            # Coordinate sequences -
            #
            altId=1
            for seqId in seqIdList:
                verList=self.__sds.getVersionIds(seqId, altId=altId, dataType="sequence", seqType='xyz')
                if len(verList) > 0:
                    ver=verList[0]                
                    seqLabel.set( seqType='xyz',seqInstId=seqId, seqAltId=altId, seqVersion=ver)
                    idXyzSeq=seqLabel.pack()
                    defSelectList.append(idXyzSeq)                            
                    
        return defSelectList
            
    
