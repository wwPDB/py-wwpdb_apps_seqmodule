##
# File:    ReferenceSequenceDataUpdate.py
# Date:    27-Mar-2013
# Updates:
##
"""
Utilities for adding out-of-band sequence and feature data to the sequence data store.

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.08"

import sys,traceback, time, os


from wwpdb.apps.seqmodule.io.SequenceDataStore              import SequenceDataStore
from wwpdb.apps.seqmodule.io.ReferenceSequenceUtils         import ReferenceSequenceUtils
from wwpdb.apps.seqmodule.io.PdbxIoUtils                    import  PdbxFileIo,ReferenceSequenceIo
from wwpdb.apps.seqmodule.util.SequenceLabel                import SequenceLabel, SequenceFeature
from wwpdb.apps.seqmodule.util.SequenceAssign               import SequenceAssignArchive,SequenceAssignDepositor,ReferenceSequence
from wwpdb.apps.seqmodule.util.SequenceReferenceData        import SequenceReferenceData
from wwpdb.apps.seqmodule.align.AlignmentStatistics         import AlignmentStatistics
from wwpdb.io.locator.PathInfo                              import PathInfo
from wwpdb.apps.seqmodule.control.SequenceDataAssemble      import SequenceDataAssemble


class ReferenceSequenceDataUpdate(object):
    """ Utilities for adding out-of-band sequence and feature data to the sequence data store.

    """
    def __init__(self,reqObj=None,verbose=False,log=sys.stderr):
        self.__verbose=verbose
        self.__debug=False
        self.__reqObj=reqObj
        self.__lfh=log
        #
        self.__srd=SequenceReferenceData(verbose=self.__verbose,log=self.__lfh)
        #
        self.__setup()

    def __setup(self):
        try:
            self.__siteId      = self.__reqObj.getValue("WWPDB_SITE_ID")
            self.__identifier  = self.__reqObj.getValue("identifier")
            self.__sessionId   = self.__reqObj.getSessionId()        
            self.__sessionObj  = self.__reqObj.getSessionObj()
            self.__sessionPath = self.__sessionObj.getPath()
            self.__pI=PathInfo(siteId=self.__siteId,sessionPath=self.__sessionPath,verbose=self.__verbose,log=self.__lfh)
            self.__sds=SequenceDataStore(reqObj=self.__reqObj,verbose=self.__verbose,log=self.__lfh)        
            #
        except:
            if (self.__verbose):
                self.__lfh.write("+ReferenceSequenceDataUpdate.__setup() sessionId %s failed\n" % (self.__sessionObj.getId()))                                         
                traceback.print_exc(file=self.__lfh)


    def doUpdate(self,entityId,maxRefAlign=100):
        entityD=self.__getEntityDetails(entityId=entityId)
        entityD['ENTRY_ID']=self.__identifier

        if (self.__verbose):
            self.__lfh.write("+ReferenceSequenceDataUpdate.doUpdate() begins for entity %r\n" % entityId)

        if (self.__debug):
            self.__lfh.write("+ReferenceSequenceDataUpdate.doUpdate() Entity %r Entity dictionary = %r\n" % (entityId,entityD.items()))
        #
        # JDW choose very short sequence limits for a deliberate search within the UI.
        self.__doEntityReferenceSearch(dataSetId=self.__identifier,entityD=entityD,minSearchSequenceLengthAA=5,minSearchSequenceLengthNA=15)
        status,eRefL=self.__readReferenceSearchResults(dataSetId=self.__identifier,entityId=entityId,fileSource='session',
                                                       wfInstanceId=None,sessionCachePath=None)
        if (self.__verbose):
            self.__lfh.write("+ReferenceSequenceDataUpdate.doUpdate() Reference search returns status %r content reference dictionary length = %r\n" % (status,len(eRefL) ))
        self.__reloadReferenceSequences(entityD, eRefL)
        return True

    def doUpdateSelections(self,entityId,maxRefAlign=100):
        if (self.__verbose):
            self.__lfh.write("+ReferenceSequenceDataUpdate.doUpdateSelections() begins for entity %r\n" % entityId)
        #
        # Reset the selection list for this entity -- 
        #
        self.__resetEntitySequenceSelections(entityId)        
        #
        #
        alstat=AlignmentStatistics(reqObj=self.__reqObj,maxRefAlign=maxRefAlign,verbose=self.__verbose,log=self.__lfh)
        alstat.doUpdate()
        return True


    def __resetEntitySequenceSelections(self,entityId):
        selectIdList=self.__sds.getSelectedIds()
        instIdList=self.__sds.getGroup(entityId)
        #
        if (self.__verbose):
            self.__lfh.write("\n+ReferenceSequenceDataUpdate.__resetEntitySequenceSelections()  entity group %s instance list %r input selected sequence list %r\n" % 
                             (entityId,instIdList,selectIdList))
        oL=[]
        sLab=SequenceLabel()
        for seqId in selectIdList:
            sLab.unpack(seqId)
            seqType=sLab.getSequenceType()
            seqInstId=sLab.getSequenceInstId()
            if seqType in ['ref','auth'] and seqInstId == entityId:
                oL.append(seqId)
            elif seqType in ['xyz'] and seqInstId in instIdList:
                oL.append(seqId)
        #
        fL = [myId for myId in selectIdList if myId not in oL]
        if (self.__verbose):
            self.__lfh.write("+ReferenceSequenceDataUpdate.__resetEntitySequenceSelections() sequences filtered for entity %r : %r\n" % 
                             (entityId,fL))
        #
        sda=SequenceDataAssemble(reqObj=self.__reqObj,verbose=self.__verbose,log=self.__lfh)
        eL=sda.makeEntityDefaultSelection(entityId, sds=self.__sds)
        fL.extend(eL)
        if (self.__verbose):
            self.__lfh.write("+ReferenceSequenceDataUpdate.resetEntitySequenceSelections() new selection for entity %r : %r\n" % 
                             (entityId,eL))
            self.__lfh.write("+ReferenceSequenceDataUpdate.resetEntitySequenceSelections() saving new full selection list: %r\n" %  fL)

        self.__sds.setSelectedIds(idList=fL)
        self.__sds.serialize()
        #
        return True
        

    def dump(self):
        """  Output the depositor and archive sequence reference assignments.  
        """
        entityIdList=self.__sds.getGroupIds()

        seqAssignD=self.__sds.getReferenceAssignments()
        sA=SequenceAssignArchive(verbose=self.__verbose,log=self.__lfh)
        sA.set(seqAssignD)
        sA.printIt(log=self.__lfh)

        for entityId in entityIdList:
            nRef=sA.getReferenceCount(entityId=entityId)
            self.__lfh.write("+ReferenceSequenceDataUpdate.dump() entityId %s archive assignment count %d\n" % (entityId,nRef))
            if nRef > 0:
                refL=sA.getReferenceList(entityId=entityId)
                for ref in refL:
                    ref.printIt(self.__lfh)
        #
        depSeqAssignD=self.__sds.getDepositorReferenceAssignments()
        sADep=SequenceAssignDepositor(verbose=self.__verbose,log=self.__lfh)
        sADep.set(depSeqAssignD)
        sADep.printIt(log=self.__lfh)

        for entityId in entityIdList:
            nRef=sADep.getReferenceCount(entityId=entityId)
            self.__lfh.write("+ReferenceSequenceDataUpdate.dump() entityId %s depositor assignment count %d\n" % (entityId,nRef))
            if nRef > 0:
                refL=sADep.getReferenceList(entityId=entityId)
                for ref in refL:
                    ref.printIt(self.__lfh)


    def __getEntityDetails(self, entityId):
        """  Get entity feature data from the stored author sequence entity sequence and feature data.

             Returns a dictionary of features for the more recent sequence/feature versions.
        """
        entityD={}
        #
        if (self.__verbose):
            self.__lfh.write("+ReferenceSequenceDataUpdate.__getAuthFeatures() entityId %s\n" % entityId)
        seqFeature=SequenceFeature()        
        
        # get the author sequence identifier
        seqIds=self.__sds.getGroup(groupId=entityId)
        
        if len(seqIds)==0:
            return entityD
        #
        # JDW ## CHANGE 
        #seqId0=seqIds[0]
        seqId0=entityId
        altId=1
        partId=1
        verList=self.__sds.getVersionIds(seqId0, partId=partId, altId=altId, dataType="sequence", seqType='auth')
        if len(verList) == 0:
            return entityD
            
        sfObj=self.__sds.getFeatureObj(seqId=seqId0,seqType='auth',partId=partId,altId=altId, version=verList[0])
        entityD['POLYMER_TYPE']=sfObj.getPolymerType()
        entityD['POLYMER_LINKING_TYPE']=sfObj.getPolymerLinkingType()
        entityD['ENTITY_DESCRIPTION']=sfObj.getEntityDescription()
        #entityD['ENTITY_NAME']=sfObj.getEntityName()

        sTupL=self.__sds.getSequence(seqId=seqId0,seqType='auth',partId=partId,altId=altId, version=verList[0])
        r1L=[]
        for sTup in sTupL:
            r1L.append(self.__srd.cnv3To1(sTup[0]))

        entityD['SEQ_ENTITY_1']=''.join(r1L)
        entityD['SEQ_ENTITY_1_CAN']=entityD['SEQ_ENTITY_1']
        entityD['ENTITY_ID']=entityId
        entityD['PART_LIST']=self.__getPartList(entityId=entityId)
        entityD['IDENTIFIER']=seqId0

        return entityD

    def __getPartList(self,entityId):
        """ Return a list of dictionaries of entity part details.

        Return dictionary has the following keys - 
        
        pD['SEQ_NUM_BEG']
        pD['SEQ_NUM_END']
        pD['SOURCE_TAXID']
        pD['SEQ_PART_TYPE']
        pD['SEQ_PART_ID']

        """
        pL=[]
        seqIds=self.__sds.getGroup(groupId=entityId)
        if len(seqIds)==0:
            return pL

        seqFeature=SequenceFeature()
        ## JDW CHANGE
        #seqId0=seqIds[0]
        seqId0=entityId
        
        partIdList=self.__sds.getPartIds(seqId0, dataType="sequence", seqType="auth")
        
        for partId in partIdList:
            pD={}
            vL=self.__sds.getVersionIds(seqId0,partId=partId, altId=1, dataType="sequence", seqType="auth")
            if len(vL)>0:
                sfObj=self.__sds.getFeatureObj(seqId0,seqType="auth",partId=partId,altId=1,version=vL[0])
                pId,pSeqBegin,pSeqEnd,pType=sfObj.getAuthPartDetails()
                taxId = sfObj.getSourceTaxId()
                pD['SEQ_NUM_BEG']=pSeqBegin
                pD['SEQ_NUM_END']=pSeqEnd
                pD['SOURCE_TAXID']=taxId
                pD['SEQ_PART_TYPE']=pType
                pD['SEQ_PART_ID']=pId         
                #
                pD['SOURCE_NAME']=sfObj.getSourceOrganism()
                pD['SOURCE_STRAIN']=sfObj.getSourceStrain()
                pL.append(pD)
        #
        return pL


    def __doEntityReferenceSearch(self,dataSetId,entityD,minSearchSequenceLengthAA=12,minSearchSequenceLengthNA=50):
        """  Perform the reference sequence database search using the input entity dictionary.

             Store matching results in the local session directory.
        """
        try:
            startTime=time.clock()        
            #
            entityId=entityD['ENTITY_ID']
            rsu=ReferenceSequenceUtils(siteId=self.__siteId,verbose=self.__verbose,log=self.__lfh)
            rsio=ReferenceSequenceIo(verbose=self.__verbose,log=self.__lfh)

            polyTypeCode=entityD['POLYMER_TYPE']
            seqLen=len(entityD['SEQ_ENTITY_1_CAN'])
            #
            skip  = ( (polyTypeCode in ['SAC']) or 
                      ((polyTypeCode in ['DNA','RNA','XNA']) and  (seqLen < minSearchSequenceLengthNA)) or 
                      ((polyTypeCode in ['AA']) and  (seqLen < minSearchSequenceLengthAA)) )
            if (self.__verbose):
                self.__lfh.write("+ReferenceSequenceDataUpdate.__doReferenceSearch() search for entity id %s type %s length %d skip status %r\n" % 
                                 (entityId,polyTypeCode,seqLen,skip))
                self.__lfh.flush()
                
            if skip:
                return False

            mR=rsu.searchEntities(entityD=entityD,saveBlast=True,filePath=self.__sessionPath,filePrefix=dataSetId)
            fn=self.__pI.getReferenceSequenceFilePath(dataSetId,entityId=entityId,fileSource='session')
            rsio.writeMatchResults(entityD,outFilePath=fn,matchResults=mR)

            if (self.__verbose):
                endTime=time.clock()
                self.__lfh.write("+ReferenceSequenceDataUpdate.__doReferenceSearch()  completed at %s (%.2f seconds)\n" % 
                                 (time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime-startTime))
        except:
            self.__lfh.write("+ReferenceSequenceDataUpdate.__doReferenceSearch() - failing \n")
            traceback.print_exc(file=self.__lfh)

        return False            

    def __readReferenceSearchResults(self,dataSetId,entityId,fileSource='session',wfInstanceId=None,sessionCachePath=None):
        """  Read the reference sequence database search results for the input entity.

             Return (status,eRefL) completion status and a dictionary of reference matching details.

             Setting fileSource='session' and sessionCachePath=<cache path> allows input of saved search results.
        """
        eRefL=[]
        #
        try:
            pI=PathInfo(siteId=self.__siteId,sessionPath=self.__sessionPath,verbose=self.__verbose,log=self.__lfh)
            if ((fileSource == 'session') and (sessionCachePath is not None)):
                pI.setSessionPath(sessionCachePath)
            #
            fn=pI.getReferenceSequenceFilePath(dataSetId,entityId=entityId,wfInstanceId=wfInstanceId,fileSource=fileSource)                
            if (os.access(fn,os.R_OK)):
                rc0 = PdbxFileIo(verbose=self.__verbose,log=self.__lfh).getContainer(fn)
                rsin=ReferenceSequenceIo(dataContainer=rc0,verbose=self.__verbose,log=self.__lfh)
                eRefL=rsin.readMatchResults()
            else:
                if (self.__verbose):
                    self.__lfh.write("+ReferenceSequenceDataUpdate.__readReferenceSearchResults() - NO reference file for entity %s in path %s\n" % (entityId,fn))
                eRefL=[]
            return True,eRefL
        except:
            self.__lfh.write("+ReferenceSequenceDataUpdate.__readReferenceSearchResults() - failing \n")
            traceback.print_exc(file=self.__lfh)

        return False,eRefL        

    def __reloadReferenceSequences(self, entityD, eRefL):
        """ reload references sequences for the input entity.

            Inputs sources -- entityD, eRefL
         
        """
        #
        entityId=entityD['ENTITY_ID']
        startTime=time.clock()        
        if (self.__verbose):
            self.__lfh.write("+ReferenceSequenceDataUpdate.__reloadReferenceSequences() for entityId %r started at %s \n" % 
                             (entityId,time.strftime("%Y %m %d %H:%M:%S", time.localtime())))
        #
        seqFeature=SequenceFeature(verbose=self.__verbose)
        #
        pDList=entityD['PART_LIST']
        sId=entityD['IDENTIFIER']
        polyType=entityD['POLYMER_LINKING_TYPE']            
        polyTypeCode=self.__srd.getPolymerTypeCode(polyType)
        if (self.__verbose):
            self.__lfh.write("+ReferenceSequenceDataUpdate.__reloadReferenceSequences() sId %s entityId %r polyType %s polyTypeCode %s partlist %r\n" % (sId,entityId,polyType,polyTypeCode,pDList))
        #
        self.__sds.filterIndex(sId, dataType="sequence", seqType="ref")
        self.__sds.filterIndex(sId, dataType="feature",  seqType="ref")

        #
        for (partNo,pD) in enumerate(pDList,start=1):
            seqNumBeg=pD['SEQ_NUM_BEG']
            seqNumEnd=pD['SEQ_NUM_END']
            seqPartId=pD['SEQ_PART_ID']
            seqPartType=pD['SEQ_PART_TYPE']
            #
            altId=1
            if (self.__debug):
                self.__lfh.write("+ReferenceSequenceDataUpdate.__reloadReferenceSequences() entityId %r partNo %r partId %r type %s begin %r end %r\n" %
                                 (entityId, partNo, seqPartId, seqPartType,seqNumBeg,seqNumEnd))
            for rD  in eRefL:
                #
                if partNo == int(rD['fragment_id']):
                    myOrderId   = int(rD['id'])
                    iBegin      = int(rD['hitFrom'])
                    iEnd        = int(rD['hitTo'])                
                    seqS=rD['subject']
                    sTup3L= self.__srd.cnv1To3ListIdx(seqS,iBegin,polyTypeCode)
                    self.__sds.setSequence(sTup3L,  sId, 'ref',  partId=partNo, altId=altId, version=1)

                    if (self.__debug):
                        self.__lfh.write("+ReferenceSequenceDataUpdate.__reloadReferenceSequences() entityId %r partNo %r dbbegin %r  dbend %r seq %s\n" %
                                         (entityId, partNo, iBegin, iEnd, seqS))
                    seqFeature.clear()
                    # disambiguate the organism and strain data -
                    if rD['db_name'] in ['SP','TR','UNP']:
                        org,strain=seqFeature.decodeUniProtSourceOrganism(rD['source_scientific'])
                        seqFeature.setSource(organism=org, strain=strain,taxid=rD['taxonomy_id'],commonName=rD['source_common'])
                    else:
                        seqFeature.setSource(organism=rD['source_scientific'], taxid=rD['taxonomy_id'],commonName=rD['source_common'])
                    #
                    seqFeature.setId(dbName=rD['db_name'],dbCode=rD['db_code'], dbAccession=rD['db_accession'],dbIsoform=rD['db_isoform'])
                    seqFeature.setRefSeqNames(proteinName=rD['name'],synonyms=rD['synonyms'], geneName=rD['gene'])
                    seqFeature.setRefSeqDetails(enzymeClass=rD['ec'], description=rD['db_description'], comments=rD['comments'], keywords=rD['keyword'])
                    seqFeature.setItem('REF_MATCH_BEGIN',iBegin)
                    seqFeature.setItem('REF_MATCH_END',iEnd)
                    seqFeature.setItem('ORG_ORDER_ID', myOrderId)
                    seqFeature.setPolymerType(polyTypeCode)      
                    seqFeature.setRefSortOrder(sortIndex=rD['sort_order'], sortMetric=rD['sort_metric'])
                    #
                    seqFeature.setItem('AUTH_REF_SEQ_SIM_BLAST',rD['seq_sim'])
                    #
                    seqFeature.setAuthPartDetails(partNo,seqNumBeg,seqNumEnd,seqPartType)
                    self.__sds.setFeature( seqFeature.get(),  sId, 'ref',  partId=partNo, altId=altId, version=1)

                    altId+=1
                    #

        self.__sds.serialize()
        if (self.__verbose):
            endTime=time.clock()
            self.__lfh.write("+ReferenceSequenceDataUpdate.__reloadReferenceSequences()  for entity %r completed at %s (%.2f seconds)\n" % 
                             (entityId, time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                              endTime-startTime))


        return True
                
if __name__ == '__main__':
    pass
