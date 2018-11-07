##
# File:    SequenceDataUpdate.py
# Date:    18-Mar-2013
#
# Updates:
#  11-Apr-2013 jdw  adjust the computation of default aligned regions of input reference sequence
#  16-Sep-2013 jdw  add preliminary makeEntityReviewForm
#  11-Nov-2013 jdw  add part support to makeEntityReviewForm
#  14-Nov-2013 jdw  refactor methods grabbing selections
#  25-Nov-2013 jdw  more host org attributes in form
#  19-Jan-2014 jdw  "HOST_ORG_CELL_LINE"
#  12-Feb-2014 jdw  expand the annotations recovered from reference sequence database entries on reload.
#                   shoudl correspond to initial load now -
# 
# 3-entity example -
#   /service/sequence_editor/new_session/wf?classID=AnnMod&identifier=D_1000000000&filesource=wf-archive&instance=W_000
# chimera
#   /service/sequence_editor/new_session/wf?classID=AnnMod&identifier=D_1000000001&filesource=wf-archive&instance=W_000
##
"""
Utilities for adding out-of-band sequence and feature data to the sequence data store.

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.08"

import sys,traceback 

from wwpdb.apps.seqmodule.io.SequenceDataStore              import SequenceDataStore
from wwpdb.apps.seqmodule.io.ReferenceSequenceUtils         import ReferenceSequenceUtils
from wwpdb.apps.seqmodule.io.TaxonomyUtils                  import TaxonomyUtils
from wwpdb.apps.seqmodule.util.SequenceLabel                import SequenceLabel, SequenceFeature, SequenceFeatureMap
from wwpdb.apps.seqmodule.util.SequenceAssign               import SequenceAssignArchive,SequenceAssignDepositor,ReferenceSequence
from wwpdb.apps.seqmodule.util.SequenceReferenceData        import SequenceReferenceData
from wwpdb.apps.seqmodule.align.AlignmentStatistics         import AlignmentStatistics
from wwpdb.utils.rcsb.UtilDataStore                         import UtilDataStore

from wwpdb.utils.pair_align.wrapper.libPairwiseAlignPackage import PairwiseAlign

class SequenceDataUpdate(object):
    """ Utilities for adding out-of-band sequence and feature data to the sequence data store.

    """
    def __init__(self,reqObj=None,verbose=False,log=sys.stderr):
        self.__verbose=verbose
        self.__debug=False
        self.__reqObj=reqObj
        self.__lfh=log
        self.__defaultInsertSortMetric=100000
        #
        self.__srd=SequenceReferenceData(verbose=self.__verbose,log=self.__lfh)
        self.__gapSymbol=self.__srd.getGapSymbol()
        self.__maxRefAlign=100
        #
        self.__entityReviewItemList=["ENTITY_DESCRIPTION","ENTITY_SYNONYMS","SOURCE_ORGANISM","SOURCE_COMMON_NAME",
                                     "SOURCE_TAXID","SOURCE_GENE_NAME","SOURCE_STRAIN","SOURCE_VARIANT","ENTITY_ENZYME_CLASS","ENTITY_FRAGMENT_DETAILS",
                                     "HOST_ORG_SOURCE","HOST_ORG_COMMON_NAME","HOST_ORG_STRAIN",
                                     "HOST_ORG_TAXID","HOST_ORG_VECTOR","HOST_ORG_VECTOR_TYPE","HOST_ORG_PLASMID", "HOST_ORG_CELL_LINE","SOURCE_METHOD",
                                     "CURRENT_AUTH_SELECT_ID","CURRENT_REF_SELECT_ID"]

        self.__setup()

    def __setup(self):
        try:
            self.__placeHolderValue="click-to-edit"
            self.__siteId      = self.__reqObj.getValue("WWPDB_SITE_ID")
            self.__sessionId   = self.__reqObj.getSessionId()        
            self.__sessionObj  = self.__reqObj.getSessionObj()
            self.__sessionPath = self.__sessionObj.getPath()
            #
            self.__selectIdList=self.__reqObj.getSummarySelectList()
            self.__sds=SequenceDataStore(reqObj=self.__reqObj,verbose=self.__verbose,log=self.__lfh)
            #

            #
            #
        except:
            if (self.__verbose):
                self.__lfh.write("+SequenceDataUpdate.__setup() sessionId %s failed\n" % (self.__sessionObj.getId()))


    def __setNewRefId(self,refId):
        uds=UtilDataStore(reqObj=self.__reqObj,verbose=self.__verbose,log=self.__lfh)
        uds.set('new-ref-seq-id',refId)
        uds.serialize()

    def updateAuthEntityDetails(self,selectIdList=None):
        """ Update current auth/sample entity attribute details based on current selected reference sequences -
        """
        if selectIdList is not None:
            self.__selectIdList=selectIdList
        #
        self.__lfh.write("+SequenceDataUpdate.updateAuthEntityDetails() starting with selectIdList %r\n" % selectIdList)
        entityIdList=self.__sds.getGroupIds()
        for entityId in entityIdList:
            instanceIds=self.__sds.getGroup(groupId=entityId)
            if len(instanceIds)==0:
                continue
            instanceId=instanceIds[0]
            partIdList=self.__sds.getPartIds(instanceId, dataType="sequence", seqType="auth")
            #
            self.__lfh.write("+SequenceDataUpdate.updateAuthEntityDetails() entityId %r instanceId %r partIdList %r\n" % (entityId,instanceId,partIdList))
            sfm=SequenceFeatureMap(verbose=self.__verbose,log=self.__lfh)
            for partId in partIdList:
                authSelectId,authSL,authFdObj=self.__getCurrentAuthSelection(entityId,partId=partId)
                seqType,seqInstId,seqPartId,seqAltId,seqVersion=authSL.get()
                #
                #
                refSelectId,refSL,refFdObj=self.__getCurrentRefSelection(entityId,partId=partId)
                #
                # filter special cases -
                #
                curRefId=authFdObj.getCurrentRefSelectId()
                hasEdit=authFdObj.getManualEditStatus()
                self.__lfh.write("+SequenceDataUpdate.updateAuthEntityDetails() entityId %r partId %r refSelectId %r curRefId %r hasEdit %r\n" % 
                                 (entityId,partId,refSelectId,curRefId,hasEdit))
                if curRefId == refSelectId and hasEdit:
                    self.__lfh.write("+SequenceDataUpdate.updateAuthEntityDetails() skip update entityId %r partId %r curRefId %r hasEdit %r\n" % (entityId,partId,curRefId,hasEdit))
                    continue
                selfRefTarget='selfref_%s_%d' % (entityId,partId)
                if refSelectId.startswith(selfRefTarget):
                    # use original auth values -
                    sfm.updateAuthOrig(authFdObj.get())
                else:
                    sfm.updateAuth(authFdObj.get(),refFdObj.get())
                self.__sds.setFeature(authFdObj.get(),seqInstId,'auth',partId=partId,altId=seqAltId,version=seqVersion)

        self.__sds.serialize()
        return True
        
    def __getCurrentAuthSelection(self,entityId,partId=1):
        """ Search selections for author sequence part 1. in order to establish the 
            the selected version.   Then return the objects associated with the input
            partId.  
        """
        try:
            seqIdList=self.__sds.getGroup(entityId)
            instanceId=seqIdList[0]
            sL=SequenceLabel()
            for sId in self.__selectIdList:
                if sId.startswith('auth'):
                    sL.unpack(sId)
                    seqType,seqInstId,seqPartId,seqAltId,seqVersion=sL.get()
                    if instanceId == seqInstId and seqPartId == 1:
                        fdObj=self.__sds.getFeatureObj(seqInstId,seqType="auth",partId=partId,altId=seqAltId,version=seqVersion)
                        self.__lfh.write("+SequenceDataUpdate._getCurrentAuthSelection() returns: entity %r instance %r partId %r altId %r version %r\n" % 
                                         (entityId,seqInstId,seqPartId,seqAltId,seqVersion))
                        return sId,sL,fdObj
        except:
            self.__lfh.write("+SequenceDataUpdate._getCurrentAuthSelection() failed for selectList %r entityId %r partId %r\n" % (self.__selectIdList,entityId,partId))
            traceback.print_exc(file=self.__lfh)                        

        self.__lfh.write("+SequenceDataUpdate._getCurrentAuthSelection() no return for entityId %r partId %r\n" %  (entityId,partId))
        return None,None,None

    def __getCurrentRefSelection(self,entityId,partId=1):
        try:
            seqIdList=self.__sds.getGroup(entityId)
            instanceId=seqIdList[0]
            sL=SequenceLabel()
            selfRefTarget='selfref_%s_%d' % (entityId,partId)
            for sId in self.__selectIdList:
                if sId.startswith('ref'):
                    sL.unpack(sId)
                    seqType,seqInstId,seqPartId,seqAltId,seqVersion=sL.get()
                    if instanceId == seqInstId and partId == seqPartId:
                        fObj=self.__sds.getFeatureObj(seqInstId,seqType="ref",partId=seqPartId,altId=seqAltId,version=seqVersion)
                        self.__lfh.write("+SequenceDataUpdate._getCurrentRefSelection() returns: entity %r instance %r partId %r altId %r version %r\n" % 
                                         (entityId,seqInstId,seqPartId,seqAltId,seqVersion))
                        return sId,sL,fObj
                elif sId.startswith(selfRefTarget):
                    #JDW JDW
                    sL.set(seqType='ref', seqInstId=instanceId, seqPartId=partId, seqAltId=1, seqVersion=1)
                    fObj=SequenceFeature()
                    return sId,sL,fObj
        except:
            self.__lfh.write("+SequenceDataUpdate._getCurrentRefSelection() failed for selectList %r entityId %s partId %r\n" % (self.__selectIdList,entityId,partId))
            traceback.print_exc(file=self.__lfh)                        

        self.__lfh.write("+SequenceDataUpdate._getCurrentRefSelection() no return for entityId %r instanceId %r partId %r selectIdList %r\n" %  (entityId,instanceId,partId,self.__selectIdList))
        return None,None,None


    def __getCurrentRefDetails(self,entityId,partId=1):
        seqIdList=self.__sds.getGroup(entityId)
        seqId=seqIdList[0]
        self.__lfh.write("+SequenceDataUpdate.__getCurrentRefDetails() entityId %r partId %r seqId %r\n" % (entityId,partId,seqId))
        self.__lfh.write("+SequenceDataUpdate.__getCurrentRefDetails() selectIdList %r\n" % (self.__selectIdList))
        sL=SequenceLabel()
        for sId in self.__selectIdList:
            if sId.startswith('ref'):
                sL.unpack(sId)
                seqType,seqInstId,seqPartId,seqAltId,seqVersion=sL.get()
                self.__lfh.write("+SequenceDataUpdate.__getCurrentRefDetails() testing seqInstId %r seqPartId %r\n" % (seqInstId,seqPartId))               
                if seqId == seqInstId and partId == seqPartId:
                    fObj=self.__sds.getFeatureObj(seqInstId,seqType="ref",partId=seqPartId,altId=seqAltId,version=seqVersion)
                    fD=fObj.get()
                    return fD['DB_NAME'],fD['DB_CODE'],fD['DB_ACCESSION'],fD['REF_MATCH_BEGIN'],fD['REF_MATCH_END']

        return None,None,None,None,None

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
            self.__lfh.write("+SequenceDataUpdate.dump() entityId %s archive assignment count %d\n" % (entityId,nRef))
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
            self.__lfh.write("+SequenceDataUpdate.dump() entityId %s depositor assignment count %d\n" % (entityId,nRef))
            if nRef > 0:
                refL=sADep.getReferenceList(entityId=entityId)
                for ref in refL:
                    ref.printIt(self.__lfh)


    def fetchArchiveReferenceSequences(self):
        entityIdList=self.__sds.getGroupIds()

        seqAssignD=self.__sds.getReferenceAssignments()
        sA=SequenceAssignArchive(verbose=self.__verbose,log=self.__lfh)
        sA.set(seqAssignD)

        for entityId in entityIdList:
            nRef=sA.getReferenceCount(entityId=entityId)
            self.__lfh.write("+SequenceDataUpdate.fetchArchiveReferenceSequences() entityId %s archive assignment count %d\n" % (entityId,nRef))
            if nRef > 0:
                partD=self.__getPartDetails(entityId=entityId)
                refL=sA.getReferenceList(entityId=entityId)
                for ref in refL:
                    dbName,dbCode,dbAccession=ref.getDbReference()
                    seq_align_beg,seq_align_end=ref.getSeqAlignRange()
                    #
                    # Match the residue ranges to entity parts  -
                    #
                    pId=self.__assignPart(partD=partD,seqBegin=seq_align_beg,seqEnd=seq_align_end)
                    if pId > 0:
                        self.addReferenceSequence(entityId=entityId,partId=pId,dbName=dbName,dbAccession=dbAccession,refSeqBeg=seq_align_beg,refSeqEnd=seq_align_end)
                    else:
                        self.__lfh.write("+SequenceDataUpdate.fetchArchiveReferenceSequences() entityId %s  no part assigned for residue range %d - %d\n" 
                                         % (entityId,seq_align_beg,seq_align_end))                    
        if self.__debug:
            self.__sds.dump(self.__lfh)
        return True
        
    
    def addReferenceSequence(self,entityId, partId, dbName,dbAccession,refSeqBeg=None,refSeqEnd=None):
        fD=self.__getAuthFeatures(entityId=entityId,partId=partId)
        seqFeature=SequenceFeature()
        seqFeature.set(fD)
        #
        rso=self.__fetchReferenceSequence(dbName,dbAccession)
        packedSeqLabel=self.__loadReferenceSequence(entityId=entityId, partId=partId, refSeqObj=rso,refSeqBeg=refSeqBeg,refSeqEnd=refSeqEnd)
        if self.__debug:
            self.__sds.dump(self.__lfh)
        self.__sds.serialize()
        # update alignment stats - 
        alstat=AlignmentStatistics(reqObj=self.__reqObj,maxRefAlign=self.__maxRefAlign,verbose=self.__verbose,log=self.__lfh)
        alstat.doUpdate()

        return packedSeqLabel
            
        

    def __fetchReferenceSequence(self,dbName,dbAccession):
        """  Return a reference sequence object loading the content dictionary for the input reference sequence  - 

        """
        rsu=ReferenceSequenceUtils(siteId=self.__siteId,verbose=self.__verbose,log=self.__lfh)
        #
        if self.__verbose:
            self.__lfh.write("+SequenceDataUpdate.fetchReferenceSequence() fetch reference sequence  %s %s\n" % (dbName,dbAccession))            
        if dbName in ['UNP']:
            dt=rsu.fetchUniProt([dbAccession])
            if dt.has_key(dbAccession):
                d=dt[dbAccession]
            else:
                d={}
            if dbAccession[0] in ['P','Q','O']:
                d['db_name']='SP'                
            else:
                d['db_name']='TR'
        elif dbName in ['GB']:
            d=rsu.fetchNcbiGi(dbAccession)
            d['db_assession']=dbAccession
            d['db_name']=dbName
        else:
            d={}

        if self.__verbose:
            self.__lfh.write("\n+SequenceDataUpdate.fetchReferenceSequence() content dictionary for  %s %s\n" % (dbName,dbAccession))
            self.__lfh.write("\n+SequenceDataUpdate.fetchReferenceSequence() content dictionary keys  %s\n" % d.keys())
            for k,v in d.items():
                self.__lfh.write(" + Key %s value %r\n" % (k,v[:512]))    
        
        rso=ReferenceSequence()
        rso.set(d)
        return rso


    def __assignPart(self,partD,seqBegin,seqEnd):
        """  Assign the entity partId to the 
        """
        pMatch=-1
        for partId,pTup in partD.items():
            if ((seqBegin == pTup[1]) and (seqEnd == pTup[2])):
                pMatch=partId

        return pMatch

    def __setPartType(self,pList,pType):
        oL=[]
        pSelectList=['false' for p in pList]
        if pType is not None and len(pType) > 1:
            pListU=[p.upper() for p in pList]
            pTypeU=pType.upper()
            try:
                idx = pListU.index(pTypeU)
                pSelectList[idx]='true'
            except ValueError:
                idx = -1 
        tL=[]
        for pt,psel in zip(pList,pSelectList):
            tL.append('{"value":"%s","label":"%s","selected":%s}' % (pt,pt,psel))
        oL.append('[')
        oL.append( ','.join(tL) )
        oL.append(']')

        return ''.join(oL)
                      
                
    def makeSourceEditForm(self,entityId,entryId=''):
        """
            <div class="ief" data-ief-edittype="select" 
             data-ief-selectvalues="[{"value":"1","label":"Presentation Label","selected":true},{"value":"2","label":"Label 2","selected":false}]">

        """
        top_form_template='''
        <div id="sectaxonomy">
        <h3>Taxonomy Data Form for Entry %s Entity %s</h3>        
        <form name="formtaxonomy" id="formtaxonomy" action="/service/sequence_editor/respond_form/taxonomy" method="post" class="taxonomy_ajaxform">
            <input type="hidden" name="sessionid" value="%s" />
            <input type="hidden" name="entityid" value="%s" />
            <input type="hidden" name="numparts" value="%d" />
            <table>
            <tr> 
               <th>Part Id</th>
               <th>Taxonomy Id</th>
               <th>Seq Begin</th>
               <th>Seq End</th>
               <th>Part Type</th>
             </tr>
        '''
        #
        bottom_form_template='''
            </table>
        <br class="clearfloat" />
        <div class="width50 fltlft">Search Sequence Database: <input type="checkbox" name="seq_search_op" id="seq_search_op" /></div>
        <div class="width50 fltlft"></div><input type="submit" name="submit" value="Submit edits" class="disableonclick submitparentform"  /></div>
        <br class="clearfloat" />
            <!-- <input type="reset" name="reset" value="Reset" /> -->
        </form>
        </div>
        '''
        #
        optList=['N-terminal tag','C-terminal tag','Biological sequence','Linker']
        optSel=['false' for opt in optList]
        #
        #                 partD[pId]=(pId,pSeqBegin,pSeqEnd,pType,taxId)        
        partD=self.__getPartDetails(entityId)
        self.__lfh.write("+SequenceDataUpdate.buildSourceEditForm() part data\n")
        for k,v in partD.items():
            self.__lfh.write(" part %r  data:  %r\n" % (k,v))
 
        partIdList=partD.keys()
        partIdList.sort(lambda x, y: int(x) - int(y))        
        oL=[]
        oL.append(top_form_template % (entryId,entityId,self.__sessionId,entityId,len(partIdList)))
        for partId in partIdList:
            pIdT,seqBeg,seqEnd,pType,taxIdT=partD[partId]
            oL.append('<tr>')
            oL.append('<td><span id="p_%d_partid"   class="ief">%d</span></td>' % (partId,partId))
            if ((taxIdT is not None) and (len(taxIdT) > 0)):
                taxId = taxIdT  
                oL.append('<td><span id="p_%d_taxid"    class="ief">%s</span></td>' % (partId,taxId))
            else:
                taxId=self.__placeHolderValue
                oL.append('<td><span id="p_%d_taxid"    class="ief greyedout">%s</span></td>' % (partId,taxId))                
            oL.append('<td><span id="p_%d_seqbegin" class="ief">%s</span></td>' % (partId,seqBeg))
            oL.append('<td><span id="p_%d_seqend"   class="ief">%s</span></td>' % (partId,seqEnd))
            jTxt=self.__setPartType(optList,pType)
            oL.append('<td><span id="p_%d_seqtype"  class="ief" data-ief-edittype="select" data-ief-selectvalues=\'%s\'>%s</span>' % (partId,jTxt,pType))
            oL.append('</tr>')

        partId=int(partIdList[-1])+1
        taxId=self.__placeHolderValue
        seqBeg=self.__placeHolderValue
        seqEnd=self.__placeHolderValue
        pType=self.__placeHolderValue

        oL.append('<tr>')
        oL.append('<td><span id="p_%d_partid"   class="ief">%d</span></td>' % (partId,partId))
        oL.append('<td><span id="p_%d_taxid"    class="ief greyedout">%s</span></td>' % (partId,taxId))
        oL.append('<td><span id="p_%d_seqbegin" class="ief greyedout">%s</span></td>' % (partId,seqBeg))
        oL.append('<td><span id="p_%d_seqend"   class="ief greyedout">%s</span></td>' % (partId,seqEnd))
        jTxt=self.__setPartType(optList,pType)
        oL.append('<td><span id="p_%d_seqtype"  class="ief greyedout" data-ief-edittype="select" data-ief-selectvalues=\'%s\'>%s</span></td>' % (partId,jTxt,pType))
        oL.append('</tr>')

        partId+=1
        oL.append('<tr>')
        oL.append('<td><span id="p_%d_partid"   class="ief">%d</span></td>' % (partId,partId))
        oL.append('<td><span id="p_%d_taxid"    class="ief greyedout">%s</span></td>' % (partId,taxId))
        oL.append('<td><span id="p_%d_seqbegin" class="ief greyedout">%s</span></td>' % (partId,seqBeg))
        oL.append('<td><span id="p_%d_seqend"   class="ief greyedout">%s</span></td>' % (partId,seqEnd))
        jTxt=self.__setPartType(optList,pType)
        oL.append('<td><span id="p_%d_seqtype"  class="ief greyedout" data-ief-edittype="select" data-ief-selectvalues=\'%s\'>%s</span></td>' % (partId,jTxt,pType))
        oL.append('</tr>')

        partId+=1
        oL.append('<tr>')
        oL.append('<td><span id="p_%d_partid"   class="ief">%d</span></td>' % (partId,partId))
        oL.append('<td><span id="p_%d_taxid"    class="ief greyedout">%s</span></td>' % (partId,taxId))
        oL.append('<td><span id="p_%d_seqbegin" class="ief greyedout">%s</span></td>' % (partId,seqBeg))
        oL.append('<td><span id="p_%d_seqend"   class="ief greyedout">%s</span></td>' % (partId,seqEnd))
        jTxt=self.__setPartType(optList,pType)
        oL.append('<td><span id="p_%d_seqtype"  class="ief greyedout" data-ief-edittype="select" data-ief-selectvalues=\'%s\'>%s</span></td>' % (partId,jTxt,pType))
        oL.append('</tr>')

        partId+=1
        oL.append('<tr>')
        oL.append('<td><span id="p_%d_partid"   class="ief">%d</span></td>' % (partId,partId))
        oL.append('<td><span id="p_%d_taxid"    class="ief greyedout">%s</span></td>' % (partId,taxId))
        oL.append('<td><span id="p_%d_seqbegin" class="ief greyedout">%s</span></td>' % (partId,seqBeg))
        oL.append('<td><span id="p_%d_seqend"   class="ief greyedout">%s</span></td>' % (partId,seqEnd))
        jTxt=self.__setPartType(optList,pType)
        oL.append('<td><span id="p_%d_seqtype"  class="ief greyedout" data-ief-edittype="select" data-ief-selectvalues=\'%s\'>%s</span></td>' % (partId,jTxt,pType))
        oL.append('</tr>')

        partId+=1
        oL.append('<tr>')
        oL.append('<td><span id="p_%d_partid"   class="ief">%d</span></td>' % (partId,partId))
        oL.append('<td><span id="p_%d_taxid"    class="ief greyedout">%s</span></td>' % (partId,taxId))
        oL.append('<td><span id="p_%d_seqbegin" class="ief greyedout">%s</span></td>' % (partId,seqBeg))
        oL.append('<td><span id="p_%d_seqend"   class="ief greyedout">%s</span></td>' % (partId,seqEnd))
        jTxt=self.__setPartType(optList,pType)
        oL.append('<td><span id="p_%d_seqtype"  class="ief greyedout" data-ief-edittype="select" data-ief-selectvalues=\'%s\'>%s</span></td>' % (partId,jTxt,pType))
        oL.append('</tr>')

        partId+=1
        oL.append('<tr>')
        oL.append('<td><span id="p_%d_partid"   class="ief">%d</span></td>' % (partId,partId))
        oL.append('<td><span id="p_%d_taxid"    class="ief greyedout">%s</span></td>' % (partId,taxId))
        oL.append('<td><span id="p_%d_seqbegin" class="ief greyedout">%s</span></td>' % (partId,seqBeg))
        oL.append('<td><span id="p_%d_seqend"   class="ief greyedout">%s</span></td>' % (partId,seqEnd))
        jTxt=self.__setPartType(optList,pType)
        oL.append('<td><span id="p_%d_seqtype"  class="ief greyedout" data-ief-edittype="select" data-ief-selectvalues=\'%s\'>%s</span></td>' % (partId,jTxt,pType))
        oL.append('</tr>')

        partId+=1
        oL.append('<tr>')
        oL.append('<td><span id="p_%d_partid"   class="ief">%d</span></td>' % (partId,partId))
        oL.append('<td><span id="p_%d_taxid"    class="ief greyedout">%s</span></td>' % (partId,taxId))
        oL.append('<td><span id="p_%d_seqbegin" class="ief greyedout">%s</span></td>' % (partId,seqBeg))
        oL.append('<td><span id="p_%d_seqend"   class="ief greyedout">%s</span></td>' % (partId,seqEnd))
        jTxt=self.__setPartType(optList,pType)
        oL.append('<td><span id="p_%d_seqtype"  class="ief greyedout" data-ief-edittype="select" data-ief-selectvalues=\'%s\'>%s</span></td>' % (partId,jTxt,pType))
        oL.append('</tr>')

        oL.append(bottom_form_template)
        
        return '\n'.join(oL)
        #

    def sourceEditFormResponder(self):
        """  Update the polymer entity data store using user provided entity part, source and taxonomy content.
        
             Form data encoded in the input request object -- 
        """
        try:
            numParts=int(str(self.__reqObj.getValue("numparts")))
            entityId=self.__reqObj.getValue("entityid")
            #
            pD={}
            for partId in range(1,numParts+8):
                taxId=self.__reqObj.getValue("p_%d_taxid" % partId )
                seqBegin=self.__reqObj.getValue("p_%d_seqbegin" % partId )
                seqEnd=self.__reqObj.getValue("p_%d_seqend" % partId )
                seqPartType=self.__reqObj.getValue("p_%d_seqtype" % partId )
                if partId > numParts:
                    pD[partId]=(partId,str(seqBegin),str(seqEnd),str(seqPartType),str(taxId))
                else:
                    pD[partId]=(partId,int(seqBegin),int(seqEnd),str(seqPartType),str(taxId))

            #
            pOrgD=self.__getPartDetails(entityId=entityId,numExtra=3)
            isChanged=False
            for partId in range(1,numParts+8):
                if pOrgD[partId] != pD[partId] and partId <= numParts:
                    self.__lfh.write("+SequenceDataUpdate.sourceEditFormResponder() source differs at partId %d current %r next %r\n" 
                                     % (partId,pOrgD[partId],pD[partId]))
                    isChanged=True
                    break
                elif pD[partId][1] != self.__placeHolderValue and partId > numParts:
                    isChanged=True
                    break
            if isChanged:
                if (self.__verbose):
                    self.__lfh.write("+SequenceDataUpdate.sourceEditFormResponder() source data for entity %s updated with %r\n" % (entityId,pD))
                self.__updateAuthPartDetails(entityId=entityId,partD=pD)
            else:
                if (self.__verbose):
                    self.__lfh.write("+SequenceDataUpdate.sourceEditFormResponder() form values unchanged\n")                
        except:
            self.__lfh.write("+SequenceDataUpdate.sourceEditResponder() failing\n")
            traceback.print_exc(file=self.__lfh)            
            return False

        if (self.__debug):
            self.__sds.dump(self.__lfh)

        return True


    def seqDbRefFormResponder(self):
        """  Update the sequence data store using data sequence database reference
        
             Form data encoded in the input request object -- 

             Return the packed sequence label of the added reference sequence or None
        """
        packedSeqLabel=None
        try:
            partId=int(str(self.__reqObj.getValue("partid")))
            entityId=self.__reqObj.getValue("entityid")
            dbName=self.__reqObj.getValue("dbname")
            dbAccession=self.__reqObj.getValue("dbaccession")
            dbSeqBegin=self.__reqObj.getValue("dbseqbegin")
            dbSeqEnd=self.__reqObj.getValue("dbseqend")
            if ((dbName == self.__placeHolderValue) or (dbAccession == self.__placeHolderValue)):
                return packedSeqLabel
            if dbSeqBegin == self.__placeHolderValue:
                dbSeqBegin = None
            else:
                try:
                    dbSeqBegin = int(str(dbSeqBegin))
                except:
                    dbSeqBegin = None                    

            if dbSeqEnd == self.__placeHolderValue:
                dbSeqEnd = None
            else:
                try:
                    dbSeqEnd = int(str(dbSeqEnd))
                except:
                    dbSeqEnd = None
                    

            packedSeqLabel=self.addReferenceSequence(entityId, partId, dbName,dbAccession,refSeqBeg=dbSeqBegin,refSeqEnd=dbSeqEnd)
            self.__reqObj.setNewRefId(packedSeqLabel)

        except:
            self.__lfh.write("+SequenceDataUpdate.seqDbRefResponder() failing\n")
            traceback.print_exc(file=self.__lfh)            

        return packedSeqLabel


    def makeSeqdbrefEditForm(self,entityId,partId=1,entryId=''):
        """    Return a preliminary form to input sequence database references. 

        """
        #
        form_template='''
        <div id="sectseqdbref">
        <h3>Reference Sequence Database Data Form for Entry %(entryid)s Entity %(entityid)s Part %(partid)s</h3>
        <form name="formseqdbref_%(partid)s" id="formseqdbref_%(partid)s" action="/service/sequence_editor/respond_form/seqdbref" method="post" class="auth_ajaxform">
            <input type="hidden" name="sessionid" value="%(sessionid)s" />
            <input type="hidden" name="partid" value="%(partid)s" />
            <input type="hidden" name="entityid" value="%(entityid)s" />            
            <table>
            <tr> 
               <th>Database name</th>
               <th>Accession Code</th>
               <th>Seq Begin</th>
               <th>Seq End</th>
             </tr>
            <tr>
            <td><span id="dbname"      class="ief %(dbname_css)s" data-ief-edittype="select" data-ief-selectvalues='[{"value":"UNP","label":"UNP","selected":false},{"value":"GB","label":"GB","selected":false}]'>%(dbname)s</span></td>
            <td><span id="dbaccession" class="ief %(dbaccession_css)s">%(dbaccession)s</span></td>
            <td><span id="dbseqbegin"  class="ief %(dbseqbegin_css)s">%(dbseqbegin)s</span></td>
            <td><span id="dbseqend"    class="ief %(dbseqend_css)s">%(dbseqend)s</span></td>
            </tr>
            </table>
            <input type="submit" name="submit" value="Submit" class="disableonclick submitparentform" />
            <!-- <input type="reset" name="reset" value="Reset" /> --> 
        </form>
       </div>
        '''
        #
        dbName,dbCode,dbAccession,dbSeqBegin,dbSeqEnd=self.__getCurrentRefDetails(entityId,partId=int(partId))
        #

        pD={}
        pD['entryid']=entryId
        pD['sessionid']=self.__sessionId
        pD['partid']=partId
        pD['entityid']=entityId        
        #
        pD['dbname']=dbName           
        pD['dbaccession']=dbAccession 
        pD['dbseqbegin']=dbSeqBegin   
        pD['dbseqend']=dbSeqEnd       
        #
        itemList=['dbname','dbaccession','dbseqbegin','dbseqend']
        for item in itemList:
            kyCss=item+"_css"
            if pD[item] is None or len(str(pD[item])) < 1:
                pD[item]=self.__placeHolderValue
                pD[kyCss]="greyedout"
            else:
                pD[kyCss]=""


        rD = {}        
        rD['htmlcontent']=(form_template % pD)
        #
        return rD


    def __getAuthRefAssignments(self,entityId):
        """  Get author provided reference assignments.   Return a lists of reference database accessions.

             Note that this is currently being do a the level of 
             entity as the entity parts are not explicit in input data at this time. 2013-11-14
        """
        depSeqAssignD=self.__sds.getDepositorReferenceAssignments()
        sADep=SequenceAssignDepositor(verbose=self.__verbose,log=self.__lfh)
        sADep.set(depSeqAssignD)
        refSeqAssignL=sADep.getReferenceList(entityId)

        if (self.__debug):
            self.__lfh.write("+SequenceDataUpdate.__getAuthRefAssignments() for entityId %r reference count is %d\n" % (entityId, sADep.getReferenceCount(entityId)))
            sADep.printIt(log=self.__lfh)
            for ii,rsa in enumerate(refSeqAssignL):
                self.__lfh.write("+SequenceDataUpdate.__getAuthRefAssignments() depositor reference  %d\n" % (ii+1))
                rsa.printIt(self.__lfh)
        #
        authDbNameList=[]
        authDbAccessionList=[]
        authDbCodeList=[]
        for rsa in refSeqAssignL:
            dbName,dbCode,dbAccession=rsa.getDbReference()
            if dbName not in ['.','?']:
                authDbNameList.append(dbName)            
            else:
                authDbNameList.append('')
            if dbAccession not in ['.','?']:                
                authDbAccessionList.append(dbAccession)
            else:
                authDbAccessionList.append('')                
            if dbCode not in ['.','?']:                
                authDbCodeList.append(dbCode)
            else:
                authDbCodeList.append('')                
        #
        if (self.__verbose):
            self.__lfh.write("+SequenceDataUpdate.__getAuthRefAssignments()  dbNames %r dbAccessions %r dbCodes %r\n" % (authDbNameList,authDbAccessionList,authDbCodeList))

        return authDbNameList,authDbAccessionList,authDbCodeList



    def makeEntityReviewForm(self,entityId,entryId=''):
        """  Wrapper for form generation entity detail input. Scope of form is for each entity including all parts -- 

        """
        rD = {}        
        rD['htmlcontent']=''
        seqIds=self.__sds.getGroup(groupId=entityId)
        if len(seqIds)==0:
            return rD
        seqId0=seqIds[0]
        curPartIdList=self.__sds.getPartIds(seqId0, dataType="sequence", seqType="auth")
        
        rL=[]
        rL.append('<div id="sectentityreview">')
        for partId in curPartIdList:
            rL.append(self.__makeEntityReviewFormPart(entityId,partId,entryId=entryId,numParts=len(curPartIdList)))
        
        rL.append('</div>')
        rD['htmlcontent']=''.join(rL)
        
        return rD


    def __makeEntityReviewFormPart(self,entityId,partId=1,entryId='',numParts=1):
        """    Return an edit form for to review entity annotations.
        """
        #
        form_template='''

        <h3>Entity Review Form for Entry %(entryid)s Entity %(entityid)s Part %(partid)s</h3>

        <form name="formentityreview_%(entityid)s_%(partid)s" id="formentityreview_%(entityid)s_%(partid)s" action="/service/sequence_editor/respond_form/entityreview" method="post" class="review_ajaxform">
            <input type="submit" name="submit" value="%(savebuttontext)s" class="disableonclick submitparentform fltrgt" />
            <input type="hidden" name="sessionid" value="%(sessionid)s" />
            <input type="hidden" name="partid" value="%(partid)s" />
            <input type="hidden" name="entityid" value="%(entityid)s" />  
            <input type="hidden" name="instid" value="%(instid)s" />  
            <input type="hidden" name="altid" value="%(altid)s" />  
            <input type="hidden" name="versionid" value="%(versionid)s" />  

            <input type="hidden" name="CURRENT_AUTH_SELECT_ID" value="%(CURRENT_AUTH_SELECT_ID)s" />          
            <input type="hidden" name="CURRENT_REF_SELECT_ID" value="%(CURRENT_REF_SELECT_ID)s" />          
            <br class="clearfloat" />
            <table>
            <tr> 
               <th>Item </th>
               <th>Current Value</th>
               <th>Reference DB</th>
               <th>Author Provided</th>
            </tr>

             <tr><td>DB Name</td>
                 <td></td>
                 <td class="bgcolcurrent"><span id="REF_DB_NAME" class="my-cell-static">%(REF_DB_NAME)s</span></td>                
                 <td class="bgcolauthor"><span id="AUTH_DB_NAME_LIST" class="my-cell-static">%(AUTH_DB_NAME_LIST)s</span></td>
             </tr>

             <tr>
                <td>DB Code </td>
                 <td></td>
                <td class="bgcolcurrent"><span id="REF_DB_CODE" class="my-cell-static">%(REF_DB_CODE)s</span></td>                
                <td class="bgcolauthor"><span id="AUTH_DB_CODE_LIST" class="my-cell-static">%(AUTH_DB_CODE_LIST)s</span></td>
             </tr>

             <tr>
                <td>DB Accession </td>
                <td></td>
                <td class="bgcolcurrent"><span id="REF_DB_ACCESSION" class="my-cell-static">%(REF_DB_ACCESSION)s</span></td> 
                <td class="bgcolauthor"><span id="AUTH_DB_ACCESSION_LIST" class="my-cell-static">%(AUTH_DB_ACCESSION_LIST)s</span></td>
             </tr>

             <!-- START  HERE -->
             <tr>
                <td>Name/Description </td>
                <td class="bgcolcurrent"><span id="ENTITY_DESCRIPTION" class="ief my-editable-cell %(ENTITY_DESCRIPTION_CSS)s">%(ENTITY_DESCRIPTION)s</span></td>
                <td class="bgcolcurrent"><span id="REF_DB_MOLECULE_NAME" class="my-cell-static">%(REF_DB_MOLECULE_NAME)s</span></td>                
                <td class="bgcolauthor"><span id="ENTITY_DESCRIPTION_ORIG" class="my-cell-static">%(ENTITY_DESCRIPTION_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Synonyms </td>
                <td class="bgcolcurrent"><span id="ENTITY_SYNONYMS" class="ief my-editable-cell %(ENTITY_SYNONYMS_CSS)s">%(ENTITY_SYNONYMS)s</span></td>
                <td class="bgcolcurrent"><span id="REF_DB_MOLECULE_SYNONYMS" class="my-cell-static">%(REF_DB_MOLECULE_SYNONYMS)s</span></td>                
                <td class="bgcolauthor"><span id="ENTITY_SYNONYMS_ORIG" class="my-cell-static">%(ENTITY_SYNONYMS_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Source Organism </td>
                <td class="bgcolcurrent"><span id="SOURCE_ORGANISM" class="ief my-editable-cell %(SOURCE_ORGANISM_CSS)s">%(SOURCE_ORGANISM)s</span></td>
                <td class="bgcolcurrent"><span id="REF_DB_SOURCE_ORGANISM" class="my-cell-static">%(REF_DB_SOURCE_ORGANISM)s</span></td>                
                <td class="bgcolauthor"><span id="SOURCE_ORGANISM_ORIG" class="my-cell-static">%(SOURCE_ORGANISM_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Source Common Name </td>
                <td class="bgcolcurrent"><span id="SOURCE_COMMON_NAME" class="ief my-editable-cell %(SOURCE_COMMON_NAME_CSS)s">%(SOURCE_COMMON_NAME)s</span></td>
                <td class="bgcolcurrent"><span id="REF_DB_SOURCE_COMMON_NAME" class="my-cell-static">%(REF_DB_SOURCE_COMMON_NAME)s</span></td>                
                <td class="bgcolauthor"><span id="SOURCE_COMMON_NAME_ORIG" class="my-cell-static">%(SOURCE_COMMON_NAME_ORIG)s</span></td>
             </tr>


             <tr>
                <td>Taxonomy ID</td>
                <td class="bgcolcurrent"><span id="SOURCE_TAXID" class="ief my-editable-cell %(SOURCE_TAXID_CSS)s">%(SOURCE_TAXID)s</span></td>
                <td class="bgcolcurrent"><span id="REF_DB_SOURCE_TAXID" class="my-cell-static">%(REF_DB_SOURCE_TAXID)s</span></td>                
                <td class="bgcolauthor"><span id="SOURCE_TAXID_ORIG" class="my-cell-static">%(SOURCE_TAXID_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Gene Name</td>
                <td class="bgcolcurrent"><span id="SOURCE_GENE_NAME" class="ief my-editable-cell %(SOURCE_GENE_NAME_CSS)s">%(SOURCE_GENE_NAME)s</span></td>
                <td class="bgcolcurrent"><span id="REF_DB_GENE_NAME" class="my-cell-static">%(REF_DB_GENE_NAME)s</span></td>                
                <td class="bgcolauthor"><span id="SOURCE_GENE_NAME_ORIG" class="my-cell-static">%(SOURCE_GENE_NAME_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Strain Name</td>
                <td class="bgcolcurrent"><span id="SOURCE_STRAIN"        class="ief my-editable-cell %(SOURCE_STRAIN_CSS)s">%(SOURCE_STRAIN)s</span></td>
                <td class="bgcolcurrent"><span id="REF_DB_SOURCE_STRAIN" class="my-cell-static">%(REF_DB_SOURCE_STRAIN)s</span></td>                
                <td class="bgcolauthor"><span id="SOURCE_STRAIN_ORIG"   class="my-cell-static">%(SOURCE_STRAIN_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Source Variant</td>
                <td class="bgcolcurrent"><span id="SOURCE_VARIANT"        class="ief my-editable-cell %(SOURCE_VARIANT_CSS)s">%(SOURCE_VARIANT)s</span></td>
                <td></td>             
                <td class="bgcolauthor"><span id="SOURCE_VARIANT_ORIG"   class="my-cell-static">%(SOURCE_VARIANT_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Enzyme Classification</td>
                <td class="bgcolcurrent"><span id="ENTITY_ENZYME_CLASS"        class="ief my-editable-cell %(ENTITY_ENZYME_CLASS_CSS)s">%(ENTITY_ENZYME_CLASS)s</span></td>
                <td class="bgcolcurrent"><span id="REF_DB_ENZYME_CLASS" class="my-cell-static">%(REF_DB_ENZYME_CLASS)s</span></td>                
                <td class="bgcolauthor"><span id="ENTITY_ENZYME_CLASS_ORIG"   class="my-cell-static">%(ENTITY_ENZYME_CLASS_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Fragment Details</td>  
                <td class="bgcolcurrent"><span id="ENTITY_FRAGMENT_DETAILS" class="ief my-editable-cell %(ENTITY_FRAGMENT_DETAILS_CSS)s">%(ENTITY_FRAGMENT_DETAILS)s</span></td>
                <td></td>             
                <td class="bgcolauthor"><span id="ENTITY_FRAGMENT_DETAILS_ORIG" class="my-cell-static">%(ENTITY_FRAGMENT_DETAILS_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Host Organism</td>
                <td class="bgcolcurrent"><span id="HOST_ORG_SOURCE" class="ief my-editable-cell %(HOST_ORG_SOURCE_CSS)s">%(HOST_ORG_SOURCE)s</span></td>
                <td></td>             
                <td class="bgcolauthor"><span id="HOST_ORG_SOURCE_ORIG" class="my-cell-static">%(HOST_ORG_SOURCE_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Host Common Name</td>
                <td class="bgcolcurrent"><span id="HOST_ORG_COMMON_NAME" class="ief my-editable-cell %(HOST_ORG_COMMON_NAME_CSS)s">%(HOST_ORG_COMMON_NAME)s</span></td>
                <td></td>             
                <td class="bgcolauthor"><span id="HOST_ORG_COMMON_NAME_ORIG" class="my-cell-static">%(HOST_ORG_COMMON_NAME_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Host Organism Strain</td>
                <td class="bgcolcurrent"><span id="HOST_ORG_STRAIN" class="ief my-editable-cell %(HOST_ORG_STRAIN_CSS)s">%(HOST_ORG_STRAIN)s</span></td>
                <td></td>             
                <td class="bgcolauthor"><span id="HOST_ORG_STRAIN_ORIG" class="my-cell-static">%(HOST_ORG_STRAIN_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Host Organism TaxID</td>
                <td class="bgcolcurrent"><span id="HOST_ORG_TAXID" class="ief my-editable-cell %(HOST_ORG_TAXID_CSS)s">%(HOST_ORG_TAXID)s</span></td>
                <td></td>             
                <td class="bgcolauthor"><span id="HOST_ORG_TAXID_ORIG" class="my-cell-static">%(HOST_ORG_TAXID_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Host Organism Vector</td>
                <td class="bgcolcurrent"><span id="HOST_ORG_VECTOR" class="ief my-editable-cell %(HOST_ORG_VECTOR_CSS)s">%(HOST_ORG_VECTOR)s</span></td>
                <td></td>             
                <td class="bgcolauthor"><span id="HOST_ORG_VECTOR_ORIG" class="my-cell-static">%(HOST_ORG_VECTOR_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Host Organism Vector Type</td>
                <td class="bgcolcurrent"><span id="HOST_ORG_VECTOR_TYPE" class="ief my-editable-cell %(HOST_ORG_VECTOR_TYPE_CSS)s">%(HOST_ORG_VECTOR_TYPE)s</span></td>
                <td></td>             
                <td class="bgcolauthor"><span id="HOST_ORG_VECTOR_TYPE_ORIG" class="my-cell-static">%(HOST_ORG_VECTOR_TYPE_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Host Organism Plasmid</td>
                <td class="bgcolcurrent"><span id="HOST_ORG_PLASMID" class="ief my-editable-cell %(HOST_ORG_PLASMID_CSS)s">%(HOST_ORG_PLASMID)s</span></td>
                <td></td>             
                <td class="bgcolauthor"><span id="HOST_ORG_PLASMID_ORIG" class="my-cell-static">%(HOST_ORG_PLASMID_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Host Organism Cell Line</td>
                <td class="bgcolcurrent"><span id="HOST_ORG_CELL_LINE" class="ief my-editable-cell %(HOST_ORG_CELL_LINE_CSS)s">%(HOST_ORG_CELL_LINE)s</span></td>
                <td></td>             
                <td class="bgcolauthor"><span id="HOST_ORG_CELL_LINE_ORIG" class="my-cell-static">%(HOST_ORG_CELL_LINE_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Source Method</td>
                <td class="bgcolcurrent"><span id="SOURCE_METHOD" class="ief my-editable-cell" data-ief-edittype="select"  data-ief-selectvalues='%(SOURCE_METHOD_SELECT)s'>%(SOURCE_METHOD)s</span></td>
                <td></td>             
                <td class="bgcolauthor"><span id="SOURCE_METHOD_ORIG" class="my-cell-static">%(SOURCE_METHOD_ORIG)s</span></td>
             <tr>
            </table>
            <input type="submit" name="submit" value="%(savebuttontext)s" class="disableonclick submitparentform fltrgt" />
            <!-- <input type="reset" name="reset" value="Reset" /> --> 
        </form>
        '''
        #
        rL=[]       
        #
        authSelectId,authSL,authFdObj=self.__getCurrentAuthSelection(entityId,partId=partId)
        refSelectId,refSL,refFdObj  =self.__getCurrentRefSelection(entityId,partId=partId)

        if authSL is None or authFdObj is None:
            return rL

        authSeqType,authSeqInstId,authSeqPartId,authSeqAltId,authSeqVersion=authSL.get()

        
        authFD=authFdObj.get()
        refFD=refFdObj.get()

        # This information is consolidated by entity here -
        authDbNameList,authDbAccessionList,authDbCodeList= self.__getAuthRefAssignments(entityId)

        pD={}
        pD['AUTH_DB_NAME_LIST']=','.join(authDbNameList)
        pD['AUTH_DB_CODE_LIST']=','.join(authDbCodeList)
        pD['AUTH_DB_ACCESSION_LIST']=','.join(authDbAccessionList)

        pD['REF_DB_NAME']=refFD['DB_NAME']
        pD['REF_DB_CODE']=refFD['DB_CODE']
        pD['REF_DB_ACCESSION']=refFD['DB_ACCESSION']
        pD['REF_DB_MOLECULE_NAME']=refFD['DB_MOLECULE_NAME']
        pD['REF_DB_MOLECULE_SYNONYMS']=refFD['DB_MOLECULE_SYNONYMS']
        pD['REF_DB_SOURCE_ORGANISM']=refFD['SOURCE_ORGANISM']
        pD['REF_DB_SOURCE_COMMON_NAME']=refFD['SOURCE_COMMON_NAME']
        pD['REF_DB_SOURCE_TAXID']=refFD['SOURCE_TAXID']
        pD['REF_DB_SOURCE_STRAIN']=refFD['SOURCE_STRAIN']
        pD['REF_DB_GENE_NAME']=refFD['DB_GENE_NAME']
        pD['REF_DB_ENZYME_CLASS']=refFD['DB_MOLECULE_EC']
        
        #                     
        if authFD['SOURCE_METHOD'] == "MAN":
            tup=('true','false','false')
        elif authFD['SOURCE_METHOD'] == "NAT":
            tup=('false','true','false')
        elif  authFD['SOURCE_METHOD'] == "SYN":
            tup=('false','false','true')
        else:
            tup=('false','false','false')

        pD['SOURCE_METHOD_SELECT']='[{"value":"MAN","label":"MAN","selected":%s},{"value":"NAT","label":"NAT","selected":%s},{"value":"SYN","label":"SYN","selected":%s}]' % tup

        #
        pD['entryid']=entryId
        pD['sessionid']=self.__sessionId
        pD['partid']=partId
        pD['entityid']=entityId  
        pD['instid']=authSeqInstId
        pD['altid']=authSeqAltId
        pD['versionid']=authSeqVersion  
        if numParts > 1:
            pD['savebuttontext']='Save Edits for Entity '+ str(entityId) + ' Part ' + str(partId) + '  (save each part separately)'
        else:
            pD['savebuttontext']='Save Edits for Entity '+ str(entityId) 
        #

        #
        for item in self.__entityReviewItemList:
            kyCss=item+"_CSS"
            if authFD[item] is None or len(authFD[item]) < 1:
                authFD[item]=self.__placeHolderValue
                pD[kyCss]="greyedout"
            else:
                pD[kyCss]=""

        pD.update(authFD)                
        pD['CURRENT_AUTH_SELECT_ID']=authSelectId
        pD['CURRENT_REF_SELECT_ID']=refSelectId


        #
        if (self.__debug):
            self.__lfh.write("+SequenceDataUpdate.makeEntityReview() reference sequence feature dictionary contents:\n")
            for k in sorted(refFD.keys()):
                v=refFD[k]
                self.__lfh.write(" ++++    %-40s   --  %s\n" % (k,v))                
                self.__lfh.write("+SequenceDataUpdate.makeEntityReview() author sequence feature dictionary contents:\n")
                for k in sorted(pD.keys()):
                    v=pD[k]
                    self.__lfh.write(" ++++    %-40s   --  %s\n" % (k,v))                
        rL=(form_template % pD)
        #
        return rL

    def entityReviewFormResponder(self):
        """  Update the entity annotation from form data - 
        
             Form data encoded in the input request object -- 

             Return True
        """
        try:

            partId=int(str(self.__reqObj.getValue("partid")))
            entityId=self.__reqObj.getValue("entityid")
            instId=self.__reqObj.getValue("instid")
            altId=self.__reqObj.getValue("altid")
            versionId=self.__reqObj.getValue("versionid")
            #
            authFdObj=self.__sds.getFeatureObj(instId,seqType='auth',partId=partId,altId=altId,version=versionId)

            fDUpd={}
            for ky in self.__entityReviewItemList:
                fDUpd[ky]=self.__reqObj.getValue(ky)
                if fDUpd[ky]==self.__placeHolderValue:
                    fDUpd[ky]=''
                if self.__verbose:
                    self.__lfh.write("+SequenceDataUpdate.entityReviewResponder() ky=%40s   -- value=%r\n" % (ky,fDUpd[ky]) )
            #
            # mark this sequence feature object as edited - to screen for overwrites -
            #
            fDUpd['HAS_MANUAL_EDIT']=True
            #
            authFdObj.set(fDUpd,resetAll=False)
            #
            self.__sds.setFeatureObj(authFdObj,instId,seqType='auth',partId=partId,altId=altId,version=versionId)

            if (self.__verbose):
                self.__lfh.write("+SequenceDataUpdate.entityReviewFormResponder() updating AUTH sequence %s alt %r part %r version %r\n" % (instId,altId,partId,versionId))
                authFdObj.printIt(log=self.__lfh)
            #
            self.__sds.serialize()
        except:
            self.__lfh.write("+SequenceDataUpdate.entityReviewResponder() failing\n")
            traceback.print_exc(file=self.__lfh)            

        return True


    def __getAuthSeqComment(self,idx, partD, seqBegMin, seqEndMax):
        """      Set default annotation between part -- and clean up any other comments -- 
                 Return initial annotation assignment for residue position 'idx' based on the 
                 the input part definition OR an out-of-part assignment of a terminal expression 
                 tag or an internal linker. 
                 
                 input    idx   sequence position in author/sample sequence (1-N)
                          partD[partId]=(id,seqBeg,seqEnd,partTypeComment,...)
 
        """
        # Are we in defined entity part ?
        for pId,rl in partD.items():
            if idx < rl[1] or idx > rl[2]:
                continue
            else:
                return rl[3]

        # not in part -- 
        if idx < seqBegMin or idx > seqEndMax:
            return 'expression tag'
        else:
            return 'linker'

            
    def __updateAuthPartDetails(self,entityId,partD):
        #
        # create new version of all parts -- 
        #
        if (self.__verbose):
            self.__lfh.write("\n\n+SequenceDataUpdate.__updateAuthPartDetails() Starting for entity %s and part dictionary %r\n" % (entityId,partD))
        #
        seqBegMin=100
        seqEndMax=0
        pD={}
        for tId, pTup in partD.items():
            pId,seqBeg,seqEnd,seqPartType,taxId=pTup
            if (self.__verbose):
                self.__lfh.write("+SequenceDataUpdate.__updateAuthPartDetails() Entity %s part %d new assignment form values %s\n" % (entityId,pId,pTup))
            #
            if ((seqPartType == self.__placeHolderValue) or  (seqBeg == self.__placeHolderValue) or (seqEnd == self.__placeHolderValue)):
                continue
            seqNumBeg=int(seqBeg)
            seqNumEnd=int(seqEnd)
            seqBegMin=min(seqBegMin, seqNumBeg)
            seqEndMax=max(seqEndMax, seqNumEnd)
            spt=''
            if  seqPartType in ['biological sequence']:
                spt=''
            elif seqPartType in ['n-terminal tag','c-terminal tag']:
                spt='expression tag'
            elif seqPartType in ['linker']:
                spt='linker'
            pD[pId]=(pId,seqNumBeg,seqNumEnd,spt)
                
            
        #
        if (self.__verbose):
            self.__lfh.write("\n\n+SequenceDataUpdate.__updateAuthPartDetails() entity %s seqBegMin %d   seqEndMax %d\n" % (entityId,seqBegMin,seqEndMax))

        tU=TaxonomyUtils(siteId=self.__siteId,verbose=self.__verbose,log=self.__lfh)

        seqIds=self.__sds.getGroup(groupId=entityId)
        if len(seqIds)==0:
            return False

        seqFeature=SequenceFeature()
        seqId0=seqIds[0]
        curPartIdList=self.__sds.getPartIds(seqId0, dataType="sequence", seqType="auth")

        # base versioning on the base partId = 1
        # 
        vL=self.__sds.getVersionIds(seqId0,partId=1, altId=1, dataType="sequence", seqType="auth")
        if len(vL)<1:
            return False
        # target version for update.    
        #lastVer=int(str(vL[0]))
        curVer=int(str(vL[0]))
        nextVer=curVer+1
        # JDW      authSelectId,authSL,authFdObj=self.__getCurrentAuthSelection(entityId,partId=partId)
        #          curSeqType,curInstId,curPartId,curAltId,curVer=authSL.get()
        #
        for partId, pTup in partD.items():
            authSelectId,authSL,authFdObj=self.__getCurrentAuthSelection(entityId,partId=partId)
            curSeqType,curInstId,curPartId,curAltId,curVer=authSL.get()
            if partId in curPartIdList:
                # Existing part ---
                # fObj=self.__sds.getFeatureObj(seqId0,seqType="auth",partId=partId,altId=1,version=curVer)
                taxIdOrg=authFdObj.getSourceTaxId()
                # values - 
                pId,seqBeg,seqEnd,seqPartType,taxId=pTup
                if (self.__verbose):
                    self.__lfh.write("+SequenceDataUpdate.__updateAuthPartDetails() Taxonomy assignments for part %d org taxid %s new taxid %s\n" % (partId,taxIdOrg,taxId))
                if taxIdOrg != taxId:
                    nL=tU.lookUpSource(taxId=taxId)                   
                    if len(nL) > 0:
                        authFdObj.setSource(organism=nL[0],taxid=taxId)
                    else:
                        authFdObj.setTaxId(taxid=taxId)
                authFdObj.setAuthPartDetails(partId,seqBeg,seqEnd,seqPartType)
                self.__sds.setFeatureObj(authFdObj,seqId0,seqType='auth',partId=partId,altId=1,version=nextVer)
                # Update the sequence...
                seqIdx=self.__sds.getSequence(seqId0,seqType="auth",partId=partId,altId=1,version=curVer)
                ####
                ###  JDW reset default comment state ---
                sTupL=[]
                for sT in seqIdx:
                    comment=self.__getAuthSeqComment(sT[3], pD,seqBegMin,seqEndMax)
                    sTupL.append( (sT[0],sT[1],comment,sT[3]))
                ####
                #self.__sds.setSequence(seqIdx,seqId0,seqType='auth',partId=partId,altId=1,version=nextVer)
                self.__sds.setSequence(sTupL,seqId0,seqType='auth',partId=partId,altId=1,version=nextVer)
                if (self.__verbose):
                    self.__lfh.write("+SequenceDataUpdate.__updateAuthPartDetails() updating part %d from version %d to %d\n" % (partId,curVer,nextVer))  
            else:
                # New part --- 
                pId,seqBeg,seqEnd,seqPartType,taxId=pTup
                if (self.__verbose):
                    self.__lfh.write("+SequenceDataUpdate.__updateAuthPartDetails() Entity %s part %d new assignment form values %s\n" % (entityId,partId,pTup))
                #
                if ((seqPartType == self.__placeHolderValue) or  (seqBeg == self.__placeHolderValue) or (seqEnd == self.__placeHolderValue)):
                    continue

                # 
                if ((taxId is None) or (taxId == self.__placeHolderValue)):
                    taxId=''
                    if (self.__verbose):
                        self.__lfh.write("+SequenceDataUpdate.__updateAuthPartDetails() Handle missing taxid for entity %s part %s\n" % (entityId,pId))
                #
                # Update the sequence ...
                #
                seqIdx=self.__sds.getSequence(seqId0,seqType="auth",partId=1,altId=1,version=curVer)

                ###  JDW reset default comment state ---
                sTupL=[]
                for sT in seqIdx:
                    comment=self.__getAuthSeqComment(sT[3], pD,seqBegMin,seqEndMax)
                    sTupL.append( (sT[0],sT[1],comment,sT[3]))

                #self.__sds.setSequence(seqIdx,seqId0,seqType='auth',partId=partId,altId=1,version=nextVer)
                #### JDW
                self.__sds.setSequence(sTupL,seqId0,seqType='auth',partId=partId,altId=1,version=nextVer)
                #
                # Update features ... Using data from the principal part.
                #
                # fObj=SequenceFeature()
                authSelectId1,authSL1,fObj=self.__getCurrentAuthSelection(entityId,partId=1)
                #
                nL=tU.lookUpSource(taxId=taxId)                   
                if len(nL) > 0:
                    fObj.setSource(organism=nL[0],taxid=taxId)
                else:
                    fObj.setTaxId(taxid=taxId)
                fObj.setAuthPartDetails(partId,seqBeg,seqEnd,seqPartType)
                #
                #  JDW also need to update any other entity-level details --- 
                #
                self.__sds.setFeatureObj(fObj,seqId0,seqType='auth',partId=partId,altId=1,version=nextVer)
                if (self.__verbose):
                    self.__lfh.write("+SequenceDataUpdate.sourceEditFormResponder() entity %s new part %d with version %d\n" % (entityId,partId,nextVer))
        self.__sds.serialize()
        return True


    def __getPartDetails(self,entityId,numExtra=0):
        """ Return dictionary of part boundaries and types ---   pD[pId]=(pId,pSeqBegin,pSeqEnd,pType,taxId)
        """

        pD={}
        seqIds=self.__sds.getGroup(groupId=entityId)
        if len(seqIds)==0:
            return pD

        seqFeature=SequenceFeature()
        seqId0=seqIds[0]

        polymerTypeCode='AA'
        #
        partIdList=self.__sds.getPartIds(seqId0, dataType="sequence", seqType="auth")
        
        for partId in partIdList:
            tD={}
            vL=self.__sds.getVersionIds(seqId0,partId=partId, altId=1, dataType="sequence", seqType="auth")
            if len(vL)>0:
                pfD=self.__sds.getFeature(seqId0,seqType="auth",partId=partId,altId=1,version=vL[0])
                seqFeature.set(pfD)
                pId,pSeqBegin,pSeqEnd,pType=seqFeature.getAuthPartDetails()
                polymerTypeCode=seqFeature.getPolymerType()
                taxId = seqFeature.getSourceTaxId()
                pD[pId]=(pId,pSeqBegin,pSeqEnd,pType,taxId)
                lastPart=pId

        if numExtra > 0:
            pv=self.__placeHolderValue
            pId=lastPart
            for ii in range(0,numExtra):
                pId+=1
                pD[pId]=(pId,pv,pv,pv,pv)
        #
        seqAuthIdx=self.__sds.getSequence(seqId=seqId0,seqType='auth',partId=1, altId=1,version=vL[0])
        r3List=[]
        for sTup in seqAuthIdx:
            (r3,sIdx,comment,idx)=sTup
            r3List.append(r3)
        self.__srd.cnvList3to1WithMods(r3List)
        #
        return pD

    def __writeFasta(self, filePath,sequence,comment="myquery"):
        num_per_line = 60
        l = len(sequence) / num_per_line
        x = len(sequence) % num_per_line
        m = l
        if x:
            m = l + 1

        seq = '>'+str(comment).strip()+'\n'
        for i in range(m):
            n = num_per_line
            if i == l:
                n = x
            seq += sequence[i*num_per_line:i*num_per_line+n]
            if i != (m - 1):
                seq += '\n'
        try:
            ofh=open(filePath,'w')
            ofh.write(seq)
            ofh.close()
            return True
        except:
            if (self.__verbose):
                self.__lfh.write("+RcsbDpUtility.__writeFasta() failed for path %s\n" % filePath)
                traceback.print_exc(file=self.__lfh)                    

        return False

    def __getAuthFeatures(self, entityId, partId):
        """  Get feature data for the author sequence entity sequence.

             Returns a dictionary of features for the more recent sequence version.
        """
        featureD={}
        #
        if (self.__verbose):
            self.__lfh.write("+SequenceDataUpdate.__getAuthFeatures() entityId %s partId %s\n" % (entityId,partId))
        seqFeature=SequenceFeature()        
        
        # get the author sequence identifier
        seqIds=self.__sds.getGroup(groupId=entityId)
        
        if len(seqIds)==0:
            return featureD

        #
        # Collect the feature dictionaries for the parts of each sequence version.
        seqId0=seqIds[0]
        altId=1
        verList=self.__sds.getVersionIds(seqId0, partId=partId, altId=altId, dataType="sequence", seqType='auth')
        if len(verList) == 0:
            return featureD
            
        
        featureD=self.__sds.getFeature(seqId=seqId0,seqType='auth',partId=partId,altId=altId, version=verList[0])

        if (self.__debug):
            kys=featureD.keys()
            for k in sorted( kys ):
                v=featureD[k]
                self.__lfh.write("+SequenceDataUpdate.__getAuthFeatures() %-35s : %s\n" % (k,v))

        return featureD

    def __loadReferenceSequence(self, entityId, partId, refSeqObj,refSeqBeg=None,refSeqEnd=None):
        """ Do what is needed to load the input reference sequence into the sequence data store  - 

            Return: the packed sequence label identifier for the loaded sequence or None

        """
        #
        if (self.__verbose):
            self.__lfh.write("\n+SequenceDataUpdate.__loadReferenceSequence() entityId %s partId %d\n" % (entityId,partId))
        #
        #
        pA=PairwiseAlign()
        pA.setVerbose(self.__verbose)
        #
        seqFeature=SequenceFeature(self.__verbose)                
        seqIdList=self.__sds.getGroup(groupId=entityId)            
        if len(seqIdList) < 1:
            if (self.__verbose):
                self.__lfh.write("+SequenceDataUpdate.__loadReferenceSequence() entity  group %s is empty\n" % entityId)
            return None


        seqId0=seqIdList[0]
        verList=self.__sds.getVersionIds(seqId0, partId=partId, altId=1, dataType="sequence", seqType='auth')
        if (len(verList) < 1 ):
            if (self.__verbose):
                self.__lfh.write("+SequenceDataUpdate.__loadReferenceSequence() entity %s partId %d seqId %s has no sequence data\n" %  (entityId,partId,seqId0)) 
            return None

        verLatest=verList[0]
        seqFeature=self.__sds.getFeatureObj(seqId=seqId0,seqType='auth',partId=partId, altId=1,version=verLatest)
        (pId,seqBeg,seqEnd,seqPartType)=seqFeature.getAuthPartDetails()
        polyTypeCode=seqFeature.getPolymerType()

        if (self.__verbose):
            self.__lfh.write("+SequenceDataUpdate.__loadReferenceSequence() partId %d seqBegin %d seqEnd %d pId %d type %s\n" %  (partId,seqBeg,seqEnd,pId,polyTypeCode))

        #
        altIdList=self.__sds.getAlternativeIds(seqId0, dataType="sequence", seqType="ref",partId=partId)
        lenAltIdList=len(altIdList) 
        #
        if lenAltIdList < 1:
            nextAltId = 1
        else:
            nextAltId = int(altIdList[0]) + 1


        sL=SequenceLabel()
        sL.set(seqType='ref', seqInstId=seqId0, seqPartId=partId, seqAltId=nextAltId, seqVersion=1)
        nextRefLabel=sL.pack()

        if (self.__verbose):
            # JDW -  create new reference label -
            self.__lfh.write("+SequenceDataUpdate.__loadReferenceSequence() partId %d reference sequence count %d next altId %d label %s\n" %  (partId,len(altIdList),nextAltId,nextRefLabel))
        
        #
        # Do we know the boundaries for the new reference sequence?
        #
        if refSeqBeg is None and refSeqEnd is None:
            if self.__verbose:
                self.__lfh.write("+SequenceDataUpdate.__loadReferenceSequence() no valid boundaries detected -  calculating these from alignment\n")

            # need to do alignment select the useful range of the sequence.
            # Set the reference sequence according to part data -- 
            seqAuthIdx=self.__sds.getSequence(seqId=seqId0,seqType='auth',partId=partId, altId=1,version=verLatest)
            r3L=[]
            for tup in seqAuthIdx[seqBeg-1:seqEnd]:
                r3L.append(tup[0])
            pA.setReferenceSequence(r3L, 'auth:' + seqId0 + '_P' + str(partId))
            authSeqLen=len(seqAuthIdx[seqBeg-1:seqEnd])

            testRefIdx=refSeqObj.getSequenceWithIndex(polyTypeCode=polyTypeCode,seqBegin=refSeqBeg,seqEnd=refSeqEnd)
            r3L=[]
            for tup in testRefIdx:
                r3L.append(tup[0])
            pA.addTestSequence(r3L,'ref:'+str(nextAltId)+ '_P' + str(partId) )
            pA.doAlign()
            aL=pA.getAlignment('ref:'+str(nextAltId)+ '_P' + str(partId) )
            alignLength=len(aL)
            #
            # where both are not gap
            idxL=[]
            iRef=1
            for ii,aTup in enumerate(aL,start=1):
                if (self.__verbose):
                    self.__lfh.write("+SequenceDataUpdate.__loadReferenceSequence() alIdx %d  refIdx %d ref %r auth %r\n" % (ii,iRef,aTup[0],aTup[1]))
                if aTup[0] != self.__gapSymbol:                    
                    idxL.append(iRef)
                if (aTup[1] != self.__gapSymbol):
                    iRef+=1

            refSeqBeg=min(idxL)
            refSeqEnd=max(idxL)
            if self.__verbose:
                self.__lfh.write("+SequenceDataUpdate.__loadReferenceSequence() partId %d type %s using aligned refSeqBeg %d refSeqEnd %d\n" %  
                                 (partId,polyTypeCode,refSeqBeg,refSeqEnd))            
            numMatch=0
            numMatchGaps=0
            for aTup in aL:
                if aTup[0] == aTup[1]:
                    numMatch+=1
                if aTup[0] == aTup[1] or aTup[1] == self.__gapSymbol:
                    numMatchGaps+=1                        
            if (self.__verbose):
                for ii,aTup in enumerate(aL,start=1):
                    self.__lfh.write(" -----   %d   %r\n" % (ii,aTup))

        elif refSeqBeg is not None and refSeqEnd is None:
            refSeqEnd=refSeqObj.getSequenceLength()
        #
        # Interpret the range information relative to the full length sequence stored in refSeqObj()
        #
        # Load sequence - 
        #
        seqRefIdx=refSeqObj.getSequenceWithIndex(polyTypeCode=polyTypeCode,seqBegin=refSeqBeg,seqEnd=refSeqEnd)
        if self.__verbose:
            self.__lfh.write("+SequenceDataUpdate.__loadReferenceSequence() partId %d type %s refSeqBeg %d refSeqEnd %d\n" %  (partId,polyTypeCode,refSeqBeg,refSeqEnd))
            self.__lfh.write("+SequenceDataUpdate.__loadReferenceSequence() seqRefIdx %r\n" %  seqRefIdx)
                
        self.__sds.setSequence(seqRefIdx,  seqId0, 'ref',  partId=partId, altId=nextAltId, version=1)
        #
        # Load feature data
        #
        seqFeature.clear()
        #
        dbName,dbCode,dbAccession=refSeqObj.getDatabaseInfo()
        seqFeature.setId(dbName=dbName,dbCode=dbCode, dbAccession=dbAccession)
        #
        taxId=refSeqObj.getTaxId()
        orgName=refSeqObj.getSourceName()
        strain=refSeqObj.getSourceStrain()
        sourceCommonName=refSeqObj.getSourceCommonName()
        seqFeature.setSource(organism=orgName, strain=strain,taxid=taxId,commonName=sourceCommonName)
        #
        proteinName=refSeqObj.getName()
        synonyms=refSeqObj.getSynonyms()
        geneName=refSeqObj.getGeneName()
        seqFeature.setRefSeqNames(proteinName=proteinName, synonyms=synonyms, geneName=geneName)
        ec=refSeqObj.getEnzymeClass()
        keywords=refSeqObj.getKeywords()
        comments=refSeqObj.getComments()
        seqFeature.setRefSeqDetails(enzymeClass=ec, description='', comments=comments, keywords=keywords)
        #
        seqFeature.setAuthRefAlignRange(refMatchBegin=refSeqBeg,refMatchEnd=refSeqEnd)
        seqFeature.setPolymerType(polyTypeCode)      
        #
        # jdw add missing attributes -- 
        
        # less useful --
        seqFeature.setItem('ORG_ORDER_ID', 0)
        sOrder=lenAltIdList+1
        seqFeature.setRefSortOrder(sortIndex=sOrder, sortMetric=self.__defaultInsertSortMetric)
        #
        #seqFeature.setItem('AUTH_REF_SEQ_SIM_BLAST',rD['seq_sim'])
        #
        seqFeature.setAuthPartDetails(partId,seqBeg,seqEnd,seqPartType)
        self.__sds.setFeature( seqFeature.get(),  seqId0, 'ref',  partId=partId, altId=nextAltId, version=1)
        #
        return nextRefLabel



                
if __name__ == '__main__':
    pass
