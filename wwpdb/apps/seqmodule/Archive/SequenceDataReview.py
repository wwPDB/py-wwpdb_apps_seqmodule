##
# File:    SequenceDataReview.py
# Date:    18-Sep-2013
#
# Updates:
# 
##
"""
Prepare entity level sequence annotation review presentation with selected editable annotations.

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
from wwpdb.apps.seqmodule.util.SequenceLabel                import SequenceLabel, SequenceFeature
from wwpdb.apps.seqmodule.util.SequenceAssign               import SequenceAssignArchive,SequenceAssignDepositor,ReferenceSequence
from wwpdb.apps.seqmodule.util.SequenceReferenceData        import SequenceReferenceData

class SequenceDataReview(object):
    """ Prepare entity level sequence annotation review presentation with selected editable annotations.

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
        except:
            if (self.__verbose):
                self.__lfh.write("+SequenceDataReview.__setup() sessionId %s failed\n" % (self.__sessionObj.getId()))                                         

    def __getInstanceList(self,entityId='1'):
        seqIdList=self.__sds.getGroup(entityId)
        if (self.__verbose):
            self.__lfh.write("+SequenceDataReview.__getInstanceList() entityId %s instance list %s\n" % (entityId,seqIdList))
        return seqIdList
        
    def __getFeatureObj(self,seqType='ref'):
        sL=SequenceLabel()
        for sId in self.__selectIdList:
            if sId.startswith(seqType):
                sL.unpack(sId)
                seqType0,seqInstId,seqPartId,seqAltId,seqVersion=sL.get()
                fObj=self.__sds.getFeatureObj(seqInstId,seqType=seqType,partId=seqPartId,altId=seqAltId,version=seqVersion)
                return fObj

    def __getDepositorAssignments(self):
        """  Get the depositor and prior sequence reference assignments.  
        """
        entityIdList=self.__sds.getGroupIds()

        seqAssignD=self.__sds.getReferenceAssignments()
        sA=SequenceAssignArchive(verbose=self.__verbose,log=self.__lfh)
        sA.set(seqAssignD)
        if (self.__debug):
            sA.printIt(log=self.__lfh)
            for entityId in entityIdList:
                nRef=sA.getReferenceCount(entityId=entityId)
                self.__lfh.write("+SequenceDataReview.dump() entityId %s archive assignment count %d\n" % (entityId,nRef))
                if nRef > 0:
                    refL=sA.getReferenceList(entityId=entityId)
                    for ref in refL:
                        ref.printIt(self.__lfh)
        #
        depSeqAssignD=self.__sds.getDepositorReferenceAssignments()
        sADep=SequenceAssignDepositor(verbose=self.__verbose,log=self.__lfh)
        sADep.set(depSeqAssignD)
        if (self.__debug):
            sADep.printIt(log=self.__lfh)
            for entityId in entityIdList:
                nRef=sADep.getReferenceCount(entityId=entityId)
                self.__lfh.write("+SequenceDataReview.dump() entityId %s depositor assignment count %d\n" % (entityId,nRef))
                if nRef > 0:
                    refL=sADep.getReferenceList(entityId=entityId)
                    for ref in refL:
                        ref.printIt(self.__lfh)
        return sA,sADep



    def makeEntityReviewForm(self,entityId,entryId=''):
        """    Return a preliminary form to review entity annotations.

        """
        #
        form_template='''
        <div id="sectentityreview">
        <h3>Entity Review Form for Entry %(entryid)s Entity %(entityid)s Part %(partid)s</h3>
        <form name="formentityreview_%(partid)s" id="formentityreview_%(partid)s" action="/service/sequence_editor/respond_form/entityreview" method="post" class="review_ajaxform">
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
            <td><span id="dbname"      class="ief greyedout" data-ief-edittype="select" 
                                       data-ief-selectvalues='[{"value":"UNP","label":"UNP","selected":false},{"value":"GB","label":"GB","selected":false}]'>%(dbname)s</span></td>
            <td><span id="dbaccession" class="ief greyedout">%(dbaccession)s</span></td>
            <td><span id="dbseqbegin"  class="ief greyedout">%(dbseqbegin)s</span></td>
            <td><span id="dbseqend"    class="ief greyedout">%(dbseqend)s</span></td>
            </tr>

            </table>
            <input type="submit" name="submit" value="Submit" class="disableonclick submitparentform" />
            <!-- <input type="reset" name="reset" value="Reset" /> --> 
        </form>
       </div>
        '''
        #
        sA,sADep=self.__getDepositorAssignments()
        authFObj=self.__getFeatureObj(seqType='auth')
        refFObj=self.__getFeatureObj(seqType='ref')
        instanceList=self.__getInstanceList(entityId=entityId)

        #
        
        dbName,dbAccession,dbSeqBegin,dbSeqEnd=self.__getCurrentRefDetails()
        partId=1
        pD={}
        pD['entryid']=entryId
        pD['sessionid']=self.__sessionId
        pD['partid']=partId
        pD['entityid']=entityId        
        pD['dbname']=dbName           if dbName      is not None else self.__placeHolderValue
        pD['dbaccession']=dbAccession if dbAccession is not None else self.__placeHolderValue
        pD['dbseqbegin']=dbSeqBegin   if dbSeqBegin  is not None else self.__placeHolderValue
        pD['dbseqend']=dbSeqEnd       if dbSeqEnd    is not None else self.__placeHolderValue
        rD = {}        
        rD['htmlcontent']=(form_template % pD)
        #
        return rD

    def entityReviewResponder(self):
        """  Update the entity annotation from form data - 
        
             Form data encoded in the input request object -- 

             Return True
        """
        try:
            partId=int(str(self.__reqObj.getValue("partid")))
        except:
            self.__lfh.write("+SequenceDataReview.entityReviewResponder() failing\n")
            traceback.print_exc(file=self.__lfh)            

        return True
                
if __name__ == '__main__':
    pass
