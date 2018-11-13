##
# File:    SummaryView.py
# Date:    19-Jan-2010
# Updates:
# 09-Feb-2010  jdw Adopt SequenceFeature() to access to feature dictionaries.
# 12-Feb-2010  jdw Move example data loader to the SequenceDataImportExample class.
# 14-Feb-2010  jdw Add statistics updater after primary sequence data load.
# 20-Apr-2010  jdw Ported to module seqmodule
# 27-Apr-2010  jdw Add row-level data dictionary to the summary data object.
# 02-May-2010  jdw Add SequenceSelection()
# 09-Aug-2010  rps __loadSummary() --> Highest numbered version of coordinate sequence to be selected by default in "SELECT" column of UI
# 08-Jan-2013  jdw Select higher version of author sequence in reload operations
# 03-Mar-2013  jdw Refactor -
# 05-Mar-2013  jdw Move default sequence selection from SequenceSelection()
#                  Add default self-reference selections
#                  Implement stored sorting order for reference sequeneces --
#                  Limit reference sequence to maxRefAlign --
# 12-Mar-2013  jdw adjust the assembly of author sequence details -
# 13-Apr-2013  jdw change self reference identifiers
# 15-Nov-2013  jdw major overhaul to support expanded author content display and entity review
# 12-Dec-2013  jdw make default alignment selections.
# 19-Dec-2013  jdw select any added or edited reference sequence by default
# 15-May-2014  jdw add status field for instance id with xyz sequence groups
# 22-May-2014  jdw add method to provide entry details --
# 30-Aug-2017  zf  change self.__reqObj.getSummaryAlignList & self.__reqObj.getSummarySelectList to use latest UI input values
##
"""
Controlling class for the production of data for the summary sequence view.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.08"

import sys
import traceback

from wwpdb.apps.seqmodule.io.SequenceDataStore import SequenceDataStore
from wwpdb.apps.seqmodule.util.SequenceLabel import SequenceLabel, SequenceFeature
from wwpdb.apps.seqmodule.util.SequenceFeatureDepict import SequenceFeatureDepict

from wwpdb.apps.seqmodule.util.SequenceAssign import SequenceAssignArchive, SequenceAssignDepositor, ReferenceSequence
from wwpdb.apps.seqmodule.update.UpdatePolymerEntitySourceDetails import UpdatePolymerEntitySourceDetails
from wwpdb.apps.seqmodule.control.SequenceDataAssemble import SequenceDataAssemble

from wwpdb.io.misc.FormatOut import FormatOut


class SummaryView(object):
    """ Controlling class for the production of data for the summary sequence view.

        Supported operations:

         load             = load summary from current sequence data store using default sequence selection.
         reload           = reload summary using current sequence data and current user selections.

    """

    def __init__(self, reqObj=None, maxRefAlign=100, verbose=False, log=sys.stderr):
        self.__verbose = verbose
        self.__debug = False
        self.__reqObj = reqObj
        self.__lfh = log
        self.__maxRefAlign = maxRefAlign
        #
        # placeholders for sequence identifiers picked on the summary page as selected and/or to be aligned.
        self.__summarySeqAlignList = []
        #
        self.__natureSourceTaxIds = {}
        self.__setup()

    def __setup(self):
        try:
            self.__sessionObj = self.__reqObj.getSessionObj()
            self.__sessionPath = self.__sessionObj.getPath()
            self.__sds = SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            self.__sdu = UpdatePolymerEntitySourceDetails(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            #
            # Selections coming from the the web request --
            #
            self.__summarySeqAlignList = self.__reqObj.getSummaryAlignList(usingRevAllAlignIds=False)
            self.__summarySeqSelectList = self.__reqObj.getSummarySelectList(usingRevSlectIds=False)
            if (self.__verbose):
                self.__lfh.write("+SummaryView.__setup() request input align list  %r\n" % (self.__summarySeqAlignList))
                self.__lfh.write("+SummaryView._setup()  request input select list %r\n" % (self.__summarySeqSelectList))
            #
            # Reset any newly loaded reference sequence id as it is included above --
            self.__reqObj.setNewRefId('')
            #
            # check for the selection of a new reference id --
            #
            if False:
                newId = self.__reqObj.getValue("new-ref-seq-id")
                if len(newId) > 0:
                    swapId = self.__findId(newId, self.__summarySeqSelectList)
                    if swapId is not None:
                        self.__summarySeqSelectList.remove(swapId)
                        self.__summarySeqSelectList.append(newId)
                        self.__reqObj.setNewRefId('')
        except:
            if (self.__verbose):
                self.__lfh.write("+SummaryView.__setup() sessionId %s failed\n" % (self.__sessionObj.getId()))
                traceback.print_exc(file=self.__lfh)

    def __finish(self):
        try:
            # Persist the current summary selection list --
            self.__reqObj.setSummarySelectList(self.__summarySeqSelectList)
            self.__reqObj.setSummaryAlignList(self.__summarySeqAlignList)
            self.__sds.setSelectedIds(idList=self.__summarySeqSelectList)
            self.__sds.serialize()
        except:
            if (self.__verbose):
                self.__lfh.write("+SummaryView.__finish() sessionId %s failed\n" % (self.__sessionObj.getId()))
                traceback.print_exc(file=self.__lfh)

    def __findId(self, qId, idList):
        """  find matching '[self]ref_id_part_' in idList.
        """
        try:
            fL = qId.split('_')
            ss = '_'.join(fL[:3])
            ss1 = 'selfref' + '_'.join(fL[1:3])
            for tId in idList:
                if tId.startswith(ss) or tId.startswith(ss1):
                    return tId
        except:
            pass

        return None

    def getEntryDetails(self, kyList=['STRUCT_TITLE', 'CITATION_TITLE', 'PDB_ID']):
        eD = {}
        for ky in kyList:
            eD[ky] = self.__sds.getEntryDetail(detailKey=ky)
        #
        return eD

    def loadSummary(self, operation='load'):
        """Create the data structure to populate the HTML pages containing alignment summary --

           Options for operation:
           + load   - Performs a default selection of the best matching reference sequences then
                      creates the summary data structure.
           + reload - reloads current data from the sequence and alignment data stores.
        """
        if (self.__verbose):
            self.__lfh.write("+SummaryView.loadSummary() operation %s : sessionId %s\n" % (operation, self.__sessionObj.getId()))
            self.__lfh.write("+SummaryView.loadSummary() request input align list  %r\n" % (self.__summarySeqAlignList))
            self.__lfh.write("+SummaryView.loadSummary() request input select list %r\n" % (self.__summarySeqSelectList))
        if (operation == "load"):
            # Using coming from persisted data only here
            #  -- ignore user selections from web context as these may not exist or may need to be overridden --
            self.__summarySeqSelectList = self.__sds.getSelectedIds()
            if len(self.__summarySeqSelectList) < 1:
                self.__lfh.write("+SummaryView.loadSummary() recomputing default select list\n")
                sda = SequenceDataAssemble(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
                self.__summarySeqSelectList = sda.makeDefaultSelection(sds=self.__sds)

            #
        self.__summarySeqAlignList = self.__summarySeqSelectList

        self.__sdu.updateAuthEntityDetails(selectIdList=self.__summarySeqSelectList)
        self.__sds = SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        return self.__loadSummary(op=operation)

    def getGroupIdList(self):
        return self.__sds.getGroupIds()

    def __loadSummary(self, op='load'):
        """ Assemble the data for sequence summary view using current contents of the session sequence data store.

            Returns:  summaryDataObj [<entity/groupId>]  ['auth'|'xyz'] ->  {data dict}
                                                         ['ref']        -> [(partId,{data dict}),(partId,{data dict}),...]
                                   {data dict}  with keys

                                   dT['ROW_IDS']       =rowIdList
                                   dT['ROW_STATUS']    =rowStatusList
                                   dT['ROW_DATA_DICT'] =rowDataDictList
                                   dT['SELF_REFERENCE_FLAG']=Boolean
        """
        #
        if (self.__verbose):
            self.__lfh.write("+SummaryView.__loadSummary() with sessionId %s\n" % self.__sessionObj.getId())

        #
        seqLabel = SequenceLabel()
        seqFeature = SequenceFeature()
        summaryDataObj = {}

        gIdList = self.__sds.getGroupIds()

        if (self.__verbose):
            self.__lfh.write("+SummaryView.__loadSummary() group list is %r\n" % gIdList)
            self.__lfh.write("+SummaryView.__loadSummary() starting selection list %r\n" % self.__summarySeqSelectList)
            self.__lfh.write("+SummaryView.__loadSummary() starting alignment list %r\n" % self.__summarySeqAlignList)

        for gId in gIdList:
            #
            summaryDataObj[gId] = {}
            #
            dT = self.__buildAuthSection(groupId=gId, op=op)
            summaryDataObj[gId]['auth'] = dT

            rL = self.__buildReferenceSection(groupId=gId)
            summaryDataObj[gId]['ref'] = rL

            summaryDataObj[gId]['xyz'] = self.__buildCoordinateSection(groupId=gId)
        #
        self.__finish()
        #
        warningMsg = ''
        if len(self.__natureSourceTaxIds) > 1:
            warningMsg += "Entry contains multiple natural sources:<br />\n"
            for k,v in self.__natureSourceTaxIds.items():
                if len(v) > 1:
                    warningMsg += "Entities '" + "', '".join(v) + "' have "
                else:
                    warningMsg += "Entity '" + "', '".join(v) + "' has "
                #
                warningMsg += "source taxonomy Id '" + k + "'.<br />\n"
            #
        #
        return summaryDataObj,warningMsg

    def __getAuthFeatures(self, entityId, seqId0, partIdList):
        """  Consolidate the feature data for the author sequences entity sequence respecting
             potential mutiple source details.

             Returns a dictionary of features for each sequence version.
        """
        if (self.__verbose):
            self.__lfh.write("+SummaryView.__getAuthFeatures() entityId %r seqId0 %r partIdList %r\n" % (entityId, seqId0, partIdList))

        spStr = "<br />"
        seqFeature = SequenceFeature()
        #
        #
        depSeqAssignD = self.__sds.getDepositorReferenceAssignments()
        sADep = SequenceAssignDepositor(verbose=self.__verbose, log=self.__lfh)
        sADep.set(depSeqAssignD)
        refSeqAssignL = sADep.getReferenceList(entityId)

        if (self.__debug):
            self.__lfh.write("+SummaryView.__getAuthFeatures() for entityId %r reference count is %d\n" % (entityId, sADep.getReferenceCount(entityId)))
            sADep.printIt(log=self.__lfh)
            for ii, rsa in enumerate(refSeqAssignL):
                self.__lfh.write("+SummaryView.__getAuthFeatures() depositor reference  %d\n" % (ii + 1))
                rsa.printIt(self.__lfh)
        #
        depText = ''
        for rsa in refSeqAssignL:
            dbName, dbCode, dbAccession = rsa.getDbReference()
            if ((len(dbAccession) > 0 and (dbAccession not in ['.', '?'])) or (len(dbCode) > 0 and (dbCode not in ['.', '?']))):
                tRef = "Depositor reference: %s %s %s" % (dbName, dbAccession, dbCode)
                sp = spStr if len(depText) > 0 else ''
                depText += (sp + tRef)

            refDetails = rsa.getDetails()
            if ((len(refDetails) > 0 and refDetails not in ['.', '?'])):
                tDetails = "Depositor details: %s" % refDetails
                sp = spStr if len(depText) > 0 else ''
                depText += (sp + tDetails)
        #
        if (self.__verbose):
            self.__lfh.write("+SummaryView.__getAuthFeatures() depositor info %s\n" % depText)
        #
        # Collect the feature dictionaries for the parts of each sequence version.
        #
        featureD = {}
        altId = 1
        for partId in partIdList:
            verList = self.__sds.getVersionIds(seqId0, partId=partId, altId=altId, dataType="sequence", seqType='auth')
            for ver in verList:
                seqAuthFD = self.__sds.getFeature(seqId=seqId0, seqType='auth', partId=partId, altId=altId, version=ver)
                if ver not in featureD:
                    featureD[ver] = []
                featureD[ver].append(seqAuthFD)

        if (self.__debug):
            for k, v in featureD.items():
                self.__lfh.write("+SummaryView.__getAuthFeatures() ver %d  fD %r\n" % (k, v))

        #
        # prepare marked-up details strings for each version consolidating part details.
        #  detailsD[ver] = d2   d2[source]
        #
        detailsTemplate = "Part %s: %s (%4s - %4s)"
        detailsD = {}
        for ver, fDList in featureD.items():
            d2 = {}
            for k in ['taxId', 'details', 'sourceAndStrain', 'description']:
                d2[k] = ''
            for partId, fD in enumerate(fDList, start=1):
                seqFeature.set(fD)
                #
                tDescription = seqFeature.getEntityDescription()
                tEc = seqFeature.getEntityEnzymeClass()
                tFrag = seqFeature.getEntityFragmentDetails()
                tMutation = seqFeature.getEntityMutationDetails()
                tDetails = seqFeature.getEntityDetails()
                ##
                tName = seqFeature.getEntitySynonyms()
                d2['description'] = ''
                if len(tName) > 1:
                    d2['description'] = "Name: " + tName + spStr
                if len(tDescription) > 1:
                    d2['description'] += "Description: " + tDescription + spStr
                if len(tEc) > 1:
                    d2['description'] += "EC: " + tEc + spStr
                if len(tFrag) > 1:
                    d2['description'] += "Fragment: " + tFrag + spStr
                if len(tMutation) > 1:
                    d2['description'] += "Mutation: " + tMutation + spStr
                if len(tDetails) > 1:
                    d2['description'] += "Entity details: " + tDetails + spStr

                if len(d2['description']) > 0:
                    d2['description'] += depText

                sp = '' if partId == 1 else spStr
                taxId = seqFeature.getSourceTaxId()
                d2['taxId'] += sp + taxId

                #
                pId, seqNumBeg, seqNumEnd, seqPartType = seqFeature.getAuthPartDetails()
                d2['details'] += sp + detailsTemplate % (partId, seqPartType, seqNumBeg, seqNumEnd)

                authOrg = seqFeature.getSourceOrganism()
                authStrain = seqFeature.getSourceStrain()
                if len(authStrain) > 1:
                    d2['sourceAndStrain'] += sp + authOrg + '/' + authStrain
                else:
                    d2['sourceAndStrain'] += sp + authOrg
                #
            detailsD[ver] = d2
        #
        if self.__debug:
            self.__lfh.write("SummaryView.__getAuthFeatures() detailsD %r\n" % detailsD.items())
            for k, v in detailsD.items():
                self.__lfh.write("+SummaryView.__getAuthFeatures() ver %d  detailsD %r\n" % (k, v))
        #
        return detailsD

    def __buildAuthSection(self, groupId=None, op='load'):
        """ Assemble the data content for the author sequence summary view.

            Returns: summaryDataObject)

            ** Always show a single full sequence for each entity and store multi-part information as details

        """

        seqLabel = SequenceLabel()
        seqFeature = SequenceFeature()

        taxD = {}
        rowIdList = []
        rowDataList = []
        rowStatusList = []
        rowDataDictList = []
        #
        seqIdList = self.__sds.getGroup(groupId)
        #
        # Each author entity sequence is shown once and is labeled by its first instance in the group.
        #
        #  ###JDW  CHANGE
        # seqId0=seqIdList[0]
        seqId0 = groupId

        if self.__verbose:
            self.__lfh.write("SummaryView.__buildAuthSection() groupId %r op %s \n" % (groupId, op))
        #
        altId = 1
        partIdList = self.__sds.getPartIds(seqId0, dataType='sequence', seqType='auth')
        partId0 = partIdList[0]
        #detailsD=self.__getAuthFeatures(groupId, seqId0, partIdList)
        detailsD = self.__getAuthFeaturesAll(groupId, seqId0, partIdList)

        verList = self.__sds.getVersionIds(seqId0, partId=partId0, altId=altId, dataType="sequence", seqType='auth')

        for ver in verList:
            seqAuth = self.__sds.getSequence(seqId=seqId0, seqType='auth', partId=partId0, altId=altId, version=ver)
            seqLabel.set(seqType='auth', seqInstId=seqId0, seqPartId=partId0, seqAltId=altId, seqVersion=ver)
            seqAuthId = seqLabel.pack()
            rowIdList.append(seqAuthId)
            #
            isSelected = seqAuthId in self.__summarySeqSelectList
            isAligned = seqAuthId in self.__summarySeqAlignList
            #
            if op == 'reload':
                isSelected = (ver == max(verList))
                isAligned = (ver == max(verList))
            #
            rowStatusList.append((isSelected, isAligned))
            #
            rowDataDict = {}
            rowDataDict['ROW_VERSION'] = ver
            rowDataDict['ROW_SEQ_LENGTH'] = len(seqAuth)
            rowDataDict['ROW_DETAIL_STRING'] = ''
            rowDataDict.update(detailsD[ver])
            rowDataDictList.append(rowDataDict)

        #
        dT = {}
        dT['ROW_IDS'] = rowIdList
        dT['ROW_STATUS'] = rowStatusList
        dT['ROW_DATA_DICT'] = rowDataDictList
        dT['ENTITY_NUM_PARTS'] = len(partIdList)

        return dT

    def __buildCoordinateSection(self, groupId=None):
        """ Assemble the data content for the coordinate sequence summary view.

            Returns: summaryDataObject

        """
        seqLabel = SequenceLabel()
        seqFeature = SequenceFeature()

        #
        # For each coordinate sequence -
        #
        seqIdList = self.__sds.getGroup(groupId)

        rowIdList = []
        rowStatusList = []
        rowDataDictList = []
        #
        partId = 1
        altId = 1
        for seqId in seqIdList:
            verList = self.__sds.getVersionIds(seqId=seqId, partId=partId, altId=altId, dataType="sequence", seqType='xyz')
            if len(verList) < 1:
                continue
            maxVrsnNum = verList[0]
            seqLabel.set(seqType='xyz', seqInstId=seqId, seqPartId=partId, seqAltId=altId, seqVersion=maxVrsnNum)
            maxSeqXyzId = seqLabel.pack()

            # select the highest version sequence
            if not (maxSeqXyzId in self.__summarySeqSelectList):
                self.__summarySeqSelectList.append(maxSeqXyzId)
                self.__summarySeqAlignList.append(maxSeqXyzId)

            for ver in verList:
                seqXyz = self.__sds.getSequence(seqId, 'xyz', partId=partId, altId=altId, version=ver)
                seqXyzFD = self.__sds.getFeature(seqId, 'xyz', partId=partId, altId=altId, version=ver)
                seqLabel.set(seqType='xyz', seqInstId=seqId, seqPartId=partId, seqAltId=altId, seqVersion=ver)
                seqXyzId = seqLabel.pack()
                #
                seqFeature.set(seqXyzFD)
                sfd = SequenceFeatureDepict(sfObj=seqFeature, verbose=self.__verbose, log=self.__lfh)
                detailString = sfd.markupXyzFeatures()
                #
                rowLabel = seqId
                rowIdList.append(seqXyzId)
                #
                # unselect any selected lower version for this sequence -
                if (ver < maxVrsnNum) and (seqXyzId in self.__summarySeqSelectList):
                    self.__summarySeqSelectList.remove(seqXyzId)
                    if (seqXyzId in self.__summarySeqAlignList):
                        self.__summarySeqAlignList.remove(seqXyzId)
                #
                isSelected = seqXyzId in self.__summarySeqSelectList
                isAligned = seqXyzId in self.__summarySeqAlignList
                # JDW add instance here to provide a selection group
                if (ver == maxVrsnNum):
                    rowStatusList.append((isSelected, isAligned, seqId, 'maxver'))
                else:
                    rowStatusList.append((isSelected, isAligned, seqId, ''))
                #
                rowDataDict = {}
                rowDataDict['ROW_ID_CODE'] = rowLabel
                rowDataDict['ROW_VERSION'] = ver
                rowDataDict['ROW_SEQ_LENGTH'] = len(seqXyz)
                rowDataDict['ROW_DETAIL_STRING'] = detailString
                rowDataDict.update(seqXyzFD)
                rowDataDictList.append(rowDataDict)
        #
        dT = {}
        dT['ROW_IDS'] = rowIdList
        dT['ROW_STATUS'] = rowStatusList
        dT['ROW_DATA_DICT'] = rowDataDictList
        return dT

    def __buildReferenceSection(self, groupId=None):
        """ Assemble the data content for the reference sequence summary view.

            Returns: [(partId,summaryDataObject),]

                   where:
                      partId is the integer identifier for pieces of a multipart entity sequence
                      summaryDataObject contains the 'ref' sequence summary data for the input groupId.
        """
        rL = []
        seqLabel = SequenceLabel()
        seqFeature = SequenceFeature()
        #
        # use the leading sequence Id in the group as the id for the reference -

        seqIdList = self.__sds.getGroup(groupId)
        #
        # JDW ### CHANGE
        # seqIdRef=seqIdList[0]
        seqIdRef = groupId

        #
        partIdList = self.__sds.getPartIds(seqIdRef, dataType='sequence', seqType='auth')
        for partId in partIdList:
            #
            authVerList = self.__sds.getVersionIds(seqId=seqIdRef, partId=partId, altId=1, dataType="feature", seqType='auth')
            skipPart = False
            if (len(authVerList) > 0):
                authFObj = self.__sds.getFeatureObj(seqIdRef, 'auth', partId=partId, altId=1, version=authVerList[0])
                sourceInfo = authFObj.getEntitySourceMethod().upper()
                TaxId = authFObj.getSourceTaxId()
                if (sourceInfo == 'NAT') and TaxId:
                     if self.__natureSourceTaxIds.has_key(TaxId):
                         self.__natureSourceTaxIds[TaxId].append(groupId)
                     else:
                         self.__natureSourceTaxIds[TaxId] = [ groupId ]
                     #
                #
                authPartId, authSeqBeg, authSeqEnd, authSeqPartType = authFObj.getAuthPartDetails()
                skipPart = len(authSeqPartType) > 1 and authSeqPartType.lower() != 'biological sequence'
            if skipPart:
                self.__lfh.write("+SummaryView.__buildReferenceSection() entity %s skipping part %s type %s\n" % (groupId, partId, authSeqPartType))
                continue

            rowIdList = []
            rowDataList = []
            rowStatusList = []
            rowDataDictList = []
            # List of reference sequences for this group only for the leading sequence -
            #
            altIdList = self.__sds.getAlternativeIds(seqIdRef, dataType="sequence", seqType="ref", partId=partId)

            srId = 'selfref_' + str(groupId) + '_' + str(partId)
            if srId in self.__summarySeqSelectList:
                if (self.__verbose):
                    self.__lfh.write("+SummaryView.__buildReferenceSection() using self-refereence for entity %s part %s\n" % (groupId, partId))
                selfRefFlag = True
            else:
                selfRefFlag = False
            #
            for altId in altIdList[:self.__maxRefAlign]:
                #
                verList = self.__sds.getVersionIds(seqId=seqIdRef, partId=partId, altId=altId, dataType="sequence", seqType='ref')
                for ver in verList:
                    #
                    seqRefFD = self.__sds.getFeature(seqIdRef, 'ref', partId=partId, altId=altId, version=ver)
                    seqFeature.set(seqRefFD)
                    seqLabel.set(seqType='ref', seqInstId=seqIdRef, seqPartId=partId, seqAltId=altId, seqVersion=ver)
                    seqRefId = seqLabel.pack()

                    isSelected = seqRefId in self.__summarySeqSelectList
                    isAligned = seqRefId in self.__summarySeqAlignList
                    rowStatusList.append((isSelected, isAligned))
                    #
                    sfd = SequenceFeatureDepict(sfObj=seqFeature, verbose=self.__verbose, log=self.__lfh)
                    detailString = sfd.markupReferenceAlignmentFeatures()
                    featureString = sfd.markupReferenceFeatures()
                    rowLabel = sfd.markupDatabaseReferenceWithUrl(altId)
                    rowIdList.append(seqRefId)

                    lengthRefSeq = seqFeature.getMatchLength()
                    authSimWithGaps = seqFeature.getAuthRefSimWithGaps()
                    #
                    #
                    rowDataDict = {}
                    rowDataDict['ROW_ID_CODE'] = rowLabel
                    rowDataDict['ROW_VERSION'] = ver
                    rowDataDict['ROW_SEQ_LENGTH'] = lengthRefSeq
                    rowDataDict['ROW_AUTH_REF_SIM'] = "%.3f" % authSimWithGaps
                    rowDataDict['ROW_DETAIL_STRING'] = detailString
                    rowDataDict['ROW_FEATURE_STRING'] = featureString
                    rowDataDict['ROW_IS_SELECTED'] = isSelected
                    rowDataDict['ROW_IS_ALIGNED'] = isAligned
                    rowDataDict.update(seqRefFD)
                    rowDataDictList.append(rowDataDict)
            #
            dT = {}
            dT['ROW_IDS'] = rowIdList
            dT['ROW_STATUS'] = rowStatusList
            dT['ROW_DATA_DICT'] = rowDataDictList
            dT['SELF_REFERENCE_FLAG'] = selfRefFlag
            rL.append((partId, dT))
        #
        return rL

    def __getAuthFeaturesAll(self, entityId, seqId0, partIdList):
        """  Consolidate the feature data for the author sequences entity sequence respecting
             potential mutiple source details.

             Original and current features are assembled.

             Returns a dictionaries of features for each sequence version.

             Details[ver][partId]['current']={'partdetails','sourceAndStrain','description','hostorg'}
             Details[ver][partId]['author']={'partdetails','sourceAndStrain','description','hostorg'}
        """
        if (self.__verbose):
            self.__lfh.write("+SummaryView.__getAuthFeaturesAll() entityId %r seqId0 %r partIdList %r\n" % (entityId, seqId0, partIdList))

        spStr = "<br />"
        #
        authRefAssignText = self.__markupAuthorAssignments(entityId)
        #
        #
        altId = 1
        detailsD = {}
        for partId in partIdList:
            verList = self.__sds.getVersionIds(seqId0, partId=partId, altId=altId, dataType="sequence", seqType='auth')
            for ver in verList:
                if ver not in detailsD:
                    detailsD[ver] = {}

                sfObj = self.__sds.getFeatureObj(seqId=seqId0, seqType='auth', partId=partId, altId=altId, version=ver)
                sfd = SequenceFeatureDepict(sfObj=sfObj, verbose=self.__verbose, log=self.__lfh)

                if partId not in detailsD[ver]:
                    detailsD[ver][partId] = {}

                detailsD[ver][partId]['current'] = sfd.markupCurrentEntityDetails()
                detailsD[ver][partId]['author'] = sfd.markupAuthorEntityDetails()

                if len(authRefAssignText) > 0:
                    detailsD[ver][partId]['author']['description'] += authRefAssignText

        for ver in detailsD.keys():
            np = 0
            for partId in partIdList:
                if partId in detailsD[ver]:
                    np += 1
            detailsD[ver]['numparts'] = np
        if self.__debug:
            self.__lfh.write("+SummaryView.__getAuthFeaturesAll() detailsD entity %r seqId %r partIdList %r\n" % (entityId, seqId0, partIdList))
            fOut = FormatOut()
            fOut.autoFormat("Summary Entity Details", detailsD, 3, 3)
            fOut.writeStream(self.__lfh)
            fOut.clear()
        #
        return detailsD

    def __markupAuthorAssignments(self, entityId):
        """  Markup the author provided reference assignments -
        """
        spStr = "<br />"
        authRefAssignText = ''
        depSeqAssignD = self.__sds.getDepositorReferenceAssignments()
        sADep = SequenceAssignDepositor(verbose=self.__verbose, log=self.__lfh)
        sADep.set(depSeqAssignD)
        refSeqAssignL = sADep.getReferenceList(entityId)

        if (self.__debug):
            self.__lfh.write("+SummaryView.__markupAuthorAssignments() for entityId %r author reference count is %d\n" % (entityId, sADep.getReferenceCount(entityId)))
            sADep.printIt(log=self.__lfh)
            for ii, rsa in enumerate(refSeqAssignL):
                self.__lfh.write("+SummaryView.__markupAuthorAssignments() depositor reference  %d\n" % (ii + 1))
                rsa.printIt(self.__lfh)
        #

        for rsa in refSeqAssignL:
            dbName, dbCode, dbAccession = rsa.getDbReference()
            if ((len(dbAccession) > 0 and (dbAccession not in ['.', '?'])) or (len(dbCode) > 0 and (dbCode not in ['.', '?']))):
                tRef = "Depositor reference: %s %s %s" % (dbName, dbAccession, dbCode)
                sp = spStr if len(authRefAssignText) > 0 else ''
                authRefAssignText += (sp + tRef)

            refDetails = rsa.getDetails()
            if ((len(refDetails) > 0 and refDetails not in ['.', '?'])):
                tDetails = "Depositor details: %s" % refDetails
                sp = spStr if len(authRefAssignText) > 0 else ''
                authRefAssignText += (sp + tDetails)
        #
        if (self.__verbose):
            self.__lfh.write("+SummaryView.__getAuthFeatures() depositor reference assignments for entity %s : %s\n" % (entityId, authRefAssignText))
        #
        return authRefAssignText


if __name__ == '__main__':
    pass
