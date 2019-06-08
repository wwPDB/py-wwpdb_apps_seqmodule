##
# File:    UpdatePolymerEntityReference.py
# Date:    22-Mar-2014
#
# Updates:
# 04-July-2014  jdw Fix GB reference sequence update
# 29-July-2014  jdw Extend for sequence isoforms
##
"""
Utilities for updating reference sequence of matching residue ranges

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.08"

import os
import os.path
import sys
import traceback

from wwpdb.apps.seqmodule.io.SequenceDataStore import SequenceDataStore
from wwpdb.apps.seqmodule.io.ReferenceSequenceUtils import ReferenceSequenceUtils
from wwpdb.apps.seqmodule.util.SequenceLabel import SequenceLabel, SequenceFeature, SequenceFeatureMap
from wwpdb.apps.seqmodule.util.SequenceAssign import SequenceAssignArchive, SequenceAssignDepositor, ReferenceSequence
from wwpdb.apps.seqmodule.util.SequenceReferenceData import SequenceReferenceData
from wwpdb.apps.seqmodule.align.AlignmentStatistics import AlignmentStatistics
from wwpdb.utils.align.alignlib import PairwiseAlign


class UpdatePolymerEntityReference(object):
    """ Utilities for updating reference sequence of matching residue ranges

    """

    def __init__(self, reqObj=None, verbose=False, log=sys.stderr):
        self.__verbose = verbose
        self.__debug = False
        self.__reqObj = reqObj
        self.__lfh = log
        self.__defaultInsertSortMetric = 100000
        #
        self.__srd = SequenceReferenceData(verbose=self.__verbose, log=self.__lfh)
        self.__gapSymbol = self.__srd.getGapSymbol()
        self.__maxRefAlign = 100
        #
        self.__seqFetchError = ''
        #
        self.__setup()

    def __setup(self):
        try:
            self.__placeHolderValue = "click-to-edit"
            self.__siteId = self.__reqObj.getValue("WWPDB_SITE_ID")
            self.__sessionId = self.__reqObj.getSessionId()
            self.__sessionObj = self.__reqObj.getSessionObj()
            self.__sessionPath = self.__sessionObj.getPath()
            #
            self.__selectIdList = self.__reqObj.getSummarySelectList()
            self.__sds = SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            #
        except:
            if (self.__verbose):
                self.__lfh.write("+UpdatePolymerEntityReference.__setup() sessionId %s failed\n" % (self.__sessionObj.getId()))

    def seqDbRefFormResponder(self):
        """  Update the sequence data store using data sequence database reference

             Form data encoded in the input request object --

             Return the packed sequence label of the added reference sequence or None
        """
        packedSeqLabel = None
        try:
            partId = int(str(self.__reqObj.getValue("partid")))
            entityId = self.__reqObj.getValue("entityid")
            dbName = self.__reqObj.getValue("dbname")
            dbAccession = self.__reqObj.getValue("dbaccession")

            dbSeqBegin = self.__reqObj.getValue("dbseqbegin")
            dbSeqEnd = self.__reqObj.getValue("dbseqend")
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

            dbIsoform = ''
            if dbName in ['UNP', 'SP', 'TR']:
                tL = dbAccession.split('-')
                if len(tL) > 1:
                    dbIsoform = dbAccession
                    dbAccession = tL[0]

            packedSeqLabel = self.addReferenceSequence(entityId, partId, dbName, dbAccession, dbIsoform=dbIsoform, refSeqBeg=dbSeqBegin, refSeqEnd=dbSeqEnd)
            if packedSeqLabel:
                self.__reqObj.setNewRefId(packedSeqLabel)
            #
        except:
            self.__lfh.write("+UpdatePolymerEntityReference.seqDbRefResponder() failing\n")
            traceback.print_exc(file=self.__lfh)

        return packedSeqLabel

    def makeSeqdbrefEditForm(self, entityId, partId=1, entryId=''):
        """    Return a preliminary form to input sequence database references.

        """
        #
        form_template = '''
        <div id="sectseqdbref">
        <h3>Reference Sequence Database Data Form for Entry %(entryid)s Entity %(entityid)s Part %(partid)s</h3>
        <form name="formseqdbref_%(partid)s" id="formseqdbref_%(partid)s" action="/service/sequence_editor/respond_form/seqdbref" method="post" class="auth_ajaxform">
            <input type="hidden" name="sessionid" value="%(sessionid)s" />
            <input type="hidden" name="partid" value="%(partid)s" />
            <input type="hidden" name="entityid" value="%(entityid)s" />
            <table>
            <tr>
               <th>Reference Resource</th>
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
        dbName, dbCode, dbAccession, dbIsoform, dbSeqBegin, dbSeqEnd = self.__getCurrentRefDetails(entityId, partId=int(partId))
        #

        pD = {}
        pD['entryid'] = entryId
        pD['sessionid'] = self.__sessionId
        pD['partid'] = partId
        pD['entityid'] = entityId
        #
        pD['dbname'] = self.__srd.convertDbNameToResource(dbName)
        #
        if dbIsoform is not None and len(dbIsoform) > 0:
            pD['dbaccession'] = dbIsoform
        else:
            pD['dbaccession'] = dbAccession

        pD['dbseqbegin'] = dbSeqBegin
        pD['dbseqend'] = dbSeqEnd
        #
        itemList = ['dbname', 'dbaccession', 'dbseqbegin', 'dbseqend']
        for item in itemList:
            kyCss = item + "_css"
            if pD[item] is None or len(str(pD[item])) < 1:
                pD[item] = self.__placeHolderValue
                pD[kyCss] = "greyedout"
            else:
                pD[kyCss] = ""

        rD = {}
        rD['htmlcontent'] = (form_template % pD)
        #
        return rD

    def __getCurrentRefSelection(self, entityId, partId=1):
        try:
            seqIdList = self.__sds.getGroup(entityId)
            # JDW CHANGE
            # instanceId=seqIdList[0]
            instanceId = entityId
            sL = SequenceLabel()
            selfRefTarget = 'selfref_%s_%d' % (entityId, partId)
            for sId in self.__selectIdList:
                if sId.startswith('ref'):
                    sL.unpack(sId)
                    seqType, seqInstId, seqPartId, seqAltId, seqVersion = sL.get()
                    if instanceId == seqInstId and partId == seqPartId:
                        fObj = self.__sds.getFeatureObj(seqInstId, seqType="ref", partId=seqPartId, altId=seqAltId, version=seqVersion)
                        self.__lfh.write("+UpdatePolymerEntityReference._getCurrentRefSelection() returns: entity %r instance %r partId %r altId %r version %r\n" %
                                         (entityId, seqInstId, seqPartId, seqAltId, seqVersion))
                        return sId, sL, fObj
                elif sId.startswith(selfRefTarget):
                    # JDW JDW
                    sL.set(seqType='ref', seqInstId=instanceId, seqPartId=partId, seqAltId=1, seqVersion=1)
                    fObj = SequenceFeature()
                    return sId, sL, fObj
        except:
            self.__lfh.write("+UpdatePolymerEntityReference._getCurrentRefSelection() failed for selectList %r entityId %s partId %r\n" % (self.__selectIdList, entityId, partId))
            traceback.print_exc(file=self.__lfh)

        self.__lfh.write(
            "+UpdatePolymerEntityReference._getCurrentRefSelection() no return for entityId %r instanceId %r partId %r selectIdList %r\n" %
            (entityId, instanceId, partId, self.__selectIdList))
        return None, None, None

    def __getCurrentRefDetails(self, entityId, partId=1):
        seqIdList = self.__sds.getGroup(entityId)
        # JDW CHANGE
        # seqId=seqIdList[0]
        seqId = entityId
        self.__lfh.write("+UpdatePolymerEntityReference.__getCurrentRefDetails() entityId %r partId %r seqId %r\n" % (entityId, partId, seqId))
        self.__lfh.write("+UpdatePolymerEntityReference.__getCurrentRefDetails() selectIdList %r\n" % (self.__selectIdList))
        sL = SequenceLabel()
        for sId in self.__selectIdList:
            if sId.startswith('ref'):
                sL.unpack(sId)
                seqType, seqInstId, seqPartId, seqAltId, seqVersion = sL.get()
                self.__lfh.write("+UpdatePolymerEntityReference.__getCurrentRefDetails() testing seqInstId %r seqPartId %r\n" % (seqInstId, seqPartId))
                if seqId == seqInstId and partId == seqPartId:
                    fObj = self.__sds.getFeatureObj(seqInstId, seqType="ref", partId=seqPartId, altId=seqAltId, version=seqVersion)
                    fD = fObj.get()
                    # return  self.__srd.convertDbNameToResource(fD['DB_NAME']),fD['DB_CODE'],fD['DB_ACCESSION'],fD['REF_MATCH_BEGIN'],fD['REF_MATCH_END']
                    return fD['DB_NAME'], fD['DB_CODE'], fD['DB_ACCESSION'], fD['DB_ISOFORM'], fD['REF_MATCH_BEGIN'], fD['REF_MATCH_END']

        return None, None, None, None, None, None

    def fetchArchiveReferenceSequences(self):
        entityIdList = self.__sds.getGroupIds()

        seqAssignD = self.__sds.getReferenceAssignments()
        sA = SequenceAssignArchive(verbose=self.__verbose, log=self.__lfh)
        sA.set(seqAssignD)

        for entityId in entityIdList:
            nRef = sA.getReferenceCount(entityId=entityId)
            self.__lfh.write("+UpdatePolymerEntityReference.fetchArchiveReferenceSequences() entityId %s archive assignment count %d\n" % (entityId, nRef))
            if nRef > 0:
                partD = self.__getPartDetails(entityId=entityId)
                refL = sA.getReferenceList(entityId=entityId)
                for ref in refL:
                    dbName, dbCode, dbAccession = ref.getDbReference()
                    seq_align_beg, seq_align_end = ref.getSeqAlignRange()
                    #
                    # Match the residue ranges to entity parts  -
                    #
                    pId = self.__assignPart(partD=partD, seqBegin=seq_align_beg, seqEnd=seq_align_end)
                    if pId > 0:
                        self.addReferenceSequence(entityId=entityId, partId=pId, dbName=dbName, dbAccession=dbAccession, refSeqBeg=seq_align_beg, refSeqEnd=seq_align_end)
                    else:
                        self.__lfh.write("+UpdatePolymerEntityReference.fetchArchiveReferenceSequences() entityId %s  no part assigned for residue range %d - %d\n"
                                         % (entityId, seq_align_beg, seq_align_end))
        if self.__debug:
            self.__sds.dump(self.__lfh)
        return True

    def addReferenceSequence(self, entityId, partId, dbName, dbAccession, dbIsoform=None, refSeqBeg=None, refSeqEnd=None):
        if (self.__verbose):
            self.__lfh.write("+UpdatePolymerEntityReference.addReferenceSequences() starting with entityId %r partId %r dbName %r dbAccession %r dbIsoform %r refSeqBeg %r refSeqEnd %r\n" %
                             (entityId, partId, dbName, dbAccession, dbIsoform, refSeqBeg, refSeqEnd))
        fD = self.__getAuthFeatures(entityId=entityId, partId=partId)
        seqFeature = SequenceFeature()
        seqFeature.set(fD)
        #
        rso = self.__fetchReferenceSequence(dbName, dbAccession, dbIsoform)
        if self.__seqFetchError:
            return None
        #
        if (refSeqBeg is not None) and (refSeqBeg < 1):
            if self.__seqFetchError:
                self.__seqFetchError += "\n"
            #
            self.__seqFetchError += "Invalid SEQ BEGIN number: " + str(refSeqBeg) + "."
        #
        if refSeqEnd and (refSeqEnd > rso.getSequenceLength()):
            if self.__seqFetchError:
                self.__seqFetchError += "\n"
            #
            accCode = dbAccession
            if dbIsoform:
                accCode = dbIsoform
            #
            self.__seqFetchError += "Invalid SEQ END number: " + str(refSeqEnd) + " > " + str(rso.getSequenceLength()) + " ( sequence length of " + accCode + " )."
        #
        if self.__seqFetchError:
            return None
        #
        packedSeqLabel = self.__loadReferenceSequence(entityId=entityId, partId=partId, refSeqObj=rso, refSeqBeg=refSeqBeg, refSeqEnd=refSeqEnd)
        if self.__debug:
            self.__sds.dump(self.__lfh)
        self.__sds.serialize()
        #
        # update alignment stats -
        #
        alstat = AlignmentStatistics(reqObj=self.__reqObj, maxRefAlign=self.__maxRefAlign, verbose=self.__verbose, log=self.__lfh)
        alstat.doUpdate()

        if (self.__verbose):
            self.__lfh.write("+UpdatePolymerEntityReference.addReferenceSequences() completed returning sequence label %r\n" % packedSeqLabel)
        return packedSeqLabel

    def getSeqFetchError(self):
        return self.__seqFetchError

    def __fetchReferenceSequence(self, dbName, dbAccession, dbIsoform):
        """  Return a reference sequence object loading the content dictionary for the input reference sequence  -

        """
        rsu = ReferenceSequenceUtils(siteId=self.__siteId, verbose=self.__verbose, log=self.__lfh)
        #
        dbResource = self.__srd.convertDbNameToResource(dbName)
        if self.__verbose:
            self.__lfh.write("\n+UpdatePolymerEntityReference.fetchReferenceSequence() fetch reference sequence: dbName %s accession %s isoform %s\n" %
                             (dbName, dbAccession, dbIsoform))

        if dbResource in ['UNP']:
            idCode = dbAccession
            if dbIsoform is not None and len(dbIsoform) > 0:
                idCode = dbIsoform
            if self.__verbose:
                self.__lfh.write("\n+UpdatePolymerEntityReference.fetchReferenceSequence() UNP search for idCode %s \n" % idCode)
            dt = rsu.fetchUniProt([idCode])
            if self.__verbose:
                self.__lfh.write("\n+UpdatePolymerEntityReference.fetchReferenceSequence() after UNP search for idCode %s return %r\n" % (idCode, dt.items()))
            if idCode in dt:
                d = dt[idCode]
            elif dbAccession in dt:
                d = dt[dbAccession]
            else:
                d = {}

            if d and ('db_name' not in d):
                # guess --
                if dbAccession[0] in ['P', 'Q', 'O']:
                    d['db_name'] = 'SP'
                else:
                    d['db_name'] = 'TR'
            # JDW JDW
            # if dbIsoform is not None:
            #    d['db_isoform']=dbIsoform
            # else:
            #    d['db_isoform']=''
        elif dbResource in ['GB', 'DBJ', 'EMB']:
            xmlPath = None
            if self.__debug:
                xmlPath = os.path.join(self.__sessionPath, dbAccession + '.xml')
            d = rsu.fetchNcbiGi(dbAccession, xmlPath=xmlPath)
            d['db_accession'] = dbAccession
            d['db_name'] = dbName
        else:
            d = {}

        if not d:
            accCode = dbAccession
            if dbIsoform:
                accCode = dbIsoform
            #
            self.__seqFetchError = "Fetch reference sequence [ dbName=" + dbName + ", Accession=" + accCode + "] failed."
        #

        if self.__verbose:
            self.__lfh.write(
                "+UpdatePolymerEntityReference.fetchReferenceSequence() content dictionary for  dbResrource %s accession %s isoform %s\n" %
                (dbResource, dbAccession, dbIsoform))
            self.__lfh.write("+UpdatePolymerEntityReference.fetchReferenceSequence() content dictionary keys  %s\n" % d.keys())
            for k, v in d.items():
                self.__lfh.write(" + Key %s value %r\n" % (k, v[:512]))

        rso = ReferenceSequence(verbose=self.__verbose, log=self.__lfh)
        rso.clear()
        rso.set(d)
        # rso.dump()
        #
        return rso

    def __assignPart(self, partD, seqBegin, seqEnd):
        """  Assign the entity partId to the
        """
        pMatch = -1
        for partId, pTup in partD.items():
            if ((seqBegin == pTup[1]) and (seqEnd == pTup[2])):
                pMatch = partId

        return pMatch

    def __getAuthRefAssignments(self, entityId):
        """  Get author provided reference assignments.   Return a lists of reference database accessions.

             Note that this is currently being do a the level of
             entity as the entity parts are not explicit in input data at this time. 2013-11-14
        """
        depSeqAssignD = self.__sds.getDepositorReferenceAssignments()
        sADep = SequenceAssignDepositor(verbose=self.__verbose, log=self.__lfh)
        sADep.set(depSeqAssignD)
        refSeqAssignL = sADep.getReferenceList(entityId)

        if (self.__debug):
            self.__lfh.write("+UpdatePolymerEntityReference.__getAuthRefAssignments() for entityId %r reference count is %d\n" % (entityId, sADep.getReferenceCount(entityId)))
            sADep.printIt(log=self.__lfh)
            for ii, rsa in enumerate(refSeqAssignL):
                self.__lfh.write("+UpdatePolymerEntityReference.__getAuthRefAssignments() depositor reference  %d\n" % (ii + 1))
                rsa.printIt(self.__lfh)
        #
        authDbNameList = []
        authDbAccessionList = []
        authDbCodeList = []
        for rsa in refSeqAssignL:
            dbName, dbCode, dbAccession = rsa.getDbReference()
            if dbName not in ['.', '?']:
                authDbNameList.append(dbName)
            else:
                authDbNameList.append('')
            if dbAccession not in ['.', '?']:
                authDbAccessionList.append(dbAccession)
            else:
                authDbAccessionList.append('')
            if dbCode not in ['.', '?']:
                authDbCodeList.append(dbCode)
            else:
                authDbCodeList.append('')
        #
        if (self.__verbose):
            self.__lfh.write(
                "+UpdatePolymerEntityReference.__getAuthRefAssignments()  dbNames %r dbAccessions %r dbCodes %r\n" %
                (authDbNameList, authDbAccessionList, authDbCodeList))

        return authDbNameList, authDbAccessionList, authDbCodeList

    def __getPartDetails(self, entityId, numExtra=0):
        """ Return dictionary of part boundaries and types ---   pD[pId]=(pId,pSeqBegin,pSeqEnd,pType,taxId)
        """

        pD = {}
        seqIds = self.__sds.getGroup(groupId=entityId)
        # if len(seqIds)==0:
        #    return pD

        seqFeature = SequenceFeature()
        # JDW CHANGE
        # seqId0=seqIds[0]
        seqId0 = entityId

        polymerTypeCode = 'AA'
        #
        partIdList = self.__sds.getPartIds(seqId0, dataType="sequence", seqType="auth")

        for partId in partIdList:
            tD = {}
            vL = self.__sds.getVersionIds(seqId0, partId=partId, altId=1, dataType="sequence", seqType="auth")
            if len(vL) > 0:
                pfD = self.__sds.getFeature(seqId0, seqType="auth", partId=partId, altId=1, version=vL[0])
                seqFeature.set(pfD)
                pId, pSeqBegin, pSeqEnd, pType = seqFeature.getAuthPartDetails()
                polymerTypeCode = seqFeature.getPolymerType()
                taxId = seqFeature.getSourceTaxId()
                pD[pId] = (pId, pSeqBegin, pSeqEnd, pType, taxId)
                lastPart = pId

        if numExtra > 0:
            pv = self.__placeHolderValue
            pId = lastPart
            for ii in range(0, numExtra):
                pId += 1
                pD[pId] = (pId, pv, pv, pv, pv)
        #
        seqAuthIdx = self.__sds.getSequence(seqId=seqId0, seqType='auth', partId=1, altId=1, version=vL[0])
        r3List = []
        for sTup in seqAuthIdx:
            (r3, sIdx, comment, idx) = sTup
            r3List.append(r3)
        self.__srd.cnvList3to1WithMods(r3List)
        #
        return pD

    def __getAuthFeatures(self, entityId, partId):
        """  Get feature data for the author sequence entity sequence.

             Returns a dictionary of features for the more recent sequence version.
        """
        featureD = {}
        #
        if (self.__verbose):
            self.__lfh.write("+UpdatePolymerEntityReference.__getAuthFeatures() entityId %s partId %s\n" % (entityId, partId))
        seqFeature = SequenceFeature()

        # get the author sequence identifier
        seqIds = self.__sds.getGroup(groupId=entityId)

        # JDW CHANGE
        # if len(seqIds)==0:
        #    return featureD

        #
        # Collect the feature dictionaries for the parts of each sequence version.
        # seqId0=seqIds[0]
        seqId0 = entityId
        #
        altId = 1
        verList = self.__sds.getVersionIds(seqId0, partId=partId, altId=altId, dataType="sequence", seqType='auth')
        if len(verList) == 0:
            return featureD

        featureD = self.__sds.getFeature(seqId=seqId0, seqType='auth', partId=partId, altId=altId, version=verList[0])

        if (self.__debug):
            kys = featureD.keys()
            for k in sorted(kys):
                v = featureD[k]
                self.__lfh.write("+UpdatePolymerEntityReference.__getAuthFeatures() %-35s : %s\n" % (k, v))

        return featureD

    def __loadReferenceSequence(self, entityId, partId, refSeqObj, refSeqBeg=None, refSeqEnd=None):
        """ Do what is needed to load the input reference sequence into the sequence data store  -

            Return: the packed sequence label identifier for the loaded sequence or None

        """
        #
        if (self.__verbose):
            self.__lfh.write(
                "\n+UpdatePolymerEntityReference.__loadReferenceSequence() entityId %s partId %r refSeqBeg %r refSeqEnd %r \n" %
                (entityId, partId, refSeqBeg, refSeqEnd))
        #
        #
        #
        seqFeature = SequenceFeature(self.__verbose)
        seqIdList = self.__sds.getGroup(groupId=entityId)

        seqId0 = entityId
        verList = self.__sds.getVersionIds(seqId0, partId=partId, altId=1, dataType="sequence", seqType='auth')
        if (len(verList) < 1):
            if (self.__verbose):
                self.__lfh.write("+UpdatePolymerEntityReference.__loadReferenceSequence() entity %s partId %d seqId %s has no sequence data\n" % (entityId, partId, seqId0))
            return None

        verLatest = verList[0]
        seqFeature = self.__sds.getFeatureObj(seqId=seqId0, seqType='auth', partId=partId, altId=1, version=verLatest)
        (pId, seqBeg, seqEnd, seqPartType) = seqFeature.getAuthPartDetails()
        polyTypeCode = seqFeature.getPolymerType()

        if (self.__verbose):
            self.__lfh.write(
                "+UpdatePolymerEntityReference.__loadReferenceSequence() read auth sequence for entity %r partId %d seqBegin %d seqEnd %d pId %d type %s\n" %
                (entityId, partId, seqBeg, seqEnd, pId, polyTypeCode))

        #
        altIdList = self.__sds.getAlternativeIds(seqId0, dataType="sequence", seqType="ref", partId=partId)
        lenAltIdList = len(altIdList)
        #
        if lenAltIdList < 1:
            nextAltId = 1
        else:
            nextAltId = int(altIdList[0]) + 1

        sL = SequenceLabel()
        sL.set(seqType='ref', seqInstId=seqId0, seqPartId=partId, seqAltId=nextAltId, seqVersion=1)
        nextRefLabel = sL.pack()

        if (self.__verbose):
            self.__lfh.write("+UpdatePolymerEntityReference.__loadReferenceSequence() create new ref sequence label for entity %r partId %r reference sequence count %d next altId %d label %s\n" %
                             (entityId, partId, len(altIdList), nextAltId, nextRefLabel))

        #
        # Do we know the boundaries for the new reference sequence?
        #
        if refSeqBeg is None and refSeqEnd is None:

            if self.__verbose:
                self.__lfh.write("+UpdatePolymerEntityReference.__loadReferenceSequence() NO valid boundaries detected -  calculating these from alignment\n")

            pA = PairwiseAlign()
            pA.setVerbose(self.__verbose)

            # need to do alignment to select the useful range of the sequence.
            # Set the alignment reference as the 'auth' sequence with the current part id
            seqAuthIdx = self.__sds.getSequence(seqId=seqId0, seqType='auth', partId=partId, altId=1, version=verLatest)
            r3L = []
            for tup in seqAuthIdx[seqBeg - 1:seqEnd]:
                r3L.append(tup[0])
            pA.setReferenceSequence(r3L, 'auth:' + seqId0 + '_P' + str(partId))
            authSeqLen = len(seqAuthIdx[seqBeg - 1:seqEnd])

            testRefIdx = refSeqObj.getSequenceWithIndex(polyTypeCode=polyTypeCode, seqBegin=refSeqBeg, seqEnd=refSeqEnd)
            r3L = []
            for tup in testRefIdx:
                r3L.append(tup[0])
            pA.addTestSequence(r3L, 'ref:' + str(nextAltId) + '_P' + str(partId))
            pA.doAlign()
            aL = pA.getAlignment('ref:' + str(nextAltId) + '_P' + str(partId))
            alignLength = len(aL)
            #
            # where both are not gap
            #
            idxL = []
            iRef = 1
            for ii, aTup in enumerate(aL, start=1):
                if (self.__verbose):
                    self.__lfh.write("+UpdatePolymerEntityReference.__loadReferenceSequence() alIdx %d  refIdx %d ref %r auth %r\n" % (ii, iRef, aTup[0], aTup[1]))
                if aTup[0] != self.__gapSymbol:
                    idxL.append(iRef)
                if (aTup[1] != self.__gapSymbol):
                    iRef += 1

            refSeqBeg = min(idxL)
            refSeqEnd = max(idxL)
            if self.__verbose:
                self.__lfh.write("+UpdatePolymerEntityReference.__loadReferenceSequence() partId %d type %s using aligned refSeqBeg %d refSeqEnd %d\n" %
                                 (partId, polyTypeCode, refSeqBeg, refSeqEnd))
            numMatch = 0
            numMatchGaps = 0
            for aTup in aL:
                if aTup[0] == aTup[1]:
                    numMatch += 1
                if aTup[0] == aTup[1] or aTup[1] == self.__gapSymbol:
                    numMatchGaps += 1
            if (self.__debug):
                for ii, aTup in enumerate(aL, start=1):
                    self.__lfh.write(" -----   %d   %r\n" % (ii, aTup))

        elif refSeqBeg is not None and refSeqEnd is None:
            refSeqEnd = refSeqObj.getSequenceLength()
        # JDW
        elif refSeqBeg is None and refSeqEnd is not None:
            refSeqBeg = 1

        # ---------- Using either supplied or computed ranges load the new reference sequence ----------
        # Interpret the range information relative to the full length sequence stored in refSeqObj()
        #
        seqRefIdx = refSeqObj.getSequenceWithIndex(polyTypeCode=polyTypeCode, seqBegin=refSeqBeg, seqEnd=refSeqEnd)
        if self.__verbose:
            self.__lfh.write(
                "+UpdatePolymerEntityReference.__loadReferenceSequence() entity %r partId %r type %s refSeqBeg %r refSeqEnd %r\n" %
                (entityId, partId, polyTypeCode, refSeqBeg, refSeqEnd))
        if self.__debug:
            self.__lfh.write("+UpdatePolymerEntityReference.__loadReferenceSequence() reference sequence with index - seqRefIdx %r\n" % seqRefIdx)

        self.__sds.setSequence(seqRefIdx, seqId0, 'ref', partId=partId, altId=nextAltId, version=1)
        #
        # Load feature data
        #
        seqFeature.clear()
        #
        dbName, dbCode, dbAccession, dbIsoform = refSeqObj.getDatabaseInfo()
        if (self.__verbose):
            self.__lfh.write(
                "+UpdatePolymerEntityReference.__loadReferenceSequence() loading reference sequence with dbName %r  dbCode %r dbAccession %r\n" %
                (dbName, dbCode, dbAccession))
        #
        seqFeature.setId(dbName=dbName, dbCode=dbCode, dbAccession=dbAccession, dbIsoform=dbIsoform)
        #
        taxId = refSeqObj.getTaxId()
        orgName = refSeqObj.getSourceName()
        strain = refSeqObj.getSourceStrain()
        sourceCommonName = refSeqObj.getSourceCommonName()
        seqFeature.setSource(organism=orgName, strain=strain, taxid=taxId, commonName=sourceCommonName)
        #
        proteinName = refSeqObj.getName()
        synonyms = refSeqObj.getSynonyms()
        geneName = refSeqObj.getGeneName()
        seqFeature.setRefSeqNames(proteinName=proteinName, synonyms=synonyms, geneName=geneName)
        ec = refSeqObj.getEnzymeClass()
        keywords = refSeqObj.getKeywords()
        comments = refSeqObj.getComments()
        seqFeature.setRefSeqDetails(enzymeClass=ec, description='', comments=comments, keywords=keywords)
        #
        seqFeature.setAuthRefAlignRange(refMatchBegin=refSeqBeg, refMatchEnd=refSeqEnd)
        seqFeature.setPolymerType(polyTypeCode)
        #
        dbIsoformDescription = refSeqObj.getDbIsoformDescription()
        seqFeature.setDbIsoformDescription(description=dbIsoformDescription)
        #
        # jdw add missing attributes --

        # less useful --
        seqFeature.setItem('ORG_ORDER_ID', 0)
        sOrder = lenAltIdList + 1
        seqFeature.setRefSortOrder(sortIndex=sOrder, sortMetric=self.__defaultInsertSortMetric)
        #
        seqFeature.setAuthPartDetails(partId, seqBeg, seqEnd, seqPartType)
        self.__sds.setFeature(seqFeature.get(), seqId0, 'ref', partId=partId, altId=nextAltId, version=1)
        #
        return nextRefLabel


if __name__ == '__main__':
    pass
