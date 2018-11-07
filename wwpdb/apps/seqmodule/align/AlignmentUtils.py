##
# File:    AlignmentUtils.py
# Date:    14-April-2013
#
# Updates:
# 15-Apr-2013  jdw  refactored from AlignmentView.py
# 16-Apr-2013  jdw  fix dual gap index reference error
# 16-Apr-2013  jdw  adjust terminal index conditions for conflict annotation
# 17-Apr-2013  jdw  add save and read methods for the current/revised alignment id list
# 17-Apr-2013  jdw  again adjusting last sequence index for conflict assignment -
# 23-Apr-2013  jdw  add dictionary detection of residue modifications.
# 21-May-2013  jdw  rewrite  __annotateAlignmentWithDefaultConflicts()
# 12-Jun-2013  jdw  Set default annotation for any internal conflicts lacking annotation to "conflict"
# 30-Nov-2013  jdw  filter gaps prior to alignment by default.
# 12-Dec-2013  jdw  change MET css to avoid MET/gap highlighting
#  6-Jan-2014  jdw   'initiating methionine'
# 02-Feb-2014  jdw  avoid default updates of existing annotations
# 09-Feb-2014  jdw  handle prior 'conflicts' in assignment of default annotations.
# 20-Mar-2014  jdw  Gap handling revised - now preserving any gaps in input sequences in the alignment prep step
# 21-Mar-2014  jdw  in sequence prep filter leading gaps from reference sequence if sequence is type auth
#                   also ignore trailing gaps in test sequnences --
# 14-May-2014  jdw  simplify the alignment storage io  add AlignDataStore()
# 22-May-2014  jdw  completed overhaul of alignment data flow -
#  5-Jun-2014  jdw  renumber auth sequence in the alignment after any edits to this sequence -
#  3-Jul-2014  jdw  adjust the order of alignment to auth, ref, xyz  -- ref before xyz always --
#  4-Jul-2014  jdw  make sequence alignment annotations are preserved in individual sequences -
# 10-Sep-2014  jdw  add method __updateAuthPartDetails() to revise author sequence part boundaries after alignment.
#                   part boundary revision tested for initial alignments __doAlignment().
# 16-Oct-2014  jdw fixed update of part boundaries, alignment length consistency following edits and filtering
#                  empty alignment records following deletions.
# 22-Oct-2014  jdw Add explicit mapping of author sequence part boudaries to the alignment coordinate system.
# 28-Oct-2014  jdw Adjust clear prior annotatoins -
# 30-Nov-2014  jdw Fix unitialized sequence length in gap filtering
# 26-Aug-2015  jdw Fix sequence bounds to correspond to alignment on call to __annotateAlignmentWithDefaultConflicts()
# 30-Aug-2017  zf  Add fragment range for ref sequence using addTestSequenceWithLinkAndRange() function
#                  add __getAuthDefinedMutations() function to get author defined mutation information from _entity.pdbx_mutation
#
##
##
"""
Shared alignment methods -

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.09"

import os
import sys
import copy
import itertools
import traceback
from operator import itemgetter

from wwpdb.apps.seqmodule.io.SequenceEditStore import getEditStoreFilename,SequenceEditStore
from wwpdb.apps.seqmodule.io.SequenceDataStore import SequenceDataStore
from wwpdb.apps.seqmodule.io.AlignDataStore import AlignDataStore
from wwpdb.apps.seqmodule.util.SequenceLabel import SequenceLabel, ResidueLabel, SequenceFeature
from wwpdb.apps.seqmodule.util.SequenceReferenceData import SequenceReferenceData
from wwpdb.apps.seqmodule.io.ChemCompUtils import ChemCompUtils
from wwpdb.utils.pair_align.wrapper.libPairwiseAlignPackage import PairwiseAlign
from wwpdb.utils.rcsb.PathInfo import PathInfo

def IsCompatible(comment, newResName):
    compatibilityMap = { "ala-gly-like" : [ "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", \
                                            "MET", "MSE", "PHE", "PRO", "PYL", "SEC", "SER", "THR", "TRP", "TYR", "UNK", "VAL" ], \
                         "c-gamma-like" : [ "ARG", "ASN", "ASP", "GLN", "GLU", "HIS", "ILE", "LEU", "LYS", "MET", "MSE", "PHE", \
                                            "PYL", "TRP", "THR", "TYR", "VAL" ], \
                         "c-delta-like" : [ "ARG", "GLN", "GLU", "HIS", "ILE", "LYS", "PHE", "PYL", "TRP", "TYR" ], \
                         "ile-val-like" : [ "ILE", "VAL" ], \
                         "ser-thr-like" : [ "SER", "THR" ], \
                         "phe-tyr-like" : [ "PHE", "TYR" ], \
                         "lys-pyl-like" : [ "LYS", "PYL" ] }

    #
    if not comment:
        return False
    #
    for patternType,allowList in compatibilityMap.items():
        if (comment.find(patternType) != -1) and (newResName in allowList):
             return True
        #
    #
    return False

class AlignmentUtils(object):

    """ Sequence alignment utility methods.

    """

    def __init__(self, reqObj=None, verbose=False, log=sys.stderr):
        self.__verbose = verbose
        self.__debug = True
        self.__lfh = log
        self.__reqObj = reqObj

        self.__operation = None
        self.__sessionObj = None
        self.__alignGroupId = None
        self.__excludedPartIdList = []
        ##
        self.__srd = SequenceReferenceData(verbose=self.__verbose, log=self.__lfh)
        self.__gapSymbol = self.__srd.getGapSymbol()
        self.__sds = None
        #
        self.__setup()

    def __setup(self):
        try:
            self.__sessionObj = self.__reqObj.getSessionObj()
            self.__sessionPath = self.__sessionObj.getPath()
            self.__operation = self.__reqObj.getAlignmentOp()
            self.__siteId = self.__reqObj.getValue('WWPDB_SITE_ID')
            self.__identifier = self.__reqObj.getValue("identifier")
            self.__ccU = ChemCompUtils(siteId=self.__siteId, verbose=self.__verbose, log=self.__lfh)
            # --  -----------------------------------------------------------------------------
            # Part boundaries in author sequence coordinate system -
            self.__authPartD = {}
            # Part boundaries mapped to the aligned author coordinate system -
            self.__authPartAlignD = {}
            # List of part id's order in increasing begining sequence boundary
            self.__authPartIdList = []
            # linkages between sequence parts in author sequence coordinate system --
            self.__authLinkerList = []
            # Mininum and maximum author sequence indices - author sequence coordinate system
            self.__authSeqBounds = None
            # Mininum and maximum author sequence indices - mapped to author aligned coordinate system
            self.__authSeqAlignBounds = None
            #
            self.__pI = PathInfo(siteId=self.__siteId, sessionPath=self.__sessionPath, verbose=self.__verbose, log=self.__lfh)
            if (self.__verbose):
                self.__lfh.write("+AlignmentUtils.__setup() sessionId %s\n" % self.__sessionObj.getId())
                self.__lfh.write("+AlignmentUtils.__setup() operation  %s\n" % self.__operation)
        except:
            if (self.__verbose):
                self.__lfh.write("+AlignmentUtils.__setup() failed\n")
                traceback.print_exc(file=self.__lfh)

    def doAlignment(self, inpAlignIdList=None):
        """ Top-level
        """
        self.__alignGroupId, refDef, seqList = self.__prepareSequenceData(inpAlignIdList=inpAlignIdList)

        return self.__doAlignment(refDef=refDef, seqList=seqList)

    def getAlignmentGroupId(self):
        return self.__alignGroupId

    def getSelectedAnnotations(self, inpAlignIdList=[]):
        """ Get selected annotations for display in alignment view display --
        """
        aD = {}
        for ky in ['mutation', 'mutationOrig', 'description', 'descriptionOrig']:
            aD[ky] = ''

        if self.__sds is None:
            self.__sds = SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            #
        authId = None
        for seqId in inpAlignIdList:
            if seqId.startswith('auth'):
                authId = seqId
                break
        if authId is None:
            return aD

        authLabel = SequenceLabel()
        authLabel.unpack(authId)
        (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = authLabel.get()
        authFObj = self.__sds.getFeatureObj(seqId=seqInstId, seqType=seqType, partId=seqPartId, altId=seqAltId, version=seqVersion)

        aD['mutation'] = authFObj.getEntityMutationDetails()
        aD['mutationOrig'] = authFObj.getEntityMutationDetailsOrig()
        aD['description'] = authFObj.getEntityDescription()
        aD['descriptionOrig'] = authFObj.getEntityDescriptionOrig()
        aD['source_method'] = authFObj.getEntitySourceMethod()
        return aD

    def getEntryDetails(self, kyList=['STRUCT_TITLE', 'CITATION_TITLE', 'PDB_ID']):
        eD = {}
        if self.__sds is None:
            self.__sds = SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        for ky in kyList:
            eD[ky] = self.__sds.getEntryDetail(detailKey=ky)
        #
        return eD

    def __prepareSequenceData(self, inpAlignIdList=None):
        """
        From the list of input sequences in inpAlignIdList establish the alignment reference and
        the associated test sequences for this alignment.

        Populate the following internal copies of details for the reference and test sequences:

        Returns --

        For the reference sequence -

        refSeq   = [3-letter-code, ...]
        refSeqId = (refId,refLabel)
        refSeqIdx =[(compId,str(indx),'comment',index in seq),(),...]
        refSeqFD={}   feature dictionary

        returns this as refDef (refId refLabelObj, ref3L, refSeqIdx, refSeqFD)

        For the test sequences in the alignment are stored as -

        seqList= [[seqId,seqLabelObj,seq3L,seqLWithIndex,seqFD]]

        seqId = sequence identifier (ie. auth_A_1_1_1)
        seqlabelObj = From class SequenceLabel()
        seq3L = [compId,compId,compId,...]
        seqWithIndex = sequence as stored in sequence data store below -
        seqFD  = sequence feature dictionary as from class SequenceFeature()

        Sequence data store storage model - (seqWithIndex)

        'auth' [(comp_id,  str(index in seq),         'comment'  index in sequence), (), ...]
        'ref'  [(comp_id,  str(offset in full dbseq), 'comment', index in sequence), (), ...]
        'xyz'  [(comp_id,  str(auth_seq_id),          'comment', index in sequence ), (), ...]

        """
        #
        if (self.__verbose):
            self.__lfh.write("\n\n+AlignmentUtils.__prepareSequenceData() starting sessionId %s\n" % self.__sessionObj.getId())
            self.__lfh.write("+AlignmentUtils.__prepareSequenceData() with operation:  %s\n" % self.__operation)
        #
        # Get alignment reference and ordered list of sequences to align.
        #
        refId, alignIdList = self.__getAlignmentReference(alignIdList=inpAlignIdList)
        if (self.__verbose):
            self.__lfh.write("+AlignmentUtils.__prepareSequenceData() reference sequence Id (refId) %s\n" % refId)
        if (self.__debug):
            self.__lfh.write("+AlignmentUtils.__prepareSequenceData() test align list length   %d\n" % len(alignIdList))
            for id in alignIdList:
                self.__lfh.write("+AlignmentUtils.__prepareSequenceData() test align Id %s\n" % id)
        #
        #

        #
        if self.__sds is None:
            self.__sds = SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        if (self.__debug):
            self.__sds.dump(self.__lfh)

        refLabel = SequenceLabel()
        refLabel.unpack(refId)
        (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = refLabel.get()
        #
        # get the group id for the reference sequence -
        # JDW
        # alignGroupId=self.__sds.getGroupId(seqInstId)
        alignGroupId = seqInstId
        refSeqIdx = self.__sds.getSequence(seqId=seqInstId, seqType=seqType, partId=seqPartId, altId=seqAltId, version=seqVersion)
        refSeqFD = self.__sds.getFeature(seqId=seqInstId, seqType=seqType, partId=seqPartId, altId=seqAltId, version=seqVersion)

        if (self.__verbose):
            self.__lfh.write(
                "+AlignmentUtils.__prepareSequenceData() Reference id %s type %s group id %s sequence length %d\n" %
                (seqInstId, seqType, alignGroupId, len(refSeqIdx)))
        #
        refSeq = []
        tRefSeqIdx = []
        leadingFlag = True
        for sPos in refSeqIdx:
            # JDW filter only leading gaps from autn ref sequence
            # self.__lfh.write("+AlignmentUtils.__prepareSequenceData() leadingFlag %r seqType %s residue %s\n" %  (leadingFlag,seqType,sPos[0]))
            if (leadingFlag and (seqType == 'auth') and (sPos[0] == self.__gapSymbol)):
                #self.__lfh.write("+AlignmentUtils.__prepareSequenceData() SKIPPED leadingFlag %r seqType %s residue %s\n" %  (leadingFlag,seqType,sPos[0]))
                continue
            else:
                leadingFlag = False
            if ((seqType != 'auth') and (sPos[0] == self.__gapSymbol)):
                continue
            refSeq.append(sPos[0])
            tRefSeqIdx.append(sPos)

        if (self.__verbose):
            self.__lfh.write("+AlignmentUtils.__prepareSequenceData() reference sequence Id (refId) %s\n" % refId)

        if (self.__debug):
            self.__lfh.write("+AlignmentUtils.__prepareSequenceData() align id list length   %d\n" % len(alignIdList))
            for id in alignIdList:
                self.__lfh.write("      + Alignment test sequence Id %s\n" % id)
        #
        preserveGapFlag = True
        seqList = []
        for id in alignIdList:
            aLabel = SequenceLabel()
            aLabel.unpack(id)
            (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = aLabel.get()
            seqIdx = self.__sds.getSequence(seqId=seqInstId, seqType=seqType, partId=seqPartId, altId=seqAltId, version=seqVersion)
            seqFD = self.__sds.getFeature(seqId=seqInstId, seqType=seqType, partId=seqPartId, altId=seqAltId, version=seqVersion)
            seq3L = []
            tSeqIdx = []
            typicalLink = []
            leadingFlag = True
            #
            # -- ignore trailing gaps in test sequences ---
            iLast = len(seqIdx) - 1
            for i in range(len(seqIdx) - 1, -1, -1):
                sPos = seqIdx[i]
                if (sPos[0] != self.__gapSymbol):
                    iLast = i
                    break
            #
            for sPos in seqIdx[:iLast + 1]:
                # JDW Gap handling revised  2014-Mar-20
                # Now preserving any gaps in input sequences in the alignment prep step !
                if sPos[0] != self.__gapSymbol:
                    seq3L.append(sPos[0])
                    tSeqIdx.append(sPos)
                    typicalLink.append(0 if ('long_begin' in sPos[2]) else 1)
                elif (preserveGapFlag):
                    # still discard leading gaps -
                    if (leadingFlag and (sPos[0] == self.__gapSymbol)):
                        continue
                    else:
                        seq3L.append(sPos[0])
                        tSeqIdx.append(sPos)
                        typicalLink.append(1)
                        leadingFlag = False

            seqList.append([id, aLabel, seq3L, tSeqIdx, seqFD, typicalLink])
            if (self.__debug):
                if preserveGapFlag:
                    self.__lfh.write("+AlignmentUtils.__prepareSequenceData() Saved id %s seqType %s  poly type %s partId %r ungapped length %d gapped length %d\n" %
                                     (id, seqType, seqFD['POLYMER_TYPE'], seqPartId, len(seq3L), len(seqIdx)))
                else:
                    self.__lfh.write("+AlignmentUtils.__prepareSequenceData() PRESERVE GAPS %r Saved id %s seqType %s  poly type %s partId %r gapped length %d\n" %
                                     (preserveGapFlag, id, seqType, seqFD['POLYMER_TYPE'], seqPartId, len(seq3L)))

        #
        if (self.__verbose):
            self.__lfh.write("+AlignmentUtils.__prepareSequenceData() align list length  %d\n" % len(seqList))
        #
        refDef = (refId, refLabel, refSeq, tRefSeqIdx, refSeqFD)

        # Get any partitioning details for the reference entity
        self.__buildAuthPartDetails(refId)

        if self.__debug:
            self.__dumpAlignSeqInput(io=self.__lfh, refDef=refDef, seqList=seqList)

        if (self.__verbose):
            self.__lfh.write("+AlignmentUtils.__prepareSequenceData() completed.\n")
        #
        return (alignGroupId, refDef, seqList)
     
    def __doAlignment(self, refDef, seqList):
        """ Align selected test sequences with the reference sequence.  Produce a list
            of aligned sequences annotated with index information from the input  sequences.

            Input data is prepared by __getSequenceData() -

            refDef (refId refLabelObj, ref3L, refSeqIdx, refSeqFD)

            seqList= [[seqId,seqLabelObj,seq3L,seqLWithIndex,seqFD,typicalLink]]

            On returns - Output Data storage model for aligned sequence list -

            alignSeqList [[Sequence identifier (ie. auth_A_1_1_1 as used in sequence data store),
                           SequenceLabel() object,
                           Aligned sequence with index details,
                           Conflict flag list (bool),
                           Feature dictionary for input sequence],,,]

        where aligned sequences with index is a list of tuples of the form -

        (one-letter-code, 3-letter-code, original residue index, position in sequence, position in alignment, comment)

        """
        #
        numParts = len(self.__authPartIdList)
        #
        seqFeature = SequenceFeature()
        #
        refSeqId = refDef[0]
        refSeqLabel = refDef[1]
        refSeq3 = refDef[2]
        refSeqIdx = refDef[3]
        refSeqFD = refDef[4]
        refSeqType = refSeqLabel.getSequenceType()
        seqFeature.set(refSeqFD)
        self.__lfh.write("\n\n+AlignmentUtils.doAlignment() Starting new alignment with reference sequence refId %s with %d test sequences\n" % (refSeqId, len(seqList)))
        #
        authDefinedMutations = self.__getAuthDefinedMutations(seqFeature.getEntityMutationDetailsOrig())
        #
        pA = PairwiseAlign()
        pA.setVerbose(self.__verbose)
        #   set reference sequence -> ( Sequence(e.g.  [compId,compId,...]), sequence ID(e.g. ref_A_1_1_1) )
        pA.setReferenceSequence(refSeq3, refSeqId)
        for aSeq in seqList:
            authLabel = SequenceLabel()
            authLabel.unpack(aSeq[0])
            (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = authLabel.get()
            #   set test sequences -> ( Sequence(e.g.  [compId,compId,...]), sequence ID(e.g. auth_A_1_1_1) )
            if (refSeqType == 'auth') and (seqType == 'ref'):
                #pA.addTestSequenceWithLinkAndRange(aSeq[2], aSeq[0], aSeq[5], self.__authPartD[seqPartId][1], self.__authPartD[seqPartId][2])
                pA.addTestSequenceWithLinkAndRange(aSeq[2], aSeq[0], [], self.__authPartD[seqPartId][1], self.__authPartD[seqPartId][2])
            else:
                pA.addTestSequenceWithLink(aSeq[2], aSeq[0], aSeq[5])
            #
        #
        #self.__writeAlignSeqInfo(refSeqId, refSeq3, refSeqType, seqList)
        #
        #myFails = pA.doAlignConsensus()
        myFails = pA.doMultipleAlign()
        if (self.__verbose):
            if myFails > 0:
                self.__lfh.write("+AlignmentUtils.__doAlignment() Consensus fail count = %d\n" % len(myFails))
        #
        # Post process the aligned sequences -
        #
        # Get the reference sequence the leading sequence in the alignment -
        #
        aSeq0 = seqList[0]
        aL0 = pA.getAlignment(aSeq0[0])

        alignRefSeq = []
        for aTup in aL0:
            alignRefSeq.append(aTup[0])

        alignRefSeqIdx = self.__annotateAlignmentWithIndex(alignRefSeq, refSeqIdx, clearSelectedComments=True)
        lenRefAlignment = len(alignRefSeqIdx)

        if (self.__verbose):
            self.__lfh.write("+AlignmentUtils.__doAlignment() First     alignment length for reference %s is %d\n"
                             % (refSeqId, len(alignRefSeq)))
            self.__lfh.write("+AlignmentUtils.__doAlignment() Annotated alignment length for reference %s is %d\n"
                             % (refSeqId, len(alignRefSeqIdx)))

        # ----------------
        # Create and return alignment data structure -
        #
        alignSeqList = []
        #
        # Set the default conflict state to false
        conflictRefL = [(0, '')] * len(aL0)

        # Assign annotations to regions outside of sequence parts --
        if (refSeqType == 'auth'):
            self.__mapAuthPartDetails(alignRefSeqIdx)
            if (self.__verbose):
                self.__lfh.write("+AlignmentUtils.__doAlignment() adding out-of-part annotations to the reference sequence %s\n" % refSeqId)
            self.__annotateAuthOutOfPartComment(alignRefSeqIdx, conflictRefL)

        #  Reference sequence is always first on the alignSeqList=[]
        #
        alignSeqList.append([refSeqId, refSeqLabel, alignRefSeqIdx, conflictRefL, refSeqFD])

        #  For each input test sequence in the alignment -
        for aSeq in seqList:
            testSeqId = aSeq[0]
            testSeqLab = aSeq[1]
            testSeqIdx = aSeq[3]
            testSeqFD = aSeq[4]
            testSeqType = testSeqLab.getSequenceType()
            #
            aL = pA.getAlignment(testSeqId)

            if (self.__verbose):
                self.__lfh.write("\n\n+AlignmentUtils.doAlignment() for reference id %s test id %s alignment length %d\n" % (refSeqId, testSeqId, len(aL)))
            #
            if testSeqType in ['auth', 'ref']:
                seqFeature.set(testSeqFD)
                # JDW old way
                # testPartId,testSeqBeg,testSeqEnd,testSeqPartType=seqFeature.getAuthPartDetails()
                testPartId, testSeqPartType = seqFeature.getPartInfo()
                # JDW get current part boundaries from internal data structure -
                testSeqBeg = self.__authPartAlignD[testPartId][1]
                testSeqEnd = self.__authPartAlignD[testPartId][2]
                if (self.__verbose):
                    self.__lfh.write("+AlignmentUtils.doAlignment() for testId %s testPartId %r testSeqBeg %d  testSeqEnd %d\n" %
                                     (testSeqId, testPartId, testSeqBeg, testSeqEnd))
            #

            alignTestSeq = []
            for aTup in aL:
                alignTestSeq.append(aTup[1])
            alignTestSeqIdx = self.__annotateAlignmentWithIndex(alignTestSeq, testSeqIdx)
            lenTestAlignment = len(alignTestSeqIdx)
            if (self.__debug):
                self.__formatAlignedPair(self.__lfh, refSeqId, alignRefSeqIdx, testSeqId, alignTestSeqIdx)
            #
            if ((refSeqType in ['auth']) and (testSeqType in ['ref'])):
                for ii in range(testSeqBeg - 1, testSeqEnd):
                    rT = alignRefSeqIdx[ii]
                    tT = alignTestSeqIdx[ii]
                    if ((rT[1] == tT[1]) and (len(rT[5]) > 0)):
                        if (self.__verbose):
                            self.__lfh.write("+AlignmentUtils.doAlignment() clear annotation at position %d  = %s\n" % (ii + 1, rT[5]))
                        alignRefSeqIdx[ii] = (rT[0], rT[1], rT[2], rT[3], rT[4], '')
                        alignTestSeqIdx[ii] = (tT[0], tT[1], tT[2], tT[3], tT[4], '')

            alignSeq = []
            conflictL = [(0, '')] * len(aL)
            numConflicts = 0
            for ii, aTup in enumerate(aL):
                #
                # get the index in the reference sequence
                refTupIdx = alignRefSeqIdx[ii]

                #
                # JDW  2013-11-16 filter conflicts that are outside part range
                # Skip conflicts if we are outside of the reference sequence index range for our sequence part.
                # if (((refSeqType in ['auth']) and (testSeqType in ['ref'])) and ((refIdx < refSeqBeg) or (refIdx > refSeqEnd))):
                #    continue
                #
                refCompId = aTup[0]
                testCompId = aTup[1]

                if ((refSeqType in ['auth']) and (refCompId == self.__gapSymbol) and (testCompId != self.__gapSymbol)):
                    conflictL[ii] = (1, '')
                    conflictRefL[ii] = (1, '')
                    numConflicts += 1
                    continue
                #
                #   Skip dual gaps
                if ((refCompId == self.__gapSymbol) and (testCompId == self.__gapSymbol)):
                    continue
                #
                # We can only use the sequence array index as the residue index may have undefined values for deletions
                # refIdx=refTupIdx[3]+1
                refIdx = ii + 1
                # --------------------------------------------------
                #      --  Conflict assignments performed here --
                #

                # Skip conflicts if we are outside of the sequence index range for our sequence part.
                if ((numParts > 1) and (((refSeqType in ['auth']) and (testSeqType in ['ref'])) and ((refIdx < testSeqBeg) or (refIdx > testSeqEnd)))):
                    continue

                isConflict, refConflict, testConflict = self.__assignConflict(refSeqType, refCompId, testSeqType, testCompId, ii, len(aL))
                if isConflict:
                    # assign conflict to proper sequence part
                    conflictL[ii] = testConflict
                    conflictRefL[ii] = refConflict
                    numConflicts += 1
                # else:
                #    conflictRefL[ii]=(0,'')

            if (self.__verbose):
                self.__lfh.write("+AlignmentUtils.doAlignment() %s conflicts %d length alignment %d length indexed test sequence %d\n"
                                 % (testSeqId, numConflicts, len(alignTestSeqIdx), len(testSeqIdx)))

            if ((refSeqType in ['auth']) and (testSeqType in ['ref'])):
                #self.__annotateAlignmentWithDefaultConflicts(alignRefDbSeqIdx=alignTestSeqIdx, alignAuthSeqIdx=alignRefSeqIdx,authSeqBeg=testSeqBeg,authSeqEnd=testSeqEnd,authSeqPartId=testPartId)
                # JDW -  Exchanged the sequence boundaries here for the reference sequence -- the above was wrong -
                #
                # Need to get authSeqBeg & End for the current sequence part for testPartId -- in the alignment coordinate system!!
                #
                refSeqPartBeg = self.__authPartAlignD[testPartId][1]
                refSeqPartEnd = self.__authPartAlignD[testPartId][2]
                self.__annotateAlignmentWithDefaultConflicts(alignRefDbSeqIdx=alignTestSeqIdx, alignAuthSeqIdx=alignRefSeqIdx, authSeqBeg=refSeqPartBeg,
                                                             authSeqEnd=refSeqPartEnd, authSeqPartId=testPartId, authDefinedMutations=authDefinedMutations)

            alignSeqList.append([testSeqId, testSeqLab, alignTestSeqIdx, conflictL, testSeqFD])

        #
        self.__clearExcludedConflicts(alignSeqList)
        #
        if (self.__verbose):
            for aSeq in seqList:
                self.__lfh.write("+AlignmentUtils.doAlignment() %s alignment length %d reference length %d\n" %
                                 (testSeqId, len(aSeq[2]), lenRefAlignment))

            numConflicts = 0
            for conflictTup in alignSeqList[0][3]:
                if conflictTup[0] != 0:
                    numConflicts += 1
            self.__lfh.write("+AlignmentUtils.doAlignment() Total conflicts of aligned sequences with reference %d\n" % numConflicts)

        return alignSeqList

    def filterSequenceIdList(self, seqIdList=None, groupId=None):
        """  Select the sequence id in the input alignment list corresponding to the input align group/entity id.
        """
        if self.__sds is None:
            self.__sds = SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)

        oL = []
        if groupId is None or seqIdList is None or len(seqIdList) < 1:
            return ([])
        #
        instIdList = self.__sds.getGroup(groupId)
        #
        if (self.__verbose):
            self.__lfh.write("\n+AlignmentUtils.filterSequenceIdList()  align/entity group %s instance list %r input sequence list %r\n" %
                             (groupId, instIdList, seqIdList))
        sLab = SequenceLabel()
        for seqId in seqIdList:
            sLab.unpack(seqId)
            seqType = sLab.getSequenceType()
            seqInstId = sLab.getSequenceInstId()
            if seqType in ['ref', 'auth'] and seqInstId == groupId:
                oL.append(seqId)
            elif seqType in ['xyz'] and seqInstId in instIdList:
                oL.append(seqId)

        if (self.__verbose):
            self.__lfh.write("\n+AlignmentUtils.filterSequenceIdList() returning filtered list %r\n" % oL)

        return oL

    def getSelfReferencePartIdList(self, seqIdList=None, groupId=None):
        """  Select the input sequence id list for group/entity id self referenced parts.
        """
        oL = []
        if groupId is None or seqIdList is None or len(seqIdList) < 1:
            return ([])
        #
        sLab = SequenceLabel()
        for seqId in seqIdList:
            sLab.unpack(seqId)
            seqType = sLab.getSequenceType()
            seqInstId = sLab.getSequenceInstId()
            seqPartId = sLab.getSequencePartId()
            if seqType in ['selfref'] and seqInstId == groupId:
                oL.append(seqPartId)

        if (self.__verbose):
            self.__lfh.write("\n+AlignmentUtils.getSelfReferencePartIdList() returning self-reference part id list %r\n" % oL)

        return oL

    def setExcludedPartIdList(self, partIdList):
        self.__excludedPartIdList = partIdList

    def __getAlignmentReference(self, alignIdList=None):
        """
        Select the alignment reference from the input list of sequences identifiers (e.g. ref_A_1_1_1).

        Return the tuple of the seqId of the reference and the list of test sequences

        """

        if alignIdList is None or len(alignIdList) < 2:
            return (None, [])
        #
        if (self.__verbose):
            self.__lfh.write("\n+AlignmentUtils.__getAlignmentReference() Alignment ID size - %d\n" % len(alignIdList))
            for id in alignIdList:
                self.__lfh.write("   %s\n" % id)

        # take first auth sequence -
        refId = None
        for id in alignIdList:
            if id.startswith("auth"):
                refId = id
                break
        # next try the first ref sequence
        if refId is None:
            for id in alignIdList:
                if id.startswith("ref"):
                    refId = id
                    break
        # finally - take the first of what is selected -
        if refId is None:
            refId = alignIdList[0]
        #
        # JDW 2014-July-03
        # Order the remainder:  auth, ref, xyz
        #
        oL = []
        for seqType in ['auth', 'ref', 'xyz']:
            for id in alignIdList:
                if id != refId and id.startswith(seqType):
                    oL.append(id)
        if (self.__verbose):
            self.__lfh.write("\n+AlignmentUtils.__getAlignmentReference() Reference ID  %s\n" % refId)
            for id in oL:
                self.__lfh.write(" Test Id  %s\n" % id)

        return refId, oL

    def __annotateAlignmentWithIndex(self, s3L, s3WithIndexL, clearSelectedComments=False):
        """
        This internal method adds  one-letter code, original residue index, position in the original sequence, position in
        alignment and annotation comment -

        On input:
                    s3L - is the aligned sequence as a list of 3-letter-codes or gap-symbols.

           s3WithIndexL - original sequence with indexing and annotation details as stored in the sequence data store -

                          Stored indexed sequence have storage model -

                          (3-letter-code, original residue index, comment/details, alignment index )

        Returns a list of tuples containing -

        (one-letter-code, 3-letter-code, original residue index, position in sequence (0 based), position in alignment (0-based), comment)

        var detailList = {'engineered mutation': 'engineered mutation','cloning artifact': 'cloning artifact',
                      'variant':'variant','expression tag':'expression tag', 'insertion':'insertion',
                      'deletion':'deletion', 'microheterogeneity':'microheterogeneity', 'chromophore':'chromophore',
                      'linker':'linker', 'conflict':'conflict', 'acetylation':'acetylation', 'amidation':'amidation',
                      'initiating methionine':'initiating methionine', 'modified residue': 'modified residue',
        'microheterogeneity/modified residue':'microheterogeneity/modified residue',};
        """
        ir = 0
        annotL = []
        for idx in range(0, len(s3L)):

            if ((ir < len(s3WithIndexL)) and (s3WithIndexL[ir][0] == s3L[idx])):
                if (clearSelectedComments and (s3WithIndexL[ir][2] in ['linker', 'expression tag', 'conflict', 'cloning artifact', 'insertion', 'engineered mutation'])):
                    comment = ''
                else:
                    comment = s3WithIndexL[ir][2]
                annotL.append((self.__srd.cnv3To1(s3L[idx]), s3L[idx], s3WithIndexL[ir][1], ir, idx, comment))
                ir += 1
            elif s3L[idx] == self.__gapSymbol:
                annotL.append((self.__gapSymbol, self.__gapSymbol, '', '', idx, ''))
            else:
                self.__lfh.write("+AlignmentUtils.__annotateAlignmentWithIndex() error at index %d %s %d %d\n "
                                 % (idx, s3L[idx], ir, len(s3WithIndexL)))

        return annotL

    def __assignConflict(self, refSeqType, refCompId, testSeqType, testCompId, idxAlign, alignLength):
        """
           On input:
                      refSeqType,testSeqType = ['auth'|'ref'|'xyz']
                      refCompId,testCompId   =  component 3-letter-codes

             JDW not used
                      idxAlign    = index in alignment
                      alignLength = length of alignment

           Returns:  isConflict, refConflictTup (code,correction), testConflictTup (code,correction)

           Methionine --  Any reference methione is candidate for met/mse -


        """
        refConflict = (0, '')
        testConflict = (0, '')
        isConflict = False

        if ((testSeqType.lower() in ['xyz']) and (testCompId == self.__gapSymbol) and (refCompId == 'MET')):
            return True, (9, 'MSE'), testConflict

        if ((testSeqType.lower() in ['xyz']) and (testCompId == self.__gapSymbol)):
            return isConflict, refConflict, testConflict

        if ((refSeqType in ['auth']) and (testSeqType in ['ref']) and (refCompId != testCompId) and (testCompId == self.__gapSymbol)):
            refConflict, testConflict = self.__assignConflictType(refCompId, refSeqType, testCompId, testSeqType)
            isConflict = True
        elif ((refCompId != testCompId) and (refCompId == self.__gapSymbol)):
            refConflict, testConflict = self.__assignConflictType(refCompId, refSeqType, testCompId, testSeqType)
            isConflict = True
        elif (refCompId != testCompId) and not ((refCompId == self.__gapSymbol) or (testCompId == self.__gapSymbol)):
            refConflict, testConflict = self.__assignConflictType(refCompId, refSeqType, testCompId, testSeqType)
            isConflict = True
        else:
            pass
        return isConflict, refConflict, testConflict

    def __isMatched(self, conflictTup):
        return (conflictTup[0] == 0)

    def __assignConflictType(self, r3Ref, refSeqType, r3Test, testSeqType):
        """ Assign conflict type as  -

            0 - None
            1 - non-specific in ref sequence
            2 - non-specific in test sequence

            3 - r3Ref = not gap r3Test=gap
            4 - r3Ref = gap     r3Test=not gap

            5 - r3Ref = GLU  r3Test = GLN
            6 - r3Ref - ASP  r3Test = ASN

            7 - r3Ref = Any  r3Test = ALA/GLY
            8 - r3Ref = MET  r3Test = MSE

            and return a tuple of conflict integer type and
            the likely correction.  ref,test   (code,r3),(code,r3)

        """
        if (r3Ref == r3Test):
            return ((0, r3Ref), (0, r3Ref))
        elif (r3Ref == 'MET' and r3Test == self.__gapSymbol):
            return ((9, 'MSE'), (0, r3Test))
        elif (r3Ref != self.__gapSymbol and r3Test == self.__gapSymbol):
            return ((1, r3Ref), (3, r3Ref))
        elif (r3Ref == self.__gapSymbol and r3Test != self.__gapSymbol):
            return ((1, r3Ref), (4, r3Test))
        elif (r3Ref == 'GLU' and r3Test == 'GLN'):
            if (testSeqType == 'xyz'):
                return ((1, r3Ref), (5, r3Ref))
            else:
                return ((1, r3Ref), (10, r3Test))
            #
        elif (r3Ref == 'ASP' and r3Test == 'ASN'):
            if (testSeqType == 'xyz'):
                return ((1, r3Ref), (6, r3Ref))
            else:
                return ((1, r3Ref), (10, r3Test))
            #
        elif (r3Test == 'ALA' or r3Test == 'GLY'):
            if (testSeqType == 'xyz'):
                return ((1, r3Ref), (7, r3Ref))
            else:
                return ((1, r3Ref), (10, r3Test))
            #
        elif (r3Ref == 'MET' and r3Test == 'MSE'):
            return ((8, r3Test), (1, r3Test))
        elif (r3Ref == 'MSE' and r3Test == 'MET'):
            return ((1, r3Ref), (10, r3Test))
        elif (testSeqType == 'xyz'):
            return ((1, r3Ref), (2, r3Ref))
        #
        return ((1, r3Ref), (10, r3Test))

    def __XassignConflictType(self, r3Ref, r3Test):
        """ Assign conflict type as  -

            0 - None
            1 - non-specific in ref sequence
            2 - non-specific in test sequence

            3 - r3Ref = not gap r3Test=gap
            4 - r3Ref = gap     r3Test=not gap

            5 - r3Ref = GLU  r3Test = GLN
            6 - r3Ref - ASP  r3Test = ASN

            7 - r3Ref = Any  r3Test = ALA/GLY
            8 - r3Ref = MET  r3Test = MSE

            and return a tuple of conflict integer type and
            the likely correction.  (code,r3)

        """
        if (r3Ref == r3Test):
            return (0, r3Ref)
        elif (r3Ref != self.__gapSymbol and r3Test == self.__gapSymbol):
            return (3, r3Ref)
        elif (r3Ref == self.__gapSymbol and r3Test != self.__gapSymbol):
            return (4, r3Test)
        elif (r3Ref == 'GLU' and r3Test == 'GLN'):
            return (5, 'GLU')
        elif (r3Ref == 'ASP' and r3Test == 'ASN'):
            return (6, 'ASP')
        elif (r3Test == 'ALA' or r3Test == 'GLY'):
            return (7, r3Ref)
        elif (r3Ref == 'MET' and r3Test == 'MSE'):
            return (8, 'MSE')

        return (2, r3Ref)

    def __annotateAlignmentWithDefaultConflicts(self, alignRefDbSeqIdx, alignAuthSeqIdx, authSeqBeg, authSeqEnd, authSeqPartId, authDefinedMutations):
        """
        Add default 'expression tag annotation comments' for leading or trailing gap residues, and
            comments for any conflicts which can be attributed to a residue modification.

        On input the aligned sequences are -
               alignAuthSeqIdx[]  contains the aligned author/sample sequence
               alignRefDbSeqIdx[] contains one of the aligned references seqeuence database sequences

        The sequence positions (1-indexed) define the boundaries of the author/sample sequence to annotate.

        The input aligned sequence lists contain the following index details -

        (one-letter-code, 3-letter-code, original label residue index, position in sequence, position in alignment, comment)

        On return --

        Input sequences are updated in place with annotations

        Revised - JDW - 21-May-2013 - make the annotation assignment scan span the full alignment.

        """
        lastPart = True
        if len(self.__authPartIdList) > 0:
            lastPart = (authSeqPartId == self.__authPartIdList[-1])

        firstPart = True
        if len(self.__authPartIdList) > 0:
            firstPart = (authSeqPartId == self.__authPartIdList[0])

        lenAlign = len(alignAuthSeqIdx)

        seqBeg = authSeqBeg
        seqEnd = authSeqEnd

        if (self.__verbose):
            self.__lfh.write("\n+AlignmentUtils.__annotateAlignmentWithDefaultConflicts() STARTS with seqPartId %d (first %r last %r) authSeqBeg %d authSeqEnd %d aligned sequence lengths auth %d refdb %d\n" %
                             (authSeqPartId, firstPart, lastPart, authSeqBeg, authSeqEnd, len(alignRefDbSeqIdx), len(alignAuthSeqIdx)))

        #
        # Remove prior annotations for conflicts that are not present in current alignment -
        #   --   JDW from the author sequence ONLY for the first part  --
        #
        for ii in range(0, lenAlign):
            authTup = alignAuthSeqIdx[ii]
            refDbTup = alignRefDbSeqIdx[ii]
            #
            authCompId = authTup[1]
            refCompId = refDbTup[1]
            #
            if (authCompId == refCompId):
                if firstPart:
                    sList = list(authTup)
                    sList[5] = ''
                    alignAuthSeqIdx[ii] = tuple(sList)

                sList = list(refDbTup)
                sList[5] = ''
                alignRefDbSeqIdx[ii] = tuple(sList)

            # clear prior conflicts
            if ((authCompId != refCompId) and (refCompId == self.__gapSymbol) and (authTup[5] == 'conflict')):
                if firstPart:
                    sList = list(authTup)
                    sList[5] = ''
                    alignAuthSeqIdx[ii] = tuple(sList)

                sList = list(refDbTup)
                sList[5] = ''
                alignRefDbSeqIdx[ii] = tuple(sList)

        #  Look for trailing gaps in the reference sequence
        #  jdw - note this must run to the end of the alignment --
        for ii in range(lenAlign - 1, -1, -1):
            authTup = alignAuthSeqIdx[ii]
            if (authTup[1] == self.__gapSymbol):
                continue

            # authIdx=authTup[3]+1
            authIdx = ii + 1
            if ((authIdx < seqBeg) or (authIdx > seqEnd)):
                continue
            refDbTup = alignRefDbSeqIdx[ii]
            if refDbTup[1] != self.__gapSymbol:
                break

            # if ((refDbTup[1] == self.__gapSymbol) and (len(authTup[5]) < 1)):
            if ((refDbTup[1] == self.__gapSymbol)):
                if lastPart:
                    cType = 'expression tag'
                else:
                    cType = 'linker'

                # do not overwrite existing annotations
                if (len(authTup[5]) < 1):
                    sList = list(authTup)
                    sList[5] = cType
                    alignAuthSeqIdx[ii] = tuple(sList)
                if (len(refDbTup[5]) < 1):
                    sList = list(refDbTup)
                    sList[5] = cType
                    alignRefDbSeqIdx[ii] = tuple(sList)

        #  Look for leading gaps in the reference sequence
        # jdw - note this must run from the beginning of the alignment --
        for ii in range(0, lenAlign):

            authTup = alignAuthSeqIdx[ii]
            if (authTup[1] == self.__gapSymbol):
                continue
            # authIdx=authTup[3]+1
            authIdx = ii + 1
            if ((authIdx < seqBeg) or (authIdx > seqEnd)):
                continue

            refDbTup = alignRefDbSeqIdx[ii]
            if refDbTup[1] != self.__gapSymbol:
                break

            if ((refDbTup[1] == self.__gapSymbol)):
                if firstPart:
                    if ((ii == 0) and (authTup[1] in ['MET', 'MSE'])):
                        cType = 'initiating methionine'
                    else:
                        cType = 'expression tag'
                else:
                    cType = 'linker'
                if (len(authTup[5]) < 1):
                    sList = list(authTup)
                    sList[5] = cType
                    alignAuthSeqIdx[ii] = tuple(sList)
                if (len(refDbTup[5]) < 1):
                    sList = list(refDbTup)
                    sList[5] = cType
                    alignRefDbSeqIdx[ii] = tuple(sList)

        # Look for internal deletions in author sequence and author residue modifications -
        #
        # scan status ( before, in, last, after)
        scanStatus = 'before'

        for ii in range(0, lenAlign):
            authTup = alignAuthSeqIdx[ii]
            refDbTup = alignRefDbSeqIdx[ii]
            #
            authCompId = authTup[1]
            # authSeqPos=authTup[3]
            authSeqPos = ii + 1
            refCompId = refDbTup[1]
            #
            if ((scanStatus == 'before') and (authCompId != self.__gapSymbol) and (authSeqPos >= seqBeg)):
                scanStatus = 'in'

            if (scanStatus == 'last'):
                scanStatus = 'after'

            if ((scanStatus == 'in') and (authCompId != self.__gapSymbol) and (authSeqPos == seqEnd)):
                scanStatus = 'last'

            if ((authCompId == self.__gapSymbol) and (authCompId != refCompId) and (len(authTup[5]) < 1)):
                sList = list(authTup)
                sList[5] = 'deletion'
                alignAuthSeqIdx[ii] = tuple(sList)
                sList = list(refDbTup)
                sList[5] = 'deletion'
                alignRefDbSeqIdx[ii] = tuple(sList)
                continue

            #
            # if (self.__debug):
            # if authCompId != refCompId:
            #    self.__lfh.write("AlignmentUtils.__annotateAlignmentWithDefaultConflicts() ii %r scanStatus %s authSeqPos %r seqBeg %r seqEnd %r authCompId %r refCompId %r\n" %
            #                     (ii,scanStatus,authSeqPos,seqBeg,seqEnd,authCompId,refCompId) )
            #

            if scanStatus not in ['in', 'last']:
                continue

            if self.__ccU.isModificationOf(refCompId, authCompId):
                sList = list(authTup)
                tS = sList[5]
                if sList[5] is not None and sList[5] in ['microheterogeneity']:
                    tS = 'microheterogeneity/modified residue'
                elif sList[5] is not None and sList[5] in ['microheterogeneity/modified residue']:
                    pass
                else:
                    tS = 'modified residue'
                sList[5] = tS
                alignAuthSeqIdx[ii] = tuple(sList)
                #
                sList = list(refDbTup)
                sList[5] = tS
                alignRefDbSeqIdx[ii] = tuple(sList)
            #
            #  Handle any internal insertions lacking annotation -
            elif ((authCompId != refCompId) and (refCompId == self.__gapSymbol) and (len(authTup[5]) < 1)):
                sList = list(authTup)
                sList[5] = 'insertion'
                alignAuthSeqIdx[ii] = tuple(sList)
                sList = list(refDbTup)
                sList[5] = 'insertion'
                alignRefDbSeqIdx[ii] = tuple(sList)

            #
            #  Handle any internal conflicts lacking annotation -
            elif ((authCompId != refCompId) and (len(authTup[5]) < 1)):
                val = 'conflict'
                mut1 = refDbTup[0] + refDbTup[2] + authTup[0]
                mut2 = refDbTup[0] + authTup[2] + authTup[0]
                if (mut1 in authDefinedMutations) or (mut2 in authDefinedMutations):
                    val = 'engineered mutation'
                #
                sList = list(authTup)
                sList[5] = val
                alignAuthSeqIdx[ii] = tuple(sList)
                sList = list(refDbTup)
                sList[5] = val
                alignRefDbSeqIdx[ii] = tuple(sList)

    #
    #   Alignment data I/O --

    def saveAlignmentIdList(self, alignIdList):
        try:
            idListPath = os.path.join(self.__sessionPath, "align-id-list.txt")
            ofh = open(idListPath, 'w')
            ofh.write("%s\n" % ','.join(alignIdList))
            ofh.close()
            return True
        except:
            if (self.__verbose):
                self.__lfh.write("+AlignmentUtils.__saveAlignmentList() failed for list = %r\n" % alignIdList)
                traceback.print_exc(file=self.__lfh)

    def readAlignmentIdList(self):
        """ Return the list of sequence ids in the current alignment.
        """
        try:
            idListPath = os.path.join(self.__sessionPath, "align-id-list.txt")
            ifh = open(idListPath, 'r')
            tS = str(ifh.read()).strip()
            ifh.close()
            return tS.split(',')
        except:
            traceback.print_exc(file=self.__lfh)
            return []

    def storeAlignmentData(self, alignSeqList=None, identifier='anon', entityId='1'):
        """  Store the alignment data in the input alignSeqList = [[],...] in the alignment data store -

        """
        #
        fn = self.__getAlignmentStoreFilename(identifier=identifier, entityId=entityId)
        ads = AlignDataStore(reqObj=self.__reqObj, fileName=fn, verbose=self.__verbose)
        ads.set(alignSeqList)
        ads.serialize()

    def getAlignmentData(self, identifier='anon', entityId='1'):
        """  Restore the current alignment from the alignment data store -

        """
        #
        if (self.__verbose):
            self.__lfh.write("\n+AlignmentUtils.getAlignmentData() checking for stored alignment data for data set %s entity %r\n" % (identifier, entityId))

        fn = self.__getAlignmentStoreFilename(identifier=identifier, entityId=entityId)
        ads = AlignDataStore(reqObj=self.__reqObj, fileName=fn, verbose=self.__verbose, log=self.__lfh)
        alignIdList = ads.getAlignIdList()
        alignSeqList = ads.get()
        #
        if (self.__verbose):
            self.__lfh.write("+AlignmentUtils.getAlignmentData() used stored alignment in file %s\n" % fn)
            self.__lfh.write("+AlignmentUtils.getAlignmentData() entity id  %s\n" % entityId)
            self.__lfh.write("+AlignmentUtils.getAlignmentData() stored alignment length %d\n" % len(alignIdList))
            self.__lfh.write("+AlignmentUtils.getAlignmentData() stored alignment ids    %r\n" % alignIdList)
        #
        return alignIdList, alignSeqList

    def __getAlignmentStoreFilename(self, identifier, entityId='1', fileSource='session'):
        """ Add a tag to the alignment data filename to differentiate multiple active alignment views. -
        """
        fn = os.path.basename(self.__pI.getSequenceAlignFilePath(dataSetId=identifier, entityId=entityId, fileSource=fileSource, versionId="latest"))
        return fn

    def formatAlignment(self, io, alignSeqList):
        """
        """
        alignLength = len(alignSeqList[0][2])
        io.write("\n+AlignmentUtils.__formatAlignment() - alignment list data - length = %d\n" % alignLength)
        io.write("                 Key: one-letter-code (3-letter-code) label residue index (position in sequence, position in alignment) conflict-comment\n")
        for aPos in range(0, alignLength):
            io.write("%5d:  " % aPos)
            for aTup in alignSeqList:
                id = aTup[0]
                sLabel = aTup[1]
                aL = aTup[2]
                conflictL = aTup[3]
                # (one-letter-code, 3-letter-code, original label residue index, position in sequence,position in alignment )
                rT = aL[aPos]
                io.write("%s(%3s) %5r(%5r,%5r)(%r-%-25r) -+- " % (rT[0], rT[1], str(rT[2]), str(rT[3]), str(rT[4]), conflictL[aPos][0], str(rT[5])))
            io.write("\n")
        io.flush()

    def __formatAlignedPair(self, io, refSeqId, alignRefSeqIdx, testSeqId, alignTestSeqIdx):
        """
        """
        alignLength = len(alignRefSeqIdx)
        io.write("+AlignmentUtils.__formatAlignedPair() - ref = %s  test = %s  alignment list data - length = %d\n" % (refSeqId, testSeqId, alignLength))
        io.write("                 Key: one-letter-code (3-letter-code) label residue index (position in sequence, position in alignment) conflict-comment\n")
        for aPos in range(0, alignLength):
            io.write("%5d:  " % aPos)
            # (one-letter-code, 3-letter-code, original label residue index, position in sequence,position in alignment )
            rT = alignRefSeqIdx[aPos]
            io.write("%s(%3s) %4s(%4r,%4r)(%-25s) -+- " % (rT[0], rT[1], rT[2], rT[3], rT[4], str(rT[5])))
            rT = alignTestSeqIdx[aPos]
            io.write("%s(%3s) %4s(%4r,%4r)(%-25s) -+- " % (rT[0], rT[1], rT[2], rT[3], rT[4], str(rT[5])))
            io.write("\n")

    def __writeAlignSeqInfo(self, refSeqId, refSeq3, refSeqType, seqList):
        testSeqPath = os.path.join(self.__sessionPath, "testSeq.h")
        ofh = open(testSeqPath, 'w')
        self.__writeSeq(ofh, refSeqId, refSeq3, '"')
        for aSeq in seqList:
            authLabel = SequenceLabel()
            authLabel.unpack(aSeq[0])
            (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = authLabel.get()
            self.__writeSeq(ofh, aSeq[0], aSeq[2], '"')
            if (refSeqType != 'auth') or (seqType != 'ref'):
                self.__writeSeq(ofh, aSeq[0] + "_link", aSeq[5], '')
            #
        #
        ofh.write("\n\ntypedef struct {\n")
        ofh.write("       const char*  seq_id;\n")
        ofh.write("       const char** seq;\n")
        ofh.write("       const int*   link;\n")
        ofh.write("       int   len_seq;\n")
        ofh.write("       int   len_link;\n")
        ofh.write("       int   range_start;\n")
        ofh.write("       int   range_end;\n")
        ofh.write("} SEQS;")
        ofh.write("\n\n#define NUM_SEQS  %d\n" % (len(seqList) + 1))
        ofh.write("\n\nstatic const SEQS seqs_def[NUM_SEQS] = {\n")
        ofh.write('       { "%s", %s, NULL, %d, 0, 0, 0 },\n' % ( refSeqId, refSeqId, len(refSeq3) ))
        for aSeq in seqList:
            authLabel = SequenceLabel()
            authLabel.unpack(aSeq[0])
            (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = authLabel.get()
            if (refSeqType == 'auth') and (seqType == 'ref'):
                ofh.write('       { "%s", %s, NULL, %d, 0, %d, %d },\n' % ( aSeq[0], aSeq[0], len(aSeq[2]), self.__authPartD[seqPartId][1], self.__authPartD[seqPartId][2]))
            else:
                ofh.write('       { "%s",  %s, %s, %d, %d, 0, 0 },\n' % ( aSeq[0], aSeq[0], aSeq[0] + "_link", len(aSeq[2]), len(aSeq[5]) ))
            #
        #
        ofh.write("};\n\n\n")
        ofh.close()

    def __writeSeq(self, io, Id, aSeq, quote):
        """
        """
        io.write("\n\n")
        if quote == '"':
            io.write("static const char *%s[%d] = {\n" % (Id, len(aSeq)))
        else:
            io.write("static const int %s[%d] = {\n" % (Id, len(aSeq)))
        #
        count = 0
        first = True
        for s in aSeq:
            if not first:
                io.write(", ")
            #
            if count == 20:
                io.write("\n")
                count = 0
            #
            io.write("%s%s%s" % ( quote, str(s), quote ))
            first = False
            count += 1
           
        #
        io.write("\n};\n")

    def __dumpAlignSeqInput(self, io, refDef, seqList):
        """ Output the original sequence lists and the resulting alignment lists -
        """
        refSeqIdx = refDef[3]
        refSeqId = refDef[0]
        io.write("\n+AlignmentUtils.__dumpAlignSeqInput() -------------------------------------------------------------------------\n")
        io.write("+AlignmentUtils.__dumpAlignSeqInput Reference sequence() %s full length %d\n" % (refSeqId, len(refSeqIdx)))
        for tup in refSeqIdx:
            io.write("   + %8s %8s %8s \n" % (tup[0], tup[1], tup[2]))

        for aSeq in seqList:
            io.write("+AlignmentUtils.__dumpAlignSeqInput() ------------------------------\n")
            io.write("+AlignmentUtils.__dumpAlignSeqInput() Test sequence %s full length %d\n" % (aSeq[0], len(aSeq[3])))
            for tup in aSeq[3]:
                io.write("   + %8s %8s %8s \n" % (tup[0], tup[1], tup[2]))
        io.flush()

    ##
    def transferComments(self, alignSeqList):
        """ For the input aligned sequence list propogate any comment details to corresponding
            sequences in the sequence data store.

            On input:

            alignSeqList [[Sequence identifier (ie. auth_A_1_1_1 as used in sequence data store),
                           SequenceLabel() object,
                           Aligned sequence with index details,
                           Conflict flag list (bool),
                           Feature dictionary for input sequence],,,]


            Each aligned sequences with index details has the storage model -

            (one-letter-code, 3-letter-code, original label residue index, position in sequence, position in alignment, comment )

            Stored sequences in the sequence data store have the simpler storage model -

             (3-letter-code, original residue index, comment/details, alignment index )

            Annotations in alignment comment field <conflict annotations> are transfered to the comment
            field of the associated sequence in the sequence data store.

        """
        sLabel = SequenceLabel()
        if self.__sds is None:
            self.__sds = SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        #
        # For the sequence in the alignment -
        alignRef = True
        for alignTup in alignSeqList:
            seqId = alignTup[0]
            alSeqIdx = alignTup[2]
            if (self.__verbose):
                self.__lfh.write("+AlignmentUtils.transferComments() updating comments for seqId = %r length %d\n" % (seqId, len(alSeqIdx)))
            #
            # Recover the sequence from the data store -
            sLabel.unpack(seqId)
            (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = sLabel.get()
            #
            seqIdxNext = []
            #
            # JDW oct-22 -- Review preserving gaps -
            # if seqType in ['xyz'] or ((seqType in ['auth']) and alignRef) :
            if seqType in ['xyz']:
                # preserve gaps in coordinate  -
                for aTup in alSeqIdx:
                    seqIdxNext.append((aTup[1], aTup[2], aTup[5], aTup[4]))
            else:
                #
                for aTup in alSeqIdx:
                    if aTup[1] != self.__gapSymbol:
                        seqIdxNext.append((aTup[1], aTup[2], aTup[5], aTup[4]))

            if self.__verbose:
                self.__lfh.write("+AlignmentUtils.transferComments() updated sequence %s length %d version %s\n"
                                 % (seqId, len(seqIdxNext), seqVersion))
            if self.__debug:
                for aTup in seqIdxNext:
                    if aTup[2] is not None and len(aTup[2]) > 0:
                        self.__lfh.write("      +++transferComments() %s %s %s %s \n" % ((aTup[0], aTup[1], aTup[3], aTup[2])))

            # Update revised auth sequence in the sequence  data store -
            if seqType == "auth":
                # JDW - new get parts for only the relevant version -----
                partIdList = self.__sds.getPartIdsForVersion(seqInstId, seqType=seqType, altId=seqAltId, version=seqVersion)
                if (self.__verbose):
                    self.__lfh.write("+AlignmentUtils.transferComments() part id list %r\n" % partIdList)
                # All of the parts of author sequence are complete and common -
                for pId in partIdList:
                    if (self.__verbose):
                        self.__lfh.write("+AlignmentUtils.transferComments() store sequence for part id  %r\n" % pId)
                    self.__sds.setSequence(seqIdxNext, seqId=seqInstId, seqType=seqType, partId=pId, altId=seqAltId, version=seqVersion)
            else:
                self.__sds.setSequence(seqIdxNext, seqId=seqInstId, seqType=seqType, partId=seqPartId, altId=seqAltId, version=seqVersion)

            #
            alignRef = False

        self.__sds.serialize()
        #

    def __PREVannotateAlignmentWithDefaultConflicts(self, alignRefDbSeqIdx, alignAuthSeqIdx, authSeqBeg, authSeqEnd):
        """
        Add default 'comments' for expression tags for any leading or trailing gap residues.

        On input the aligned sequences are -
               alignAuthSeqIdx[]  contains the aligned author/sample sequence
               alignRefDbSeqIdx[] contains one of the aligned references seqeuence database sequences

        The sequence positions (1-indexed) define the boundaries of the author/sample sequence to annotate.

        The input aligned sequence lists contain the following index details -

        (one-letter-code, 3-letter-code, original label residue index, position in sequence, position in alignment, comment)

        On return --

        Input sequences are updated in place with annotations

        """
        seqBeg = authSeqBeg
        seqEnd = authSeqEnd
        if (self.__verbose):
            self.__lfh.write("AlignmentUtils.__annotateAlignmentWithDefaultConflicts() alignment lengths auth %d refdb %d authSeqBeg %d authSeqEnd %d\n" %
                             (len(alignRefDbSeqIdx), len(alignAuthSeqIdx), authSeqBeg, authSeqEnd))

        # jdw - note this must run to the end of the alignment --
        for ii in range(len(alignRefDbSeqIdx) - 1, -1, -1):
            authTup = alignAuthSeqIdx[ii]
            if (authTup[1] == self.__gapSymbol):
                continue

            authIdx = authTup[3] + 1
            if ((authIdx < seqBeg) or (authIdx > seqEnd)):
                continue
            refDbTup = alignRefDbSeqIdx[ii]
            if refDbTup[1] != self.__gapSymbol:
                break

            if ((refDbTup[1] == self.__gapSymbol) and (len(authTup[5]) < 1)):
                sList = list(authTup)
                sList[5] = 'expression tag'
                alignAuthSeqIdx[ii] = tuple(sList)
                sList = list(refDbTup)
                sList[5] = 'expression tag'
                alignRefDbSeqIdx[ii] = tuple(sList)

        for ii in range(seqBeg - 1, seqEnd):
            authTup = alignAuthSeqIdx[ii]
            if (authTup[1] == self.__gapSymbol):
                continue
            authIdx = authTup[3] + 1
            if ((authIdx < seqBeg) or (authIdx > seqEnd)):
                continue

            refDbTup = alignRefDbSeqIdx[ii]
            if refDbTup[1] != self.__gapSymbol:
                break

            if ((refDbTup[1] == self.__gapSymbol) and (len(authTup[5]) < 1)):
                sList = list(authTup)
                sList[5] = 'expression tag'
                alignAuthSeqIdx[ii] = tuple(sList)
                sList = list(refDbTup)
                sList[5] = 'expression tag'
                alignRefDbSeqIdx[ii] = tuple(sList)

        for ii in range(seqBeg - 1, seqEnd):
            authTup = alignAuthSeqIdx[ii]
            refDbTup = alignRefDbSeqIdx[ii]

            # moved this as deletions will not have an index
            # authIdx=authTup[3]+1
            # if ((authIdx < seqBeg) or (authIdx > seqEnd)):
            #    continue

            if ((authTup[1] == self.__gapSymbol) and (len(authTup[5]) < 1)):
                sList = list(authTup)
                sList[5] = 'deletion'
                alignAuthSeqIdx[ii] = tuple(sList)
                sList = list(refDbTup)
                sList[5] = 'deletion'
                alignRefDbSeqIdx[ii] = tuple(sList)
                continue

            authIdx = authTup[3] + 1
            #
            if (self.__debug):
                self.__lfh.write("+AlignmentUtils.__annotateAlignmentWithDefaultConflicts() ii %r authIdx %r seqBeg %r seqEnd %r refCompId %r authCompId %r\n" %
                                 (ii, authIdx, seqBeg, seqEnd, refDbTup[1], authTup[1]))
            #

            if ((authIdx < seqBeg) or (authIdx > seqEnd)):
                continue

            if self.__ccU.isModificationOf(refDbTup[1], authTup[1]):
                # if ((authTup[1] == 'MSE') and (refDbTup[1] == 'MET')):
                sList = list(authTup)
                sList[5] = 'modified residue'
                alignAuthSeqIdx[ii] = tuple(sList)
                sList = list(refDbTup)
                sList[5] = 'modified residue'
                alignRefDbSeqIdx[ii] = tuple(sList)
            #

    def __annotateAuthOutOfPartComment(self, alignAuthSeqIdx, conflictL):
        """      Return annotation assignments for out-of-part regions of author sequence.
                 (e.g. out-of-part  terminal expression tags or internal linkers)

                input   authSeqwithIndex and ConflictL

                 returns --   True for success or false otherwise
        """
        if (self.__verbose):
            self.__lfh.write(
                "+AlignmentUtils.__annotateAuthOutOfPartComment() starting with aligned sequence length %d and bounds %r\n" %
                (len(alignAuthSeqIdx), self.__authSeqAlignBounds))
        if self.__authSeqAlignBounds is None:
            return False
        #
        for idx in range(1, len(alignAuthSeqIdx) + 1):
            inPart = False
            for partId in self.__authPartIdList:
                if ((idx >= self.__authPartAlignD[partId][1]) and (idx <= self.__authPartAlignD[partId][2])):
                    inPart = True
                    break
            if inPart:
                continue
            # not in part --
            sL = list(alignAuthSeqIdx[idx - 1])
            if idx < self.__authSeqAlignBounds[0] or idx > self.__authSeqAlignBounds[1]:
                sL[5] = 'expression tag'
            else:
                sL[5] = 'linker'

            if (self.__verbose):
                self.__lfh.write("+AlignmentUtils.__annotateAuthOutOfPartComment() assigning %s at position %d\n" % (sL[5], idx))
            alignAuthSeqIdx[idx - 1] = tuple(sL)
            conflictL[idx - 1] = (1, 'out-of-part')
        #
        return True

    def __getClosestKey(self, iVal, mD):
        """  Find the key in the input dictionary with minimum absolute difference
             to the target input value.
        """
        retVal = None
        try:
            iL = mD.keys()
            minV = 100000
            for i in iL:
                iDiff = abs(iVal - i)
                if iDiff < minV:
                    retVal = i
                    minV = iDiff
        except:
            traceback.print_exc(file=self.__lfh)
        return retVal

    #

    def __mapAuthPartDetails(self, authSeqIdx):
        """  Update the internal data structures containing 'auth' entity part details for input sample sequence after
             alignment.   This is necessary if the length of the sample sequence is extended owing to deletions relative
             to the reference sequence or due to other edits --

             authSeqIdx, the aligned sequence records have the structure:

             [(one-letter-code, 3-letter-code, original residue index, position in sequence, position in alignment, comment),...]

        """
        #
        if (self.__verbose):
            self.__lfh.write("\n+AlignmentUtils.__mapAuthPartDetails() Starting with alignment length %d\n" % len(authSeqIdx))
        #
        # Create a mapping dictionary -- residue index -> alignment position 1-N
        #
        mD = {}
        for ii, seqTup in enumerate(authSeqIdx):
            if seqTup[0] != self.__gapSymbol:
                mD[int(seqTup[2])] = ii + 1
        #
        self.__authPartAlignD = {}
        #
        for pId in self.__authPartD.keys():
            pTup = self.__authPartD[pId]
            (authPartId, authSeqBeg, authSeqEnd, authSeqPartType) = pTup
            # include check for missing mapping due to potenitial residue deletion
            if authSeqBeg not in mD:
                seqBeg = self.__getClosestKey(authSeqBeg, mD)
            else:
                seqBeg = mD[authSeqBeg]

            if authSeqEnd not in mD:
                seqEnd = self.__getClosestKey(authSeqEnd, mD)
            else:
                seqEnd = mD[authSeqEnd]

            self.__authPartAlignD[pId] = (authPartId, seqBeg, seqEnd, authSeqPartType)

        #
        # Get the partId list ordered by the begining sequence position.   (aligned coordinate system)
        #
        self.__authSeqAlignBounds = None
        pTupL = self.__authPartAlignD.values()
        pTupL.sort(key=itemgetter(1))
        authSeqPosMin = pTupL[0][1]
        authSeqPosMax = pTupL[-1][2]
        self.__authSeqAlignBounds = (authSeqPosMin, authSeqPosMax)
        if (self.__verbose):
            self.__lfh.write("+AlignmentUtils.__mapAuthPartDetails() part dictionary parts %r\n" % self.__authPartAlignD.items())
            self.__lfh.write("+AlignmentUtils.__mapAuthPartDetails() auth sequence minimum position %r maximum position %r\n" %
                             (self.__authSeqAlignBounds[0], self.__authSeqAlignBounds[1]))
        #

    def __updateAuthPartDetails(self, authSeqIdx):
        """  Update the internal data structures containing 'auth' entity part details for input sample sequence after
             alignment.   This is necessary if the length of the sample sequence is extended owing to deletions relative
             to the reference sequence or due to other edits --

             authSeqIdx, the aligned sequence records have the structure:

             [(one-letter-code, 3-letter-code, original residue index, position in sequence, position in alignment, comment),...]

        """
        #
        if (self.__verbose):
            self.__lfh.write("\n+AlignmentUtils.__updateAuthPartDetails() Starting with alignment length %d\n" % len(authSeqIdx))

        #
        #  Remap sequence endpoints in case alignment edits have changed the sequence length ---
        #
        self.__adjustPartBounds(authSeqIdx)
        #
        # Create a mapping dictionary -- prior residue index -> new position in aligned sequence
        mD = {}
        for ii, seqTup in enumerate(authSeqIdx):
            if seqTup[0] != self.__gapSymbol:
                # mD[int(seqTup[2])]=ii+1
                mD[int(seqTup[2])] = int(seqTup[2])
        #

        #
        for pId in self.__authPartD.keys():
            pTup = self.__authPartD[pId]
            (authPartId, authSeqBeg, authSeqEnd, authSeqPartType) = pTup
            # include check for missing mapping due to potenitial residue deletion
            if authSeqBeg not in mD:
                seqBeg = self.__getClosestKey(authSeqBeg, mD)
            else:
                seqBeg = mD[authSeqBeg]

            if authSeqEnd not in mD:
                seqEnd = self.__getClosestKey(authSeqEnd, mD)
            else:
                seqEnd = mD[authSeqEnd]

            self.__authPartD[pId] = (authPartId, seqBeg, seqEnd, authSeqPartType)

        self.__authPartDetailsWorker()

    def getAuthPartDetails(self, authSeqId):
        self.__buildAuthPartDetails(authSeqId)
        return self.__authPartD

    def __buildAuthPartDetails(self, authSeqId):
        """  Create the internal data structures containing 'auth' entity part details for input sample sequence.

             **** boundaries are coordinate system of the author sequence 1-N ungapped  ****

             Update:

             self.__authPartD[partId]=(int(authPartId),int(authSeqBeg),int(authSeqEnd),str(authSeqPartType.lower()))
             self.__authPartIdList=[] ordered by sequence beginning position -
             self.__authLinkerList=[]
             self.__authSeqBounds=(minPos,maxPos) of author sequence parts - (1-N) ungapped
        """
        self.__authPartD = {}

        if authSeqId.startswith('auth'):
            if (self.__verbose):
                self.__lfh.write("\n+AlignmentUtils.__buildAuthPartDetails() create part dictionary for author reference sequence %s\n" % authSeqId)

            if self.__sds is None:
                self.__sds = SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)

            authLabel = SequenceLabel()
            authLabel.unpack(authSeqId)
            (authSeqType, authSeqInstId, authSeqPartId, authSeqAltId, authSeqVersion) = authLabel.get()

            #
            partIdList = self.__sds.getPartIdsForVersion(authSeqInstId, dataType='sequence', seqType='auth', altId=authSeqAltId, version=authSeqVersion)
            #authFObj=self.__sds.getFeatureObj(authSeqInstId,'auth', partId=authSeqPartId, altId=authSeqAltId,version=authSeqVersion)

            for partId in partIdList:
                authFObj = self.__sds.getFeatureObj(authSeqInstId, 'auth', partId=partId, altId=authSeqAltId, version=authSeqVersion)
                authPartId, authSeqBeg, authSeqEnd, authSeqPartType = authFObj.getAuthPartDetails()
                # JDW
                self.__authPartD[int(partId)] = (int(authPartId), int(authSeqBeg), int(authSeqEnd), str(authSeqPartType.lower()))

            self.__authPartDetailsWorker()

    def __authPartDetailsWorker(self):
        """ Worker method to build internal data structures of author sequence parts ...
        """
        self.__authPartIdList = []
        self.__authLinkerList = []
        self.__authSeqBounds = None

        #
        # Get the partId list ordered by the begining sequence position.
        #
        pTupL = self.__authPartD.values()
        pTupL.sort(key=itemgetter(1))
        for pTup in pTupL:
            self.__authPartIdList.append(pTup[0])

        ##
        try:
            authSeqPosMin = pTupL[0][1]
            authSeqPosMax = pTupL[-1][2]
            self.__authSeqBounds = (authSeqPosMin, authSeqPosMax)
        except:
            self.__authSeqBounds = None
        #
        self.__authLinkerList = []
        if len(self.__authPartIdList) > 1:
            pTup = self.__authPartD[self.__authPartIdList[0]]
            seqEnd = pTup[2]
            for partId in self.__authPartIdList[1:]:
                pTup = self.__authPartD[partId]
                if ((pTup[1] - seqEnd) > 1):
                    for ii in range(seqEnd + 1, pTup[1]):
                        self.__authLinkerList.append(ii)

        if (self.__verbose):
            self.__lfh.write("+AlignmentUtils.__authPartDetailsWorker() part dictionary parts %r\n" % self.__authPartD.items())
            self.__lfh.write("+AlignmentUtils.__authPartDetailsWorker() ordered part list  %r\n" % self.__authPartIdList)
            self.__lfh.write("+AlignmentUtils.__authPartDetailsWorker() authLinkerList %r\n" % (self.__authLinkerList))
            self.__lfh.write("+AlignmentUtils.__authPartDetailsWorker() auth sequence minimum position %r maximum position %r\n" %
                             (self.__authSeqBounds[0], self.__authSeqBounds[1]))

    def __adjustPartBounds(self, authSeqIdx):
        """   Adjust __authPartD to correspond with the current sequence endpoints in the edited aligned
              sequence.

              Input is the current aligned author sequence -
        """
        #
        if (self.__verbose):
            self.__lfh.write("\n+AlignmentUtils.__adjustPartBounds() alignment length %d\n" % len(authSeqIdx))

        minIdx = 10000000
        maxIdx = -100
        for seqTup in authSeqIdx:
            if seqTup[0] != self.__gapSymbol:
                if int(seqTup[2]) < minIdx:
                    minIdx = int(seqTup[2])
                if int(seqTup[2]) > maxIdx:
                    maxIdx = int(seqTup[2])
        #
        # Get the partId list ordered by the begining sequence position.
        #
        pTupL = self.__authPartD.values()
        pTupL.sort(key=itemgetter(1))

        if minIdx != pTupL[0][1]:
            pId = pTupL[0][0]
            seqBeg = pTupL[0][1]
            seqEnd = pTupL[0][2]
            pType = pTupL[0][3]
            self.__authPartD[pId] = (pId, minIdx, seqEnd, pType)

        if maxIdx != pTupL[-1][2]:
            pId = pTupL[-1][0]
            seqBeg = pTupL[-1][1]
            seqEnd = pTupL[-1][2]
            pType = pTupL[-1][3]
            self.__authPartD[pId] = (pId, seqBeg, maxIdx, pType)

        if (self.__verbose):
            self.__lfh.write("+AlignmentUtils.__adjustPartBounds() auth sequence part dictionary parts %r\n" % self.__authPartD.items())
            self.__lfh.write("+AlignmentUtils.__adjustPartBounds() auth sequence minimum position %r maximum position %r\n" % (minIdx, maxIdx))

    def reAnnotateAlignment(self, alignSeqList):
        """ reAnnotate the input alignment

            refDef (refId refLabelObj, ref3L, refSeqIdx, refSeqFD)

            seqList= [[seqId,seqLabelObj,seq3L,seqLWithIndex,seqFD,typicalLink]]

            alignSeqList [[Sequence identifier (ie. auth_A_1_1_1 as used in sequence data store),
                           SequenceLabel() object,
                           Aligned sequence with index details,
                           Conflict flag list (bool),
                           Feature dictionary for input sequence],,,]

        where aligned sequences with index is a list of tuples of the form -

        (one-letter-code, 3-letter-code, original residue index, position in sequence, position in alignment, comment)

        operates on input alignSeqList in place --

        """
        #
        seqFeature = SequenceFeature()
        #
        # Reference sequence is the first sequence on the align list --
        #
        alignRefSeq = alignSeqList[0]
        refSeqId = alignRefSeq[0]
        refSeqLabel = alignRefSeq[1]
        refSeqType = refSeqLabel.getSequenceType()
        alignRefSeqIdx = alignRefSeq[2]
        #conflictRefL   =alignRefSeq[3]
        refSeqFD = alignRefSeq[4]
        seqFeature.set(refSeqFD)
        lenRefAlignment = len(alignRefSeqIdx)
        #
        authDefinedMutations = self.__getAuthDefinedMutations(seqFeature.getEntityMutationDetailsOrig())

        # Changing the default state to false
        conflictRefL = [(0, '')] * len(alignRefSeqIdx)
        #
        refPartId, refSeqPartType = seqFeature.getPartInfo()
        self.__buildAuthPartDetails(refSeqId)
        numParts = len(self.__authPartIdList)

        if (self.__verbose):
            self.__lfh.write("\n\n+AlignmentUtils.reAnnotateAlignment() reannotating alignment with refId %s refPartId %r of numParts %d\n" % (refSeqId, refPartId, numParts))
        #
        # Assign annotations to regions outside of sequence parts --
        if (refSeqType == 'auth'):
            self.__mapAuthPartDetails(alignRefSeqIdx)
            if (self.__verbose):
                self.__lfh.write("+AlignmentUtils.__reAnnotateAlignment() adding out-of-part annotations to the reference sequence %s\n" % refSeqId)
            self.__annotateAuthOutOfPartComment(alignRefSeqIdx, conflictRefL)
            if (False):
                # JDW clear annotations
                for testPartId in self.__authPartIdList:
                    testSeqBeg = self.__authPartAlignD[testPartId][1]
                    testSeqEnd = self.__authPartAlignD[testPartId][2]
                    for ii in range(testSeqBeg - 1, testSeqEnd):
                        rT = alignRefSeqIdx[ii]
                        if len(rT[5]) > 0:
                            if (self.__verbose):
                                self.__lfh.write("+AlignmentUtils.__reAnnotateAlignment() clear annotation at position %d  = %s\n" % (ii + 1, rT[5]))
                        alignRefSeqIdx[ii] = (rT[0], rT[1], rT[2], rT[3], rT[4], '')

        #  For each input test sequence in the alignment -
        for aSeq in alignSeqList[1:]:
            testSeqId = aSeq[0]
            testSeqLab = aSeq[1]
            alignTestSeqIdx = aSeq[2]
            conflictL = aSeq[3]
            testSeqFD = aSeq[4]
            testSeqType = testSeqLab.getSequenceType()
            # ###
            if (self.__verbose):
                self.__lfh.write("\n\n+AlignmentUtils.reAnnotateAlignment() for reference id %s test id %s alignment length %d\n" %
                                 (refSeqId, testSeqId, len(alignTestSeqIdx)))
            #
            if testSeqType in ['auth', 'ref']:
                seqFeature.set(testSeqFD)
                testPartId, testSeqPartType = seqFeature.getPartInfo()
                # JDW get current part boundaries from internal data structure -
                testSeqBeg = self.__authPartAlignD[testPartId][1]
                testSeqEnd = self.__authPartAlignD[testPartId][2]
                if (self.__verbose):
                    self.__lfh.write("+AlignmentUtils.reAnnotateAlignment() for testId %s testPartId %r testSeqBeg %d  testSeqEnd %d\n" %
                                     (testSeqId, testPartId, testSeqBeg, testSeqEnd))
                #
            if testSeqType in ['ref']:
                for ii in range(testSeqBeg - 1, testSeqEnd):
                    rT = alignRefSeqIdx[ii]
                    tT = alignTestSeqIdx[ii]
                    if ((rT[1] == tT[1]) and (len(rT[5]) > 0)):
                        if (self.__verbose):
                            self.__lfh.write("+AlignmentUtils.__reAnnotateAlignment() clear annotation at position %d  = %s\n" % (ii + 1, rT[5]))
                        alignRefSeqIdx[ii] = (rT[0], rT[1], rT[2], rT[3], rT[4], '')
                        alignTestSeqIdx[ii] = (tT[0], tT[1], tT[2], tT[3], tT[4], '')
            #
            lenTestAlignment = len(alignTestSeqIdx)
            if (self.__debug):
                self.__formatAlignedPair(self.__lfh, refSeqId, alignRefSeqIdx, testSeqId, alignTestSeqIdx)
            #

            conflictL = [(0, '')] * lenRefAlignment
            numConflicts = 0

            for ii in range(0, len(alignRefSeqIdx)):
                # get the index in the reference sequence
                refTupIdx = alignRefSeqIdx[ii]
                testTupIdx = alignTestSeqIdx[ii]

                refCompId = refTupIdx[1]
                testCompId = testTupIdx[1]

                if ((refSeqType in ['auth']) and (refCompId == self.__gapSymbol) and (testCompId != self.__gapSymbol)):
                    conflictL[ii] = (1, '')
                    conflictRefL[ii] = (1, '')
                    numConflicts += 1
                    continue
                #
                #   Skip dual gaps
                if ((refCompId == self.__gapSymbol) and (testCompId == self.__gapSymbol)):
                    continue

                # We can only use the sequence array index here,  as the residue index may have undefined values for deletions!
                # refIdx=refTupIdx[3]+1
                # JDW -- Now using part mapping to the alignment coordinate system
                refIdx = ii + 1

                # ------------------------------------------------------
                #       --  Conflict assignments performed here --
                #

                # Skip conflicts if we are outside of the sequence index range for our sequence part.
                if ((numParts > 1) and (((refSeqType in ['auth']) and (testSeqType in ['ref'])) and ((refIdx < testSeqBeg) or (refIdx > testSeqEnd)))):
                    continue

                isConflict, refConflict, testConflict = self.__assignConflict(refSeqType, refCompId, testSeqType, testCompId, ii, lenRefAlignment)
                if isConflict:
                    # assign conflict to proper sequence part
                    conflictL[ii] = testConflict
                    conflictRefL[ii] = refConflict
                    numConflicts += 1
                    if self.__verbose and testSeqType in ['auth', 'ref']:
                        self.__lfh.write("+AlignmentUtils.reAnnotateAlignment() conflict %d align position %d ref index %d between %s and %s (refconflict %r testconflict %r\n" %
                                         (numConflicts, ii, refIdx, refCompId, testCompId, refConflict, testConflict))
                        myInPart = not ((refIdx < testSeqBeg) or (refIdx > testSeqEnd))
                        self.__lfh.write("+AlignmentUtils.reAnnotateAlignment() conflict %d align position %d refIdx %r testSeqBeg %r testSeqEnd %r inpart %r\n" %
                                         (numConflicts, ii, refIdx, testSeqBeg, testSeqEnd, myInPart))
            if (self.__verbose):
                self.__lfh.write("+AlignmentUtils.reAnnotateAlignment() %s conflicts %d length alignment %d\n"
                                 % (testSeqId, numConflicts, len(alignTestSeqIdx)))

            if ((refSeqType in ['auth']) and (testSeqType in ['ref'])):
                #
                #refSeqPartBeg = self.__authPartD[testPartId][1]
                #refSeqPartEnd = self.__authPartD[testPartId][2]
                #
                # Need to get authSeqBeg & End for the current sequence part for testPartId -- in the alignment coordinate system!!
                #
                refSeqPartBeg = self.__authPartAlignD[testPartId][1]
                refSeqPartEnd = self.__authPartAlignD[testPartId][2]
                self.__annotateAlignmentWithDefaultConflicts(alignRefDbSeqIdx=alignTestSeqIdx, alignAuthSeqIdx=alignRefSeqIdx, authSeqBeg=refSeqPartBeg,
                                                             authSeqEnd=refSeqPartEnd, authSeqPartId=testPartId, authDefinedMutations=authDefinedMutations)

            aSeq[3] = conflictL

        alignRefSeq[3] = conflictRefL
        self.__clearExcludedConflicts(alignSeqList)
        return True

    ##
    def __getAlignedSequenceIndex(self, id, alignSeqList):
        for idx in range(0, len(alignSeqList)):
            if (id == alignSeqList[idx][0]):
                return idx
        return -1

    def __getAlignedSequence(self, id, alignSeqList):
        for idx in range(0, len(alignSeqList)):
            if (id == alignSeqList[idx][0]):
                return alignSeqList[idx][2]
        return []

    ##
    def updateSequences(self, edSeqIdList, alignIdList, alignSeqList):
        """ For the input list of sequences (edSeqIdList) propogate the alignment edits to new versions of the
            sequences in the sequence data store.

            On input:

            alignSeqList [[Sequence identifier (ie. auth_A_1_1_1 as used in sequence data store),
                           SequenceLabel() object,
                           Aligned sequence with index details,
                           Conflict flag list (bool),
                           Feature dictionary for input sequence],,,]


            Each aligned sequences with index details has the storage model -

            (one-letter-code, 3-letter-code, original label residue index, position in sequence, position in alignment, comment )

            Stored sequences in the sequence data store have the simpler storage model -

             (3-letter-code, original residue index, comment/details, alignment index )

            Returns-  a revised list of sequences to be aligned in the next iteration.
                      being updated with any sequence ids updated with new versions in this module.

            ****

            Annotations in alignment comment field <conflict annotations> are transfered to the comment
            field of the associated sequence in the sequence data store.

        """
        sLabel = SequenceLabel()
        sFeature = SequenceFeature()
        if self.__sds is None:
            self.__sds = SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        #
        modSeqIdList = []
        retAlignIdList = copy.deepcopy(alignIdList)
        #
        # Update the part details -
        #
        tId = alignIdList[0]
        sLabel.unpack(tId)
        (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = sLabel.get()
        if seqType == 'auth':
            if (self.__verbose):
                self.__lfh.write("\n\n+AlignmentUtils.updateSequences() build part dictionary for seqId = %r\n" % tId)
            self.__buildAuthPartDetails(tId)
            alSeqIdx = self.__getAlignedSequence(tId, alignSeqList)
            # COME BACK
            self.__updateAuthPartDetails(alSeqIdx)

        #
        for seqId in edSeqIdList:
            if (self.__verbose):
                self.__lfh.write("+AlignmentUtils.updateSequences() Starts by applying edits to  seqId = %r (min %r max %r)\n" %
                                 (seqId, self.__authSeqBounds[0], self.__authSeqBounds[1]))
            #
            idxAl = self.__getAlignedSequenceIndex(seqId, alignSeqList)
            if (idxAl < 0):
                continue
            #
            # Edited sequence from alignment -
            #
            alSeqIdx = self.__getAlignedSequence(seqId, alignSeqList)
            #
            # Recover the sequence and feature dictionary from the data store -
            #
            sLabel.unpack(seqId)
            (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = sLabel.get()
            seqIdx = self.__sds.getSequence(seqId=seqInstId, seqType=seqType, partId=seqPartId, altId=seqAltId, version=seqVersion)
            #
            #
            seqVersionNext = seqVersion + 1
            seqIdxNext = []
            #
            if seqType in ['auth']:
                # **** JDW  Discard gaps in the auth sequence
                # **** NOTE that prior comments are copied over  --
                minAl = -1
                maxAl = -1
                for ii, aTup in enumerate(alSeqIdx):
                    # JDW JDW 2014-Oct-21
                    if aTup[1] == self.__gapSymbol:
                        continue
                    if ((len(aTup[2]) > 0) and (int(aTup[2]) == self.__authSeqBounds[0])):
                        minAl = ii
                    if ((len(aTup[2]) > 0) and (int(aTup[2]) == self.__authSeqBounds[1])):
                        maxAl = ii
                    # if self.__verbose:
                    #    self.__lfh.write("+AlignmentUtils.updateSequences() ii %d %r  aTup[2] %r minAl %d maxAl %d\n" % (ii,aTup[0],aTup[2],minAl,maxAl))
                if (self.__verbose):
                    self.__lfh.write("+AlignmentUtils.updateSequences() seq %s alignment index min %r max %r\n" % (seqId, minAl, maxAl))

                for ii, aTup in enumerate(alSeqIdx):
                    if ((aTup[1] != self.__gapSymbol) and (ii >= minAl) and (ii <= maxAl)):
                        seqIdxNext.append((aTup[1], aTup[2], aTup[5], aTup[4]))
            elif seqType in ['xyz']:
                # **** Preserve gaps in coordinate sequences -
                #      comments are preserved from prior sequence --
                for aTup in alSeqIdx:
                    seqIdxNext.append((aTup[1], aTup[2], aTup[5], aTup[4]))
            else:
                #
                # **** Gaps are filtered from reference sequences -
                for aTup in alSeqIdx:
                    if aTup[1] != self.__gapSymbol:
                        seqIdxNext.append((aTup[1], aTup[2], aTup[5], aTup[4]))

            if self.__verbose:
                self.__lfh.write("+AlignmentUtils.updateSequences() updated sequence %s length %d new version %s\n"
                                 % (seqId, len(seqIdxNext), seqVersionNext))
                # for atup in seqIdxNext:
                #    self.__lfh.write("      +++__updateSequences() %s %s %s %s \n" % ((aTup[1],aTup[2],aTup[5],aTup[4])))

            #  Add new sequence version to data store -
            #
            partIdList = self.__sds.getPartIds(seqInstId, seqType=seqType)
            if (self.__verbose):
                self.__lfh.write("+AlignmentUtils.updateSequences() part id list %r\n" % partIdList)
                try:
                    self.__lfh.write("+AlignmentUtils.updateSequences() authPartD    %r\n" % self.__authPartD.items())
                except:
                    self.__lfh.write("+AlignmentUtils.updateSequences() authPartD    %r\n" % self.__authPartD)

            if seqType == "auth":
                idx = 1
                sPos = 1
                tSeqIdxNext = []
                for sTup in seqIdxNext:
                    tSeqIdxNext.append((sTup[0], str(sPos), sTup[2], idx))
                    if (sTup[0] != self.__gapSymbol):
                        sPos += 1
                    idx += 1

                if (self.__verbose):
                    self.__lfh.write("+AlignmentUtils.updateSequences() auth sequence length %d  idx-1 %d orig length %d\n" %
                                     (len(seqIdxNext), idx - 1, len(seqIdx)))
                #
                # All of the parts of author sequence are complete and common -
                #
                for pId in partIdList:
                    pTup = self.__authPartD[pId]
                    if (self.__verbose):
                        self.__lfh.write("+AlignmentUtils.updateSequences() store %s sequence for part id  %r prior length %d current length %d part min %d part max %d\n" %
                                         (seqType, pId, len(seqIdx), len(tSeqIdxNext), pTup[1], pTup[2]))

                    if ((pId == partIdList[-1]) and (len(tSeqIdxNext) != pTup[2])):
                        # if last part -
                        self.__authPartD[pId] = (pTup[0], pTup[1], len(tSeqIdxNext), pTup[3])
                        if (self.__verbose):
                            self.__lfh.write("+AlignmentUtils.updateSequences() part id %r sequence length extended to %r\n" % (pId, len(tSeqIdxNext)))

                    self.__sds.setSequence(tSeqIdxNext, seqId=seqInstId, seqType=seqType, partId=pId, altId=seqAltId, version=seqVersionNext)
            else:
                self.__sds.setSequence(seqIdxNext, seqId=seqInstId, seqType=seqType, partId=seqPartId, altId=seqAltId, version=seqVersionNext)

            #
            # Update the feature dictionary -
            #

            #
            if seqType == "auth":
                for pId in partIdList:
                    # Each author sequence part has a different feature dictionary -
                    if (self.__verbose):
                        self.__lfh.write("+AlignmentUtils.updateSequences() update stored features for instance/entity %s part id  %r\n" % (seqInstId, pId))
                    seqFD = self.__sds.getFeature(seqId=seqInstId, seqType=seqType, partId=pId, altId=seqAltId, version=seqVersion)
                    sFeature.set(seqFD)
                    sFeature.clearAlignDetails()
                    #
                    #  JDW JDW Sep 10, 2014  -- Revised to update all sequence parts
                    # if ((pId == partIdList[-1]) and (len(seqIdxNext) !=  len(seqIdx))):
                    if (True):
                        #  JDW Use the internal updated part boundaries --
                        # jpartId,jSeqPartBeg,jSeqPartEnd,jSeqPartType=sFeature.getAuthPartDetails()
                        #
                        jPartId, jSeqPartType = sFeature.getPartInfo()
                        jSeqPartBeg = self.__authPartD[pId][1]
                        jSeqPartEnd = self.__authPartD[pId][2]
                        #sFeature.setAuthPartDetails(pId,jSeqPartBeg, len(seqIdxNext), jSeqPartType)
                        #
                        sFeature.setAuthPartDetails(pId, jSeqPartBeg, jSeqPartEnd, jSeqPartType)

                        # Update related reference features with the revised length
                        altIdList = self.__sds.getAlternativeIds(seqInstId, dataType="sequence", seqType="ref", partId=pId)
                        for altId in altIdList:
                            fObj = self.__sds.getFeatureObj(seqId=seqInstId, seqType='ref', partId=pId, altId=altId, version=seqVersion)
                            # kPartId,kSeqPartBeg,kSeqPartEnd,kSeqPartType=fObj.getAuthPartDetails()
                            #fObj.setAuthPartDetails(kPartId,kSeqPartBeg, len(seqIdxNext), kSeqPartType)
                            fObj.setAuthPartDetails(pId, jSeqPartBeg, jSeqPartEnd, jSeqPartType)
                            self.__sds.setFeature(fObj.get(), seqId=seqInstId, seqType='ref', partId=pId, altId=altId, version=seqVersion)
                    if (self.__debug):
                        sFeature.printIt(log=self.__lfh)
                    self.__sds.setFeature(sFeature.get(), seqId=seqInstId, seqType=seqType, partId=pId, altId=seqAltId, version=seqVersionNext)
            else:
                seqFD = self.__sds.getFeature(seqId=seqInstId, seqType=seqType, partId=seqPartId, altId=seqAltId, version=seqVersion)
                sFeature.set(seqFD)
                sFeature.clearAlignDetails()
                self.__sds.setFeature(sFeature.get(), seqId=seqInstId, seqType=seqType, partId=seqPartId, altId=seqAltId, version=seqVersionNext)
            #
            # Create new version of the sequence Id and exchange this on the alignId list.
            #
            sLabel.set(seqType=seqType, seqInstId=seqInstId, seqPartId=seqPartId, seqAltId=seqAltId, seqVersion=seqVersionNext)
            seqIdNext = sLabel.pack()
            retAlignIdList.remove(seqId)
            retAlignIdList.append(seqIdNext)
            modSeqIdList.append((seqId, seqIdNext))

        self.__sds.serialize()
        #
        return retAlignIdList, modSeqIdList
    ##

    def updateAlignment(self, alignSeqList, entityId=''):
        """  Apply any stored edits to the current alignment in the input sequence alignment list.

             Update the aligned sequences stored as -

             alignSeqList [[Sequence Id,  SequenceLabel() object, aligned sequence w/ index, conflict list, sequence feature dictionary],,,]

             Returns - edIdList[] list id's for edited sequences -

             Note -- Currently removing edits after application -

        """
        #
        realignFlag = False
        edIdList = []
        #
        if (self.__verbose):
            self.__lfh.write("+AlignmentUtils.updateAlignment() operation %s starting\n" % self.__operation)
        #
        esfn = getEditStoreFilename(entityId)
        ses = SequenceEditStore(sessionObj=self.__sessionObj, fileName=esfn, verbose=self.__verbose, log=self.__lfh)
        edObjList = ses.getList()
        if (self.__verbose):
            self.__lfh.write("+AlignmentUtils.updateAlignment() edit store %s edit list length %d\n" % (esfn, len(edObjList)))
        #
        if (len(edObjList) < 1):
            self.__normalizeAlignmentLength(alignSeqList)
            return edIdList, realignFlag
        #
        resLabel = ResidueLabel()
        #
        # Apply edits to aligned sequences in memory
        #
        # store applied edits for deletion -
        deleteObjList = []
        #
        # Find the stored edits associated with each of the aligned sequences and apply these,
        # maintain a list of applied edits so that these may be deleted from the edit store
        # after application.
        #
        for aSeq in alignSeqList:
            id = aSeq[0]
            aLabel = aSeq[1]
            aLabel.unpack(id)
            (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = aLabel.get()
            seqWithIdx = aSeq[2]

            if (self.__verbose):
                self.__lfh.write("+AlignmentUtils.__updateAlignment() searching edits for : seqType %s InstId %s partId %d altId %d seqVersion %d \n"
                                 % (seqType, seqInstId, seqPartId, seqAltId, seqVersion))
            #
            # Make a list of edit objects for the current aligned sequence.
            #
            eL = []
            for edObj in edObjList:
                tId = edObj.getTargetElementId()
                resLabel.unpack(tId)

                (seqTypeEd, seqInstIdEd, seqPartIdEd, seqAltIdEd, seqVersionEd) = resLabel.getSeq()
                #
                # self.__lfh.write("+AlignmentUtils.__updateAlignment() check edit for seqType %s InstI %s AltId %d seqVersion %d \n"
                #                 %  (seqTypeEd,seqInstIdEd,int(seqAltIdEd),int(seqVersionEd)))

                if ((seqType == seqTypeEd) and
                        (seqInstId == seqInstIdEd) and
                        (int(seqPartId) == int(seqPartIdEd)) and
                        (int(seqAltId) == int(seqAltIdEd)) and
                        (int(seqVersion) == int(seqVersionEd))):
                    if (self.__verbose):
                        self.__lfh.write("+AlignmentUtils.__updateAlignment() matched this edit for sequence tId %s\n" % tId)
                    eL.append(edObj)
            #
            if (self.__verbose):
                self.__lfh.write("+AlignmentUtils.__updateAlignment() matched %d edits for sequence id %s\n" % (len(eL), id))
                self.__lfh.flush()

            #
            # Apply edits
            #
            if len(eL) > 0:
                if (self.__verbose):
                    self.__lfh.write("+AlignmentUtils.__updateAlignment() for sequence id %s aligned sequence length %d the edit list length is %d\n" %
                                     (id, len(seqWithIdx), len(eL)))
                    self.__lfh.flush()
                aSeq[2], mFlag = self.__applyEditsToAlignments(eL, seqWithIdx)
                if mFlag:
                    realignFlag = True
                a3L = []
                for sPos in aSeq[2]:
                    a3L.append(sPos[1])
                if (self.__verbose):
                    self.__lfh.write("+AlignmentUtils.__updateAlignment() after applying edits on sequence id %s the test sequence length is %d\n" % (id, len(a3L)))
                #
                # Renumber - update auth index for any insertions --  entity_poly_seq.num and sequence position index must be unique and exist
                #                                                     for any non gap residues --  Care here for data type --
                #
                if seqType in ['auth']:
                    ii = 1
                    for iPos in range(0, len(aSeq[2])):
                        if aSeq[2][iPos][0] != self.__gapSymbol:
                            aSeq[2][iPos] = (aSeq[2][iPos][0], aSeq[2][iPos][1], str(ii), ii - 1, aSeq[2][iPos][4], aSeq[2][iPos][5])
                            ii += 1

                #
                deleteObjList.extend(eL)

                edIdList.append(id)

            self.__lfh.flush()

        ses.deleteEditList(deleteObjList)
        #
        self.__normalizeAlignmentLength(alignSeqList)
        #
        return edIdList, realignFlag

    def __normalizeAlignmentLength(self, alignSeqList):
        """
        Extend the aligned length of each sequence to correspond to the maximun aligned sequence length (after edits/inserts):

            alignSeqList [[Sequence identifier (ie. auth_A_1_1_1 as used in sequence data store),
                           SequenceLabel() object,
                           Aligned sequence with index details,
                           Conflict flag list (bool),
                           Feature dictionary for input sequence],,,]

        seqWithIndex[]     aligned sequence as a list of tuples of  -
                            (one-letter-code, 3-letter-code, original label residue index, position in sequence, position in alignment, comment )

        """
        maxL = -1
        for aSeq in alignSeqList:
            iLen = len(aSeq[2])
            maxL = max(maxL, iLen)
        #
        if (self.__verbose):
            self.__lfh.write("+AlignmentUtils.__normalizeAlignmentLength() maximum sequence length is %r\n" % maxL)
        for aSeq in alignSeqList:
            if len(aSeq[2]) < maxL:
                for iPos in range(len(aSeq[2]), maxL):
                    sTup = ('.', '.', '', '', iPos, '')
                    aSeq[2].append(sTup)
                for iPos in range(len(aSeq[2]), maxL):
                    aSeq[3].append((0, ''))
        #
        # Filter out empty slots accross the whole alignment arising from deletions ---
        #
        okL = [False] * maxL
        for iPos in range(0, maxL):
            for aSeq in alignSeqList:
                sTup = aSeq[2][iPos]
                if sTup[0] != self.__gapSymbol:
                    okL[iPos] = True
                    break
        #
        #
        iCount = okL.count(False)
        if (self.__verbose):
            self.__lfh.write("+AlignmentUtils.__normalizeAlignmentLength() removing %d empty records\n" % iCount)
        #
        if iCount > 0:
            for aSeq in alignSeqList:
                aSeqNew = list(itertools.compress(aSeq[2], okL))
                aSeq[2] = []
                for ii, sTup in enumerate(aSeqNew):
                    t = list(sTup)
                    t[4] = ii
                    aSeq[2].append(tuple(t))

    def __applyEditsToAlignments(self, editList, seqWithIndexList):
        """
        Apply edits to the input aligned sequence  - return the edited indexed sequence -

        editList[]          list of edit objects to be applied to indexed input sequence

        seqWithIndexList[]  aligned sequence as a list of tuples of  -
                            (one-letter-code, 3-letter-code, original label residue index, position in sequence, position in alignment, comment )

        Return edited sequence as a list of tuples containing -

           (one-letter-code, 3-letter-code, original label residue index, position in sequence, position in alignment, comment)
        """
        realignFlag = False
        resLabel = ResidueLabel()
        #
        # create a list of deletions
        dList = []
        for edObj in editList:
            tId = edObj.getTargetElementId()
            resLabel.unpack(tId)
            iPos = resLabel.getAlignmentIndex()
            if (edObj.getEditType() == 'delete'):
                dList.append(iPos)
        #
        if (self.__verbose):
            self.__lfh.write("+AlignmentUtils.__applyEditsToAlignments() delete list is %r\n" % dList)
        if len(dList) > 0:
            realignFlag = True

        for edObj in editList:
            tId = edObj.getTargetElementId()
            edType = edObj.getEditType()
            edValue = edObj.getValueNew()
            resLabel.unpack(tId)
            try:
                iPos = int(resLabel.getAlignmentIndex())
            except:
                continue

            if (self.__verbose):
                self.__lfh.write("+AlignmentUtils.__applyEditsToAlignments() operation %s at alignment index %d with %r\n" %
                                 (edType, iPos, edValue))

            if (edType == 'replace' or edType == 'insert'):
                edv = " ".join(edValue)
                seqWithIndexList[iPos] = (seqWithIndexList[iPos][0], edv, seqWithIndexList[iPos][2],
                                          seqWithIndexList[iPos][3], iPos, seqWithIndexList[iPos][5])
            elif (edType == 'replaceid'):
                try:
                    newId = edObj.getNewElementId()
                    dstLab = ResidueLabel()
                    dstLab.unpack(newId)
                    dstAlignPos = dstLab.getAlignmentIndex()
                    # make this an integer always --  JDW THIS IS CREATING PROBLEMS IN HAVE INSERTION CODES MAY 29,2014
                    #                    dstSeqPos  =int(dstLab.getSequenceIndex())
                    dstSeqPos = dstLab.getSequenceIndex()
                    dstLblInd = dstLab.getResidueLabelIndex()
                    dstResType = dstLab.getResidueType()
                    edv = " ".join(edValue)
                    oneLc = self.__srd.cnv3To1(edv)
                    if (self.__verbose):
                        self.__lfh.write("+AlignmentUtils.__applyEditsToAlignments() newId %s lbl index %s seq position %s edv %s olc %s\n" %
                                         (newId, dstLblInd, dstSeqPos, edv, oneLc))
                        self.__lfh.flush()
                    # Don't propagate comments on gap positions or deletions
                    if ((oneLc != self.__gapSymbol) and (seqWithIndexList[iPos][5] not in ['deletion'])):
                        seqWithIndexList[iPos] = (oneLc, edv, dstLblInd, dstSeqPos, iPos, seqWithIndexList[iPos][5])
                    else:
                        seqWithIndexList[iPos] = (oneLc, edv, dstLblInd, dstSeqPos, iPos, '')
                except:
                    traceback.print_exc(file=self.__lfh)
                    self.__lfh.flush()

            elif (edType == 'details'):
                seqWithIndexList[iPos] = (seqWithIndexList[iPos][0], seqWithIndexList[iPos][1], seqWithIndexList[iPos][2],
                                          seqWithIndexList[iPos][3], iPos, edValue)

        #
        self.__updateMissingIndices(seqWithIndexList)

        # Build the edited sequence in seqT[]
        seqT = []
        insC = "ZABCDEFGHIJKLMNOPQRSTUVWXY"
        for idx in range(0, len(seqWithIndexList)):
            # Skip over alignment positions marked for deletion -
            #
            if idx not in dList:
                if (len(seqWithIndexList[idx][1].strip().split()) > 1):
                    realignFlag = True
                    # handle insertions  -
                    tL = seqWithIndexList[idx][1].strip().split()
                    seqT.append((self.__srd.cnv3To1(tL[0]), tL[0], seqWithIndexList[idx][2], seqWithIndexList[idx][3], idx, seqWithIndexList[idx][5]))
                    for it in range(1, len(tL)):
                        jj = it % len(insC)
                        seqT.append((self.__srd.cnv3To1(tL[it]), tL[it], str(seqWithIndexList[idx][2]) + insC[jj], seqWithIndexList[idx][3], idx, seqWithIndexList[idx][5]))
                else:
                    seqT.append((self.__srd.cnv3To1(seqWithIndexList[idx][1]), seqWithIndexList[idx][1],
                                 seqWithIndexList[idx][2], seqWithIndexList[idx][3], idx, seqWithIndexList[idx][5]))
        #
        # for tup in seqT:
        #    self.__lfh.write("+AlignmentUtils.__applyEdits() - ed %s - index %s\n" % (tup[0],tup[1]))
        #
        if (self.__verbose):
            self.__lfh.write("+AlignmentUtils.__applyEditsToAlignments() returns sequence length %d realignFlag %r\n" % (len(seqT), realignFlag))

        return seqT, realignFlag

    def __updateMissingIndices(self, alSeq):
        """  Scan the aligned sequence for missing index values making reasonable substitutions
             based on surronding context.

            Aligned sequences have storage model -

            (one-letter-code, 3-letter-code, original label residue index, position in sequence, position in alignment, comment )
        """
        # scan in reverse
        alSeq.reverse()
        for ii in range(1, len(alSeq)):
            prevTup = alSeq[ii - 1]
            curTup = alSeq[ii]
            if ((curTup[1] != self.__gapSymbol) and (not curTup[2]) and prevTup[2]):
                v = int(str(prevTup[2])) - 1
                alSeq[ii] = (curTup[0], curTup[1], str(v), curTup[3], curTup[4], curTup[5])
                if self.__debug:
                    self.__lfh.write("   ++++++++++  + __updateMissingIndices() udpdating missing index %r %r \n" % (curTup[1], v))
            else:
                # if self.__debug:
                #    self.__lfh.write("    + __updateMissingIndices() ii %d curTup[1] %r curTup[2] %r prevTup[1] %r prevTup[2] %r not gap %r is missing %r\n" %
                #                     (ii,curTup[1],curTup[2],prevTup[1],prevTup[2], curTup[1] != self.__gapSymbol, not curTup[2]) )
                pass

        # scan forward -
        alSeq.reverse()
        for ii in range(1, len(alSeq)):
            prevTup = alSeq[ii - 1]
            curTup = alSeq[ii]
            if ((curTup[1] != self.__gapSymbol) and (not curTup[2]) and prevTup[2]):
                v = int(str(prevTup[2])) + 1
                alSeq[ii] = (curTup[0], curTup[1], str(v), curTup[3], curTup[4], curTup[5])
                if self.__debug:
                    self.__lfh.write("   ++++++++++  + __updateMissingIndices() udpdating missing index %r %r \n" % (curTup[1], v))

    def __clearExcludedConflicts(self, alignSeqList):
        """ For the input aligned sequences clear conflicts and comment details in excluded entity parts.

            On input:

            alignSeqList [[Sequence identifier (ie. auth_A_1_1_1 as used in sequence data store),
                           SequenceLabel() object,
                           Aligned sequence with index details,
                           Conflict flag list (bool),
                           Feature dictionary for input sequence],,,]


            Each aligned sequences with index details has the storage model -

            (one-letter-code, 3-letter-code, original label residue index, position in sequence, position in alignment, comment )

        """
        if ((len(self.__excludedPartIdList) < 1) or (alignSeqList is None) or (len(alignSeqList) < 1)):
            return False

        sLabel = SequenceLabel()

        alignRef = True
        for alignTup in alignSeqList:
            seqId = alignTup[0]
            alSeqIdx = alignTup[2]
            cList = alignTup[3]
            sLabel.unpack(seqId)
            (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = sLabel.get()
            if (self.__verbose):
                self.__lfh.write("+AlignmentUtils.__clearExcludedConflicts() clearing updating comments for seqId = %r length %d\n" % (seqId, len(alSeqIdx)))
            #
            if ((seqType in ['auth']) and alignRef):
                for pId in self.__excludedPartIdList:
                    refSeqPartBeg = self.__authPartAlignD[pId][1]
                    refSeqPartEnd = self.__authPartAlignD[pId][2]

                    for aPos in range(0, len(alSeqIdx)):
                        aTup = alSeqIdx[aPos]
                        # refIdx=aTup[3]+1
                        refIdx = aPos + 1
                        if ((refIdx < refSeqPartBeg) or (refIdx > refSeqPartEnd)):
                            continue
                        cList[aPos] = (0, '')
            #
            alignRef = False

        return True

    def __XXXgetAuthSeqOutOfPartComment(self, idx):
        """      Return annotation assignment for residue position 'idx' based on the
                 the input part definition dictionary for out-of-part  terminal expression
                 tags or internal linkers.

                 input    idx   sequence position in author/sample sequence (1-N)

                 returns --   either '', 'expression tag', or 'linker'.
        """
        if self.__authSeqBounds is None:
            return ''
        for partId in self.__authPartIdList:
            if idx >= self.__authPartD[partId][1] or idx <= self.__authPartD[partId][2]:
                # return blank if in part --
                return ''

        # not in part --
        if idx < self.__authSeqBounds[0] or idx > self.__authSeqBounds[1]:
            return 'expression tag'
        else:
            return 'linker'

    def __getAuthDefinedMutations(self, entityMutationDetails):
        """ Extract all author input mutation information (e.g. V178A ) from _entity.pdbx_mutation
        """
        authDefinedMutations = []
        if not entityMutationDetails:
            return authDefinedMutations
        #
        details = entityMutationDetails.upper().replace(',', ' ')
        detailList = details.split(' ')
        for val in detailList:
            if len(val) < 3:
                continue
            #
            hasDigit = False
            allDigit = True
            for i in range(1, len(val) - 1):
                if val[i].isdigit():
                    hasDigit = True
                else:
                    allDigit = False
                #
            #
            if hasDigit and allDigit and val[0].isalpha() and val[len(val) - 1].isalpha():
                authDefinedMutations.append(val)
                self.__lfh.write("+AlignmentUtils.__getAuthDefinedMutations insert %s\n" % val)
            else:
                self.__lfh.write("+AlignmentUtils.__getAuthDefinedMutations skip %s\n" % val)
            #
        #
        return authDefinedMutations
