##
# File:  AlignmentStatistics.py
# Date:  12-Feb-2010
#
# Updates:
# 20-Apr-2010 jdw Ported to module seqmodule.
# 10-Aug-2010 RPS: doUpdate() now logging via self.__lfh.write() instead of sys.stderr.write()
# 18-Dec-2012 jdw  align full list of candidates ---
# 04-Mar-2013 jdw  Overhaul to support the entity polymer parts.
#  4-Dec-2013 jdw  Add linkage details for auth/xyz alignments
##
"""
Update alignment statistics for sequences in the sequence data store.


Methods are provided to compute pairwise alignment statistics between current
latest version of the author sequences and the latest version of each coordinate
and reference sequences.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.07"


import sys
from wwpdb.apps.seqmodule.io.SequenceDataStore import SequenceDataStore
from wwpdb.apps.seqmodule.util.SequenceReferenceData import SequenceReferenceData
from wwpdb.apps.seqmodule.util.SequenceLabel import SequenceFeature

#
from wwpdb.utils.align.alignlib import PairwiseAlign  # pylint: disable=no-name-in-module


class AlignmentStatistics(object):
    """
    Update alignment statistics for sequences in the sequence data store.

     Methods are provided to compute pairwise alignment stats between current
    latest version of the author sequence and the latest version of each coordinate
    and reference sequence.

    """

    def __init__(self, reqObj=None, maxRefAlign=100, verbose=False, log=sys.stderr):
        self.__reqObj = reqObj
        self.__verbose = verbose
        self.__debug = False
        self.__lfh = log
        self.__maxRefAlign = maxRefAlign
        #
        self.__sessionObj = self.__reqObj.getSessionObj()
        self.__sessionPath = self.__sessionObj.getPath()

        #
        self.__sessionPath = "."
        #
        self.__srd = SequenceReferenceData(verbose=self.__verbose, log=self.__lfh)
        self.__gapSymbol = self.__srd.getGapSymbol()
        #
        self.__ok = True
        self.__setup()

    def __setup(self):
        """ """
        try:
            self.__sessionPath = self.__sessionObj.getPath()
            self.__sds = SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)

        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+AlignmentStatistics.__setup() - using session path %s failing \n" % self.__sessionPath)
            self.__ok = False

    #
    def doUpdate(self):
        if not self.__ok:
            if self.__verbose:
                self.__lfh.write("+AlignmentStatistics.doUpdate()  sessionId %s failed\n" % self.__sessionObj.getId())
            return False

        if self.__verbose:
            self.__lfh.write("\n+AlignmentStatistics.doUpdate()  STARTING for sessionId %s\n" % (self.__sessionObj.getId()))

        seqFeature = SequenceFeature(self.__verbose)
        rD = self.__updateAuthXyzAlignments()
        for kTup, vTup in rD.items():
            fD = self.__sds.getFeature(seqId=kTup[0], seqType=kTup[1], partId=kTup[2], altId=kTup[3], version=kTup[4])
            seqFeature.clear()
            seqFeature.set(fD)
            sim = float(vTup[2]) / float(vTup[1])
            simG = float(vTup[3]) / float(vTup[1])
            seqFeature.setAuthXyzAlignDetails(seqLen=vTup[0], alignLen=vTup[1], seqSim=sim, seqSimWithGaps=simG)
            seqFeature.setAuthXyzAlignRange(seqBegin=vTup[4], seqEnd=vTup[5])
            self.__sds.setFeature(seqFeature.get(), seqId=kTup[0], seqType=kTup[1], partId=kTup[2], altId=kTup[3], version=kTup[4])

        rD = self.__updateAuthRefAlignments()
        for kTup, vTup in rD.items():
            fD = self.__sds.getFeature(seqId=kTup[0], seqType=kTup[1], partId=kTup[2], altId=kTup[3], version=kTup[4])
            seqFeature.clear()
            seqFeature.set(fD)
            sim = float(vTup[2]) / float(vTup[1])
            simG = float(vTup[3]) / float(vTup[1])
            seqFeature.setAuthRefAlignDetails(seqLen=vTup[0], alignLen=vTup[1], seqSim=sim, seqSimWithGaps=simG)
            self.__sds.setFeature(seqFeature.get(), seqId=kTup[0], seqType=kTup[1], partId=kTup[2], altId=kTup[3], version=kTup[4])

        self.__sds.serialize()

        if self.__verbose:
            self.__lfh.write("+AlignmentStatistics.doUpdate()  COMPLETED\n\n")
        return True

    #
    def __updateAuthXyzAlignments(self):
        """Update alignment statistics for author/entity sequence and  each coordinate sequence.

        Author and coordinate sequences are full length matches independent of any sequence patitioning.
        """
        #
        if self.__verbose:
            self.__lfh.write("\n+AlignmentStatistics.__updateAuthXyzAlignments()  sessionId %s\n" % self.__sessionObj.getId())
        #
        # For each entity sequence  - compute pairwise alignment with the instance chains.
        #
        refSeqLen = 0
        altId = 1
        partId = 1
        resultD = {}
        pA = PairwiseAlign()
        pA.setVerbose(self.__verbose)
        #
        # Get the alignment reference sequence (e.g. auth sequence with highest version number)
        #
        groupIdList = self.__sds.getGroupIds()
        if self.__verbose:
            self.__lfh.write("+AlignmentStatistics.__updateAuthXyzAlignments() group list %r\n" % groupIdList)

        for gId in groupIdList:
            seqIdList = self.__sds.getGroup(gId)
            if self.__verbose:
                self.__lfh.write("+AlignmentStatistics.__updateAuthXyzAlignments() sequence instances in group %s: %r\n" % (gId, seqIdList))
            #
            # JDW CHANGE
            if len(seqIdList) < 1:
                if self.__verbose:
                    self.__lfh.write("+AlignmentStatistics.__updateAuthXyzAlignments()  group %s is empty\n" % gId)
                continue
            #
            # seqId0=seqIdList[0]
            seqId0 = gId
            #
            # Show the author sequence(s) for only the first sequence in the group (all versions)
            #
            verList = self.__sds.getVersionIds(seqId0, partId=partId, altId=altId, dataType="sequence", seqType="auth")
            if len(verList) < 1:
                if self.__verbose:
                    self.__lfh.write("+AlignmentStatistics.__updateAuthXyzAlignments()  group %s no auth sequences for id %s\n" % (gId, seqId0))
                continue
            #
            #
            verLatest = verList[0]
            seqAuthIdx = self.__sds.getSequence(seqId=seqId0, seqType="auth", partId=partId, altId=altId, version=verLatest)
            #
            if self.__verbose:
                self.__lfh.write(
                    "+AlignmentStatistics.__updateAuthXyzAlignments() auth sequence entity group %s instance %s version %s length %d\n" % (gId, seqId0, verLatest, len(seqAuthIdx))
                )
            #
            refSeqLen = len(seqAuthIdx)
            #
            r3L = []
            for tup in seqAuthIdx:
                r3L.append(str(tup[0]))
            pA.setReferenceSequence(r3L, "auth:" + str(seqId0))

            # -----
            # Get the instance sequences from the coordinate data.
            #
            # Stored sequences in the sequence data store have the storage model -
            #
            # (3-letter-code, original residue index, comment/details, alignment index )
            #
            sLen = {}
            sRange = {}
            seqIdExistList = []
            for seqId in seqIdList:
                verList = self.__sds.getVersionIds(seqId=seqId, partId=partId, altId=1, dataType="sequence", seqType="xyz")
                if len(verList) < 1:
                    if self.__verbose:
                        self.__lfh.write("+AlignmentStatistics.__updateAuthXyzAlignments()  group %s no xyz sequences for id %s\n" % (gId, seqId))
                    continue
                seqIdExistList.append(seqId)
                verLatest = verList[0]
                seqXyzIdx = self.__sds.getSequence(seqId, "xyz", partId=partId, altId=1, version=verLatest)

                if self.__verbose:
                    self.__lfh.write(
                        "+AlignmentStatistics.__updateAuthXyzAlignments() xyz  sequence entity group %s instance %s version %s length %d\n"
                        % (gId, seqId, verLatest, len(seqXyzIdx))
                    )
                sBeg = seqXyzIdx[0][1]
                sEnd = seqXyzIdx[len(seqXyzIdx) - 1][1]
                sRange[seqId] = (sBeg, sEnd)
                r3L = []
                typicalLink = []
                for tup in seqXyzIdx:
                    # Filter any gaps prior to alignment  --
                    if tup[0] != self.__gapSymbol:
                        r3L.append(str(tup[0]))
                        typicalLink.append(0 if ("long_begin" in tup[2]) else 1)

                sLen[seqId] = len(r3L)
                pA.addTestSequenceWithLink(r3L, "xyz:" + str(seqId), typicalLink)
            #
            pA.doAlign()
            #

            #
            for seqId in seqIdList:
                if seqId not in seqIdExistList:
                    continue
                verList = self.__sds.getVersionIds(seqId=seqId, partId=partId, altId=altId, dataType="sequence", seqType="xyz")
                verLatest = verList[0]
                #
                aL = pA.getAlignment("xyz:" + str(seqId))
                alignLength = len(aL)
                numMatch = 0
                numMatchGaps = 0
                sBeg, sEnd = sRange[seqId]
                for aTup in aL:
                    if aTup[0] == aTup[1]:
                        numMatch += 1
                    if aTup[0] == aTup[1] or aTup[1] == self.__gapSymbol:
                        numMatchGaps += 1

                resultD[(seqId, "xyz", partId, altId, verLatest)] = (sLen[seqId], alignLength, numMatch, numMatchGaps, sBeg, sEnd)
        #
        if self.__verbose:
            self.__lfh.write("\n+AlignmentStatistics.__updateAuthXyzfAlignments()  Alignment summary for author (length=%d) and coordinate sequences\n" % refSeqLen)
        if self.__debug:
            for kTup, v in resultD.items():
                self.__lfh.write(
                    "  -  xyz seq  partId %d alt %d ver %d seqId %s seq length %d align length %d matches %d matches w/ gaps %d sBeg %s sEnd %s\n"
                    % (kTup[2], kTup[3], kTup[4], kTup[0], v[0], v[1], v[2], v[3], v[4], v[5])
                )
        #

        return resultD

    # ##
    def __updateAuthRefAlignments(self):
        """Update alignment statistics for author/entity sequence and each reference sequence."""
        #
        if self.__verbose:
            self.__lfh.write("\n+AlignmentStatistics.__updateAuthRefAlignments()  sessionId %s\n" % (self.__sessionObj.getId()))
        #
        # For each entity sequence  - compute pairwise alignment with the each reference sequence.
        #
        resultD = {}
        pA = PairwiseAlign()
        pA.setVerbose(self.__verbose)
        #
        seqFeature = SequenceFeature(self.__verbose)
        #
        # Get the alignment reference sequence (e.g. auth sequence with highest version number)
        #
        groupIdList = self.__sds.getGroupIds()
        for gId in groupIdList:
            seqIdList = self.__sds.getGroup(gId)

            if len(seqIdList) < 1:
                if self.__verbose:
                    self.__lfh.write("+AlignmentStatistics.__updateAuthRefAlignments() group %s is empty\n" % gId)
                # continue
            # JDW CHANGE
            # seqId0=seqIdList[0]
            seqId0 = gId
            partIdList = self.__sds.getPartIds(seqId0, dataType="sequence", seqType="auth")
            for partId in partIdList:

                verList = self.__sds.getVersionIds(seqId0, partId=partId, altId=1, dataType="sequence", seqType="auth")
                if len(verList) < 1:
                    if self.__verbose:
                        self.__lfh.write("+AlignmentStatistics.__updateAuthXyzAlignments()  group %s no auth sequences for id %s\n" % (gId, seqId0))
                    continue

                verLatest = verList[0]
                seqAuthIdx = self.__sds.getSequence(seqId=seqId0, seqType="auth", partId=partId, altId=1, version=verLatest)
                fD = self.__sds.getFeature(seqId=seqId0, seqType="auth", partId=partId, altId=1, version=verLatest)
                seqFeature.clear()
                seqFeature.set(fD)
                (pId, seqBeg, seqEnd, _seqPartType) = seqFeature.getAuthPartDetails()

                if self.__debug:
                    self.__lfh.write("+AlignmentStatistics.__updateAuthXyzAlignments()  partId %d seqBegin %d seqEnd %d pId %d\n" % (partId, seqBeg, seqEnd, pId))

                r3L = []
                for tup in seqAuthIdx[seqBeg - 1 : seqEnd]:
                    r3L.append(str(tup[0]))
                pA.setReferenceSequence(r3L, "auth:" + str(seqId0) + "_P" + str(partId))
                refSeqLen = len(seqAuthIdx[seqBeg - 1 : seqEnd])

                # Get the sequence of each reference sequence in this group -
                #
                # List of reference sequences for this group only for the leading sequence -
                #
                sLen = {}
                altIdList = self.__sds.getAlternativeIds(seqId0, dataType="sequence", seqType="ref", partId=partId)
                if self.__debug:
                    self.__lfh.write("+AlignmentStatistics.__updateAuthXyzAlignments()  partId %d altIdList %r\n" % (partId, altIdList))

                if self.__verbose:
                    self.__lfh.write("+AlignmentStatistics.__updateAuthXyzAlignments()  entity group %s partId %d length altIds %d\n" % (gId, partId, len(altIdList)))

                #
                # Limit the number of reference sequences -
                #
                shrtAltIdList = altIdList[: self.__maxRefAlign]
                #
                for altId in shrtAltIdList:
                    verList = self.__sds.getVersionIds(seqId0, partId=partId, altId=altId, dataType="sequence", seqType="ref")
                    verLatest = verList[0]
                    seqRefIdx = self.__sds.getSequence(seqId=seqId0, seqType="ref", partId=partId, altId=altId, version=verLatest)
                    r3L = []
                    for tup in seqRefIdx:
                        r3L.append(str(tup[0]))
                    sLen[altId] = len(r3L)
                    pA.addTestSequence(r3L, "ref:" + str(altId) + "_P" + str(partId))

                #
                pA.doAlign()

                for altId in shrtAltIdList:
                    verList = self.__sds.getVersionIds(seqId0, partId=partId, altId=altId, dataType="sequence", seqType="ref")
                    if self.__debug:
                        self.__lfh.write("+AlignmentStatistics.__updateAuthXyzAlignments()  partId %d altId %d versions %r\n" % (partId, altId, verList))

                    verLatest = verList[0]
                    #
                    aL = pA.getAlignment("ref:" + str(altId) + "_P" + str(partId))
                    alignLength = len(aL)
                    numMatch = 0
                    numMatchGaps = 0
                    for aTup in aL:
                        if aTup[0] == aTup[1]:
                            numMatch += 1
                        if aTup[0] == aTup[1] or aTup[1] == self.__gapSymbol:
                            numMatchGaps += 1

                    resultD[(seqId0, "ref", partId, altId, verLatest)] = (sLen[altId], alignLength, numMatch, numMatchGaps)

                #
                if self.__verbose:
                    self.__lfh.write("\n+AlignmentStatistics.__updateAuthRefAlignments()  Alignment summary author sequence (length=%d) with reference sequences\n" % refSeqLen)
                if self.__debug:
                    for kTup, v in resultD.items():
                        if kTup[2] == partId:
                            self.__lfh.write(
                                "  -  ref seq  partId %2d alt %3d ver %d chainId %s length %5d align length %5d matches %5d matches w/ gaps %5d\n"
                                % (kTup[2], kTup[3], kTup[4], kTup[0], v[0], v[1], v[2], v[3])
                            )
        return resultD
