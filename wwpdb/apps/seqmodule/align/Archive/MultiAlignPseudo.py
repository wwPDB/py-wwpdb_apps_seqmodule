##
# File:    MultiAlignPseudo.py
# Date:    19-Jan-2010
#
# Updates:
#
# 20-Apr-2010 jdw Ported to module seqmodule.
# 20-Apr-2010 jdw Move rendering methods to AlignmentViewDepiction.
# 18-May-2010 jdw Add support for 'replaceid' type edits
# 27-Jul-2010 RPS Added support for accommodating different ordering of sequence types as per user preferences
# 11-Aug-2010 RPS Updated doAlignment() and getAlignmentData() to exclude gap conflicts from list of conflicts
# 08-Jan-2013 jdw include leading and trailing conflicts of any type.
# 27-Feb-2013 jdw replace ALignmentDataStore with new SequenceDataStore -
# 07-Mar-2013 jdw add sequence part support --
# 10-Mar-2013 jdw gut and refactor -- minimize global references  - fix conflict assignment
# 12-Mar-2013 jdw add method to save sequence ids for the current alignment
# 03-Apr-2013 jdw Rollup of many changes related to managing indexing and annotating conflicts.
# 19-Apr-2013 jdw Refactor with AlignmentUtils()
# 10-Dec-2013 jdw Renumber auth sequence after any edits.
# 19-May-2014 jdw Reorganize -
# 22-May-2014 jdw added export operation option to recover data for exporting assignments and mapping.
# 22-May-2014 jdw restored missing code
#  6-Jun-2014 jdw clear annotations in any self-referenced entities -
#  4-Jul-2014 jdw transfer comments after reannotation
# 30-Aug-2017 zf  change operation value from 'load' to 'loadandstore' for runSelected()
##
"""
Controlling class for the production of pseudo multiple-sequence alignments.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.09"

import os
import sys
import copy
import traceback

from wwpdb.apps.seqmodule.align.AlignmentStatistics import AlignmentStatistics
from wwpdb.apps.seqmodule.align.AlignmentUtils import AlignmentUtils

from wwpdb.apps.seqmodule.util.SequenceLabel import SequenceLabel

from wwpdb.apps.seqmodule.io.SequenceDataStore import SequenceDataStore
#


class MultiAlignPseudo(object):

    """ Controlling class for the production of pseduo multiple-sequence alignments.

        Supported operations:

         - load:                align the input sequence list, store and render the alignment.
         - realign:             apply the current set of stored edits and recompute and store the alignment.

    """

    def __init__(self, reqObj=None, verbose=False, log=sys.stderr):
        self.__verbose = verbose
        self.__debug = True
        self.__lfh = log
        self.__reqObj = reqObj

        #
        self.__alignOnReload = False
        #
        self.__operation = None
        self.__sessionObj = None
        #
        self.__sds = None
        self.__aU = AlignmentUtils(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        #
        self.__sessionObj = self.__reqObj.getSessionObj()
        self.__sessionPath = self.__sessionObj.getPath()
        self.__aD = {}
        self.__eD = {}
        #
        # List of tuples of revised sequences in the current alignment context.  [(id_org,id_rev),(id_org,id_rev),...]
        self.__modSeqIdList = []

    def runSelected(self, identifier):
        """ Create initial alignments using default selections --
        """

        if self.__verbose:
            self.__lfh.write("\n+MultiAlignPseudo.runSelected() STARTING for %r\n" % identifier)
        if self.__sds is None:
            self.__sds = SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        # if (self.__debug):
        #   self.__sds.dump(self.__lfh)
        #
        gIdList = self.__sds.getGroupIds()
        selectedIdList = self.__sds.getSelectedIds()
        #
        if self.__verbose:
            self.__lfh.write("+MultiAlignPseudo.runSelected() using stored gIdList %r\n" % gIdList)
            self.__lfh.write("+MultiAlignPseudo.runSelected() using stored selectedIdlist %r\n" % selectedIdList)

        for gId in gIdList:
            aL = self.run(operation='loadandstore', identifier=identifier, alignIdList=selectedIdList, alignGroupId=gId, selectedIdList=selectedIdList)

    def run(self, operation='load', identifier=None, alignIdList=None, alignGroupId=None, selectedIdList=None):
        """ Build or load the the alignment of input list of sequence identifiers -

        Supported operations:

         - load:                align the input sequence list and store the alignment.

         - realign:             apply the current set of stored edits to the current alignment,
                                also to corresponding sequences in the sequence data store,
                                and then recompute and store the alignment.

         Returns - alignment data list appropriate for rendering by the MultiAlignPseudoDepiciton() class.

             alignSeqList [[Sequence Id  from sequence data store (ie. auth_1_1_1_1),
                            SequenceLabel() object,
                            aligned sequence with index details,
                            conflict flag list (bool),
                            feature dictionary for input sequence],,,]

             the first sequence in the list is the reference sequence followed by the aligned test sequences -

        Each aligned sequence (w/index) contains a list of  -

        (one-letter-code, 3-letter-code, original label residue index, position in sequence, position in alignment, comment)

        where the comment contains one of  -

        ['engineered mutation','cloning artifact','variant','expression tag','insertion','deletion','microheterogeneity','chromophore',
           'linker','conflict','acetylation','amidation', 'initiating methionine']

         """
        if self.__verbose:
            self.__lfh.write("\n+MultiAlignPseudo.run() STARTING with operation %s for identifier %r entity %r\n" % (operation, identifier, alignGroupId))

        self.__alignGroupId = alignGroupId
        self.__selectedIdList = selectedIdList
        self.__identifier = identifier

        # if ((operation in  ["re-load", "realign"]) and (alignIdList is None or len(alignIdList) < 2):

        alignIdListF = self.__aU.filterSequenceIdList(seqIdList=alignIdList, groupId=alignGroupId)
        excludedPartIdList = self.__aU.getSelfReferencePartIdList(seqIdList=selectedIdList, groupId=alignGroupId)
        self.__aU.setExcludedPartIdList(excludedPartIdList)

        if self.__verbose:
            self.__lfh.write("+MultiAlignPseudo.run() input align sequence id list %r\n" % alignIdList)
            self.__lfh.write("+MultiAlignPseudo.run() filtered align sequence id list %r\n" % alignIdListF)
            self.__lfh.write("+MultiAlignPseudo.run() current selected sequence id list %r\n" % selectedIdList)
            self.__lfh.write("+MultiAlignPseudo.run() excluded part id list %r\n" % excludedPartIdList)
        #

        if operation in ["re-load", "realign"]:
            return self.__updateAndLoad(identifier=identifier, alignGroupId=alignGroupId, inputAlignIdList=alignIdListF)

        elif ((operation == "load-old") or (operation == "loadandstore")):
            return self.__load(identifier=identifier, alignGroupId=alignGroupId, inputAlignIdList=alignIdListF)

        elif ((operation == "load") or (operation == "loadandstore")):
            return self.__loadFromStore(identifier=identifier, alignGroupId=alignGroupId, inputAlignIdList=alignIdListF)

        elif ((operation == "export")):
            return self.__export(identifier=identifier, alignGroupId=alignGroupId, inputAlignIdList=alignIdListF)

        else:
            return []

    def getAlignGroupId(self):
        return self.__alignGroupId

    def getSelectedAnnotations(self):
        return self.__aD

    def getEntryDetails(self):
        return self.__eD

    def __loadFromStore(self, identifier, alignGroupId, inputAlignIdList):
        """
        """
        #     --  Check for a stored alignment for this input sequence list --
        #
        alignIdList, alignSeqList = self.__aU.getAlignmentData(identifier=identifier, entityId=alignGroupId)

        if self.__verbose:
            self.__lfh.write("\n+MultiAlignPseudo.__loadFromStore() starting for identifier %r entityGroup %r inputAlignIdList %r stored alignIdList %r\n"
                             % (identifier, alignGroupId, inputAlignIdList, alignIdList))
        #
        # Compare the input id list with return list -
        #
        if ((len(alignIdList) > 0) and (set(inputAlignIdList) == set(alignIdList))):
            if (self.__verbose):
                self.__lfh.write("+MultiAlignPseudo.__loadFromStore() USING EXISTING STORED alignment for sequence list: %r\n" % alignIdList)
            self.__aD = self.__aU.getSelectedAnnotations(inpAlignIdList=alignIdList)
            self.__eD = self.__aU.getEntryDetails()
            return alignSeqList
        else:
            if (self.__verbose):
                self.__lfh.write("+MultiAlignPseudo.__loadFromStore() RECALCULATING alignment for %r\n" % inputAlignIdList)
            return self.__load(identifier=identifier, alignGroupId=alignGroupId, inputAlignIdList=inputAlignIdList)

    def __load(self, identifier, alignGroupId, inputAlignIdList):
        """ Using the input sequence list  - get the sequences from the sequence store,
            assign the reference sequence, align the reference with the rest, and store the
            alignment in a data store with the input alignment tag.

            Return:  an sequence alignment list -

             alignSeqList [[Sequence Id  (ie. auth_A_1_1_1 as used in sequence data store),
                            SequenceLabel() object,
                            aligned sequence with index details,
                            conflict flag list (bool),
                            feature dictionary for input sequence],,,]

             the first sequence in the list is the reference sequence followed by the aligned test sequences -

            Aligned sequences with indices have storage model -

            (one-letter-code, 3-letter-code, original label residue index, position in sequence, position in alignment, comment )

        """
        #
        alignIdList = copy.deepcopy(inputAlignIdList)
        if (self.__verbose):
            self.__lfh.write("\n\n+MultiAlignPseudo.__load() Starting alignment for entity group  %s and id list %r\n" % (alignGroupId, inputAlignIdList))

        alignSeqList = self.__aU.doAlignment(inpAlignIdList=alignIdList)
        self.__aU.transferComments(alignSeqList)

        self.__aU.storeAlignmentData(alignSeqList=alignSeqList, identifier=identifier, entityId=alignGroupId)
        self.__aU.saveAlignmentIdList(alignIdList=alignIdList)
        self.__aD = self.__aU.getSelectedAnnotations(inpAlignIdList=alignIdList)
        self.__eD = self.__aU.getEntryDetails()

        if (self.__debug):
            self.__aU.formatAlignment(io=self.__lfh, alignSeqList=alignSeqList)

        if (self.__verbose):
            self.__lfh.write("+MultiAlignPseudo.__load() Completed for alignment for entity group %s with ids %r\n\n" %
                             (alignGroupId, alignIdList))
        return alignSeqList

    def __updateAndLoad(self, identifier, alignGroupId, inputAlignIdList):
        """  Apply any pending edits to the sequences within the current alignment, then apply those
             edits to sequences in the current sequence data store creating new versions for any edited sequences,
             and ONLY IF SPECIFIED realign the edited sequences.

             Return the revised sequence alignment list from the __load() method.
        """
        self.__modSeqIdList = []
        inputAlignIdList = copy.deepcopy(inputAlignIdList)
        if (self.__verbose):
            self.__lfh.write("\n\n+MultiAlignPseudo.__updateAndLoad() Starting with input alignment entity group %s with id list %r\n" %
                             (alignGroupId, inputAlignIdList))

        alignIdList, alignSeqList = self.__aU.getAlignmentData(identifier=identifier, entityId=alignGroupId)

        if (self.__verbose):
            self.__lfh.write("+MultiAlignPseudo.__updateAndLoad() Using input alignment entity group %s with stored id list %r\n" %
                             (alignGroupId, alignIdList))
        #
        # apply any pending edits to the current alignment --
        #
        edSeqIdList, multiInsertFlag = self.__aU.updateAlignment(alignSeqList, entityId=alignGroupId)
        if (self.__verbose):
            self.__lfh.write("+MultiAlignPseudo.__updateAndLoad() multiInsertFlag %r sequences edited in this alignment %d\n" %
                             (multiInsertFlag, len(edSeqIdList)))
        #
        #
        if (len(edSeqIdList) > 0):
            revAlignIdList, modSeqIdList = self.__aU.updateSequences(edSeqIdList, alignIdList, alignSeqList)
            #
            if (self.__alignOnReload):
                alstat = AlignmentStatistics(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
                alstat.doUpdate()
        else:
            revAlignIdList = alignIdList
            modSeqIdList = []

        self.__modSeqIdList = modSeqIdList
        if (self.__verbose):
            self.__lfh.write("+MultiAlignPseudo.__updateAndLoad() revised  alignment id list %r\n" % revAlignIdList)
            self.__lfh.write("+MultiAlignPseudo.__updateAndLoad() modified sequence identifiers %r\n" % modSeqIdList)
        #
        # Update any revised sequence identifiers in the alignment  +++
        #
        for (seqIdPrev, seqIdNext) in modSeqIdList:
            for ii in range(0, len(alignSeqList)):
                aTup = alignSeqList[ii]
                if aTup[0] == seqIdPrev:
                    aNew = list(aTup)
                    sLab = SequenceLabel()
                    aNew[0] = seqIdNext
                    sLab.unpack(seqIdNext)
                    aNew[1] = sLab
                    alignSeqList[ii] = aNew
        #
        #
        #
        if (self.__alignOnReload or multiInsertFlag):
            # JDW -- minority path now -- avoid the realignment --
            if (self.__verbose):
                self.__lfh.write("+MultiAlignPseudo.__updateAndLoad() re-alignment starting with id list %r\n" % revAlignIdList)
            oL = self.__load(identifier=identifier, alignGroupId=alignGroupId, inputAlignIdList=revAlignIdList)
        else:
            # JDW -- preserve the edited alignment -- July 4, 2014 JDW invert order
            self.__aU.reAnnotateAlignment(alignSeqList)
            self.__aU.transferComments(alignSeqList)
            oL = alignSeqList

        #
        if (self.__verbose):
            self.__lfh.write("+MultiAlignPseudo.__updateAndLoad() revised  alignment id list %r\n" % revAlignIdList)
        #
        if (self.__debug):
            self.__aU.formatAlignment(io=self.__lfh, alignSeqList=oL)

        self.__aU.storeAlignmentData(alignSeqList=oL, identifier=identifier, entityId=alignGroupId)
        self.__aU.saveAlignmentIdList(alignIdList=revAlignIdList)
        self.__aD = self.__aU.getSelectedAnnotations(inpAlignIdList=revAlignIdList)
        self.__eD = self.__aU.getEntryDetails()

        if (self.__verbose):
            self.__lfh.write("+MultiAlignPseudo.__updateAndLoad() Completed for alignment enity group %s with ids %r\n\n" %
                             (alignGroupId, revAlignIdList))
        return oL

    def updateSequenceIdList(self, seqIdList):
        oL = []
        for seqId in seqIdList:
            tId = seqId
            for modTup in self.__modSeqIdList:
                if tId == modTup[0]:
                    tId = modTup[1]
            oL.append(tId)
        #
        return oL

    def __export(self, identifier, alignGroupId, inputAlignIdList):
        """  Recover any existing alignnment for the input sequence id list.  Compute a new alignment if none exists.

             Return the a saved sequence alignment list or the recalculate this using the __load() method.
        """
        if (self.__verbose):
            self.__lfh.write("\n\n+MultiAlignPseudo.__export() Starting with input alignment entity group %s ids %r\n" %
                             (alignGroupId, inputAlignIdList))
        #
        # Check for a stored alignment -
        #
        alignIdList, alignSeqList = self.__aU.getAlignmentData(identifier=identifier, entityId=alignGroupId)
        #
        # compare the input id list with return list -
        #
        if ((len(alignIdList) > 0) and set(inputAlignIdList).issubset(alignIdList)):
            # The stored alignment will do -- reorder only --
            pass
        else:
            # Calculate the alignment and store result -
            alignSeqList = self.__load(identifier=identifier, alignGroupId=alignGroupId, inputAlignIdList=inputAlignIdList)
            self.__aU.storeAlignmentData(alignSeqList=alignSeqList, identifier=identifier, entityId=alignGroupId)

        #
        #
        if (self.__debug):
            self.__aU.formatAlignment(io=self.__lfh, alignSeqList=alignSeqList)

        #
        # reorder/filter elements --
        oL = []
        for id in inputAlignIdList:
            for atup in alignSeqList:
                if atup[0] == id:
                    # oL.append([atup[0],atup[2],atup[3]])
                    oL.append(atup)
                    break
        #
        if (self.__verbose):
            self.__lfh.write("+MultiAlignPseudo.__export() Completed for alignment enity group %s with ids %r\n\n" %
                             (alignGroupId, inputAlignIdList))
        return oL
