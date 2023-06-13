##
# File:    UpdatePolymerEntityPartitions.py
# Date:    22-Mar-2014
#
# Updates:
#   18-Apr-2014  jdw   add option for enity consolidation --
#    7-Aug-2015  jdw   upper case input one-letter-code sequence -
#    7-Sep-2017  zf    modified __updatePolymerEntityPartitions() to remove extra incorrect fragment assignment(s)
#    7-Aug-2021  zf    add option for changing _entity_poly.type
#
##
"""
Utilities for adding out-of-band sequence and feature data to the sequence data store.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.08"

import sys
import traceback
import copy

from wwpdb.apps.seqmodule.io.SequenceDataStore import SequenceDataStore
from wwpdb.apps.seqmodule.io.TaxonomyUtils import TaxonomyUtils
from wwpdb.apps.seqmodule.util.SequenceLabel import SequenceLabel, SequenceFeature
from wwpdb.apps.seqmodule.util.SequenceReferenceData import SequenceReferenceData


class UpdatePolymerEntityPartitions(object):

    """Utilities for updating polymer entity sequence and sequence partitioning and taxonomy assignment."""

    def __init__(self, reqObj=None, verbose=False, log=sys.stderr):
        self.__verbose = verbose
        self.__reqObj = reqObj
        self.__lfh = log
        #
        self.__srd = SequenceReferenceData(verbose=self.__verbose, log=self.__lfh)

        self.__setup()

    def __setup(self):
        try:
            self.__placeHolderValue = "click-to-edit"
            self.__siteId = self.__reqObj.getValue("WWPDB_SITE_ID")
            self.__sessionId = self.__reqObj.getSessionId()
            self.__sessionObj = self.__reqObj.getSessionObj()
            #
            self.__selectIdList = self.__reqObj.getSummarySelectList()
            self.__sds = SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            #
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+UpdatePolymerEntityPartitions.__setup() sessionId %s failed\n" % (self.__sessionObj.getId()))

    def makePolymerEntityPartEditForm(self, entityId, entryId=""):
        """
        <div class="ief" data-ief-edittype="select"
         data-ief-selectvalues="[{"value":"1","label":"Presentation Label","selected":true},{"value":"2","label":"Label 2","selected":false}]">

        <!-- #### merge option ####
        <table>
           <tr><td>Add instances from entity: </td>
               <td><span id="merge_instance_1"  class="ief greyedout" data-ief-edittype="select" data-ief-selectvalues='%s'>%s</span></td>
           </tr>
        </table>
         -->

        """
        top_form_template = """
        <div id="sectaxonomy">
        <h3>Sequence Parition and Taxonomy Data Form for Entry %s Entity %s</h3>
        <form name="formtaxonomy" id="formtaxonomy" action="/service/sequence_editor/respond_form/taxonomy" method="post" class="taxonomy_ajaxform">
            <input type="hidden" name="sessionid" value="%s" />
            <input type="hidden" name="entityid" value="%s" />
            <input type="hidden" name="numparts" value="%d" />

            <table>
            <tr><th>Entity Sequence</th><th>Entity Type</th></tr>
            <tr>
            <td style="text-align: left; font-family: monospace; white-space: pre;"><span id="entity_seq_1"  data-ief-edittype="textarea-mono" class="ief">%s</span></td>
            <td><span id="entity_type_1" class="%s" data-ief-edittype="select" data-ief-selectvalues='%s'>%s</span></td>
            </tr>
            </table>
            <br />
            <br />
            <table>
            <tbody id="seq_partition_table">
            <tr>
               <th>Part Id</th>
               <th>Taxonomy Id</th>
               <th>Seq Begin</th>
               <th>Seq End</th>
               <th>Part Type</th>
             </tr>
        """
        #
        bottom_form_template = """
            </tbody>
            </table>
        <br class="clearfloat" />
        <div class="width50 fltlft">Search Sequence Database: &nbsp;&nbsp; <input type="checkbox" name="seq_search_op" id="seq_search_op" /></div>
        <div class="width50 fltrgt"> Add instances from entity:&nbsp;&nbsp;
             <span id="merge_instance_1"  class="ief greyedout" data-ief-edittype="select" data-ief-selectvalues='%s'>%s</span></div>
        <br />
        <br />

        <div class="width50 fltlft"><input type="submit" name="submit" value="Submit edits" class="disableonclick submitparentform"  /></div>
        <div class="fltrgt"><input type="button" id="add_row_button" name="add_row_button" value="Add rows"/></div>
        <br class="clearfloat" />
            <!-- <input type="reset" name="reset" value="Reset" /> -->
            <input type="hidden" name="total_numparts" id="total_numparts" value="%s" />
            <input type="hidden" id="seq_length" value="%d" />
        </form>
        </div>
        """
        #
        # optList=['Biological sequence','Linker','N-terminal tag','C-terminal tag']
        optList = ["", "Biological sequence"]
        # optSel = ['false' for opt in optList]
        #
        entityTypeList = [
            "polypeptide(L)",
            "polypeptide(D)",
            "polydeoxyribonucleotide",
            "polyribonucleotide",
            "polydeoxyribonucleotide/polyribonucleotide hybrid",
            "cyclic-pseudo-peptide",
            "peptide nucleic acid",
            "other",
        ]
        #
        seq_length, seq1, partD, entityType = self.__getEntityPartDetails(entityId)
        if self.__verbose:
            self.__lfh.write("+UpdatePolymerEntityPartitions.makePolymerEntityPartEditForm() entity partition starting data\n")
            self.__lfh.write(" Sequence:\n%s\n" % seq1)
            for k, v in partD.items():
                self.__lfh.write(" part %r  data:  %r\n" % (k, v))
            #
        #
        partIdList = list(partD.keys())
        partIdList.sort(key=int)
        #
        entityDisplayClass = "ief"
        if entityType == self.__placeHolderValue:
            entityDisplayClass = "ief greyedout"
        #
        entityTypeTxt = self.__formatSelectList(entityTypeList, entityType)
        oL = []
        oL.append(top_form_template % (entryId, entityId, self.__sessionId, entityId, len(partIdList), seq1, entityDisplayClass, entityTypeTxt, entityType))
        #
        for partId in partIdList:
            _pIdT, seqBeg, seqEnd, pType, taxIdT = partD[partId]
            oL.append("<tr>")
            oL.append('<td><span id="p_%d_partid">%d</span><input type="hidden" name="p_%d_partid" value="%d" /></td>' % (partId, partId, partId, partId))
            if (taxIdT is not None) and (len(taxIdT) > 0):
                taxId = taxIdT
                oL.append('<td><span id="p_%d_taxid"    class="ief">%s</span></td>' % (partId, taxId))
            else:
                taxId = self.__placeHolderValue
                oL.append('<td><span id="p_%d_taxid"    class="ief greyedout">%s</span></td>' % (partId, taxId))
            oL.append('<td><span id="p_%d_seqbegin" class="ief">%s</span></td>' % (partId, seqBeg))
            oL.append('<td><span id="p_%d_seqend"   class="ief">%s</span></td>' % (partId, seqEnd))
            jTxt = self.__formatSelectList(optList, pType)
            oL.append('<td><span id="p_%d_seqtype"  class="ief" data-ief-edittype="select" data-ief-selectvalues=\'%s\'>%s</span></td>' % (partId, jTxt, pType))
            oL.append("</tr>")

        partId = int(partIdList[-1])
        taxId = self.__placeHolderValue
        seqBeg = self.__placeHolderValue
        seqEnd = self.__placeHolderValue
        pType = self.__placeHolderValue

        for _i in range(1, 8):
            partId += 1
            oL.append("<tr>")
            oL.append('<td><span id="p_%d_partid">%d</span><input type="hidden" name="p_%d_partid" value="%d" /></td>' % (partId, partId, partId, partId))
            oL.append('<td><span id="p_%d_taxid"    class="ief greyedout">%s</span></td>' % (partId, taxId))
            oL.append('<td><span id="p_%d_seqbegin" class="ief greyedout">%s</span></td>' % (partId, seqBeg))
            oL.append('<td><span id="p_%d_seqend"   class="ief greyedout">%s</span></td>' % (partId, seqEnd))
            jTxt = self.__formatSelectList(optList, pType)
            oL.append('<td><span id="p_%d_seqtype"  class="ief greyedout" data-ief-edittype="select" data-ief-selectvalues=\'%s\'>%s</span></td>' % (partId, jTxt, pType))
            oL.append("</tr>")
        #

        eIdList = self.__sds.getGroupIds()
        myList = []
        myList.append(self.__placeHolderValue)
        myList.extend([eId for eId in sorted(eIdList) if eId != entityId])
        jTxt = self.__formatSelectList(myList, "None")
        oL.append(bottom_form_template % (jTxt, self.__placeHolderValue, partId, seq_length))

        return "\n".join(oL)
        #

    def polymerEntityPartEditFormResponder(self):
        """Update the polymer entity data store using user provided entity part, source and taxonomy content.

        Form data encoded in the input request object --
        """
        try:
            numParts = int(str(self.__reqObj.getValue("numparts")))
            totalNumParts = int(str(self.__reqObj.getValue("total_numparts")))
            entityId = self.__reqObj.getValue("entityid")
            seq1 = str(self.__reqObj.getValue("entity_seq_1")).upper().strip()
            entity_type_1 = str(self.__reqObj.getValue("entity_type_1")).strip()
            #
            mergeEntityId = self.__reqObj.getValue("merge_instance_1")
            if mergeEntityId != self.__placeHolderValue:
                self.__lfh.write("+UpdatePolymerEntityPartitions.polymerEntityPartEditFormResponder() merging entity id %r\n" % mergeEntityId)
                self.__updateGroupInstances(entityId, mergeEntityId)

            #
            pD = {}
            # for partId in range(1, numParts + 8):
            for partId in range(1, totalNumParts + 1):
                taxId = self.__reqObj.getValue("p_%d_taxid" % partId)
                if (taxId is None) or (taxId == self.__placeHolderValue):
                    taxId = ""
                #
                seqBegin = self.__reqObj.getValue("p_%d_seqbegin" % partId)
                seqEnd = self.__reqObj.getValue("p_%d_seqend" % partId)
                seqPartType = self.__reqObj.getValue("p_%d_seqtype" % partId)
                if (
                    (not seqPartType)
                    or (seqPartType == self.__placeHolderValue)
                    or (not seqBegin)
                    or (seqBegin == self.__placeHolderValue)
                    or (not seqEnd)
                    or (seqEnd == self.__placeHolderValue)
                ):
                    continue
                #
                if partId > numParts:
                    pD[partId] = (partId, str(seqBegin), str(seqEnd), str(seqPartType), str(taxId))
                else:
                    pD[partId] = (partId, int(seqBegin), int(seqEnd), str(seqPartType), str(taxId))
                #
            #
            self.__lfh.write("+UpdatePolymerEntityPartitions.polymerEntityPartEditFormResponder() pD=%d\n" % len(pD))
            _seq_length, seq1Org, pOrgD, entityType = self.__getEntityPartDetails(entityId=entityId, numExtra=(totalNumParts - numParts))
            #
            if seq1 != seq1Org:
                self.__lfh.write("+UpdatePolymerEntityPartitions.polymerEntityPartEditFormResponder() sequence has changed\n")
                self.__lfh.write("+UpdatePolymerEntityPartitions.polymerEntityPartEditFormResponder() sequence org:\n%s\n" % seq1Org)
                self.__lfh.write("+UpdatePolymerEntityPartitions.polymerEntityPartEditFormResponder() sequence new:\n%s\n" % seq1)
                updateSequenceFlag = True
            else:
                updateSequenceFlag = False
                self.__lfh.write("+UpdatePolymerEntityPartitions.polymerEntityPartEditFormResponder() sequence is unchanged\n")
            #
            updatePolymerTypeFlag = False
            if entity_type_1 and (entity_type_1 != entityType):
                updatePolymerTypeFlag = True
            #
            updatePartitionFlag = False
            # for partId in range(1, numParts + 8):
            for partId in range(1, totalNumParts + 1):
                # if pOrgD[partId] != pD[partId] and partId <= numParts:
                if (partId in pOrgD) and (partId in pD) and pOrgD[partId] != pD[partId]:
                    self.__lfh.write(
                        "+UpdatePolymerEntityPartitions.polymerEntityPartEditFormResponder() source differs at partId %d current %r next %r\n" % (partId, pOrgD[partId], pD[partId])
                    )
                    updatePartitionFlag = True
                    break
                # elif pD[partId][1] != self.__placeHolderValue and partId > numParts:
                elif ((partId in pD) and (partId not in pOrgD)) or ((partId not in pD) and (partId in pOrgD)):
                    updatePartitionFlag = True
                    break
                #
            #
            if updatePartitionFlag or updatePolymerTypeFlag or updateSequenceFlag:
                if self.__verbose:
                    self.__lfh.write("+UpdatePolymerEntityPartitions.polymerEntityPartEditFormResponder() source data for entity %s updated with %r\n" % (entityId, pD))
                self.__updatePolymerEntityPartitions(
                    entityId=entityId, partD=pD, seq1=seq1, polymerType=entity_type_1, updateSequenceFlag=updateSequenceFlag, updatePolymerTypeFlag=updatePolymerTypeFlag
                )
            else:
                if self.__verbose:
                    self.__lfh.write("+UpdatePolymerEntityPartitions.polymerEntityPartEditFormResponder() form values unchanged\n")
                #
            #
        except:  # noqa: E722 pylint: disable=bare-except
            self.__lfh.write("+UpdatePolymerEntityPartitions.sourceEditResponder() failing\n")
            traceback.print_exc(file=self.__lfh)
            return False

        return True

    def __getCurrentAuthSelection(self, entityId, partId=1):
        """Search selections for author sequence part 1. in order to establish the
        the selected version.   Then return the objects associated with the input
        partId.
        """
        try:
            # seqIdList = self.__sds.getGroup(entityId)
            # JDW CHANGE
            # instanceId=seqIdList[0]
            instanceId = entityId
            sL = SequenceLabel()
            for sId in self.__selectIdList:
                if sId.startswith("auth"):
                    sL.unpack(sId)
                    _seqType, seqInstId, seqPartId, seqAltId, seqVersion = sL.get()
                    if instanceId == seqInstId and seqPartId == 1:
                        fdObj = self.__sds.getFeatureObj(seqInstId, seqType="auth", partId=partId, altId=seqAltId, version=seqVersion)
                        self.__lfh.write(
                            "+UpdatePolymerEntityPartitions._getCurrentAuthSelection() returns: entity %r instance %r partId %r altId %r version %r\n"
                            % (entityId, seqInstId, seqPartId, seqAltId, seqVersion)
                        )
                        return sId, sL, fdObj
        except:  # noqa: E722 pylint: disable=bare-except
            self.__lfh.write("+UpdatePolymerEntityPartitions._getCurrentAuthSelection() failed for selectList %r entityId %r partId %r\n" % (self.__selectIdList, entityId, partId))
            traceback.print_exc(file=self.__lfh)

        self.__lfh.write("+UpdatePolymerEntityPartitions._getCurrentAuthSelection() no return for entityId %r partId %r\n" % (entityId, partId))
        return None, None, None

    # def __assignPart(self, partD, seqBegin, seqEnd):
    #     """Assign the entity partId to the"""
    #     pMatch = -1
    #     for partId, pTup in partD.items():
    #         if (seqBegin == pTup[1]) and (seqEnd == pTup[2]):
    #             pMatch = partId

    #     return pMatch

    def __formatSelectList(self, pList, pType):
        oL = []
        pSelectList = ["false" for p in pList]
        if pType is not None and len(pType) > 1:
            pListU = [p.upper() for p in pList]
            pTypeU = pType.upper()
            try:
                idx = pListU.index(pTypeU)
                pSelectList[idx] = "true"
            except ValueError:
                idx = -1
        tL = []
        for pt, psel in zip(pList, pSelectList):
            tL.append('{"value":"%s","label":"%s","selected":%s}' % (pt, pt, psel))
        oL.append("[")
        oL.append(",".join(tL))
        oL.append("]")

        return "".join(oL)

    def __getAuthSeqComment(self, idx, partD, seqBegMin, seqEndMax):
        """Set default annotation between part -- and clean up any other comments --
        Return initial annotation assignment for residue position 'idx' based on the
        the input part definition OR an out-of-part assignment of a terminal expression
        tag or an internal linker.

        input    idx   sequence position in author/sample sequence (1-N)
                 partD[partId]=(id,seqBeg,seqEnd,partTypeComment,...)

        """
        # Are we in defined entity part ?
        for _pId, rl in partD.items():
            if idx < rl[1] or idx > rl[2]:
                continue
            else:
                return rl[3]

        # not in part --
        if idx < seqBegMin or idx > seqEndMax:
            return "expression tag"
        else:
            return "linker"

    # def __updateAuthPartDetails(self, entityId, partD):
    #     #
    #     # create new version of all parts --
    #     #
    #     if self.__verbose:
    #         self.__lfh.write("\n\n+UpdatePolymerEntityPartitions.__updateAuthPartDetails() Starting for entity %s and part dictionary %r\n" % (entityId, partD))
    #     #
    #     seqBegMin = 100
    #     seqEndMax = 0
    #     pD = {}
    #     for tId, pTup in partD.items():
    #         pId, seqBeg, seqEnd, seqPartType, taxId = pTup
    #         if self.__verbose:
    #             self.__lfh.write("+UpdatePolymerEntityPartitions.__updateAuthPartDetails() Entity %s part %d new assignment form values %s\n" % (entityId, pId, pTup))
    #         #
    #         if (seqPartType == self.__placeHolderValue) or (seqBeg == self.__placeHolderValue) or (seqEnd == self.__placeHolderValue):
    #             continue
    #         seqNumBeg = int(seqBeg)
    #         seqNumEnd = int(seqEnd)
    #         seqBegMin = min(seqBegMin, seqNumBeg)
    #         seqEndMax = max(seqEndMax, seqNumEnd)
    #         spt = ""
    #         if seqPartType in ["biological sequence"]:
    #             spt = ""
    #         elif seqPartType in ["n-terminal tag", "c-terminal tag"]:
    #             spt = "expression tag"
    #         elif seqPartType in ["linker"]:
    #             spt = "linker"
    #         pD[pId] = (pId, seqNumBeg, seqNumEnd, spt)

    #     #
    #     if self.__verbose:
    #         self.__lfh.write("\n\n+UpdatePolymerEntityPartitions.__updateAuthPartDetails() entity %s seqBegMin %d   seqEndMax %d\n" % (entityId, seqBegMin, seqEndMax))

    #     tU = TaxonomyUtils(siteId=self.__siteId, verbose=self.__verbose, log=self.__lfh)

    #     # seqIds = self.__sds.getGroup(groupId=entityId)

    #     # if len(seqIds)==0:
    #     #    return False

    #     # seqFeature = SequenceFeature()
    #     # JDW CHANGE
    #     # seqId0=seqIds[0]
    #     seqId0 = entityId
    #     curPartIdList = self.__sds.getPartIds(seqId0, dataType="sequence", seqType="auth")

    #     # base versioning on the base partId = 1
    #     #
    #     vL = self.__sds.getVersionIds(seqId0, partId=1, altId=1, dataType="sequence", seqType="auth")
    #     if len(vL) < 1:
    #         return False
    #     # target version for update.
    #     # lastVer=int(str(vL[0]))
    #     curVer = int(str(vL[0]))
    #     nextVer = curVer + 1
    #     # JDW      authSelectId,authSL,authFdObj=self.__getCurrentAuthSelection(entityId,partId=partId)
    #     #          curSeqType,curInstId,curPartId,curAltId,curVer=authSL.get()
    #     #
    #     for partId, pTup in partD.items():
    #         authSelectId, authSL, authFdObj = self.__getCurrentAuthSelection(entityId, partId=partId)
    #         curSeqType, curInstId, curPartId, curAltId, curVer = authSL.get()
    #         if partId in curPartIdList:
    #             # Existing part ---
    #             # fObj=self.__sds.getFeatureObj(seqId0,seqType="auth",partId=partId,altId=1,version=curVer)
    #             taxIdOrg = authFdObj.getSourceTaxId()
    #             # values -
    #             pId, seqBeg, seqEnd, seqPartType, taxId = pTup
    #             if self.__verbose:
    #                 self.__lfh.write(
    #                     "+UpdatePolymerEntityPartitions.__updateAuthPartDetails() Taxonomy assignments for part %d org taxid %s new taxid %s\n" % (partId, taxIdOrg, taxId)
    #                 )
    #             if taxIdOrg != taxId:
    #                 nL = tU.lookUpSource(taxId=taxId)
    #                 if len(nL) > 0:
    #                     authFdObj.setSource(organism=nL[0], taxid=taxId)
    #                 else:
    #                     authFdObj.setTaxId(taxid=taxId)
    #             authFdObj.setAuthPartDetails(partId, seqBeg, seqEnd, seqPartType)
    #             self.__sds.setFeatureObj(authFdObj, seqId0, seqType="auth", partId=partId, altId=1, version=nextVer)
    #             # Update the sequence...
    #             seqIdx = self.__sds.getSequence(seqId0, seqType="auth", partId=partId, altId=1, version=curVer)
    #             ####
    #             # JDW reset default comment state ---
    #             sTupL = []
    #             for sT in seqIdx:
    #                 comment = self.__getAuthSeqComment(sT[3], pD, seqBegMin, seqEndMax)
    #                 sTupL.append((sT[0], sT[1], comment, sT[3], sT[4], sT[5]))
    #             ####
    #             # self.__sds.setSequence(seqIdx,seqId0,seqType='auth',partId=partId,altId=1,version=nextVer)
    #             self.__sds.setSequence(sTupL, seqId0, seqType="auth", partId=partId, altId=1, version=nextVer)
    #             if self.__verbose:
    #                 self.__lfh.write("+UpdatePolymerEntityPartitions.__updateAuthPartDetails() updating part %d from version %d to %d\n" % (partId, curVer, nextVer))
    #         else:
    #             # New part ---
    #             pId, seqBeg, seqEnd, seqPartType, taxId = pTup
    #             if self.__verbose:
    #                 self.__lfh.write("+UpdatePolymerEntityPartitions.__updateAuthPartDetails() Entity %s part %d new assignment form values %s\n" % (entityId, partId, pTup))
    #             #
    #             if (seqPartType == self.__placeHolderValue) or (seqBeg == self.__placeHolderValue) or (seqEnd == self.__placeHolderValue):
    #                 continue

    #             #
    #             if (taxId is None) or (taxId == self.__placeHolderValue):
    #                 taxId = ""
    #                 if self.__verbose:
    #                     self.__lfh.write("+UpdatePolymerEntityPartitions.__updateAuthPartDetails() Handle missing taxid for entity %s part %s\n" % (entityId, pId))
    #             #
    #             # Update the sequence ...
    #             #
    #             seqIdx = self.__sds.getSequence(seqId0, seqType="auth", partId=1, altId=1, version=curVer)

    #             # JDW reset default comment state ---
    #             sTupL = []
    #             for sT in seqIdx:
    #                 comment = self.__getAuthSeqComment(sT[3], pD, seqBegMin, seqEndMax)
    #                 sTupL.append((sT[0], sT[1], comment, sT[3], sT[4], sT[5]))

    #             # self.__sds.setSequence(seqIdx,seqId0,seqType='auth',partId=partId,altId=1,version=nextVer)
    #             # JDW
    #             self.__sds.setSequence(sTupL, seqId0, seqType="auth", partId=partId, altId=1, version=nextVer)
    #             #
    #             # Update features ... Using data from the principal part.
    #             #
    #             # fObj=SequenceFeature()
    #             authSelectId1, authSL1, fObj = self.__getCurrentAuthSelection(entityId, partId=1)
    #             #
    #             nL = tU.lookUpSource(taxId=taxId)
    #             if len(nL) > 0:
    #                 fObj.setSource(organism=nL[0], taxid=taxId)
    #             else:
    #                 fObj.setTaxId(taxid=taxId)
    #             fObj.setAuthPartDetails(partId, seqBeg, seqEnd, seqPartType)
    #             #
    #             #  JDW also need to update any other entity-level details ---
    #             #
    #             self.__sds.setFeatureObj(fObj, seqId0, seqType="auth", partId=partId, altId=1, version=nextVer)
    #             if self.__verbose:
    #                 self.__lfh.write("+UpdatePolymerEntityPartitions.polymerEntityPartEditFormResponder() entity %s new part %d with version %d\n" % (entityId, partId, nextVer))
    #             #
    #         #
    #     #
    #     self.__sds.serialize()
    #     return True

    def __getEntityPartDetails(self, entityId, numExtra=0):
        """Return the entity sequence and a dictionary of part boundaries and types ---   pD[pId]=(pId,pSeqBegin,pSeqEnd,pType,taxId)"""

        pD = {}
        # seqIds = self.__sds.getGroup(groupId=entityId)
        # JDW CHANGE
        # if len(seqIds)==0:
        #    return pD

        seqFeature = SequenceFeature()
        #
        # JDW CHANGE
        # seqId0=seqIds[0]
        #
        seqId0 = entityId

        # polymerTypeCode = "AA"
        entityType = self.__placeHolderValue
        #
        partIdList = self.__sds.getPartIds(seqId0, dataType="sequence", seqType="auth")

        for partId in partIdList:
            vL = self.__sds.getVersionIds(seqId0, partId=partId, altId=1, dataType="sequence", seqType="auth")
            if len(vL) > 0:
                pfD = self.__sds.getFeature(seqId0, seqType="auth", partId=partId, altId=1, version=vL[0])
                seqFeature.set(pfD)
                pId, pSeqBegin, pSeqEnd, pType = seqFeature.getAuthPartDetails()
                # polymerTypeCode = seqFeature.getPolymerType()
                if entityType == self.__placeHolderValue:
                    polymerType = seqFeature.getPolymerLinkingType()
                    if polymerType:
                        entityType = polymerType
                    #
                #
                taxId = seqFeature.getSourceTaxId()
                pD[pId] = (pId, pSeqBegin, pSeqEnd, pType, taxId)
                lastPart = pId

        if numExtra > 0:
            pv = self.__placeHolderValue
            pId = lastPart
            for _ii in range(0, numExtra):
                pId += 1
                pD[pId] = (pId, pv, pv, pv, pv)
        #
        seqAuthIdx = self.__sds.getSequence(seqId=seqId0, seqType="auth", partId=1, altId=1, version=vL[0])
        r3List = []
        for sTup in seqAuthIdx:
            (r3, _sIdx, _comment, _idx, _r1, _org_r3) = sTup
            r3List.append(r3)
        seq1 = self.__srd.cnvList3to1WithModsFormatted(r3List, maxLine=60)
        #
        return len(r3List), seq1, pD, entityType

    def __updatePolymerEntityPartitions(self, entityId, partD, seq1=None, polymerType=None, updateSequenceFlag=False, updatePolymerTypeFlag=False):
        """Create new versions of sequence and all related parts --"""
        if self.__verbose:
            self.__lfh.write(
                "\n\n+UpdatePolymerEntityPartitions.__updatePolymerEntityPartitions() Starting for entity %s sequence update flag %r part dictionary\n%r\n"
                % (
                    entityId,
                    updateSequenceFlag,
                    partD,
                )
            )
        #
        seqBegMin = 100
        seqEndMax = 0
        nP = 0
        pD = {}
        for _tId, pTup in partD.items():
            pId, seqBeg, seqEnd, seqPartType, taxId = pTup
            if self.__verbose:
                self.__lfh.write("+UpdatePolymerEntityPartitions.__updatePolymerEntityPartitions() Entity %s part %d new assignment form values %s\n" % (entityId, pId, pTup))
            #
            if (seqPartType == self.__placeHolderValue) or (seqBeg == self.__placeHolderValue) or (seqEnd == self.__placeHolderValue):
                continue
            #
            seqNumBeg = int(seqBeg)
            seqNumEnd = int(seqEnd)
            seqBegMin = min(seqBegMin, seqNumBeg)
            seqEndMax = max(seqEndMax, seqNumEnd)
            spt = ""
            if seqPartType in ["biological sequence"]:
                spt = ""
            elif seqPartType in ["n-terminal tag", "c-terminal tag"]:
                spt = "expression tag"
            elif seqPartType in ["linker"]:
                spt = "linker"
            #
            pD[pId] = (pId, seqNumBeg, seqNumEnd, spt)
            nP += 1
        #
        if self.__verbose:
            self.__lfh.write("\n\n+UpdatePolymerEntityPartitions.__updatePolymerEntityPartitions() entity %s seqBegMin %d   seqEndMax %d\n" % (entityId, seqBegMin, seqEndMax))
        #

        tU = TaxonomyUtils(siteId=self.__siteId, verbose=self.__verbose, log=self.__lfh)

        # seqIds = self.__sds.getGroup(groupId=entityId)
        # JDW CHANGE
        # if len(seqIds)==0:
        #    return False

        # seqFeature = SequenceFeature()
        # JDW CHANGE
        # seqId0=seqIds[0]
        seqId0 = entityId
        curPartIdList = self.__sds.getPartIds(seqId0, dataType="sequence", seqType="auth")

        # Versioning of 'auth' sequence is based partId = 1
        #
        vL = self.__sds.getVersionIds(seqId0, partId=1, altId=1, dataType="sequence", seqType="auth")
        if len(vL) < 1:
            return False
        #
        curVer = int(str(vL[0]))
        nextVer = curVer + 1
        self.__lfh.write("+UpdatePolymerEntityPartitions.__updatePolymerEntityPartitions() nextVer=%d\n" % nextVer)
        #
        seqIdxNext = []
        if updateSequenceFlag:
            # First update the sequences for each part IF the sequence has changed -
            _authSelectId, authSL, authFdObj = self.__getCurrentAuthSelection(entityId, partId=1)
            polyTypeCode = authFdObj.getPolymerType()
            if updatePolymerTypeFlag:
                polyTypeCode = self.__srd.getPolymerTypeCode(polymerType)
            #
            (_r1L, r3L) = self.__srd.parseSequence(str(seq1).upper(), polyTypeCode)
            #
            #  Handle any change in sequence length at the endpoints --
            #
            seqEndMax = len(r3L)
            (pId, seqNumBeg, seqNumEnd, spt) = pD[nP]
            pD[nP] = (pId, seqNumBeg, seqEndMax, spt)
            pId, seqBeg, seqEnd, seqPartType, taxId = partD[nP]
            partD[nP] = (pId, seqBeg, seqEndMax, seqPartType, taxId)
            #
            ir = 1
            for r3 in r3L:
                comment = self.__getAuthSeqComment(ir, pD, seqBegMin, seqEndMax)
                seqIdxNext.append((r3, str(ir), comment, ir, self.__srd.cnv3To1(r3), r3))
                ir += 1
            #
            for pId in pD.keys():
                self.__sds.setSequence(seqIdxNext, seqId0, "auth", partId=pId, altId=1, version=nextVer)
            #
        #
        for partId, pTup in partD.items():
            _authSelectId, authSL, authFdObj = self.__getCurrentAuthSelection(entityId, partId=partId)
            _curSeqType, _curInstId, _curPartId, _curAltId, curVer = authSL.get()
            if partId in curPartIdList:
                # Existing part ---
                # fObj=self.__sds.getFeatureObj(seqId0,seqType="auth",partId=partId,altId=1,version=curVer)
                taxIdOrg = authFdObj.getSourceTaxId()
                # values -
                pId, seqBeg, seqEnd, seqPartType, taxId = pTup
                if self.__verbose:
                    self.__lfh.write(
                        "+UpdatePolymerEntityPartitions.__updatePolymerEntityPartitions() Taxonomy assignments for part %d org taxid %s new taxid %s\n" % (partId, taxIdOrg, taxId)
                    )
                if taxIdOrg != taxId:
                    nL = tU.lookUpSource(taxId=taxId)
                    if len(nL) > 0:
                        authFdObj.setSource(organism=nL[0], taxid=taxId)
                    else:
                        authFdObj.setTaxId(taxid=taxId)
                    #
                #
                authFdObj.setAuthPartDetails(partId, seqBeg, seqEnd, seqPartType)
                if updatePolymerTypeFlag:
                    authFdObj.setPolymerLinkingType(polymerType)
                    authFdObj.setPolymerType(self.__srd.getPolymerTypeCode(polymerType))
                #
                self.__sds.setFeatureObj(authFdObj, seqId0, seqType="auth", partId=partId, altId=1, version=nextVer)
                #
                if not updateSequenceFlag:
                    #  Update the annotation of the prior sequence ...
                    seqIdx = self.__sds.getSequence(seqId0, seqType="auth", partId=partId, altId=1, version=curVer)
                    sTupL = []
                    for sT in seqIdx:
                        comment = self.__getAuthSeqComment(sT[3], pD, seqBegMin, seqEndMax)
                        sTupL.append((sT[0], sT[1], comment, sT[3], sT[4], sT[5]))
                    self.__sds.setSequence(sTupL, seqId0, seqType="auth", partId=partId, altId=1, version=nextVer)
                if self.__verbose:
                    self.__lfh.write("+UpdatePolymerEntityPartitions.__updatePolymerEntityPartitions() updating part %d from version %d to %d\n" % (partId, curVer, nextVer))
            else:
                #  Adding a new entity part ---
                pId, seqBeg, seqEnd, seqPartType, taxId = pTup
                if self.__verbose:
                    self.__lfh.write(
                        "+UpdatePolymerEntityPartitions.__updatePolymerEntityPartitions() Entity %s part %d new assignment form values %s\n" % (entityId, partId, pTup)
                    )
                #
                if (seqPartType == self.__placeHolderValue) or (seqBeg == self.__placeHolderValue) or (seqEnd == self.__placeHolderValue):
                    continue
                #
                if (taxId is None) or (taxId == self.__placeHolderValue):
                    taxId = ""
                    if self.__verbose:
                        self.__lfh.write("+UpdatePolymerEntityPartitions.__updatePolymerEntityPartitions() Handle missing taxid for entity %s part %s\n" % (entityId, pId))
                    #

                #
                # Update the sequence ... if it has not already been done above.
                #
                if not updateSequenceFlag:
                    seqIdx = self.__sds.getSequence(seqId0, seqType="auth", partId=1, altId=1, version=curVer)
                    sTupL = []
                    for sT in seqIdx:
                        comment = self.__getAuthSeqComment(sT[3], pD, seqBegMin, seqEndMax)
                        sTupL.append((sT[0], sT[1], comment, sT[3], sT[4], sT[5]))
                    #
                    self.__sds.setSequence(sTupL, seqId0, seqType="auth", partId=partId, altId=1, version=nextVer)
                #
                # Update features ... Using data from the principal part.
                #
                # fObj=SequenceFeature()
                _authSelectId1, _authSL1, fObj = self.__getCurrentAuthSelection(entityId, partId=1)
                #
                nL = tU.lookUpSource(taxId=taxId)
                if len(nL) > 0:
                    fObj.setSource(organism=nL[0], taxid=taxId)
                else:
                    fObj.setTaxId(taxid=taxId)
                #
                fObj.setAuthPartDetails(partId, seqBeg, seqEnd, seqPartType)
                if updatePolymerTypeFlag:
                    fObj.setPolymerLinkingType(polymerType)
                    fObj.setPolymerType(self.__srd.getPolymerTypeCode(polymerType))
                #
                #  JDW also need to update any other entity-level details ---
                #
                self.__sds.setFeatureObj(fObj, seqId0, seqType="auth", partId=partId, altId=1, version=nextVer)
                if self.__verbose:
                    self.__lfh.write("+UpdatePolymerEntityPartitions.polymerEntityPartEditFormResponder() entity %s new part %d with version %d\n" % (entityId, partId, nextVer))
                #
            #
        #
        for partId in curPartIdList:
            if partId in partD.keys():
                continue
            #
            self.__sds.removePartId(seqId0, dataType="sequence", seqType="auth", partId=partId)
        #
        self.__sds.serialize()
        return True

    def __updateGroupInstances(self, entityId, mergeEntityId):
        """Move/add instances from mergeEntity to entityId"""
        try:
            currentIds = self.__sds.getGroup(groupId=entityId)
            mergeIds = self.__sds.getGroup(groupId=mergeEntityId)
            tMergeIds = copy.deepcopy(mergeIds)
            self.__lfh.write("+UpdatePolymerEntityPartitions.__updateGroupInstances() entity %r instance list:  %r\n" % (entityId, currentIds))
            self.__lfh.write("+UpdatePolymerEntityPartitions.__updateGroupInstances() merge entity %r instance list:  %r\n" % (mergeEntityId, tMergeIds))
            for mId in tMergeIds:
                if mId not in currentIds:
                    currentIds.append(mId)
                    mergeIds.remove(mId)
            self.__sds.setGroup(entityId, currentIds)
            self.__sds.setGroup(mergeEntityId, mergeIds)
            self.__sds.serialize()
            self.__lfh.write("+UpdatePolymerEntityPartitions.__updateGroupInstances() entity %r instance list updated:  %r\n" % (entityId, currentIds))
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            self.__lfh.write("+UpdatePolymerEntityPartitions.__updateGroupInstances() failed for entity %r merge entity %r\n" % (entityId, mergeEntityId))
            traceback.print_exc(file=self.__lfh)
        return False


if __name__ == "__main__":
    pass
