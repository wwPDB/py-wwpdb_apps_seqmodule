##
# File:  SequenceDataExport_v2.py
# Date:  14-Feb-2010
#
# Updates:
#   20-Apr-2010 jdw Ported to module seqmodule.
#   23-Apr-2010 jdw Added method for self-references.
#   02-May-2010 jdw Add file path to constructor, move makeDefaultSelection() to SequenceSelection class.
#                   Add handling of self-references...
#  01-Mar-2013  jdw remove SiteInterface -
#  08-Mar-2013  jdw adpot the standard file name conventions.
#  08-Mar-2013  jdw add support use of sequence parts.
#  07-Apr-2013  jdw add support for isoforms/filter self references from db references
#  23-Apr-2013  jdw add separate category for export of modified residues
#  05-May-2013  jdw filter gap residues from author/xyz mapping list
#  20-Sep-2013  jdw fix filter for mapping from auth/xyz mapping
#  14-Nov-2013  jdw add new entity details export
#  18-Dec-2013  jdw adjust conflict array lengths after alignment --
#  19-Jan-2014  jdw "HOST_ORG_CELL_LINE"
#  27-Apr-2014  jdw working on fix to source detail export -
#  28-Apr-2014  jdw verify that updateArchive=False is default --
#                   Changes will always be applied to the first model in the session  --
#
# 19-May-2014   jdw refactor to better reuse new stored alignment data -
# 22-May-2014   jdw add output conflict list
#  3-Jun-2014   jdw fix self-reference issue
# 16-Oct-2014   jdw add explicit flag for sequence heterogeneity
# 22-Oct-2014   jdw add combination microheterogeneity/modified residue
# 01-Feb-2015   jdw add host org variant
# 01-Feb-2016   jdw change file handling for model merge -
# 07-Sep-2017   zf  add __cifCatAnnoDbRef() to export annotator's editing db reference information
# 31-Oct-2018   zf  split orignial SequenceDataExport into SequenceDataExport_v2 & AlignmentExport
##
"""
Export sequence and alignment details to production RCSB pipeline.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.07"

import os
import sys
import traceback

from mmcif.api.DataCategory import DataCategory
from mmcif.api.PdbxContainers import DataContainer
from mmcif.io.PdbxWriter import PdbxWriter

#
from wwpdb.apps.seqmodule.align.AlignmentExport import AlignmentExport
from wwpdb.apps.seqmodule.io.SequenceDataStore import SequenceDataStore
from wwpdb.apps.seqmodule.util.SequenceReferenceData import SequenceReferenceData
from wwpdb.apps.seqmodule.util.SequenceLabel import SequenceLabel
from wwpdb.apps.seqmodule.util.Autodict import Autodict

#
from wwpdb.io.misc.FormatOut import FormatOut
from wwpdb.io.locator.PathInfo import PathInfo
from wwpdb.utils.dp.RcsbDpUtility import RcsbDpUtility


class SequenceDataExport(object):
    """
    This class encapsulates all of the data export operations
    of sequence and other data to the RCSB/WF data pipeline.

    Storage model - exported data is loaded from the sequence data store
                    where it is managed by the SequenceDataStore() class.
    """

    def __init__(self, reqObj=None, exportList=None, verbose=False, log=sys.stderr):
        self.__verbose = verbose
        self.__reqObj = reqObj
        self.__lfh = log
        if exportList is None:
            exportList = []
        self.__debug = False
        self.__sessionObj = None
        self.__sessionPath = "."
        #
        # This is the setting for all selected sequences to be exported -
        #
        self.__summarySeqSelectList = exportList
        self.__selfRefList = []
        self.__setup()

    def __setup(self):
        try:
            self.__sessionObj = self.__reqObj.getSessionObj()
            self.__sessionPath = self.__sessionObj.getPath()

            self.__siteId = self.__reqObj.getValue("WWPDB_SITE_ID")
            self.__identifier = self.__reqObj.getValue("identifier")

            self.__srd = SequenceReferenceData(verbose=self.__verbose, log=self.__lfh)

            self.__pI = PathInfo(siteId=self.__siteId, sessionPath=self.__sessionPath, verbose=self.__verbose, log=self.__lfh)
            self.__exportFilePathSession = self.__pI.getSequenceAssignmentFilePath(self.__identifier, fileSource="session", versionId="next")
            if self.__verbose:
                self.__lfh.write("+SequenceDataExport.__setup() session assignment file path %s\n" % self.__exportFilePathSession)
            #
            self.__sds = SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            if self.__debug:
                self.__sds.dump(self.__lfh)
            #
            if not self.__summarySeqSelectList:
                self.__summarySeqSelectList = self.__sds.getSelectedIds()
            #
            self.__selfReferenceEntityList = []
            self.__trueSelfReferenceEntityList = []
            self.__annoRefDbL = []
            self.__annoSeqAuthRefMap = {}
            self.__I = self.__makeExportIndex(self.__summarySeqSelectList)
            self.__makeSelfReferenceIndex(self.__summarySeqSelectList)
            if self.__verbose:
                self.__lfh.write("+SequenceDataExport.__setup() incoming selected id list %r\n" % (self.__summarySeqSelectList))
                self.__lfh.write("+SequenceDataExport.__setup() self reference list %r\n" % self.__selfReferenceEntityList)
            #
        except:  # noqa: E722 pylint: disable=bare-except
            self.__lfh.write("+SequenceDataExport.__setup() failed for entry id %s\n" % (self.__identifier))
            traceback.print_exc(file=self.__lfh)
        #

    def exportAssignments(self):
        """Export assign file and associated updated model data file --"""
        numConflicts = 0
        conflictList = []
        try:
            ok, numConflicts, conflictList, warningMsg = self.__exportSeqMapping()
            return ok, numConflicts, conflictList, warningMsg
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+SequenceDataExport.exportAssignments() failed for entry id %s\n" % (self.__identifier))
                traceback.print_exc(file=self.__lfh)
            #
        #
        return False, numConflicts, conflictList, ""

    def getAllEntityIdList(self):
        """Interface to get all entity Id list from SequenceDataStore"""
        return self.__sds.getGroupIds()

    def applyAssignmentsToModel(self):
        """Create an updated model file in the current session --"""
        try:
            #  ----------------
            # Changes will always be applied to the earliest model in the current session  --
            #
            for v in ["1", "2", "3", "4"]:
                pdbxPath = self.__pI.getModelPdbxFilePath(self.__identifier, fileSource="session", versionId=v)
                self.__lfh.write("\n+SequenceDataExport.applyAssignmentsToModel() verifying starting model target path %s\n" % pdbxPath)
                if os.access(pdbxPath, os.R_OK):
                    break
                #
            #
            pdbxPathNext = self.__pI.getModelPdbxFilePath(self.__identifier, fileSource="session", versionId="next")
            logPath = os.path.join(self.__sessionPath, "annot-seqmod-merge.log")
            if os.access(logPath, os.R_OK):
                os.remove(logPath)
            #
            self.__lfh.write("\n+SequenceDataExport.applyAssignmentsToModel() starting model target path %s\n" % pdbxPath)
            self.__lfh.write("+SequenceDataExport.applyAssignmentsToModel() model output path %s\n" % pdbxPathNext)
            self.__lfh.write("+SequenceDataExport.applyAssignmentsToModel() assignment file path %s\n" % self.__exportFilePathSession)

            dp = RcsbDpUtility(tmpPath=self.__sessionPath, siteId=self.__siteId, verbose=self.__verbose, log=self.__lfh)
            dp.setDebugMode()
            dp.imp(pdbxPath)
            dp.addInput(name="seqmod_assign_file_path", value=self.__exportFilePathSession, type="file")
            dp.op("annot-merge-sequence-data")
            dp.exp(pdbxPathNext)
            dp.expLog(logPath)
            # dp.cleanup()
            if not os.access(pdbxPathNext, os.R_OK):
                return False
            #
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+SequenceDataExport.applyAssignmentsToModel() model update failed for assignment file  %s\n" % self.__exportFilePathSession)
                traceback.print_exc(file=self.__lfh)
            return False
        #
        return True

    def printIt(self):
        fOut = FormatOut()
        fOut.autoFormat("Summary Selection List", self.__summarySeqSelectList, 3, 3)
        fOut.autoFormat("Self Reference List", self.__selfReferenceEntityList, 3, 3)
        fOut.autoFormat("Selection Index", self.__I, 3, 3)
        fOut.writeStream(self.__lfh)
        fOut.clear()

    def __makeExportIndex(self, selectIdList):
        """Return an index dictionary for the input sequence id list -

        The sequence labels  have the format:

        seqType + "_" + seqInstId + "_" + seqPartId + "_" + seqAltId  + "_" + seqVersion

        or, for the special case of a self reference --

        "selfref"  "_"  entityId  "_"  partId
        """
        if self.__verbose:
            self.__lfh.write("+SequenceDataExport.__makeExportIndex() selected sequence id list %r\n" % selectIdList)
        #
        I = Autodict()  # noqa: E741
        try:
            for seqId in selectIdList:
                tup = seqId.strip().split("_")
                # Skip self reference ids
                if str(tup[0]).strip() == "selfref":
                    self.__selfRefList.append(seqId)
                    continue
                #
                I[tup[0]][tup[1]][int(tup[2])] = (int(tup[3]), int(tup[4]))
            #
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+SequenceDataExport.__makeExportIndex() failed \n")
                traceback.print_exc(file=self.__lfh)
            #
        #
        return I

    def __makeSelfReferenceIndex(self, selectIdList):
        """Return the list self referenced entities for the input sequence id list -

        The sequence labels  have the format:

        seqType + "_" + seqInstId + "_" + seqPartId + "_" + seqAltId  + "_" + seqVersion

        or, for the special case of a self reference --

        "selfref"  "_"  entityId  "_"  partId
        """
        if self.__verbose:
            self.__lfh.write("+SequenceDataExport.__makeSelfReferenceIndex() selected sequence id list %r\n" % selectIdList)
        #
        try:
            for seqId in selectIdList:
                tup = seqId.strip().split("_")
                if str(tup[0]).strip() == "selfref":
                    if self.__verbose:
                        self.__lfh.write("+SequenceDataExport.__makeSelfReferenceIndex() - assigning reference for sequence entity %s part %s (%s)\n" % (tup[1], tup[2], tup[0]))
                    #
                    seqId = str(tup[1]).strip()
                    partId = str(tup[2]).strip()
                    self.__selfReferenceEntityList.append((seqId, partId))
                    foundAnnInputRefInfo = False
                    if (seqId in self.__I["auth"]) and (int(partId) in self.__I["auth"][seqId]):
                        altId = self.__I["auth"][seqId][int(partId)][0]
                        verId = self.__I["auth"][seqId][int(partId)][1]
                        seqAuthFD = self.__sds.getFeature(seqId=seqId, seqType="auth", partId=int(partId), altId=altId, version=verId)
                        row = []
                        row.append(len(self.__annoRefDbL) + 1)
                        row.append(seqId)
                        for item in ("ANNO_EDIT_DB_NAME", "ANNO_EDIT_DB_CODE", "ANNO_EDIT_DB_ACCESSION", "ANNO_EDIT_DB_ALIGN_BEGIN", "ANNO_EDIT_DB_ALIGN_END"):
                            if (item in seqAuthFD) and seqAuthFD[item]:
                                row.append(seqAuthFD[item])
                            else:
                                break
                            #
                        #
                        if len(row) == 7:
                            seqAuthRefMap = []
                            if (seqId in self.__annoSeqAuthRefMap) and self.__annoSeqAuthRefMap[seqId]:
                                seqAuthRefMap = self.__annoSeqAuthRefMap[seqId]
                            else:
                                seqAuthIdx = self.__sds.getSequence(seqId=seqId, seqType="auth", partId=int(partId), altId=altId, version=verId)
                                idx = 1
                                for seqAuth in seqAuthIdx:
                                    seqAuthRefMap.append([".", seqId, seqAuth[0], str(idx), seqId, ".", ".", "."])
                                    idx += 1
                                #
                            #
                            if seqAuthRefMap:
                                authFObj = self.__sds.getFeatureObj(seqId, "auth", partId=int(partId), altId=altId, version=verId)
                                _authPartId, authSeqBeg, authSeqEnd, _authSeqPartType = authFObj.getAuthPartDetails()
                                partBeg = int(authSeqBeg)
                                partEnd = int(authSeqEnd)
                                foundMap = False
                                startNum = int(row[5])
                                for seqAuth in seqAuthRefMap:
                                    idx = int(seqAuth[3])
                                    if (idx >= partBeg) and (idx <= partEnd):
                                        if idx > partBeg:
                                            startNum += 1
                                        #
                                        seqAuth[5] = seqAuth[2]
                                        seqAuth[6] = str(startNum)
                                        seqAuth[7] = str(partId)
                                        foundMap = True
                                    #
                                #
                                if foundMap:
                                    row.append(partId)
                                    self.__annoRefDbL.append(row)
                                    foundAnnInputRefInfo = True
                                    self.__annoSeqAuthRefMap[seqId] = seqAuthRefMap
                                #
                            #
                        #
                    #
                    if not foundAnnInputRefInfo:
                        self.__trueSelfReferenceEntityList.append((seqId, partId))
                    #
                #
            #
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+SequenceDataExport.__makeSelfReferenceIndex() failed \n")
                traceback.print_exc(file=self.__lfh)
            #
        #

    def __exportSeqMapping(self):
        """Export worker function -"""
        #
        # Get entry ids
        depDataSetId = self.__sds.getEntryDetail("DEPOSITION_DATA_SET_ID")
        # pdbId = self.__sds.getEntryDetail("PDB_ID")
        #
        # Get the entity group list -
        gIdList = self.__sds.getGroupIds()
        if len(gIdList) < 1:
            if self.__verbose:
                self.__lfh.write("+SequenceDataExport.__exportSeqMapping() for data set %s group list is empty\n" % depDataSetId)
            #
        #
        if self.__verbose:
            self.__lfh.write("+SequenceDataExport.__exportSeqMapping() number of groups %d depDataSetId %s\n" % (len(gIdList), depDataSetId))
        #

        seqLabel = SequenceLabel()
        refFeatureD = {}
        rptRefL = []
        rptXyzL = []
        rptCommentL = []
        rptCommentModL = []
        rptDeleteL = []
        allRefSeqIdxD = {}
        numConflicts = 0
        conflictList = []
        warningMsg = ""
        natureSourceTaxIds = {}
        #
        polyAlaCaseList = []
        for gId in gIdList:
            #
            # Get seqIds in group -
            #
            seqIdList = self.__sds.getGroup(gId)
            if len(seqIdList) < 1:
                if self.__verbose:
                    self.__lfh.write("+SequenceDataExport.__exportSeqMapping() skipping entity group %s which has no instances\n" % (gId))
                #
                continue
            #
            if self.__verbose:
                self.__lfh.write("+SequenceDataExport.__exportSeqMapping() group %s sequence list %r\n" % (gId, seqIdList))
            #
            # Author sequence -
            #
            sourceInfo = ""
            idAuthSeq = ""
            first = True
            if gId in self.__I["auth"]:
                if self.__verbose:
                    self.__lfh.write("+SequenceDataExport.__exportSeqMapping() entity group %s author sequence s list %r\n" % (gId, list(self.__I["auth"][gId].items())))
                #
                for partId, (altId, ver) in self.__I["auth"][gId].items():
                    fD = self.__sds.getFeature(gId, seqType="auth", partId=partId, altId=altId, version=ver)
                    if ("SOURCE_METHOD" in fD) and fD["SOURCE_METHOD"].upper() == "NAT":
                        sourceInfo = "NAT"
                        if ("SOURCE_TAXID" in fD) and fD["SOURCE_TAXID"]:
                            if fD["SOURCE_TAXID"] in natureSourceTaxIds:
                                if gId not in natureSourceTaxIds[fD["SOURCE_TAXID"]]:
                                    natureSourceTaxIds[fD["SOURCE_TAXID"]].append(gId)
                                #
                            else:
                                natureSourceTaxIds[fD["SOURCE_TAXID"]] = [gId]
                            #
                        #
                    #
                    seqLabel.set(seqType="auth", seqInstId=gId, seqPartId=partId, seqAltId=altId, seqVersion=ver)
                    idAuthSeq = seqLabel.pack()
                    if self.__verbose:
                        self.__lfh.write("+SequenceDataExport.__exportSeqMapping() from entity group %s exporting author sequence %s\n" % (gId, idAuthSeq))
                    #
                    if first and ("POLYMER_TYPE" in fD) and (fD["POLYMER_TYPE"] == "AA"):
                        seqAuth = self.__sds.getSequence(seqId=gId, seqType="auth", partId=partId, altId=altId, version=ver)
                        polyALA_assignment = self.__checkPolyAlaAssignment(seqAuth)
                        if polyALA_assignment > 0:
                            polyAlaCaseList.append((gId, polyALA_assignment))
                        #
                    #
                    first = False
                #
            #
            if not idAuthSeq:
                # report missing author sequence
                if self.__verbose:
                    self.__lfh.write("+SequenceDataExport.__exportSeqMapping() missing author sequence from entity group %s\n" % gId)
                #
                continue
            #
            # Reference sequences -
            #
            idListRef = []
            if gId in self.__I["ref"]:
                for partId, (altId, ver) in self.__I["ref"][gId].items():
                    seqRefIdx = self.__sds.getSequence(seqId=gId, seqType="ref", partId=partId, altId=altId, version=ver)
                    seqRefFD = self.__sds.getFeature(seqId=gId, seqType="ref", partId=partId, altId=altId, version=ver)
                    seqLabel.set(seqType="ref", seqInstId=gId, seqPartId=partId, seqAltId=altId, seqVersion=ver)
                    idRefSeq = seqLabel.pack()
                    idListRef.append(idRefSeq)
                    #
                    refFeatureD[(gId, partId)] = (idRefSeq, altId, seqRefFD, partId)
                    allRefSeqIdxD[idRefSeq] = seqRefIdx
                    #
                    if self.__verbose:
                        self.__lfh.write("+SequenceDataExport.__exportSeqMapping() from entity group %s part %s exporting reference sequence %s\n" % (gId, partId, idRefSeq))
                    #
                #
            else:
                # report missing reference sequence
                if self.__verbose:
                    self.__lfh.write("+SequenceDataExport.__exportSeqMapping() missing reference sequence from group %s\n" % gId)
                #
            #
            # Coordinate sequences -
            #
            idListXyz = []
            for seqId in seqIdList:
                if seqId in self.__I["xyz"]:
                    for partId, (altId, ver) in self.__I["xyz"][seqId].items():
                        seqLabel.set(seqType="xyz", seqInstId=seqId, seqPartId=partId, seqAltId=altId, seqVersion=ver)
                        idXyzSeq = seqLabel.pack()
                        idListXyz.append(idXyzSeq)
                    #
                else:
                    # report missing coordinate sequence
                    pass
                #
            #
            if (len(idListRef) < 1) and (len(idListXyz) < 1):
                continue
            #
            alignExport = AlignmentExport(reqObj=self.__reqObj, entityId=gId, pathInfo=self.__pI, seqDataStore=self.__sds, verbose=self.__verbose, log=self.__lfh)
            localRptRefL, localRptCommentL, localRptCommentModL, localRptDeleteL, localRptXyzL, message, numC = alignExport.doExport(
                idAuthSeq, idListRef, idListXyz, self.__selfRefList, sourceInfo
            )
            if self.__verbose:
                errMsg = alignExport.getErrorMessage()
                if errMsg:
                    self.__lfh.write("\nError message for entity '%s'\n%s\n" % (gId, errMsg))
                #
            #
            if len(localRptRefL) > 0:
                rptRefL.extend(localRptRefL)
            #
            if len(localRptCommentL) > 0:
                rptCommentL.extend(localRptCommentL)
            #
            if len(localRptCommentModL) > 0:
                rptCommentModL.extend(localRptCommentModL)
            #
            if len(localRptDeleteL) > 0:
                rptDeleteL.extend(localRptDeleteL)
            #
            if len(localRptXyzL) > 0:
                rptXyzL.extend(localRptXyzL)
            #
            if message:
                warningMsg += message
            #
            if numC > 0:
                numConflicts += numC
                conflictList.append((gId, numC))
            #
        #  END of iteration over  --- gId  ----
        #
        if len(natureSourceTaxIds) > 1:
            warningMsg += "Entry contains multiple natural sources:<br />\n"
            for k, v in natureSourceTaxIds.items():
                if len(v) > 1:
                    warningMsg += "Entities '" + "', '".join(v) + "' have "
                else:
                    warningMsg += "Entity '" + "', '".join(v) + "' has "
                #
                warningMsg += "source taxonomy Id '" + k + "'.<br />\n"
            #
        #
        if polyAlaCaseList:
            for polyAlaCaseTupL in polyAlaCaseList:
                if warningMsg:
                    warningMsg += "<br />\n"
                #
                warningMsg += "Entity [" + polyAlaCaseTupL[0] + "]"
                if polyAlaCaseTupL[1] == 1:
                    warningMsg += " is composed only of poly-ALA.<br />\n"
                else:
                    warningMsg += " has stretches of poly-ALA.<br />\n"
                #
            #
        # ref FeatureD ---
        rptRefDbL = self.__refDdReport(refFeatureD)
        # rptDeleteL = self.__deletionReport(refFeatureD, allRefSeqIdxD)

        try:
            # Make a local copy of the mapping file in the session directory and then copy the file as needed.
            #
            blockList = []
            ofh = open(self.__exportFilePathSession, "w")
            dataBlock = DataContainer(depDataSetId)
            if len(self.__trueSelfReferenceEntityList) > 0:
                dataBlock.append(self.__cifCatSelfReference(self.__trueSelfReferenceEntityList))
            #
            if len(rptRefDbL) > 0:
                dataBlock.append(
                    self.__writeGeneralCifCategory(
                        "pdbx_seqtool_db_ref", ("ordinal", "entity_id", "db_name", "db_code", "db_accession", "db_isoform", "match_begin", "match_end", "entity_part_id"), rptRefDbL
                    )
                )
            #
            if len(self.__annoRefDbL) > 0:
                dataBlock.append(
                    self.__writeGeneralCifCategory(
                        "pdbx_seqtool_anno_db_ref", ("ordinal", "entity_id", "db_name", "db_code", "db_accession", "match_begin", "match_end", "entity_part_id"), self.__annoRefDbL
                    )
                )
                #
                for k, seqAuthRefMap in self.__annoSeqAuthRefMap.items():
                    for seqAuth in seqAuthRefMap:
                        seqAuth[0] = str(len(rptRefL) + 1)
                        rptRefL.append(tuple(seqAuth))
                    #
                #
            #
            dataBlock.append(self.__getEntitySourceCategory())
            #
            if len(rptCommentL) > 0:
                dataBlock.append(
                    self.__writeGeneralCifCategory("pdbx_seqtool_mapping_comment", ("ordinal", "entity_id", "entity_mon_id", "entity_seq_num", "comment"), rptCommentL)
                )
            #
            if len(rptCommentModL) > 0:
                dataBlock.append(
                    self.__writeGeneralCifCategory(
                        "pdbx_seqtool_mapping_modification_comment", ("ordinal", "entity_id", "entity_mon_id", "entity_seq_num", "comment"), rptCommentModL
                    )
                )
            #
            if len(rptDeleteL) > 0:
                dataBlock.append(
                    self.__writeGeneralCifCategory(
                        "pdbx_seqtool_ref_deletions", ("ordinal", "entity_id", "part_id", "db_name", "db_accession", "db_isoform", "ref_mon_id", "ref_mon_num"), rptDeleteL
                    )
                )
            #
            if len(rptRefL) > 0:
                dataBlock.append(
                    self.__writeGeneralCifCategory(
                        "pdbx_seqtool_mapping_ref",
                        ("ordinal", "entity_id", "entity_mon_id", "entity_seq_num", "auth_asym_id", "ref_mon_id", "ref_mon_num", "entity_part_id"),
                        rptRefL,
                    )
                )
            #
            if len(rptXyzL) > 0:
                dataBlock.append(
                    self.__writeGeneralCifCategory(
                        "pdbx_seqtool_mapping_xyz",
                        ("ordinal", "entity_id", "entity_mon_id", "entity_seq_num", "auth_asym_id", "pdb_mon_id", "pdb_mon_num", "pdb_ins_code", "hetero_flag", "org_mon_id"),
                        rptXyzL,
                    )
                )
            #
            blockList.append(dataBlock)
            pdbxW = PdbxWriter(ofh)
            pdbxW.write(blockList)
            ofh.close()
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+SequenceDataExport.__exportSeqMapping() saving export mapping file failed %s\n" % self.__exportFilePathSession)
                traceback.print_exc(file=self.__lfh)
            return False, numConflicts, conflictList, warningMsg

        return True, numConflicts, conflictList, warningMsg

    def __checkPolyAlaAssignment(self, seqAuth):
        """Check if the sequence contains 10 or more consecutive ALA residues"""
        has_consecutive_ALA = False
        count = 0
        for seqTupL in seqAuth:
            if seqTupL[0] == "ALA":
                count += 1
            else:
                if count > 9:
                    has_consecutive_ALA = True
                #
                count = 0
            #
        #
        if count > 9:
            has_consecutive_ALA = True
        #
        if count == len(seqAuth):
            return 1
        elif has_consecutive_ALA:
            return 2
        #
        return 0

    def __refDdReport(self, refFeatureD):
        """Return a list of content rows for reference sequence database accession report --"""
        rptRefDbL = []
        irow = 1
        for (gId, partId), tup in refFeatureD.items():
            if (str(gId), str(partId)) in self.__selfReferenceEntityList:
                continue
            #
            # sId = tup[0]
            fD = tup[2]
            rptRefDbL.append(
                [
                    irow,
                    gId,
                    self.__srd.convertDbNameToResource(fD["DB_NAME"]),
                    fD["DB_CODE"],
                    fD["DB_ACCESSION"],
                    fD["DB_ISOFORM"],
                    fD["REF_MATCH_BEGIN"],
                    fD["REF_MATCH_END"],
                    partId,
                ]
            )
            irow += 1
        #
        return rptRefDbL

    # def __deletionReport(self, refFeatureD, allRefSeqIdxD):
    #     """Return a list of content rows for the sequence deletion report."""
    #     rptDeleteL = []
    #     idel = 1
    #     for (gId, partId), tup in refFeatureD.items():
    #         if (str(gId), str(partId)) in self.__selfReferenceEntityList:
    #             continue
    #         #
    #         sId = tup[0]
    #         fD = tup[2]
    #         if sId in allRefSeqIdxD:
    #             if self.__verbose:
    #                 self.__lfh.write("+SequenceDataExport.__exportSeqMapping() searching deletions in seqId %s\n" % sId)
    #             #
    #             refSeqIdx = allRefSeqIdxD[sId]
    #             for sTup in refSeqIdx:
    #                 if sTup[2] in ["Deletion", "deletion"]:
    #                     if self.__verbose:
    #                         self.__lfh.write("+SequenceDataExport.__exportSeqMapping() position %s comment %s\n" % (sTup[1], sTup[2]))
    #                     #
    #                     rptDeleteL.append([idel, gId, partId, self.__srd.convertDbNameToResource(fD["DB_NAME"]), fD["DB_ACCESSION"], fD["DB_ISOFORM"], sTup[0], sTup[1]])
    #                     idel += 1
    #                 #
    #             #
    #         else:
    #             self.__lfh.write("+SequenceDataExport.__exportSeqMapping() delete scan missing reference sequence for seqId %s keys() %r\n" % (sId, allRefSeqIdxD.keys()))
    #         #
    #     #
    #     return rptDeleteL

    def __cifCatSelfReference(self, entityList):
        """Update entities for self reference  ---"""
        aCat = DataCategory("pdbx_seqtool_self_ref")
        aCat.appendAttribute("entity_id")
        aCat.appendAttribute("entity_part_id")
        for eid, partId in entityList:
            aCat.append([eid, partId])
        return aCat

    def __writeGeneralCifCategory(self, categoryName, categoryItems, dataList):
        """Write pdbx_seqtool_db_ref, pdbx_seqtool_anno_db_ref, pdbx_seqtool_ref_deletions, pdbx_seqtool_mapping_ref, pdbx_seqtool_mapping_comment,
        pdbx_seqtool_mapping_modification_comment, pdbx_seqtool_mapping_xyz categories
        """
        if len(dataList) < 1:
            return None
        #
        aCat = DataCategory(categoryName)
        for itemName in categoryItems:
            aCat.appendAttribute(itemName)
        #
        rowCount = 0
        for row in dataList:
            rowCount += 1
            row[0] = str(rowCount)
            aCat.append(row)
        #
        return aCat

    def __getEntitySourceCategory(self):
        """Return category object of source details."""
        itemTupList = [
            ("ORDINAL", "ordinal"),
            ("ENTITY_ID", "entity_id"),
            ("AUTH_SEQ_PART_ID", "seq_part_id"),
            ("AUTH_SEQ_NUM_BEGIN", "seq_part_beg"),
            ("AUTH_SEQ_NUM_END", "seq_part_end"),
            ("AUTH_SEQ_PART_TYPE", "seq_part_type"),
            ("ENTITY_DESCRIPTION", "entity_description"),
            ("ENTITY_SYNONYMS", "entity_synonyms"),
            ("SOURCE_GENE_NAME", "gene_name"),
            ("SOURCE_TAXID", "taxonomy_id"),
            ("SOURCE_ORGANISM", "source_scientific_name"),
            ("SOURCE_STRAIN", "source_strain"),
            ("SOURCE_COMMON_NAME", "source_common_name"),
            ("SOURCE_VARIANT", "variant"),
            ("POLYMER_LINKING_TYPE", "entity_polymer_type"),
            ("ENTITY_ENZYME_CLASS", "entity_enzyme_class"),
            ("ENTITY_FRAGMENT_DETAILS", "entity_fragment_details"),
            ("ENTITY_MUTATION_DETAILS", "entity_mutation_details"),
            ("ENTITY_DETAILS", "entity_details"),
            ("SOURCE_METHOD", "source_method"),
            ("HOST_ORG_SOURCE", "host_org_scientific_name"),
            ("HOST_ORG_STRAIN", "host_org_strain"),
            ("HOST_ORG_TAXID", "host_org_taxonomy_id"),
            ("HOST_ORG_VECTOR", "host_org_vector"),
            ("HOST_ORG_VECTOR_TYPE", "host_org_vector_type"),
            ("HOST_ORG_PLASMID", "host_org_plasmid"),
            ("HOST_ORG_COMMON_NAME", "host_org_common_name"),
            ("HOST_ORG_CELL_LINE", "host_org_cell_line"),
            ("HOST_ORG_VARIANT", "host_org_variant"),
        ]

        if self.__verbose:
            self.__lfh.write("+SequenceDataExport.__getSourceDetails() starting\n")
            # self.printIt()
        #
        rowList = []
        entityIdList = self.__sds.getGroupIds()
        iRow = 1
        for entityId in entityIdList:
            seqIds = self.__sds.getGroup(entityId)
            if len(seqIds) < 1:
                if self.__verbose:
                    self.__lfh.write("+SequenceDataExport.__getSourceDetails() entityId %s is empty\n" % (entityId))
                #
                continue
            #
            seqId0 = entityId
            partIdList = self.__sds.getPartIds(seqId0, dataType="sequence", seqType="auth")
            #
            for partId in partIdList:
                if self.__verbose:
                    self.__lfh.write("+SequenceDataExport.__getSourceDetails() entityId %r partId %r  __I  %r\n" % (entityId, partId, self.__I["auth"][seqId0][partId]))
                # JDW
                # altId,versionId=self.__I["auth"][seqId0][partId]
                altId, versionId = self.__I["auth"][seqId0][1]
                fD = self.__sds.getFeature(seqId0, seqType="auth", partId=partId, altId=altId, version=versionId)
                d = {}
                for itemTup in itemTupList:
                    ky = itemTup[0]
                    if ky in fD:
                        d[ky] = fD[ky]
                    else:
                        d[ky] = ""
                    #
                #
                d["ORDINAL"] = iRow
                d["ENTITY_ID"] = entityId
                d["AUTH_SEQ_PART_ID"] = partId
                rowList.append(d)
                iRow += 1
                #
            #
        #
        aCat = DataCategory("pdbx_seqtool_entity_details")
        for itemTup in itemTupList:
            attribName = itemTup[1]
            aCat.appendAttribute(attribName)
        #
        for rD in rowList:
            vL = []
            for itemTup in itemTupList:
                ky = itemTup[0]
                vL.append(rD[ky])
            #
            aCat.append(vL)
        #
        return aCat
