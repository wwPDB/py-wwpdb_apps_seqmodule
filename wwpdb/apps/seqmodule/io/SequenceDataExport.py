##
# File:  SequenceDataExport.py
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
#  19-Jan-2014  jdw 'HOST_ORG_CELL_LINE'
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
##
"""
Export sequence and alignment details to production RCSB pipeline.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.07"

import sys
import os
import string
import shutil
import traceback

from wwpdb.apps.seqmodule.io.SequenceDataStore import SequenceDataStore
from wwpdb.apps.seqmodule.util.SequenceReferenceData import SequenceReferenceData
from wwpdb.apps.seqmodule.util.SequenceLabel import SequenceLabel, SequenceFeature
from wwpdb.apps.seqmodule.util.Autodict import Autodict
# new
from wwpdb.apps.seqmodule.align.MultiAlignPseudo import MultiAlignPseudo
#
from wwpdb.io.locator.PathInfo import PathInfo
from wwpdb.utils.dp.RcsbDpUtility import RcsbDpUtility
from wwpdb.io.misc.FormatOut import FormatOut
from wwpdb.utils.align.alignlib import PairwiseAlign
#
from mmcif.io.PdbxWriter import PdbxWriter
from mmcif.api.PdbxContainers import *
#
#
#


class SequenceDataExport(object):

    """
     This class encapsulates all of the data export operations
     of sequence and other data to the RCSB/WF data pipeline.

     Storage model - exported data is loaded from the sequence data store
                     where it is managed by the SequenceDataStore() class.

    """

    def __init__(self, reqObj=None, exportList=[], verbose=False, log=sys.stderr):
        self.__verbose = verbose
        self.__reqObj = reqObj
        self.__lfh = log
        self.__debug = False
        self.__sessionObj = None
        self.__sessionPath = '.'
        #
        # This is the setting for all selected sequences to be exported -
        #
        self.__summarySeqSelectList = exportList
        #
        # if (len(self.__summarySeqSelectList) == 0):
        #    self.__summarySeqSelectList=self.__reqObj.getSummarySelectList()
        #
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
            if (self.__verbose):
                self.__lfh.write("+SequenceDataExport.__setup() session assignment file path %s\n" % self.__exportFilePathSession)

            self.__sds = SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            if (self.__debug):
                self.__sds.dump(self.__lfh)
            #

            self.__selfReferenceEntityList = []
            self.__trueSelfReferenceEntityList = []
            self.__annoRefDbL = []
            self.__annoSeqAuthRefMap = {}
            self.__I = self.makeExportIndex(self.__summarySeqSelectList)
            self.makeSelfReferenceIndex(self.__summarySeqSelectList)
            if (self.__verbose):
                self.__lfh.write("+SequenceDataExport.__setup() incoming selected id list %r\n" % (self.__summarySeqSelectList))
                self.__lfh.write("+SequenceDataExport.__setup() self reference list %r\n" % self.__selfReferenceEntityList)
            #
        except:
            self.__lfh.write("+SequenceDataExport.__setup() failed for entry id %s\n" % (self.__identifier))
            traceback.print_exc(file=self.__lfh)

    def exportAssignments(self):
        """ Export assign file and associated updated model data file --
        """
        numConflicts = 0
        conflictList = []
        try:
            ok, numConflicts, conflictList, warningMsg = self.__exportSeqMapping()
            return ok, numConflicts, conflictList, warningMsg
        except:
            if (self.__verbose):
                self.__lfh.write("+SequenceDataExport.exportAssignments() failed for entry id %s\n" % (self.__identifier))
                traceback.print_exc(file=self.__lfh)
        return False, numConflicts, conflictList, ''

    def applyAssignmentsToModel(self):
        """  Create an updated model file in the current session --
        """
        try:
            #  ----------------
            # Changes will always be applied to the earliest model in the current session  --
            #
            for v in ["1", "2", "3", "4"]:
                pdbxPath = self.__pI.getModelPdbxFilePath(self.__identifier, fileSource="session", versionId=v)
                self.__lfh.write("\n+SequenceDataExport.applyAssignmentsToModel() verifying starting model target path %s\n" % pdbxPath)
                if os.access(pdbxPath, os.R_OK):
                    break
            pdbxPathNext = self.__pI.getModelPdbxFilePath(self.__identifier, fileSource="session", versionId="next")
            logPath = os.path.join(self.__sessionPath, "annot-seqmod-merge.log")

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
        except:
            if (self.__verbose):
                self.__lfh.write("+SequenceDataExport.applyAssignmentsToModel() model update failed for assignment file  %s\n"
                                 % self.__exportFilePathSession)
                traceback.print_exc(file=self.__lfh)
            return False
        #

        return True

    def __updateAssignmentsArchive(self, identifier, exportFilePathSession):
        """  Copy the assignment file from session to archive storage --
        """
        exportFilePathArchive = self.__pI.getSequenceAssignmentFilePath(identifier, fileSource="archive", versionId="next")
        try:
            if (self.__verbose):
                self.__lfh.write("+SequenceDataExport.__updateArchive() copying assignment file to archive %s\n" % exportFilePathArchive)
            shutil.copyfile(exportFilePathSession, exportFilePathArchive)
            return True
        except:
            if (self.__verbose):
                self.__lfh.write("+SequenceDataExport.__updateArchive() copying mapping file failed for: %s\n" % exportFilePathArchive)
        return False

    def makeExportIndex(self, selectIdList):
        """ Return an index dictionary for the input sequence id list -

            The sequence labels  have the format:

            seqType + '_' + seqInstId + '_' + seqPartId + '_' + seqAltId  + '_' + seqVersion

            or, for the special case of a self reference --

            'selfref'  '_'  entityId  '_'  partId
        """
        if (self.__verbose):
            self.__lfh.write("+SequenceDataExport.makeExportIndex() selected sequence id list %r\n" % selectIdList)
        I = Autodict()
        try:
            for seqId in selectIdList:
                tup = seqId.strip().split('_')
                # Skip self reference ids
                if str(tup[0]).strip() == 'selfref':
                    continue
                I[tup[0]][tup[1]][int(tup[2])] = (int(tup[3]), int(tup[4]))
                #
        except:
            if (self.__verbose):
                self.__lfh.write("+SequenceDataExport.makeExportIndex() failed \n")
                traceback.print_exc(file=self.__lfh)

        return I

    def makeSelfReferenceIndex(self, selectIdList):
        """ Return the list self referenced entities for the input sequence id list -

            The sequence labels  have the format:

            seqType + '_' + seqInstId + '_' + seqPartId + '_' + seqAltId  + '_' + seqVersion

            or, for the special case of a self reference --

            'selfref'  '_'  entityId  '_'  partId
        """

        if (self.__verbose):
            self.__lfh.write("+SequenceDataExport.makeSelfReferenceIndex() selected sequence id list %r\n" % selectIdList)

        try:
            for seqId in selectIdList:
                tup = seqId.strip().split('_')
                if str(tup[0]).strip() == 'selfref':
                    if (self.__verbose):
                        self.__lfh.write("+SequenceDataExport.makeSelfReferenceIndex() - assigning reference for sequence entity %s part %s (%s)\n" %
                                         (tup[1], tup[2], tup[0]))
                    seqId = str(tup[1]).strip()
                    partId = str(tup[2]).strip()
                    self.__selfReferenceEntityList.append((seqId, partId))
                    foundAnnInputRefInfo = False
                    if self.__I['auth'].has_key(seqId) and self.__I['auth'][seqId].has_key(int(partId)):
                        altId = self.__I['auth'][seqId][int(partId)][0]
                        verId = self.__I['auth'][seqId][int(partId)][1]
                        seqAuthFD = self.__sds.getFeature(seqId=seqId, seqType='auth', partId=int(partId), altId=altId, version=verId)
                        row = []
                        row.append(len(self.__annoRefDbL) + 1)
                        row.append(seqId)
                        for item in ( 'ANNO_EDIT_DB_NAME', 'ANNO_EDIT_DB_CODE', 'ANNO_EDIT_DB_ACCESSION', 'ANNO_EDIT_DB_ALIGN_BEGIN', 'ANNO_EDIT_DB_ALIGN_END' ):
                            if seqAuthFD.has_key(item) and seqAuthFD[item]:
                                row.append(seqAuthFD[item])
                            else:
                                break
                            #
                        #
                        if len(row) == 7:
                            seqAuthRefMap = []
                            if self.__annoSeqAuthRefMap.has_key(seqId) and self.__annoSeqAuthRefMap[seqId]:
                                seqAuthRefMap = self.__annoSeqAuthRefMap[seqId]
                            else:
                                seqAuthIdx = self.__sds.getSequence(seqId=seqId,seqType='auth',partId=int(partId),altId=altId,version=verId)
                                idx = 1
                                for seqAuth in seqAuthIdx:
                                    seqAuthRefMap.append([ '.', seqId, seqAuth[0], str(idx), seqId, '.', '.', '.' ])
                                    idx += 1
                                #
                            #
                            if seqAuthRefMap:
                                authFObj = self.__sds.getFeatureObj(seqId, 'auth', partId=int(partId),altId=altId,version=verId)
                                authPartId, authSeqBeg, authSeqEnd, authSeqPartType = authFObj.getAuthPartDetails()
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
        except:
            if (self.__verbose):
                self.__lfh.write("+SequenceDataExport.makeSelfReferenceIndex() failed \n")
                traceback.print_exc(file=self.__lfh)
            #
        #

    def getFilePath(self):
        return self.__exportFilePathSession

    def printIt(self):
        fOut = FormatOut()
        fOut.autoFormat("Summary Selection List", self.__summarySeqSelectList, 3, 3)
        fOut.autoFormat("Self Reference List", self.__selfReferenceEntityList, 3, 3)
        fOut.autoFormat("Selection Index", self.__I, 3, 3)
        fOut.writeStream(self.__lfh)
        fOut.clear()

    def __exportSeqMapping(self):
        """  Export worker function -
        """
        #
        # Get entry ids
        depDataSetId = self.__sds.getEntryDetail('DEPOSITION_DATA_SET_ID')
        pdbId = self.__sds.getEntryDetail('PDB_ID')
        #
        # Get the entity group list -
        #
        gIdList = self.__sds.getGroupIds()
        if len(gIdList) < 1:
            if (self.__verbose):
                self.__lfh.write("+SequenceDataExport.__exportSeqMapping() for data set %s group list is empty\n" % depDataSetId)
        #
        if (self.__verbose):
            self.__lfh.write("+SequenceDataExport.__exportSeqMapping() number of groups %d depDataSetId %s\n" % (len(gIdList), depDataSetId))

        seqLabel = SequenceLabel()
        seqLenD = {}
        refFeatureD = {}
        rptRefL = []
        rptXyzL = []
        rptCommentL = []
        rptCommentModL = []
        allRefSeqIdxD = {}
        numConflicts = 0
        conflictList = []
        warningMsg = ''
        natureSourceTaxIds = {}
        #
        for gId in gIdList:
            #
            # Get seqIds in group -
            #
            seqIdList = self.__sds.getGroup(gId)
            flagRefOk = True
            if len(seqIdList) < 1:
                if (self.__verbose):
                    self.__lfh.write("+SequenceDataExport.__exportSeqMapping() skipping entity group %s which has no instances\n" % (gId))
                continue
            #
            if (self.__verbose):
                self.__lfh.write("+SequenceDataExport.__exportSeqMapping() group %s sequence list %r\n" % (gId, seqIdList))
            #
            # JDW CHANGE
            # seqId0=seqIdList[0]
            seqId0 = gId

            #
            # Author sequence -
            #
            sourceInfo = ''
            if seqId0 in self.__I['auth']:
                if (self.__verbose):
                    self.__lfh.write("+SequenceDataExport.__exportSeqMapping() entity group %s author sequence %s list %r\n" % (gId, seqId0, self.__I['auth'][seqId0].items()))
                for partId, (altId, ver) in self.__I['auth'][seqId0].items():
                    fD = self.__sds.getFeature(seqId0, seqType="auth", partId=partId, altId=altId, version=ver)
                    if fD.has_key('SOURCE_METHOD') and fD['SOURCE_METHOD'].upper() == 'NAT':
                        sourceInfo = 'NAT'
                        if fD.has_key('SOURCE_TAXID') and fD['SOURCE_TAXID']:
                            if natureSourceTaxIds.has_key(fD['SOURCE_TAXID']):
                                if gId not in natureSourceTaxIds[fD['SOURCE_TAXID']]:
                                    natureSourceTaxIds[fD['SOURCE_TAXID']].append(gId)
                                #
                            else:
                                natureSourceTaxIds[fD['SOURCE_TAXID']] = [ gId ]
                            #
                        #
                    #
                    # seqAuthIdx=self.__sds.getSequence(seqId=seqId0,seqType='auth',partId=partId,altId=altId,version=ver)
                    seqLabel.set(seqType='auth', seqInstId=seqId0, seqPartId=partId, seqAltId=altId, seqVersion=ver)
                    idAuthSeq = seqLabel.pack()
                    # seqLenD[idAuthSeq]=len(seqAuthIdx)
                    if (self.__verbose):
                        self.__lfh.write("+SequenceDataExport.__exportSeqMapping() from entity group %s exporting author sequence %s\n" %
                                         (gId, idAuthSeq))
            else:
                # report missing author sequence
                if (self.__verbose):
                    self.__lfh.write("+SequenceDataExport.__exportSeqMapping() missing author sequence from entity group %s seqId %s\n" %
                                     (gId, seqId0))
                continue

            #
            # Reference sequences -
            #
            # sfObj=SequenceFeature()
            # refSeqIdxD={}
            idListRef = []
            if seqId0 in self.__I['ref']:
                for partId, (altId, ver) in self.__I['ref'][seqId0].items():

                    seqRefIdx = self.__sds.getSequence(seqId=seqId0, seqType='ref', partId=partId, altId=altId, version=ver)
                    seqRefFD = self.__sds.getFeature(seqId=seqId0, seqType='ref', partId=partId, altId=altId, version=ver)
                    seqLabel.set(seqType='ref', seqInstId=seqId0, seqPartId=partId, seqAltId=altId, seqVersion=ver)
                    idRefSeq = seqLabel.pack()
                    idListRef.append(idRefSeq)
                    #
                    # seqLenD[idRefSeq]=len(seqRefIdx)
                    # sfObj.set(seqRefFD)
                    # orgName,strainName=sfObj.decodeUniProtSourceName()
                    # taxId=sfObj.getSourceTaxId()
                    # dbName=sfObj.getRefDatabaseName()
                    #refFeatureD[(gId,partId)]=(idRefSeq,altId,seqRefFD,partId, dbName, orgName,strainName,taxId)

                    refFeatureD[(gId, partId)] = (idRefSeq, altId, seqRefFD, partId)
                    # refSeqIdxD[idRefSeq]=seqRefIdx
                    allRefSeqIdxD[idRefSeq] = seqRefIdx
                    #
                    if (self.__verbose):
                        self.__lfh.write("+SequenceDataExport.__exportSeqMapping() from entity group %s part %s exporting reference sequence %s\n" %
                                         (gId, partId, idRefSeq))
                #
            else:
                # report missing reference sequence
                if (self.__verbose):
                    self.__lfh.write("+SequenceDataExport.__exportSeqMapping() missing reference sequence from group %s seqId %s\n" %
                                     (gId, seqId0))
                flagRefOk = False

            #
            # Coordinate sequences -
            #
            idListXyz = []
            # seqXyzIdxD={}
            for seqId in seqIdList:
                if seqId in self.__I['xyz']:
                    for partId, (altId, ver) in self.__I['xyz'][seqId].items():
                        seqLabel.set(seqType='xyz', seqInstId=seqId, seqPartId=partId, seqAltId=altId, seqVersion=ver)
                        idXyzSeq = seqLabel.pack()
                        idListXyz.append(idXyzSeq)
                        #
                        #seqXyzIdxD[idXyzSeq]=self.__sds.getSequence(seqId=seqId,seqType='xyz',partId=partId, altId=altId,version=ver)
                        # seqLenD[idXyzSeq]=len(seqXyzIdxD[idXyzSeq])
                        # if (self.__verbose):
                        #    self.__lfh.write("+SequenceDataExport.__exportSeqMapping() from group %s exporting reference sequence %s length %d\n" %
                        #                     (gId,idXyzSeq,seqLenD[idXyzSeq]))
                else:
                    # report missing coordinate sequence
                    pass

            #
            if flagRefOk:
                #
                alignIdList = []
                alignIdList.append(idAuthSeq)
                # alignIdList.extend(refSeqIdxD.keys())
                alignIdList.extend(idListRef)
                map = MultiAlignPseudo(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
                alignRef = map.run(operation='export', identifier=self.__identifier, alignIdList=alignIdList, alignGroupId=gId, selectedIdList=self.__summarySeqSelectList)

                #alignRef = self.__align(seqAuthIdx,idAuthSeq,refSeqIdxD)
                #
                #rptRefL  = self.__alignRefReport3(gId,alignRef,seqRefFD,rptRefL)
                rptRefL = self.__alignRefReport(gId, alignRef, rptRefL)
                rptCommentL,eelCommentL = self.__authCommentReport(gId, alignRef, rptCommentL)
                if len(eelCommentL) > 0 and sourceInfo == 'NAT':
                    warningMsg += "Entity '" + gId + "' with natural source has '" + "', '".join(eelCommentL) + "' sequence annotation information.<br />\n"
                #
                rptCommentModL = self.__authCommentModReport(gId, alignRef, rptCommentModL)
            #
            if len(idListXyz) > 0:
                alignIdList = []
                alignIdList.append(idAuthSeq)
                # alignIdList.extend(seqXyzIdxD.keys())
                alignIdList.extend(idListXyz)
                map = MultiAlignPseudo(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
                alignXyz = map.run(operation='export', identifier=self.__identifier, alignIdList=alignIdList, alignGroupId=gId, selectedIdList=self.__summarySeqSelectList)
                #alignXyz = self.__align(seqAuthIdx,idAuthSeq,seqXyzIdxD)

                rptXyzL, numC = self.__alignReport(gId, alignXyz, rptXyzL)
                if numC > 0:
                    numConflicts += numC
                    conflictList.append((gId, numC))
                # bail out
            #
            #  END of iteration over  --- gId  ----
        #
        if len(natureSourceTaxIds) > 1:
            warningMsg += "Entry contains multiple natural sources:<br />\n"
            for k,v in natureSourceTaxIds.items():
                if len(v) > 1:
                    warningMsg += "Entities '" + "', '".join(v) + "' have "
                else:
                    warningMsg += "Entity '" + "', '".join(v) + "' has "
                #
                warningMsg += "source taxonomy Id '" + k + "'.<br />\n"
            #
        #
        ##
        ##
        # rptSource=self.__getSourceDetails(refFeatureD)
        #
        # ref FeatureD ---
        #
        rptRefDbL = self.__refDdReport(refFeatureD)
        rptDeleteL = self.__deletionReport(refFeatureD, allRefSeqIdxD)

        try:
            # Make a local copy of the mapping file in the session directory and then copy the file as needed.
            #
            blockList = []
            ofh = open(self.__exportFilePathSession, "w")
            dataBlock = DataContainer(depDataSetId)
            if (len(self.__trueSelfReferenceEntityList) > 0):
                dataBlock.append(self.__cifCatSelfReference(self.__trueSelfReferenceEntityList))
            #
            dataBlock.append(self.__cifCatDbRef(rptRefDbL))
            if (len(self.__annoRefDbL) > 0):
                dataBlock.append(self.__cifCatAnnoDbRef(self.__annoRefDbL))
                for k,seqAuthRefMap in self.__annoSeqAuthRefMap.items():
                    for seqAuth in seqAuthRefMap:
                        seqAuth[0] = str(len(rptRefL) + 1)
                        rptRefL.append(tuple(seqAuth))
                    #
                #
            #
            dataBlock.append(self.__getEntitySourceCategory())
            #
            dataBlock.append(self.__cifCatAuthComment(rptCommentL))
            #
            dataBlock.append(self.__cifCatAuthCommentMod(rptCommentModL))
            if (len(rptDeleteL) > 0):
                dataBlock.append(self.__cifCatDbRefDeletions(rptDeleteL))
            #
            # dataBlock.append(self.__cifCatSourceDetails(rptSource))

            dataBlock.append(self.__cifCatRefMapping(rptRefL))
            dataBlock.append(self.__cifCatXyzMapping(rptXyzL))

            blockList.append(dataBlock)
            pdbxW = PdbxWriter(ofh)
            pdbxW.write(blockList)
            ofh.close()
        except:
            if (self.__verbose):
                self.__lfh.write("+SequenceDataExport.__exportSeqMapping() saving export mapping file failed %s\n" % self.__exportFilePathSession)
                traceback.print_exc(file=self.__lfh)
            return False, numConflicts, conflictList, warningMsg

        return True, numConflicts, conflictList, warningMsg

    def __cifCatDbRef(self, rptRefDbL):
        aCat = DataCategory("pdbx_seqtool_db_ref")
        aCat.appendAttribute("ordinal")
        aCat.appendAttribute("entity_id")
        aCat.appendAttribute("db_name")
        aCat.appendAttribute("db_code")
        aCat.appendAttribute("db_accession")
        aCat.appendAttribute("db_isoform")
        aCat.appendAttribute("match_begin")
        aCat.appendAttribute("match_end")
        aCat.appendAttribute("entity_part_id")
        for row in rptRefDbL:
            aCat.append(row)
        return aCat

    def __cifCatAnnoDbRef(self, annoRefDbL):
        aCat = DataCategory("pdbx_seqtool_anno_db_ref")
        aCat.appendAttribute("ordinal")
        aCat.appendAttribute("entity_id")
        aCat.appendAttribute("db_name")
        aCat.appendAttribute("db_code")
        aCat.appendAttribute("db_accession")
        aCat.appendAttribute("match_begin")
        aCat.appendAttribute("match_end")
        aCat.appendAttribute("entity_part_id")
        for row in annoRefDbL:
            aCat.append(row)
        return aCat

    def __cifCatDbRefDeletions(self, rptDeleteL):
        aCat = DataCategory("pdbx_seqtool_ref_deletions")
        aCat.appendAttribute("ordinal")
        aCat.appendAttribute("entity_id")
        aCat.appendAttribute("part_id")
        aCat.appendAttribute("db_name")
        aCat.appendAttribute("db_accession")
        aCat.appendAttribute("db_isoform")
        aCat.appendAttribute("ref_mon_id")
        aCat.appendAttribute("ref_mon_num")
        for row in rptDeleteL:
            aCat.append(row)
        return aCat

    def __cifCatSelfReference(self, entityList):
        """ Update entities for self reference  ---
        """
        aCat = DataCategory("pdbx_seqtool_self_ref")
        aCat.appendAttribute("entity_id")
        aCat.appendAttribute("entity_part_id")
        for id, partId in entityList:
            aCat.append([id, partId])
        return aCat

    def __cifCatRefMapping(self, rptRefL):
        aCat = DataCategory("pdbx_seqtool_mapping_ref")
        aCat.appendAttribute("ordinal")
        aCat.appendAttribute("entity_id")
        aCat.appendAttribute("entity_mon_id")
        aCat.appendAttribute("entity_seq_num")
        aCat.appendAttribute("auth_asym_id")
        aCat.appendAttribute("ref_mon_id")
        aCat.appendAttribute("ref_mon_num")
        aCat.appendAttribute("entity_part_id")
        for row in rptRefL:
            aCat.append(row)
        return aCat

    def __cifCatAuthComment(self, rptRefL):
        aCat = DataCategory("pdbx_seqtool_mapping_comment")
        aCat.appendAttribute("ordinal")
        aCat.appendAttribute("entity_id")
        aCat.appendAttribute("entity_mon_id")
        aCat.appendAttribute("entity_seq_num")
        # aCat.appendAttribute("auth_asym_id")
        aCat.appendAttribute("comment")
        for row in rptRefL:
            aCat.append(row)
        return aCat

    def __cifCatAuthCommentMod(self, rptRefL):
        aCat = DataCategory("pdbx_seqtool_mapping_modification_comment")
        aCat.appendAttribute("ordinal")
        aCat.appendAttribute("entity_id")
        aCat.appendAttribute("entity_mon_id")
        aCat.appendAttribute("entity_seq_num")
        # aCat.appendAttribute("auth_asym_id")
        aCat.appendAttribute("comment")
        for row in rptRefL:
            aCat.append(row)
        return aCat

    def __cifCatXyzMapping(self, rptXyzL):
        aCat = DataCategory("pdbx_seqtool_mapping_xyz")
        aCat.appendAttribute("ordinal")
        aCat.appendAttribute("entity_id")
        aCat.appendAttribute("entity_mon_id")
        aCat.appendAttribute("entity_seq_num")
        aCat.appendAttribute("auth_asym_id")
        aCat.appendAttribute("pdb_mon_id")
        aCat.appendAttribute("pdb_mon_num")
        aCat.appendAttribute("pdb_ins_code")
        aCat.appendAttribute("hetero_flag")
        for row in rptXyzL:
            aCat.append(row)
        return aCat

    def __refDdReport(self, refFeatureD):
        """  Return a list of content rows for reference sequence database accession report --
        """
        rptRefDbL = []
        irow = 1
        for (gId, partId), tup in refFeatureD.items():
            if ((str(gId), str(partId)) in self.__selfReferenceEntityList):
                continue
            sId = tup[0]
            fD = tup[2]
            rptRefDbL.append(
                (irow,
                 gId,
                 self.__srd.convertDbNameToResource(
                     fD['DB_NAME']),
                    fD['DB_CODE'],
                    fD['DB_ACCESSION'],
                    fD['DB_ISOFORM'],
                    fD['REF_MATCH_BEGIN'],
                    fD['REF_MATCH_END'],
                    partId))
            irow += 1

        return rptRefDbL

    def __deletionReport(self, refFeatureD, allRefSeqIdxD):
        """  Return a list of content rows for the sequence deletion report.
        """
        rptDeleteL = []
        idel = 1
        for (gId, partId), tup in refFeatureD.items():
            if ((str(gId), str(partId)) in self.__selfReferenceEntityList):
                continue
            sId = tup[0]
            fD = tup[2]
            if sId in allRefSeqIdxD:
                if (self.__verbose):
                    self.__lfh.write("+SequenceDataExport.__exportSeqMapping() searching deletions in seqId %s\n" % sId)
                refSeqIdx = allRefSeqIdxD[sId]
                for sTup in refSeqIdx:
                    if sTup[2] in ['Deletion', 'deletion']:
                        if (self.__verbose):
                            self.__lfh.write("+SequenceDataExport.__exportSeqMapping() position %s comment %s\n" % (sTup[1], sTup[2]))
                        rptDeleteL.append((idel, gId, partId, self.__srd.convertDbNameToResource(fD['DB_NAME']), fD['DB_ACCESSION'], fD['DB_ISOFORM'], sTup[0], sTup[1]))
                        idel += 1
            else:
                self.__lfh.write("+SequenceDataExport.__exportSeqMapping() delete scan missing reference sequence for seqId %s keys() %r\n" % (sId, allRefSeqIdxD.keys()))

        return rptDeleteL

    #
    def __alignReport(self, gId, align, oList):
        """  Using a standard alignment ---  generate mapping for Sample and coordinate alignment
        """
        gapSymbol = self.__srd.getGapSymbol()
        numConflicts = 0
        # refLab=align[0][0]
        # sLabel=SequenceLabel()
        # sLabel.unpack(refLab)
        #(seqType,seqInstId,seqPartId,seqAltId,seqVersion)=sLabel.get()
        #
        (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = align[0][1].get()
        refSeq = align[0][2]
        if (self.__debug):
            self.__lfh.write("len refSeq %r\n" % len(refSeq))
        #
        for alTest in align[1:]:
            # sLabel.unpack(alTest[0])
            tstSeq = alTest[2]
            if (self.__debug):
                self.__lfh.write("len tstSeq %r\n" % len(tstSeq))
            (seqTypeT, seqInstIdT, seqPartIdT, seqAltIdT, seqVersionT) = alTest[1].get()
            if (self.__debug):
                self.__lfh.write("gId %r seqTypeT %r seqInstIdT %r seqPartIdT %r seqAltIdT %r seqVersionT %r\n" % (gId, seqTypeT, seqInstIdT, seqPartIdT, seqAltIdT, seqVersionT))
                self.__lfh.write("refSeq %r\n" % len(refSeq))
                #self.__lfh.write("tstSeq %r\n" % tstSeq)
            for idx in range(0, len(refSeq)):
                if str(refSeq[idx][1]) in ['.', '?']:
                    continue
                if str(refSeq[idx][2]) in ['.', '?']:
                    continue

                if (self.__debug):
                    self.__lfh.write("idx %r tstSeq %r\n" % (idx, tstSeq[idx]))

                    # if ( (type =='xyz') and (len(rT[5]) > 0) and (rT[5].find("hetero") != -1 )):
                    #    cssPosClassBg += " bgxyzhetero  "
                    #    ii=rT[5].find("hetero-")
                    #    idS+='_'+rT[5][ii+7:]

                if str(tstSeq[idx][1]) in ['.', '?', '']:
                    # missing test residue
                    tstResId = '.'
                    resIdx = '.'
                    insCode = '.'
                    #
                    # continue
                else:
                    tstResId = str(tstSeq[idx][1])
                    tstResIdx = tstSeq[idx][2]
                    #self.__lfh.write(" testResId %r tstResIdx %r number test %r\n" % (tstResId,tstResIdx,str(tstResIdx).isdigit()))
                    if str(tstSeq[idx][2])[-1].isdigit():
                        insCode = '.'
                        resIdx = str(tstSeq[idx][2])
                    else:
                        insCode = str(tstSeq[idx][2])[-1]
                        resIdx = str(tstSeq[idx][2])[:-1]

                if ((len(tstSeq[idx][5]) > 0) and (tstSeq[idx][5].find("hetero") != -1)):
                    ii = tstSeq[idx][5].find("hetero-")
                    tS = tstSeq[idx][5][ii + 7:]
                    hL = tS.split(':')
                    for h in hL:
                        irow = len(oList) + 1
                        oList.append((str(irow), str(gId), str(refSeq[idx][1]), str(refSeq[idx][2]),
                                      str(seqInstIdT), h, resIdx, insCode, 'Y'))
                else:
                    if ((tstResId is not gapSymbol) and (str(refSeq[idx][1]) != tstResId)):
                        numConflicts += 1
                    irow = len(oList) + 1
                    oList.append((str(irow), str(gId), str(refSeq[idx][1]), str(refSeq[idx][2]),
                                  str(seqInstIdT), tstResId, resIdx, insCode, 'N'))
        #
        return oList, numConflicts

    def __alignRefReport(self, gId, align, oList):
        """  Using a standard alignment list --
        """
        gapSymbol = self.__srd.getGapSymbol()
        # refLab=align[0][0]
        # sLabel=SequenceLabel()
        # sLabel.unpack(refLab)

        (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = align[0][1].get()
        refSeq = align[0][2]
        #
        for idx in range(0, len(refSeq)):
            entity_mon_id = str(refSeq[idx][1])
            entity_seq_num = str(refSeq[idx][2])
            comment = ''
            if (len(refSeq[idx][5]) > 0):
                comment = str(refSeq[idx][5])
            #
            if (comment == 'chromophore') and ((not entity_mon_id) or entity_mon_id == '.') and ((not entity_seq_num) or entity_seq_num == '.'):
                start = idx - 2
                if start < 0:
                    start = 0
                #
                end = idx + 3
                if end > len(refSeq):
                    end = len(refSeq)
                #
                for idx1 in range(start, end):
                    if (len(refSeq[idx1][5]) > 0) and (str(refSeq[idx1][5]) == 'chromophore') and str(refSeq[idx1][1]) and \
                       (str(refSeq[idx1][1]) != '.') and str(refSeq[idx1][2]) and (str(refSeq[idx1][2]) != '.'):
                        entity_mon_id = str(refSeq[idx1][1])
                        entity_seq_num = str(refSeq[idx1][2])
                        break
                    #
                #
            #
            irow = len(oList) + 1
            tstCompId = gapSymbol
            tstSeqNum = '.'
            tstPartId = '.'
            for alTest in align[1:]:
                # sLabel.unpack(alTest[0])
                #(seqTypeT,seqInstIdT,seqPartIdT,seqAltIdT,seqVersionT)=sLabel.get()
                (seqTypeT, seqInstIdT, seqPartIdT, seqAltIdT, seqVersionT) = alTest[1].get()
                tstSeq = alTest[2]
                if (tstSeq[idx][1] != gapSymbol):
                    tstCompId = tstSeq[idx][1]
                    tstSeqNum = tstSeq[idx][2]
                    tstPartId = seqPartIdT
            oList.append((str(irow), str(gId), entity_mon_id, entity_seq_num, str(seqInstIdT), str(tstCompId), str(tstSeqNum), str(tstPartId)))
        #
        return oList

    def __authCommentModReport(self, gId, align, oList):
        """
        """
        # refLab=align[0][0]
        # sLabel=SequenceLabel()
        # sLabel.unpack(refLab)
        (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = align[0][1].get()
        refSeq = align[0][2]
        #
        for idx in range(0, len(refSeq)):
            if (len(refSeq[idx][5]) > 0):
                comment = str(refSeq[idx][5])
                if comment in ["modified residue", "microheterogeneity/modified residue"]:
                    irow = len(oList) + 1
                    oList.append((str(irow), str(gId), str(refSeq[idx][1]), str(refSeq[idx][2]), 'modified residue'))
            #
        return oList

    def __authCommentReport(self, gId, align, oList):
        """
        """
        eelCommentL = []
        # refLab=align[0][0]
        # sLabel=SequenceLabel()
        # sLabel.unpack(refLab)
        (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = align[0][1].get()
        refSeq = align[0][2]
        #
        for idx in range(0, len(refSeq)):
            if (len(refSeq[idx][5]) > 0):
                comment = str(refSeq[idx][5])
                for cType in ( 'engineered mutation', 'expression tag', 'linker' ):
                    if (comment.find(cType) != -1) and (cType not in eelCommentL):
                        eelCommentL.append(cType)
                    #
                #
                if comment in ["microheterogeneity/modified residue"]:
                    comment = "microheterogeneity"
                if comment in ["modified residue", "deletion"]:
                    continue
                if (comment == 'chromophore') and ((not str(refSeq[idx][1])) or str(refSeq[idx][1]) == '.') and \
                   ((not str(refSeq[idx][2])) or str(refSeq[idx][2]) == '.'):
                    continue
                irow = len(oList) + 1
                oList.append((str(irow), str(gId), str(refSeq[idx][1]), str(refSeq[idx][2]), comment))

            #
        #
        return oList,eelCommentL

    def __getEntitySourceCategory(self):
        """ Return category object of source details.
        """
        itemTupList = [('ORDINAL', 'ordinal'),
                       ('ENTITY_ID', 'entity_id'),
                       ('AUTH_SEQ_PART_ID', 'seq_part_id'),
                       ('AUTH_SEQ_NUM_BEGIN', 'seq_part_beg'),
                       ('AUTH_SEQ_NUM_END', 'seq_part_end'),
                       ('AUTH_SEQ_PART_TYPE', 'seq_part_type'),
                       ('ENTITY_DESCRIPTION', 'entity_description'),
                       ('ENTITY_SYNONYMS', 'entity_synonyms'),
                       ('SOURCE_GENE_NAME', 'gene_name'),
                       ('SOURCE_TAXID', 'taxonomy_id'),
                       ('SOURCE_ORGANISM', 'source_scientific_name'),
                       ('SOURCE_STRAIN', 'source_strain'),
                       ('SOURCE_COMMON_NAME', 'source_common_name'),
                       ('SOURCE_VARIANT', 'variant'),
                       ('ENTITY_ENZYME_CLASS', 'entity_enzyme_class'),
                       ('ENTITY_FRAGMENT_DETAILS', 'entity_fragment_details'),
                       ('ENTITY_MUTATION_DETAILS', 'entity_mutation_details'),
                       ('ENTITY_DETAILS', 'entity_details'),
                       ('SOURCE_METHOD', 'source_method'),
                       ('HOST_ORG_SOURCE', 'host_org_scientific_name'),
                       ('HOST_ORG_STRAIN', 'host_org_strain'),
                       ('HOST_ORG_TAXID', 'host_org_taxonomy_id'),
                       ('HOST_ORG_VECTOR', 'host_org_vector'),
                       ('HOST_ORG_VECTOR_TYPE', 'host_org_vector_type'),
                       ('HOST_ORG_PLASMID', 'host_org_plasmid'),
                       ('HOST_ORG_COMMON_NAME', 'host_org_common_name'),
                       ('HOST_ORG_CELL_LINE', 'host_org_cell_line'),
                       ('HOST_ORG_VARIANT', 'host_org_variant'),
                       ]

        if (self.__verbose):
            self.__lfh.write("+SequenceDataExport.__getSourceDetails() starting\n")
            # self.printIt()
        #
        rowList = []
        entityIdList = self.__sds.getGroupIds()
        iRow = 1
        for entityId in entityIdList:
            seqIds = self.__sds.getGroup(entityId)
            if len(seqIds) < 1:
                if (self.__verbose):
                    self.__lfh.write("+SequenceDataExport.__getSourceDetails() entityId %s is empty\n" % (entityId))
                continue

            # seqId0=seqIds[0]
            seqId0 = entityId
            partIdList = self.__sds.getPartIds(seqId0, dataType="sequence", seqType="auth")
            #
            # JDW JDW
            for partId in partIdList:
                if (self.__verbose):
                    self.__lfh.write("+SequenceDataExport.__getSourceDetails() entityId %r partId %r  __I  %r\n" % (entityId, partId, self.__I['auth'][seqId0][partId]))
                # JDW
                # altId,versionId=self.__I['auth'][seqId0][partId]
                altId, versionId = self.__I['auth'][seqId0][1]
                fD = self.__sds.getFeature(seqId0, seqType="auth", partId=partId, altId=altId, version=versionId)
                d = {}
                for itemTup in itemTupList:
                    ky = itemTup[0]
                    if ky in fD:
                        d[ky] = fD[ky]
                    else:
                        d[ky] = ''
                d['ORDINAL'] = iRow
                d['ENTITY_ID'] = entityId
                d['AUTH_SEQ_PART_ID'] = partId
                rowList.append(d)
                iRow += 1
                #
        #
        aCat = DataCategory("pdbx_seqtool_entity_details")
        for itemTup in itemTupList:
            attribName = itemTup[1]
            aCat.appendAttribute(attribName)
        for rD in rowList:
            vL = []
            for itemTup in itemTupList:
                ky = itemTup[0]
                vL.append(rD[ky])
            aCat.append(vL)
        return aCat

    ##
    # OBSOLETE CODE
    ##
    def dump(self, io):
        """ Output the original sequence lists and the resulting alignment lists -
        """
        io.write("------------------------------\n")
        io.write("Sequence %s\n" % self.__refSeqId[0])
        for tup in self.__refSeqIdx:
            io.write("   + %8s %8s %8s \n" % (tup[0], tup[1], tup[2]))

        for aSeq in self.__seqList:
            io.write("------------------------------\n")
            io.write("Sequence %s\n" % aSeq[0])
            for tup in aSeq[3]:

                io.write("   + %8s %8s %8s \n" % (tup[0], tup[1], tup[2]))

        io.write("------------------------------\n")
        io.write("ALIGNED SEQUENCES\n")
        self.__formatAlignment(io)

    def __XformatAlignment(self, io):
        alignLength = len(self.__alignSeqList[0][2])
        io.write("\n+SequenceDataExport.__formatAlignment() - Dumping alignment list data - %d\n" % alignLength)
        for aPos in range(0, alignLength):
            io.write("%5d:  " % aPos)
            for aTup in self.__alignSeqList:
                id = aTup[0]
                sLabel = aTup[1]
                aL = aTup[2]
                conflictL = aTup[3]
                # (one-letter-code, 3-letter-code, original label residue index, position in sequence )
                rT = aL[aPos]
                io.write("%s(%3s) %5s(%5s)(%r) - " % (rT[0], rT[1], str(rT[2]), str(rT[3]), conflictL[aPos][0]))
            io.write("\n")
        io.flush()

    def __XannotateAlignmentWithIndex(self, s3L, s3WithIndexL):
        """
        Input is the aligned sequence as a list of 3-letter-codes or gap-symbols.

        Add one-letter code and residue index in the original coordinate sequence and
        returns a list of tuples containing -

        (one-letter-code, 3-letter-code, original label residue index, position in sequence, position in alignment, comment)
        """

        gapSymbol = self.__srd.getGapSymbol()
        #
        ir = 0
        annotL = []
        for idx in range(0, len(s3L)):

            if ((ir < len(s3WithIndexL)) and (s3WithIndexL[ir][0] == s3L[idx])):
                # self.__lfh.write("orgidx-> %s ir %d s3l-> %s  idx %d len s3L %d len s3WithIndexL %d\n" %
                #                 (s3WithIndexL[ir][0],ir,s3L[idx],idx,len(s3L),len(s3WithIndexL)))
                annotL.append((self.__srd.cnv3To1(s3L[idx]), s3L[idx], s3WithIndexL[ir][1], ir, idx, s3WithIndexL[ir][2]))
                ir += 1
            elif s3L[idx] == gapSymbol:
                annotL.append((gapSymbol, gapSymbol, '', '', idx, ''))
            else:
                self.__lfh.write("+SequenceDataExport.__annotateAlignmentWithIndex() error at index %d %s %d %d\n "
                                 % (idx, s3L[idx], ir, len(s3WithIndexL)))
        return annotL

    def __Xalign(self, seqRefIdx, idRef, testSeqD):
        """ Align the test sequences against the reference sequence.  Produce a list
            of aligned sequences annotated with index information from the input
            sequences.

            Return data storage model for aligned sequence list -

            alignSeqList [[Index from sequence data store,  aligned sequence with index details, conflict flag list (bool)],,,,]
        """
        #
        if (self.__verbose):
            self.__lfh.write("+SequenceDataExport.__align() Beginning  alignment with reference %s\n" % idRef)

        pA = PairwiseAlign()
        pA.setVerbose(self.__verbose)
        r3L = []
        for tup in seqRefIdx:
            r3L.append(tup[0])
        pA.setReferenceSequence(r3L, idRef)

        idList = []
        for id, seqIdx in testSeqD.items():
            r3L = []
            typicalLink = []
            for tup in seqIdx:
                r3L.append(tup[0])
                typicalLink.append(0 if ('long_begin' in tup[2]) else 1)
            if (self.__debug):
                #self.__lfh.write("%s (%d) -> %r\n" % (id,len(r3L),r3L))
                pass
            pA.addTestSequenceWithLink(r3L, id, typicalLink)
            idList.append(id)
        #
        pA.doAlign()
        #        if (self.__verbose):
        #            pA.wrAlignmentFull(ostream(self.__lfh))
        #            self.__lfh.flush()
        #

        # Post process the aligned sequences -
        #
        # Get the reference sequence from the first alignment -
        #
        aL0 = pA.getAlignment(idList[0])
        alignRefSeq = []
        for tup in aL0:
            alignRefSeq.append(tup[0])

        alignRefSeqIdx = self.__annotateAlignmentWithIndex(alignRefSeq, seqRefIdx)
        lenRefAlignment = len(alignRefSeqIdx)
        if (self.__verbose):
            self.__lfh.write("+SequenceDataExport.__align() First     alignment length for reference %s is %d\n"
                             % (idRef, len(alignRefSeq)))
            self.__lfh.write("+SequenceDataExport.__align() Annotated alignment length for reference %s is %d\n"
                             % (idRef, len(alignRefSeqIdx)))
        llen = len(aL0)
        for id, seqIdx in testSeqD.items():
            aL = pA.getAlignment(id)
            llen = max(llen, len(aL))
        #
        alignSeqList = []
        # conflictRefL=[(0,'')]*len(aL0)
        conflictRefL = [(0, '')] * llen
        alignSeqList.append([idRef, alignRefSeqIdx, conflictRefL])

        for id, seqIdx in testSeqD.items():
            aL = pA.getAlignment(id)
            alignSeq = []
            conflictL = [(0, '')] * len(aL)
            if (self.__verbose):
                self.__lfh.write("+SequenceDataExport.__align() reference %s and test %s length %d\n"
                                 % (idRef, id, len(aL)))
            ii = 0
            numConflicts = 0
            for aTup in aL:
                alignSeq.append(aTup[1])
                if aTup[0] != aTup[1]:
                    conflictL[ii] = 1
                    conflictRefL[ii] = (1, aTup[0])
                    numConflicts += 1
                ii += 1

            if (self.__verbose):
                self.__lfh.write("+SequenceDataExport.__align() %s conflicts %d length alignment %d length indexed sequence %d\n"
                                 % (id, numConflicts, len(alignSeq), len(seqIdx)))

            alignSeqIdx = self.__annotateAlignmentWithIndex(alignSeq, seqIdx)
            alignSeqList.append([id, alignSeqIdx, conflictL])

        if (self.__verbose):
            numConflicts = 0
            for conflictTup in alignSeqList[0][2]:
                if conflictTup[0] != 0:
                    numConflicts += 1
            self.__lfh.write("+SequenceDataExport.__align() Leaving with total conflicts between test and reference sequence = %d\n" % numConflicts)

        return alignSeqList

    def __XalignRefReportOld(self, gId, align, seqRefFD, oList):
        """
        """
        refLab = align[0][0]
        sLabel = SequenceLabel()
        sLabel.unpack(refLab)
        (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = sLabel.get()
        refSeq = align[0][1]
        #
        for alTest in align[1:]:
            sLabel.unpack(alTest[0])
            tstSeq = alTest[1]
            (seqTypeT, seqInstIdT, seqPartIdT, seqAltIdT, seqVersionT) = sLabel.get()
            for idx in range(0, len(refSeq)):
                irow = len(oList) + 1
                oList.append((str(irow), str(gId), str(refSeq[idx][1]), str(refSeq[idx][2]), str(seqInstIdT),
                              str(tstSeq[idx][1]), str(tstSeq[idx][2]), str(seqPartIdT)))
        #
        return oList

    def __XalignReport3(self, gId, align, oList):
        """  Sample and coordinate alignment
        """
        gapSymbol = self.__srd.getGapSymbol()
        numConflicts = 0
        refLab = align[0][0]
        sLabel = SequenceLabel()
        sLabel.unpack(refLab)
        (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = sLabel.get()
        refSeq = align[0][1]
        if (self.__debug):
            self.__lfh.write("len refSeq %r\n" % len(refSeq))
        #
        for alTest in align[1:]:
            sLabel.unpack(alTest[0])
            tstSeq = alTest[1]
            if (self.__debug):
                self.__lfh.write("len tstSeq %r\n" % len(tstSeq))
            (seqTypeT, seqInstIdT, seqPartIdT, seqAltIdT, seqVersionT) = sLabel.get()
            if (self.__debug):
                self.__lfh.write("gId %r seqTypeT %r seqInstIdT %r seqPartIdT %r seqAltIdT %r seqVersionT %r\n" % (gId, seqTypeT, seqInstIdT, seqPartIdT, seqAltIdT, seqVersionT))
                self.__lfh.write("refSeq %r\n" % len(refSeq))
                #self.__lfh.write("tstSeq %r\n" % tstSeq)
            for idx in range(0, len(refSeq)):
                if str(refSeq[idx][1]) in ['.', '?']:
                    continue
                if str(refSeq[idx][2]) in ['.', '?']:
                    continue

                if (self.__debug):
                    self.__lfh.write("idx %r tstSeq %r\n" % (idx, tstSeq[idx]))

                if str(tstSeq[idx][1]) in ['.', '?', '']:
                    # missing test residue
                    tstResId = '.'
                    resIdx = '.'
                    insCode = '.'
                    #
                    # continue
                else:
                    tstResId = str(tstSeq[idx][1])
                    tstResIdx = tstSeq[idx][2]
                    #self.__lfh.write(" testResId %r tstResIdx %r number test %r\n" % (tstResId,tstResIdx,str(tstResIdx).isdigit()))
                    if str(tstSeq[idx][2])[-1].isdigit():
                        insCode = '.'
                        resIdx = str(tstSeq[idx][2])
                    else:
                        insCode = str(tstSeq[idx][2])[-1]
                        resIdx = str(tstSeq[idx][2])[:-1]

                irow = len(oList) + 1
                # oList.append(( str(irow), str(gId), str(refSeq[idx][1]), str(refSeq[idx][2]),
                #               str(seqInstIdT), str(tstSeq[idx][1]), str(tstSeq[idx][2])))
                if ((tstResId is not gapSymbol) and (str(refSeq[idx][1]) != tstResId)):
                    numConflicts += 1
                oList.append((str(irow), str(gId), str(refSeq[idx][1]), str(refSeq[idx][2]),
                              str(seqInstIdT), tstResId, resIdx, insCode))
        #
        return oList, numConflicts

    def __XalignRefReport3(self, gId, align, seqRefFD, oList):
        """
        """
        gapSymbol = self.__srd.getGapSymbol()
        refLab = align[0][0]
        sLabel = SequenceLabel()
        sLabel.unpack(refLab)
        (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = sLabel.get()
        refSeq = align[0][1]
        #
        for idx in range(0, len(refSeq)):
            irow = len(oList) + 1
            tstCompId = gapSymbol
            tstSeqNum = '.'
            tstPartId = '.'
            for alTest in align[1:]:
                sLabel.unpack(alTest[0])
                (seqTypeT, seqInstIdT, seqPartIdT, seqAltIdT, seqVersionT) = sLabel.get()
                tstSeq = alTest[1]
                if (tstSeq[idx][1] != gapSymbol):
                    tstCompId = tstSeq[idx][1]
                    tstSeqNum = tstSeq[idx][2]
                    tstPartId = seqPartIdT
            oList.append((str(irow), str(gId), str(refSeq[idx][1]), str(refSeq[idx][2]), str(seqInstIdT),
                          str(tstCompId), str(tstSeqNum), str(tstPartId)))
        #
        return oList

    def __XauthCommentModRepor3(self, gId, align, oList):
        """
        """
        refLab = align[0][0]
        sLabel = SequenceLabel()
        sLabel.unpack(refLab)
        (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = sLabel.get()
        refSeq = align[0][1]
        #
        for idx in range(0, len(refSeq)):
            if (len(refSeq[idx][5]) > 0):
                comment = str(refSeq[idx][5])
                if comment != "modified residue":
                    continue
                irow = len(oList) + 1
                oList.append((str(irow), str(gId), str(refSeq[idx][1]), str(refSeq[idx][2]), str(seqInstId), comment))
            #
        return oList

    def __XauthCommentReport3(self, gId, align, oList):
        """
        """
        refLab = align[0][0]
        sLabel = SequenceLabel()
        sLabel.unpack(refLab)
        (seqType, seqInstId, seqPartId, seqAltId, seqVersion) = sLabel.get()
        refSeq = align[0][1]
        #
        for idx in range(0, len(refSeq)):
            if (len(refSeq[idx][5]) > 0):
                comment = str(refSeq[idx][5])
                if comment == "modified residue":
                    continue
                irow = len(oList) + 1
                oList.append((str(irow), str(gId), str(refSeq[idx][1]), str(refSeq[idx][2]),
                              str(seqInstId), comment))
            #
        return oList
