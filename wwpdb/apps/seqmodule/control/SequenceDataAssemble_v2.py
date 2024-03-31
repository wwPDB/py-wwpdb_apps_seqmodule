##
# File:  SequenceDataAssemble.py
# Date:  25-Feb-2013
#
# Updates:
#  26-Feb-2013 jdw  assign polymer parts to sequences and features on input
#  03-Mar-2013 jdw  create rendered table of long intra-polymer linkages
#  03-Mar-2013 jdw  upldate new part data items in sequence data store.
#  05-Mar-2013 jdw  nuke and pave
#  29-Nov-2013 jdw  save model file name in UtilDataStore() for later use by visualization app.
#   3-Dec-2013 jdw  skip PDB format data file production
#  19-Jan-2014 jdw  add 'HOST_ORG_CELL_LINE'
#  19-Jan-2014 jdw  standardize path polymer link report using getPolyLinkReportFilePath()
#  01-Feb-2015 jdw  add host org variant
#  24-Aug-2017  zf  remove __doAssembleFromArchive() & __doAssembleFromModel(), move all file copy operations into DataImporter.py
#                   modify __calcPolymerLinkages(), also modify backend C++ code to speed up file process
#                   modify __doReferenceSearch() to use previous blast results if existing & applicable
#  01-Oct-2022  zf  modified __getAuthDefinedRefD() method to include author provided "seq_align_begin", "seq_align_end", "db_align_end",
#                   "db_align_beg", & "db_seq_one_letter_code" information
#                   improved process to handle the author-provided reference IDs and new entity summary page
#  06-Jan-2024  zf  using relevant UniProt sequence from Sequence Builder to build the reference sequence information instead of BLAST searching
##
"""
Assemble sequence and other required data to start/restart the sequence editor tool.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.09"

try:
    import cPickle as pickle
except ImportError:
    import pickle as pickle
#
import json
import multiprocessing
import os
import shutil
import sys
import time
import traceback

from mmcif.api.DataCategory import DataCategory
from mmcif.api.PdbxContainers import DataContainer
from mmcif.io.PdbxWriter import PdbxWriter

from wwpdb.apps.seqmodule.align.AlignmentTools import AlignmentTools
from wwpdb.apps.seqmodule.io.AlignmentDataStore import AlignmentDataStore
from wwpdb.apps.seqmodule.io.SequenceDataExport_v2 import SequenceDataExport
from wwpdb.apps.seqmodule.link.PolymerLinkageDepict import PolymerLinkageDepict
from wwpdb.apps.seqmodule.update.UpdatePolymerEntitySourceDetails import UpdatePolymerEntitySourceDetails
from wwpdb.apps.seqmodule.util.LocalBlastSearchUtils import LocalBlastSearchUtils
from wwpdb.apps.seqmodule.util.SeqReferenceSearchUtils import SeqAnnotationSearchUtils
# from wwpdb.apps.seqmodule.util.SeqReferenceSearchUtils import SeqReferenceSearchUtils
from wwpdb.apps.seqmodule.util.SequenceLabel import SequenceLabel
from wwpdb.apps.seqmodule.util.UpdateSequenceDataStoreUtils import UpdateSequenceDataStoreUtils
from wwpdb.apps.seqmodule.util.MultiProcLimit import MultiProcLimit
from wwpdb.utils.config.ConfigInfo import ConfigInfo
from rcsb.utils.multiproc.MultiProcUtil import MultiProcUtil
from wwpdb.io.locator.PathInfo import PathInfo
from wwpdb.utils.dp.RcsbDpUtility import RcsbDpUtility
from wwpdb.utils.session.UtilDataStore import UtilDataStore


class SequenceDataAssemble(UpdateSequenceDataStoreUtils):
    """
    This class encapsulates the data assembly operations
    of model sequence and reference data matching data for the sequence
    editor tool.

    Storage model - imported data is loaded into the sequence data store
                    where it is managed by the SequenceDataStore() class.

    Methods in this class extract the author and coordinate sequence data
    and source information from PDBx CIF.

    Reference sequence data is extracted from the processed BLAST results
    (ie. wwPDB seqdb-match data files) for each polymer entity.
    """

    def __init__(self, reqObj=None, verbose=False, log=sys.stderr):
        super(SequenceDataAssemble, self).__init__(reqObj=reqObj, seqDataStore=None, verbose=verbose, log=log)
        self.__sessionObj = self._reqObj.getSessionObj()
        self.__sessionPath = self.__sessionObj.getPath()
        self.__siteId = self._reqObj.getValue("WWPDB_SITE_ID")
        self.__cI = ConfigInfo(self.__siteId)
        self.__pI = PathInfo(siteId=self.__siteId, sessionPath=self.__sessionPath, verbose=self._verbose, log=self._lfh)
        #
        self.__dataSetId = str(self._reqObj.getValue("identifier")).upper()
        self.__reRunBlastSearchFlag = False
        self.__doRefSearchFlag = False
        self.__forceBlastSearchFlag = False
        self.__resetAlignmentFlag = False
        self.__misMatchPath = self.__pI.getFilePath(self.__dataSetId, contentType="mismatch-warning", formatType="pic", fileSource="session")
        #
        # minimum sequence length to search in reference databases -
        #
        self.__minSearchSequenceLengthAA = 9
        self.__minSearchSequenceLengthNA = 50
        self.__entityIdList = []
        #
        self.__authDefinedRefD = {}
        #
        self.__autoProcessFlag = False

    def doAssemble(self, calPolyLink=True, doRefSearch=False, doAutoProcess=False):
        """Perform all steps required to construct the session sequence data store using data from within
        the input 'fileSource'.
        """
        self.__doRefSearchFlag = doRefSearch
        #
        if self._verbose:
            self._lfh.write("\n+SequenceDataAssemble.doAssemble() STARTING with depId %r calPolyLink %r doRefSearch %r doAutoProcess %r\n" %
                            (self.__dataSetId, calPolyLink, doRefSearch, doAutoProcess))
        #
        pdbxFilePath = self.__pI.getModelPdbxFilePath(self.__dataSetId, fileSource="session", versionId="1")
        polyLinkPath = self.__pI.getFilePath(self.__dataSetId, contentType="polymer-linkage-distances", formatType="json", fileSource="session")
        #
        # Calculate polymer linkage information
        #
        if calPolyLink or (not os.access(polyLinkPath, os.R_OK)):
            self.__calcPolymerLinkages(pdbxFilePath=pdbxFilePath, polyLinkPath=polyLinkPath)
        #
        if not os.access(polyLinkPath, os.R_OK):
            return False
        #
        ok, no_coord_issue, seqSimilarityInfoList, polymerLinkDistList, missingResidueInfoList, mixMseMetResidueList, entryD, entityD, instanceD, statisticsMap, depSeqAssign, \
            seqAssign = self.__readJsonObject(polyLinkPath)
        #
        if not ok:
            return False
        #
        if doAutoProcess:
            self.__autoProcessFlag = no_coord_issue
            if len(seqSimilarityInfoList) > 0:
                self.__autoProcessFlag = False
            #
        #
        self.__getAuthDefinedRefD(depSeqAssign)
        #
        entryD["PDBX_MODEL_FILE_PATH"] = pdbxFilePath
        #
        # Render polymer linkage information
        #
        polyLinkHtmlPath = self.__pI.getPolyLinkReportFilePath(self.__dataSetId, fileSource="session")
        self.__renderPolymerLinkages(polymerLinkDistList=polymerLinkDistList, polyLinkHtmlPath=polyLinkHtmlPath)
        #
        # Update SequenceDataStore object - reset SequenceDataStore object here in order to incorporate the partial annotated result(s)
        self.resetSequenceDataStore()
        #
        prevModelFilePath, prevSavedFilePathMap = self.__checkIfPartialAnnotatedResultExist()
        newOldEntityIdMap, misMatchList, notFoundMatchList = self.__getPartialAnnotatedResult(prevModelFilePath, prevSavedFilePathMap, entityD, instanceD)
        if newOldEntityIdMap:
            self.__doRefSearchFlag = False
        #
        self.__entityIdList = list(entityD.keys())
        searchEntityIdList = []
        for entityId in self.__entityIdList:
            if entityId in newOldEntityIdMap:
                continue
            #
            searchEntityIdList.append(entityId)
        #
        eRefD, ownRefD, eSSRefD = self.__doReferenceSearch(entityD, searchEntityIdList)
        #
        if self._verbose:
            self._lfh.write("+SequenceDataAssemble.doAssemble() - Reference match length    %d\n" % len(eRefD))
        #
        sortedEntityIdList = []
        try:
            intEntityIdList = []
            for entityId in self.__entityIdList:
                intEntityIdList.append(int(entityId))
            #
            intEntityIdList.sort()
            for iEntityId in intEntityIdList:
                sortedEntityIdList.append(str(iEntityId))
            #
        except:  # noqa: E722 pylint: disable=bare-except
            sortedEntityIdList = self.__entityIdList
        #
        self.__createSeqMatchFileAndMisMatchPickleFile(sortedEntityIdList, entityD, ownRefD, eSSRefD, list(newOldEntityIdMap.keys()), misMatchList,
                                                       notFoundMatchList, seqSimilarityInfoList, missingResidueInfoList, mixMseMetResidueList)
        #
        self.__updateUtilDataStore(pdbxFilePath)
        #
        groupDict = {}
        groupPartD = {}
        polymerTypeCode = {}
        skipInstanceIdList = []
        for eId, eD in entityD.items():
            # sId = eD["ENTITY_ID"]
            if eId in newOldEntityIdMap:
                skipInstanceIdList.extend(eD["INSTANCE_LIST"])
            #
            groupDict[eId] = eD["INSTANCE_LIST"]
            groupPartD[eId] = list(range(1, len(eD["PART_LIST"]) + 1))
            if self._verbose:
                self._lfh.write("+SequenceDataAssemble.doAssemble() entity %s PART_LIST %r\n" % (eId, eD["PART_LIST"]))
                self._lfh.write("+SequenceDataAssemble.doAssemble() entity %s numParts %d groupPartD  %r\n" % (eId, len(eD["PART_LIST"]), groupPartD[eId]))
            for inst in eD["INSTANCE_LIST"]:
                polymerTypeCode[inst] = eD["POLYMER_TYPE_CODE"]
                if self._verbose:
                    self._lfh.write("+SequenceDataAssemble.doAssemble()  instance %s  type %s\n" % (inst, polymerTypeCode[inst]))
                #
            #
        #
        startTime = time.time()
        if self._verbose:
            self._lfh.write(
                "+SequenceDataAssemble loadSequenceDataStore for sessionId %s started at %s \n" % (self.__sessionObj.getId(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))
            )
        #
        if ("PDB_ID" in entryD["PDB_ID"]) and (len(entryD["PDB_ID"]) == 4):
            self.setPdbId(entryD["PDB_ID"])
        #
        self.setPolymerTypeCode(polymerTypeCode)
        self.setDepositorSeqAssign(depSeqAssign)
        self.setArchiveSeqAssign(seqAssign)
        self.setCoordinateInstanceInfo(instanceD, statisticsMap, skipList=skipInstanceIdList)
        self.setMultipleEntitiesSequenceInfo(entityD, skipList=list(newOldEntityIdMap.keys()))
        self.setMultipleEntitiesRefSequenceInfo(entityD, eRefD, ownRefD, eSSRefD, skipList=list(newOldEntityIdMap.keys()))
        self.setGroup(groupDict)
        self.setGroupParts(groupPartD)
        self.setEntryDetail(entryD)
        self.setDefaultSelectedIds()
        self.saveSequenceDataStore()
        #
        if self._verbose:
            self.dumpSequenceDataStore()
        #
        self.__createDefaultAlignments()
        #
        # Export the assignment file for auto-complete sequence processing
        if self.__autoProcessFlag:
            defaultSelectedIdList = self.getSelectedIds()
            sdu = UpdatePolymerEntitySourceDetails(reqObj=self._reqObj, verbose=self._verbose, log=self._lfh)
            sdu.updateAuthEntityDetails(selectIdList=defaultSelectedIdList)
            #
            dataExport = SequenceDataExport(reqObj=self._reqObj, exportList=defaultSelectedIdList, verbose=self._verbose, log=self._lfh)
            ok, numConflicts, _conflictList, _warningMsg = dataExport.exportAssignments()
            if (not ok) or (numConflicts > 0):
                assignFilePath = self.__pI.getSequenceAssignmentFilePath(self.__dataSetId, fileSource="session", versionId="latest")
                if assignFilePath and os.access(assignFilePath, os.R_OK):
                    os.remove(assignFilePath)
                #
            #
        #
        if self._verbose:
            endTime = time.time()
            self._lfh.write(
                "+SequenceDataAssemble loadSequenceDataStore completed at %s (%.2f seconds)\n" % (time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
            )
            self._lfh.write("+SequenceDataAssemble.doAssemble() COMPLETED for %r\n\n" % self.__dataSetId)
        #
        return True

    def reRunSearch(self):
        """Re-run reference sequence search for selected entities"""
        selectedEntityIdList = str(self._reqObj.getValue("entityids")).strip().split(",")
        if not selectedEntityIdList:
            return True
        #
        self.__reRunBlastSearchFlag = True
        self._entityAlignInfoMap = {}
        self.__reRunSearch(selectedEntityIdList)
        self.__updateSelections(selectedEntityIdList)
        self.saveSequenceDataStore()
        return True

    def doUpdateSelections(self):
        """ """
        entityId = str(self._reqObj.getValue("activegroupid")).strip()
        if not entityId:
            return True
        #
        self.__doRefSearchFlag = False
        self._entityAlignInfoMap = {}
        #
        seqSearchOp = str(self._reqObj.getValue("seq_search_op"))
        if seqSearchOp == "on":
            flag = False
            withRefInfo = str(self._reqObj.getValue("withref_info"))
            if withRefInfo == "yes":
                flag = True
            #
            self.__minSearchSequenceLengthAA = 5
            self.__minSearchSequenceLengthNA = 15
            self.__reRunSearch([entityId], withRefFlag=flag)
        #
        self.__updateSelections([entityId])
        self.saveSequenceDataStore()
        return True

    def getEntityIdList(self):
        return self.__entityIdList

    def setLogHandle(self, log=sys.stderr):
        """Reset the stream for logging output."""
        try:
            self._lfh = log
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            return False

    def __calcPolymerLinkages(self, pdbxFilePath=None, polyLinkPath=None):
        #
        # Calculate intra-polymer linkage distances -
        #
        try:
            dp = RcsbDpUtility(tmpPath=self.__sessionPath, siteId=self.__siteId, verbose=self._verbose, log=self._lfh)
            dp.imp(pdbxFilePath)
            dp.op("annot-poly-link-dist-json")
            dp.exp(polyLinkPath)
            dp.cleanup()
            if self._verbose:
                self._lfh.write("+SequenceDataAssemble.__calcPolymerLinkages() - Polymer link distance file copied to: %s\n" % polyLinkPath)
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            if self._verbose:
                self._lfh.write("+SequenceDataAssemble.__calcPolymerLinkages() - Polymer link distance calculation failed for file %s\n" % pdbxFilePath)
                traceback.print_exc(file=self._lfh)
            #
        #
        return False

    def __readJsonObject(self, jsonFilePath):
        """ """
        try:
            no_coord_issue = True
            seqSimilarityInfoList = []
            polymerLinkDistList = []
            missingResidueInfoList = []
            mixMseMetResidueList = []
            entryD = {}
            entityD = {}
            instanceD = {}
            statisticsMap = {}
            depSeqAssign = {}
            seqAssign = {}
            with open(jsonFilePath) as DATA:
                jsonObj = json.load(DATA)
                if ("HAS_SEQ_COOR_ISSUE" in jsonObj) and jsonObj["HAS_SEQ_COOR_ISSUE"]:
                    no_coord_issue = False
                #
                if ("SEQ_ITENTITY_INFO" in jsonObj) and jsonObj["SEQ_ITENTITY_INFO"]:
                    seqSimilarityInfoList = jsonObj["SEQ_ITENTITY_INFO"]
                #
                if ("linkageL" in jsonObj) and jsonObj["linkageL"]:
                    polymerLinkDistList = jsonObj["linkageL"]
                #
                if ("entryD" in jsonObj) and jsonObj["entryD"]:
                    entryD = jsonObj["entryD"]
                #
                #
                if ("entityD" in jsonObj) and jsonObj["entityD"]:
                    entityD = jsonObj["entityD"]
                #
                if ("instanceD" in jsonObj) and jsonObj["instanceD"]:
                    instanceD = jsonObj["instanceD"]
                #
                if ("alignStatisticsD" in jsonObj) and jsonObj["alignStatisticsD"]:
                    statisticsMap = jsonObj["alignStatisticsD"]
                #
                if ("depSeqAssign" in jsonObj) and jsonObj["depSeqAssign"]:
                    depSeqAssign = jsonObj["depSeqAssign"]
                    for key, valD in depSeqAssign.items():
                        if ("ref_list" not in valD) or (not valD["ref_list"]):
                            continue
                        #
                        for refD in valD["ref_list"]:
                            for item in ("db_name", "db_accession", "db_code", "seq_align_begin", "seq_align_end", "db_align_end",
                                         "db_align_beg", "db_seq_one_letter_code"):
                                if (item in refD) and refD[item]:
                                    refD[item] = refD[item].strip().upper()
                                #
                            #
                        #
                    #
                #
                if ("annSeqAssign" in jsonObj) and jsonObj["annSeqAssign"]:
                    seqAssign = jsonObj["annSeqAssign"]
                #
            #
            for key, valD in entityD.items():
                if ("SEQ_TUP_LIST" not in valD) or (not valD["SEQ_TUP_LIST"]):
                    continue
                #
                hasMSE = False
                hasMET = False
                for seqTupl in valD["SEQ_TUP_LIST"]:
                    if seqTupl[0] == "MSE":
                        hasMSE = True
                    elif seqTupl[0] == "MET":
                        hasMET = True
                    #
                #
                if hasMSE and hasMET:
                    mixMseMetResidueList.append("Entity " + key + " contains both MET and MSE.")
                #
                if ("INSTANCE_LIST" not in valD) or (not valD["INSTANCE_LIST"]):
                    continue
                #
                for instId in valD["INSTANCE_LIST"]:
                    if instId not in instanceD:
                        continue
                    #
                    if (10 * len(instanceD[instId])) > (7 * len(valD["SEQ_TUP_LIST"])):
                        continue
                    #
                    percent = float((len(valD["SEQ_TUP_LIST"]) - len(instanceD[instId])) * 100) / float(len(valD["SEQ_TUP_LIST"]))
                    missingResidueInfoList.append("%.1f" % percent + "% residues of chain '" + instId + "' (%d/%d residues) are missing in coordinates."
                                                  % (len(valD["SEQ_TUP_LIST"]) - len(instanceD[instId]), len(valD["SEQ_TUP_LIST"])))
                #
            #
            return True, no_coord_issue, seqSimilarityInfoList, polymerLinkDistList, missingResidueInfoList, mixMseMetResidueList, entryD, entityD, instanceD, \
                statisticsMap, depSeqAssign, seqAssign
        except:  # noqa: E722 pylint: disable=bare-except
            if self._verbose:
                traceback.print_exc(file=self._lfh)
            #
        #
        return False, False, [], [], [], [], {}, {}, {}, {}, {}, {}

    def __getAuthDefinedRefD(self, depSeqAssign):
        """Get author provided reference information from pdbx_struct_ref_seq_depositor_info category"""
        if not depSeqAssign:
            return
        #
        for key, valD in depSeqAssign.items():
            if ("ref_list" not in valD) or (not valD["ref_list"]):
                continue
            #
            entity_id = key.strip()
            codeList = []
            refList = []
            for refD in valD["ref_list"]:
                data = []
                for item in ("db_name", "db_accession", "db_code", "seq_align_begin", "seq_align_end", "db_align_end",
                             "db_align_beg", "db_seq_one_letter_code"):
                    if (item in refD) and refD[item]:
                        data.append(refD[item].strip().upper())
                    else:
                        data.append("")
                    #
                #
                if (data[0] == "") or ((data[1] == "") and (data[2] == "")):
                    continue
                #
                for idx in (1, 2):
                    val = data[idx]
                    if val and (val not in codeList):
                        codeList.append(val)
                    #
                #
                refList.append(data)
            #
            if len(refList) == 0:
                continue
            elif len(refList) > 1:
                # Chimera case: check the "seq_align_begin", "seq_align_end" values
                missingSeqAlignRange = False
                for refData in refList:
                    if (refData[3] == "") or (refData[4] == ""):
                        missingSeqAlignRange = True
                        break
                    #
                #
                if missingSeqAlignRange:
                    continue
                #
            #
            self.__authDefinedRefD[entity_id] = (codeList, refList)
        #

    def __renderPolymerLinkages(self, polymerLinkDistList, polyLinkHtmlPath):
        """Create HTML table rendering of the input polymer linkage distance list.

        HTML output is written to the "polyLinkHtmlPath":

        Returns True for success of False otherwise.
        """
        try:
            pld = PolymerLinkageDepict(upperBound=1.75, lowerBound=1.2, verbose=self._verbose, log=self._lfh)
            oL = pld.buildPolymerLinkageTable(polymerLinkDistList)
            ofh = open(polyLinkHtmlPath, "w")
            ofh.write("%s" % "".join(oL))
            ofh.close()
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            self._lfh.write("+SequenceDataAssemble.__renderPolymerLinkage() - length linkage list %d  html path %s\n" % (len(oL), polyLinkHtmlPath))
            traceback.print_exc(file=self._lfh)
        #
        return False

    def __checkIfPartialAnnotatedResultExist(self):
        """ """
        partialFilePath = self.__pI.getFilePath(dataSetId=self.__dataSetId, contentType="partial-seq-annotate", formatType="txt")
        if (not partialFilePath) or (not os.access(partialFilePath, os.R_OK)):
            return "", {}
        #
        fph = open(partialFilePath, "r")
        inData = fph.read()
        fph.close()
        #
        prevModelFilePath = ""
        prevSavedFilePathMap = {}
        for line in inData.split("\n"):
            tupL = line.split(" ")
            if len(tupL) == 4 and os.access(tupL[3], os.R_OK):
                if tupL[0] == "model":
                    prevModelFilePath = tupL[3]
                    continue
                #
                if tupL[0] in prevSavedFilePathMap:
                    prevSavedFilePathMap[tupL[0]][tupL[2]] = tupL
                else:
                    myD = {}
                    myD[tupL[2]] = tupL
                    prevSavedFilePathMap[tupL[0]] = myD
                #
            #
        #
        if (prevModelFilePath == "") or (not prevSavedFilePathMap) or (not os.access(prevModelFilePath, os.R_OK)):
            return "", {}
        #
        return prevModelFilePath, prevSavedFilePathMap

    def __getPartialAnnotatedResult(self, prevModelFilePath, prevSavedFilePathMap, currEntityD, currInstanceD):
        """ """
        if prevModelFilePath == "":
            return {}, [], []
        #
        if ("seq-data-stats" not in prevSavedFilePathMap) or ("1" not in prevSavedFilePathMap["seq-data-stats"]):
            return {}, [], []
        #
        tmpPolyLinkPath = os.path.join(self.__sessionPath, self.__dataSetId + "-poly-link.json")
        self.__calcPolymerLinkages(pdbxFilePath=prevModelFilePath, polyLinkPath=tmpPolyLinkPath)
        #
        ok, _no_coord_issue, _seqSimilarityInfoList, _polymerLinkDistList, _missingResidueInfoList, _mixMseMetResidueList, _entryD, prevEntityD, prevInstanceD, _statisticsMap, \
            _depSeqAssign, _seqAssign = self.__readJsonObject(tmpPolyLinkPath)
        #
        if not ok:
            return {}, [], []
        #
        newOldEntityIdMap = {}
        reverseMap = {}
        for newId, newD in currEntityD.items():
            for oldId, oldD in prevEntityD.items():
                sameSeqFlag, sameEntityFlag = self.__isSameEntity(newD, oldD)
                if not sameSeqFlag:
                    continue
                #
                # Copying blast search result pickle file
                #
                if ("seqdb-match" in prevSavedFilePathMap) and (oldId in prevSavedFilePathMap["seqdb-match"]):
                    tupL = prevSavedFilePathMap["seqdb-match"][oldId]
                    importFilePath = self.__pI.getFilePath(dataSetId=self.__dataSetId, contentType=tupL[0], formatType=tupL[1], fileSource="session", partNumber=newId)
                    shutil.copyfile(tupL[3], importFilePath)
                #
                if not sameEntityFlag:
                    break
                #
                foundCoordSeqMisMatch = False
                for cId in newD["INSTANCE_LIST"]:
                    if (cId not in prevInstanceD) or (cId not in currInstanceD):
                        foundCoordSeqMisMatch = True
                        break
                    #
                    if not self.__isSameSequence(currInstanceD[cId], prevInstanceD[cId]):
                        foundCoordSeqMisMatch = True
                        break
                    #
                #
                if foundCoordSeqMisMatch:
                    break
                #
                # Copy sequence alignment pickle file
                #
                if ("seq-align-data" in prevSavedFilePathMap) and (oldId in prevSavedFilePathMap["seq-align-data"]):
                    tupL = prevSavedFilePathMap["seq-align-data"][oldId]
                    importFilePath = self.__pI.getFilePath(dataSetId=self.__dataSetId, contentType=tupL[0], formatType=tupL[1], fileSource="session", partNumber=newId)
                    shutil.copyfile(tupL[3], importFilePath)
                    if newId != oldId:
                        alignUtil = AlignmentDataStore(reqObj=self._reqObj, entityId=newId, pathInfo=self.__pI, verbose=self._verbose, log=self._lfh)
                        alignUtil.updateEntityId(oldId)
                        alignUtil.serialize()
                    #
                    newOldEntityIdMap[newId] = oldId
                    reverseMap[oldId] = newId
                #
                break
            #
        #
        misMatchList = []
        notFoundMatchList = []
        #
        if newOldEntityIdMap:
            self.mergePartialAnnotatedResult(newOldEntityIdMap, prevSavedFilePathMap["seq-data-stats"]["1"][3])
            if ("mismatch-warning" in prevSavedFilePathMap) and ("1" in prevSavedFilePathMap["mismatch-warning"]):
                tupL = prevSavedFilePathMap["mismatch-warning"]["1"]
                try:
                    fb = open(tupL[3], "rb")
                    warningD = pickle.load(fb)
                    fb.close()
                    #
                    if "mismatch" in warningD:
                        for entityId in warningD["mismatch"]:
                            if entityId in reverseMap:
                                misMatchList.append(reverseMap[entityId])
                            #
                        #
                    #
                    if "not_found_existing_match" in warningD:
                        for entityId in warningD["not_found_existing_match"]:
                            if entityId in reverseMap:
                                notFoundMatchList.append(reverseMap[entityId])
                            #
                        #
                    #
                except:  # noqa: E722 pylint: disable=bare-except
                    if self._verbose:
                        traceback.print_exc(file=self._lfh)
                    #
                #
            #
        #
        return newOldEntityIdMap, misMatchList, notFoundMatchList

    def __isSameEntity(self, currED, prevED):
        """ """
        sameSeqFlag = False
        if ("SEQ_TUP_LIST" in currED) and ("SEQ_TUP_LIST" in prevED):
            sameSeqFlag = self.__isSameSequence(currED["SEQ_TUP_LIST"], prevED["SEQ_TUP_LIST"])
        #
        if not sameSeqFlag:
            return sameSeqFlag, sameSeqFlag
        #
        if ("INSTANCE_LIST" not in currED) or ("INSTANCE_LIST" not in prevED):
            return sameSeqFlag, False
        #
        currInstList = sorted(set(currED["INSTANCE_LIST"]))
        prevInstList = sorted(set(prevED["INSTANCE_LIST"]))
        if currInstList != prevInstList:
            return sameSeqFlag, False
        #
        return sameSeqFlag, True

    def __isSameSequence(self, seqTupList1, seqTupList2):
        """ """
        if len(seqTupList1) != len(seqTupList2):
            return False
        #
        for i in range(0, len(seqTupList1)):
            if seqTupList1[i][0] != seqTupList2[i][0]:
                return False
            #
        #
        return True

    def __doReferenceSearch(self, entityD, entityIdList, withRefFlag=False):
        """Perform the reference sequence database search using the input entity dictionary.
        Store matching results in the local session directory.
        """
        try:
            startTime = time.time()
            #
            entityTupList = []
            for eId in entityIdList:
                if eId not in entityD:
                    continue
                #
                eD = entityD[eId]
                sId = eD["ENTITY_ID"]
                polyTypeCode = eD["POLYMER_TYPE_CODE"]
                seqLen = len(eD["SEQ_ENTITY_1_CAN"])
                #
                skip = (
                    (polyTypeCode in ["SAC"])
                    or ((polyTypeCode in ["DNA", "RNA", "XNA"]) and (seqLen < self.__minSearchSequenceLengthNA))
                    or ((polyTypeCode in ["AA"]) and (seqLen < self.__minSearchSequenceLengthAA))
                )
                #
                NA_flag = False
                if polyTypeCode == "DNA":
                    NA_flag = True
                #
                if skip and (polyTypeCode == "AA"):
                    oneLetterCodeSeq = str(eD["SEQ_ENTITY_1_CAN"]).strip().upper()
                    countALA = 0
                    for oneLetterCode in oneLetterCodeSeq:
                        if oneLetterCode == "A":
                            countALA += 1
                        #
                    #
                    # Stop autoProcess for poly-ALA Seq per ticket# DAOTHER-8256
                    if (10 * countALA) > (9 * seqLen):
                        self.__autoProcessFlag = False
                    #
                #
                if self._verbose:
                    self._lfh.write("+SequenceDataAssemble.__doReferenceSearch() search for entity id %s type %s length %d skip status %r\n" % (eId, polyTypeCode, seqLen, skip))
                    self._lfh.flush()
                #
                if sId in self.__authDefinedRefD:
                    entityTupList.append((sId, eD, self.__authDefinedRefD[sId], skip, NA_flag))
                else:
                    entityTupList.append((sId, eD, (), skip, NA_flag))
                #
            #
            if not entityTupList:
                if self._verbose:
                    endTime = time.time()
                    self._lfh.write(
                        "+SequenceDataAssemble.__doReferenceSearch()  completed at %s (%.2f seconds)\n"
                        % (time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
                    )
                #
                return {}, {}, {}
            #
            if withRefFlag and (len(entityTupList) == 1):
                seqSearchUtil = LocalBlastSearchUtils(siteId=self.__siteId, sessionPath=self.__sessionPath, pathInfo=self.__pI,
                                                      doRefSearchFlag=self.__doRefSearchFlag, verbose=self._verbose, log=self._lfh)
                predefinedRefD = seqSearchUtil.getPredefinedRefSequence(dataSetId=self.__dataSetId, entityD=entityTupList[0][1])
                if predefinedRefD:
                    return predefinedRefD, {}, {}
                #
            #
            ownRefD, otherRefD = self.__runSameSeqAnnotationSearch(entityTupList)
            #
            blastList = []
            skipList = []
            naList = []
            if self.__reRunBlastSearchFlag:
                blastList = entityTupList
                self.__autoProcessFlag = False
            else:
                for entityTup in entityTupList:
                    hasOwnRef = False
                    if (entityTup[0] in ownRefD) and ("auto_match_status" in ownRefD[entityTup[0]]) and ownRefD[entityTup[0]]["auto_match_status"]:
                        hasOwnRef = True
                    #
                    hasSameSeqRef = False
                    if (entityTup[0] in otherRefD) and ("auto_match_status" in otherRefD[entityTup[0]]) and otherRefD[entityTup[0]]["auto_match_status"]:
                        hasSameSeqRef = True
                    #
                    if hasOwnRef:
                        continue
                    elif hasSameSeqRef:
                        if self.__forceBlastSearchFlag and (not entityTup[3]):
                            blastList.append(entityTup)
                        #
                        continue
                    #

                    if entityTup[3]:
                        skipList.append(entityTup)
                    #
                    if (len(entityTup[2]) > 0) or (not entityTup[3]):
                        blastList.append(entityTup)
                        if not entityTup[3]:
                            if entityTup[4]:
                                naList.append(entityTup[0])
#                           else:
#                               self.__autoProcessFlag = False
                            #
                        #
                    #
                #
            #
            refD = {}
            if blastList:
                autoMatchStatus, refD = self.__runSeqBlastSearch(blastList)
                if (not autoMatchStatus) or (not refD):
                    self.__autoProcessFlag = False
                #
            #
            if refD:
                for k, _v in refD.items():
                    if k in naList:
                        self.__autoProcessFlag = False
                    #
                #
            #
            if self.__autoProcessFlag:
                for entityTup in skipList:
                    if entityTup[0] in refD:
                        continue
                    #
                    #                   selfRefmap = {}
                    #                   selfRefmap["1"] = True
                    #                   ownRefD[entityTup[0]] = selfRefmap
                    ownRefD[entityTup[0]] = {"1": {"selfref": True}}
                #
            if self._verbose:
                endTime = time.time()
                self._lfh.write(
                    "+SequenceDataAssemble.__doReferenceSearch()  completed at %s (%.2f seconds)\n" % (time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
                )
            #
            return refD, ownRefD, otherRefD
        except:  # noqa: E722 pylint: disable=bare-except
            if self._verbose:
                self._lfh.write("+SequenceDataAssemble.__doReferenceSearch() - failing \n")
                traceback.print_exc(file=self._lfh)
            #
        #
        return {}, {}, {}

    def __runSameSeqAnnotationSearch(self, entityTupList):
        """Search same sequence annotation information from processed entries"""
        try:
            startTime = time.time()
            #
            if self.__cI.get("USE_COMPUTE_CLUSTER"):
                numProc = len(entityTupList)
            else:
                numProc = int(multiprocessing.cpu_count() / 2)
            mpu = MultiProcUtil(verbose=True)
            mpu.set(workerObj=self, workerMethod="runMultiSameSeqAnnotationSearches")
            mpu.setWorkingDir(self.__sessionPath)
            #
            _ok, _failList, retLists, _diagList = mpu.runMulti(dataList=entityTupList, numProc=numProc, numResults=1)
            #
            ownRefD = {}
            otherRefD = {}
            for tupList in retLists:
                for tup in tupList:
                    if tup[1]:
                        ownRefD[tup[0]] = tup[1]
                    #
                    if tup[2]:
                        otherRefD[tup[0]] = tup[2]
                    #
                #
            #
            if self._verbose:
                endTime = time.time()
                self._lfh.write(
                    "+SequenceDataAssemble.__runSameSeqAnnotationSearch()  completed at %s (%.2f seconds)\n"
                    % (time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
                )
            #
            return ownRefD, otherRefD
        except:  # noqa: E722 pylint: disable=bare-except
            if self._verbose:
                self._lfh.write("+SequenceDataAssemble.__runSameSeqAnnotationSearch() - failing \n")
                traceback.print_exc(file=self._lfh)
            #
        #
        return {}, {}

    def runMultiSameSeqAnnotationSearches(self, dataList, procName, optionsD, workingDir):  # pylint: disable=unused-argument
        """Multiple same sequence annotation search processing API"""
        rList = []
        for tupL in dataList:
            seqSearchUtil = SeqAnnotationSearchUtils(siteId=self.__siteId, sessionPath=self.__sessionPath, pathInfo=self.__pI, verbose=self._verbose, log=self._lfh)
            selfRefD, sameSeqRefD = seqSearchUtil.getSameSeqRefInfo(self.__dataSetId, tupL[1], tupL[2])
            rList.append((tupL[0], selfRefD, sameSeqRefD))
        #
        return rList, rList, []

    def __runSeqBlastSearch(self, entityTupList):
        """Perform blast reference sequence search"""
        try:
            startTime = time.time()
            #
            if self.__cI.get("USE_COMPUTE_CLUSTER"):
                numProc = len(entityTupList)
            else:
                numProc = int(multiprocessing.cpu_count() / 2)
            mpu = MultiProcUtil(verbose=True)
            mpu.set(workerObj=self, workerMethod="runMultiBlastReferenceSearches")
            mpu.setWorkingDir(self.__sessionPath)
            #
            # Setup limits for NCBI requets
            apikey = self.__cI.get("NCBI_API_KEY")
            apirate = self.__cI.get("NCBI_API_RATE")
            if apikey:
                if apirate:
                    rate = int(apirate)
                else:
                    rate = 5
                #
            else:
                rate = 1
            #
            mpl = MultiProcLimit(rate)
            mpu.setOptions({"ncbilock": mpl})
            #
            _ok, _failList, retLists, _diagList = mpu.runMulti(dataList=entityTupList, numProc=numProc, numResults=1)
            #
            autoMatchStatus = True
            refD = {}
            for tupList in retLists:
                for tup in tupList:
                    if tup[1]:
                        refD[tup[0]] = tup[1]
                        if (len(tup) != 3) or (not tup[2]):
                            autoMatchStatus = False
                        #
                    else:
                        autoMatchStatus = False
                    #
                #
            #
            if not refD:
                autoMatchStatus = False
            #
            if self._verbose:
                endTime = time.time()
                self._lfh.write(
                    "+SequenceDataAssemble.__runSeqBlastSearch()  completed at %s (%.2f seconds)\n" % (time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
                )
            #
            return autoMatchStatus, refD
        except:  # noqa: E722 pylint: disable=bare-except
            if self._verbose:
                self._lfh.write("+SequenceDataAssemble.__runSeqBlastSearch() - failing \n")
                traceback.print_exc(file=self._lfh)
            #
        #
        return False, {}

    def runMultiBlastReferenceSearches(self, dataList, procName, optionsD, workingDir):  # pylint: disable=unused-argument
        """Multiple blast reference sequence search processing API"""
        ncbilock = optionsD.get("ncbilock", None)
        rList = []
        for tupL in dataList:
            seqSearchUtil = LocalBlastSearchUtils(
                siteId=self.__siteId,
                sessionPath=self.__sessionPath,
                pathInfo=self.__pI,
                doRefSearchFlag=self.__doRefSearchFlag,
                verbose=self._verbose,
                log=self._lfh,
                ncbilock=ncbilock,
            )
            autoMatchStatus, eRefD = seqSearchUtil.searchSeqReference(dataSetId=self.__dataSetId, entityD=tupL[1], authRefList=tupL[2])
            rList.append((tupL[0], eRefD, autoMatchStatus))
            #
        #
        return rList, rList, []

    def __createSeqMatchFileAndMisMatchPickleFile(self, sortedEntityIdList, entityD, ownRefD, eSSRefD, prevEntityIdList, prevMisMatchList,
                                                  prevNotFoundMatchList, seqSimilarityInfoList, missingResidueInfoList, mixMseMetResidueList):
        """Create D_xxxxxxxxxx_seqmatch_P1.cif.Vx file"""
        annList = []
        misMatchList = []
        notFoundMatchList = []
        for entityId in sortedEntityIdList:
            if entityId in prevEntityIdList:
                if entityId in prevMisMatchList:
                    misMatchList.append(entityId)
                #
                if entityId in prevNotFoundMatchList:
                    notFoundMatchList.append(entityId)
                #
                continue
            #
            annTuple = []
            annTuple.append(entityId)
            if (
                (entityId in ownRefD)
                and ("1" in ownRefD[entityId])
                and ("REF_ENTRY_ID" in ownRefD[entityId]["1"])
                and ("REF_PDB_ID" in ownRefD[entityId]["1"])
                and ("REF_ENTRY_ENTITY_ID" in ownRefD[entityId]["1"])
                and ("db_name" in ownRefD[entityId]["1"])
                and ("db_code" in ownRefD[entityId]["1"])
            ):
                annTuple.append("Y")
                annTuple.append(ownRefD[entityId]["1"]["REF_ENTRY_ID"])
                annTuple.append(ownRefD[entityId]["1"]["REF_PDB_ID"])
                annTuple.append(ownRefD[entityId]["1"]["REF_ENTRY_ENTITY_ID"])
                annTuple.append(ownRefD[entityId]["1"]["db_name"])
                annTuple.append(ownRefD[entityId]["1"]["db_code"])
            elif (
                entityId in eSSRefD
                and ("1" in eSSRefD[entityId])
                and ("REF_ENTRY_ID" in eSSRefD[entityId]["1"])
                and ("REF_PDB_ID" in eSSRefD[entityId]["1"])
                and ("REF_ENTRY_ENTITY_ID" in eSSRefD[entityId]["1"])
                and ("db_name" in eSSRefD[entityId]["1"])
                and ("db_code" in eSSRefD[entityId]["1"])
            ):
                annTuple.append("Y")
                annTuple.append(eSSRefD[entityId]["1"]["REF_ENTRY_ID"])
                annTuple.append(eSSRefD[entityId]["1"]["REF_PDB_ID"])
                annTuple.append(eSSRefD[entityId]["1"]["REF_ENTRY_ENTITY_ID"])
                annTuple.append(eSSRefD[entityId]["1"]["db_name"])
                annTuple.append(eSSRefD[entityId]["1"]["db_code"])
            else:
                annTuple.append("N")
                annTuple.append(".")
                annTuple.append(".")
                annTuple.append(".")
                annTuple.append(".")
                annTuple.append(".")
                notFoundMatchList.append(entityId)
            #
            if "SEQ_COORD_MISMATCH" in entityD[entityId]:
                annTuple.append("Y")
                annTuple.append("\n".join(entityD[entityId]["SEQ_COORD_MISMATCH"]))
                misMatchList.append(entityId)
                #
            else:
                annTuple.append("N")
                annTuple.append(".")
            #
            annList.append(annTuple)
        #
        warningMsgMap = {}
        if len(seqSimilarityInfoList) > 0:
            warningMsgMap["seq_warning_info"] = "</br>\n".join(seqSimilarityInfoList)
        #
        warningMsgMap["mismatch"] = misMatchList
        warningMsgMap["not_found_existing_match"] = notFoundMatchList
        warningMsgMap["missing_residue"] = missingResidueInfoList
        warningMsgMap["mix_mse_met"] = mixMseMetResidueList
        #
        fb = open(self.__misMatchPath, "wb")
        pickle.dump(warningMsgMap, fb)
        fb.close()
        #
        if prevEntityIdList:
            return
        #
        try:
            matchPath = self.__pI.getFilePath(self.__dataSetId, contentType="seqmatch", formatType="pdbx", versionId="next")
            ofh = open(matchPath, "w")
            #
            dataBlock = DataContainer(self.__dataSetId)
            #
            aCat = DataCategory("entry")
            aCat.appendAttribute("id")
            aCat.setValue(self.__dataSetId, "id", 0)
            dataBlock.append(aCat)
            #
            dataBlock.append(
                self.__writeGeneralCifCategory(
                    "pdbx_sequence_annotation",
                    (
                        "entity_id",
                        "entity_matched",
                        "matched_entry_id",
                        "matched_pdb_id",
                        "matched_entity_id",
                        "reference_db_name",
                        "reference_db_code",
                        "coord_seq_mismatch",
                        "coord_seq_mismatch_message",
                    ),
                    annList,
                )
            )
            #
            blockList = []
            blockList.append(dataBlock)
            pdbxW = PdbxWriter(ofh)
            pdbxW.write(blockList)
            ofh.close()
        except:  # noqa: E722 pylint: disable=bare-except
            if self._verbose:
                traceback.print_exc(file=self._lfh)
            #
        #

    def __writeGeneralCifCategory(self, categoryName, categoryItems, dataList):
        """Write cif category`"""
        aCat = DataCategory(categoryName)
        for itemName in categoryItems:
            aCat.appendAttribute(itemName)
        #
        for row in dataList:
            aCat.append(row)
        #
        return aCat

    def __updateUtilDataStore(self, pdbxFilePath):
        uds = UtilDataStore(reqObj=self._reqObj, verbose=self._verbose, log=self._lfh)
        pdbxFileName = os.path.basename(pdbxFilePath)
        uds.set("modelfilename", pdbxFileName)
        uds.serialize()
        #
        # Make an extra copy of the CIF file with a standard file extension  --  (may not be needed)
        #
        dstPath = os.path.join(self.__sessionPath, self.__dataSetId + ".cif")
        try:
            shutil.copyfile(pdbxFilePath, dstPath)
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self._lfh)
        #

    def __createDefaultAlignments(self):
        """Generate all default selected sequence alignments for all entities"""
        if not self._entityAlignInfoMap:
            return
        #
        #
        # try:
        #     numProc = multiprocessing.cpu_count() / 2
        #     mpu = MultiProcUtil(verbose = True)
        #     mpu.set(workerObj = self, workerMethod = "multiCreateAlignmentProcess")
        #     mpu.setWorkingDir(self.__sessionPath)
        #     ok,failList,retLists,diagList = mpu.runMulti(dataList = self._entityAlignInfoMap.values(), numProc = numProc, numResults = 1)
        #     #
        # except:  # noqa: E722 pylint: disable=bare-except
        #     if (self._verbose):
        #         traceback.print_exc(file=self._lfh)
        #     #
        # #
        #
        for dataD in self._entityAlignInfoMap.values():
            self.__createDefaultAlignmentObject(dataD)
        #

    def multiCreateAlignmentProcess(self, dataList, procName, optionsD, workingDir):  # pylint: disable=unused-argument
        """Miltiple auth vs coord vs ref sequences alignment processing API"""
        rList = []
        for dataD in dataList:
            self.__createDefaultAlignmentObject(dataD)
            rList.append(dataD["id"])
        #
        return rList, rList, []

    def __createDefaultAlignmentObject(self, dataD):
        """Generate a default selected sequence alignment for a entity"""
        alignUtil = AlignmentTools(
            reqObj=self._reqObj, entityId=dataD["id"], pathInfo=self.__pI, seqDataStore=self.getSequenceDataStoreObj(), deserializeFlag=False, verbose=self._verbose, log=self._lfh
        )
        if self.__resetAlignmentFlag:
            alignUtil.resetInputAlignData(alignD=dataD)
        else:
            alignUtil.setInputAlignData(alignD=dataD)
        #

    def __reRunSearch(self, selectedEntityIdList, withRefFlag=False):
        """ """
        self.__doRefSearchFlag = True
        self.__forceBlastSearchFlag = True
        #
        startTime = time.time()
        if self._verbose:
            self._lfh.write(
                "+SequenceDataAssemble.__reRunSearch() for sessionId %s started at %s \n" % (self.__sessionObj.getId(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))
            )
        #
        self.__getAuthDefinedRefD(self.getDepositorSeqAssign())
        #
        entityD = {}
        for entityId in selectedEntityIdList:
            authLabel, eD = self.getEntityDetails(entityId)
            if eD:
                entityD[entityId] = eD
                self.__setEntityAlignInfoMap(authLabel, eD)
            #
        #
        eRefD, ownRefD, eSSRefD = self.__doReferenceSearch(entityD, list(entityD.keys()), withRefFlag=withRefFlag)
        self.setMultipleEntitiesRefSequenceInfo(entityD, eRefD, ownRefD, eSSRefD)
        #
        self.__resetAlignmentFlag = True
        self.__createDefaultAlignments()
        #
        if self._verbose:
            endTime = time.time()
            self._lfh.write("+SequenceDataAssemble.__reRunSearch() completed at %s (%.2f seconds)\n" % (time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime))
        #

    def __updateSelections(self, selectedEntityIdList):
        """Update the sequence selections"""
        if not self._entityAlignInfoMap:
            self.__getCurrentEntityAlignInfoMap(selectedEntityIdList)
        #
        removeList = []
        for entityId in selectedEntityIdList:
            if (entityId not in self._entityAlignInfoMap) or ("alignids" not in self._entityAlignInfoMap[entityId]):
                continue
            #
            removeList.append("auth_" + entityId)
            removeList.append("ref_" + entityId)
            removeList.append("selfref_" + entityId)
            #
            pdbChainIdList = self.getGroup(entityId)
            for chainId in pdbChainIdList:
                removeList.append("xyz_" + chainId)
            #
        #
        updatedSelectIdList = []
        selectIdList = self.getSelectedIds()
        for seqId in selectIdList:
            tL = str(seqId).split("_")
            key = tL[0] + "_" + tL[1]
            if key in removeList:
                continue
            #
            updatedSelectIdList.append(seqId)
        #
        for entityId, entityAlignInfo in self._entityAlignInfoMap.items():
            for seqId in entityAlignInfo["alignids"]:
                updatedSelectIdList.append(seqId)
            #
        #
        self.setSelectedIds(updatedSelectIdList)

    def __setEntityAlignInfoMap(self, authLabel, eD):
        """Set self._entityAlignInfoMap for coordinate parts"""
        sId = eD["ENTITY_ID"]
        #
        alignids = []
        xyz_label = []
        alignids.append(authLabel)
        xyz_label.append(authLabel)
        #
        entityAlignInfo = {}
        entityAlignInfo["id"] = sId
        entityAlignInfo["auth_label"] = authLabel
        #
        seqLabel = SequenceLabel(verbose=self._verbose)
        for cId in eD["INSTANCE_LIST"]:
            seqLabel.set(seqType="xyz", seqInstId=cId, seqPartId=1, seqAltId=1, seqVersion=1)
            alignids.append(seqLabel.pack())
            xyz_label.append(seqLabel.pack())
        #
        entityAlignInfo["alignids"] = alignids
        entityAlignInfo["xyz_label"] = xyz_label
        #
        partInfo = {}
        pDList = eD["PART_LIST"]
        for (partNo, pD) in enumerate(pDList, start=1):
            partInfo[partNo] = (int(pD["SEQ_NUM_BEG"]), int(pD["SEQ_NUM_END"]))
        #
        entityAlignInfo["part_info"] = partInfo
        self._entityAlignInfoMap[sId] = entityAlignInfo

    def __getCurrentEntityAlignInfoMap(self, selectedEntityIdList):
        """ """
        for entityId in selectedEntityIdList:
            authId = self.getLatestVersionSeqId(seqType="auth", seqInstId=entityId)
            if not authId:
                continue
            #
            alignids = []
            alignids.append(authId)
            #
            pdbChainIdList = self.getGroup(entityId)
            for chainId in pdbChainIdList:
                xyzId = self.getLatestVersionSeqId(seqType="xyz", seqInstId=chainId)
                if xyzId:
                    alignids.append(xyzId)
                #
            #
            partIdList = self.getPartIdList("auth", entityId)
            for partId in partIdList:
                altIdList = self.getAlternativeIdList("ref", entityId, partId)
                if len(altIdList) == 0:
                    alignids.append("selfref_" + str(entityId) + "_" + str(partId))
                    continue
                #
                refId = self.getLatestVersionSeqId(seqType="ref", seqInstId=entityId, seqPartId=partId, seqAltId=altIdList[0])
                if refId:
                    alignids.append(refId)
                else:
                    alignids.append("selfref_" + str(entityId) + "_" + str(partId))
                #
            #
            entityAlignInfo = {}
            entityAlignInfo["alignids"] = alignids
            self._entityAlignInfoMap[entityId] = entityAlignInfo
        #
