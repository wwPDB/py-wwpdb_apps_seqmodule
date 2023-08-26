##
# File:  UpdateSequenceDataStoreUtils.py
# Date:  23-Nov-2018
#
# Updates:
#  03-Oct-2022 zf  improve the handling for author provided reference sequence information
#
##
"""
Wrap class to update SequenceDataStore storage

"""
__docformat__ = "restructuredtext en"
__author__ = "Zukang Feng"
__email__ = "zfeng@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.09"

import copy
import os
import sys
import traceback

from wwpdb.apps.seqmodule.io.SequenceDataStore import SequenceDataStore
from wwpdb.apps.seqmodule.util.SequenceLabel import SequenceLabel, SequenceFeature
from wwpdb.io.file.mmCIFUtil import mmCIFUtil
from wwpdb.io.locator.ChemRefPathInfo import ChemRefPathInfo


class UpdateSequenceDataStoreUtils(object):
    """ """

    def __init__(self, reqObj=None, seqDataStore=None, verbose=False, log=sys.stderr):
        """ """
        self._reqObj = reqObj
        self.__sds = seqDataStore
        self._verbose = verbose
        self._lfh = log
        #
        self.__dataSetId = str(self._reqObj.getValue("identifier")).upper()
        self.__pdbId = self.__dataSetId
        #
        self.__seqLabel = SequenceLabel(verbose=self._verbose)
        self.__seqFeature = SequenceFeature(verbose=self._verbose, log=self._lfh)
        #
        self.__polymerTypeCode = {}
        self.__defSelList = []
        self._entityAlignInfoMap = {}
        self.__chainIDSeqLabelMap = {}

    def __checkSequenceDataStore(self):
        """ """
        if not self.__sds:
            self.__sds = SequenceDataStore(reqObj=self._reqObj, verbose=self._verbose, log=self._lfh)
        #

    def resetSequenceDataStore(self):
        """ """
        self.__checkSequenceDataStore()
        self.__sds.reset()

    def saveSequenceDataStore(self):
        """ """
        self.__checkSequenceDataStore()
        self.__sds.serialize()

    def getSequenceDataStoreObj(self):
        """ """
        self.__checkSequenceDataStore()
        return self.__sds

    def dumpSequenceDataStore(self):
        """ """
        self.__checkSequenceDataStore()
        self.__sds.dump(self._lfh)

    def getSequence(self, seqType, seqInstId, seqPartId, seqAltId, seqVersion):
        """Get sequence from SequenceDataStore object"""
        try:
            self.__checkSequenceDataStore()
            return self.__sds.getSequence(seqId=seqInstId, seqType=seqType, partId=seqPartId, altId=seqAltId, version=seqVersion)
        except:  # noqa: E722 pylint: disable=bare-except
            if self._verbose:
                traceback.print_exc(file=self._lfh)
            #
        #
        return []

    def setSequence(self, seqList, seqType, seqInstId, seqPartId=1, seqAltId=1, seqVersion=1):
        """Insert sequence to SequenceDataStore object"""
        self.__checkSequenceDataStore()
        self.__sds.setSequence(seqList, seqId=seqInstId, seqType=seqType, partId=seqPartId, altId=seqAltId, version=seqVersion)

    def getFeature(self, seqType, seqInstId, seqPartId, seqAltId, seqVersion):
        """Get feature from SequenceDataStore object"""
        try:
            self.__checkSequenceDataStore()
            return self.__sds.getFeature(seqId=seqInstId, seqType=seqType, partId=seqPartId, altId=seqAltId, version=seqVersion)
        except:  # noqa: E722 pylint: disable=bare-except
            if self._verbose:
                traceback.print_exc(file=self._lfh)
            #
        #
        return {}

    def setFeature(self, featureD, seqType, seqInstId, seqPartId=1, seqAltId=1, seqVersion=1):
        """Insert feature to SequenceDataStore object"""
        self.__checkSequenceDataStore()
        self.__sds.setFeature(featureD, seqId=seqInstId, seqType=seqType, partId=seqPartId, altId=seqAltId, version=seqVersion)

    def getGroup(self, entityId):
        """Get PDB chain ID group from SequenceDataStore object"""
        self.__checkSequenceDataStore()
        return self.__sds.getGroup(entityId)

    def setGroup(self, groupDict):
        """Set entity ID vs. PDB chain IDs map"""
        self.__checkSequenceDataStore()
        for k, v in groupDict.items():
            self.__sds.setGroup(k, v)
        #

    def getPartIdList(self, seqType, entityId):
        """Get pard Id list from SequenceDataStore object"""
        self.__checkSequenceDataStore()
        return self.__sds.getPartIds(entityId, dataType="sequence", seqType=seqType)

    def getVersionIdList(self, seqType, seqInstId, seqPartId, seqAltId):
        """Get version Id list from SequenceDataStore object"""
        self.__checkSequenceDataStore()
        return self.__sds.getVersionIds(seqInstId, partId=seqPartId, altId=seqAltId, dataType="sequence", seqType=seqType)

    def getNextVersionNumber(self, seqType, seqInstId, seqPartId, seqAltId):
        """Get next version number"""
        existingVersList = self.getVersionIdList(seqType, seqInstId, seqPartId, seqAltId)
        if len(existingVersList) == 0:
            return 1
        else:
            return int(existingVersList[0]) + 1
        #

    def getLatestVersionSeqId(self, seqType="auth", seqInstId="1", seqPartId=1, seqAltId=1):
        """Get latest version of auth/xyz sequence Id"""
        existingVersList = self.getVersionIdList(seqType, seqInstId, seqPartId, seqAltId)
        if len(existingVersList) > 0:
            if len(existingVersList) > 1:
                existingVersList.sort(reverse=True)
            #
            self.__seqLabel.set(seqType=seqType, seqInstId=seqInstId, seqPartId=seqPartId, seqAltId=seqAltId, seqVersion=existingVersList[0])
            return self.__seqLabel.pack()
        #
        return ""

    def getAlternativeIdList(self, seqType, seqInstId, seqPartId):
        """Get alternative Id list from SequenceDataStore object"""
        self.__checkSequenceDataStore()
        return self.__sds.getAlternativeIds(seqInstId, dataType="sequence", seqType=seqType, partId=seqPartId)

    def getNextAlternativeNumber(self, seqType, seqInstId, seqPartId):
        """Get next alternative number"""
        altIdList = self.getAlternativeIdList(seqType, seqInstId, seqPartId)
        if len(altIdList) == 0:
            return 1
        else:
            return int(altIdList[0]) + 1
        #

    def getEntryDetail(self, key):
        """Get entry details from SequenceDataStore object"""
        self.__checkSequenceDataStore()
        return self.__sds.getEntryDetail(detailKey=key)

    def setEntryDetail(self, entryD):
        """Set entry details"""
        self.__checkSequenceDataStore()
        for k, v in entryD.items():
            self.__sds.setEntryDetail(k, v)
        #

    def getEntityDetails(self, entityId):
        """Get entity feature data from the stored author sequence entity sequence and feature data.
        Returns a dictionary of features for the more recent sequence/feature versions.
        """
        entityD = {}
        verList = self.getVersionIdList("auth", entityId, 1, 1)
        if len(verList) == 0:
            return "", entityD
        #
        siteId = str(self._reqObj.getValue("WWPDB_SITE_ID"))
        #
        sTupL = self.getSequence("auth", entityId, 1, 1, verList[0])
        r1L = []
        r1L_CAN = []
        r1L_SEARCH = []
        for sTup in sTupL:
            r1L_CAN.append(sTup[4])
            if sTup[4] == "X":
                r1L.append("(" + sTup[0] + ")")
                #
                ccCode = ""
                ccId = sTup[0].upper()
                crpi = ChemRefPathInfo(siteId=siteId, verbose=self._verbose, log=self._lfh)
                ccFilePath = crpi.getFilePath(ccId, "CC")
                if ccFilePath and os.access(ccFilePath, os.F_OK):
                    cifObj = mmCIFUtil(filePath=ccFilePath)
                    ccCode = cifObj.GetSingleValue("chem_comp", "one_letter_code")
                #
                if ccCode:
                    r1L_SEARCH.append(ccCode)
                else:
                    r1L_SEARCH.append(sTup[4])
                #
            else:
                r1L.append(sTup[4])
                r1L_SEARCH.append(sTup[4])
            #
        #
        self.__seqFeature.set(self.getFeature("auth", entityId, 1, 1, verList[0]))
        #
        entityD["POLYMER_TYPE_CODE"] = self.__seqFeature.getPolymerType()
        entityD["POLYMER_LINKING_TYPE"] = self.__seqFeature.getPolymerLinkingType()
        entityD["ENTITY_DESCRIPTION"] = self.__seqFeature.getEntityDescription()
        entityD["SEQ_ENTITY_1"] = "".join(r1L)
        entityD["SEQ_ENTITY_1_CAN"] = "".join(r1L_CAN)
        entityD["SEQ_ENTITY_1_SEARCH"] = "".join(r1L_SEARCH)
        entityD["ENTITY_ID"] = entityId
        entityD["PART_LIST"] = self.__getPartList(entityId=entityId)
        entityD["IDENTIFIER"] = entityId
        entityD["ENTRY_ID"] = self.__dataSetId
        entityD["INSTANCE_LIST"] = self.getGroup(entityId)
        #
        self.__seqLabel.set(seqType="auth", seqInstId=entityId, seqPartId=1, seqAltId=1, seqVersion=verList[0])
        return self.__seqLabel.pack(), entityD

    def setPdbId(self, pdbId):
        """ """
        self.__pdbId = pdbId

    def setGroupParts(self, groupPartD):
        """ """
        self.__checkSequenceDataStore()
        for k, v in groupPartD.items():
            self.__sds.setGroupParts(k, v)
        #

    def setPolymerTypeCode(self, polymerTypeCode):
        """ """
        self.__polymerTypeCode = polymerTypeCode

    def setDepositorSeqAssign(self, depSeqAssign):
        """ """
        self.__checkSequenceDataStore()
        self.__sds.setDepositorReferenceAssignments(assignD=depSeqAssign)

    def getDepositorSeqAssign(self):
        """ """
        self.__checkSequenceDataStore()
        return self.__sds.getDepositorReferenceAssignments()

    def setArchiveSeqAssign(self, seqAssign):
        """ """
        self.__checkSequenceDataStore()
        self.__sds.setReferenceAssignments(assignD=seqAssign)

    def setSelectedIds(self, selectedIdList):
        """ """
        self.__checkSequenceDataStore()
        self.__sds.setSelectedIds(idList=selectedIdList)

    def getSelectedIds(self):
        """ """
        self.__checkSequenceDataStore()
        return self.__sds.getSelectedIds()

    def setCoordinateInstanceInfo(self, instanceD, statisticsMap, skipList=None):
        """Load coordinate sequences & features  (coordinate sequences have no parts i.e. partId=1)"""
        if skipList is None:
            skipList = []
        for cId, S3L in instanceD.items():
            if cId in skipList:
                continue
            #
            # Only store sequences if we know the polymer type details -
            if cId not in self.__polymerTypeCode:
                continue
            #
            self.__seqLabel.set(seqType="xyz", seqInstId=cId, seqPartId=1, seqAltId=1, seqVersion=1)
            self.__defSelList.append(self.__seqLabel.pack())
            self.__chainIDSeqLabelMap[cId] = self.__seqLabel.pack()
            #
            sTupL = []
            for sList in S3L:
                # ( 0: 3_L_code, 1: auth_numbering, 2: comment, 3: sequential_index (starting 1), 4: 1_L_code, 5: Orig_3_L_codem 6: linkage_flag )
                sTupL.append((sList[0], sList[1], sList[2], int(sList[3]), sList[4], sList[0], sList[5]))
            #
            self.setSequence(sTupL, "xyz", cId)
            self.__seqFeature.clear()
            self.__seqFeature.setId(dbName="PDB", dbCode=self.__pdbId, dbAccession=self.__pdbId)
            self.__seqFeature.setPolymerType(self.__polymerTypeCode[cId])
            if len(S3L) > 0:
                self.__seqFeature.setAuthXyzAlignRange(seqBegin=S3L[0][1], seqEnd=S3L[len(S3L) - 1][1])
            #
            if cId in statisticsMap:
                self.__seqFeature.setAuthXyzAlignDetails(
                    seqLen=len(S3L), alignLen=int(statisticsMap[cId][0]), seqSim=float(statisticsMap[cId][1]), seqSimWithGaps=float(statisticsMap[cId][2])
                )
            #
            self.setFeature(self.__seqFeature.get(), "xyz", cId)
        #

    def setMultipleEntitiesSequenceInfo(self, entityD, skipList=None):
        """Load author/entity sequence and features for multiple entities"""
        if skipList is None:
            skipList = []
        for eId, eD in entityD.items():
            if eId in skipList:
                continue
            #
            self.setSingleEntitySequenceInfo(eId, eD)
        #

    def setSingleEntitySequenceInfo(self, eId, eD):  # pylint: disable=unused-argument
        """Load author/entity sequence and features for single entity - for multipart entity, distinguishing sequence features are stored for each part."""
        sId = eD["ENTITY_ID"]
        #
        self.__seqLabel.set(seqType="auth", seqInstId=sId, seqPartId=1, seqAltId=1, seqVersion=1)
        self.__defSelList.append(self.__seqLabel.pack())
        #
        alignids = []
        xyz_label = []
        alignids.append(self.__seqLabel.pack())
        xyz_label.append(self.__seqLabel.pack())
        #
        entityAlignInfo = {}
        entityAlignInfo["id"] = sId
        entityAlignInfo["auth_label"] = self.__seqLabel.pack()
        for cId in eD["INSTANCE_LIST"]:
            if cId in self.__chainIDSeqLabelMap:
                alignids.append(self.__chainIDSeqLabelMap[cId])
                xyz_label.append(self.__chainIDSeqLabelMap[cId])
            #
        #
        entityAlignInfo["alignids"] = alignids
        entityAlignInfo["xyz_label"] = xyz_label
        entityAlignInfo["auth_coord_align_index"] = eD["ALIGN_INDEX"]
        #
        sTupL = []
        for sList in eD["SEQ_TUP_LIST"]:
            # ( 0: 3_L_code, 1: sequential_numbering (starting 1), 2: comment, 3: sequential_index (starting 1), 4: 1_L_code, 5: Orig_3_L_code )
            sTupL.append((sList[0], sList[1], sList[2], int(sList[1]), sList[3], sList[0]))
        #
        partInfo = {}
        pDList = eD["PART_LIST"]
        for (partNo, pD) in enumerate(pDList, start=1):
            partInfo[partNo] = (int(pD["SEQ_NUM_BEG"]), int(pD["SEQ_NUM_END"]))
            #
            self.setSequence(sTupL, "auth", sId, seqPartId=partNo)
            self.setFeature(self.getAuthFeatureObj(self.__pdbId, eD, partNo, pD).get(), "auth", sId, seqPartId=partNo)
        #
        entityAlignInfo["part_info"] = partInfo
        self._entityAlignInfoMap[sId] = entityAlignInfo

    def setMultipleEntitiesRefSequenceInfo(self, entityD, eRefD, ownRefD, eSSRefD, skipList=None):
        """Load reference sequence and features for multiple entities"""
        if skipList is None:
            skipList = []
        for eId, eD in entityD.items():
            if eId in skipList:
                continue
            #
            self.setSingleEntityRefSequenceInfo(eId, eD, eRefD, ownRefD, eSSRefD)
        #

    def setSingleEntityRefSequenceInfo(self, eId, eD, eRefD, ownRefD, eSSRefD):
        """Load reference sequence(s) and features for single entity"""
        ref_label = {}
        ref_align_index = {}
        sId = eD["ENTITY_ID"]
        for (partNo, pD) in enumerate(eD["PART_LIST"], start=1):
            rList = []
            autoList = []
            selfRefFlag = False
            if (eId in eRefD) and eRefD[eId] and (str(partNo) in eRefD[eId]) and eRefD[eId][str(partNo)]:
                rList = eRefD[eId][str(partNo)]
#               for refD in eRefD[eId][str(partNo)]:
#                   if ("auto_match_status" in refD) and refD["auto_match_status"]:
#                       autoList.append(refD)
#                   else:
#                       rList.append(refD)
#                   #
#               #
            #
            if (eId in eSSRefD) and eSSRefD[eId] and (str(partNo) in eSSRefD[eId]) and eSSRefD[eId][str(partNo)] and ("selfref" not in eSSRefD[eId][str(partNo)]):
                # rList.append(eSSRefD[eId][str(partNo)])
                if ("auto_match_status" in eSSRefD[eId][str(partNo)]) and eSSRefD[eId][str(partNo)]["auto_match_status"]:
                    autoList.append(eSSRefD[eId][str(partNo)])
                else:
                    rList.append(eSSRefD[eId][str(partNo)])
                #
            #
            if (eId in ownRefD) and ownRefD[eId] and (str(partNo) in ownRefD[eId]) and ownRefD[eId][str(partNo)]:
                if ("selfref" in ownRefD[eId][str(partNo)]) and ownRefD[eId][str(partNo)]["selfref"]:
                    selfRefFlag = True
                else:
                    # rList.append(ownRefD[eId][str(partNo)])
                    if ("auto_match_status" in ownRefD[eId][str(partNo)]) and ownRefD[eId][str(partNo)]["auto_match_status"]:
                        autoList.append(ownRefD[eId][str(partNo)])
                    else:
                        rList.append(ownRefD[eId][str(partNo)])
                    #
                #
            #
            # put auto_match_status=True references at the bottom
            if len(autoList) > 0:
                rList.extend(autoList)
            #
            if (len(rList) > 0) and (not selfRefFlag):
                if (("statistics" in rList[-1]) and (rList[-1]["statistics"][2] > 0.899999)) or (len(autoList) > 0):
                    self.__seqLabel.set(seqType="ref", seqInstId=sId, seqPartId=partNo, seqAltId=len(rList), seqVersion=1)
                    self.__defSelList.append(self.__seqLabel.pack())
                    ref_label[partNo] = [self.__seqLabel.pack()]
                    self._entityAlignInfoMap[sId]["alignids"].append(self.__seqLabel.pack())
                else:
                    self.__defSelList.append("selfref_" + str(sId) + "_" + str(partNo))
                    self._entityAlignInfoMap[sId]["alignids"].append("selfref_" + str(sId) + "_" + str(partNo))
                #
            else:
                self.__defSelList.append("selfref_" + str(sId) + "_" + str(partNo))
                self._entityAlignInfoMap[sId]["alignids"].append("selfref_" + str(sId) + "_" + str(partNo))
                if not rList:
                    continue
                #
            #
            for (altId, rD) in enumerate(rList, start=1):
                if "alignment" in rD:
                    self.__seqLabel.set(seqType="ref", seqInstId=sId, seqPartId=partNo, seqAltId=altId, seqVersion=1)
                    ref_align_index[self.__seqLabel.pack()] = rD["alignment"]
                #
                self.setSequence(rD["seq_tup_list"], "ref", sId, seqPartId=partNo, seqAltId=altId)
                self.setFeature(
                    self.getRefFeatureObj(eD["POLYMER_TYPE_CODE"], partNo, int(pD["SEQ_NUM_BEG"]), int(pD["SEQ_NUM_END"]), pD["SEQ_PART_TYPE"], rD).get(),
                    "ref",
                    sId,
                    seqPartId=partNo,
                    seqAltId=altId,
                )
            #
        #
        self._entityAlignInfoMap[sId]["ref_label"] = ref_label
        self._entityAlignInfoMap[sId]["auth_ref_align_index"] = ref_align_index

    def setDefaultSelectedIds(self):
        """ """
        self.setSelectedIds(self.__defSelList)

    def getAuthFeatureObj(self, pdbId, eD, partNo, pD):
        """ """
        self.__seqFeature.clear()
        self.__seqFeature.setPolymerType(eD["POLYMER_TYPE_CODE"])
        self.__seqFeature.setPolymerLinkingType(eD["POLYMER_LINKING_TYPE"])
        self.__seqFeature.setId(dbName="PDB", dbCode=pdbId, dbAccession=pdbId)
        #
        self.__seqFeature.setSource(
            organism=pD["SOURCE_NAME"],
            strain=pD["SOURCE_STRAIN"],
            taxid=pD["SOURCE_TAXID"],
            gene=pD["SOURCE_GENE_NAME"],
            method=eD["SOURCE_METHOD"],
            commonName=pD["SOURCE_COMMON_NAME"],
            variant=pD["SOURCE_VARIANT"],
        )
        self.__seqFeature.setSourceOrig(
            organism=pD["SOURCE_NAME"],
            strain=pD["SOURCE_STRAIN"],
            taxid=pD["SOURCE_TAXID"],
            gene=pD["SOURCE_GENE_NAME"],
            method=eD["SOURCE_METHOD"],
            commonName=pD["SOURCE_COMMON_NAME"],
            variant=pD["SOURCE_VARIANT"],
        )
        #
        seqNumBeg = int(pD["SEQ_NUM_BEG"])
        seqNumEnd = int(pD["SEQ_NUM_END"])
        seqPartType = pD["SEQ_PART_TYPE"]
        self.__seqFeature.setAuthPartDetails(partNo, seqNumBeg, seqNumEnd, seqPartType)
        self.__seqFeature.setAuthPartDetailsOrig(partNo, seqNumBeg, seqNumEnd, seqPartType)
        #
        self.__seqFeature.setEntityDescription(description=eD["ENTITY_DESCRIPTION"])
        self.__seqFeature.setEntityDescriptionOrig(description=eD["ENTITY_DESCRIPTION"])
        #
        self.__seqFeature.setEntitySynonyms(synonyms=eD["ENTITY_SYNONYMS"])
        self.__seqFeature.setEntitySynonymsOrig(synonyms=eD["ENTITY_SYNONYMS"])
        #
        self.__seqFeature.setEntityEnzymeClass(ec=eD["ENZYME_CLASS"])
        self.__seqFeature.setEntityEnzymeClassOrig(ec=eD["ENZYME_CLASS"])
        #
        self.__seqFeature.setEntityFragmentDetails(details=eD["FRAGMENT_DETAILS"])
        self.__seqFeature.setEntityFragmentDetailsOrig(details=eD["FRAGMENT_DETAILS"])
        #
        self.__seqFeature.setEntityMutationDetails(details=eD["MUTATION_DETAILS"])
        self.__seqFeature.setEntityMutationDetailsOrig(details=eD["MUTATION_DETAILS"])
        #
        self.__seqFeature.setEntityDetails(details=eD["ENTITY_DETAILS"])
        self.__seqFeature.setEntityDetailsOrig(details=eD["ENTITY_DETAILS"])
        #
        self.__seqFeature.setHostOrgDetails(
            source=pD["HOST_ORG_SOURCE"],
            strain=pD["HOST_ORG_STRAIN"],
            taxid=pD["HOST_ORG_TAXID"],
            vector=pD["HOST_ORG_VECTOR"],
            vectorType=pD["HOST_ORG_VECTOR_TYPE"],
            plasmid=pD["HOST_ORG_PLASMID"],
            commonName=pD["HOST_ORG_COMMON_NAME"],
            cellLine=pD["HOST_ORG_CELL_LINE"],
            variant=pD["HOST_ORG_VARIANT"],
        )
        self.__seqFeature.setHostOrgDetailsOrig(
            source=pD["HOST_ORG_SOURCE"],
            strain=pD["HOST_ORG_STRAIN"],
            taxid=pD["HOST_ORG_TAXID"],
            vector=pD["HOST_ORG_VECTOR"],
            vectorType=pD["HOST_ORG_VECTOR_TYPE"],
            plasmid=pD["HOST_ORG_PLASMID"],
            commonName=pD["HOST_ORG_COMMON_NAME"],
            cellLine=pD["HOST_ORG_CELL_LINE"],
            variant=pD["HOST_ORG_VARIANT"],
        )
        #
        return self.__seqFeature

    def getRefFeatureObj(self, polymerType, partNo, seqNumBeg, seqNumEnd, seqPartType, rD):
        """ """
        self.__seqFeature.clear()
        # disambiguate the organism and strain data -
        if ("strain" in rD) and (len(rD["strain"]) > 0):
            self.__seqFeature.setSource(organism=rD["source_scientific"], strain=rD["strain"], taxid=rD["taxonomy_id"], commonName=rD["source_common"])
        elif rD["db_name"] in ["SP", "TR", "UNP"]:
            org, strain = self.__seqFeature.decodeUniProtSourceOrganism(rD["source_scientific"])
            self.__seqFeature.setSource(organism=org, strain=strain, taxid=rD["taxonomy_id"], commonName=rD["source_common"])
        else:
            self.__seqFeature.setSource(organism=rD["source_scientific"], taxid=rD["taxonomy_id"], commonName=rD["source_common"])
        #
        self.__seqFeature.setId(dbName=rD["db_name"], dbCode=rD["db_code"], dbAccession=rD["db_accession"], dbIsoform=rD["db_isoform"])
        self.__seqFeature.setDbIsoformDescription(description=rD["db_isoform_description"])

        self.__seqFeature.setRefSeqNames(proteinName=rD["name"], synonyms=rD["synonyms"], geneName=rD["gene"])
        self.__seqFeature.setRefSeqDetails(enzymeClass=rD["ec"], description=rD["db_description"], comments=rD["comments"], keywords=rD["keyword"])
        if ("variant" in rD) and rD["variant"]:
            self.__seqFeature.setRefSeqVariant(variant=rD["variant"])
        #
        if ("entity_fragment_details" in rD) and rD["entity_fragment_details"]:
            self.__seqFeature.setRefSeqFragmentDetails(fragmentDetails=rD["entity_fragment_details"])
        #
        if "hitFrom" in rD:
            self.__seqFeature.setItem("REF_MATCH_BEGIN", int(rD["hitFrom"]))
        #
        if "hitTo" in rD:
            self.__seqFeature.setItem("REF_MATCH_END", int(rD["hitTo"]))
        #
        if "id" in rD:
            self.__seqFeature.setItem("ORG_ORDER_ID", int(rD["id"]))
        #
        if "db_length" in rD:
            self.__seqFeature.setItem("FULL_LENGTH", int(rD["db_length"]))
        #
        self.__seqFeature.setPolymerType(polymerType)
        if ("sort_order" in rD) and ("sort_metric" in rD):
            self.__seqFeature.setRefSortOrder(sortIndex=rD["sort_order"], sortMetric=rD["sort_metric"])
        #
        if "seq_sim" in rD:
            self.__seqFeature.setItem("AUTH_REF_SEQ_SIM_BLAST", rD["seq_sim"])
        #
        self.__seqFeature.setAuthPartDetails(partNo, seqNumBeg, seqNumEnd, seqPartType)
        if "statistics" in rD:
            self.__seqFeature.setAuthRefAlignDetails(seqLen=len(rD["seq_tup_list"]), alignLen=rD["statistics"][0], seqSim=rD["statistics"][1], seqSimWithGaps=rD["statistics"][2])
        #
        # For 100% sequence match previous processed entry
        #
        for item in ("REF_ENTRY_ID", "REF_ENTRY_ENTITY_ID", "REF_ENTRY_STATUS", "REF_ENTRY_ANN"):
            if item in rD:
                self.__seqFeature.setItem(item, rD[item])
            #
        #
        if "author_provided_id" in rD:
            self.__seqFeature.setItem("IS_AUTH_PROVIDED_ID", rD["author_provided_id"])
        #
        return self.__seqFeature

    def mergePartialAnnotatedResult(self, newOldEntityIdMap, dataStoreFile):
        """Merge previous partial SequenceDataStore object values into current SequenceDataStore object"""
        if (not newOldEntityIdMap) or (not dataStoreFile) or (not os.access(dataStoreFile, os.R_OK)):
            return
        #
        prevSDS = SequenceDataStore(reqObj=self._reqObj, fileName=dataStoreFile, verbose=self._verbose, log=self._lfh)
        #
        oldNewSeqIdMap = {}
        for newId, oldId in newOldEntityIdMap.items():
            oldNewSeqIdMap["auth_" + oldId] = newId
            oldNewSeqIdMap["ref_" + oldId] = newId
            oldNewSeqIdMap["selfref_" + oldId] = newId
            self.__copyData(prevSDS, "auth", newId, oldId)
            self.__copyData(prevSDS, "ref", newId, oldId)
            chainIdList = prevSDS.getGroup(oldId)
            for chainId in chainIdList:
                oldNewSeqIdMap["xyz_" + chainId] = chainId
                self.__copyData(prevSDS, "xyz", chainId, chainId)
            #
        #
        # Update default selected Ids
        #
        prevSelectedIdList = prevSDS.getSelectedIds()
        for seqId in prevSelectedIdList:
            tL = str(seqId).split("_")
            key = tL[0] + "_" + tL[1]
            if key not in oldNewSeqIdMap:
                continue
            #
            tL[1] = oldNewSeqIdMap[key]
            self.__defSelList.append("_".join(tL))
        #

    def __getPartList(self, entityId):
        """Return a list of dictionaries of entity part details.
        Return dictionary has the following keys -

        pD["SEQ_NUM_BEG"]
        pD["SEQ_NUM_END"]
        pD["SOURCE_TAXID"]
        pD["SEQ_PART_TYPE"]
        pD["SEQ_PART_ID"]
        """
        pL = []
        #
        partIdList = self.getPartIdList("auth", entityId)
        for partId in partIdList:
            verList = self.getVersionIdList("auth", entityId, partId, 1)
            if len(verList) == 0:
                continue
            #
            self.__seqFeature.set(self.getFeature("auth", entityId, partId, 1, verList[0]))
            pId, pSeqBegin, pSeqEnd, pType = self.__seqFeature.getAuthPartDetails()
            #
            pD = {}
            pD["SEQ_NUM_BEG"] = pSeqBegin
            pD["SEQ_NUM_END"] = pSeqEnd
            pD["SEQ_PART_TYPE"] = pType
            pD["SEQ_PART_ID"] = pId
            pD["SOURCE_TAXID"] = self.__seqFeature.getSourceTaxId()
            pL.append(pD)
        #
        return pL

    def __copyData(self, prevSDS, seqType, newInstId, oldInstId):
        """Copy sequence & feature data from previous SequenceDataStore object into current SequenceDataStore object"""
        partIdList = prevSDS.getPartIds(oldInstId, dataType="sequence", seqType=seqType)
        for partId in partIdList:
            altIdList = prevSDS.getAlternativeIds(oldInstId, dataType="sequence", seqType=seqType, partId=partId)
            for altId in altIdList:
                versionList = prevSDS.getVersionIds(oldInstId, partId=partId, altId=altId, dataType="sequence", seqType=seqType)
                for version in versionList:
                    sequence = prevSDS.getSequence(oldInstId, seqType, partId=partId, altId=altId, version=version)
                    if sequence:
                        self.setSequence(copy.deepcopy(sequence), seqType, newInstId, seqPartId=partId, seqAltId=altId, seqVersion=version)
                    #
                    feature = prevSDS.getFeature(oldInstId, seqType, partId=partId, altId=altId, version=version)
                    if feature:
                        self.setFeature(copy.deepcopy(feature), seqType, newInstId, seqPartId=partId, seqAltId=altId, seqVersion=version)
                    #
                #
            #
        #
