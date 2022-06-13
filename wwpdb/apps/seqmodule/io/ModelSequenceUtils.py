##
# File:  ModelSequenceUtils.py
# Date:  21-Feb-2013
#
# Updates:
#   23-Feb-2013 jdw  generalized source feature information to support multi-part entities
#   27-Feb-2013 jdw  use instance nomenclature uniformly.  Add sequence separated sequence
#                    identifier to entityD structure
##
"""
Assemble sequence details from the model coordinate data.
"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.09"

import sys

from wwpdb.apps.seqmodule.io.PdbxIoUtils import ModelFileIo
from wwpdb.apps.seqmodule.util.SequenceReferenceData import SequenceReferenceData


class ModelSequenceUtils(object):
    """
    Assemble sequence details from the model coordinate data.
    """

    def __init__(self, dataContainer=None, verbose=False, log=sys.stderr):
        self.__verbose = verbose
        # self.__debug = False
        self.__lfh = log
        self.__sdf = ModelFileIo(dataContainer=dataContainer, verbose=self.__verbose, log=self.__lfh)
        self.__cName = dataContainer.getName()
        #
        self.__entityD = {}
        #

    def getEntitySequenceDetails(self):
        """
                    Extract the author entity sequence data
                    and source information from the current CIF file.

                    Return dictionary:
                       Dictionary with keys ->
                                IDENTIFIER             = the sequence identifier presented for this entity (e.g. first auth_asym_id)
                                ENTITY_ID              = copy of the entity_id
                                INSTANCE_LIST = []
                                POLYMER_LINKING_TYPE   = polypeptide(L)("_entity_poly.type")
                                SEQ_ENTITY_1           = PNFSGNWKI...  ("_entity_poly.pdbx_seq_one_letter_code")
                                SOURCE_METHOD          = MAN           ("_entity.src_method")
                                ENTITY_DESCRIPTION      = _entity.pdbx_description
                                ENTITY_SYNONYMS         = from entity_name_common
                                POLYMER_REFSEQ_DB_NAME = db_name, author provided from entity_poly
                                POLYMER_REFSEQ_DB_ID   = db_id, author provided from entity_poly

                                ENZYME_CLASS           = EC code
                                FRAGMENT_DETAILS       = fragment information
                                ENTITY_DETAILS         = entity detail information
                                MUTATION_DETAILS       = mutation information
        mod

                                PART_LIST              = list of dictionaries with keys:
                                         'SOURCE_NAME'   scientific name
                                         'SOURCE_TAXID'  NCBI Taxonomy ID
                                         'SOURCE_STRAIN' strain name
                                         'SEQ_NUM_BEG'   begining sequence offset (1-sequence length )
                                         'SEQ_NUM_END'   ending sequence offset
                                         'SEQ_PART_TYPE' 'N-terminal tag|C-terminal tag'|'Biological sequence'|'Linker'
                                         'SEQ_PART_ID'   1,2,3...

                                         'SOURCE_GENE_NAME'   gene name
                                         'HOST_ORG_SOURCE'    host Org Source organism scientific name
                                         'HOST_ORG_VECTOR'    host Org Vector type
                                         'HOST_ORG_STRAIN'    host Org Strain name
                                         'HOST_ORG_TAXID'     host Org NCBI Tax Id
                                         'HOST_ORG_PLASMID'   host Org Plasmid name

        """
        if self.__verbose:
            self.__lfh.write("+ModelSequenceUtils.getEntitySequenceDetails() from file %s\n" % self.__cName)
        self.__entityD = self.__getEntityDetails()
        return self.__entityD

    def getSequenceAssignmentDetails(self):
        """Return dictionaries of reference assignment details provided by depositor and assigned by the archive."""
        depSeqAssign = {}
        polyEntityList = self.__sdf.getEntityPolyList()
        for eId in polyEntityList:
            refList, difD = self.__sdf.getDepositorSeqDbRefDetails(entityId=eId)
            depSeqAssign[eId] = {}
            depSeqAssign[eId]["ref_list"] = refList
            depSeqAssign[eId]["ref_dif_dict"] = difD

        seqAssign = {}
        for eId in polyEntityList:
            refList, alignD, difD = self.__sdf.getArchiveSeqDbRefDetails(entityId=eId)
            seqAssign[eId] = {}
            seqAssign[eId]["ref_list"] = refList
            seqAssign[eId]["ref_align_dict"] = alignD
            seqAssign[eId]["ref_dif_dict"] = difD

        return depSeqAssign, seqAssign

    def getCoordinateSequenceDetails(self, useAtomSite=True, linkInstD=None):
        """
         Extract coordinate sequence data from data model.

         Return dictionary:

        instD{} =[(res3, orig auth index + ins code (str), comment/details,  align index placeholder,linkOutlierFlag),(),..]

        as a dictionary with chain id key containing sequences stored as a list of tuples.

        Sequence tuple - (3-letter-code, original auth residue index w/ins.code (str), comment/details, alignment index, linkOutlierFlag(bool))

        """
        if self.__verbose:
            self.__lfh.write("+ModelSequenceUtils.getCoordinateSequenceDetails() use atom_site %r from file %s\n" % (useAtomSite, self.__cName))

        instD = {}
        if useAtomSite:
            #
            nAtoms = self.__sdf.getAtomSiteTableLength()
            if nAtoms > 1:
                if self.__verbose:
                    self.__lfh.write("+ModelSequenceUtils.__getCoordinateSequenceDetails() atom site table length = %d\n" % nAtoms)
                    instD = self.__sdf.getSequenceFeaturesFromAtomSite(linkInstD=linkInstD)
                    if self.__verbose:
                        for k, v in instD.items():
                            self.__lfh.write("+ModelSequenceUtils.__getCoordinateSequenceDetails() chain %s sequence length %d\n" % (k, len(v)))
                else:
                    # no coordinates found
                    if self.__verbose:
                        self.__lfh.write("+ModelSequenceUtils.__getCoordinateSequenceDetails() no coordinate data to process.\n")
        else:
            # Check for residue sequence data in category - pdbx_poly_seq_scheme
            if self.__sdf.getResidueTableLength() > 2:

                instD = self.__getCoordinateSequencesFromMappingTable()
                if self.__verbose:
                    self.__lfh.write("+ModelSequenceUtils.__getCoordinateDetails() sequences from residue table for %d chains\n" % len(instD))
            else:
                if self.__verbose:
                    self.__lfh.write("+ModelSequenceUtils.__getCoordinateDetails() no coordinate mapping data to process.\n")
        #
        return instD

    def __getCoordinateSequencesFromMappingTable(self):
        """Get coordinate sequence from the pdbx_poly_seq_scheme residue mapping table.

        Return dictionary:  {chain/auth_asym_id} -> [(comp_id, auth_seq_id, '', component index in chain), (), ...]
        """
        if len(self.__entityD) == 0:
            self.__entityD = self.__getEntityDetails()

        instD = {}
        for _id, eD in self.__entityD.items():
            instIdList = eD["INSTANCE_LIST"]
            for instId in instIdList:
                sL = self.__sdf.getCoordinateSequenceList(instId)
                instD[instId] = sL
        return instD

    def __getEntityDetails(self):
        """Internal working method to assemble dictionary of entity features.

        Dictionary with keys ->
                 IDENTIFIER             = the sequence identifier presented for this entity (e.g. first auth_asym_id)
                 ENTITY_ID              = copy of the entity_id
                 INSTANCE_LIST = []
                 POLYMER_LINKING_TYPE   = polypeptide(L)("_entity_poly.type")
                 SEQ_ENTITY_1           = PNFSGNWKI...  ("_entity_poly.pdbx_seq_one_letter_code")
                 SOURCE_METHOD          = MAN           ("_entity.src_method")
                 ENTITY_DESCRIPTION      = _entity.pdbx_description
                 ENTITY_SYNONYMS         = from entity_name_common
                 POLYMER_REFSEQ_DB_NAME = db_name, author provided from entity_poly
                 POLYMER_REFSEQ_DB_ID   = db_id, author provided from entity_poly

                 ENZYME_CLASS           = EC code
                 FRAGMENT_DETAILS       = fragment information
                 ENTITY_DETAILS         = entity detail information
                 MUTATION_DETAILS       = mutation information


                 PART_LIST              = list of dictionaries with keys:
                          'SOURCE_NAME'   scientific name
                          'SOURCE_TAXID'  NCBI Taxonomy ID
                          'SOURCE_STRAIN' strain name
                          'SEQ_NUM_BEG'   begining sequence offset (1-sequence length )
                          'SEQ_NUM_END'   ending sequence offset
                          'SEQ_PART_TYPE' 'N-terminal tag|C-terminal tag'|'Biological sequence'|'Linker'
                          'SEQ_PART_ID'   1,2,3...

                          'SOURCE_GENE_NAME'   gene name
                          'HOST_ORG_SOURCE'    host Org Source organism scientific name
                          'HOST_ORG_VECTOR'    host Org Vector type
                          'HOST_ORG_STRAIN'    host Org Strain name
                          'HOST_ORG_TAXID'     host Org NCBI Tax Id
                          'HOST_ORG_PLASMID'   host Org Plasmid name
                          'HOST_ORG_CELL_LINE' host Org cell line

        """
        sdr = SequenceReferenceData(self.__verbose, self.__lfh)
        entityD = {}
        # pdbId = self.__sdf.getDbCode("PDB")
        entryId = self.__sdf.getEntryId()
        cName = self.__sdf.getContainerName()

        polyEntityList = self.__sdf.getEntityPolyList()
        for entityId in polyEntityList:

            if self.__verbose:
                self.__lfh.write("+ModelSequenceUtils.__getEntityDetails() entity id= %s\n" % entityId)
            sD = {}
            sD["ENTITY_ID"] = entityId
            sD["ENTRY_ID"] = entryId if len(entryId) > 0 else cName
            instList = self.__sdf.getPdbChainIdList(entityId)
            sD["INSTANCE_LIST"] = instList
            if self.__verbose:
                self.__lfh.write("+ModelSequenceUtils.__getEntityDetails() entity id= %s chain list %r\n" % (entityId, sD["INSTANCE_LIST"]))
            # if len(instList) > 0:
            #     sId = instList[0]
            # else:
            #     sId = "E_" + str(entityId)

            # sD['IDENTIFIER']=sId
            #
            # changed this to entityId
            sD["IDENTIFIER"] = entityId
            sD["ENTITY_DESCRIPTION"] = self.__sdf.getEntityDescription(entityId)
            sD["ENTITY_SYNONYMS"] = self.__sdf.getEntityName(entityId)
            sD["POLYMER_LINKING_TYPE"] = self.__sdf.getPolymerEntityType(entityId)
            #
            sD["SOURCE_METHOD"] = str(self.__sdf.getSourceMethod(entityId)).upper()

            sD["ENZYME_CLASS"] = str(self.__sdf.getEntityEnzymeClass(entityId))
            sD["FRAGMENT_DETAILS"] = str(self.__sdf.getEntityFragmentDetails(entityId))
            sD["ENTITY_DETAILS"] = str(self.__sdf.getEntityDetails(entityId))
            sD["MUTATION_DETAILS"] = str(self.__sdf.getEntityMutationDetails(entityId))

            if self.__verbose:
                self.__lfh.write(
                    "+ModelSequenceUtils.__getEntityDetails() entity id= %s source method %s linking type %s\n" % (entityId, sD["SOURCE_METHOD"], sD["POLYMER_LINKING_TYPE"])
                )

            #
            sD["SEQ_ENTITY_1"] = self.__sdf.getSequence(entityId)
            r1L, _r3L = sdr.parseSequence(sD["SEQ_ENTITY_1"], sdr.getPolymerTypeCode(sD["POLYMER_LINKING_TYPE"]))
            sD["SEQ_ENTITY_1_CAN"] = "".join(r1L)

            if self.__verbose:
                self.__lfh.write("+ModelSequenceUtils.__getEntityDetails() entity id= %s sequence %s\n" % (entityId, sD["SEQ_ENTITY_1"]))
                self.__lfh.write("+ModelSequenceUtils.__getEntityDetails() entity id= %s cannonical sequence %s\n" % (entityId, sD["SEQ_ENTITY_1_CAN"]))
                self.__lfh.write("+ModelSequenceUtils.__getEntityDetails() entity id= %s cannonical sequence length %d\n" % (entityId, len(sD["SEQ_ENTITY_1_CAN"])))
            seqLength = len(sD["SEQ_ENTITY_1_CAN"])
            #
            # The feature list will include dictionary of source org, taxid, strain and part residue boundaries
            #
            sD["PART_LIST"] = []
            #
            catList = ["entity_src_gen", "entity_src_nat", "pdbx_entity_src_syn"]
            for cat in catList:
                sD["PART_LIST"].extend(self.__sdf.getSourceDetailsList(entityId, sourceCategoryName=cat, seqLength=seqLength))
                if len(sD["PART_LIST"]) > 0:
                    break
            #
            # if there is no source information then assign some defaults for the partitioning.
            if len(sD["PART_LIST"]) < 1:
                sD["PART_LIST"] = self.__sdf.assignSourceDefaultList(entityId, seqLength=seqLength)
            #
            if self.__verbose:
                self.__lfh.write("+ModelSequenceUtils.__getEntityDetails() entity id= %s part list %r\n" % (entityId, sD["PART_LIST"]))
            #
            # Check missing polymer type assignment
            #
            if len(sD["POLYMER_LINKING_TYPE"]) < 2:
                sD["POLYMER_LINKING_TYPE"] = sdr.guessPolymerType(sD["SEQ_ENTITY_1_CAN"])
            #
            if self.__verbose:
                self.__lfh.write("+ModelSequenceUtils.__getEntityDetails() entity id= %s polymer linking type %s\n" % (entityId, sD["POLYMER_LINKING_TYPE"]))

            entityD[entityId] = sD

        return entityD
