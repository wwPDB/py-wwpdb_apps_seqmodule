##
# File:  PdbxIoUtils.py
# Date:  21-Feb-2013
#
# Updates:
#
# 24-Feb-2013 jdw  Add method to output of sequence match results.
# 25-Feb-2013 jdw  Add attributes for 'sort_order' and 'sort_metric'
# 27-Feb-2013 jdw
# 05-Mar-2013 jdw  add getSequenceFeaturesFromAtomSite()
# 08-Mar-2013 jdw  update feature calculation getSequenceFeaturesFromAtomSite()
# 22-Mar-2013 jdw  fix last residue problem
#  1-Apr-2013 jdw  add entity sequence index to link reader.
#  2-Apr-2013 jdw  add model methods getDepositorAssemblyDetails() getDepositorAssemblyDetailsRcsb()
#  6-Apr-2013 jdw  Avoid stripping (nnn) from all sequences in removeArtifact()
#  7-Apr-2013 jdw  add isoform support
# 14-Jun-2013 jdw  include insertion code in PDB residue index.
#  1-Dec-2013 jdw  Add linkage outlier information to the coordinate sequence as an additional part
#                  of the comment element.  Residues at the beginning/end of a long linkage are separately
#                  labeled (long_begin|long_end)
# 19-Jan-2014 jdw  add                 d['HOST_ORG_CELL_LINE']=hostOrgCellLine
#
#  9-Feb-2014 jdw  add HOST_ORG_VECTOR_TYPE and SOURCE_VARIANT and orig analogs
# 18-Feb-2014 jdw  look for label_comp_id if auth_comp_id is missing.
#  5-May-2014 jdw  fix edge condition with short and similarly numbered chains.
# 22-May-2014 jdw  add title methods
# 01-Feb-2015 jdw  add host org variant
# 01-Jul-2015 jdw  add getAssemblyDetails() and  __getAttributeDictList()
# 02-Feb-2017 ep   add getDepositorAssemblyEvidence(), getDepositorStructOperList(), getDepositorAssemblyGen(),
#                  getDepositorAssemblyClassification() and __getDepositorDetails()
#                  for common method support
# 24-Aug-2017 zf   add PolymerInstanceIo
##
"""
Utility methods for accessing entity, chain anad sequence details from
deposited data files.

Content specific classes operate on PDBx container object input.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.07"


import sys
import re
import traceback

from mmcif.io.IoAdapterPy import IoAdapterPy
from wwpdb.apps.seqmodule.util.SequenceReferenceData import SequenceReferenceData
from mmcif.api.PdbxContainers import DataContainer
from mmcif.api.DataCategory import DataCategory


class PdbxFileIo(object):

    """Read PDBx data files and package content as PDBx container object or container object list
    Write PDBx data using source PDBx container object list source content.

    PDBx container object represents
    """

    def __init__(self, ioObj=IoAdapterPy(), verbose=True, log=sys.stderr):
        """Input processing can be performed using either native Python or C++ Io libraries
        by choosing the appropriate input adapter.
        """

        self.__ioObj = ioObj
        self.__verbose = verbose
        self.__lfh = log

    def getContainer(self, fPath, index=0):
        try:
            cList = self.__ioObj.readFile(fPath)
            return cList[index]
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
            return None

    def getContainerList(self, fPath):
        try:
            cList = self.__ioObj.readFile(fPath)
            return cList
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
            return []

    def writeContainerList(self, fPath, containerList=None):
        try:
            return self.__ioObj.writeFile(fPath, containerList)
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
            return False


class PolymerLinkageIo(object):

    """
    Assemble data from polymer linkage distance data files -

    """

    def __init__(self, dataContainer=None, verbose=True, log=sys.stderr):  # pylint: disable=unused-argument
        self.__currentContainer = dataContainer
        # self.__verbose = verbose
        # self.__lfh = log

    def getPolymerLinkDistances(self):
        """Returns a dictionary of linkage distance information of atypical linkages.  d > 1.7  or d < 1.2

        Example -

        data_RCSB056215
        #
        loop_
        _pdbx_polymer_linkage_distance.id
        _pdbx_polymer_linkage_distance.PDB_model_num
        _pdbx_polymer_linkage_distance.auth_asym_id_1
        _pdbx_polymer_linkage_distance.auth_comp_id_1
        _pdbx_polymer_linkage_distance.auth_seq_id_1
        _pdbx_polymer_linkage_distance.label_seq_id_1
        _pdbx_polymer_linkage_distance.PDB_ins_code_1
        _pdbx_polymer_linkage_distance.auth_asym_id_2
        _pdbx_polymer_linkage_distance.auth_comp_id_2
        _pdbx_polymer_linkage_distance.auth_seq_id_2
        _pdbx_polymer_linkage_distance.label_seq_id_2
        _pdbx_polymer_linkage_distance.PDB_ins_code_2
        _pdbx_polymer_linkage_distance.dist
        1    1 B VAL 7  7   ? B LYS 8  8   ? 1.338
        2    1 B LYS 8  8   ? B GLU 9  9   ? 1.329
        3    1 B GLU 9  9   ? B LEU 10 10  ? 1.334
        4    1 B LEU 10 10  ? B LEU 11 11  ? 1.326
        5    1 B LEU 11 11  ? B GLU 12 12  ? 1.329
        6    1 B GLU 12 12  ? B ALA 13 13  ? 1.328
        7    1 B ALA 13 13  ? B GLY 14 14  ? 1.327
        8    1 B GLY 14 14  ? B VAL 15 15  ? 1.329
        9    1 B VAL 15 15  ? B HIS 16 16  ? 1.342
        10   1 B HIS 16 16  ? B PHE 17 17  ? 1.329

        """
        #
        if not self.__currentContainer.exists("pdbx_polymer_linkage_distance"):
            return []

        refTable = self.__currentContainer.getObj("pdbx_polymer_linkage_distance")
        nRows = refTable.getRowCount()
        #
        colNames = list(refTable.getAttributeList())

        myList = [
            "id",
            "PDB_model_num",
            "auth_asym_id_1",
            "auth_comp_id_1",
            "auth_seq_id_1",
            "label_seq_id_1",
            "PDB_ins_code_1",
            "auth_asym_id_2",
            "auth_comp_id_2",
            "auth_seq_id_2",
            "label_seq_id_2",
            "PDB_ins_code_2",
            "dist",
        ]

        rList = []
        for iRow in range(0, nRows):
            rD = {}
            row = refTable.getRow(iRow)
            for col in myList:
                if col in colNames:
                    val = str(row[colNames.index(col)])
                    if val is None:
                        val = ""
                    elif (val == ".") or (val == "?"):
                        val = ""
                    rD[col] = val
                else:
                    rD[col] = ""
            rList.append(rD)
        return rList


class PolymerInstanceIo(object):

    """
    Assemble data from polymer linkage distance data files -

    """

    def __init__(self, dataContainer=None, verbose=True, log=sys.stderr):  # pylint: disable=unused-argument
        self.__currentContainer = dataContainer
        # self.__verbose = verbose
        # self.__lfh = log

    def getPolymerInstances(self):
        """Returns a list of polymer instances

        Example -

        data_RCSB001900
        #
        loop_
        _pdbx_polymer_instance.auth_asym_id
        _pdbx_polymer_instance.auth_comp_id
        _pdbx_polymer_instance.auth_seq_id
        _pdbx_polymer_instance.comment
        _pdbx_polymer_instance.index
        A LYS 5   ?                           1
        A ASN 6   ?                           2
        A ILE 7   ?                           3
        A VAL 8   ?                           4
        A PHE 9   ?                           5
        A ILE 10  ?                           6
        A GLY 11  ?                           7
        A PHE 12  ?                           8
        A MSE 13  "mean_occ=0.667 disordered" 9
        A GLY 14  ?                           10
        A SER 15  ?                           11
        A GLY 16  ?                           12
        """
        #
        if not self.__currentContainer.exists("pdbx_polymer_instance"):
            return {}

        refTable = self.__currentContainer.getObj("pdbx_polymer_instance")
        nRows = refTable.getRowCount()
        #
        colNames = list(refTable.getAttributeList())

        myList = ["auth_asym_id", "auth_comp_id", "auth_seq_id", "comment", "index"]
        for col in myList:
            if col not in colNames:
                return {}
            #
        #

        rDic = {}
        for iRow in range(0, nRows):
            rD = {}
            row = refTable.getRow(iRow)
            for col in myList:
                val = str(row[colNames.index(col)])
                if val is None:
                    val = ""
                elif (val == ".") or (val == "?"):
                    val = ""
                if col == "index":
                    rD[col] = int(val)
                else:
                    rD[col] = val
                #
            #
            if rD["auth_asym_id"] in rDic:
                rDic[rD["auth_asym_id"]].append((rD["auth_comp_id"], rD["auth_seq_id"], rD["comment"], rD["index"]))
            else:
                rDic[rD["auth_asym_id"]] = [(rD["auth_comp_id"], rD["auth_seq_id"], rD["comment"], rD["index"])]
            #
        #
        return rDic


class ReferenceSequenceIo(object):

    """
    Assemble data from reference sequence matching data files

    """

    def __init__(self, dataContainer=None, verbose=True, log=sys.stderr):
        #
        self.__currentContainer = dataContainer
        self.__verbose = verbose
        self.__lfh = log
        self.__matchEntityAttribNameList = [
            "id",
            "fragment_id",
            "beg_seq_num",
            "end_seq_num",
            "sort_order",
            "sort_metric",
            "db_name",
            "db_code",
            "db_accession",
            "db_isoform",
            "match_length",
            "queryFrom",
            "queryTo",
            "hitFrom",
            "hitTo",
            "identity",
            "positive",
            "gaps",
            "alignLen",
            "query",
            "subject",
            "midline",
            "query_length",
            "name",
            "source_scientific",
            "source_common",
            "taxonomy_id",
            "gene",
            "synonyms",
            "comments",
            "keyword",
            "ec",
            "db_length",
            "db_description",
            "db_isoform_description",
        ]

    def readMatchResults(self):
        """Returns a dictionary of reference sequence details."""
        if not self.__currentContainer.exists("match_entity"):
            return []

        rD = {}
        refTable = self.__currentContainer.getObj("match_entity")
        nRows = refTable.getRowCount()
        #
        colNames = list(refTable.getAttributeList())

        myIntList = ["match_length", "queryFrom", "queryTo", "hitFrom", "hitTo", "identity", "positive", "gaps", "alignLen", "taxonomy_id", "query_length", "db_length"]

        rList = []
        for iRow in range(0, nRows):
            ok = True
            rD = {}
            for ky in self.__matchEntityAttribNameList:
                rD[ky] = ""
            row = refTable.getRow(iRow)
            for col in self.__matchEntityAttribNameList:
                if col in colNames:
                    if col in myIntList:
                        try:
                            rD[col] = int(str(row[colNames.index(col)]))
                        except:  # noqa: E722 pylint: disable=bare-except
                            if col in ["taxonomy_id"]:
                                rD[col] = 0
                                self.__lfh.write(
                                    "+ReferenceSequenceIo.readMatchResults() integer read error row %d at %d colname %s val %s\n"
                                    % (iRow, nRows, col, str(row[colNames.index(col)]))
                                )
                            else:
                                self.__lfh.write(
                                    "+ReferenceSequenceIo.readMatchResults() integer read error skipping row %d at %d colname %s val %s\n"
                                    % (iRow, nRows, col, str(row[colNames.index(col)]))
                                )

                                ok = False
                    elif col in ["query", "subject", "midline"]:
                        rD[col] = str(row[colNames.index(col)]).upper()
                    else:
                        if str(row[colNames.index(col)]) in [".", "?"]:
                            rD[col] = ""
                        else:
                            rD[col] = str(row[colNames.index(col)])
            if ok:
                try:
                    seq_sim = float(rD["identity"]) / float(rD["alignLen"])
                    # self.__lfh.write("+ReferenceSequenceIo.readMatchResults() row %d of %d identity %r alignlen %r seq_sim %f\n"
                    # % (iRow,nRows,rD['identity'],rD['alignLen'],seq_sim))
                except:  # noqa: E722 pylint: disable=bare-except
                    self.__lfh.write("+ReferenceSequenceIo.readMatchResults() failed on row %d of %d identity %r alignlen %r\n" % (iRow, nRows, rD["identity"], rD["alignLen"]))
                    seq_sim = 0.0
                rD["seq_sim"] = seq_sim
                # placeholder
                rD["source_string_distance"] = 100
                rD["taxid_string_distance"] = 100
                rList.append(rD)
        self.__lfh.write("+ReferenceSequenceIo.readMatchResults() nrows  %d return list len %d\n" % (nRows, len(rList)))

        return rList

    def writeMatchResults(self, entityD, outFilePath, matchResults=None):
        """Export the reference sequence matching results for an entity including sub-parts."""
        if matchResults is None:
            return

        entityId = entityD["ENTITY_ID"]
        entryId = entityD["ENTRY_ID"]
        seqFull = entityD["SEQ_ENTITY_1"]

        curContainer = DataContainer("match_entity")
        #
        #  Table = info
        t = DataCategory("info")
        t.appendAttribute("struct_id")
        t.appendAttribute("entity_id")
        t.appendAttribute("sequence")
        t.appendAttribute("fragment_count")
        #
        t.setValue(entryId, "struct_id", 0)
        t.setValue(entityId, "entity_id", 0)
        t.setValue(self.__formatSequence(seqFull), "sequence", 0)

        t.setValue(len(matchResults), "fragment_count", 0)
        curContainer.append(t)

        #
        # Table = match_entity
        #
        t = DataCategory("match_entity")
        for item in self.__matchEntityAttribNameList:
            t.appendAttribute(item)

        #
        # Table = org_sequence
        #
        t1 = DataCategory("org_sequence")
        t1.appendAttribute("id")
        t1.appendAttribute("sequence")

        iRow = 0
        for hitList in matchResults:
            self.__lfh.write("+writeMatchResults() Length of hit list %d\n" % len(hitList))
            for hit in hitList:
                # self.__lfh.write("+writeMatchResults() iRow %d hit dictionary %r\n" % (iRow,hit.items()))
                t.setValue(str(iRow + 1), "id", iRow)
                for attrib in self.__matchEntityAttribNameList:
                    if attrib in hit:
                        t.setValue(str(hit[attrib]), attrib, iRow)
                t1.setValue(str(iRow + 1), "id", iRow)
                if "sequence" in hit:
                    seq = re.sub("[\t \n]", "", hit["sequence"])
                    t1.setValue(self.__formatSequence(str(seq)), "sequence", iRow)

                iRow = iRow + 1
            #
        #

        curContainer.append(t)
        curContainer.append(t1)

        myContainerList = []
        myContainerList.append(curContainer)
        pf = PdbxFileIo(verbose=self.__verbose, log=self.__lfh)
        pf.writeContainerList(outFilePath, myContainerList)

    def __formatSequence(self, sequence):
        num_per_line = 60
        l = len(sequence) / num_per_line  # noqa: E741
        x = len(sequence) % num_per_line
        m = l
        if x:
            m = l + 1

        seq = ""
        for i in range(m):
            n = num_per_line
            if i == l:
                n = x
            seq += sequence[i * num_per_line : i * num_per_line + n]
            if i != (m - 1):
                seq += "\n"

        return seq


class ModelFileIo(object):

    """
    Assemble sample and coordinate sequence details from model coordinate data file.

    """

    def __init__(self, dataContainer=None, verbose=True, log=sys.stderr):
        self.__currentContainer = dataContainer
        self.__verbose = verbose
        self.__lfh = log
        #
        self.__polymerEntityChainDict = {}
        self.__chainPolymerEntityDict = {}
        self.__buildPolymerEntityChainDict()

    def getContainerName(self):
        return self.__currentContainer.getName()

    def getEntryId(self):
        """Return _entry.id"""
        try:
            catObj = self.__currentContainer.getObj("entry")
            return catObj.getValue("id", 0)
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getStructTitle(self):
        """Return _struct.title"""
        try:
            catObj = self.__currentContainer.getObj("struct")
            return catObj.getValue("title", 0)
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getCitationTitle(self, id="primary"):  # pylint: disable=redefined-builtin
        """Return _citation.title  where id = primary"""
        try:
            catObj = self.__currentContainer.getObj("citation")
            vL1 = catObj.selectValuesWhere("title", id, "id")
            v1 = self.__firstOrDefault(vL1, default="")
            return v1
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getPolyPeptideLEntityCount(self):
        """Get the number of entities of type polypeptide(L)"""
        return self.getPolymerEntityCount(type="polypeptide(L)")

    def getPolyPeptideDEntityCount(self):
        """Get the number of entities of type polypeptide(D)"""
        return self.getPolymerEntityCount(type="polypeptide(D)")

    def getEntityCount(self, type):  # pylint: disable=redefined-builtin
        """Returns the integer count of entities of the input 'type'
        (polymer, non-polymer, macrolide, water)
        """
        if type in SequenceReferenceData._entityTypes:  # pylint: disable=protected-access
            catObj = self.__currentContainer.getObj("entity")
            indices = catObj.selectIndices(type, "type")
            return len(indices)
        else:
            return 0

    def getPolymerEntityCount(self, type):  # pylint: disable=redefined-builtin
        """Returns the integer count of polymer entities of the input 'type'
        ++ allowed types are reference._polymerEntityTypes -
        """
        if type in SequenceReferenceData._polymerEntityTypes:  # pylint: disable=protected-access
            catObj = self.__currentContainer.getObj("entity_poly")
            indices = catObj.selectIndices(type, "type")
            return len(indices)
        else:
            return 0

    def __isEmptyValue(self, val):
        if (val is None) or (len(val) == 0) or (val in [".", "?"]):
            return True
        else:
            return False

    def __firstOrDefault(self, valList, default=""):
        if len(valList) > 0 and not self.__isEmptyValue(valList[0]):
            return valList[0]
        else:
            return default

    def getPolymerEntityType(self, entityId):
        catObj = self.__currentContainer.getObj("entity_poly")
        vals = catObj.selectValuesWhere("type", entityId, "entity_id")
        return self.__firstOrDefault(vals, default="")

    def getPolymerEntityDbInfo(self, entityId):
        catObj = self.__currentContainer.getObj("entity_poly")
        vL1 = catObj.selectValuesWhere("pdbx_seq_db_name", entityId, "entity_id")
        vL2 = catObj.selectValuesWhere("pdbx_seq_db_id", entityId, "entity_id")
        v1 = self.__firstOrDefault(vL1, default="")
        v2 = self.__firstOrDefault(vL2, default="")
        return v1, v2

    def getPolymerEntityList(self, type=None):  # pylint: disable=redefined-builtin
        """Returns a list of polymer entity id's  of the input 'type'
        type is an entity type (all, polymer, non-polymer,  any)  or
            one of the polymer entity types.
        """
        try:
            if type in ["any", "all", "polymer"]:
                tType = "polymer"
                catObj = self.__currentContainer.getObj("entity")
                eList = catObj.selectValuesWhere("id", tType, "type")
                return eList
            elif type in SequenceReferenceData._polymerEntityTypes:  # pylint: disable=protected-access
                catObj = self.__currentContainer.getObj("entity_poly")
                eList = catObj.selectValuesWhere("entity_id", type, "type")
            else:
                return []
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+PdbxIoUtils.getPolymerEntityList() WARNING - Likely missing entity or entity_poly data categories\n")

        return []

    def getDbCode(self, dbId):
        """Return the database code for the input database id/name"""
        try:
            catObj = self.__currentContainer.getObj("database_2")
            vals = catObj.selectValuesWhere("database_code", dbId, "database_id")
            return self.__firstOrDefault(vals, default="NOID")
        except:  # noqa: E722 pylint: disable=bare-except
            return "NOID"

    def getSequence(self, entityId):
        """Return one-letter-code sequence for the input entity."""
        try:
            catObj = self.__currentContainer.getObj("entity_poly")
            if catObj.hasAttribute("pdbx_seq_one_letter_code"):
                vals = catObj.selectValuesWhere("pdbx_seq_one_letter_code", entityId, "entity_id")
            elif catObj.hasAttribute("ndb_seq_one_letter_code"):
                vals = catObj.selectValuesWhere("ndb_seq_one_letter_code", entityId, "entity_id")
            else:
                vals = []
            oneLetterCodeSeq = self.__firstOrDefault(vals, default="")
            return self.___removeArtifacts(oneLetterCodeSeq)
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
            return ""

    def ___removeArtifacts(self, oneLetterCodeSeq):
        seq = oneLetterCodeSeq.upper()
        seq = re.sub("[\t \n]", "", seq)
        # JDW Try this for now
        seq = re.sub("\\(PYL\\)", "O", seq)
        seq = re.sub("\\(SEC\\)", "U", seq)
        #
        # seq = re.sub('\(MSE\)', 'M', seq)
        # seq = re.sub('\([A-Z]{2,3}\)', 'X', seq)
        return seq

    def getSourceMethod(self, entityId):
        """Return source method for the input entity -"""
        try:
            catObj = self.__currentContainer.getObj("entity")
            vals = catObj.selectValuesWhere("src_method", entityId, "id")
            return self.__firstOrDefault(vals, default="")
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getEntityEnzymeClass(self, entityId):
        """Return any entity E.C. assignment"""
        try:
            catObj = self.__currentContainer.getObj("entity")
            vals = catObj.selectValuesWhere("pdbx_ec", entityId, "id")
            return self.__firstOrDefault(vals, default="")
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getEntityFragmentDetails(self, entityId):
        """Return any entity fragment details ."""
        try:
            catObj = self.__currentContainer.getObj("entity")
            vals = catObj.selectValuesWhere("pdbx_fragment", entityId, "id")
            return self.__firstOrDefault(vals, default="")
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getEntityMutationDetails(self, entityId):
        """Return any entity details ."""
        try:
            catObj = self.__currentContainer.getObj("entity")
            vals = catObj.selectValuesWhere("pdbx_mutation", entityId, "id")
            return self.__firstOrDefault(vals, default="")
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getEntityDetails(self, entityId):
        """Return any entity details ."""
        try:
            catObj = self.__currentContainer.getObj("entity")
            vals = catObj.selectValuesWhere("details", entityId, "id")
            return self.__firstOrDefault(vals, default="")
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getEntityDescription(self, entityId):
        """
        Return a the _entity.pdbx_description or empty string.
        """
        try:
            catObj = self.__currentContainer.getObj("entity")
            if catObj.hasAttribute("pdbx_description"):
                vals = catObj.selectValuesWhere("pdbx_description", entityId, "id")
            elif catObj.hasAttribute("ndb_chain_id"):
                vals = catObj.selectValuesWhere("ndb_description", entityId, "id")
            else:
                vals = []
            return self.__firstOrDefault(vals, "")
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getEntityName(self, entityId):
        """
        Return the _entity.name or empty string.
        """
        try:
            catObj = self.__currentContainer.getObj("entity_name_com")
            if catObj.hasAttribute("name"):
                vals = catObj.selectValuesWhere("name", entityId, "entity_id")
            else:
                vals = []
            return self.__firstOrDefault(vals, "")
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getResidueTableLength(self):
        if self.__currentContainer.exists("pdbx_poly_seq_scheme"):
            tab = self.__currentContainer.getObj("pdbx_poly_seq_scheme")
            return tab.getRowCount()
        else:
            return 0

    def getPdbChainIdList(self, entityId):
        catObj = self.__currentContainer.getObj("entity_poly")
        if catObj.hasAttribute("pdbx_strand_id"):
            vals = catObj.selectValuesWhere("pdbx_strand_id", entityId, "entity_id")
        elif catObj.hasAttribute("ndb_chain_id"):
            vals = catObj.selectValuesWhere("ndb_chain_id", entityId, "entity_id")
        else:
            vals = []
        st = self.__firstOrDefault(vals, "")

        tList = []
        if st is not None and len(st) > 0:
            if (len(st) > 1) and (st.count(",") > 0):
                tList = st.split(",")
            elif (len(st) > 1) and (st.count(" ") > 0):
                tList = st.split()
            else:
                tList = st.split(",")
        rList = []
        for ch in tList:
            if len(ch) == 0 or ch in [".", "?"]:
                continue
            rList.append(str(ch).strip())
        return rList

    def getSequence3AlignList(self, chainId):
        """
        Return a list of lists containing (PDB Monomer id (xyz), PDBx mononer(entity), Ent. Poly Seq Num,  PDB residue Num)
        for the aligned entity and coordinate sequences (expressed as 3-letter-codes).

        """
        aList = ["pdb_mon_id", "mon_id", "seq_id", "pdb_seq_num"]
        aa3List = []
        try:
            tabO = self.__currentContainer.getObj("pdbx_poly_seq_scheme")
            aa3List = tabO.selectValueListWhere(aList, chainId, "pdb_strand_id")
        except:  # noqa: E722 pylint: disable=bare-except
            return aa3List

        return aa3List

    def residueMapTableExists(self):
        if self.__currentContainer.exists("pdbx_poly_seq_scheme"):
            return True
        else:
            return False

    def getCoordinateSequenceList(self, chainId):
        """
        Return a list of tuples containing (PDB Monomer 3-letter-code, PDB residue num, '', Ent. Poly Seq Num, )
        for the aligned entity and coordinate sequences (expressed as 3-letter-codes).

        """
        aa3List = []
        try:
            catObj = self.__currentContainer.getObj("pdbx_poly_seq_scheme")
            indices = catObj.selectIndices(chainId, "pdb_strand_id")
            aa3List = []
            for ii in indices:
                monId = catObj.getValue("pdb_mon_id", ii)
                if (len(monId) > 0) and (monId != "?") and (monId != "."):
                    a = catObj.getValue("pdb_mon_id", ii)
                    b = catObj.getValue("pdb_seq_num", ii)
                    c = catObj.getValue("seq_id", ii)
                    aa3List.append((a, b, "", c))
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)

        return aa3List

    def __getFirstValueFromList(self, attributeNameList, catObj=None, rowIndex=0):
        """Return the value from the first non null attribute found in the input attribute list
        from the input category object/rowIndex.
        """
        try:
            for at in attributeNameList:
                if catObj.hasAttribute(at):
                    val = catObj.getValue(at, rowIndex)
                    if not self.__isEmptyValue(val):
                        return val
            return ""
        except:  # noqa: E722 pylint: disable=bare-except
            return None

    def getSourceDetailsList(self, entityId, sourceCategoryName="entity_src_gen", seqLength=0):
        """Return a list of dictionaries containing source feature details.

               'SOURCE_NAME'   scientific name
        'SOURCE_COMMON_NAME'   common name
               'SOURCE_TAXID'  NCBI Taxonomy ID
               'SOURCE_STRAIN' strain name
               'SOURCE_VARIANT' variant

               'SOURCE_GENE_NAME'   gene name
               'HOST_ORG_SOURCE'    host Org Source organism scientific name
               'HOST_ORG_VECTOR'    host Org Vector
          'HOST_ORG_VECTOR_TYPE'    host Org Vector type
               'HOST_ORG_STRAIN'    host Org Strain name
               'HOST_ORG_TAXID'     host Org NCBI Tax Id
               'HOST_ORG_PLASMID'   host Org Plasmid name
           'HOST_ORG_COMMON_NAME'   host Org common name
           'HOST_ORG_VARIANT'       host Org variant


               'SEQ_PART_ID'   _entity_src_gen.pdbx_src_id  = partId
               'SEQ_PART_TYPE' _entity_src_gen.pdbx_seq_type =  'N-terminal tag|C-terminal tag'|'Biological sequence'|'Linker'

                Regions of the sequence pertaining to this part -

               'SEQ_NUM_BEG'   begining sequence offset (1-sequence length )
               'SEQ_NUM_END'   ending sequence offset

                *_ORIG'        copy of the above




        """
        rL = []
        if not self.__currentContainer.exists(sourceCategoryName):
            if self.__verbose:
                self.__lfh.write("+PdbxIoUtils.getSourceDetailsList() No category %s\n" % sourceCategoryName)
            return rL
        try:
            catObj = self.__currentContainer.getObj(sourceCategoryName)
            indices = catObj.selectIndices(entityId, "entity_id")
            if self.__verbose:
                self.__lfh.write("+PdbxIoUtils.getSourceDetailsList() found %d matches for entity %s in category %s \n" % (len(indices), entityId, sourceCategoryName))
            for ii in indices:
                d = {}

                scientificName = ""
                taxId = ""
                strain = ""
                beg = str(1)
                end = str(seqLength)
                entitySeqType = "Biological sequence"
                entityPartId = str(1)
                details = ""
                #
                scientificName = self.__getFirstValueFromList(["pdbx_gene_src_scientific_name", "pdbx_organism_scientific", "organism_scientific"], catObj, ii)
                strain = self.__getFirstValueFromList(
                    [
                        "gene_src_strain",
                        "strain",
                    ],
                    catObj,
                    ii,
                )
                taxId = self.__getFirstValueFromList(["pdbx_gene_src_ncbi_taxonomy_id", "pdbx_ncbi_taxonomy_id", "ncbi_taxonomy_id"], catObj, ii)
                commonName = self.__getFirstValueFromList(["gene_src_common_name", "common_name", "organism_common_name"], catObj, ii)
                variant = self.__getFirstValueFromList(["pdbx_gene_src_variant", "pdbx_variant"], catObj, ii)
                #
                geneName = self.__getValueOrDefault(catObj, "pdbx_gene_src_gene", ii, "")
                hostOrgSource = self.__getValueOrDefault(catObj, "pdbx_host_org_scientific_name", ii, "")
                hostOrgVector = self.__getValueOrDefault(catObj, "pdbx_host_org_vector", ii, "")
                hostOrgVectorType = self.__getValueOrDefault(catObj, "pdbx_host_org_vector_type", ii, "")
                hostOrgStrain = self.__getValueOrDefault(catObj, "pdbx_host_org_strain", ii, "")
                hostOrgTaxId = self.__getValueOrDefault(catObj, "pdbx_host_org_ncbi_taxonomy_id", ii, "")
                hostOrgPlasmidName = self.__getValueOrDefault(catObj, "plasmid_name", ii, "")
                hostOrgCommonName = self.__getValueOrDefault(catObj, "host_org_common_name", ii, "")
                hostOrgCellLine = self.__getValueOrDefault(catObj, "pdbx_host_org_cell_line", ii, "")
                hostOrgVariant = self.__getValueOrDefault(catObj, "pdbx_host_org_variant", ii, "")

                beg = self.__getValueOrDefault(catObj, "pdbx_beg_seq_num", ii, beg)
                end = self.__getValueOrDefault(catObj, "pdbx_end_seq_num", ii, end)
                #
                entitySeqType = self.__getValueOrDefault(catObj, "pdbx_seq_type", ii, entitySeqType)
                entityPartId = self.__getValueOrDefault(catObj, "pdbx_src_id", ii, entityPartId)
                details = self.__getValueOrDefault(catObj, "details", ii, details)

                d["SOURCE_NAME"] = scientificName
                d["SOURCE_COMMON_NAME"] = commonName
                d["SOURCE_TAXID"] = taxId
                d["SOURCE_STRAIN"] = strain
                d["SOURCE_DETAILS"] = details
                d["SOURCE_GENE_NAME"] = geneName
                d["SOURCE_VARIANT"] = variant

                d["HOST_ORG_SOURCE"] = hostOrgSource
                d["HOST_ORG_VECTOR"] = hostOrgVector
                d["HOST_ORG_VECTOR_TYPE"] = hostOrgVectorType
                d["HOST_ORG_STRAIN"] = hostOrgStrain
                d["HOST_ORG_TAXID"] = hostOrgTaxId
                d["HOST_ORG_PLASMID"] = hostOrgPlasmidName
                d["HOST_ORG_COMMON_NAME"] = hostOrgCommonName
                d["HOST_ORG_CELL_LINE"] = hostOrgCellLine
                d["HOST_ORG_VARIANT"] = hostOrgVariant

                d["SEQ_NUM_BEG"] = int(str(beg))
                d["SEQ_NUM_END"] = int(str(end))
                #
                d["SOURCE_NAME_ORIG"] = scientificName
                d["SOURCE_COMMON_NAME_ORIG"] = commonName
                d["SOURCE_TAXID_ORIG"] = taxId
                d["SOURCE_STRAIN_ORIG"] = strain
                d["SOURCE_VARIANT_ORIG"] = variant

                d["SEQ_NUM_BEG_ORIG"] = int(str(beg))
                d["SEQ_NUM_END_ORIG"] = int(str(end))

                d["SOURCE_GENE_NAME_ORIG"] = geneName
                d["HOST_ORG_SOURCE_ORIG"] = hostOrgSource
                d["HOST_ORG_VECTOR_ORIG"] = hostOrgVector
                d["HOST_ORG_VECTOR_TYPE_ORIG"] = hostOrgVectorType
                d["HOST_ORG_STRAIN_ORIG"] = hostOrgStrain
                d["HOST_ORG_TAXID_ORIG"] = hostOrgTaxId
                d["HOST_ORG_PLASMID_ORIG"] = hostOrgPlasmidName
                d["HOST_ORG_COMMON_NAME_ORIG"] = hostOrgCommonName
                d["HOST_ORG_VARIANT_ORIG"] = hostOrgVariant

                d["SEQ_PART_TYPE"] = entitySeqType
                d["SEQ_PART_ID"] = int(str(entityPartId))
                rL.append(d)
            #
            # reset the part id to 1 if there is only 1
            if len(rL) == 1:
                rL[0]["SEQ_PART_ID"] = 1
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)

        return rL

    def __getValueOrDefault(self, catObj, attrib, row, default):
        """Returns default value if attribute does not exist or row number out of range"""
        # API changed for lower code. getValueOrDefault used to return default for missing attribute.
        # Now raises exception.
        if catObj.hasAttribute(attrib):
            return catObj.getValueOrDefault(attrib, row, default)
        else:
            return default

    def assignSourceDefaultList(self, entityId, seqLength=0):  # pylint: disable=unused-argument
        """Assign default dictionaries containing source feature details.

        'SOURCE_NAME'   scientific name   =blank
        'SOURCE_TAXID'  NCBI Taxonomy ID  =blank
        'SOURCE_STRAIN' strain name       =blank

        'SEQ_NUM_BEG'   begining sequence offset (1-sequence length )
        'SEQ_NUM_END'   ending sequence offset

         *_ORIG'        copy of the above

        'SEQ_PART_TYPE' _entity_src_gen.pdbx_seq_type =  'Biological sequence'
        'SEQ_PART_ID'   _entity_src_gen.pdbx_src_id  = 1


        """
        rL = []
        try:
            d = {}
            beg = str(1)
            end = str(seqLength)
            entitySeqType = "Biological sequence"
            entityPartId = str(1)
            details = ""
            d["SOURCE_NAME"] = ""
            d["SOURCE_TAXID"] = ""
            d["SOURCE_STRAIN"] = ""
            d["SOURCE_DETAILS"] = details
            d["SOURCE_COMMON_NAME"] = ""
            d["SOURCE_VARIANT"] = ""
            d["SEQ_NUM_BEG"] = int(str(beg))
            d["SEQ_NUM_END"] = int(str(end))

            d["SOURCE_GENE_NAME"] = ""
            d["SOURCE_VARIANT"] = ""
            d["HOST_ORG_SOURCE"] = ""
            d["HOST_ORG_VECTOR"] = ""
            d["HOST_ORG_VECTOR_TYPE"] = ""
            d["HOST_ORG_STRAIN"] = ""
            d["HOST_ORG_TAXID"] = ""
            d["HOST_ORG_PLASMID"] = ""
            d["HOST_ORG_COMMON_NAME"] = ""
            d["HOST_ORG_CELL_LINE"] = ""
            d["HOST_ORG_VARIANT"] = ""

            d["SOURCE_GENE_NAME_ORIG"] = ""
            d["HOST_ORG_SOURCE_ORIG"] = ""
            d["HOST_ORG_VECTOR_ORIG"] = ""
            d["HOST_ORG_VECTOR_TYPE_ORIG"] = ""
            d["HOST_ORG_STRAIN_ORIG"] = ""
            d["HOST_ORG_TAXID_ORIG"] = ""
            d["HOST_ORG_PLASMID_ORIG"] = ""
            d["HOST_ORG_COMMON_NAME_ORIG"] = ""
            d["HOST_ORG_CELL_LINE_ORIG"] = ""
            d["HOST_ORG_VARIANT_ORIG"] = ""

            #
            d["SOURCE_NAME_ORIG"] = ""
            d["SOURCE_TAXID_ORIG"] = ""
            d["SOURCE_STRAIN_ORIG"] = ""
            d["SOURCE_COMMON_NAME_ORIG"] = ""
            d["SOURCE_VARIANT_ORIG"] = ""

            d["SEQ_NUM_BEG_ORIG"] = int(str(beg))
            d["SEQ_NUM_END_ORIG"] = int(str(end))
            d["SEQ_PART_TYPE"] = entitySeqType
            d["SEQ_PART_ID"] = int(str(entityPartId))
            rL.append(d)
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)

        return rL

    def getDepositorSeqDbRefDetails(self, entityId):
        """Returns details of depositor provided sequence database references and confilict annotation.

        Reads categories:  pdbx_struct_ref_seq_depositor_info and pdbx_struct_ref_seq_dif_depositor_info

        For pdbx_struct_ref_seq_depositor_info returns a list of dictionaries with keys -

        ['ref_id', 'entity_id', 'db_code','db_name',   'db_accession', 'db_align_beg', 'db_align_end', 'seq_align_begin', 'seq_align_end','seq_one_letter_code', 'details']

        For pdbx_struct_ref_seq_dif_depositor_info returns a dictionary with key [ref_id] corresponding to the alignment details for the input entityId.

        ['ref_id' ,'auth_mon_id' ,'auth_seq_id','db_code' ,'db_accession' ,'db_mon_id' ,'db_seq_id' ,'annotation' ,'ordinal' ]

        """
        refList = []
        difD = {}
        refIdList = []
        try:
            categoryName = "pdbx_struct_ref_seq_depositor_info"
            catObj = self.__currentContainer.getObj(categoryName)
            if catObj is None:
                return refList, difD
            if entityId is not None:
                indices = catObj.selectIndices(entityId, "entity_id")
            else:
                indices = []

            if self.__verbose:
                self.__lfh.write("+PdbxIoUtils.getDepositorSeqDbRefDetails() found %d matches for entity %s in category %s \n" % (len(indices), entityId, categoryName))

            for ii in indices:
                d1 = {}
                #
                attribList = [
                    "ref_id",
                    "entity_id",
                    "db_code",
                    "db_name",
                    "db_accession",
                    "db_align_beg",
                    "db_align_end",
                    "seq_align_begin",
                    "seq_align_end",
                    "seq_one_letter_code",
                    "details",
                ]

                for attribName in attribList:
                    if catObj.hasAttribute(attribName):
                        d1[attribName] = catObj.getValue(attribName, ii)
                    else:
                        d1[attribName] = None
                if "ref_id" in d1:
                    refIdList.append(d1["ref_id"])
                refList.append(d1)

            refIdList = list(set(refIdList))
            categoryName = "pdbx_struct_ref_seq_dif_depositor_info"
            catObj = self.__currentContainer.getObj(categoryName)
            if catObj is None:
                return refList, difD

            for refId in refIdList:
                difList = []
                indices = catObj.selectIndices(refId, "ref_id")
                if self.__verbose:
                    self.__lfh.write("+PdbxIoUtils.getDepositorSeqDbRefDetails() found %d matches for refId %s in category %s \n" % (len(indices), refId, categoryName))
                for ii in indices:
                    d2 = {}
                    #
                    attribList = ["ref_id", "auth_mon_id", "auth_seq_id", "db_code", "db_accession", "db_mon_id", "db_seq_id", "annotation", "ordinal"]
                    for attribName in attribList:
                        if catObj.hasAttribute(attribName):
                            d2[attribName] = catObj.getValue(attribName, ii)
                        else:
                            d2[attribName] = None
                    difList.append(d2)

                difD[refId] = difList
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+PdbxIoUtils.getDepositorSeqDbRefDetails() failing for entityId %s\n" % entityId)
                traceback.print_exc(file=self.__lfh)
        #
        return refList, difD

    def getArchiveSeqDbRefDetails(self, entityId):
        """Returns details of depositor provided sequence database references and confilict annotation.

        Reads categories:  struct_ref, struct_ref_seq and struct_ref_seq_dif -

        category struct_ref -
          ['id', 'db_name', 'db_code', 'pdbx_db_accession', 'entity_id', 'pdbx_seq_one_letter_code', 'pdbx_align_begin', 'biol_id']

        category struct_ref_seq -

          ['align_id', 'ref_id','pdbx_PDB_id_code', 'pdbx_strand_id','seq_align_beg','pdbx_seq_align_beg_ins_code', 'seq_align_end', 'pdbx_seq_align_end_ins_code', 'pdbx_db_accession',
            'db_align_beg','db_align_end', 'pdbx_auth_seq_align_beg', 'pdbx_auth_seq_align_end']

        category struct_ref_seq_dif -

        ['align_id', 'pdbx_pdb_id_code', 'mon_id', 'pdbx_pdb_strand_id', 'seq_num', 'pdbx_pdb_ins_code', 'pdbx_seq_db_name',
        'pdbx_seq_db_accession_code', 'db_mon_id', 'pdbx_seq_db_seq_num',
        'details', 'pdbx_auth_seq_num', 'pdbx_ordinal']

        """
        refList = []
        alignD = {}
        difD = {}
        refIdList = []
        try:
            categoryName = "struct_ref"
            catObj = self.__currentContainer.getObj(categoryName)
            if catObj is None:
                return refList, alignD, difD
            indices = catObj.selectIndices(entityId, "entity_id")
            if self.__verbose:
                self.__lfh.write("+PdbxIoUtils.getArchiveSeqDbRefDetails() found %d matches for entity %s in category %s \n" % (len(indices), entityId, categoryName))
            for ii in indices:
                d1 = {}
                #
                attribList = ["id", "db_name", "db_code", "pdbx_db_accession", "entity_id", "pdbx_seq_one_letter_code", "pdbx_align_begin", "biol_id"]
                for attribName in attribList:
                    if catObj.hasAttribute(attribName):
                        d1[attribName] = catObj.getValue(attribName, ii)
                    else:
                        d1[attribName] = None
                if "id" in d1:
                    refIdList.append(d1["id"])
                refList.append(d1)

            refIdList = list(set(refIdList))

            categoryName = "struct_ref_seq"
            catObj = self.__currentContainer.getObj(categoryName)
            if catObj is None:
                return refList, alignD, difD

            #
            #  Get the list of alignments --
            #
            alignIdList = []
            for refId in refIdList:
                alignList = []
                indices = catObj.selectIndices(refId, "ref_id")
                if self.__verbose:
                    self.__lfh.write("+PdbxIoUtils.getArchiveSeqDbRefDetails() found %d matches for refId %s in category %s \n" % (len(indices), refId, categoryName))
                for ii in indices:
                    d2 = {}
                    #
                    attribList = [
                        "align_id",
                        "ref_id",
                        "pdbx_PDB_id_code",
                        "pdbx_strand_id",
                        "seq_align_beg",
                        "pdbx_seq_align_beg_ins_code",
                        "seq_align_end",
                        "pdbx_seq_align_end_ins_code",
                        "pdbx_db_accession",
                        "db_align_beg",
                        "db_align_end",
                        "pdbx_auth_seq_align_beg",
                        "pdbx_auth_seq_align_end",
                    ]
                    for attribName in attribList:
                        if catObj.hasAttribute(attribName):
                            d2[attribName] = catObj.getValue(attribName, ii)
                        else:
                            d2[attribName] = None
                    if "align_id" in d2:
                        alignIdList.append(d2["align_id"])
                    alignList.append(d2)
                alignD[refId] = alignList

            #
            #  Get the list of residue conflicts in the alignment list --
            #
            categoryName = "struct_ref_seq_dif"
            catObj = self.__currentContainer.getObj(categoryName)
            if catObj is None:
                return refList, alignD, difD

            alignIdList = list(set(alignIdList))

            for alignId in alignIdList:
                difList = []
                indices = catObj.selectIndices(alignId, "align_id")
                if self.__verbose:
                    self.__lfh.write("+PdbxIoUtils.getArchiveSeqDbRefDetails() found %d matches for alignId %s in category %s \n" % (len(indices), alignId, categoryName))
                for ii in indices:
                    d2 = {}
                    #
                    attribList = [
                        "align_id",
                        "pdbx_pdb_id_code",
                        "mon_id",
                        "pdbx_pdb_strand_id",
                        "seq_num",
                        "pdbx_pdb_ins_code",
                        "pdbx_seq_db_name",
                        "pdbx_seq_db_accession_code",
                        "db_mon_id",
                        "pdbx_seq_db_seq_num",
                        "details",
                        "pdbx_auth_seq_num",
                        "pdbx_ordinal",
                    ]
                    for attribName in attribList:
                        if catObj.hasAttribute(attribName):
                            d2[attribName] = catObj.getValue(attribName, ii)
                        else:
                            d2[attribName] = None
                    difList.append(d2)

                difD[alignId] = difList

        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+PdbxIoUtils.getArchiveSeqDbRefDetails() failing for entityId %s\n" % entityId)
                traceback.print_exc(file=self.__lfh)
        #
        return refList, alignD, difD

    def getEntityPolyList(self):
        """Returns a list of polymer entity ids"""
        if not self.__currentContainer.exists("entity_poly"):
            return []
        #
        catObj = self.__currentContainer.getObj("entity_poly")
        nRows = catObj.getRowCount()

        eList = []
        for ii in range(0, nRows):
            eId = catObj.getValue("entity_id", ii)
            eList.append(eId)
        return eList

    def getSequence1Xyz(self, chainId):
        """
        Get the one-letter-code sequence from the deposited coordinate data
        """
        if not self.__currentContainer.exists("pdbx_poly_seq_scheme"):
            return ""
        #
        catObj = self.__currentContainer.getObj("pdbx_poly_seq_scheme")
        indices = catObj.selectIndices(chainId, "pdb_strand_id")
        mon3List = []
        for ii in indices:
            mon3List.append(catObj.getValue("pdb_mon_id", ii))

        # translate to one-letter-codes --
        sq = ""
        for mon in mon3List:
            if mon in SequenceReferenceData._monDict3:  # pylint: disable=protected-access
                sq += SequenceReferenceData._monDict3[mon]  # pylint: disable=protected-access
            elif (mon in [".", "?"]) or (len(mon) == 0):
                sq += "-"
            else:
                sq += "X"

        return sq

    ##

    def getMisalignCount(self, chainId):
        """
        Return the count of misaligned residues in maxit sequence alignment.

        """
        catObj = self.__currentContainer.getObj("pdbx_poly_seq_scheme")
        indices = catObj.selectIndices(chainId, "pdb_id")
        mCount = 0
        for ii in indices:
            pdb_mon_id = catObj.getValue("pdb_mon_id", ii)
            mon_id = catObj("mon_id", ii)
            if pdb_mon_id != mon_id:
                mCount += 1
        return mCount

    def getAtomSiteTableLength(self):
        if self.__currentContainer.exists("atom_site"):
            catObj = self.__currentContainer.getObj("atom_site")
            return catObj.getRowCount()
        else:
            return 0

    def getEncapsulatedCoordinates(self):
        """_pdbx_original_pdb_coordinates.coord_section"""
        try:
            row = [""]  # Keep pylint happy.
            encapCoordTable = [""]  # For pylint
            if self.__currentContainer.exists("pdbx_original_pdb_coordinates"):
                encapCoordTable = self.__currentContainer.getObj("pdbx_original_pdb_coordinates")
            elif self.__currentContainer.exists("ndb_original_pdb_coordinates"):
                encapCoordTable = self.__currentContainer.getObj("ndb_original_pdb_coordinates")
            elif self.__currentContainer.exists("ndb_original_ndb_coordinates"):
                encapCoordTable = self.__currentContainer.getObj("ndb_original_ndb_coordinates")
            row = encapCoordTable.getRow(0)
        except:  # noqa: E722 pylint: disable=bare-except
            row = [""]

        return row[0]

    #
    #
    def getPolymerEntityChainDict(self):
        if len(self.__polymerEntityChainDict) == 0:
            self.__buildPolymerEntityChainDict()
        return self.__polymerEntityChainDict

    def getChainPolymerEntityDict(self):
        if len(self.__chainPolymerEntityDict) == 0:
            self.__buildPolymerEntityChainDict()
        return self.__chainPolymerEntityDict

    def __buildPolymerEntityChainDict(self):
        """Build entity chain mapping information --  Chain details must be provided"""
        self.__polymerEntityChainDict = {}
        pEntityList = self.getPolymerEntityList("all")
        for eId in pEntityList:
            tL = self.getPdbChainIdList(eId)
            if len(tL) > 0:
                self.__polymerEntityChainDict[eId] = tL
        #
        self.__chainPolymerEntityDict = {}
        for eId, cList in self.__polymerEntityChainDict.items():
            for cId in cList:
                self.__chainPolymerEntityDict[cId] = eId

    def getEntitySequence1Auth(self, kd):
        """Get the dictoinary of author provided one-letter-code sequences corresponding
        to the input keyword arguments:

        kd  -->   entityId=1|2...
                  chainId=A|B...

        example:  getEntitySequenceAuth1(chainId="B")
                  getEntitySequenceAuth1(entityId='1')
        """
        sq = None
        for k, v in kd.items():
            eId = None
            if k == "entityId":
                eId = v
            elif k == "chainId":
                if v in self.__chainPolymerEntityDict:
                    eId = self.__chainPolymerEntityDict[v]

            if eId is not None:
                try:
                    catObj = self.__currentContainer.getObj("entity_poly")
                    indices = catObj.selectIndices(eId, "entity_id")
                    sqR = catObj.getValue("pdbx_seq_one_letter_code", indices[0])
                    sq = ""
                    for cc in sqR:
                        if cc.isspace():
                            continue
                        sq += cc
                    return sq
                except:  # noqa: E722 pylint: disable=bare-except
                    pass

        return sq

    def getEntityCanSequence1Auth(self, kd):
        """Get the dictionary of author provided cannonical one-letter-code sequences corresponding
        to the input keyword arguments:
        kd  -->   entityId=1|2...
                  chainId=A|B...
        example:  getEntityCanSequenceAuth1(chainId="B")
                  getEntityCanSequenceAuth1(entityId='1',entityId='2')
        """
        sq = None
        for k, v in kd.items():
            eId = None
            if k == "entityId":
                eId = v
            elif k == "chainId":
                if v in self.__chainPolymerEntityDict:
                    eId = self.__chainPolymerEntityDict[v]

            if eId is not None:
                catObj = self.__currentContainer.getObj("entity_poly")
                indices = catObj.selectIndices(eId, "entity_id")
                sqR = catObj.getValue("pdbx_seq_one_letter_code_can", indices[0])
                sq = ""
                for cc in sqR:
                    if cc.isspace():
                        continue
                    sq += cc
                return sq
        return sq

    def getAssemblyDetails(self):
        """
        #
        loop_
        _pdbx_struct_assembly.id
        _pdbx_struct_assembly.details
        _pdbx_struct_assembly.method_details
        _pdbx_struct_assembly.oligomeric_count
        1 author_defined_assembly   ?    3
        2 software_defined_assembly PISA 4
        #
        loop_
        _pdbx_struct_assembly_gen.assembly_id
        _pdbx_struct_assembly_gen.oper_expression
        _pdbx_struct_assembly_gen.asym_id_list
        1 1       A,B,C,D,E,F,G,H,I,J,K
        1 2       A,B,C,D,E,F,G,H,I,J,K
        1 3       A,B,C,D,E,F,G,H,I,J,K
        2 1,4,5,6 A,B,C,D,E,F,G,H,I,J,K
        #
        loop_
        _pdbx_struct_assembly_prop.biol_id
        _pdbx_struct_assembly_prop.type
        _pdbx_struct_assembly_prop.value
        _pdbx_struct_assembly_prop.details
        2 "SSA (A^2)"  58580 ?
        2 "ABSA (A^2)" 31160 ?
        2 MORE         -337  ?
        #
        loop_
        _pdbx_struct_oper_list.id
        _pdbx_struct_oper_list.type
        _pdbx_struct_oper_list.name
        _pdbx_struct_oper_list.matrix[1][1]
        _pdbx_struct_oper_list.matrix[1][2]
        _pdbx_struct_oper_list.matrix[1][3]
        _pdbx_struct_oper_list.vector[1]
        _pdbx_struct_oper_list.matrix[2][1]
        _pdbx_struct_oper_list.matrix[2][2]
        _pdbx_struct_oper_list.matrix[2][3]
        _pdbx_struct_oper_list.vector[2]
        _pdbx_struct_oper_list.matrix[3][1]
        _pdbx_struct_oper_list.matrix[3][2]
        _pdbx_struct_oper_list.matrix[3][3]
        _pdbx_struct_oper_list.vector[3]
        1 "identity operation"         1_555 1.0000000000  0.0000000000 0.0000000000 0.0000000000   0.0000000000 1.0000000000  0.0000000000 0.0000000000  0.0000000000 0.0000000000 1.0000000000  0.0000000000  # noqa: E501
        2 "crystal symmetry operation" 2_566 -1.0000000000 0.0000000000 0.0000000000 0.0000000000   0.0000000000 -1.0000000000 0.0000000000 95.9710000000 0.0000000000 0.0000000000 1.0000000000  137.3490000000  # noqa: E501
        3 "crystal symmetry operation" 2_656 -1.0000000000 0.0000000000 0.0000000000 80.9760000000  0.0000000000 -1.0000000000 0.0000000000 0.0000000000  0.0000000000 0.0000000000 1.0000000000  137.3490000000  # noqa: E501
        4 "crystal symmetry operation" 2_765 -1.0000000000 0.0000000000 0.0000000000 161.9520000000 0.0000000000 -1.0000000000 0.0000000000 95.9710000000 0.0000000000 0.0000000000 1.0000000000  0.0000000000  # noqa: E501
        5 "crystal symmetry operation" 3_757 -1.0000000000 0.0000000000 0.0000000000 161.9520000000 0.0000000000 1.0000000000  0.0000000000 0.0000000000  0.0000000000 0.0000000000 -1.0000000000 274.6980000000  # noqa: E501
        6 "crystal symmetry operation" 4_567 1.0000000000  0.0000000000 0.0000000000 0.0000000000   0.0000000000 -1.0000000000 0.0000000000 95.9710000000 0.0000000000 0.0000000000 -1.0000000000 274.6980000000  # noqa: E501
        #
        """
        assemL = []
        assemGenL = []
        assemOpL = []
        if not self.__currentContainer.exists("pdbx_struct_assembly"):
            return assemL, assemGenL, assemOpL
        #
        catObj = self.__currentContainer.getObj("pdbx_struct_assembly")
        myList = ["id", "details", "method_details", "oligomeric_count"]
        assemL = self.__getAttributeDictList(catObj=catObj, attributeList=myList)
        #
        catObj = self.__currentContainer.getObj("pdbx_struct_assembly_gen")
        myList = ["assembly_id", "oper_expression", "asym_id_list"]
        assemGenL = self.__getAttributeDictList(catObj=catObj, attributeList=myList)

        catObj = self.__currentContainer.getObj("pdbx_struct_oper_list")
        myList = ["id", "type", "name"]
        assemOpL = self.__getAttributeDictList(catObj=catObj, attributeList=myList)

        return assemL, assemGenL, assemOpL

    def __getAttributeDictList(self, catObj, attributeList):
        #
        rList = []
        try:
            colNames = list(catObj.getAttributeList())
            nRows = catObj.getRowCount()
            for iRow in range(0, nRows):
                rD = {}
                row = catObj.getRow(iRow)
                for col in attributeList:
                    if col in colNames:
                        val = str(row[colNames.index(col)])
                        if val is None:
                            val = ""
                        elif (val == ".") or (val == "?"):
                            val = ""
                        rD[col] = val
                    else:
                        rD[col] = ""
                rList.append(rD)
        except:  # noqa: E722 pylint: disable=bare-except
            self.__lfh.write("PdbxIoUtils.__getAttributeDictList - failed ")
            traceback.print_exc(file=self.__lfh)

        return rList

    def __getDepositorDetails(self, tableName, myList):
        """Returns a dictionary of assembly details using a list of attributes"""

        if not self.__currentContainer.exists(tableName):
            return []

        catObj = self.__currentContainer.getObj(tableName)

        return self.__getAttributeDictList(catObj, myList)

    def getDepositorAssemblyDetails(self):
        """Returns a dictionary of assembly details provided at deposition."""
        myList = ["id", "details", "matrix_flag", "method_details", "oligomeric_count", "oligomeric_details", "upload_file_name"]

        return self.__getDepositorDetails("pdbx_struct_assembly_depositor_info", myList)

    def getDepositorAssemblyDetailsRcsb(self):
        """Returns a dictionary of assembly details provided at deposition (for current system)"""
        myList = ["id", "details", "rcsb_description", "method_details", "pdbx_aggregation_state", "pdbx_assembly_method", "pdbx_formula_weight", "pdbx_formula_weight_method"]

        return self.__getDepositorDetails("struct_biol", myList)

    def getDepositorAssemblyGen(self):
        """Returns a dictionary of assembly details provided at deposition (for current system)"""
        myList = ["id", "asym_id_list", "assembly_id", "oper_expression", "full_matrices", "at_unit_matrix", "chain_id_list", "all_chains", "helical_rotation", "helical_rise"]

        return self.__getDepositorDetails("pdbx_struct_assembly_gen_depositor_info", myList)

    def getDepositorStructOperList(self):
        """Returns a dictionary of _pdbx_struct_oper_list_depositor_info."""
        myList = [
            "id",
            "name",
            "symmetry_operation",
            "type",
            "matrix[1][1]",
            "matrix[1][2]",
            "matrix[1][3]",
            "matrix[2][1]",
            "matrix[2][2]",
            "matrix[2][3]",
            "matrix[3][1]",
            "matrix[3][2]",
            "matrix[3][3]",
            "vector[1]",
            "vector[2]",
            "vector[3]",
        ]

        return self.__getDepositorDetails("pdbx_struct_oper_list_depositor_info", myList)

    def getDepositorAssemblyEvidence(self):
        """Returns a dictionary of _pdbx_struct_oper_list_depositor_info."""
        myList = ["id", "assembly_id", "experimental_support", "details"]

        return self.__getDepositorDetails("pdbx_struct_assembly_auth_evidence", myList)

    def getDepositorAssemblyClassification(self):
        """Returns a dictionary of _pdbx_struct_assembly_auth_classification"""
        myList = ["assembly_id", "reason_for_interest"]

        return self.__getDepositorDetails("pdbx_struct_assembly_auth_classification", myList)

    def __getLeastCommmonComp(self, compIdList):
        srd = SequenceReferenceData(verbose=self.__verbose, log=self.__lfh)
        try:
            for compId in compIdList:
                if not srd.isStandard3(compId):
                    return compId
            return compIdList[0]
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            return "UNK"

    def getSequenceFeaturesFromAtomSite(self, linkInstD=None):
        """Extract the polymer sequence from the PDBx atom_site category for the first model.

        Returns -

        chD{} =[(a3, orig auth index + ins code (str), comment/details,  align index placeholder),(),..]

        as a dictionary with chain id key containing sequences stored as a list of tuples.

        Each tuple - (3-letter-code, original auth residue index w/ins.code (str), comment/details, alignment index))

        On input the link instance dictionary as structure -

         linkIsntD[authAsymId,compId,authSeqId+insCode] = (compId,authSeqId+insCode,float(dist),seqId,'begin|end',longFlag)
        """

        if not self.__currentContainer.exists("atom_site"):
            return {}

        linkD = linkInstD if linkInstD is not None else {}

        catObj = self.__currentContainer.getObj("atom_site")
        nRows = catObj.getRowCount()
        #
        colNames = list(catObj.getAttributeList())
        # i0 = colNames.index("group_PDB")
        i1 = colNames.index("auth_asym_id")

        if "auth_comp_id" in colNames:
            i2 = colNames.index("auth_comp_id")
        else:
            i2 = colNames.index("label_comp_id")

        i3 = colNames.index("auth_seq_id")
        #
        i3L = colNames.index("label_seq_id")
        #
        i4 = 0  # For pylint - should never happen
        if "pdbx_PDB_model_num" in colNames:
            i4 = colNames.index("pdbx_PDB_model_num")
        elif "ndb_model" in colNames:
            i4 = colNames.index("ndb_model")
        #
        i8 = 0  # Should never happen - keep pylint happy
        if "pdbx_PDB_ins_code" in colNames:
            i8 = colNames.index("pdbx_PDB_ins_code")
        elif "ndb_ins_code" in colNames:
            i8 = colNames.index("ndb_ins_code")

        #
        i5 = colNames.index("occupancy")
        i6 = colNames.index("label_alt_id")
        occCount = 0
        occSum = 0.0
        disorderFlag = False
        #
        compHeteroFlag = False
        compHeteroL = []

        #
        chD = {}
        rList = []
        idx = 1
        row = catObj.getRow(0)
        # atGroup = row[i0]
        chId_cur = row[i1]
        compId_cur = row[i2]
        seqId_cur = row[i3]
        modelId_cur = row[i4]
        #
        occupancy = row[i5]
        altId = row[i6]
        occCount += 1
        occSum += float(str(occupancy))
        if len(altId) > 0 and altId not in [".", "?"]:
            disorderFlag = True
        #
        insCode_cur = row[i8]

        #
        for iRow in range(1, nRows):
            row = catObj.getRow(iRow)
            # atGroup = row[i0]
            seqIdL = row[i3L]
            if seqIdL in [".", "?"]:
                continue
            #
            # restore the atom group filter 18-Sep-2013 to handle abberant data sets.
            # if atGroup != 'ATOM':
            #     continue
            chId = row[i1]
            compId = row[i2]
            seqId = row[i3]
            modelId = row[i4]
            insCode = row[i8]
            #
            # Only read the first model ---
            #
            if modelId != modelId_cur:
                break

            if (compId != compId_cur) and (seqId == seqId_cur) and (insCode == insCode_cur) and (chId == chId_cur):
                #             -- Special case of sequence heterogeneity --
                compHeteroFlag = True
                if compId_cur not in compHeteroL:
                    compHeteroL.append(compId_cur)
                if compId not in compHeteroL:
                    compHeteroL.append(compId)
                compId_cur = compId
            ##
            if (compId == compId_cur) and (seqId == seqId_cur) and (insCode == insCode_cur) and (chId == chId_cur):
                occupancy = row[i5]
                altId = row[i6]
                occCount += 1
                occSum += float(str(occupancy))
                if len(altId) > 0 and altId not in [".", "?"]:
                    disorderFlag = True
                continue
            else:
                #  Author residue numbers include appended insertion codes --
                seqIdInsCode = seqId_cur
                if len(insCode_cur) > 0 and insCode_cur not in [".", "?"]:
                    seqIdInsCode = seqId_cur + insCode_cur

                # Flag disorder
                disorderText = "disordered" if disorderFlag else ""

                # Occupancy outliers
                if occCount > 0:
                    meanOcc = occSum / occCount
                    occText = "" if meanOcc > 0.75 else "mean_occ=%.3f" % meanOcc
                else:
                    if self.__verbose:
                        self.__lfh.write("  ++ occCount %d occSum %f idx %d disorderText %s \n" % (occCount, occSum, idx, disorderText))
                        self.__lfh.write("  ++     chId %s compId %s seqId %s modelId %s\n" % (chId, compId, seqId, modelId))
                    occText = ""
                #
                # Polymer linkage outliers ---
                #
                tTup = (chId, compId_cur, seqIdInsCode)
                if tTup in linkD:
                    #  This captures link=x.xx for the outlier
                    linkText = "link=%.2f" % linkD[tTup][2]
                    #
                    # Flag the first partner in a long linkage for sequence alignment
                    if linkD[tTup][4] == "begin" and linkD[tTup][5]:
                        linkText += ",long_begin"
                    elif linkD[tTup][4] == "end" and linkD[tTup][5]:
                        linkText += ",long_end"
                else:
                    linkText = ""
                #
                if compHeteroFlag:
                    lcId = self.__getLeastCommmonComp(compHeteroL)
                    heteroText = "hetero-" + ":".join(compHeteroL)
                    comment = str(occText + "  " + disorderText + " " + linkText + " " + heteroText).strip()
                    rList.append((lcId, seqIdInsCode, comment, idx))
                else:
                    comment = str(occText + "  " + disorderText + " " + linkText).strip()
                    rList.append((compId_cur, seqIdInsCode, comment, idx))

                idx += 1
                compId_cur = compId
                seqId_cur = seqId
                occCount = 1
                occSum = float(str(occupancy))
                disorderFlag = False
                insCode_cur = insCode
                #
                compHeteroFlag = False
                compHeteroL = []

            #
            if chId != chId_cur:
                if chId_cur not in chD:
                    chD[chId_cur] = rList
                rList = []
                idx = 1
                chId_cur = chId

        # assign the last residue and chain
        if chId_cur not in chD:
            disorderText = "disordered" if disorderFlag else ""
            if occCount > 0:
                meanOcc = occSum / occCount
                occText = "" if meanOcc > 0.75 else "mean_occ=%.3f" % meanOcc
            else:
                self.__lfh.write("  ++ occCount %d occSum %f idx %d disorderText %s \n" % (occCount, occSum, idx, disorderText))
                self.__lfh.write("  ++     chId %s compId %s seqId %s modelId %s\n" % (chId, compId, seqId, modelId))
                occText = ""

            seqIdInsCode = seqId_cur
            if len(insCode_cur) > 0 and insCode_cur not in [".", "?"]:
                seqIdInsCode = seqId_cur + insCode_cur

            tTup = (chId, compId_cur, seqIdInsCode)
            if tTup in linkD:
                linkText = "link=%.2f" % linkD[tTup][2]
                # Flag the first partner in a long linkage for sequence alignment
                if linkD[tTup][4] == "begin" and linkD[tTup][5]:
                    linkText += ",long_begin"
                elif linkD[tTup][4] == "end" and linkD[tTup][5]:
                    linkText += ",long_end"

            else:
                linkText = ""

            if compHeteroFlag:
                lcId = self.__getLeastCommmonComp(compHeteroL)
                heteroText = "hetero-" + ":".join(compHeteroL)
                comment = str(occText + "  " + disorderText + " " + linkText + " " + heteroText).strip()
                rList.append((lcId, seqIdInsCode, comment, idx))
            else:
                comment = str(occText + "  " + disorderText + " " + linkText).strip()
                rList.append((compId_cur, seqIdInsCode, comment, idx))

            chD[chId_cur] = rList
        #
        return chD


if __name__ == "__main__":
    pass
