# There are some "global" variables and cannot determine if we "own" them
# pylint: disable=attribute-defined-outside-init
##
# File:    SequenceLabel.py
# Date:    16-Dec-2009
#
# Updates:
#  15-Jan-2010 jdw Refactor
#  18-Jan-2010 jdw Add sequence label
#  18-Jan-2010 jdw Add sequence version identifier
#  19-Jan-2010 jdw Add alternative sequence identifier
#  12-Feb-2010 jdw Add class SequenceFeature()
#  20-Apr-2010 jdw Port to module seqmodule.
#  02-May-2010 jdw Add accessors for source features.
#  27-Feb-2013 jdw Add support for polymer parts.
#                  Improve type specific attribute handling.
#  03-Mar-2013 jdw add method setAuthPartDetails()
#  03-Mar-2013 jdw move markup methods to a SequenceFeatureDepict()
#  04-Mar-2013 jdw add data items for capture ordering information
#  28-Mar-2013 jdw add feature POLYMER_LINKING_TYPE
#  19-Apr-2013 jdw strip surrounding whitepsace when decoding strain.
#  18-Sep-2013 jdw fill in missing placeholders for original values,
#                  add additional reference sequence details, and additonal
#                  sample source details.
#  04-Nov-2013 jdw update type of ref_sort_metric
#  06-Nov-2013 jdw add source and host org common names
#  11-Nov-2013 jdw update accessors for author provided values
#  22-Nov-2013 jdw add mapping to reset to original author provided values
#   5-Dec-2013 jdw add operator methods to SequenceLabel() class
#  19-Jan-2014 jdw add host org cell line methods
#  09-Feb-2014 jdw added SOURCE_VARIANT, HOST_ORG_VECTOR_TYPE and orig analogs
#  20-Feb-2014 jdw updated regex for deconvuluting source and strain.
#  29-Jul-2014 jdw add 'DB_ISOFORM_DESCRIPTION',
#  10-Sep-2014 jdw add  method getPartInfo()
#  25-Nov-2014 jdw return source + strain from __getstrain()
#  01-Feb-2015 jdw add host org variant
#  07-Sep-2017 zf  add ALIGN_LENGTH, ANNO_EDIT_DB_NAME, ANNO_EDIT_DB_CODE, ANNO_EDIT_DB_ACCESSION, ANNO_EDIT_DB_ALIGN_BEGIN,
#                  ANNO_EDIT_DB_ALIGN_END
#                  updated updateAuth() to keep original author's molecular name for DNA/RNA entities if existing
#  02-Sep-2020 zf  Excluded 'Uncharacterized protein' molecule name from Uniprot
#  03-Oct-2022 zf  Excluded "SOURCE_STRAIN" from updateAuth() method.
#                  add IS_AUTH_PROVIDED_ID
##
"""
Containers for sequence and residue labels/features used as identifiers and classifiers
of aligned sequence data.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.07"

import sys
import re
import traceback

from wwpdb.apps.seqmodule.io.TaxonomyDbUtils import TaxonomyDbUtils


class SequenceFeatureMap(object):
    def __init__(self, verbose=False, log=sys.stderr):
        self.__verbose = verbose
        self.__lfh = log
        self.__debug = False

    def updateAuth(self, authFD, refFD):
        """Map editable reference values to current author values --"""
        mapTupList = [
            ("ENTITY_DESCRIPTION", "DB_MOLECULE_NAME"),
            ("ENTITY_SYNONYMS", "DB_MOLECULE_SYNONYMS"),
            ("SOURCE_GENE_NAME", "DB_GENE_NAME"),
            # ('SOURCE_TAXID', 'SOURCE_TAXID'),
            # ('SOURCE_ORGANISM', 'SOURCE_ORGANISM'),
            # ("SOURCE_STRAIN", "SOURCE_STRAIN"),
            # ('SOURCE_COMMON_NAME', 'SOURCE_COMMON_NAME'),
            ("ENTITY_ENZYME_CLASS", "DB_MOLECULE_EC"),
        ]
        if self.__debug:
            for mapTup in mapTupList:
                self.__lfh.write("++++Before Mapping %30s : %-40r  %30s:  %r\n" % (mapTup[0], authFD[mapTup[0]], mapTup[1], refFD[mapTup[0]]))
        #
        if ("SOURCE_TAXID" in authFD) and authFD["SOURCE_TAXID"]:
            taxInfoUtil = TaxonomyDbUtils(verbose=self.__verbose, log=self.__lfh)
            scientific_name, common_name = taxInfoUtil.getTaxonomyNames(authFD["SOURCE_TAXID"])
            if scientific_name:
                authFD["SOURCE_ORGANISM"] = scientific_name
                authFD["SOURCE_COMMON_NAME"] = common_name
            #
        #
        for mapTup in mapTupList:
            if mapTup[0] == "ENTITY_DESCRIPTION":
                updateFlag = False
                # Keep the original author provided name for DNA/RNA entities
                if ("DB_NAME" in refFD) and (refFD["DB_NAME"] in ["GB", "DBJ", "EMB", "EMBL", "REF"]):
                    if ((mapTup[0] not in authFD) or (not authFD[mapTup[0]])) and (mapTup[1] in refFD) and len(refFD[mapTup[1]]) > 1:
                        authFD[mapTup[0]] = refFD[mapTup[1]]
                        updateFlag = True
                    #
                elif refFD[mapTup[1]] is not None and len(refFD[mapTup[1]]) > 1:
                    if (refFD[mapTup[1]].strip().upper() != "UNCHARACTERIZED PROTEIN") and (refFD[mapTup[1]].strip().upper() != "PREDICTED PROTEIN"):
                        authFD[mapTup[0]] = refFD[mapTup[1]]
                        updateFlag = True
                    #
                #
                if not updateFlag:
                    if (authFD["ENTITY_DESCRIPTION_ORIG"] is not None) and (len(authFD["ENTITY_DESCRIPTION_ORIG"]) > 1):
                        authFD[mapTup[0]] = authFD["ENTITY_DESCRIPTION_ORIG"]
                    else:
                        authFD[mapTup[0]] = ""
                    #
                #
            elif refFD[mapTup[1]] is not None and len(refFD[mapTup[1]]) > 1:
                authFD[mapTup[0]] = refFD[mapTup[1]]
            else:
                authFD[mapTup[0]] = ""
            #
        #
        if ("REF_SEQ_FRAGMENT_DETAILS" in refFD) and (len(refFD["REF_SEQ_FRAGMENT_DETAILS"]) > 1):
            authFD["ENTITY_FRAGMENT_DETAILS"] = refFD["REF_SEQ_FRAGMENT_DETAILS"]
        #
        if self.__debug:
            for mapTup in mapTupList:
                self.__lfh.write("++++After %30s  %r\n" % (mapTup[0], authFD[mapTup[0]]))
        #
        return True

    def updateAuthOrig(self, authFD):
        """Map editable original author values to current author values --"""
        mapTupList = [
            ("ENTITY_DESCRIPTION", "ENTITY_DESCRIPTION_ORIG"),
            ("ENTITY_SYNONYMS", "ENTITY_SYNONYMS_ORIG"),
            ("SOURCE_GENE_NAME", "SOURCE_GENE_NAME_ORIG"),
            ("SOURCE_TAXID", "SOURCE_TAXID_ORIG"),
            ("SOURCE_ORGANISM", "SOURCE_ORGANISM_ORIG"),
            ("SOURCE_STRAIN", "SOURCE_STRAIN_ORIG"),
            ("SOURCE_COMMON_NAME", "SOURCE_COMMON_NAME_ORIG"),
            ("ENTITY_ENZYME_CLASS", "ENTITY_ENZYME_CLASS_ORIG"),
        ]
        if self.__debug:
            for mapTup in mapTupList:
                self.__lfh.write("++++Before Mapping %30s : %-40r  %30s:  %r\n" % (mapTup[0], authFD[mapTup[0]], mapTup[1], authFD[mapTup[0]]))
        #
        for mapTup in mapTupList:
            if authFD[mapTup[1]] is not None and len(authFD[mapTup[1]]) > 1:
                authFD[mapTup[0]] = authFD[mapTup[1]]
            else:
                authFD[mapTup[0]] = ""
        #
        if self.__debug:
            for mapTup in mapTupList:
                self.__lfh.write("++++After %30s  %r\n" % (mapTup[0], authFD[mapTup[0]]))
        #
        return True


class SequenceFeature(object):

    """Manages access to sequence feature data -

    Storage model is dictionary of key value pairs -

    Content include database identifiers, source details, and sequence comparison statistics.
    """

    def __init__(self, verbose=False, log=sys.stderr):
        self.__verbose = verbose
        self.__lfh = log
        self.__debug = False
        #
        self.__fListStr = [
            "DB_NAME",
            "DB_CODE",
            "DB_ACCESSION",
            "DB_ISOFORM",
            "DB_ISOFORM_DESCRIPTION",
            "ANNO_EDIT_DB_NAME",
            "ANNO_EDIT_DB_CODE",
            "ANNO_EDIT_DB_ACCESSION",
            "ANNO_EDIT_DB_ALIGN_BEGIN",
            "ANNO_EDIT_DB_ALIGN_END",
            "SOURCE_ORGANISM",
            "SOURCE_STRAIN",
            "SOURCE_TAXID",
            "SOURCE_GENE_NAME",
            "SOURCE_VARIANT",
            "SOURCE_METHOD",
            "POLYMER_TYPE",
            "POLYMER_LINKING_TYPE",
            "AUTH_SEQ_PART_TYPE",
            "ENTITY_SYNONYMS",
            "ENTITY_DESCRIPTION",
            "ENTITY_ENZYME_CLASS",
            "ENTITY_FRAGMENT_DETAILS",
            "SOURCE_ORGANISM_ORIG",
            "SOURCE_STRAIN_ORIG",
            "SOURCE_TAXID_ORIG",
            "SOURCE_GENE_NAME_ORIG",
            "SOURCE_VARIANT_ORIG",
            "SOURCE_METHOD_ORIG",
            "SOURCE_COMMON_NAME",
            "SOURCE_COMMON_NAME_ORIG",
            "ENTITY_SYNONYMS_ORIG",
            "ENTITY_DESCRIPTION_ORIG",
            "ENTITY_ENZYME_CLASS_ORIG",
            "ENTITY_FRAGMENT_DETAILS_ORIG",
            "ENTITY_MUTATION_DETAILS_ORIG",
            "ENTITY_MUTATION_DETAILS",
            "ENTITY_DETAILS_ORIG",
            "ENTITY_DETAILS",
            "AUTH_SEQ_PART_TYPE_ORIG",
            "HOST_ORG_SOURCE",
            "HOST_ORG_STRAIN",
            "HOST_ORG_TAXID",
            "HOST_ORG_VECTOR",
            "HOST_ORG_VECTOR_TYPE",
            "HOST_ORG_PLASMID",
            "HOST_ORG_COMMON_NAME",
            "HOST_ORG_CELL_LINE",
            "HOST_ORG_SOURCE_ORIG",
            "HOST_ORG_STRAIN_ORIG",
            "HOST_ORG_TAXID_ORIG",
            "HOST_ORG_VECTOR_ORIG",
            "HOST_ORG_VECTOR_TYPE_ORIG",
            "HOST_ORG_VARIANT",
            "HOST_ORG_VARIANT_ORIG",
            "HOST_ORG_PLASMID_ORIG",
            "HOST_ORG_COMMON_NAME_ORIG",
            "HOST_ORG_CELL_LINE_ORIG",
            "DB_MOLECULE_NAME",
            "DB_MOLECULE_SYNONYMS",
            "DB_GENE_NAME",
            "DB_MOLECULE_EC",
            "DB_MOLECULE_DESCRIPTION",
            "DB_MOLECULE_COMMENTS",
            "DB_MOLECULE_KEYWORDS",
            "REF_SEQ_FRAGMENT_DETAILS",
            "REF_ENTRY_ID",
            "REF_ENTRY_ENTITY_ID",
            "REF_ENTRY_STATUS",
            "REF_ENTRY_ANN",
            "AUTH_XYZ_SEQ_BEGIN",
            "AUTH_XYZ_SEQ_END",
            "CURRENT_AUTH_SELECT_ID",
            "CURRENT_REF_SELECT_ID",
        ]
        #
        self.__fListInt = [
            "FULL_LENGTH",
            "ALIGN_LENGTH",
            "MATCH_LENGTH",
            "REF_MATCH_BEGIN",
            "REF_MATCH_END",
            "ORG_ORDER_ID",
            "AUTH_SEQ_NUM_BEGIN",
            "AUTH_SEQ_NUM_END",
            "AUTH_SEQ_PART_ID",
            "AUTH_SEQ_NUM_BEGIN_ORIG",
            "AUTH_SEQ_NUM_END_ORIG",
            "AUTH_SEQ_PART_ID",
            "REF_SORT_ORDER_INDEX",
            "AUTH_SEQ_PART_ID_ORIG",
        ]

        self.__fListFloat = ["AUTH_XYZ_SEQ_SIM", "AUTH_XYZ_SEQ_SIM_WITH_GAPS", "AUTH_REF_SEQ_SIM", "AUTH_REF_SEQ_SIM_WITH_GAPS", "AUTH_REF_SEQ_SIM_BLAST", "REF_SORT_METRIC"]
        self.__fListBool = ["HAS_MANUAL_EDIT", "IS_AUTH_PROVIDED_ID"]

        # self.__fIndex={k:k for k in  self.__fListStr + self.__fListInt + self.__fListBool}
        self.__fD = {}
        self.__reset()

    def __reset(self):
        self.__fD = {}
        for ff in self.__fListStr:
            self.__fD[ff] = ""
        for ff in self.__fListInt:
            self.__fD[ff] = 0
        for ff in self.__fListFloat:
            self.__fD[ff] = 0.0
        for ff in self.__fListBool:
            self.__fD[ff] = False

    def printIt(self, log=sys.stderr):
        log.write("+SequenceFeature() - Feature Dictionary\n")
        for k in sorted(self.__fD.keys()):
            v = self.__fD[k]
            log.write("  Key=%-40s  : %r\n" % (k, v))

    def clear(self):
        self.__reset()

    def get(self):
        return self.__fD

    def set(self, featureDictionary, resetAll=True):
        if resetAll:
            self.__reset()
        for k, v in featureDictionary.items():
            if k in self.__fListStr:
                self.__fD[k] = str(v)
            elif k in self.__fListInt:
                self.__fD[k] = int(str(v))
            elif k in self.__fListFloat:
                self.__fD[k] = float(str(v))
            elif k in self.__fListBool:
                try:
                    self.__fD[k] = bool(v)
                except:  # noqa: E722 pylint: disable=bare-except
                    self.__fD[k] = False

    def setItem(self, ky, val):
        if ky in self.__fListStr:
            self.__fD[ky] = str(val)
        elif ky in self.__fListInt:
            self.__fD[ky] = int(str(val))
        elif ky in self.__fListFloat:
            self.__fD[ky] = float(str(val))
        elif ky in self.__fListBool:
            try:
                self.__fD[ky] = bool(val)
            except:  # noqa: E722 pylint: disable=bare-except
                self.__fD[ky] = False

    def getItem(self, ky):
        if ky in self.__fD:
            return self.__fD[ky]
        else:
            return ""

    def setPolymerType(self, type):  # pylint: disable=redefined-builtin
        self.__fD["POLYMER_TYPE"] = type

    def getPolymerType(self):
        """One of - AA, NA, RNA, DNA, SAC, XNA (i.e.hybrid)"""
        if "POLYMER_TYPE" in self.__fD:
            return self.__fD["POLYMER_TYPE"]
        else:
            return ""

    def setPolymerLinkingType(self, linkingType):
        self.__fD["POLYMER_LINKING_TYPE"] = linkingType

    def getPolymerLinkingType(self):
        """One of the Pdbx full polymer linking types ..."""
        if "POLYMER_LINKING_TYPE" in self.__fD:
            return self.__fD["POLYMER_LINKING_TYPE"]
        else:
            return ""

    def setId(self, dbName="", dbCode="", dbAccession="", dbIsoform=""):
        self.__fD["DB_NAME"] = str(dbName)
        self.__fD["DB_CODE"] = str(dbCode)
        self.__fD["DB_ACCESSION"] = str(dbAccession)
        self.__fD["DB_ISOFORM"] = str(dbIsoform)

    #
    def setDbIsoformDescription(self, description=""):
        self.__fD["DB_ISOFORM_DESCRIPTION"] = str(description)

    #

    def getDbIsoformDescription(self):
        return self.__fD["DB_ISOFORM_DESCRIPTION"]

    def setEntitySynonyms(self, synonyms=""):
        self.__fD["ENTITY_SYNONYMS"] = str(synonyms)

    def setEntityEnzymeClass(self, ec):
        self.__fD["ENTITY_ENZYME_CLASS"] = str(ec)

    def setEntityDescription(self, description=""):
        self.__fD["ENTITY_DESCRIPTION"] = str(description)

    #

    def setEntityDescriptionOrig(self, description=""):
        self.__fD["ENTITY_DESCRIPTION_ORIG"] = str(description)

    def setEntitySynonymsOrig(self, synonyms=""):
        self.__fD["ENTITY_SYNONYMS_ORIG"] = str(synonyms)

    def setEntityEnzymeClassOrig(self, ec):
        self.__fD["ENTITY_ENZYME_CLASS_ORIG"] = str(ec)

    def setEntityFragmentDetails(self, details):
        self.__fD["ENTITY_FRAGMENT_DETAILS"] = str(details)

    def setEntityFragmentDetailsOrig(self, details):
        self.__fD["ENTITY_FRAGMENT_DETAILS_ORIG"] = str(details)

    def setEntityMutationDetails(self, details):
        self.__fD["ENTITY_MUTATION_DETAILS"] = str(details)

    def setEntityMutationDetailsOrig(self, details):
        self.__fD["ENTITY_MUTATION_DETAILS_ORIG"] = str(details)

    def setEntityDetails(self, details):
        self.__fD["ENTITY_DETAILS"] = str(details)

    def setEntityDetailsOrig(self, details):
        self.__fD["ENTITY_DETAILS_ORIG"] = str(details)

    #
    def getEntityDescription(self):
        return self.__fD["ENTITY_DESCRIPTION"]

    def getEntitySynonyms(self):
        return self.__fD["ENTITY_SYNONYMS"]

    def getEntityEnzymeClass(self):
        return self.__fD["ENTITY_ENZYME_CLASS"]

    #
    def getEntityDescriptionOrig(self):
        return self.__fD["ENTITY_DESCRIPTION_ORIG"]

    def getEntitySynonymsOrig(self):
        return self.__fD["ENTITY_SYNONYMS_ORIG"]

    def getEntityEnzymeClassOrig(self):
        return self.__fD["ENTITY_ENZYME_CLASS_ORIG"]

    def getEntityFragmentDetails(self):
        return self.__fD["ENTITY_FRAGMENT_DETAILS"]

    def getEntityFragmentDetailsOrig(self):
        return self.__fD["ENTITY_FRAGMENT_DETAILS_ORIG"]

    def getEntityMutationDetails(self):
        return self.__fD["ENTITY_MUTATION_DETAILS"]

    def getEntityMutationDetailsOrig(self):
        return self.__fD["ENTITY_MUTATION_DETAILS_ORIG"]

    def getEntityDetails(self):
        return self.__fD["ENTITY_DETAILS"]

    def getEntityDetailsOrig(self):
        return self.__fD["ENTITY_DETAILS_ORIG"]

    def getEntitySourceMethodOrig(self):
        return self.__fD["SOURCE_METHOD_ORIG"]

    def getEntitySourceMethod(self):
        return self.__fD["SOURCE_METHOD"]

    def getSourceGeneName(self):
        return self.__fD["SOURCE_GENE_NAME"]

    def getSourceGeneNameOrig(self):
        return self.__fD["SOURCE_GENE_NAME_ORIG"]

    #
    def setTaxId(self, taxid=""):
        self.__fD["SOURCE_TAXID"] = str(taxid)

    def setSource(self, organism="", strain="", taxid="", gene="", method="MAN", commonName="", variant=""):
        self.__fD["SOURCE_ORGANISM"] = str(organism)
        self.__fD["SOURCE_STRAIN"] = str(strain)
        self.__fD["SOURCE_TAXID"] = str(taxid)
        self.__fD["SOURCE_GENE_NAME"] = str(gene)
        self.__fD["SOURCE_METHOD"] = str(method)
        self.__fD["SOURCE_COMMON_NAME"] = str(commonName)
        self.__fD["SOURCE_VARIANT"] = str(variant)

    def setSourceOrig(self, organism="", strain="", taxid="", gene="", method="MAN", commonName="", variant=""):
        self.__fD["SOURCE_ORGANISM_ORIG"] = str(organism)
        self.__fD["SOURCE_STRAIN_ORIG"] = str(strain)
        self.__fD["SOURCE_TAXID_ORIG"] = str(taxid)
        self.__fD["SOURCE_GENE_NAME_ORIG"] = str(gene)
        self.__fD["SOURCE_METHOD_ORIG"] = str(method)
        self.__fD["SOURCE_COMMON_NAME_ORIG"] = str(commonName)
        self.__fD["SOURCE_VARIANT_ORIG"] = str(variant)

    def setHostOrgDetails(self, source="", strain="", taxid="", vector="", vectorType="", plasmid="", commonName="", cellLine="", variant=""):
        self.__fD["HOST_ORG_SOURCE"] = str(source)
        self.__fD["HOST_ORG_STRAIN"] = str(strain)
        self.__fD["HOST_ORG_TAXID"] = str(taxid)
        self.__fD["HOST_ORG_VECTOR"] = str(vector)
        self.__fD["HOST_ORG_VECTOR_TYPE"] = str(vectorType)
        self.__fD["HOST_ORG_PLASMID"] = str(plasmid)
        self.__fD["HOST_ORG_COMMON_NAME"] = str(commonName)
        self.__fD["HOST_ORG_CELL_LINE"] = str(cellLine)
        self.__fD["HOST_ORG_VARIANT"] = str(variant)

    def setHostOrgDetailsOrig(self, source="", strain="", taxid="", vector="", vectorType="", plasmid="", commonName="", cellLine="", variant=""):
        self.__fD["HOST_ORG_SOURCE_ORIG"] = str(source)
        self.__fD["HOST_ORG_STRAIN_ORIG"] = str(strain)
        self.__fD["HOST_ORG_TAXID_ORIG"] = str(taxid)
        self.__fD["HOST_ORG_VECTOR_ORIG"] = str(vector)
        self.__fD["HOST_ORG_VECTOR_TYPE_ORIG"] = str(vectorType)
        self.__fD["HOST_ORG_PLASMID_ORIG"] = str(plasmid)
        self.__fD["HOST_ORG_COMMON_NAME_ORIG"] = str(commonName)
        self.__fD["HOST_ORG_CELL_LINE_ORIG"] = str(cellLine)
        self.__fD["HOST_ORG_VARIANT_ORIG"] = str(variant)

    def setAuthPartDetails(self, partId, seqNumBegin, seqNumEnd, seqPartType=""):
        self.__fD["AUTH_SEQ_NUM_BEGIN"] = int(seqNumBegin)
        self.__fD["AUTH_SEQ_NUM_END"] = int(seqNumEnd)
        self.__fD["AUTH_SEQ_PART_ID"] = int(partId)
        self.__fD["AUTH_SEQ_PART_TYPE"] = str(seqPartType)

    def setAuthPartDetailsOrig(self, partId, seqNumBegin, seqNumEnd, seqPartType=""):
        self.__fD["AUTH_SEQ_NUM_BEGIN_ORIG"] = int(seqNumBegin)
        self.__fD["AUTH_SEQ_NUM_END_ORIG"] = int(seqNumEnd)
        self.__fD["AUTH_SEQ_PART_ID_ORIG"] = int(partId)
        self.__fD["AUTH_SEQ_PART_TYPE_ORIG"] = str(seqPartType)

    def setRefSeqNames(self, proteinName="", synonyms="", geneName=""):
        self.__fD["DB_MOLECULE_NAME"] = str(proteinName)
        self.__fD["DB_MOLECULE_SYNONYMS"] = str(synonyms)
        self.__fD["DB_GENE_NAME"] = str(geneName)
        if self.__debug:
            self.__lfh.write(
                "+SequenceFeature.setRefSeqName() accession %r name %r synonyms %r gene %r\n"
                % (self.__fD["DB_ACCESSION"], self.__fD["DB_MOLECULE_NAME"], self.__fD["DB_MOLECULE_SYNONYMS"], self.__fD["DB_GENE_NAME"])
            )

    def getRefSeqNames(self):
        try:
            return self.__fD["DB_MOLECULE_NAME"], self.__fD["DB_MOLECULE_SYNONYMS"], self.__fD["DB_GENE_NAME"]
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+SequenceFeature.getRefSeqName() failed\n")
        return "", "", ""

    #

    def setRefSeqDetails(self, enzymeClass="", description="", comments="", keywords=""):
        self.__fD["DB_MOLECULE_EC"] = str(enzymeClass)
        self.__fD["DB_MOLECULE_DESCRIPTION"] = str(description)
        self.__fD["DB_MOLECULE_COMMENTS"] = str(comments)
        self.__fD["DB_MOLECULE_KEYWORDS"] = str(keywords)
        if self.__debug:
            self.__lfh.write(
                "+SequenceFeature.setRefSeqDetails() accession %r EC %r descriptions %r comments %r keywords %r\n"
                % (
                    self.__fD["DB_ACCESSION"],
                    self.__fD["DB_MOLECULE_EC"],
                    self.__fD["DB_MOLECULE_DESCRIPTION"],
                    self.__fD["DB_MOLECULE_COMMENTS"],
                    self.__fD["DB_MOLECULE_KEYWORDS"],
                )
            )

    def getRefSeqDetails(self):
        try:
            return self.__fD["DB_MOLECULE_EC"], self.__fD["DB_MOLECULE_DESCRIPTION"], self.__fD["DB_MOLECULE_COMMENTS"], self.__fD["DB_MOLECULE_KEYWORDS"]
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+SequenceFeature.getRefSeqDetails() failed\n")
        return "", "", "", ""

    def setRefSeqVariant(self, variant=""):
        if variant:
            self.__fD["REF_SEQ_VARIANT"] = str(variant)
        #

    def getRefSeqVariant(self):
        if "REF_SEQ_VARIANT" in self.__fD:
            return self.__fD["REF_SEQ_VARIANT"]
        #
        return ""

    def setRefSeqFragmentDetails(self, fragmentDetails=""):
        if fragmentDetails:
            self.__fD["REF_SEQ_FRAGMENT_DETAILS"] = str(fragmentDetails)
        #

    def getRefSeqDepositorInfo(self):
        return {}

    def setRefSortOrder(self, sortIndex, sortMetric):
        self.__fD["REF_SORT_ORDER_INDEX"] = int(float(sortIndex))
        self.__fD["REF_SORT_METRIC"] = float(sortMetric)

    def getRefSortOrderIndex(self):
        try:
            return self.__fD["REF_SORT_ORDER_INDEX"]
        except:  # noqa: E722 pylint: disable=bare-except
            return 0

    def getRefSortMetric(self):
        try:
            return self.__fD["REF_SORT_METRIC"]
        except:  # noqa: E722 pylint: disable=bare-except
            return 0.0

    def getPartInfo(self):
        try:
            return (self.__fD["AUTH_SEQ_PART_ID"], self.__fD["AUTH_SEQ_PART_TYPE"])
        except:  # noqa: E722 pylint: disable=bare-except
            return (None, None)

    def getAuthPartDetails(self):
        try:
            return (self.__fD["AUTH_SEQ_PART_ID"], self.__fD["AUTH_SEQ_NUM_BEGIN"], self.__fD["AUTH_SEQ_NUM_END"], self.__fD["AUTH_SEQ_PART_TYPE"])
        except:  # noqa: E722 pylint: disable=bare-except
            return (None, None, None, None)

    def getAuthPartDetailsOrig(self):
        try:
            return (self.__fD["AUTH_SEQ_PART_ID_ORIG"], self.__fD["AUTH_SEQ_NUM_BEGIN_ORIG"], self.__fD["AUTH_SEQ_NUM_END_ORIG"], self.__fD["AUTH_SEQ_PART_TYPE_ORIG"])
        except:  # noqa: E722 pylint: disable=bare-except
            return (None, None, None, None)

    def getRefDatabaseName(self):
        if "DB_NAME" in self.__fD:
            return self.__fD["DB_NAME"]
        else:
            return ""

    def getRefIsoform(self):
        if "DB_ISOFORM" in self.__fD:
            return self.__fD["DB_ISOFORM"]
        else:
            return ""

    def getSourceOrganism(self):
        if "SOURCE_ORGANISM" in self.__fD:
            return self.__fD["SOURCE_ORGANISM"]
        else:
            return ""

    def getSourceTaxId(self):
        if "SOURCE_TAXID" in self.__fD:
            return self.__fD["SOURCE_TAXID"]
        else:
            return ""

    def getSourceStrain(self):
        if "SOURCE_STRAIN" in self.__fD:
            return self.__fD["SOURCE_STRAIN"]
        else:
            return ""

    def getSourceCommonName(self):
        if "SOURCE_COMMON_NAME" in self.__fD:
            return self.__fD["SOURCE_COMMON_NAME"]
        else:
            return ""

    def getSourceVariant(self):
        if "SOURCE_VARIANT" in self.__fD:
            return self.__fD["SOURCE_VARIANT"]
        else:
            return ""

    def getSourceOrganismOrig(self):
        if "SOURCE_ORGANISM_ORIG" in self.__fD:
            return self.__fD["SOURCE_ORGANISM_ORIG"]
        else:
            return ""

    def getSourceTaxIdOrig(self):
        if "SOURCE_TAXID_ORIG" in self.__fD:
            return self.__fD["SOURCE_TAXID_ORIG"]
        else:
            return ""

    def getSourceStrainOrig(self):
        if "SOURCE_STRAIN_ORIG" in self.__fD:
            return self.__fD["SOURCE_STRAIN_ORIG"]
        else:
            return ""

    def getSourceCommonNameOrig(self):
        if "SOURCE_COMMON_NAME_ORIG" in self.__fD:
            return self.__fD["SOURCE_COMMON_NAME_ORIG"]
        else:
            return ""

    def getSourceVariantOrig(self):
        if "SOURCE_VARIANT_ORIG" in self.__fD:
            return self.__fD["SOURCE_VARIANT_ORIG"]
        else:
            return ""

    def getHostOrgSourceOrganism(self):
        try:
            return self.__fD["HOST_ORG_SOURCE"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getHostOrgStrain(self):
        try:
            return self.__fD["HOST_ORG_STRAIN"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getHostOrgTaxId(self):
        try:
            return self.__fD["HOST_ORG_TAXID"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getHostOrgVector(self):
        try:
            return self.__fD["HOST_ORG_VECTOR"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getHostOrgVectorType(self):
        try:
            return self.__fD["HOST_ORG_VECTOR_TYPE"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getHostOrgPlasmid(self):
        try:
            return self.__fD["HOST_ORG_PLASMID"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getHostOrgCommonName(self):
        try:
            return self.__fD["HOST_ORG_COMMON_NAME"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getHostOrgCellLine(self):
        try:
            return self.__fD["HOST_ORG_CELL_LINE"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getHostOrgVariant(self):
        try:
            return self.__fD["HOST_ORG_VARIANT"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getHostOrgSourceOrganismOrig(self):
        try:
            return self.__fD["HOST_ORG_SOURCE_ORIG"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getHostOrgStrainOrig(self):
        try:
            return self.__fD["HOST_ORG_STRAIN_ORIG"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getHostOrgTaxIdOrig(self):
        try:
            return self.__fD["HOST_ORG_TAXID_ORIG"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getHostOrgVectorOrig(self):
        try:
            return self.__fD["HOST_ORG_VECTOR_ORIG"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getHostOrgVectorTypeOrig(self):
        try:
            return self.__fD["HOST_ORG_VECTOR_TYPE_ORIG"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getHostOrgPlasmidOrig(self):
        try:
            return self.__fD["HOST_ORG_PLASMID_ORIG"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getHostOrgCommonNameOrig(self):
        try:
            return self.__fD["HOST_ORG_COMMON_NAME_ORIG"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getHostOrgCellLineOrig(self):
        try:
            return self.__fD["HOST_ORG_CELL_LINE_ORIG"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getHostOrgVariantOrig(self):
        try:
            return self.__fD["HOST_ORG_VARIANT_ORIG"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    #
    def setManualEditStatus(self, status):
        if status:
            self.__fD["HAS_MANUAL_EDIT"] = True
        else:
            self.__fD["HAS_MANUAL_EDIT"] = False

    def getManualEditStatus(self):
        try:
            return self.__fD["HAS_MANUAL_EDIT"]
        except:  # noqa: E722 pylint: disable=bare-except
            return False

    def setCurentRefSelectId(self, sId):
        self.__fD["CURRENT_REF_SELECT_ID"] = str(sId)

    def getCurrentRefSelectId(self):
        try:
            return self.__fD["CURRENT_REF_SELECT_ID"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    #

    def setCurentAuthSelectId(self, sId):
        self.__fD["CURRENT_AUTH_SELECT_ID"] = str(sId)

    def getCurrentAuthSelectId(self):
        try:
            return self.__fD["CURRENT_AUTH_SELECT_ID"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    #

    def getMatchLength(self):
        if "MATCH_LENGTH" in self.__fD:
            return self.__fD["MATCH_LENGTH"]
        else:
            return 0

    def getAuthRefSimWithGaps(self):
        if "AUTH_REF_SEQ_SIM_WITH_GAPS" in self.__fD:
            return self.__fD["AUTH_REF_SEQ_SIM_WITH_GAPS"]
        else:
            return 0.0

    def setAuthXyzAlignDetails(self, seqLen=0, alignLen=0, seqSim=0.0, seqSimWithGaps=0.0):
        self.__fD["MATCH_LENGTH"] = int(seqLen)
        self.__fD["ALIGN_LENGTH"] = int(alignLen)
        self.__fD["AUTH_XYZ_SEQ_SIM"] = float(seqSim)
        self.__fD["AUTH_XYZ_SEQ_SIM_WITH_GAPS"] = float(seqSimWithGaps)

    def setAuthXyzAlignRange(self, seqBegin="", seqEnd=""):
        """These can potentially include trailing insertion codes so they are treated as strings."""
        self.__fD["AUTH_XYZ_SEQ_BEGIN"] = str(seqBegin)
        self.__fD["AUTH_XYZ_SEQ_END"] = str(seqEnd)

    def getAuthXyzAlignRange(self):
        """These can potentially include trailing insertion codes so they are treated as strings."""
        return self.__fD["AUTH_XYZ_SEQ_BEGIN"], self.__fD["AUTH_XYZ_SEQ_END"]

    def setAuthRefAlignDetails(self, seqLen=0, alignLen=0, seqSim=0.0, seqSimWithGaps=0.0):
        self.__fD["MATCH_LENGTH"] = int(str(seqLen))
        self.__fD["ALIGN_LENGTH"] = int(str(alignLen))
        self.__fD["AUTH_REF_SEQ_SIM"] = float(str(seqSim))
        self.__fD["AUTH_REF_SEQ_SIM_WITH_GAPS"] = float(str(seqSimWithGaps))

    def setAuthRefAlignRange(self, refMatchBegin=0, refMatchEnd=0):
        self.__fD["REF_MATCH_BEGIN"] = refMatchBegin
        self.__fD["REF_MATCH_END"] = refMatchEnd

    def clearAlignDetails(self):
        self.__fD["REF_MATCH_BEGIN"] = 0
        self.__fD["REF_MATCH_END"] = 0
        self.__fD["MATCH_LENGTH"] = 0
        self.__fD["ALIGN_LENGTH"] = 0
        self.__fD["AUTH_XYZ_SEQ_SIM"] = 0.0
        self.__fD["AUTH_XYZ_SEQ_SIM_WITH_GAPS"] = 0.0
        self.__fD["AUTH_REF_SEQ_SIM"] = 0.0
        self.__fD["AUTH_REF_SEQ_SIM_WITH_GAPS"] = 0.0

    def decodeUniProtSourceName(self):
        """Return source and strain -"""
        if "SOURCE_ORGANISM" in self.__fD:
            return self.__getStrain(sourceName=self.__fD["SOURCE_ORGANISM"])
        else:
            return "", ""

    def decodeUniProtSourceOrganism(self, sourceName):
        """Return source and strain -"""
        try:
            return self.__getStrain(sourceName=sourceName)
        except:  # noqa: E722 pylint: disable=bare-except
            pass
        return "", ""

    def __getStrain(self, sourceName=""):
        myRegex = r"\((.*)\)"
        try:
            # self.__lfh.write("+SequenceFeature.decodeUniProtSourceName() decoding %r\n" % sourceName)
            m = re.findall(myRegex, sourceName)
            if len(m) > 0:
                # newSourceName = re.sub(myRegex, "", sourceName)
                if str(m[0]).startswith("strain"):
                    newStrain = str(m[0])[6:]
                else:
                    newStrain = str(m[0])

                if (len(newStrain) < 1) or (newStrain.upper() == "NONE"):
                    newStrain = ""
                # return str(newSourceName).strip(),str(newStrain).strip()
                return str(sourceName).strip(), str(newStrain).strip()
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+SequenceFeature.decodeUniProtSourceName() failed \n")
                traceback.print_exc(file=self.__lfh)

        return sourceName, ""


class ResidueLabel(object):

    """
    Container for residue labels within sequence alignments -

    The residue label includes the following information:

    Sequence type (e.g. ref,auth,coordinate(xyz),)
    Sequence instance identifier    (e.g. PDB chain identifier, PDBx asym_id,)
    Sequence version identifier     (e.g. integer 1, 2, ...  where 1 is the starting version)
    Sequence alternative identifier (e.g. integer 1, 2, ...  where 1 is the primary case)
    Residue code - (e.g. amino acid 3-letter code)
    Residue label index  - (PDB format residue index)
    Sequence index  - residue position in sequence zero-based index
    Alignment index - residue position in alignment (zero-based index)
    Type - AA, RNA, DNA, SA
    Sequence polymer part   (e.g. integer 1,2... where 1 is the starting version)

    The packed residue label has the format:

    seqType + '_' + seqInstId + '_' + seqAltId + '_' + seqVersion + '_' + residueCode3 + '_'
            + residueLabelIndex + '_' +  seqIndex + '_' + alignIndex + '_' + residueType + '_' seqPartId

    """

    def __init__(self, verbose=False):  # pylint: disable=unused-argument
        # self.__versbose = verbose
        self.__reset()

    def __reset(self):
        self.__seqType = ""
        self.__seqInstId = ""
        self.__seqAltId = 1
        self.__seqVersion = 1
        self.__residueCode3 = ""
        self.__residueLabelIndex = ""
        self.__alignIndex = 0
        self.__seqIndex = 0
        self.__residueType = "AA"
        self.__seqPartId = 1

    def set(self, seqType="ref", seqInstId="", seqAltId=1, seqVersion=1, residueCode3="", residueLabelIndex=0, alignIndex=0, seqIndex=0, residueType="AA", seqPartId=1):
        self.__seqType = seqType
        self.__seqInstId = seqInstId
        self.__seqAltId = int(seqAltId)
        self.__seqVersion = int(seqVersion)
        self.__residueCode3 = residueCode3
        self.__residueLabelIndex = residueLabelIndex
        self.__seqIndex = seqIndex
        self.__alignIndex = int(alignIndex)
        self.__residueType = residueType
        self.__seqPartId = seqPartId

    def getSequenceType(self):
        return self.__seqType

    def setSequenceType(self, seqType):
        self.__seqType = seqType

    def getResidueType(self):
        return self.__residueType

    def setResidueType(self, residueType):
        self.__residueType = residueType

    def getSequenceInstId(self):
        return self.__seqInstId

    def setSequenceInstId(self, seqInstId):
        self.__seqInstId = seqInstId

    def getSequenceAltId(self):
        return int(self.__seqAltId)

    def setSequenceAltId(self, seqAltId):
        self.__seqAltId = int(seqAltId)

    def getResidueCode3(self):
        return self.__residueCode3

    def setResidueCode3(self, residueCode3):
        self.__residueCode3 = residueCode3

    def getResidueLabelIndex(self):
        return self.__residueLabelIndex

    def setResidueLabelIndex(self, residueLabelIndex):
        self.__residueLabelIndex = residueLabelIndex

    def getAlignmentIndex(self):
        return int(self.__alignIndex)

    def setAlignmentIndex(self, alignIndex):
        self.__alignIndex = int(alignIndex)

    def getSequenceIndex(self):
        return self.__seqIndex

    def setSequenceIndex(self, seqIndex):
        self.__seqIndex = seqIndex

    def getSequenceVersion(self):
        return int(self.__seqVersion)

    def setSequenceVersion(self, seqVersion):
        self.__seqVersion = int(seqVersion)

    def getSequencePartId(self):
        return int(self.__seqPartId)

    def setSequencePartId(self, seqPartId):
        self.__seqPartId = int(seqPartId)

    def printIt(self, ofh):
        ofh.write("\nResidue Label Contents:\n")
        ofh.write("  Sequence type         %s\n" % self.__seqType)
        ofh.write("  Sequence instance Id  %s\n" % self.__seqInstId)
        ofh.write("  Sequence alternative  %s\n" % str(self.__seqAltId))
        ofh.write("  Sequence version  Id  %d\n" % str(self.__seqVersion))
        ofh.write("  Residue 3-letter code %s\n" % self.__residueCode3)
        ofh.write("  Residue label index   %s\n" % self.__residueLabelIndex)
        ofh.write("  Sequence  index       %6s\n" % self.__seqIndex)
        ofh.write("  Alignment index       %6d\n" % int(self.__alignIndex))
        ofh.write("  Residue type          %6s\n" % str(self.__residueType))
        ofh.write("  Sequence part Id      %6s\n" % str(self.__seqPartId))

    def pack(self):
        """
        Return a string identifier containing residue label details.

        """
        tCode = self.__residueCode3
        if tCode == ".":
            tCode = "-"
        idS = (
            self.__seqType
            + "_"
            + self.__seqInstId
            + "_"
            + str(self.__seqAltId)
            + "_"
            + str(self.__seqVersion)
            + "_"
            + tCode
            + "_"
            + str(self.__residueLabelIndex)
            + "_"
            + str(self.__seqIndex)
            + "_"
            + str(self.__alignIndex)
            + "_"
            + str(self.__residueType)
            + "_"
            + str(self.__seqPartId)
        )

        return idS

    def unpack(self, residueLabelString):
        """
        Update the internal residue details with the data stored in the input residue label.
        """
        self.__reset()
        try:
            idL = residueLabelString.split("_")
            self.__seqType = idL[0]
            self.__seqInstId = idL[1]
            self.__seqAltId = idL[2]
            self.__seqVersion = idL[3]
            self.__residueCode3 = idL[4]
            if self.__residueCode3 == "-":
                self.__residueCode3 = "."
            self.__residueLabelIndex = idL[5]
            self.__seqIndex = idL[6]
            self.__alignIndex = idL[7]
            self.__residueType = idL[8]
            self.__seqPartId = idL[9]

            return True
        except:  # noqa: E722 pylint: disable=bare-except
            return False

    def getSeq(self):
        return (self.__seqType, self.__seqInstId, self.__seqPartId, self.__seqAltId, self.__seqVersion)

    def getDict(self):
        """
        Return a dictionary of the elements of the residue label details.
        """
        idD = {}
        idD["TYPE"] = self.__seqType
        idD["CHAIN_ID"] = self.__seqInstId
        idD["ALTERNATIVE"] = self.__seqAltId
        idD["VERSION"] = self.__seqVersion
        idD["CODE3"] = self.__residueCode3
        idD["INDEX_LABEL"] = self.__residueLabelIndex
        idD["INDEX_ALIGNMENT"] = self.__alignIndex
        idD["INDEX_SEQUENCE"] = self.__seqIndex
        idD["RESIDUE_TYPE"] = self.__residueType
        idD["PART_ID"] = self.__residueType

        return idD


class SequenceLabel(object):

    """
    Container for sequence labels within sequence alignment depictions -

    The sequence label includes the following information:

    Sequence type                 (e.g. ref,auth,coordinate(xyz),)
    Sequence instance identifier  (e.g. entity or asym_id identifier, PDBx asym_id,)
    Sequence part identifier      (Identifies parts of multi-source entities 1,2,3 ...  1 is default)
    Alternative identifier        (Identifies alternatives via integer 0, 1, 2, ...  where 0 is primary case)
    Version identifier            (e.g. integer 0, 1, 2, ...  where 0 is the starting version)

    The packed sequence label has the format:

    seqType + '_' + seqInstId + '_' + seqPartId + '_' + seqAltId  + '_' + seqVersion

    """

    def __init__(self, verbose=False):  # pylint: disable=unused-argument
        # self.__versbose = verbose
        self.__reset()
        self.__typeOrder = ["auth", "xyz", "ref"]

    def __reset(self):
        self.seqType = ""
        self.seqInstId = ""
        self.seqPartId = 1
        self.seqAltId = 1
        self.seqVersion = 1

    def __eq__(self, other):
        if isinstance(other, SequenceLabel):
            return self.get() == other.get()
        return NotImplemented

    def __ne__(self, other):
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def __str__(self):
        return "Type: %5s Id: %4s Part: %2d Alt: %4d Version: %3d" % self.get()

    def __repr__(self):
        return "Type: %5s Id: %4s Part: %2d Alt: %4d Version: %3d" % self.get()

    def __lt__(self, other):
        if isinstance(other, SequenceLabel):
            if self.seqType != other.seqType:
                try:
                    return self.__typeOrder.index(self.seqType) < self.__typeOrder.index(other.seqType)
                except:  # noqa: E722 pylint: disable=bare-except
                    # ignore failures
                    return False
            else:
                return (self.seqVersion > other.seqVersion) or ((self.seqInstId < other.seqInstId) and (self.seqPartId < other.seqPartId) and (self.seqAltId > other.seqAltId))
        return NotImplemented

    def set(self, seqType="ref", seqInstId="", seqPartId=1, seqAltId=1, seqVersion=1):
        self.seqType = seqType if seqType in self.__typeOrder else "unknown"
        self.seqInstId = str(seqInstId)
        self.seqPartId = int(seqPartId)
        self.seqAltId = int(seqAltId)
        self.seqVersion = int(seqVersion)

    def get(self):
        return (self.seqType, self.seqInstId, int(self.seqPartId), int(self.seqAltId), int(self.seqVersion))

    def getSequenceType(self):
        return self.seqType

    def setSequenceType(self, seqType):
        self.seqType = seqType
        return seqType in self.__typeOrder

    def getSequenceInstId(self):
        return self.seqInstId

    def setSequenceInstId(self, seqInstId):
        self.seqInstId = str(seqInstId)

    def getSequencePartId(self):
        return self.seqPartId

    def setSequencePartId(self, seqPartId):
        self.seqInstId = int(seqPartId)

    def getSequenceAlternativeId(self):
        return int(self.seqAltId)

    def setSequenceAlternativeId(self, seqAltId):
        self.seqAltId = int(seqAltId)

    def getSequenceVersion(self):
        return int(self.seqVersion)

    def setSequenceVersion(self, seqVersion):
        self.seqVersion = int(seqVersion)

    def printIt(self, ofh):
        ofh.write("\n+SequenceLabel.printIt() Sequence Label Contents:\n")
        ofh.write("  Sequence type            %s\n" % self.seqType)
        ofh.write("  Sequence instance    Id  %s\n" % self.seqInstId)
        ofh.write("  Sequence part        Id  %s\n" % str(self.seqPartId))
        ofh.write("  Sequence alternative Id  %s\n" % str(self.seqAltId))
        ofh.write("  Sequence version     Id  %s\n" % str(self.seqVersion))

    def pack(self):
        """
        Return a string identifier containing sequence label details.

        """
        idS = self.seqType + "_" + self.seqInstId + "_" + str(self.seqPartId) + "_" + str(self.seqAltId) + "_" + str(self.seqVersion)
        return idS

    def unpack(self, sequenceLabelString):
        """
        Update the internal sequence details with the data stored in the input sequence label.
        """
        self.__reset()
        try:
            idL = sequenceLabelString.split("_")
            self.seqType = idL[0]
            self.seqInstId = idL[1]
            self.seqPartId = int(idL[2])
            self.seqAltId = int(idL[3])
            self.seqVersion = int(idL[4])
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            return False

    def getDict(self):
        """
        Return a dictionary of the elements of the sequence label details.
        """
        idD = {}
        idD["TYPE"] = self.seqType
        idD["INSTANCE_ID"] = self.seqInstId
        idD["PART_ID"] = self.seqPartId
        idD["ALTERNATIVE_ID"] = self.seqAltId
        idD["VERSION"] = self.seqVersion
        return idD


class SequenceLabelUtils(object):

    """
    Utilitity functions on sequence labels.
    """

    def __init__(self, verbose=False, log=sys.stderr):  # pylint: disable=unused-argument
        # self.__verbose = verbose
        # self.__lfh = log
        pass

    def getAlignGroupId(self, alignIdList):
        sLabel = SequenceLabel()
        for alignId in alignIdList:
            sLabel.unpack(alignId)
            (seqType, seqInstId, _seqPartId, _seqAltId, _seqVersion) = sLabel.get()
            if seqType in ["auth", "ref"]:
                return seqInstId
