##
# File:    SequenceFeatureDepict.py
# Date:    03-Mar-2013
#
# Updates:
# 10-Nov-2013  -- Add entity detail dictionary markup
# 19-Oct-2022  zf add author provided reference sequence information
#
##
"""
Convenience methods to markup sequence features.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.07"

import sys


class SequenceFeatureDepict(object):
    """Methods to markup sequence feature details."""

    def __init__(self, sfObj=None, verbose=False, log=sys.stderr):  # pylint: disable=unused-argument
        """sfObj is a sequence feature object -"""
        self.__sfObj = sfObj
        # self.__verbose = verbose
        # self.__lfh = log
        #
        self.__fD = self.__sfObj.get()

    def markupDatabaseReferenceWithUrl(self, seqAltId=0):
        return self.__markupDatabaseReferenceWithUrl(
            dbName=self.__fD["DB_NAME"],
            dbAccession=self.__fD["DB_ACCESSION"],
            dbIsoForm=self.__fD["DB_ISOFORM"],
            seqAltId=seqAltId,
            entryId=self.__fD["REF_ENTRY_ID"],
            entityId=self.__fD["REF_ENTRY_ENTITY_ID"],
            statusCode=self.__fD["REF_ENTRY_STATUS"],
            annInitial=self.__fD["REF_ENTRY_ANN"],
            isAuthProvidedId=self.__fD["IS_AUTH_PROVIDED_ID"]
        )

    def __markupDatabaseReferenceWithUrl(
        self, dbName, dbAccession, dbIsoForm="", seqAltId=0, entryId="", entityId="", statusCode="", annInitial="", isAuthProvidedId=False
    ):  # pylint: disable=unused-argument

        style = ""
        span_start = ""
        span_end = ""
        if isAuthProvidedId:
            style = 'style="color:green;"'
            span_start = '<span style="color:green;">'
            span_end = '</span>'
        #
        displayCode = dbAccession
        if len(dbIsoForm) > 0:
            displayCode = dbIsoForm
        #
        if dbName in ["UNP", "SP", "TR"]:
            lab = '<a href="http://www.uniprot.org/uniprot/%s" target="blank"><span %s class="nowrap">%s</span></a>' \
                  % (dbAccession, style, dbName + ":" + displayCode)
#       elif dbName in ["GENBANK", "GB"]:
#           lab = dbName + ":" + displayCode
        else:
            lab = span_start + dbName + ":" + displayCode + span_end
        #
        if entryId:
            lab += "<br/><br/>From:<br/>DepID:" + entryId + "<br/>EntityID:" + entityId + "<br/>Status:" + statusCode + "<br/>Ann:" + annInitial
        #
        return lab

    # def markupAuthorFeatures(self):
    #     return self.__markupAlignmentFeatures(dbSourceOrg=self.__fD["SOURCE_ORGANISM"], taxId=self.__fD["SOURCE_TAXID"])

    def markupReferenceSimilarttFeatures(self):
        return self.__markupSimilarttFeatures(seqSimWithGaps=self.__fD["AUTH_REF_SEQ_SIM_WITH_GAPS"], seqSim=self.__fD["AUTH_REF_SEQ_SIM"])

    def markupReferenceAlignmentFeatures(self):
        return self.__markupAlignmentFeatures(
            refSeqFullLength=self.__fD["FULL_LENGTH"],
            alignLength=self.__fD["ALIGN_LENGTH"],
            seqSim=self.__fD["AUTH_REF_SEQ_SIM"],
            seqSimWithGaps=self.__fD["AUTH_REF_SEQ_SIM_WITH_GAPS"],
            alignBegin=self.__fD["REF_MATCH_BEGIN"],
            alignEnd=self.__fD["REF_MATCH_END"],
        )

    def __markupSimilarttFeatures(self, seqSimWithGaps=0.0, seqSim=0.0):
        dL = []
        #
        if seqSimWithGaps > 0.001:
            tS = '<span class="detailkey">w/ gaps: </span><span class="detailvalue">%6.3f</span><br />' % float(seqSimWithGaps)
            dL.append(tS)

        if seqSim > 0.001:
            tS = '<span class="detailkey">w/o gaps: </span><span class="detailvalue">%6.3f</span>' % float(seqSim)
            dL.append(tS)

        return "\n".join(dL)

    def __markupAlignmentFeatures(
        self,  # pylint: disable=unused-argument
        refSeqFullLength=0,
        alignLength=0,
        seqSim=0.0,  # pylint: disable=unused-argument
        seqSimWithGaps=0.0,  # pylint: disable=unused-argument
        alignBegin=0,
        alignEnd=0,
    ):

        dL = []

        if refSeqFullLength > 0:
            tS = '<span class="detailkey">Full sequence length: </span><span class="detailvalue">%d</span><br />' % int(refSeqFullLength)
            dL.append(tS)

        # """
        # if ( seqSim > 0.001):
        #     tS='<span class="detailkey">Identity (w/o gaps): </span><span class="detailvalue">%6.3f</span><br />' % float(seqSim)
        #     dL.append(tS)

        # if ( seqSimWithGaps > 0.001):
        #     tS='<span class="detailkey">Identity (w/ gaps): </span><span class="detailvalue">%6.3f</span><br />' % float(seqSimWithGaps)
        #     dL.append(tS)
        # """

        if alignLength > 0:
            tS = '<span class="detailkey">Align length (w/ author sequence): </span><span class="detailvalue">%d</span><br />' % int(alignLength)
            dL.append(tS)

        if alignEnd > 0:
            tS = '<span class="detailkey">Align range: </span><span class="detailvalue">%d-%d</span><br />' % (int(alignBegin), int(alignEnd))
            dL.append(tS)

        return "\n".join(dL)

    def markupReferenceFeatures(self):
        return self.__markupReferenceFeatures(
            sourceOrg=self.__fD["SOURCE_ORGANISM"],
            geneName=self.__fD["DB_GENE_NAME"],
            ec=self.__fD["DB_MOLECULE_EC"],
            synonyms=self.__fD["DB_MOLECULE_SYNONYMS"],
            moleculeName=self.__fD["DB_MOLECULE_NAME"],
            orgCommonName=self.__fD["SOURCE_COMMON_NAME"],
            strain=self.__fD["SOURCE_STRAIN"],
            isoformDescription=self.__fD["DB_ISOFORM_DESCRIPTION"],
        )

    def __markupReferenceFeatures(self, sourceOrg="", geneName="", ec="", synonyms="", moleculeName="", orgCommonName="", strain="", isoformDescription=""):
        dL = []
        if len(moleculeName) > 0:
            tS = '<span class="detailkey">Name: </span><span class="detailvalue">%s</span><br />' % moleculeName
            dL.append(tS)

        if len(geneName) > 0:
            tS = '<span class="detailkey">Gene: </span><span class="detailvalue">%s</span><br />' % geneName
            dL.append(tS)

        if len(synonyms) > 0:
            tS = '<span class="detailkey">Synonyms: </span><span class="detailvalue">%s</span><br />' % synonyms
            dL.append(tS)

        if len(sourceOrg) > 0:
            tS = '<span class="detailkey">Organism: </span><span class="detailvalue">%s</span><br />' % sourceOrg
            dL.append(tS)

        if len(orgCommonName) > 0:
            tS = '<span class="detailkey">Org Common: </span><span class="detailvalue">%s</span><br />' % orgCommonName
            dL.append(tS)

        if len(strain) > 0:
            tS = '<span class="detailkey">Strain: </span><span class="detailvalue">%s</span><br />' % strain
            dL.append(tS)

        if len(ec) > 0:
            tS = '<span class="detailkey">EC: </span><span class="detailvalue">%s</span><br />' % ec
            dL.append(tS)

        if len(isoformDescription) > 0:
            tS = '<span class="detailkey">Isoform: </span><span class="detailvalue">%s</span><br />' % isoformDescription
            dL.append(tS)

        return "\n".join(dL)

    def markupXyzFeatures(self):
        return self.__markupXyzFeatures(alignLength=self.__fD["ALIGN_LENGTH"], seqSim=self.__fD["AUTH_XYZ_SEQ_SIM"], seqSimWithGaps=self.__fD["AUTH_XYZ_SEQ_SIM_WITH_GAPS"])

    def __markupXyzFeatures(self, alignLength=0, seqSim=0.0, seqSimWithGaps=0.0):

        dL = []
        if alignLength > 0:
            tS = '<span class="detailkey">Align length (w/ author sequence): </span><span class="detailvalue">%d</span><br />' % int(alignLength)
            dL.append(tS)

        if seqSim > 0.001:
            tS = '<span class="detailkey">Identity (w/o gaps):</span><span class="detailvalue">%6.3f </span><br />' % float(seqSim)
            dL.append(tS)

        if seqSimWithGaps > 0.001:
            tS = '<span class="detailkey">Identity (w/ gaps):</span><span class="detailvalue">%6.3f </span><br />' % float(seqSimWithGaps)
            dL.append(tS)

        return "\n".join(dL)

    def markupCurrentEntityDetails(self):
        """Markup the current assignments for entity/source details from the input feature dictionary"""
        detailsTemplate = "Part %s:<br /> %s <br />(%4s - %4s)"
        spStr = "<br />"
        mD = {}
        for k in ["partdetails", "sourceAndStrain", "description", "hostorg"]:
            mD[k] = ""

        tDescription = self.__sfObj.getEntityDescription()
        tEc = self.__sfObj.getEntityEnzymeClass()
        tFrag = self.__sfObj.getEntityFragmentDetails()
        tMutation = self.__sfObj.getEntityMutationDetails()
        tDetails = self.__sfObj.getEntityDetails()
        tName = self.__sfObj.getEntitySynonyms()
        tMethod = self.__sfObj.getEntitySourceMethod()

        mD["description"] = ""
        if len(tDescription) > 1:
            mD["description"] += "<b>Name:&nbsp;</b> " + tDescription + spStr
        if len(tName) > 1:
            mD["description"] += "<b>Synonyms:&nbsp;</b> " + tName + spStr
        if len(tEc) > 1:
            mD["description"] += "<b>EC:&nbsp;</b> " + tEc + spStr
        if len(tFrag) > 1:
            mD["description"] += "<b>Fragment:&nbsp;</b> " + tFrag + spStr
        if len(tMutation) > 1:
            mD["description"] += "<b>Mutation:&nbsp;</b> " + tMutation + spStr
        if len(tDetails) > 1:
            mD["description"] += "<b>Entity details:&nbsp;</b> " + tDetails + spStr
        if len(tMethod) > 1:
            mD["description"] += "<b>Source method:&nbsp;</b> " + tMethod + spStr

        seqPartId, seqNumBeg, seqNumEnd, seqPartType = self.__sfObj.getAuthPartDetails()
        mD["partdetails"] = detailsTemplate % (seqPartId, seqPartType, seqNumBeg, seqNumEnd)

        authOrg = self.__sfObj.getSourceOrganism()
        authStrain = self.__sfObj.getSourceStrain()
        authCommonOrg = self.__sfObj.getSourceCommonName()
        taxId = self.__sfObj.getSourceTaxId()
        authGeneName = self.__sfObj.getSourceGeneName()
        authVariant = self.__sfObj.getSourceVariant()

        if len(authOrg) > 1:
            mD["sourceAndStrain"] += "<b>Organism:&nbsp;</b> " + authOrg + spStr
        if len(authCommonOrg) > 1:
            mD["sourceAndStrain"] += "<b>Common name:&nbsp;</b> " + authCommonOrg + spStr
        if len(authStrain) > 1:
            mD["sourceAndStrain"] += "<b>Strain:&nbsp;</b> " + authStrain + spStr
        if len(taxId) > 1:
            mD["sourceAndStrain"] += "<b>TaxID:&nbsp;</b> " + taxId + spStr
        if len(authGeneName) > 1:
            mD["sourceAndStrain"] += "<b>Gene:&nbsp;</b> " + authGeneName + spStr
        if len(authVariant) > 1:
            mD["sourceAndStrain"] += "<b>Variant:&nbsp;</b> " + authVariant + spStr

        #
        #  Host organism details
        hostOrgName = self.__sfObj.getHostOrgSourceOrganism()
        if len(hostOrgName) > 1:
            mD["hostorg"] += "<b>Name:&nbsp;</b> " + hostOrgName + spStr
        hostOrgCommonName = self.__sfObj.getHostOrgCommonName()
        if len(hostOrgCommonName) > 1:
            mD["hostorg"] += "<b>Common Name:&nbsp;</b> " + hostOrgCommonName + spStr
        hostOrgStrain = self.__sfObj.getHostOrgStrain()
        if len(hostOrgStrain) > 1:
            mD["hostorg"] += "<b>Strain:&nbsp;</b> " + hostOrgStrain + spStr
        hostOrgTaxId = self.__sfObj.getHostOrgTaxId()
        if len(hostOrgTaxId) > 1:
            mD["hostorg"] += "<b>TaxId:&nbsp;</b> " + hostOrgTaxId + spStr

        hostOrgVector = self.__sfObj.getHostOrgVector()
        if len(hostOrgVector) > 1:
            mD["hostorg"] += "<b>Vector:&nbsp;</b> " + hostOrgVector + spStr
        hostOrgVectorType = self.__sfObj.getHostOrgVectorType()
        if len(hostOrgVectorType) > 1:
            mD["hostorg"] += "<b>Vector type:&nbsp;</b> " + hostOrgVectorType + spStr

        hostOrgPlasmid = self.__sfObj.getHostOrgPlasmid()
        if len(hostOrgPlasmid) > 1:
            mD["hostorg"] += "<b>Plasmid:&nbsp;</b> " + hostOrgPlasmid + spStr
        return mD

    def markupAuthorEntityDetails(self):
        """Markup the author provided entity/source details from the input feature dictionary..."""
        spStr = "<br />"
        detailsTemplate = "Part %s:<br /> %s <br />(%4s - %4s)"
        mD = {}
        for k in ["partdetails", "sourceAndStrain", "description", "hostorg"]:
            mD[k] = ""

        tDescription = self.__sfObj.getEntityDescriptionOrig()
        tEc = self.__sfObj.getEntityEnzymeClassOrig()
        tFrag = self.__sfObj.getEntityFragmentDetailsOrig()
        tMutation = self.__sfObj.getEntityMutationDetailsOrig()
        tDetails = self.__sfObj.getEntityDetailsOrig()
        tName = self.__sfObj.getEntitySynonymsOrig()
        tMethod = self.__sfObj.getEntitySourceMethodOrig()

        mD["description"] = ""
        if len(tDescription) > 1:
            mD["description"] += "<b>Name:&nbsp;</b> " + tDescription + spStr
        if len(tName) > 1:
            mD["description"] += "<b>Synonyms:&nbsp;</b> " + tName + spStr
        if len(tEc) > 1:
            mD["description"] += "<b>EC:&nbsp;</b> " + tEc + spStr
        if len(tFrag) > 1:
            mD["description"] += "<b>Fragment:&nbsp;</b> " + tFrag + spStr
        if len(tMutation) > 1:
            mD["description"] += "<b>Mutation:&nbsp;</b> " + tMutation + spStr
        if len(tDetails) > 1:
            mD["description"] += "<b>Entity details:&nbsp;</b> " + tDetails + spStr
        if len(tMethod) > 1:
            mD["description"] += "<b>Source method:&nbsp;</b> " + tMethod + spStr

        seqPartId, seqNumBeg, seqNumEnd, seqPartType = self.__sfObj.getAuthPartDetailsOrig()
        mD["partdetails"] = detailsTemplate % (seqPartId, seqPartType, seqNumBeg, seqNumEnd)

        authOrg = self.__sfObj.getSourceOrganismOrig()
        authStrain = self.__sfObj.getSourceStrainOrig()
        authCommonOrg = self.__sfObj.getSourceCommonNameOrig()
        taxId = self.__sfObj.getSourceTaxIdOrig()
        authGeneName = self.__sfObj.getSourceGeneNameOrig()
        authVariant = self.__sfObj.getSourceVariantOrig()

        if len(authOrg) > 1:
            mD["sourceAndStrain"] += "<b>Organism:&nbsp;</b> " + authOrg + spStr
        if len(authCommonOrg) > 1:
            mD["sourceAndStrain"] += "<b>Common name:&nbsp;</b> " + authCommonOrg + spStr
        if len(authStrain) > 1:
            mD["sourceAndStrain"] += "<b>Strain:&nbsp;</b> " + authStrain + spStr
        if len(taxId) > 1:
            mD["sourceAndStrain"] += "<b>TaxID:&nbsp;</b> " + taxId + spStr
        if len(authGeneName) > 1:
            mD["sourceAndStrain"] += "<b>Gene:&nbsp;</b> " + authGeneName + spStr
        if len(authVariant) > 1:
            mD["sourceAndStrain"] += "<b>Variant:&nbsp;</b> " + authVariant + spStr

        #
        #  Host organism details
        hostOrgName = self.__sfObj.getHostOrgSourceOrganismOrig()
        if len(hostOrgName) > 1:
            mD["hostorg"] += "<b>Name:&nbsp;</b> " + hostOrgName + spStr
        hostOrgCommonName = self.__sfObj.getHostOrgCommonNameOrig()
        if len(hostOrgCommonName) > 1:
            mD["hostorg"] += "<b>Common Name:&nbsp;</b> " + hostOrgCommonName + spStr
        hostOrgStrain = self.__sfObj.getHostOrgStrainOrig()
        if len(hostOrgStrain) > 1:
            mD["hostorg"] += "<b>Strain:&nbsp;</b> " + hostOrgStrain + spStr
        hostOrgTaxId = self.__sfObj.getHostOrgTaxIdOrig()
        if len(hostOrgTaxId) > 1:
            mD["hostorg"] += "<b>TaxId:&nbsp;</b> " + hostOrgTaxId + spStr
        hostOrgVector = self.__sfObj.getHostOrgVectorOrig()
        if len(hostOrgVector) > 1:
            mD["hostorg"] += "<b>Vector:&nbsp;</b> " + hostOrgVector + spStr
        hostOrgVectorType = self.__sfObj.getHostOrgVectorTypeOrig()
        if len(hostOrgVectorType) > 1:
            mD["hostorg"] += "<b>Vector type:&nbsp;</b> " + hostOrgVectorType + spStr
        hostOrgPlasmid = self.__sfObj.getHostOrgPlasmidOrig()
        if len(hostOrgPlasmid) > 1:
            mD["hostorg"] += "<b>Plasmid:&nbsp;</b> " + hostOrgPlasmid + spStr

        return mD
