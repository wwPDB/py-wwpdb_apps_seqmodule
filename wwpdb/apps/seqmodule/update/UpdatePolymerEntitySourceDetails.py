##
# File:    UpdatePolymerEntitySourceDetails.py
# Date:    22-Mar-2014
#
# Updates:
# 1-Feb-2015  jdw add HOST_ORG_VARIANT
# 07-Sep-2017 zf  modified to allow annotator to manually input Norine reference information for small polymer sequence
##
"""
Utilities for updating polymer entity source attribute details.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.08"

import sys
import traceback

from wwpdb.apps.seqmodule.io.SequenceDataStore import SequenceDataStore
from wwpdb.apps.seqmodule.util.SequenceLabel import SequenceLabel, SequenceFeature, SequenceFeatureMap
from wwpdb.apps.seqmodule.util.SequenceAssign import SequenceAssignArchive, SequenceAssignDepositor


class UpdatePolymerEntitySourceDetails(object):

    """Utilities for adding out-of-band sequence and feature data to the sequence data store."""

    def __init__(self, reqObj=None, verbose=False, log=sys.stderr):
        self.__verbose = verbose
        self.__debug = False
        self.__reqObj = reqObj
        self.__lfh = log
        #
        self.__entityReviewItemList = [
            "ENTITY_DESCRIPTION",
            "ENTITY_SYNONYMS",
            "SOURCE_ORGANISM",
            "SOURCE_COMMON_NAME",
            "SOURCE_TAXID",
            "SOURCE_GENE_NAME",
            "SOURCE_STRAIN",
            "SOURCE_VARIANT",
            "ENTITY_ENZYME_CLASS",
            "ENTITY_FRAGMENT_DETAILS",
            "ENTITY_MUTATION_DETAILS",
            "ENTITY_DETAILS",
            "HOST_ORG_SOURCE",
            "HOST_ORG_COMMON_NAME",
            "HOST_ORG_STRAIN",
            "HOST_ORG_TAXID",
            "HOST_ORG_VECTOR",
            "HOST_ORG_VECTOR_TYPE",
            "HOST_ORG_PLASMID",
            "HOST_ORG_CELL_LINE",
            "HOST_ORG_VARIANT",
            "ANNO_EDIT_DB_NAME",
            "ANNO_EDIT_DB_CODE",
            "ANNO_EDIT_DB_ACCESSION",
            "ANNO_EDIT_DB_ALIGN_BEGIN",
            "ANNO_EDIT_DB_ALIGN_END",
            "SOURCE_METHOD",
            "CURRENT_AUTH_SELECT_ID",
            "CURRENT_REF_SELECT_ID",
        ]

        self.__selectIdList = []
        self.__setup()

    def __setup(self):
        try:
            self.__placeHolderValue = "click-to-edit"
            # self.__siteId = self.__reqObj.getValue("WWPDB_SITE_ID")
            self.__sessionId = self.__reqObj.getSessionId()
            self.__sessionObj = self.__reqObj.getSessionObj()
            # self.__sessionPath = self.__sessionObj.getPath()
            #
            self.__selectIdList = self.__reqObj.getSummarySelectList()
            self.__sds = SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+UpdatePolymerEntitySourceDetails.__setup() sessionId %s failed\n" % (self.__sessionObj.getId()))
            #
        #

    def updateAuthEntityDetails(self, selectIdList=None):
        """Update current auth/sample entity attribute details based on current selected reference sequences -"""
        if selectIdList is not None:
            self.__selectIdList = selectIdList
        #
        self.__lfh.write("+UpdatePolymerEntitySourceDetails.updateAuthEntityDetails() starting with selectIdList %r\n" % selectIdList)
        entityIdList = self.__sds.getGroupIds()
        for entityId in entityIdList:
            instanceIds = self.__sds.getGroup(groupId=entityId)
            if len(instanceIds) == 0:
                continue
            # JDW CHANGE
            # instanceId=instanceIds[0]
            instanceId = entityId
            partIdList = self.__sds.getPartIds(instanceId, dataType="sequence", seqType="auth")
            #
            self.__lfh.write("+UpdatePolymerEntitySourceDetails.updateAuthEntityDetails() entityId %r instanceId %r partIdList %r\n" % (entityId, instanceId, partIdList))
            sfm = SequenceFeatureMap(verbose=self.__verbose, log=self.__lfh)
            for partId in partIdList:
                _authSelectId, authSL, authFdObj = self.__getCurrentAuthSelection(entityId, partId=partId)
                _seqType, seqInstId, _seqPartId, seqAltId, seqVersion = authSL.get()
                #
                refSelectId, _refSL, refFdObj = self.__getCurrentRefSelection(entityId, partId=partId)
                #
                # filter special cases -
                #
                curRefId = authFdObj.getCurrentRefSelectId()
                hasEdit = authFdObj.getManualEditStatus()
                self.__lfh.write(
                    "+UpdatePolymerEntitySourceDetails.updateAuthEntityDetails() entityId %r partId %r refSelectId %r curRefId %r hasEdit %r\n"
                    % (entityId, partId, refSelectId, curRefId, hasEdit)
                )
                #
                if curRefId == refSelectId and hasEdit:
                    self.__lfh.write(
                        "+UpdatePolymerEntitySourceDetails.updateAuthEntityDetails() skip update entityId %r partId %r curRefId %r hasEdit %r\n"
                        % (entityId, partId, curRefId, hasEdit)
                    )
                    #
                    continue
                #
                selfRefTarget = "selfref_%s_%d" % (entityId, partId)
                if refSelectId and refSelectId.startswith(selfRefTarget):
                    # use original auth values -
                    sfm.updateAuthOrig(authFdObj.get())
                elif refFdObj:
                    sfm.updateAuth(authFdObj.get(), refFdObj.get())
                self.__sds.setFeature(authFdObj.get(), seqInstId, "auth", partId=partId, altId=seqAltId, version=seqVersion)
            #
        #
        self.__sds.serialize()
        return True

    def makeEntitySourceDetailsForm(self, entityId, entryId=""):
        """Wrapper for form generation entity detail input. Scope of form is for each entity including all parts --"""
        curPartIdList = self.__sds.getPartIds(entityId, dataType="sequence", seqType="auth")
        #
        rL = []
        rL.append('<div id="sectentityreview">')
        for partId in curPartIdList:
            rL.append(self.__makeEntityPartSourceDetailsForm(entityId, partId, entryId=entryId, numParts=len(curPartIdList)))
        #
        rL.append("</div>")
        #
        rD = {}
        rD["htmlcontent"] = "".join(rL)
        return rD

    def __makeEntityPartSourceDetailsForm(self, entityId, partId=1, entryId="", numParts=1):
        """Return an edit form for to review entity annotations."""
        #
        form_template = """

        <h3>Entity Review Form for Entry %(entryid)s Entity %(entityid)s Part %(partid)s</h3>
        <h4>Entry title: %(STRUCT_TITLE)s</h4>
        <h4>Citation title: %(CITATION_TITLE)s</h4>

        <form name="formentityreview_%(entityid)s_%(partid)s" id="formentityreview_%(entityid)s_%(partid)s"
          action="/service/sequence_editor/respond_form/entityreview" method="post" class="review_ajaxform">
            <input type="submit" name="submit" value="%(savebuttontext)s" class="disableonclick submitparentform fltrgt" />
            <input type="hidden" name="sessionid" value="%(sessionid)s" />
            <input type="hidden" name="partid" value="%(partid)s" />
            <input type="hidden" name="entityid" value="%(entityid)s" />
            <input type="hidden" name="instid" value="%(instid)s" />
            <input type="hidden" name="altid" value="%(altid)s" />
            <input type="hidden" name="versionid" value="%(versionid)s" />

            <input type="hidden" name="CURRENT_AUTH_SELECT_ID" value="%(CURRENT_AUTH_SELECT_ID)s" />
            <input type="hidden" name="CURRENT_REF_SELECT_ID" value="%(CURRENT_REF_SELECT_ID)s" />
            <br class="clearfloat" />
            <table>
            <tr>
               <th>Item </th>
               <th>Current Value</th>
               <th>Reference DB</th>
               <th>Author Provided</th>
            </tr>

             <tr><td>DB Name</td>
                 <td %(bgcolclass)s>%(db_name)s</td>
                 <td class="bgcolcurrent"><span id="REF_DB_NAME" class="my-cell-static">%(REF_DB_NAME)s</span></td>
                 <td class="bgcolauthor"><span id="AUTH_DB_NAME_LIST" class="my-cell-static">%(AUTH_DB_NAME_LIST)s</span></td>
             </tr>

             <tr>
                <td>DB Code </td>
                 <td %(bgcolclass)s>%(db_code)s</td>
                <td class="bgcolcurrent"><span id="REF_DB_CODE" class="my-cell-static">%(REF_DB_CODE)s</span></td>
                <td class="bgcolauthor"><span id="AUTH_DB_CODE_LIST" class="my-cell-static">%(AUTH_DB_CODE_LIST)s</span></td>
             </tr>

             <tr>
                <td>DB Accession </td>
                <td %(bgcolclass)s>%(db_accession)s</td>
                <td class="bgcolcurrent"><span id="REF_DB_ACCESSION" class="my-cell-static">%(REF_DB_ACCESSION)s</span></td>
                <td class="bgcolauthor"><span id="AUTH_DB_ACCESSION_LIST" class="my-cell-static">%(AUTH_DB_ACCESSION_LIST)s</span></td>
             </tr>
%(db_align_range)s

             <!-- START  HERE -->
             <tr>
                <td>Name/Description </td>
                <td class="bgcolcurrent"><span id="ENTITY_DESCRIPTION" class="ief my-editable-cell %(ENTITY_DESCRIPTION_CSS)s">%(ENTITY_DESCRIPTION)s</span></td>
                <td class="bgcolcurrent"><span id="REF_DB_MOLECULE_NAME" class="my-cell-static">%(REF_DB_MOLECULE_NAME)s</span></td>
                <td class="bgcolauthor"><span id="ENTITY_DESCRIPTION_ORIG" class="my-cell-static">%(ENTITY_DESCRIPTION_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Synonyms </td>
                <td class="bgcolcurrent"><span id="ENTITY_SYNONYMS" class="ief my-editable-cell %(ENTITY_SYNONYMS_CSS)s">%(ENTITY_SYNONYMS)s</span></td>
                <td class="bgcolcurrent"><span id="REF_DB_MOLECULE_SYNONYMS" class="my-cell-static">%(REF_DB_MOLECULE_SYNONYMS)s</span></td>
                <td class="bgcolauthor"><span id="ENTITY_SYNONYMS_ORIG" class="my-cell-static">%(ENTITY_SYNONYMS_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Source Organism </td>
                <td class="bgcolcurrent"><span id="SOURCE_ORGANISM" class="ief my-editable-cell %(SOURCE_ORGANISM_CSS)s">%(SOURCE_ORGANISM)s</span></td>
                <td class="bgcolcurrent"><span id="REF_DB_SOURCE_ORGANISM" class="my-cell-static">%(REF_DB_SOURCE_ORGANISM)s</span></td>
                <td class="bgcolauthor"><span id="SOURCE_ORGANISM_ORIG" class="my-cell-static">%(SOURCE_ORGANISM_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Source Common Name </td>
                <td class="bgcolcurrent"><span id="SOURCE_COMMON_NAME" class="ief my-editable-cell %(SOURCE_COMMON_NAME_CSS)s">%(SOURCE_COMMON_NAME)s</span></td>
                <td class="bgcolcurrent"><span id="REF_DB_SOURCE_COMMON_NAME" class="my-cell-static">%(REF_DB_SOURCE_COMMON_NAME)s</span></td>
                <td class="bgcolauthor"><span id="SOURCE_COMMON_NAME_ORIG" class="my-cell-static">%(SOURCE_COMMON_NAME_ORIG)s</span></td>
             </tr>


             <tr>
                <td>Taxonomy ID</td>
                <td class="bgcolcurrent"><span id="SOURCE_TAXID" class="ief my-editable-cell %(SOURCE_TAXID_CSS)s">%(SOURCE_TAXID)s</span></td>
                <td class="bgcolcurrent"><span id="REF_DB_SOURCE_TAXID" class="my-cell-static">%(REF_DB_SOURCE_TAXID)s</span></td>
                <td class="bgcolauthor"><span id="SOURCE_TAXID_ORIG" class="my-cell-static">%(SOURCE_TAXID_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Gene Name</td>
                <td class="bgcolcurrent"><span id="SOURCE_GENE_NAME" class="ief my-editable-cell %(SOURCE_GENE_NAME_CSS)s">%(SOURCE_GENE_NAME)s</span></td>
                <td class="bgcolcurrent"><span id="REF_DB_GENE_NAME" class="my-cell-static">%(REF_DB_GENE_NAME)s</span></td>
                <td class="bgcolauthor"><span id="SOURCE_GENE_NAME_ORIG" class="my-cell-static">%(SOURCE_GENE_NAME_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Strain Name</td>
                <td class="bgcolcurrent"><span id="SOURCE_STRAIN"        class="ief my-editable-cell %(SOURCE_STRAIN_CSS)s">%(SOURCE_STRAIN)s</span></td>
                <td class="bgcolcurrent"><span id="REF_DB_SOURCE_STRAIN" class="my-cell-static">%(REF_DB_SOURCE_STRAIN)s</span></td>
                <td class="bgcolauthor"><span id="SOURCE_STRAIN_ORIG"   class="my-cell-static">%(SOURCE_STRAIN_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Source Variant</td>
                <td class="bgcolcurrent"><span id="SOURCE_VARIANT"        class="ief my-editable-cell %(SOURCE_VARIANT_CSS)s">%(SOURCE_VARIANT)s</span></td>
                <td></td>
                <td class="bgcolauthor"><span id="SOURCE_VARIANT_ORIG"   class="my-cell-static">%(SOURCE_VARIANT_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Enzyme Classification</td>
                <td class="bgcolcurrent"><span id="ENTITY_ENZYME_CLASS"        class="ief my-editable-cell %(ENTITY_ENZYME_CLASS_CSS)s">%(ENTITY_ENZYME_CLASS)s</span></td>
                <td class="bgcolcurrent"><span id="REF_DB_ENZYME_CLASS" class="my-cell-static">%(REF_DB_ENZYME_CLASS)s</span></td>
                <td class="bgcolauthor"><span id="ENTITY_ENZYME_CLASS_ORIG"   class="my-cell-static">%(ENTITY_ENZYME_CLASS_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Fragment Details</td>
                <td class="bgcolcurrent"><span id="ENTITY_FRAGMENT_DETAILS" class="ief my-editable-cell %(ENTITY_FRAGMENT_DETAILS_CSS)s">%(ENTITY_FRAGMENT_DETAILS)s</span></td>
                <td></td>
                <td class="bgcolauthor"><span id="ENTITY_FRAGMENT_DETAILS_ORIG" class="my-cell-static">%(ENTITY_FRAGMENT_DETAILS_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Mutation Details</td>
                <td class="bgcolcurrent"><span id="ENTITY_MUTATION_DETAILS" class="ief my-editable-cell %(ENTITY_MUTATION_DETAILS_CSS)s">%(ENTITY_MUTATION_DETAILS)s</span></td>
                <td></td>
                <td class="bgcolauthor"><span id="ENTITY_MUTATION_DETAILS_ORIG" class="my-cell-static">%(ENTITY_MUTATION_DETAILS_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Entity Details</td>
                <td class="bgcolcurrent"><span id="ENTITY_DETAILS" class="ief my-editable-cell %(ENTITY_DETAILS_CSS)s">%(ENTITY_DETAILS)s</span></td>
                <td></td>
                <td class="bgcolauthor"><span id="ENTITY_DETAILS_ORIG" class="my-cell-static">%(ENTITY_DETAILS_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Host Organism</td>
                <td class="bgcolcurrent"><span id="HOST_ORG_SOURCE" class="ief my-editable-cell %(HOST_ORG_SOURCE_CSS)s">%(HOST_ORG_SOURCE)s</span></td>
                <td></td>
                <td class="bgcolauthor"><span id="HOST_ORG_SOURCE_ORIG" class="my-cell-static">%(HOST_ORG_SOURCE_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Host Common Name</td>
                <td class="bgcolcurrent"><span id="HOST_ORG_COMMON_NAME" class="ief my-editable-cell %(HOST_ORG_COMMON_NAME_CSS)s">%(HOST_ORG_COMMON_NAME)s</span></td>
                <td></td>
                <td class="bgcolauthor"><span id="HOST_ORG_COMMON_NAME_ORIG" class="my-cell-static">%(HOST_ORG_COMMON_NAME_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Host Organism Strain</td>
                <td class="bgcolcurrent"><span id="HOST_ORG_STRAIN" class="ief my-editable-cell %(HOST_ORG_STRAIN_CSS)s">%(HOST_ORG_STRAIN)s</span></td>
                <td></td>
                <td class="bgcolauthor"><span id="HOST_ORG_STRAIN_ORIG" class="my-cell-static">%(HOST_ORG_STRAIN_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Host Organism TaxID</td>
                <td class="bgcolcurrent"><span id="HOST_ORG_TAXID" class="ief my-editable-cell %(HOST_ORG_TAXID_CSS)s">%(HOST_ORG_TAXID)s</span></td>
                <td></td>
                <td class="bgcolauthor"><span id="HOST_ORG_TAXID_ORIG" class="my-cell-static">%(HOST_ORG_TAXID_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Host Organism Vector</td>
                <td class="bgcolcurrent"><span id="HOST_ORG_VECTOR" class="ief my-editable-cell %(HOST_ORG_VECTOR_CSS)s">%(HOST_ORG_VECTOR)s</span></td>
                <td></td>
                <td class="bgcolauthor"><span id="HOST_ORG_VECTOR_ORIG" class="my-cell-static">%(HOST_ORG_VECTOR_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Host Organism Vector Type</td>
                <td class="bgcolcurrent"><span id="HOST_ORG_VECTOR_TYPE" class="ief my-editable-cell %(HOST_ORG_VECTOR_TYPE_CSS)s">%(HOST_ORG_VECTOR_TYPE)s</span></td>
                <td></td>
                <td class="bgcolauthor"><span id="HOST_ORG_VECTOR_TYPE_ORIG" class="my-cell-static">%(HOST_ORG_VECTOR_TYPE_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Host Organism Plasmid</td>
                <td class="bgcolcurrent"><span id="HOST_ORG_PLASMID" class="ief my-editable-cell %(HOST_ORG_PLASMID_CSS)s">%(HOST_ORG_PLASMID)s</span></td>
                <td></td>
                <td class="bgcolauthor"><span id="HOST_ORG_PLASMID_ORIG" class="my-cell-static">%(HOST_ORG_PLASMID_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Host Organism Cell Line</td>
                <td class="bgcolcurrent"><span id="HOST_ORG_CELL_LINE" class="ief my-editable-cell %(HOST_ORG_CELL_LINE_CSS)s">%(HOST_ORG_CELL_LINE)s</span></td>
                <td></td>
                <td class="bgcolauthor"><span id="HOST_ORG_CELL_LINE_ORIG" class="my-cell-static">%(HOST_ORG_CELL_LINE_ORIG)s</span></td>
             </tr>

            <tr>
                <td>Host Organism Variant</td>
                <td class="bgcolcurrent"><span id="HOST_ORG_VARIANT" class="ief my-editable-cell %(HOST_ORG_VARIANT_CSS)s">%(HOST_ORG_VARIANT)s</span></td>
                <td></td>
                <td class="bgcolauthor"><span id="HOST_ORG_VARIANT_ORIG" class="my-cell-static">%(HOST_ORG_VARIANT_ORIG)s</span></td>
             </tr>

             <tr>
                <td>Source Method</td>
                <td class="bgcolcurrent"><span id="SOURCE_METHOD" class="ief my-editable-cell" data-ief-edittype="select"  data-ief-selectvalues='%(SOURCE_METHOD_SELECT)s'>%(SOURCE_METHOD)s</span></td>
                <td></td>
                <td class="bgcolauthor"><span id="SOURCE_METHOD_ORIG" class="my-cell-static">%(SOURCE_METHOD_ORIG)s</span></td>
             <tr>
            </table>
            <input type="submit" name="submit" value="%(savebuttontext)s" class="disableonclick submitparentform fltrgt" />
            <!-- <input type="reset" name="reset" value="Reset" /> -->
        </form>
        """
        #
        rL = []
        #
        authSelectId, authSL, authFdObj = self.__getCurrentAuthSelection(entityId, partId=partId)
        refSelectId, refSL, refFdObj = self.__getCurrentRefSelection(entityId, partId=partId)
        #
        if authSL is None or authFdObj is None:
            return rL
        #
        authSeqType, authSeqInstId, authSeqPartId, authSeqAltId, authSeqVersion = authSL.get()
        authFD = authFdObj.get()
        #
        refFD = {}
        if refFdObj:
            refFD = refFdObj.get()
        #
        norineFlag = False
        # seqLength = ""
        if refSL and (numParts == 1) and (authFD["POLYMER_LINKING_TYPE"].find("peptide") != -1):
            seqIdx = self.__sds.getSequence(seqId=authSeqInstId, seqType=authSeqType, partId=authSeqPartId, altId=authSeqAltId, version=authSeqVersion)
            if len(seqIdx) < 25:
                norineFlag = True
                # seqLength = str(len(seqIdx))
            #
        #
        # This information is consolidated by entity here -
        authDbNameList, authDbAccessionList, authDbCodeList = self.__getAuthRefAssignments(entityId)
        #
        pD = {}
        if norineFlag:
            # allow annotator to manually input Norine reference information for small polymer sequence
            #
            for item in ("bgcolclass", "db_name", "db_code", "db_accession"):
                pD[item] = ""
            #
            pD["bgcolclass"] = 'class="bgcolcurrent"'
            db_name_class = "greyedout"
            db_name_value_select = '[{"value":"","label":"","selected":true},{"value":"NOR","label":"NOR","selected":false}]'
            db_name_value_curr = self.__placeHolderValue
            if ("ANNO_EDIT_DB_NAME" in authFD) and authFD["ANNO_EDIT_DB_NAME"].upper() == "NOR":
                db_name_class = ""
                db_name_value_select = '[{"value":"","label":"","selected":false},{"value":"NOR","label":"NOR","selected":true}]'
                db_name_value_curr = authFD["ANNO_EDIT_DB_NAME"].upper()
            #
            pD["db_name"] = (
                '<span id="ANNO_EDIT_DB_NAME" class="ief my-editable-cell '
                + db_name_class
                + '" data-ief-edittype="select"  data-ief-selectvalues=\''
                + db_name_value_select
                + "'>"
                + db_name_value_curr
                + "</span>"
            )
            #
            span_template = '<span id="%s" class="ief my-editable-cell %s">%s</span>'
            for items in (("db_code", "ANNO_EDIT_DB_CODE"), ("db_accession", "ANNO_EDIT_DB_ACCESSION")):
                styleClass = "greyedout"
                value = self.__placeHolderValue
                if (items[1] in authFD) and authFD[items[1]]:
                    styleClass = ""
                    value = authFD[items[1]]
                #
                pD[items[0]] = span_template % (items[1], styleClass, value)
            #
            tr_template = """

             <tr>
                <td>%s</td>
                <td class="bgcolcurrent"><span id="%s" class="ief my-editable-cell %s">%s</span></td>
                <td></td>
                <td></td>
             </tr>
            """
            pD["db_align_range"] = ""
            for items in (("DB Align Begin", "ANNO_EDIT_DB_ALIGN_BEGIN"), ("DB Align End", "ANNO_EDIT_DB_ALIGN_END")):
                styleClass = "greyedout"
                value = self.__placeHolderValue
                if (items[1] in authFD) and authFD[items[1]]:
                    styleClass = ""
                    value = authFD[items[1]]
                #
                pD["db_align_range"] += tr_template % (items[0], items[1], styleClass, value)
            #
        else:
            for item in ("bgcolclass", "db_name", "db_code", "db_accession", "db_align_range"):
                pD[item] = ""
            #
        #
        pD["AUTH_DB_NAME_LIST"] = ",".join(authDbNameList)
        pD["AUTH_DB_CODE_LIST"] = ",".join(authDbCodeList)
        pD["AUTH_DB_ACCESSION_LIST"] = ",".join(authDbAccessionList)

        if refFD:
            pD["REF_DB_NAME"] = refFD["DB_NAME"]
            pD["REF_DB_CODE"] = refFD["DB_CODE"]
            displayCode = refFD["DB_ACCESSION"]
            if refFD["DB_ISOFORM"] is not None and len(refFD["DB_ISOFORM"]) > 1:
                displayCode = refFD["DB_ISOFORM"]
            pD["REF_DB_ACCESSION"] = displayCode
            pD["REF_DB_MOLECULE_NAME"] = refFD["DB_MOLECULE_NAME"]
            pD["REF_DB_MOLECULE_SYNONYMS"] = refFD["DB_MOLECULE_SYNONYMS"]
            pD["REF_DB_SOURCE_ORGANISM"] = refFD["SOURCE_ORGANISM"]
            pD["REF_DB_SOURCE_COMMON_NAME"] = refFD["SOURCE_COMMON_NAME"]
            pD["REF_DB_SOURCE_TAXID"] = refFD["SOURCE_TAXID"]
            pD["REF_DB_SOURCE_STRAIN"] = refFD["SOURCE_STRAIN"]
            pD["REF_DB_GENE_NAME"] = refFD["DB_GENE_NAME"]
            pD["REF_DB_ENZYME_CLASS"] = refFD["DB_MOLECULE_EC"]
        else:
            for item in (
                "REF_DB_NAME",
                "REF_DB_CODE",
                "REF_DB_ACCESSION",
                "REF_DB_MOLECULE_NAME",
                "REF_DB_MOLECULE_SYNONYMS",
                "REF_DB_SOURCE_ORGANISM",
                "REF_DB_SOURCE_COMMON_NAME",
                "REF_DB_SOURCE_TAXID",
                "REF_DB_SOURCE_STRAIN",
                "REF_DB_GENE_NAME",
                "REF_DB_ENZYME_CLASS",
            ):
                pD[item] = ""
            #
        #
        if authFD["SOURCE_METHOD"] == "MAN":
            tup = ("true", "false", "false")
        elif authFD["SOURCE_METHOD"] == "NAT":
            tup = ("false", "true", "false")
        elif authFD["SOURCE_METHOD"] == "SYN":
            tup = ("false", "false", "true")
        else:
            tup = ("false", "false", "false")

        pD["SOURCE_METHOD_SELECT"] = '[{"value":"MAN","label":"MAN","selected":%s},{"value":"NAT","label":"NAT","selected":%s},{"value":"SYN","label":"SYN","selected":%s}]' % tup

        #
        pD["entryid"] = entryId
        pD["sessionid"] = self.__sessionId
        pD["partid"] = partId
        pD["entityid"] = entityId
        pD["instid"] = authSeqInstId
        pD["altid"] = authSeqAltId
        pD["versionid"] = authSeqVersion
        if numParts > 1:
            pD["savebuttontext"] = "Save Edits for Entity " + str(entityId) + " Part " + str(partId) + "  (save each part separately)"
        else:
            pD["savebuttontext"] = "Save Edits for Entity " + str(entityId)
        #
        #
        for item in self.__entityReviewItemList:
            kyCss = item + "_CSS"
            if authFD[item] is None or len(authFD[item]) < 1:
                authFD[item] = self.__placeHolderValue
                pD[kyCss] = "greyedout"
            else:
                pD[kyCss] = ""

        pD.update(authFD)
        pD["CURRENT_AUTH_SELECT_ID"] = authSelectId
        if refSelectId:
            pD["CURRENT_REF_SELECT_ID"] = refSelectId
        #
        pD["STRUCT_TITLE"] = self.__sds.getEntryDetail("STRUCT_TITLE")
        pD["CITATION_TITLE"] = self.__sds.getEntryDetail("CITATION_TITLE")
        #
        if self.__debug and refFD:
            self.__lfh.write("+UpdatePolymerEntitySourceDetails.makeEntityReview() reference sequence feature dictionary contents:\n")
            for k in sorted(refFD.keys()):
                v = refFD[k]
                self.__lfh.write(" ++++    %-40s   --  %s\n" % (k, v))
                self.__lfh.write("+UpdatePolymerEntitySourceDetails.makeEntityReview() author sequence feature dictionary contents:\n")
                for k2 in sorted(pD.keys()):
                    v = pD[k2]
                    self.__lfh.write(" ++++    %-40s   --  %s\n" % (k2, v))
                #
            #
        #
        rL = form_template % pD
        #
        return rL

    def entitySourceDetailsFormResponder(self):
        """Update the entity annotation from form data -
        Form data encoded in the input request object --
        Return True
        """
        try:
            partId = int(str(self.__reqObj.getValue("partid")))
            # entityId = self.__reqObj.getValue("entityid")
            instId = self.__reqObj.getValue("instid")
            altId = self.__reqObj.getValue("altid")
            versionId = self.__reqObj.getValue("versionid")
            #
            authFdObj = self.__sds.getFeatureObj(instId, seqType="auth", partId=partId, altId=altId, version=versionId)

            fDUpd = {}
            for ky in self.__entityReviewItemList:
                fDUpd[ky] = self.__reqObj.getValue(ky)
                if fDUpd[ky] == self.__placeHolderValue:
                    fDUpd[ky] = ""
                #
                if self.__verbose:
                    self.__lfh.write("+UpdatePolymerEntitySourceDetails.entityReviewResponder() ky=%40s   -- value=%r\n" % (ky, fDUpd[ky]))
                #
            #
            # mark this sequence feature object as edited - to screen for overwrites -
            #
            fDUpd["HAS_MANUAL_EDIT"] = True
            #
            authFdObj.set(fDUpd, resetAll=False)
            #
            self.__sds.setFeatureObj(authFdObj, instId, seqType="auth", partId=partId, altId=altId, version=versionId)

            if self.__verbose:
                self.__lfh.write(
                    "+UpdatePolymerEntitySourceDetails.entityReviewFormResponder() updating AUTH sequence %s alt %r part %r version %r\n" % (instId, altId, partId, versionId)
                )
                authFdObj.printIt(log=self.__lfh)
            #
            self.__sds.serialize()
        except:  # noqa: E722 pylint: disable=bare-except
            self.__lfh.write("+UpdatePolymerEntitySourceDetails.entityReviewResponder() failing\n")
            traceback.print_exc(file=self.__lfh)
        #
        return True

    def __getCurrentAuthSelection(self, entityId, partId=1):
        """Search selections for author sequence part 1. in order to establish the
        the selected version.   Then return the objects associated with the input
        partId.
        """
        try:
            instanceId = entityId
            sL = SequenceLabel()
            for sId in self.__selectIdList:
                if sId.startswith("auth"):
                    sL.unpack(sId)
                    _seqType, seqInstId, seqPartId, seqAltId, seqVersion = sL.get()
                    if instanceId == seqInstId and seqPartId == 1:
                        fdObj = self.__sds.getFeatureObj(seqInstId, seqType="auth", partId=partId, altId=seqAltId, version=seqVersion)
                        self.__lfh.write(
                            "+UpdatePolymerEntitySourceDetails._getCurrentAuthSelection() returns: entity %r instance %r partId %r altId %r version %r\n"
                            % (entityId, seqInstId, seqPartId, seqAltId, seqVersion)
                        )
                        return sId, sL, fdObj
                    #
                #
            #
        except:  # noqa: E722 pylint: disable=bare-except
            self.__lfh.write(
                "+UpdatePolymerEntitySourceDetails._getCurrentAuthSelection() failed for selectList %r entityId %r partId %r\n" % (self.__selectIdList, entityId, partId)
            )
            traceback.print_exc(file=self.__lfh)
        #
        self.__lfh.write("+UpdatePolymerEntitySourceDetails._getCurrentAuthSelection() no return for entityId %r partId %r\n" % (entityId, partId))
        return None, None, None

    def __getCurrentRefSelection(self, entityId, partId=1):
        try:
            instanceId = entityId
            sL = SequenceLabel()
            selfRefTarget = "selfref_%s_%d" % (entityId, partId)
            for sId in self.__selectIdList:
                if sId.startswith("ref"):
                    sL.unpack(sId)
                    _seqType, seqInstId, seqPartId, seqAltId, seqVersion = sL.get()
                    if instanceId == seqInstId and partId == seqPartId:
                        fObj = self.__sds.getFeatureObj(seqInstId, seqType="ref", partId=seqPartId, altId=seqAltId, version=seqVersion)
                        self.__lfh.write(
                            "+UpdatePolymerEntitySourceDetails._getCurrentRefSelection() returns: entity %r instance %r partId %r altId %r version %r\n"
                            % (entityId, seqInstId, seqPartId, seqAltId, seqVersion)
                        )
                        return sId, False, fObj
                    #
                elif sId.startswith(selfRefTarget):
                    # JDW JDW
                    sL.set(seqType="ref", seqInstId=instanceId, seqPartId=partId, seqAltId=1, seqVersion=1)
                    fObj = SequenceFeature()
                    return sId, True, fObj
                #
            #
        except:  # noqa: E722 pylint: disable=bare-except
            self.__lfh.write(
                "+UpdatePolymerEntitySourceDetails._getCurrentRefSelection() failed for selectList %r entityId %s partId %r\n" % (self.__selectIdList, entityId, partId)
            )
            traceback.print_exc(file=self.__lfh)
        #
        self.__lfh.write(
            "+UpdatePolymerEntitySourceDetails._getCurrentRefSelection() no return for entityId %r instanceId %r partId %r selectIdList %r\n"
            % (entityId, instanceId, partId, self.__selectIdList)
        )
        return None, False, None

    # def __getCurrentRefDetails(self, entityId, partId=1):
    #     # seqIdList = self.__sds.getGroup(entityId)
    #     # JDW CHANGE
    #     # seqId=seqIdList[0]
    #     seqId = entityId
    #     self.__lfh.write("+UpdatePolymerEntitySourceDetails.__getCurrentRefDetails() entityId %r partId %r seqId %r\n" % (entityId, partId, seqId))
    #     self.__lfh.write("+UpdatePolymerEntitySourceDetails.__getCurrentRefDetails() selectIdList %r\n" % (self.__selectIdList))
    #     sL = SequenceLabel()
    #     for sId in self.__selectIdList:
    #         if sId.startswith("ref"):
    #             sL.unpack(sId)
    #             seqType, seqInstId, seqPartId, seqAltId, seqVersion = sL.get()
    #             self.__lfh.write("+UpdatePolymerEntitySourceDetails.__getCurrentRefDetails() testing seqInstId %r seqPartId %r\n" % (seqInstId, seqPartId))
    #             if seqId == seqInstId and partId == seqPartId:
    #                 fObj = self.__sds.getFeatureObj(seqInstId, seqType="ref", partId=seqPartId, altId=seqAltId, version=seqVersion)
    #                 fD = fObj.get()
    #                 return fD["DB_NAME"], fD["DB_CODE"], fD["DB_ACCESSION"], fD["REF_MATCH_BEGIN"], fD["REF_MATCH_END"]

    #     return None, None, None, None, None

    def __getAuthRefAssignments(self, entityId):
        """Get author provided reference assignments.   Return a lists of reference database accessions.

        Note that this is currently being do a the level of
        entity as the entity parts are not explicit in input data at this time. 2013-11-14
        """
        depSeqAssignD = self.__sds.getDepositorReferenceAssignments()
        sADep = SequenceAssignDepositor(verbose=self.__verbose, log=self.__lfh)
        sADep.set(depSeqAssignD)
        refSeqAssignL = sADep.getReferenceList(entityId)

        if self.__debug:
            self.__lfh.write("+UpdatePolymerEntitySourceDetails.__getAuthRefAssignments() for entityId %r reference count is %d\n" % (entityId, sADep.getReferenceCount(entityId)))
            sADep.printIt(log=self.__lfh)
            for ii, rsa in enumerate(refSeqAssignL):
                self.__lfh.write("+UpdatePolymerEntitySourceDetails.__getAuthRefAssignments() depositor reference  %d\n" % (ii + 1))
                rsa.printIt(self.__lfh)
        #
        authDbNameList = []
        authDbAccessionList = []
        authDbCodeList = []
        for rsa in refSeqAssignL:
            dbName, dbCode, dbAccession = rsa.getDbReference()
            if dbName not in [".", "?"]:
                authDbNameList.append(dbName)
            else:
                authDbNameList.append("")
            if dbAccession not in [".", "?"]:
                authDbAccessionList.append(dbAccession)
            else:
                authDbAccessionList.append("")
            if dbCode not in [".", "?"]:
                authDbCodeList.append(dbCode)
            else:
                authDbCodeList.append("")
        #
        if self.__verbose:
            self.__lfh.write(
                "+UpdatePolymerEntitySourceDetails.__getAuthRefAssignments()  dbNames %r dbAccessions %r dbCodes %r\n" % (authDbNameList, authDbAccessionList, authDbCodeList)
            )

        return authDbNameList, authDbAccessionList, authDbCodeList

    def dump(self):
        """Output the depositor and archive sequence reference assignments."""
        entityIdList = self.__sds.getGroupIds()

        seqAssignD = self.__sds.getReferenceAssignments()
        sA = SequenceAssignArchive(verbose=self.__verbose, log=self.__lfh)
        sA.set(seqAssignD)
        sA.printIt(log=self.__lfh)

        for entityId in entityIdList:
            nRef = sA.getReferenceCount(entityId=entityId)
            self.__lfh.write("+UpdatePolymerEntitySourceDetails.dump() entityId %s archive assignment count %d\n" % (entityId, nRef))
            if nRef > 0:
                refL = sA.getReferenceList(entityId=entityId)
                for ref in refL:
                    ref.printIt(self.__lfh)
        #
        depSeqAssignD = self.__sds.getDepositorReferenceAssignments()
        sADep = SequenceAssignDepositor(verbose=self.__verbose, log=self.__lfh)
        sADep.set(depSeqAssignD)
        sADep.printIt(log=self.__lfh)

        for entityId in entityIdList:
            nRef = sADep.getReferenceCount(entityId=entityId)
            self.__lfh.write("+UpdatePolymerEntitySourceDetails.dump() entityId %s depositor assignment count %d\n" % (entityId, nRef))
            if nRef > 0:
                refL = sADep.getReferenceList(entityId=entityId)
                for ref in refL:
                    ref.printIt(self.__lfh)


if __name__ == "__main__":
    pass
