##
# File:    SummaryViewDepiction.py
# Date:    21-Jan-2010
#
# Updates:
# 20-Apr-2010 jdw     Port to module seqmodule.
# 04-Aug-2010, RPS    Updating formatTable() to accommodate tablesorting needs of table used
#                     for displaying coordinate sequences as well as reference sequences.
# 06-Aug-2010, RPS    formatTable()-corrected logic error in if-elif-else block used for tablesorter class handling
#
# 07-Mar-2013 jdw     Reorganize and generatlize
# 08-Mar-2013 jdw     Add support for entity parts
# 12-Mar-2013 jdw     add identifier reference section for call back -
# 04-Nov-2013 jdw     make select reference first in order -
# 15-Nov-2013 jdw     revise author content rendering
# 30-Nov-2013 jdw     retain state of active entity
# 15-May-2014 jdw     Add radio groups for versions coordinate sequence instances
# 15-Jun-2014 jdw     Make coordinate sequence instances checkboxes again owing to case discrimination of name attributes
# 03-May-2016 ep      Move the more sequences button to template before the coordinates
# 17-May-2016 ep      Change the showmore to a class as broke multiple entity case of more sequences button. Same for srchmore
##
"""
The SummaryViewDepiction() class provides methods for rendering summary sequence view.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.07"


import sys


class SummaryViewDepiction(object):

    """The SummaryViewDepiction() class provides methods for rendering summary sequence view.

    Data requirements -

    For each sequence type (e.g. ref, auth, xyz) within each sequence group -

    summaryDataObj = [groupId][seqType]  = {<content dictionary>}

                     (for seqType == 'ref' [{},{},...] for each sequence part matched)

    The content dictionary containing the following keys -

                          [ROW_IDS]       =[]  (sequence identifier)
                          [ROW_STATUS]    =[]  (selected or aligned)
                          [ROW_DATA_DICT] = [{},{},...]  contains the data


    """

    def __init__(self, verbose=False, log=sys.stderr):
        self.__verbose = verbose  # pylint: disable=unused-private-member
        self.__lfh = log  # pylint: disable=unused-private-member
        #

    def isRealValue(self, val):
        if val is not None:
            if (len(str(val)) > 0) and (str(val) not in [".", "?"]):
                # real value
                return True
            else:
                return False
        else:
            return False

    def buildSummaryView(self, summaryDataObj=None, activeGroupId=None):
        """Render data in the input summaryDataObj  -"""
        if summaryDataObj is None:
            summaryDataObj = {}
        oL = []
        gIdList = list(summaryDataObj.keys())

        gIdList.sort(key=int)
        gId0 = activeGroupId if activeGroupId in gIdList else gIdList[0]

        oL.append('<span class="padtop6 fltlft"><b>Select Entity:</b></span>\n<div id="pagi" class="noprint fltlft"></div>\n<br class="clearfloat" />\n')
        for gId in gIdList:

            tabStyle = "displaynone"
            if gId == gId0:
                tabStyle = "_current"

            gObj = summaryDataObj[gId]

            oL.append('<div id="p%s" class="%s tabscount">\n' % (str(gId), tabStyle))

            oL.append("<h3>Sequence Entity %s</h3>\n" % str(gId))

            oL.append('<div class="inneraccordion">\n')

            # author section -
            if "auth" in gObj:
                # oL.append('<a id="auth_entity_%s" class="head ui-corner-all" href="#"><span class="ui-icon ui-icon-circle-arrow-s fltlft"></span> Author</a>\n' % str(gId))
                oL.append('<div class="head ui-corner-all">\n')
                oL.append('<div class="fltlft width50">\n')
                oL.append('<a id="auth_entity_%s" href="#" class="toggle">' % str(gId))
                oL.append('<span class="fltlft width20px"><span class="ui-icon ui-icon-circle-arrow-s"></span></span><span>Author</span>')
                oL.append("</a>\n")
                oL.append("</div>\n")

                oL.append('<div class="fltrgt width50 txtr">')
                oL.append('<span class="colauthor">|&nbsp;Author</span>&nbsp;|&nbsp;<span class="colcurrent">Current&nbsp;|</span>&nbsp;&nbsp;')
                oL.append('<button class="loadtaxonomy">Edit Sequence/Partition/Taxonomy</button>')
                oL.append('<button class="loadentityreview">Edit Source Details</button>')
                oL.append("</div>\n")
                oL.append('<br class="clearfloat" />\n')
                oL.append("</div>")

                oL.append("<div>\n")
                gD = gObj["auth"]
                if gD["ENTITY_NUM_PARTS"] > 1:
                    columnNameList = ["Version", "Length", "Range", "Source&nbsp;Organism", "Name/Description", "Host"]
                    rowDataNameList = ["ROW_VERSION", "ROW_SEQ_LENGTH", "partdetails", "sourceAndStrain", "description", "hostorg"]
                    oL.extend(self.formatTable(gId, "auth", columnNameList, gD["ROW_IDS"], rowDataNameList, gD["ROW_STATUS"], gD["ROW_DATA_DICT"]))
                else:
                    columnNameList = ["Version", "Length", "Source&nbsp;Organism", "Name/Description", "Host"]
                    rowDataNameList = ["ROW_VERSION", "ROW_SEQ_LENGTH", "sourceAndStrain", "description", "hostorg"]
                    oL.extend(self.formatTable(gId, "auth", columnNameList, gD["ROW_IDS"], rowDataNameList, gD["ROW_STATUS"], gD["ROW_DATA_DICT"]))

                oL.append("</div>\n")
                # end of author section

            # reference section -
            if "ref" in gObj:
                gDList = gObj["ref"]
                for (partId, gD) in gDList:
                    oL.append('<div class="head ui-corner-all">\n')
                    oL.append('<div class="fltlft width50">\n')
                    oL.append('<a id="ref_%s_part_%s" href="#" class="toggle">' % (str(gId), str(partId)))
                    oL.append('<span class="fltlft width20px"><span class="ui-icon ui-icon-circle-arrow-s"></span></span><span>Reference for part&nbsp;%s</span>' % str(partId))
                    oL.append("</a>\n")
                    oL.append("</div>\n")
                    oL.append('<div class="fltrgt width50 txtr">')
                    oL.append('<button class="loadseqdbref">Load Reference Sequence</button>')
                    oL.append("</div>\n")
                    oL.append('<br class="clearfloat" />\n')
                    oL.append("</div>")

                    oL.append('<div class="refdiv">\n')

                    # columnNameList=["ID Code","Version","Length","Source Organism","Tax ID","Identity<br />(w/ gaps)","Alignment Details"]
                    # rowDataNameList=["ROW_ID_CODE","ROW_VERSION","ROW_SEQ_LENGTH","SOURCE_ORGANISM","SOURCE_TAXID","ROW_AUTH_REF_SIM","ROW_DETAIL_STRING"]
                    # columnNameList = ["ID Code", "Version", "Length", "Molecular Details", "Tax ID", "Identity<br />(w/ gaps)", "Alignment Details"]
                    # rowDataNameList = ["ROW_ID_CODE", "ROW_VERSION", "ROW_SEQ_LENGTH", "ROW_FEATURE_STRING", "SOURCE_TAXID", "ROW_AUTH_REF_SIM", "ROW_DETAIL_STRING"]
                    columnNameList = ["ID Code", "Version", "Molecular Details", "Tax ID", "Identity", "Alignment Details"]
                    rowDataNameList = ["ROW_ID_CODE", "ROW_VERSION", "ROW_FEATURE_STRING", "SOURCE_TAXID", "ROW_AUTH_REF_SIM", "ROW_DETAIL_STRING"]  #
                    oL.extend(
                        self.formatTable(
                            gId, "ref", columnNameList, gD["ROW_IDS"], rowDataNameList, gD["ROW_STATUS"], gD["ROW_DATA_DICT"], selfRefFlag=gD["SELF_REFERENCE_FLAG"], partId=partId
                        )
                    )
                    oL.append("</div>\n")

            # end of reference section
            oL.append("<div>\n")
            oL.append('<br class="clearfloat" />\n')
            oL.append("</div>\n")

            oL.append('<div class="srchdiv ui-state-highlight ui-corner-all noprint displaynone" >\n')
            oL.append('     <input type="button" name="showmore" class="showmore fltrgt" value="More Reference Sequences" />\n')
            oL.append('     <br class="clearfloat" />\n')
            oL.append("</div>\n")

            # coordinate section -
            if "xyz" in gObj:
                # oL.append('<a class="head ui-corner-all" href="#"><span class="ui-icon ui-icon-circle-arrow-s fltlft"></span> Coordinate</a>\n')
                oL.append('<div class="head ui-corner-all">\n')
                oL.append("<div>\n")
                oL.append('<a href="#" class="toggle">')
                oL.append('<span class="fltlft width20px"><span class="ui-icon ui-icon-circle-arrow-e"></span></span><span>Coordinate</span>')
                oL.append("</a>\n")
                oL.append("</div>\n")
                oL.append("</div>")
                oL.append('<div style="display: none;">\n')
                gD = gObj["xyz"]
                columnNameList = ["Chain ID", "Version", "Length", "Alignment Details"]
                rowDataNameList = ["ROW_ID_CODE", "ROW_VERSION", "ROW_SEQ_LENGTH", "ROW_DETAIL_STRING"]
                oL.extend(self.formatTable(gId, "xyz", columnNameList, gD["ROW_IDS"], rowDataNameList, gD["ROW_STATUS"], gD["ROW_DATA_DICT"]))
                oL.append("</div>\n")
                # end of coordinate section

            # end of sequence group
            oL.append("</div>\n")
            oL.append("</div>\n")
        #
        return oL

    def formatTable(self, gId, seqType, columnNameList, rowIdList, rowDataNameList, rowStatusList, rowDataDictList, selfRefFlag=False, partId=1):
        """Render summary tables of author, coordinate and reference sequences."""
        labelAlignAllType = "align_" + str(gId) + "_" + seqType
        labelSelectWithinType = "sel_" + str(gId) + "_" + seqType + "_" + str(partId)
        labelSelectAllType = "sel_" + str(gId) + "_" + seqType + "_" + str(partId)
        oL = []

        if seqType == "ref":
            oL.append('<table class="tablesorter_ref">\n')
        elif seqType == "xyz":
            oL.append('<table class="tablesorter_xyz">\n')
        else:
            oL.append("<table>\n")
        #
        #  Table header ----
        oL.append("<thead>\n")

        if seqType == "ref":
            if selfRefFlag:
                isSelected = 'checked="checked"'
            else:
                isSelected = ""
            oL.append('<tr id="selfref_%s_%s">' % (str(gId), str(partId)))
            oL.append('<th>Select<br /><input type="radio" name="%s" class="saveeles refselection" %s /> Self</th>' % (labelSelectWithinType, isSelected))
            oL.append('<th><input type="checkbox" name="%s" id="%s" class="checkall" />' % (labelAlignAllType, labelAlignAllType))
            oL.append('<label for="%s">Align all</label></th>' % labelAlignAllType)
            for col in columnNameList:
                oL.append("<th>%s</th>" % col)

            oL.append("</tr>\n")

        elif seqType == "xyz":
            oL.append("<tr>")
            oL.append('<th><input type="checkbox" name="%s" id="%s" class="selectall" />' % (labelSelectAllType, labelSelectAllType))
            oL.append('<label for="%s">Select all</label></th>' % labelSelectAllType)
            oL.append('<th><input type="checkbox" name="%s" id="%s" class="checkall" />' % (labelAlignAllType, labelAlignAllType))
            oL.append('<label for="%s">Align all</label></th>' % labelAlignAllType)
            for col in columnNameList:
                oL.append("<th>%s</th>" % col)

            oL.append("</tr>\n")

        elif seqType == "auth":
            oL.append("<tr>")
            oL.append("<th>Select</th>")
            oL.append('<th><input type="checkbox" name="%s" id="%s" class="checkall" />' % (labelAlignAllType, labelAlignAllType))
            oL.append('<label for="%s">Align all</label></th>' % labelAlignAllType)

            for col in columnNameList:
                if col.startswith("Source"):
                    oL.append('<th class="width20">%s</th>' % col)
                elif col.startswith("Range"):
                    oL.append('<th class="width12">%s</th>' % col)
                elif col.startswith("Host"):
                    oL.append('<th class="width20">%s</th>' % col)
                elif col.startswith("Name"):
                    oL.append('<th class="width30">%s</th>' % col)
                else:
                    oL.append("<th>%s</th>" % col)
            oL.append("</tr>\n")

        oL.append("</thead>\n")

        #  Table body  ----
        oL.append("<tbody>\n")
        if seqType == "auth":
            # ---
            #     columnNameList=["Version","Length","Source Organism","Name/Description","Host"]
            #     rowDataNameList=["ROW_VERSION","ROW_SEQ_LENGTH","sourceAndStrain","description","hostorg"]
            #  Details[ver][partId]['current']={'partdetails','sourceAndStrain','description','hostorg'}
            #  Details[ver][partId]['author']={'partdetails','sourceAndStrain','description','hostorg'}
            #
            #  if parts == 1
            #  columns -  select --  align -- version -- seq-length -- Author  -->>  sourceStrain -- Name/Descrip  -- HostOrg --
            #                                                       -- Assigned -->> sourceStrain -- Name/Descrip  -- HostOrg --
            # ---
            #   ---- For each version ---
            for ir in range(0, len(rowDataDictList)):
                fD = rowDataDictList[ir]
                numParts = fD["numparts"]
                if numParts == 1:
                    partId = 1
                    oL.append('<tr id="%s" >' % rowIdList[ir])

                    if rowStatusList[ir][0]:
                        isSelected = 'checked="checked"'
                    else:
                        isSelected = ""

                    oL.append('<td rowspan=2><input type="radio" name="%s" class="saveeles" %s /></td>' % (labelSelectWithinType, isSelected))

                    if rowStatusList[ir][1]:
                        isSelected = 'checked="checked"'
                    else:
                        isSelected = ""

                    oL.append('<td rowspan=2><input type="checkbox" name="align_%s" class="aligneles" %s /></td>' % (rowIdList[ir], isSelected))

                    oL.append("<td rowspan=2>%d</td>" % fD["ROW_VERSION"])
                    oL.append("<td rowspan=2>%d</td>" % fD["ROW_SEQ_LENGTH"])

                    for val in rowDataNameList[2:]:
                        oL.append('<td class="textleft bgcolauthor">')
                        if self.isRealValue(fD[partId]["author"][val]):
                            oL.append(str(fD[partId]["author"][val]))
                        oL.append("</td>")
                    oL.append("</tr>\n")
                    #
                    oL.append("<tr>")
                    for val in rowDataNameList[2:]:
                        oL.append('<td class="textleft bgcolcurrent">')
                        if self.isRealValue(fD[partId]["current"][val]):
                            oL.append(str(fD[partId]["current"][val]))
                        oL.append("</td>")
                    oL.append("</tr>\n")

                elif numParts > 1:
                    rowSpan = 2 * numParts
                    oL.append('<tr id="%s" >' % rowIdList[ir])

                    if rowStatusList[ir][0]:
                        isSelected = 'checked="checked"'
                    else:
                        isSelected = ""

                    oL.append('<td rowspan=%d><input type="radio" name="%s" class="saveeles" %s /></td>' % (rowSpan, labelSelectWithinType, isSelected))

                    if rowStatusList[ir][1]:
                        isSelected = 'checked="checked"'
                    else:
                        isSelected = ""

                    oL.append('<td rowspan=%d><input type="checkbox" name="align_%s" class="aligneles" %s /></td>' % (rowSpan, rowIdList[ir], isSelected))
                    oL.append("<td rowspan=%d>%d</td>" % (rowSpan, fD["ROW_VERSION"]))
                    oL.append("<td rowspan=%d>%d</td>" % (rowSpan, fD["ROW_SEQ_LENGTH"]))
                    #

                    rowSpan = 2
                    for partId in range(1, numParts + 1):
                        if partId > 1:
                            oL.append("<tr>")

                        oL.append("<td rowspan=%d>%s</td>" % (rowSpan, fD[partId]["current"]["partdetails"]))

                        for val in rowDataNameList[3:]:
                            oL.append('<td class="textleft bgcolauthor">')
                            if self.isRealValue(fD[partId]["author"][val]):
                                oL.append(str(fD[partId]["author"][val]))
                            oL.append("</td>")
                        oL.append("</tr>\n")
                        #
                        oL.append("<tr>")
                        for val in rowDataNameList[3:]:
                            oL.append('<td class="textleft bgcolcurrent">')
                            if self.isRealValue(fD[partId]["current"][val]):
                                oL.append(str(fD[partId]["current"][val]))
                            oL.append("</td>")
                        oL.append("</tr>\n")

        elif seqType == "xyz":
            for ir in range(0, len(rowDataDictList)):
                oL.append('<tr id="%s" >' % rowIdList[ir])
                if rowStatusList[ir][0]:
                    isSelected = 'checked="checked"'
                else:
                    isSelected = ""
                #
                sg = rowStatusList[ir][2]
                vr = rowStatusList[ir][3]
                # oL.append('<td><input type="checkbox" name="%s" class="saveeles" %s /></td>'     %   (labelSelectWithinType,isSelected))
                # oL.append('<td><input type="radio" name="%s_%s" class="saveeles  selectablexyz %s" %s /></td>'     %   (labelSelectWithinType,sg,vr,isSelected))
                oL.append('<td><input type="checkbox" name="%s_%s" class="saveeles  selectablexyz %s" %s /></td>' % (labelSelectWithinType, sg, vr, isSelected))
                #
                if rowStatusList[ir][1]:
                    isSelected = 'checked="checked"'
                else:
                    isSelected = ""
                oL.append('<td><input type="checkbox" name="align_%s" class="aligneles" %s /></td>' % (rowIdList[ir], isSelected))

                fD = rowDataDictList[ir]
                for val in rowDataNameList:
                    oL.append("<td>%s</td>" % fD[val])
                oL.append("</tr>\n")

        elif seqType == "ref":
            for ir in range(0, len(rowDataDictList)):
                if not rowStatusList[ir][0]:
                    continue
                oL.append('<tr id="%s" >' % rowIdList[ir])
                if rowStatusList[ir][0]:
                    isSelected = 'checked="checked"'
                else:
                    isSelected = ""
                oL.append('<td><input type="radio" name="%s" class="saveeles refselection" %s /></td>' % (labelSelectWithinType, isSelected))

                if rowStatusList[ir][1]:
                    isSelected = 'checked="checked"'
                else:
                    isSelected = ""
                oL.append('<td><input type="checkbox" name="align_%s" class="aligneles" %s /></td>' % (rowIdList[ir], isSelected))

                fD = rowDataDictList[ir]
                for val in rowDataNameList:
                    if val in ["ROW_FEATURE_STRING"]:
                        oL.append('<td class="textleft">%s</td>' % fD[val])
                    else:
                        oL.append("<td>%s</td>" % fD[val])
                oL.append("</tr>\n")

            for ir in range(0, len(rowDataDictList)):
                if rowStatusList[ir][0]:
                    continue
                oL.append('<tr id="%s" >' % rowIdList[ir])
                if rowStatusList[ir][0]:
                    isSelected = 'checked="checked"'
                else:
                    isSelected = ""
                oL.append('<td><input type="radio" name="%s" class="saveeles refselection" %s /></td>' % (labelSelectWithinType, isSelected))

                if rowStatusList[ir][1]:
                    isSelected = 'checked="checked"'
                else:
                    isSelected = ""
                oL.append('<td><input type="checkbox" name="align_%s" class="aligneles" %s /></td>' % (rowIdList[ir], isSelected))

                fD = rowDataDictList[ir]
                for val in rowDataNameList:
                    if val in ["ROW_FEATURE_STRING"]:
                        oL.append('<td class="textleft">%s</td>' % fD[val])
                    else:
                        oL.append("<td>%s</td>" % fD[val])
                oL.append("</tr>\n")

        else:
            pass

        oL.append("</tbody>\n")
        oL.append("</table>\n")
        return oL
