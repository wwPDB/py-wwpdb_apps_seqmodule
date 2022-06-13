##
# File:    PolymerLinkageDepict.py
# Date:    03-May-2010
#
# Updates:
#  03-Mar-2013 jdw refactored --
##
"""
The PolymerLinkDepiction() class provides methods for rendering summary sequence view.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.07"


import sys


class PolymerLinkageDepict(object):
    """The PolymerLinkDepiction() class provides methods for rendering a table of polymer linkage distances."""

    def __init__(self, lowerBound=1.2, upperBound=1.7, verbose=False, log=sys.stderr):  # pylint: disable=unused-argument
        # self.__verbose = verbose
        # self.__lfh = log
        #
        self.__upperB = upperBound
        self.__lowerB = lowerBound
        #

    def buildPolymerLinkageTable(self, rowDictList=None):
        """Render data in the input list of dictionaries representing the polymer
        linkage data..
        """
        if rowDictList is None:
            rowDictList = []
        iCount = 0
        for rD in rowDictList:
            dist = float(str(rD["dist"]))
            if dist > self.__upperB or dist < self.__lowerB:
                iCount += 1

        if iCount > 0:
            return self.__formatPolymerLinkageTable(rowDictList)
        else:
            oL = []
            oL.append("<div>")
            oL.append("<p> No atypical polymer linkages found in this entry. </p>")
            oL.append("</div>")
            return oL

    def __formatPolymerLinkageTable(self, rowDictList=None):
        """Render data in the input list of dictionaries representing the polymer"""
        if rowDictList is None:
            rowDictList = []

        columnNameList = ["Model", "Chain 1", "Residue 1", "Chain 2", "Residue 2", "Distance"]
        oL = []

        oL.append("<table>\n")
        #
        oL.append("<thead>\n")
        oL.append("<tr>")
        for col in columnNameList:
            oL.append("<th>%s</th>" % col)
        oL.append("</tr>\n")
        oL.append("</thead>\n")
        oL.append("<tbody>\n")

        for rD in rowDictList:
            model = rD["PDB_model_num"]
            ch1 = rD["auth_asym_id_1"]
            res1 = rD["auth_comp_id_1"] + " " + rD["auth_seq_id_1"] + rD["PDB_ins_code_1"]
            ch2 = rD["auth_asym_id_2"]
            res2 = rD["auth_comp_id_2"] + " " + rD["auth_seq_id_2"] + rD["PDB_ins_code_2"]
            dist = float(str(rD["dist"]))
            if dist > self.__upperB or dist < self.__lowerB:
                oL.append("<tr>")
                oL.append("<td>%s</td>" % str(model))
                oL.append("<td>%s</td>" % str(ch1))
                oL.append("<td>%s</td>" % str(res1))
                oL.append("<td>%s</td>" % str(ch2))
                oL.append("<td>%s</td>" % str(res2))
                oL.append("<td>%8.2f</td>" % dist)
                oL.append("</tr>\n")

        oL.append("</tbody>\n")
        oL.append("</table>\n")
        return oL
