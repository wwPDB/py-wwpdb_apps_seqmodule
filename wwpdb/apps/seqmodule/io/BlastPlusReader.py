##
# File:    BlastPlusReader.py
# Date:    17-Apr-2013  JDW
#
# Adapted from  seqdb_v2/ReadNcbiBlastXml.py
#
#
# Updates:
#
#  18-Apr-2013  jdw  modify to handle current conventions in UNP fasta files comment lines
#                    Extract isoform information from comment lines --
#  19-Apr-2013  jdw  remove any identity-based filtering of output
#  24-Oct-2022  zf   added _getUniprotInfo() method to parse Uniprot information from "Hit_def" tag.
##


from xml.dom import minidom
import sys
import traceback


class BlastPlusReader(object):
    """Read Blast+ result file (xml format)  and return the list of dictionaries containing the following --

    dict['db_name']
    dict['db_code']
    dict['db_accession']
    dict['db_isoform']
    dict['identity']
    dict['positive']
    dict['gaps']
    dict['midline']
    dict['query']
    dict['queryFrom']
    dict['queryTo']
    dict['subject']
    dict['hitFrom']
    dict['hitTo']
    dict['alignLen']
    dict['match_length']
    dict['query_length']
    dict['db_length']

    """

    def __init__(self, verbose=True, log=sys.stderr):
        self.__verbose = verbose
        # self.__debug = True
        self.__lfh = log
        self._resultList = []
        self.__sequenceType = "polypeptide(L)"

    def setSequenceType(self, polyType):
        if polyType in ["polypeptide(L)", "polypeptide(D)", "polypeptide", "polyribonucleotide"]:
            self.__sequenceType = polyType
            return True
        else:
            return False

    def readFile(self, filePath):
        try:
            if self.__verbose:
                self.__lfh.write("+BlastPlusReader.readFile() reading %s\n" % filePath)
            ifh = open(filePath, "r")
            inputText = ifh.read()
            ifh.close()
            if len(inputText) > 0:
                domObj = minidom.parseString(inputText)
                self._resultList = self._parse(domObj)
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+BlastPlusReader.readFile() failed for %s\n" % filePath)
                traceback.print_exc(file=self.__lfh)

        return self._resultList

    def GetResultList(self):
        return self._resultList

    def _parse(self, domObj):
        resultList = []
        dlist = domObj.getElementsByTagName("Hit")
        if not dlist:
            return resultList

        length = None
        for node in domObj.getElementsByTagName("BlastOutput_query-len"):
            if node.nodeType != node.ELEMENT_NODE:
                continue

            if node.tagName == "BlastOutput_query-len":
                length = node.firstChild.data
                break

        for node in dlist:
            if node.nodeType != node.ELEMENT_NODE:
                continue

            alignlist = self._ProcessHitTag(node.childNodes, length)
            if alignlist:
                for align in alignlist:
                    if length:
                        align["query_length"] = length
                    resultList.append(align)

        return resultList

    def _ProcessHitTag(self, nodelist, length):
        resultList = []
        alignlist = []
        description = ""
        length = ""
        # _code = ""
        dbName = ""
        dbCode = ""
        isoForm = ""
        for node in nodelist:
            if node.nodeType != node.ELEMENT_NODE:
                continue
            #
            if node.tagName == "Hit_id":
                dbName, dbCode, accCode, isoForm = self._parseID(node.firstChild.data)
                if not accCode:
                    return resultList
                #
            elif node.tagName == "Hit_accession":
                _code = node.firstChild.data  # noqa: F841
            elif node.tagName == "Hit_def":
                description = node.firstChild.data
            elif node.tagName == "Hit_len":
                length = node.firstChild.data
            elif node.tagName == "Hit_hsps":
                plist = self._ProcessHit_hspsTag(node.childNodes, length)
                if plist:
                    for li in plist:
                        alignlist.append(li)
                    #
                #
            #
        #
        if not accCode or not alignlist:
            return resultList
        #
        desInfo = {}
        if description and (dbName in ("SP", "TR", "UNP")):
            desInfo = self._getUniprotInfo(description)
        #
        for align in alignlist:
            align["db_name"] = dbName
            align["db_code"] = dbCode
            align["db_accession"] = accCode
            align["db_isoform"] = isoForm
            align["db_description"] = description
            align["db_length"] = length
            if desInfo:
                for key, val in desInfo.items():
                    align[key] = val
                #
            #
            resultList.append(align)
        #
        return resultList

    def _parseID(self, data):
        accCode = ""
        dbCode = ""
        dbName = ""
        isoForm = ""
        dlist = data.split("|")
        if len(dlist) >= 2 and dlist[2] != "pdb":
            f0 = str(dlist[0]).upper()
            if f0 in ["TR", "SP"]:
                dbName = f0
                accCode = str(dlist[1])
                dbCode = str(dlist[2])
                isoForm = ""
                if accCode is not None and (accCode.find("-") != -1):
                    tL = accCode.split("-")
                    if len(tL) > 1 and len(tL[1]) > 0:
                        #  Isoform is the acession + '-' + variant #
                        # isoForm=tL[1]
                        isoForm = accCode
                        accCode = tL[0]
                    #
                #
            elif f0 in ["GI"]:
                dbName = str(dlist[2]).upper()
                accCode = str(dlist[1])
                dbCode = str(dlist[3])
            #
        #
        return dbName, dbCode, accCode, isoForm

    def _ProcessHit_hspsTag(self, nodelist, length):
        resultList = []
        for node in nodelist:
            if node.nodeType != node.ELEMENT_NODE:
                continue
            #
            if node.tagName != "Hsp":
                continue
            #
            adict = self._GetMatchAlignment(node.childNodes, length)
            if adict:
                resultList.append(adict)
            #
        #
        return resultList

    def _GetMatchAlignment(self, nodelist, length):  # pylint: disable=unused-argument
        rdict = {}
        hsp_num = ""
        for node in nodelist:
            if node.nodeType != node.ELEMENT_NODE:
                continue
            #
            if node.tagName == "Hsp_identity":
                rdict["identity"] = node.firstChild.data
            elif node.tagName == "Hsp_positive":
                rdict["positive"] = node.firstChild.data
            elif node.tagName == "Hsp_gaps":
                rdict["gaps"] = node.firstChild.data
            elif node.tagName == "Hsp_midline":
                if self.__sequenceType == "polyribonucleotide":
                    rdict["midline"] = node.firstChild.data.replace("T", "U")
                else:
                    rdict["midline"] = node.firstChild.data
                #
            elif node.tagName == "Hsp_qseq":
                if self.__sequenceType == "polyribonucleotide":
                    rdict["query"] = node.firstChild.data.replace("T", "U")
                else:
                    rdict["query"] = node.firstChild.data
                #
            elif node.tagName == "Hsp_query-from":
                rdict["queryFrom"] = node.firstChild.data
            elif node.tagName == "Hsp_query-to":
                rdict["queryTo"] = node.firstChild.data
            elif node.tagName == "Hsp_hseq":
                if self.__sequenceType == "polyribonucleotide":
                    rdict["subject"] = node.firstChild.data.replace("T", "U")
                else:
                    rdict["subject"] = node.firstChild.data
                #
            elif node.tagName == "Hsp_hit-from":
                rdict["hitFrom"] = node.firstChild.data
            elif node.tagName == "Hsp_hit-to":
                rdict["hitTo"] = node.firstChild.data
            elif node.tagName == "Hsp_align-len":
                rdict["alignLen"] = node.firstChild.data
                rdict["match_length"] = node.firstChild.data
            elif node.tagName == "Hsp_num":
                hsp_num = node.firstChild.data
            #
        # only take the first alignment within the hit
        if str(hsp_num) != "1":
            rdict.clear()
        #
        return rdict

    def _getUniprotInfo(self, description):
        desInfo = {}
        for tokenTupL in (("SV", ""), ("PE", ""), ("GN", "gene"), ("OX", "taxonomy_id"), ("OS", "source_scientific")):
            idx = description.find(tokenTupL[0] + "=")
            if idx < 0:
                continue
            #
            val = description[idx + 3:].strip()
            description = description[:idx].strip()
            if tokenTupL[1] and val:
                if tokenTupL[1] == "source_scientific":
                    idx1 = val.find("(strain")
                    if (val[-1] == ")") and idx1 >= 0:
                        source_scientific = val[:idx1].strip()
                        if source_scientific:
                            desInfo[tokenTupL[1]] = source_scientific
                        #
                        strain = val[idx1 + 7:-1].strip()
                        if strain:
                            desInfo["strain"] = strain
                        #
                    else:
                        desInfo[tokenTupL[1]] = val
                #
                elif tokenTupL[1] == "taxonomy_id":
                    if val.isnumeric():
                        desInfo[tokenTupL[1]] = val
                    #
                else:
                    desInfo[tokenTupL[1]] = val
                #
            #
        #
        if description:
            wList = description.split(" ")
            if (len(wList) > 3) and (wList[0] == "Isoform") and wList[1].isnumeric() and (wList[2] == "of"):
                desInfo["name"] = " ".join(wList[3:])
            else:
                desInfo["name"] = description
            #
        #
        return desInfo
