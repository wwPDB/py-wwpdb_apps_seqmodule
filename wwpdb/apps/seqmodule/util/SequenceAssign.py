##
# File:    SequenceAssign.py
# Date:    18-Mar-2013
#
# Updates:
# 19-Mar-2013 jdw  Add ReferenceSequenceAssign()
# 12-Feb-2014 jdw  expand attributes returned for ReferenceSequence()
# 29-Jul-2014 jdw  add db_isoform_description
#
##
"""
Containers for archive and depositor sequence assignment details.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.07"

import sys
import re
import traceback

# from operator import itemgetter
from wwpdb.apps.seqmodule.util.SequenceReferenceData import SequenceReferenceData


class ReferenceSequence(object):
    def __init__(self, verbose=False, log=sys.stderr):
        """Provides methods which encapsulate results obtained from fetches of sequence entries from GB and UNP.

          Selected content include:

          GB ['taxonomy_id', 'source_scientific', 'sequence']
          UNP  dictionary :

        dict['db_code']           - code
        dict['db_accession']      - first accession code
        dict['db_name']           -
        dict['db_isoform']        - isoform

        dict['name']              - protein name
        dict['keyword']           - keywords
        dict['sequence']          - sequence
        dict['comments']          - Uniprot comments
        dict['synonyms']          - protein synonyms
        dict['source_scientific'] - source scientific name
        dict['source_strain']     - source strain (derived)

        dict['taxonomy_id']       - source taxonomy ID
        dict['gene']              - gene

        dict['source_common']     - source common name
        dict['ec']                - EC number(s)



          This class adds  - source_strain if it is embedded in the scientific name  -
        """
        self.__verbose = verbose
        self.__lfh = log
        self.__D = {}

        self.__attribStr = [
            "db_name",
            "db_accession",
            "db_code",
            "db_isoform",
            "db_isoform_description",
            "name",
            "keyword",
            "sequence",
            "comments",
            "synonyms",
            "source_scientific",
            "source_strain",
            "taxonomy_id",
            "gene",
            "source_common",
            "ec",
        ]
        self.__attribInt = ["seq_length"]

    def clear(self):
        self.__D = {}
        for attrib in self.__attribStr:
            self.__D[attrib] = ""
        for attrib in self.__attribInt:
            self.__D[attrib] = 0

    def set(self, seqD):
        try:
            for k, v in seqD.items():
                if k in self.__attribStr:
                    if k in ["sequence"]:
                        self.__D[k] = re.sub("[\t \n]", "", str(v))
                    else:
                        self.__D[k] = str(v)
                elif k in self.__attribInt:
                    self.__D[k] = int(str(v))

            if (len(self.__D["source_scientific"]) > 0) and (self.__D["source_scientific"].find("(") != -1):
                s1, s2 = self.__getStrain(sourceName=self.__D["source_scientific"])
                self.__D["source_scientific"] = s1
                self.__D["source_strain"] = s2

            self.__D["seq_length"] = len(self.__D["sequence"])
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+ReferenceSequence.__set() failed %r\n" % seqD.items())
                traceback.print_exc(file=self.__lfh)
        return False

    def getSequenceLength(self):
        return self.__D["seq_length"]

    #

    def getEnzymeClass(self):
        try:
            return self.__D["ec"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getKeywords(self):
        try:
            return self.__D["keyword"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getComments(self):
        try:
            return self.__D["comment"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    #
    def getName(self):
        try:
            return self.__D["name"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getSynonyms(self):
        try:
            return self.__D["synonyms"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getGeneName(self):
        try:
            return self.__D["gene"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getDbIsoformDescription(self):
        try:
            return self.__D["db_isoform_description"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    #
    def getSourceName(self):
        try:
            return self.__D["source_scientific"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getSourceCommonName(self):
        try:
            return self.__D["source_common"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getSourceStrain(self):
        try:
            return self.__D["source_strain"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getSequence(self):
        """Get the full one-letter-code sequence"""
        try:
            return self.__D["sequence"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def XgetSequenceWithIndex(self, polyTypeCode="AA", seqBegin=1, seqEnd=None):
        try:
            sdr = SequenceReferenceData(verbose=self.__verbose, log=self.__lfh)
            if seqEnd is None or seqEnd > len(self.__D["sequence"]):
                seqEnd = len(self.__D["sequence"])
            if seqBegin is None or seqBegin < 1:
                seqBegin = 1
            return sdr.cnv1To3ListIdx(self.__D["sequence"][seqBegin - 1 : seqEnd], seqBegin, polyTypeCode)
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+ReferenceSequence.getSequenceWithIndex() failed seqBegin %r seqEnd %r\n" % (seqBegin, seqEnd))
                traceback.print_exc(file=self.__lfh)
            return []

    def getSequenceWithIndex(self, polyTypeCode="AA", seqBegin=1, seqEnd=None):
        """Convert the one-letter code sequence from the reference resource to internal indexed list
        format seqIdx=[(3-letter-code, ref-db-index, position in sequence (1-length), (),  ... ]

        If seqBegin > seqEnd the sense of the sequence must be reversed.
        """
        try:
            sdr = SequenceReferenceData(verbose=self.__verbose, log=self.__lfh)
            if seqEnd is None or seqEnd > len(self.__D["sequence"]):
                seqEnd = len(self.__D["sequence"])
            if seqBegin is None or seqBegin < 1:
                seqBegin = 1
            if seqBegin < seqEnd:
                return sdr.cnv1To3ListIdx(self.__D["sequence"][seqBegin - 1 : seqEnd], seqBegin, polyTypeCode)
            else:
                # invert ranges on slice -
                # sTupL=sdr.cnv1To3ListIdx(self.__D['sequence'][seqBegin-1:seqEnd] ,seqBegin,polyTypeCode)
                #
                s1S = self.__D["sequence"][seqEnd - 1 : seqBegin]
                c1S = sdr.compliment1NA(s1S, polyTypeCode)

                sTupL = sdr.cnv1To3ListIdx(c1S, seqEnd, polyTypeCode)
                # reverse the order in and reindex in the sequence tuple list
                rTupL = []
                idx = 1
                for sTup in sTupL[::-1]:
                    rTupL.append((sTup[0], sTup[1], sTup[2], idx))
                    idx += 1
                if self.__verbose:
                    for sTup in rTupL:
                        self.__lfh.write("+ReferenceSequence.getSequenceWithIndex() %r %r %r %r\n" % (sTup[0], sTup[1], sTup[2], sTup[3]))
                return rTupL
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+ReferenceSequence.getSequenceWithIndex() failed seqBegin %r seqEnd %r\n" % (seqBegin, seqEnd))
                traceback.print_exc(file=self.__lfh)
            return []

    def getTaxId(self):
        try:
            return self.__D["taxonomy_id"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getIsoform(self):
        try:
            return self.__D["db_isoform"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getDatabaseInfo(self):
        try:
            return self.__D["db_name"], self.__D["db_code"], self.__D["db_accession"], self.__D["db_isoform"]
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+ReferenceSequence.getDatabaseInfo() failed %r\n" % self.__D.items())
                traceback.print_exc(file=self.__lfh)
            return "", "", ""

    def __getStrain(self, sourceName=""):
        myRegex = r"\((.*)\)"
        try:
            # self.__lfh.write("+ReferenceSequence.__strainStrain() decoding %s\n" % sourceName)
            m = re.findall(myRegex, sourceName)
            if len(m) > 0:
                newSourceName = re.sub(myRegex, "", sourceName)
                if str(m[0]).startswith("strain"):
                    newStrain = str(m[0])[6:]
                else:
                    newStrain = str(m[0])

                if (len(newStrain) < 1) or (newStrain.upper() == "NONE"):
                    newStrain = ""
                return str(newSourceName).strip(), str(newStrain).strip()
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+ReferenceSequence.__getStrain() failed %s\n" % sourceName)
                traceback.print_exc(file=self.__lfh)

        return sourceName, ""

    def get(self):
        return self.__D

    def dump(self):
        for k, v in self.__D.items():
            self.__lfh.write("+ReferenceSequence.dump  key: %s value: %s\n" % (k, str(v)[:100]))


class ReferenceSequenceAssign(object):
    def __init__(self, verbose=True, log=sys.stderr):
        """Provides methods which encapsulate the details of reference sequence assignments."""
        self.__verbose = verbose
        self.__lfh = log
        self.__D = {}
        self.__attribStr = ["db_name", "db_accession", "db_code", "seq_one_letter_code", "full_seq_one_letter_code", "seq_id", "entity_id", "details"]
        self.__attribInt = ["ref_id", "db_align_end", "db_align_beg", "seq_align_begin", "seq_align_end", "auth_seq_align_begin", "auth_seq_align_end"]
        self.clear()

    def clear(self):
        self.__D = {}
        for attrib in self.__attribStr:
            self.__D[attrib] = ""
        for attrib in self.__attribInt:
            self.__D[attrib] = 0

    def set(self, seqD):
        try:
            for k, v in seqD.items():
                if (v is None) or (v == "") or (v in ["?", "."]):
                    continue
                if k in self.__attribStr:
                    self.__D[k] = str(v)
                elif k in self.__attribInt:
                    self.__D[k] = int(str(v))
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+RefenceSequenceAssign.set() - set failed on k %r v %r\n" % (k, v))
                traceback.print_exc(file=self.__lfh)
            return False
        #
        return True

    def get(self):
        return self.__D

    def getRefId(self):
        try:
            return self.__D["ref_id"]
        except:  # noqa: E722 pylint: disable=bare-except
            return 0

    def getEntityId(self):
        try:
            return self.__D["entity_id"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getDetails(self):
        try:
            return self.__D["details"]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def getSeqAlignRange(self):
        try:
            return self.__D["seq_align_begin"], self.__D["seq_align_end"]
        except:  # noqa: E722 pylint: disable=bare-except
            # traceback.print_exc(file=sys.stderr)
            return 0, 0

    def getDbAlignRange(self):
        try:
            return self.__D["db_align_beg"], self.__D["db_align_end"]
        except:  # noqa: E722 pylint: disable=bare-except
            return 0, 0

    def getDbReference(self):
        try:
            return self.__D["db_name"], self.__D["db_code"], self.__D["db_accession"]
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            return "", "", ""

    def getSequence(self):
        try:
            if "seq_one_letter_code" in self.__D:
                return self.__D["seq_one_letter_code"]
            elif "full_seq_one_letter_code" in self.__D:
                return self.__D["full_seq_one_letter_code"][self.__D["db_align_beg"] - 1 : self.__D["db_align_end"]]
            else:
                return ""
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def printIt(self, log=sys.stderr):
        log.write("\n+RefenceSequenceAssign.printIt() - entity Id    %s\n" % self.getEntityId())
        log.write("+RefenceSequenceAssign.printIt() - reference Id %s\n" % self.getRefId())
        log.write("+RefenceSequenceAssign.printIt() - entity sequence align range      %4d : %4d\n" % self.getSeqAlignRange())
        log.write("+RefenceSequenceAssign.printIt() - reference sequence align range   %4d : %4d\n" % self.getDbAlignRange())
        log.write("+RefenceSequenceAssign.printIt() - reference sequence name,code,acc %s : %s : %s \n" % self.getDbReference())
        log.write("+RefenceSequenceAssign.printIt() - reference sequence\n%s\n" % self.getSequence())
        log.write("+RefenceSequenceAssign.printIt() - reference details: %s\n" % self.getDetails())


class SequenceAssignDepositor(object):
    """Manages access to sequence assignment details provided during deposition.

    Storage model is dictionary of key value pairs -

    Content includes database identifiers, source details, and sequence comparison statistics.
    """

    def __init__(self, verbose=False, log=sys.stderr):
        self.__verbose = verbose
        self.__lfh = log
        self.__fD = {}
        self.__reset()

    def clear(self):
        self.__reset()

    def get(self):
        return self.__fD

    def set(self, assignD):
        self.__reset()
        self.__fD = assignD

    def __reset(self):
        self.__fD = {}
        #

    def getReferenceCount(self, entityId):
        count = 0
        try:
            return len(self.__fD[entityId]["ref_list"])
        except:  # noqa: E722 pylint: disable=bare-except
            pass
        return count

    def getReferenceList(self, entityId):
        """ """
        retL = []
        try:
            for d in self.__fD[entityId]["ref_list"]:
                sdm = ReferenceSequenceAssign(verbose=self.__verbose, log=self.__lfh)
                sdm.set(d)
                retL.append(sdm)
        except:  # noqa: E722 pylint: disable=bare-except
            self.__lfh.write("+SequenceAssignDepositor.getReferenceList() failing for entity %r\n" % entityId)
            traceback.print_exc(file=self.__lfh)
            return retL
        return retL

    def printIt(self, log=sys.stderr):
        """

        someList.sort(key=itemgetter('attrib_name'),reverse=True)
        """
        log.write("\n+SequenceAssignDepositor() - contents - \n")

        # for each entity
        #
        entityList = list(self.__fD.keys())
        if len(entityList) < 1:
            return

        entityList.sort()

        for eId in entityList:
            d = self.__fD[eId]
            if "ref_list" in d:
                log.write("+SequenceAssignDepositor.printIt() Depositor reference list for entity %r\n" % eId)
                for ii, rfD in enumerate(d["ref_list"]):
                    log.write(" +Reference %d\n" % (ii + 1))
                    for k, v in rfD.items():
                        log.write("     Key=%20s  : %s\n" % (k, v))

            if "ref_dif_dict" in d:
                log.write("\n+SequenceAssignDepositor.printIt() Depositor reference alignment difference dictionary for entity %r \n" % eId)
                rfD = d["ref_dif_dict"]
                for k, vL in rfD.items():
                    for ii, tD in enumerate(vL):
                        log.write(" +Sequence differences -------------------------  %d\n" % ii + 1)
                        for k, v in tD.items():
                            log.write("      Key=%30s  : %s\n" % (k, v))


class SequenceAssignArchive(object):
    """Manages access to archive sequence assignment details.

    Storage model is dictionary of key value pairs -

    Content include database identifiers, source details, and sequence comparison statistics.
    """

    def __init__(self, verbose=False, log=sys.stderr):
        self.__verbose = verbose
        self.__lfh = log
        self.__fD = {}
        self.__reset()
        #
        self.__attribMapB = {
            "pdbx_strand_id": "seq_id",
            "pdbx_db_accession": "db_accession",
            "align_id": "align_id",
            "ref_id": "ref_id",
            "pdbx_auth_seq_align_end": "auth_seq_align_end",
            "seq_align_beg": "seq_align_begin",
            "db_align_end": "db_align_end",
            "pdbx_auth_seq_align_beg": "auth_seq_align_beg",
            "pdbx_seq_align_beg_ins_code": "seq_align_beg_ins_code",
            "db_align_beg": "db_align_beg",
            "seq_align_end": "seq_align_end",
        }

        self.__attribMapA = {
            "id": "ref_id",
            "entity_id": "entity_id",
            "pdbx_db_accession": "db_accession",
            "pdbx_seq_one_letter_code": "seq_one_letter_code",
            "db_seq_one_letter_code": "seq_one_letter_code",
            "seq_one_letter_code": "seq_one_letter_code",
            "pdbx_align_begin": "seq_align_begin",
            "db_name": "db_name",
            "db_code": "db_code",
        }

    def clear(self):
        self.__reset()

    def get(self):
        return self.__fD

    def set(self, alignD):
        self.__reset()
        self.__fD = alignD

    def __reset(self):
        self.__fD = {}

    def getReferenceCount(self, entityId):
        count = 0
        try:
            return len(self.__fD[entityId]["ref_list"])
        except:  # noqa: E722 pylint: disable=bare-except
            pass
        return count

    def getReferenceList(self, entityId):
        """Get the list of sequence references assigned to this entity.

        Returns a list of
        """
        retL = []
        try:
            for td in self.__fD[entityId]["ref_list"]:
                #
                refD = {}
                for k, v in td.items():
                    if k in self.__attribMapA:
                        refD[self.__attribMapA[k]] = v

                #
                # Add the details from the first alignment if this exists.
                #
                refId = refD["ref_id"]
                if refId in self.__fD[entityId]["ref_align_dict"]:
                    aL = self.__fD[entityId]["ref_align_dict"][refId]
                    if len(aL) > 0:
                        aD = aL[0]
                        for k, v in aD.items():
                            if k in self.__attribMapB:
                                refD[self.__attribMapB[k]] = v
                else:
                    if self.__verbose:
                        self.__lfh.write("+SequenceAssignArchive.getReferenceList() no alignment data for ref_id %s\n" % refId)

                if self.__verbose:
                    self.__lfh.write("+SequenceAssignArchive.getReferenceList() refId %s  reference dictionary %r\n" % (refId, list(refD.items())))
                sdm = ReferenceSequenceAssign()
                sdm.set(refD)
                retL.append(sdm)

        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+SequenceAssignArchive - getReferenceList failing for entity %s\n" % entityId)
                traceback.print_exc(file=self.__lfh)

        return retL

    def printIt(self, log=sys.stderr):
        log.write("\n+SequenceAssignArchive() - Dictionary\n")

        # for each entity
        #
        entityList = list(self.__fD.keys())
        if len(entityList) < 1:
            return

        entityList.sort()

        for eId in entityList:
            d = self.__fD[eId]
            if "ref_list" in d:
                log.write("+SequenceAssignArchive.printIt() Archive reference list for entity %s\n" % eId)
                for _ii, rfD in enumerate(d["ref_list"]):
                    log.write(" +Reference %s\n" % rfD["id"])
                    for k, v in rfD.items():
                        log.write("     Key=%20s  : %s\n" % (k, v))

            if "ref_align_dict" in d:
                log.write("\n+SequenceAssignArchive.printIt() Archive reference alignment dictionary for entity %s \n" % eId)
                rfD = d["ref_align_dict"]
                for k, vL in rfD.items():
                    for tD in vL:
                        log.write(" +Alignment -------------------------  ref_id %s align_id %s\n" % (tD["ref_id"], tD["align_id"]))
                        for k, v in tD.items():
                            log.write("      Key=%30s  : %s\n" % (k, v))

            if "ref_dif_dict" in d:
                log.write("\n+SequenceAssignArchive.printIt() Archive reference alignment difference dictionary for entity %s \n" % eId)
                rfD = d["ref_dif_dict"]
                for k, vL in rfD.items():
                    for tD in vL:
                        log.write(" +Differences -------------------------  align_id %s\n" % tD["align_id"])
                        for k, v in tD.items():
                            log.write("      Key=%30s  : %s\n" % (k, v))
