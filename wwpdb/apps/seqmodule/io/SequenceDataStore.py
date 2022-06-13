##
# File:    SequenceDataStore.py
# Date:    14-Jan-2010
#
# Updates:
#   20-Apr-2010 jdw Ported to module seqmodule.
#   04-May-2010 jdw Add methods to capture sequence selections.
#   25-Feb-2013 jdw update documentation
#   27-Feb-2013 jdw replace autoDict with expanded external class Autodict()
#   27-Feb-2013 jdw add support for polymer entity parts -
#                   consolidate alignment data within this class
#   03-Mar-2013 jdw add group/sub-part support
#   06-Mar-2013 jdw change default to use common tool file name conventions
#                   via PathInfo()
#   07-Mar-2013 jdw Use fileName in constructor to override default name convention.
#   18-Mar-2013 jdw Add storage of depositor and archive assignment data.
#   19-Mar-2013 jdw reverse the order of alternative ids to reflect changed
#                   organization of reference sequence list.
#   20-Mar-2012 jdw add methods getFeatureObj() and setFeatureObj()
#   27-Mar-2013 jdw add filterIndex method which addresses the problem of removing
#                   index entries from the Autodict() data structure.
#    1-Dec-2013 jdw add support for storing polymer linkage distances -
#    4-Dec-2013 jdw add method to obtain groupId for an instance seqId
#   19-Dec-2013 jdw add getPartIdsForVersion(self, seqId, dataType="sequence", seqType="ref",altId=1,version=1)
#    3-Sep-2017 zf  add removePartId()
##
"""
Provide a storage interface for sequence and sequence feature data for
use by the sequence editing tool.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.07"

try:
    import cPickle as pickle
except ImportError:
    import pickle as pickle
#

import sys
import time
import pprint
import traceback

from wwpdb.apps.seqmodule.util.Autodict import Autodict
from wwpdb.io.locator.PathInfo import PathInfo
from wwpdb.apps.seqmodule.util.SequenceLabel import SequenceFeature


class SequenceDataStore(object):

    """All sequences are stored as lists and features are stored in dictionaries -
    - These objects are opaque to this class which only manages their storage.

    Objects are identified by  -

      dataType:  Supported data types are 'sequence', 'feature', 'link-distance'
      seqType:   Sequence types are: auth (author provided), xyz (coordinate),
                 ref (from reference sequence database).
      seqId:     Identifies an instance of a sequence or a feature object within a
                 sequence type.  Sequence and feature instances are versioned.  This
                 identifier serves can be a PDB ChainId or the PDBx (_struct_asym.id)
                 for the author and coordinate sequence types.
      partId:    Distinguishes continuous regions of a sequence with unique/separately
                 specified features. This identifier is used to maintain correpondences
                 between portions of an author provided sequence and related reference
                 sequences.    This is an integer value with default value 1.
      altId:     Distinguishes alternatives associated with a particular sequence
                 instance.   This may be used to support micro-hetereogeneity and
                 alternative reference sequences (default=1)
      version:   integer revision id where 1 is the original version.

      self.__I[dataType][seqType][seqId][partId][altId][version]=((dataId,timeStamp)

      self.__D[dataId] = sequence list [...] or  feature dictionary {...}

      where dataId is a concatenated identifier= seqType_dataType_seqId_partId_altId_version

      Sequence instance groups  (e.g. entity->instance,... ) identify groups of related
      sequence instances:

         Group Id/entity Id -> [seqId1, seqId2, ...]

         self.__G[groupId] = [seqId List]

      Group sub-subsequences (parts):

         self.__P[groupId/entityId] = [partId1, partId2,...]


      Sequence id selections:

      self.__S[]= [seqId label list]


      Sequence Id selections:

      self.__S[]= [seqId label list]

      Sequence Alignment Id lists:

      self.__L[]= [seqId label list]

      Entry features are stored in a dictoinary of key value pairs. The provides a
      container for storing features of the entry -

      self.__E={}

      Depositor sequence assignment and conflict annotation -

      self.__depositorAssignD={<entityId>}

      Archive sequence assignment and conflict annotation -

      self.__assignD={<entityId>}

    """

    def __init__(self, reqObj=None, fileName=None, verbose=False, log=sys.stderr):
        self.__verbose = verbose
        self.__debug = False
        self.__lfh = log
        self.__reqObj = reqObj
        self.__sessionObj = self.__reqObj.getSessionObj()
        self.__sessionPath = self.__sessionObj.getPath()
        #
        # storeType =  sequence|alignment
        #
        self.__fileName = fileName
        self.__filePath = None
        #
        self.__clear()
        #
        # self.__pickleProtocol = pickle.HIGHEST_PROTOCOL
        self.__pickleProtocol = 0
        #
        # Declar for pylint
        self.__D = {}
        self.__I = Autodict()
        self.__E = {}
        self.__G = {}
        self.__P = {}
        self.__S = []
        self.__L = []
        self.__depositorAssignD = {}
        self.__assignD = {}
        #
        self.__setup()

    def __clear(self):
        """ """
        self.__D = {}
        self.__I = Autodict()
        self.__E = {}
        self.__G = {}
        self.__P = {}
        self.__S = []
        self.__L = []
        self.__depositorAssignD = {}
        self.__assignD = {}

    def __setup(self):

        try:
            self.__siteId = self.__reqObj.getValue("WWPDB_SITE_ID")
            self.__identifier = self.__reqObj.getValue("identifier")
            self.__pI = PathInfo(siteId=self.__siteId, sessionPath=self.__sessionPath, verbose=self.__verbose, log=self.__lfh)
            if self.__fileName is not None:
                # self.__filePath = os.path.join(self.__sessionPath, self.__fileName)
                # Using full path file name instead
                self.__filePath = self.__fileName
            else:
                self.__filePath = self.__pI.getSequenceStatsFilePath(self.__identifier, fileSource="session")

            if self.__verbose:
                self.__lfh.write("+SequenceDataStore.__setup() - Starting with session id %s \n" % self.__sessionObj.getId())
                self.__lfh.write("+SequenceDataStore.__setup() - using data store path %s\n" % self.__filePath)

            self.deserialize()
        except:  # noqa: E722 pylint: disable=bare-except
            self.__lfh.write("+SequenceDataStore.__setup() - Failed to open data store for session id %s data store path %s\n" % (self.__sessionObj.getId(), self.__filePath))
            traceback.print_exc(file=self.__lfh)

    def reset(self):
        self.__D = {}
        self.__I = Autodict()
        self.__G = {}

    def getFilePath(self):
        return self.__filePath

    def serialize(self):
        try:
            fb = open(self.__filePath, "wb")
            pickle.dump(self.__E, fb, self.__pickleProtocol)
            pickle.dump(self.__G, fb, self.__pickleProtocol)
            pickle.dump(self.__P, fb, self.__pickleProtocol)
            pickle.dump(self.__I, fb, self.__pickleProtocol)
            pickle.dump(self.__D, fb, self.__pickleProtocol)
            pickle.dump(self.__S, fb, self.__pickleProtocol)
            pickle.dump(self.__L, fb, self.__pickleProtocol)
            pickle.dump(self.__depositorAssignD, fb, self.__pickleProtocol)
            pickle.dump(self.__assignD, fb, self.__pickleProtocol)
            fb.close()
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+SequenceDataStore.__serialize() - failing for %s\n" % self.__filePath)
                traceback.print_exc(file=self.__lfh)
            #
        #

    def deserialize(self):
        try:
            fb = open(self.__filePath, "rb")
            self.__E = pickle.load(fb)
            self.__G = pickle.load(fb)
            self.__P = pickle.load(fb)
            self.__I = pickle.load(fb)
            self.__D = pickle.load(fb)
            self.__S = pickle.load(fb)
            self.__L = pickle.load(fb)
            self.__depositorAssignD = pickle.load(fb)
            self.__assignD = pickle.load(fb)
            fb.close()
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__debug:
                self.__lfh.write("+SequenceDataStore.__deserialize() - failing for %s\n" % self.__filePath)
                traceback.print_exc(file=self.__lfh)
            #
            self.__clear()
        #

    def __makeId(self, dataType, seqType, seqId, partId=1, altId=1, version=1):
        rid = "%s_%s_%s_%d_%d_%d" % (dataType, seqType, seqId, int(partId), int(altId), int(version))
        return rid

    def __updateIndex(self, tid, dataType, seqType, seqId, partId=1, altId=1, version=1):
        lt = time.strftime("%Y %m %d %H:%M:%S", time.localtime())

        self.__I[dataType][seqType][seqId][int(partId)][int(altId)][int(version)] = (tid, lt)

    def setSequence(self, sL, seqId, seqType, partId=1, altId=1, version=1):
        try:
            rid = self.__makeId(dataType="sequence", seqType=seqType, seqId=seqId, partId=partId, altId=altId, version=version)
            self.__updateIndex(rid, dataType="sequence", seqType=seqType, seqId=seqId, partId=partId, altId=altId, version=version)
            if sL is None or len(sL) == 0:
                if self.__verbose:
                    self.__lfh.write(
                        "+SequenceDataStore.__setSequence() - empty object seqId %s seqType %s partId %d altId %d version %d\n" % (seqId, seqType, partId, altId, version)
                    )
            self.__D[rid] = sL
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            return False

    def getSequence(self, seqId, seqType, partId=1, altId=1, version=1):
        try:
            sid = self.__makeId(dataType="sequence", seqType=seqType, seqId=seqId, partId=partId, altId=altId, version=version)
            return self.__D[sid]
        except:  # noqa: E722 pylint: disable=bare-except
            return []

    def dumpSequence(self, seqId, seqType, partId=1, altId=1, version=1):
        try:
            sid = self.__makeId(dataType="sequence", seqType=seqType, seqId=seqId, partId=partId, altId=altId, version=version)
            sTupL = self.__D[sid]
            self.__lfh.write("+SequenceDataStore.dumpSequence() - seqId %s seqType %s partId %d altId %d version %d\n" % (seqId, seqType, partId, altId, version))
            for sTup in sTupL:
                self.__lfh.write("         ++ %5s %6s %d (%s)\n" % (sTup[0], sTup[1], sTup[3], sTup[2]))
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            return False

    def setFeature(self, fD, seqId, seqType, partId=1, altId=1, version=1):
        try:
            fid = self.__makeId(dataType="feature", seqType=seqType, seqId=seqId, partId=partId, altId=altId, version=version)
            self.__updateIndex(fid, dataType="feature", seqType=seqType, seqId=seqId, partId=partId, altId=altId, version=version)
            self.__D[fid] = fD
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            return False

    def setFeatureObj(self, sfObj, seqId, seqType, partId=1, altId=1, version=1):
        try:
            idx = self.__makeId(dataType="feature", seqType=seqType, seqId=seqId, partId=partId, altId=altId, version=version)
            self.__updateIndex(idx, dataType="feature", seqType=seqType, seqId=seqId, partId=partId, altId=altId, version=version)
            self.__D[idx] = sfObj.get()
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            self.__lfh.write("+SequenceDataStore.setFeatureObj() - failing for seqId %s seqType %s\n" % (seqId, seqType))
            traceback.print_exc(file=self.__lfh)
            return False

    def getFeature(self, seqId, seqType, partId=1, altId=1, version=1):
        try:
            fid = self.__makeId(dataType="feature", seqType=seqType, seqId=seqId, partId=partId, altId=altId, version=version)
            return self.__D[fid]
        except:  # noqa: E722 pylint: disable=bare-except
            return {}

    def getFeatureObj(self, seqId, seqType, partId=1, altId=1, version=1):
        sf = SequenceFeature()
        try:
            fid = self.__makeId(dataType="feature", seqType=seqType, seqId=seqId, partId=partId, altId=altId, version=version)
            sf.set(self.__D[fid])
            return sf
        except:  # noqa: E722 pylint: disable=bare-except
            return sf

    def getSequenceTypes(self, dataType="sequence"):
        try:
            return list(self.__I[dataType].keys())
        except:  # noqa: E722 pylint: disable=bare-except
            return []

    def getSelectedIds(self):
        try:
            return self.__S
        except:  # noqa: E722 pylint: disable=bare-except
            return []

    def addSelectedId(self, id):  # pylint: disable=redefined-builtin
        try:
            self.__S.append(id)
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            return False

    def setSelectedIds(self, idList=None):
        if idList is None:
            idList = []
        try:
            self.__S = []
            self.__S.extend(idList)
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+SequenceDataStore.setSelectedIds() - failing for %r\n" % idList)
                traceback.print_exc(file=self.__lfh)
            return False

    def getDataTypes(self):
        try:
            return list(self.__I.keys())
        except:  # noqa: E722 pylint: disable=bare-except
            return []

    def getIds(self, dataType="sequence", seqType="ref"):
        try:
            return list(self.__I[dataType][seqType].keys())
        except:  # noqa: E722 pylint: disable=bare-except
            return []

    def getPartIds(self, seqId, dataType="sequence", seqType="ref"):
        """Return the list of part ids in "ascending" order (smallest id  first)."""
        try:
            pL = sorted(self.__I[dataType][seqType][seqId].keys())
            return pL
        except:  # noqa: E722 pylint: disable=bare-except
            return []

    def getPartIdsForVersion(self, seqId, dataType="sequence", seqType="ref", altId=1, version=1):
        """Return the list of part ids in "ascending" order (smallest id  first).

        return parts for the input version only.
        """
        try:
            pL = []
            tpL = list(self.__I[dataType][seqType][seqId].keys())
            for p in tpL:
                if int(version) in self.__I[dataType][seqType][seqId][int(p)][int(altId)]:
                    pL.append(p)
            pL.sort()
            return pL
        except:  # noqa: E722 pylint: disable=bare-except
            return []

    def removePartId(self, seqId, dataType="sequence", seqType="auth", partId=2):
        """ """
        try:
            if int(partId) in self.__I[dataType][seqType][seqId]:
                del self.__I[dataType][seqType][seqId][int(partId)]
                return True
            #
        except:  # noqa: E722 pylint: disable=bare-except
            return False
        #
        return False

    def clearIndex(self, seqId, dataType="sequence", seqType="ref"):
        """Clear index below the sequence id."""
        try:
            self.__I[dataType][seqType][seqId] = {}
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            return False

    def getAlternativeIds(self, seqId, dataType="sequence", seqType="ref", partId=1):
        """Return the list of alternative ids in "descending" order (largest id  first)."""
        try:
            aL = list(self.__I[dataType][seqType][seqId][int(partId)].keys())
            aL.sort(reverse=True)
            return aL
        except:  # noqa: E722 pylint: disable=bare-except
            return []

    def getVersionIds(self, seqId, partId=1, altId=1, dataType="sequence", seqType="ref"):
        """Return the list of version ids in "descending" order (highest version first)."""
        try:
            vers = list(self.__I[dataType][seqType][seqId][int(partId)][int(altId)].keys())
            vers.sort(reverse=True)
            return vers
        except:  # noqa: E722 pylint: disable=bare-except
            return []

    def setGroup(self, groupId, seqIdList=None):
        if seqIdList is None:
            seqIdList = []
        try:
            if seqIdList is not None:
                self.__G[groupId] = seqIdList
            else:
                self.__G[groupId] = []
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            return False

    def getGroup(self, groupId):
        try:
            return self.__G[groupId]
        except:  # noqa: E722 pylint: disable=bare-except
            return []

    def getGroupId(self, seqId):
        """Return the entity/groupID corresponding to the input instance seqId. or None."""
        for groupId, idList in self.__G.items():
            if seqId in idList:
                return groupId
        return None

    def getGroupIds(self):
        try:
            keys = list(self.__G.keys())
            keys.sort(key=int)
            return keys
        except:  # noqa: E722 pylint: disable=bare-except
            return []

    def setEntryDetail(self, detailKey, detailValue):
        try:
            self.__E[detailKey] = detailValue
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            return False

    def getEntryDetail(self, detailKey):
        try:
            return self.__E[detailKey]
        except:  # noqa: E722 pylint: disable=bare-except
            return ""

    def setAlignIdList(self, alignIdList):
        """Set the list of identifiers for aligned sequences."""
        self.__L = alignIdList

    def getAlignIdList(self):
        """Return the list of identifiers for the stored aligned sequences."""
        return self.__L

    def setGroupParts(self, groupId, partIdList=None):
        if partIdList is None:
            partIdList = []
        try:
            self.__P[groupId] = partIdList
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            return False

    def getGroupParts(self, groupId):
        try:
            return self.__P[groupId]
        except:  # noqa: E722 pylint: disable=bare-except
            return []

    def getGroupPartCount(self, groupId):
        try:
            return len(self.__P[groupId])
        except:  # noqa: E722 pylint: disable=bare-except
            return 0

    def setDepositorReferenceAssignments(self, assignD=None):
        try:
            self.__depositorAssignD = assignD if assignD is not None else {}
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            return False

    def setReferenceAssignments(self, assignD=None):
        try:
            self.__assignD = assignD if assignD is not None else {}
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            return False

    def getDepositorReferenceAssignments(self):
        try:
            return self.__depositorAssignD
        except:  # noqa: E722 pylint: disable=bare-except
            return {}

    def getReferenceAssignments(self):
        try:
            return self.__assignD
        except:  # noqa: E722 pylint: disable=bare-except
            return {}

    def filterIndex(self, seqId, dataType="sequence", seqType="ref"):
        """Copy indexes -

        self.__I[dataType][seqType][seqId][partId][altId][version]=((dataId,timeStamp)

        """
        J = Autodict()
        for dataType0, v0 in self.__I.items():
            for seqType0, v1 in v0.items():
                for seqId0, v2 in v1.items():
                    for partId0, v3 in v2.items():
                        for altId0, v4 in v3.items():
                            for verId0, vtup in v4.items():
                                if (dataType == dataType0) and (seqType == seqType0) and (seqId == seqId0):
                                    continue
                                else:
                                    J[dataType0][seqType0][seqId0][partId0][altId0][verId0] = vtup

        if self.__debug:
            self.__lfh.write("\n  + %s index:\n" % dataType)
            for vtype, v0 in J[dataType].items():
                for vid, v1a in v0.items():
                    for pId, v1b in v1a.items():
                        for altId, v2 in v1b.items():
                            for ver, ival in v2.items():
                                self.__lfh.write(
                                    "   type %5s id %4s partId %d verId %2s altId %4d updated %12s seq len %10d\n" % (vtype, vid, pId, ver, altId, ival[1], len(self.__D[ival[0]]))
                                )
        self.__I = J

    def dump(self, ofh):
        """Dump indexes -

        self.__I[dataType][seqType][seqId][partId][altId][version]=((dataId,timeStamp)

        """
        ofh.write("\n+BEGIN++BEGIN++BEGIN++BEGIN++BEGIN+  SequenceDataStore.dump() +BEGIN++BEGIN++BEGIN++BEGIN++BEGIN++BEGIN++BEGIN+\n")
        ofh.write("+Contents of sequence data store: %s \n" % self.__filePath)
        #
        #
        nItems = 0
        for _type, v0 in self.__I["sequence"].items():
            for _id, v1a in v0.items():
                for pId, v1b in v1a.items():
                    for altId, v2 in v1b.items():
                        for ver, ival in v2.items():
                            nItems += 1

        ofh.write("\n  +Sequence Index contains %d items\n" % nItems)

        if self.__debug:
            ofh.write("\n  +Sequence Index:\n")
            for vtype, v0 in self.__I["sequence"].items():
                for vid, v1a in v0.items():
                    for pId, v1b in v1a.items():
                        for altId, v2 in v1b.items():
                            for ver, ival in v2.items():
                                ofh.write(
                                    "   type %5s id %4s partId %d verId %2s altId %4d updated %12s seq len %10d\n" % (vtype, vid, pId, ver, altId, ival[1], len(self.__D[ival[0]]))
                                )
        nItems = 0
        for _type, v0 in self.__I["feature"].items():
            for _id, v1a in v0.items():
                for pId, v1b in v1a.items():
                    for altId, v2 in v1b.items():
                        for ver, ival in v2.items():
                            nItems += 1

        ofh.write("\n  +Feature Index contains %d items\n" % nItems)

        if self.__debug:
            ofh.write("\n  +Feature Index:\n")
            for vtype, v0 in self.__I["feature"].items():
                for vid, v1a in v0.items():
                    for pId, v1b in v1a.items():
                        for altId, v2 in v1b.items():
                            for ver, ival in v2.items():
                                ofh.write(
                                    "   type %5s id %4s partID %2d verId %2s altId %4d updated %12s feat len %10d\n"
                                    % (vtype, vid, pId, ver, altId, ival[1], len(self.__D[ival[0]]))
                                )

        ofh.write("\n  +Sequence Group/Entity Index:\n")
        for gId, seqIdList in self.__G.items():
            ofh.write("  +Sequence group:  %5s\n" % gId)
            for sid in seqIdList:
                ofh.write("   %2s" % sid)
            ofh.write("\n")

        ofh.write("\n  +Sequence Group/Entity Subpart Index:\n")
        for gId, partIdList in self.__P.items():
            ofh.write("   +Sequence group:  %5s\n" % gId)
            for pid in partIdList:
                ofh.write("   %2s" % pid)
            ofh.write("\n")

        ofh.write("\n  +Entry Details:\n")
        for k, v in self.__E.items():
            ofh.write("    +Key:  %-35s value: %s\n" % (k, v))

        ofh.write("\n  +Selected sequences:\n")
        for k in self.__S:
            ofh.write("    +Sequence ID:  %s\n" % k)
        #
        ofh.write("\n  +Aligned sequences:\n")
        for k in self.__L:
            ofh.write("    +Sequence ID:  %s\n" % k)

        ofh.write("\n  +Depositor reference sequences assignments:\n")
        for eId, dd in self.__depositorAssignD.items():
            ofh.write("    +Sequence Group/Entity ID:  %s\n" % eId)
            for k, v in dd.items():
                ofh.write("    +Key:  %-35s length: %d\n" % (k, len(v)))

        ofh.write("\n  +Archive reference sequences assignments:\n")
        for eId, dd in self.__assignD.items():
            ofh.write("    +Sequence Group/Entity ID:  %s\n" % eId)
            for k, v in dd.items():
                ofh.write("    +Key:  %-35s length: %d\n" % (k, len(v)))
        ofh.write("\n+END++END++END++END++END++END++END++END++END++END+ SequenceDataStore.dump() +END++END++END++END++END++END++END++END+\n\n")

    def dumpData(self, ofh, seqId, seqType, dataType="sequence", partId=1):
        for altId, vOb in self.__I[dataType][seqType][seqId][partId].items():
            for ver, _ival in vOb.items():
                did = self.__makeId(dataType=dataType, seqType=seqType, seqId=seqId, partId=partId, altId=altId, version=ver)
                ofh.write("Data contents for %10s type %5s id %5s part %d alternative %2d version %2d\n" % (dataType, seqType, seqId, partId, altId, ver))
                if did in self.__D:
                    pprint.pprint(self.__D[did], stream=ofh)
                else:
                    ofh.write(" NO DATA FOUND\n")
