##
# File:    SequenceEditStore.py
# Date:    14-Jan-2010
#
# Updates:
# 26-Jan 2010 jdw Change pickel version - and storage model. Add list options.
# 14-Feb-2010 jdw Add delete and delete list operations / add filename to constructor
# 20-Apr-2010 jdw Ported to module seqmodule.
# 15-May-2010 jdw Add new element id to SequenceEdit object
# 27-Feb-2013 jdw Update documentation and call interface.
# 15-Sep-2017 zf  Move getEditStoreFilename() from AlignmentEdit.py to here
##
"""
Provide a storage interface for recording incremental edits in sequence data for
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
import os.path


def getEditStoreFilename(alignmentTag=""):
    """Add a tag to the edit store data filename to differentiate multiple active alignment views. -"""
    if len(alignmentTag) > 0:
        fn = "editDataStore-" + str(alignmentTag) + ".pic"
    else:
        fn = "editDataStore.pic"
    return fn


class SequenceEdit(object):
    """Container for the sequence data edit records --
    An edit record contains the following:
          Edit Operation id (ordering/grouping identifier),
          Edited element Id (e.g. residue label),
          Edit operation type,
          New value,
          Prior value
          ...
    """

    def __init__(self, verbose=False):  # pylint: disable=unused-argument
        # self.__verbose = verbose

        self.__editOpId = 0
        self.__editType = None
        self.__targetElementId = None
        self.__valueNew = None
        self.__valuePrevious = None
        self.__stylePrevious = None
        self.__newElementId = None

        self.__reset()

    def __reset(self):
        self.__editOpId = 0
        self.__editType = None
        self.__targetElementId = None
        self.__valueNew = None
        self.__valuePrevious = None
        self.__stylePrevious = None
        self.__newElementId = None

    def setValueNew(self, value):
        self.__valueNew = value

    def getValueNew(self):
        return self.__valueNew

    def setValuePrevious(self, value):
        self.__valuePrevious = value

    def getValuePrevious(self):
        return self.__valuePrevious

    def setStylePrevious(self, value):
        self.__stylePrevious = value

    def getStylePrevious(self):
        return self.__stylePrevious

    def setEditType(self, eType):
        self.__editType = eType

    def getEditType(self):
        return self.__editType

    def setTargetElementId(self, elementId):
        self.__targetElementId = elementId

    def getTargetElementId(self):
        return self.__targetElementId

    def setEditOpId(self, id):  # pylint:  disable=redefined-builtin
        self.__editOpId = id

    def getEditOpId(self):
        return self.__editOpId

    def setNewElementId(self, newElementId):
        self.__newElementId = newElementId

    def getNewElementId(self):
        return self.__newElementId

    def pack(self):
        return (self.__editOpId, self.__editType, self.__targetElementId, self.__valueNew, self.__valuePrevious, self.__stylePrevious, self.__newElementId)

    def unpack(self, editTuple):
        # self.__reset()
        # try:
        self.__editOpId = editTuple[0]
        self.__editType = editTuple[1]
        self.__targetElementId = editTuple[2]
        self.__valueNew = editTuple[3]
        self.__valuePrevious = editTuple[4]
        self.__stylePrevious = editTuple[5]
        self.__newElementId = editTuple[6]

    # except:
    #    pass

    def printIt(self, ofh):
        ofh.write("\nSequence Edit Object Contents:\n")
        ofh.write("  Edit operation Id %s\n" % self.__editOpId)
        ofh.write("  Edit type         %s\n" % self.__editType)
        ofh.write("  Target element Id %r\n" % self.__targetElementId)
        ofh.write("  New value         %s\n" % self.__valueNew)
        ofh.write("  Previous value    %s\n" % self.__valuePrevious)
        ofh.write("  Previous style    %s\n" % self.__stylePrevious)
        ofh.write("  New element Id    %s\n" % self.__newElementId)


class SequenceEditStore(object):
    """Store incremental edits on the sequences -

        An edit record contains the following:
              Edit Operation id (ordering/grouping identifier),
              Edited element Id (e.g. residue label),
              Edit operation type (replace,insert,delete,or details),
              New value,
              Prior value
              ...
    Storage model is a list of tuples where each tuple contains the above attributes.


    """

    def __init__(self, sessionObj, fileName="sequenceEditStore.pic", verbose=False, log=sys.stderr):
        self.__verbose = verbose
        self.__lfh = log
        self.__sessionObj = sessionObj
        self.__fileName = fileName
        #
        # List of SequenceEdit objects -
        self.__editList = []
        #

        self.__sessionPath = "."
        self.__filePath = self.__fileName

        # self.__pickleProtocol = pickle.HIGHEST_PROTOCOL
        self.__pickleProtocol = 0

        self.__setup()
        #

    def __setup(self):
        try:
            self.__sessionPath = self.__sessionObj.getPath()
            self.__filePath = os.path.join(self.__sessionPath, self.__fileName)
            if self.__verbose:
                self.__lfh.write("+SequenceEditStore.__setup - session id %s  edit store file path%s\n" % (self.__sessionObj.getId(), self.__filePath))

            self.deserialize()
            if self.__verbose:
                self.__lfh.write("+SequenceEditStore.__setup - opening edit store and the following recovering edit list length %d\n" % len(self.__editList))
                self.printIt(self.__lfh)

        except:  # noqa: E722 pylint: disable=bare-except
            self.__lfh.write("+SequenceEditStore.__setup - Failed opening edit store for session id %s  edit store file path%s\n" % (self.__sessionObj.getId(), self.__filePath))

    def length(self):
        return len(self.__editList)

    def reset(self):
        self.__editList = []
        self.serialize()

    def serialize(self):
        try:
            fb = open(self.__filePath, "wb")
            pickle.dump(self.__editList, fb, self.__pickleProtocol)
            fb.close()
        except:  # noqa: E722 pylint: disable=bare-except
            pass

    def deserialize(self):
        try:
            fb = open(self.__filePath, "rb")
            self.__editList = pickle.load(fb)
            fb.close()
            # for tup in tupList:
            #    sEd=SequenceEdit()
            #    sEd.unpack(tup)
            #    self.__editList.append(sEd)
        except:  # noqa: E722 pylint: disable=bare-except
            pass

    def storeEdit(self, seqEdObj):
        try:
            eTup = seqEdObj.pack()
            self.__editList.append(eTup)
            self.serialize()
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            return False

    def deleteEdit(self, seqEdObj):
        try:
            eTup = seqEdObj.pack()
            if eTup in self.__editList:
                self.__editList.remove(eTup)
            self.serialize()
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            return False

    def storeEditList(self, seqEdObjList):
        try:
            for seqEdObj in seqEdObjList:
                eTup = seqEdObj.pack()
                self.__editList.append(eTup)
            self.serialize()
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            return False

    def deleteEditList(self, seqEdObjList):
        try:
            for seqEdObj in seqEdObjList:
                eTup = seqEdObj.pack()
                if eTup in self.__editList:
                    self.__editList.remove(eTup)
            self.serialize()
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            return False

    def printIt(self, ofh):
        ofh.write("\nSequence Edit Store Contents\n")
        ofh.write("  Storage path:  %s\n" % self.__filePath)
        ofh.write("  Edit count:    %d\n" % len(self.__editList))
        for eTup in self.__editList:
            sEd = SequenceEdit()
            sEd.unpack(eTup)
            sEd.printIt(ofh)

    def getLastEditOp(self):
        try:
            eTup = self.__editList[-1]
            sEd = SequenceEdit()
            sEd.unpack(eTup)
            return sEd.getEditOpId()
        except:  # noqa: E722 pylint: disable=bare-except
            return 0

    #
    def get(self, opId):
        """Return the list of edit objects corresponding to the edit operation opId."""
        oL = []
        for eTup in self.__editList:
            sEd = SequenceEdit()
            sEd.unpack(eTup)
            if opId == sEd.getEditOpId():
                oL.append(sEd)
        return oL

    def getList(self):
        """Return the list of edit objects."""
        oL = []
        for eTup in self.__editList:
            sEd = SequenceEdit()
            sEd.unpack(eTup)
            oL.append(sEd)
        return oL

    def remove(self, opId):
        """Remove edit objects with the input edit operation Id from the current store
        and save the result.
        """
        newList = []
        for eTup in self.__editList:
            sEd = SequenceEdit()
            sEd.unpack(eTup)
            if opId != sEd.getEditOpId():
                newList.append(eTup)
        self.__editList = newList
        self.serialize()

    def getDetailsByTarget(self):
        dD = {}
        for eTup in self.__editList:
            sEd = SequenceEdit()
            sEd.unpack(eTup)
            if sEd.getEditType() == "details":
                dD[sEd.getTargetElementId()] = sEd.getValueNew()
        return dD
