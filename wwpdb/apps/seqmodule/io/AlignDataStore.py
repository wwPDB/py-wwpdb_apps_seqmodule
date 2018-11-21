##
# File:    AlignDataStore.py
# Date:    14-May-2014
#
# Updates:
#
##
"""
Provide a storage interface for recording sequence alignments.

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"


try:
    import cPickle as pickle
except ImportError:
    import pickle
import os
import os.path
import sys
import traceback

from wwpdb.apps.seqmodule.util.SequenceLabel import SequenceLabel


class AlignDataStore(object):
    """ Store and recover sequence alignment objects -

    Storage model is a list of lists where each list contains:

            alignSeqList [[Sequence identifier (ie. auth_A_1_1_1 as used in sequence data store),
                           SequenceLabel() object,
                           Aligned sequence with index details,
                           Conflict flag list (bool),
                           Feature dictionary for input sequence],,,]

    """

    def __init__(self, reqObj, fileName='alignDataStore.pic', verbose=False, log=sys.stderr):
        self.__verbose = verbose
        self.__lfh = log
        self.__reqObj = reqObj
        self.__sessionObj = self.__reqObj.getSessionObj()
        self.__fileName = fileName
        #
        self.__seqAlignList = []
        #

        self.__sessionPath = '.'
        self.__filePath = self.__fileName

        # self.__pickleProtocol = pickle.HIGHEST_PROTOCOL
        self.__pickleProtocol = 0

        self.__setup()
        #

    def __setup(self):
        try:
            self.__sessionPath = self.__sessionObj.getPath()
            self.__filePath = os.path.join(self.__sessionPath, self.__fileName)
            if (self.__verbose):
                self.__lfh.write("+AlignDataStore.__setup - using align store file path %s\n" % self.__filePath)

            if os.access(self.__filePath, os.R_OK):
                self.deserialize()
                if (self.__verbose):
                    self.__lfh.write("+AlignDataStore.__setup - opening align store with align list length %d\n" % len(
                        self.__seqAlignList))
            else:
                self.__lfh.write("+AlignDataStore.__setup - MISSING align store file path%s\n" % self.__filePath)

        except:
            self.__lfh.write(
                "+AlignDataStore.__setup - Failed opening align data store for session id %s  file path%s\n" %
                (self.__sessionObj.getId(), self.__filePath))
            traceback.print_exc(file=self.__lfh)

    def length(self):
        return len(self.__seqAlignList)

    def reset(self):
        self.__seqAlignList = []
        self.serialize()

    def serialize(self):
        try:
            oList = []
            fb = open(self.__filePath, 'wb')
            for ii in range(len(self.__seqAlignList)):
                t = self.__seqAlignList[ii][1].pack()
                self.__lfh.write("+AlignDataStore.serialize - saving sequence %d id %s\n" % (ii, t))
                oList.append([self.__seqAlignList[ii][0], t, self.__seqAlignList[ii][2], self.__seqAlignList[ii][3],
                              self.__seqAlignList[ii][4]])
            pickle.dump(oList, fb, self.__pickleProtocol)
            fb.close()
        except:
            pass

    def deserialize(self):
        try:
            fb = open(self.__filePath, 'rb')
            self.__seqAlignList = pickle.load(fb)
            for ii in range(len(self.__seqAlignList)):
                sLab = SequenceLabel()
                sLab.unpack(self.__seqAlignList[ii][1])
                self.__seqAlignList[ii][1] = sLab
            fb.close()
        except:
            pass

    def get(self):
        return self.__seqAlignList

    def set(self, seqAlignList):
        self.__seqAlignList = seqAlignList
        return True

    def getAlignIdList(self):
        try:
            return [a[0] for a in self.__seqAlignList]
        except:
            return []
