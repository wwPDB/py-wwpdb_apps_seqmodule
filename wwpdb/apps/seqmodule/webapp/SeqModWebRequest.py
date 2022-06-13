##
# File:    SeqModWebRequest.py
# Date:    20-Feb-2013
#          Derived for common WebRequest().
# Updated:
# 8-Jul-2014  jdw add groupids response object -
# 2-Dec-2014  jdw capture updates in selection and alignment id lists after alignment edits -- rev-allalignids,rev-selectids
# 30=Aug-2017 zf  change the behavior of getSummaryAlignList() & getSummarySelectList() related to rev-allalignids & rev-selectids.
#                 add setWarning()
##
"""
SeqModWebRequest provides containers and accessors for managing request parameter information
with extensions for the Sequence module.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.07"

import sys
from wwpdb.utils.session.WebRequest import InputRequest, ResponseContent
from wwpdb.utils.session.UtilDataStore import UtilDataStore


class SeqModInputRequest(InputRequest):
    def __init__(self, paramDict, verbose=False, log=sys.stderr):
        super(SeqModInputRequest, self).__init__(paramDict, verbose=verbose, log=log)
        self.__verbose = verbose
        self.__lfh = log

    def setNewRefId(self, refId):
        """Save in local request and persistent session store -"""
        self.setValue(refId, "new-ref-seq-id")
        uds = UtilDataStore(reqObj=self, verbose=self.__verbose, log=self.__lfh)
        uds.set("new-ref-seq-id", refId)
        uds.serialize()

    def __saveDataStore(self, key, value):
        """Save input key-value pair in persistent session store -"""
        uds = UtilDataStore(reqObj=self, verbose=self.__verbose, log=self.__lfh)
        uds.set(key, value)
        uds.serialize()

    def getDeleteList(self):
        """Return a list of tuples containing [(elementId,previousValue),(),...]"""
        oL = []
        try:
            dString = self._getStringValue("deleteselect")
            # self.__lfh.write("\n+SequenceInputRequest.getDeleteList() dString length %d : %r\n" % (len(dString),dString))
            dList = dString.split(",")
            if self.__verbose:
                self.__lfh.write("\n+SequenceInputRequest.getDeleteList() dList length %d : %r\n" % (len(dList), dList))
            for d in dList:
                # self.__lfh.write("\n+SequenceInputRequest.getDeleteList() dList %s\n" % d)
                dTup = d.split("|")
                # self.__lfh.write("\n+SequenceInputRequest.getDeleteList() dList %s %s\n" % (dTup[0],dTup[1]))
                oL.append((dTup[0].strip(), dTup[1].strip()))
        except:  # noqa: E722 pylint: disable=bare-except
            self.__lfh.write("\n+SequenceInputRequest.getDeleteList() error processing delete select list\n")
        return oL

    def getEditMenuSelectList(self):
        """Process the global edit menu select list --

        'targetId|priorVal|newValue|priorCss', 'targetId|priorVal|newValue|priorCss', 'targetId|priorVal|newValue|priorCss', ...,

        Return a list of tuples -

        [(targetId,priorValue,NewValue,priorCss),(,,),,,]

        """
        oL = []
        #        try:
        selectString = self._getStringValue("selectedids")
        rowList = selectString.split(",")
        if self.__verbose:
            self.__lfh.write("\n+SequenceInputRequest.getEditMenuSelectList() rowlist length %d : %r\n" % (len(rowList), rowList))
        for row in rowList:
            tp = row.split("|")
            oL.append((tp[0], tp[1], tp[2], tp[3]))

            #        except:
            #            self.__lfh.write("\n+SequenceInputRequest.getEditMenuSelectList(self): failed processing select list\n")
        return oL

    def getSelectList(self):
        """Return a list of tuples containing [(id,value,type),(id,value,type),(id,value),(),...]

        3 id,value pairs are returned per conflict record in the order ref, aligned, details.
        """
        oL = []
        colLab = ["reference", "aligned", "details"]
        try:
            selectString = self._getStringValue("selectedids")
            rowList = selectString.split(",")
            if self.__verbose:
                self.__lfh.write("\n+SequenceInputRequest.getSelectList() rowlist length %d : %r\n" % (len(rowList), rowList))
            for row in rowList:
                idTupList = row.split("|")
                #
                # self.__lfh.write("row = %r\n" % row)
                # self.__lfh.write("idtuplist = %r\n" % idTupList)
                # if (self.__verbose):
                #    self.__lfh.write("\n+SequenceInputRequest.getSelectList() idTupList length %d : %r\n" % (len(idTupList),idTupList))

                ii = 0
                for idTup in idTupList:
                    sTup = idTup.split("~")
                    oL.append((sTup[0], sTup[1], colLab[ii]))
                    ii += 1
                    # self.__lfh.write("idTup %r\n" % idTup)
                    # self.__lfh.write("sTup %r\n" % sTup)
                    # self.__lfh.write("length sTup %d\n" % len(sTup))
                    # self.__lfh.write("\n+SequenceInputRequest.getSelectList() id %s value %s type %s\n" % (sTup[0],sTup[1],"static"))
                    # self.__lfh.write(" sTup[0] = %s\n" % sTup[0])
                    # self.__lfh.write(" sTup[1] = %s\n" % sTup[1])
                    # self.__lfh.write(" ii = %d\n" % ii)
                    # self.__lfh.write(" colLab  = %s\n" % colLab[ii])

        except:  # noqa: E722 pylint: disable=bare-except
            self.__lfh.write("\n+SequenceInputRequest.getSelectList() error processing select list\n")

        return oL

    def getAlignmentOp(self):
        tV = self._getStringValue("operation")
        if tV in ["re-align", "reload", "re-load"]:
            return "realign"
        elif tV in ["reset"]:
            return "reset"
        elif tV in ["load"]:
            return "load"
        elif tV in ["loadandstore"]:
            return "loadandstore"
        else:
            return "load"

    def getAlignmentOrdering(self):
        """RPS, 20100727: Return ordering code as sourced from UI drop-down box via "viewalign_ordering" parameter."""
        tV = self._getStringValue("viewalign_order")
        # self.__lfh.write("\n+WebRequest.getAlignmentOrdering() got viewalign_order as: %s\n" % tV )

        if tV in ["auth-xyz-ref"]:
            return "auth-xyz-ref"
        elif tV in ["auth-ref-xyz"]:
            return "auth-ref-xyz"
        elif tV in ["xyz-auth-ref"]:
            return "xyz-auth-ref"
        elif tV in ["xyz-ref-auth"]:
            return "xyz-ref-auth"
        elif tV in ["ref-auth-xyz"]:
            return "ref-auth-xyz"
        elif tV in ["ref-xyz-auth"]:
            return "ref-xyz-auth"
        else:
            # self.__lfh.write("\n+WebRequest.getAlignmentOrdering() got nothing from front-end!!" )
            return "auth-xyz-ref"

        # return tV

    def getEditOp(self):
        tV = self._getStringValue("operation")
        if tV in ["edit"]:
            return "edit"
        elif tV in ["reset"]:
            return "reset"
        elif tV in ["undo"]:
            return "undo"
        elif tV in ["delete"]:
            return "delete"
        elif tV in ["global_edit_form"]:
            return "global_edit_form"
        elif tV in ["global_edit"]:
            return "global_edit"
        elif tV in ["global_edit_auth_seq"]:
            return "global_edit_auth_seq"
        elif tV in ["global_edit_menu", "globalmenuedit"]:
            return "global_edit_menu"
        elif tV in ["move", "edit_move"]:
            return "move"
        else:
            return "edit"

    def getSummaryOp(self):
        tV = self._getStringValue("operation")
        if tV in ["load", "reload"]:
            return "load"
        elif tV in ["loadexample"]:
            return "loadexample"
        elif tV in ["loadrcsb"]:
            return "loadrcsb"
        else:
            return "loadexample"

    #
    #
    def getSummaryAlignList(self, usingRevAllAlignIds=True):
        """Sequence ids for all entities/instances selected for alignment --"""
        tL = self._getRawValue("rev-allalignids")
        if (tL is not None) and (len(tL) > 0) and usingRevAllAlignIds:
            self.__saveDataStore("rev-allalignids", None)
            return tL
        #
        idString = self._getStringValue("allalignids")
        if self.__verbose:
            self.__lfh.write("\n+SeqModInputRequest.getSummaryAlignList() align list string: %r\n" % idString)
        if len(idString) > 0:
            if (tL is not None) and (len(tL) > 0):
                self.__saveDataStore("rev-allalignids", None)
            #
            idList = str(idString).split(",")
            return idList
        elif (tL is not None) and (len(tL) > 0):
            self.__saveDataStore("rev-allalignids", None)
            return tL
        else:
            return []
        #

    def setSummaryAlignList(self, alignIdList):
        """Sequence ids picked for alignment for all entities/instances --"""
        return self.setValue("allalignids", ",".join(alignIdList))

    def getAlignList(self):
        """Sequence ids picked from alignment for the current entity -"""
        idString = self._getStringValue("alignids")
        if len(idString) > 0:
            idList = str(idString).split(",")
            return idList
        else:
            return []
        #

    def setAlignList(self, alignIdList):
        """Sequence ids picked for alignment for all entities/instances --"""
        return self.setValue("alignids", ",".join(alignIdList))

    def __refIdEqual(self, id1, id2):
        """Are the two reference sequence id's pointing at the same entity/part"""
        try:
            fL1 = str(id1).split("_")
            fL2 = str(id2).split("_")
            if ((fL1[0] == fL2[0]) and (fL1[1] == fL2[1]) and (fL1[2] == fL2[2])) or (
                (fL1[0] in ["ref", "selfref"]) and (fL2[0] in ["ref", "selfref"]) and (fL1[1] == fL2[1]) and (fL1[2] == fL2[2])
            ):
                return True
            else:
                return False
        except:  # noqa: E722 pylint: disable=bare-except
            return False

    def setSummarySelectList(self, selectIdList):
        """Selected sequence ids for all entities/instances --"""
        return self.setValue("selectids", ",".join(selectIdList))

    def getSummarySelectList(self, usingRevSlectIds=True):
        """Sequence ids selected for all entities/instances.

        Preprocess list with any non-blank reference id stored in 'new-ref-seq-id'

        """
        tL = self._getRawValue("rev-selectids")
        if (tL is not None) and (len(tL) > 0) and usingRevSlectIds:
            self.__saveDataStore("rev-selectids", None)
            return tL

        idString = self._getStringValue("selectids")
        if self.__verbose:
            self.__lfh.write("\n+SeqModInputRequest.getSummarySelectList() select list string: %r\n" % idString)
        if len(idString) > 0:
            if (tL is not None) and (len(tL) > 0):
                self.__saveDataStore("rev-selectids", None)
            #
            idList = str(idString).split(",")
            newRefId = self._getStringValue("new-ref-seq-id")
            if len(newRefId) > 0:
                newIdList = []
                for sId in idList:
                    if self.__refIdEqual(sId, newRefId):
                        newIdList.append(newRefId)
                    else:
                        newIdList.append(sId)
                        #
                # update the current setting of 'selectids'
                self.setValue("selectids", ",".join(newIdList))
                return newIdList
            else:
                return idList
            #
        elif (tL is not None) and (len(tL) > 0):
            self.__saveDataStore("rev-selectids", None)
            return tL
        else:
            return []
        #

    def XgetAlignIdList(self):
        """Not used --"""
        idString = self._getStringValue("alignids")
        if self.__verbose:
            self.__lfh.write("\n+SequenceInputRequest() original idString %s\n" % idString)
        if (len(idString) > 0) and idString.startswith("['") and idString.endswith("']"):
            tString = idString[2:-2]
        else:
            tString = idString
        if self.__verbose:
            self.__lfh.write("\n+SequenceInputRequest() processed  idString %s\n" % tString)

        tList = tString.split(",")
        idL = []
        for tId in tList:
            idL.append(tId.strip())
        return idL

    def XgetSeqIdList(self):
        """Not used"""
        return []


#


class SeqModResponseContent(ResponseContent):
    def __init__(self, reqObj=None, verbose=False, log=sys.stderr):
        """
        Manage content items transfered from the application to the calling service.

        """
        super(SeqModResponseContent, self).__init__(reqObj, verbose=verbose, log=log)
        # self.__verbose = verbose
        # self.__lfh = log
        #
        self._kyList = [
            "title",
            "citationtitle",
            "pdbid",
            "identifier",
            "groupids",
            "alignids",
            "selectids",
            "conflictreportpath",
            "annotationreportpath",
            "aligntag",
            "viewalign_order",
            "editopid",
            "edittimestamp",
            "aligntimestamp",
            "warningtext",
        ]
        for ky in self._kyList:
            self._cD[ky] = ""
        #
        self._fgList = ["conflictreportflag", "annotationreportflag", "warningflag"]
        for fg in self._fgList:
            self._cD[fg] = False
        #

    #
    def setGroupIdList(self, gIdList):
        self._cD["groupids"] = ",".join([str(i) for i in gIdList])

    def setStructTitle(self, text):
        self._cD["title"] = text

    def setCitationTitle(self, text):
        self._cD["citationtitle"] = text

    def setPdbCode(self, code):
        self._cD["pdbid"] = code

    def setIdentifier(self, identifier):
        self._cD["identifier"] = identifier

    def setAlignIdList(self, alignIdList):
        self._cD["alignids"] = ",".join(alignIdList)

    def setSelectIdList(self, selectIdList):
        self._cD["selectids"] = ",".join(selectIdList)

    def setConflictReportPath(self, aPath):
        self._cD["conflictreportpath"] = aPath

    def setConflictReportFlag(self, aFlag):
        self._cD["conflictreportflag"] = aFlag

    def setAnnotationReportPath(self, aPath):
        self._cD["annotationreportpath"] = aPath

    def setAnnotationReportFlag(self, aFlag):
        self._cD["annotationreportflag"] = aFlag

    def setWarning(self, warnMsg=""):
        if not warnMsg:
            return
        #
        self._cD["warningtext"] = warnMsg
        self._cD["warningflag"] = True

    def setAlignTag(self, aTag):
        self._cD["aligntag"] = aTag

    def setViewAlignOrder(self, aOrder):
        """RPS, 20100727: Set code that determines display ordering of seq types for alignment interface."""
        self._cD["viewalign_order"] = aOrder

    def setEditOp(self, aOp):
        self._cD["editopid"] = aOp

    def setEditTimeStamp(self, tS):
        self._cD["edittimestamp"] = int(tS)

    def setAlignmentTimeStamp(self, tS):
        self._cD["aligntimestamp"] = int(tS)


if __name__ == "__main__":
    rC = SeqModResponseContent()
