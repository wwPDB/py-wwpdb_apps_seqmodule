##
# File:    AlignmentFrontEndUIEditor.py
# Date:    25-Jan-2010
# Updates:
#
# 27-Jan-2010 jdw - Improved handling of edits on "details"
# 14-Feb-2010 jdw - add support for tagged file names -
# 20-Apr-2010 jdw - ported to module seqmodule.
# 05-May-2010 jdw - Change returned residue value global edit.
# 09-May-2010 jdw - Add support for edit move operations.
# 15-May-2010 jdw - Refactor response and edit store operations (on-going)
# 17-May-2010 jdw - Revised undo for move
#  1-Apr-2013 jdw - adjust css classes for move operations.   Still not fully working for cases
#                   multiple correlated move/shifts.
# 15-Sep-2017 zf  - Move getEditStoreFilename() to SequenceEditStore.py
# 31-Oct-2018 zf  - Changed from AlignmentUtils to AlignmentTools and moved the functionality from self.__checkCompatibility()
#                   to AlignmentTools.checkCompatibility()
##

"""
Controlling class for managing edit operations on displayed sequence alignments.
"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.07"


import sys
import traceback

from wwpdb.apps.seqmodule.util.SequenceLabel import ResidueLabel
from wwpdb.apps.seqmodule.util.SequenceReferenceData import SequenceReferenceData
from wwpdb.apps.seqmodule.io.SequenceEditStore import getEditStoreFilename, SequenceEdit, SequenceEditStore
from wwpdb.apps.seqmodule.align.AlignmentTools import AlignmentTools


class AlignmentFrontEndUIEditor(object):
    """Controlling class for managing edit operations on displayed sequence alignments.

    Supported operations:

    -  edit :
    -  undo :
    -  delete :
    -  global_edit_menu :
    -  global_edit :
    -  move :

    """

    def __init__(self, reqObj=None, verbose=False, log=sys.stderr):
        """

        parameters:

        - reqObj :   [default=None]
        - verbose :  [default=False]
        - log :      [default=`sys.stderr`]

        """
        self.__verbose = verbose
        self.__lfh = log
        self.__reqObj = reqObj
        self.__srd = SequenceReferenceData(verbose=self.__verbose, log=self.__lfh)
        self.__gapSymbol = self.__srd.getGapSymbol()

        self.__alignmentTag = ""
        #
        self.__operation = None
        # self.__sessionObj = None
        self.__sessionId = None
        self.__cssClassRemove = "bgcolinsert bgcolconflict bgcolundo bgcolreplace bgcoldelete bgcoldna bgcolrna ovflow cf-misc-ref cf-misc-test cf-gap-test cf-gap-ref cf-glu-gln cf-asp-asn cf-ala-gly cf-rep-ALA cf-rep-ARG cf-rep-ASN cf-rep-ASP cf-rep-ASX cf-rep-CYS cf-rep-GLN cf-rep-GLU cf-rep-GLX cf-rep-GLY cf-rep-HIS cf-rep-ILE cf-rep-LEU cf-rep-LYS cf-rep-MET cf-rep-PHE cf-rep-PRO cf-rep-SER cf-rep-THR cf-rep-TRP cf-rep-TYR cf-rep-VAL cf-rep-SEC cf-rep-PYL"  # noqa: E501

        self.__setup()

    def __setup(self):
        try:
            # self.__sessionObj = self.__reqObj.getSessionObj()
            self.__operation = self.__reqObj.getEditOp()
            self.__sessionId = self.__reqObj.getSessionId()
            if self.__verbose:
                self.__lfh.write("+AlignmentFrontEndUIEditor.__setup()  session id %s\n" % self.__sessionId)
                self.__lfh.write("+AlignmentFrontEndUIEditor.__setup()  operation %s\n" % self.__operation)
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
                self.__lfh.write("+AlignmentFrontEndUIEditor.__setup()  Missing input for session id %s operation %s\n" % (self.__sessionId, self.__operation))
            #
        #

    def setAlignmentTag(self, tagValue):
        self.__alignmentTag = tagValue

    def getAlignmentTag(self):
        return self.__alignmentTag

    def edit(self):
        self.__lfh.write("+AlignmentFrontEndUIEditor __operation=%r\n" % self.__operation)
        if self.__operation == "edit":
            return self.__editAny()
        elif self.__operation == "delete":
            return self.__delete()
        elif self.__operation == "undo":
            rD = {}
            rD["editlist"] = self.__undo()
            return rD
        elif self.__operation == "global_edit":
            return self.__global_edit()
        elif self.__operation == "global_edit_auth_seq":
            return self.__global_edit_auth_seq()
        elif self.__operation == "global_edit_menu":
            rD = {}
            rD["editlist"] = self.__global_edit_menu()
            return rD
        elif self.__operation == "move":
            return self.__editMove()
        else:
            return {}
        #

    def undo(self):
        if self.__operation == "undo":
            return self.__undo()
        else:
            return {}
        #

    def delete(self):
        if self.__operation == "delete":
            return self.__delete()
        else:
            return {}
        #

    def __editMove(self):
        """Store edit arising for a drag-and-drop edit operations. These are translated into
        a pair of edit replace operations.

        """
        rDDict = {}
        edList = []
        srcResidueId = None
        dstResidueId = None
        srcVal = ""  # In case of exception
        try:
            srcResidueId = self.__reqObj.getValue("source")
            dstResidueId = self.__reqObj.getValue("destination")
            srcVal = self.__reqObj.getValue("sourceval")
            dstVal = self.__reqObj.getValue("destinationval")
            ses = self.__openEditStore()
            editOpLast = ses.getLastEditOp()
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
                self.__lfh.write(
                    "+AlignmentFrontEndUIEditor.__editMove()  Missing source %s or destination  %s  for session id %s\n" % (srcResidueId, dstResidueId, self.__sessionId)
                )
            #
            rD = self.__makeResponseDict(sid=srcResidueId, val=srcVal, val3=srcVal, editType="error", editOpId=None, newId=srcResidueId)
            if srcResidueId is not None:
                rDDict[srcResidueId] = rD
            return rDDict
        #
        # Handle NOOP -- or if move is NOT to a gap position --
        #
        if srcResidueId == dstResidueId:
            if self.__verbose:
                self.__lfh.write("+AlignmentMove.__editMove()  Edit(NOOP) for srcResidueId %s dstResidueId %s\n" % (srcResidueId, dstResidueId))
            #
            editType = "noop"
            rD = self.__makeResponseDict(sid=srcResidueId, val=srcVal, val3=srcVal, editType=editType, editOpId=editOpLast, newId=srcResidueId)
            rDDict[srcResidueId] = rD
            rD = self.__makeResponseDict(sid=dstResidueId, val=dstVal, val3=dstVal, editType=editType, editOpId=editOpLast, newId=dstResidueId)
            rDDict[srcResidueId] = rD
            return rDDict
        #
        if self.__verbose:
            self.__lfh.write(
                "+AlignmentFrontEndUIEditor.__editMove()  Source residue Id %s destination residue Id  %s  for session id %s\n" % (srcResidueId, dstResidueId, self.__sessionId)
            )
        #
        # Now decompose the residue ids into components -
        #
        srcLab = ResidueLabel()
        srcLab.unpack(srcResidueId)
        #
        # srcAlignPos = srcLab.getAlignmentIndex()
        srcSeqPos = srcLab.getSequenceIndex()
        srcLblInd = srcLab.getResidueLabelIndex()
        srcResType = srcLab.getResidueType()
        srcCode3 = srcLab.getResidueCode3()
        # deal with confusion in returned residue codes
        if len(srcVal) < len(srcCode3):
            srcCode1 = self.__srd.cnv3To1(srcCode3)
        else:
            srcCode3 = srcVal
            srcCode1 = self.__srd.cnv3To1(srcCode3)
        #
        dstLab = ResidueLabel()
        dstLab.unpack(dstResidueId)
        seqType = dstLab.getSequenceType()
        # dstAlignPos = dstLab.getAlignmentIndex()
        # dstSeqPos = dstLab.getSequenceIndex()
        # dstLblInd = dstLab.getResidueLabelIndex()
        # dstResType = dstLab.getResidueType()
        # dstCode3 = dstLab.getResidueCode3()
        #
        # Check if the moving is allowed
        #
        if (srcCode3 != "UNK") and ((seqType == "auth") or (seqType == "xyz")):
            aU = AlignmentTools(reqObj=self.__reqObj, entityId=self.__alignmentTag, verbose=self.__verbose, log=self.__lfh)
            errorMsg = aU.checkResidueMovingCompatibility(srcResidueId, dstResidueId)
            if errorMsg:
                return self.__makeResponseDict(sid=srcResidueId, val=srcCode1, val3=srcCode3, editType="error", editOpId=None, newId=None, errorText=errorMsg)
            #
        #
        # One edit applies to source residue -- which becomes a gap -- and retains its position in alignment.
        #
        newSrcVal = self.__gapSymbol
        srcLab.setResidueCode3(newSrcVal)
        srcLab.setResidueLabelIndex("")
        srcLab.setSequenceIndex("")
        srcIdNew = srcLab.pack()
        #
        # The second edit applies to the destination residue which assumes the identity of the source residue.
        #
        dstLab.setResidueCode3(srcCode3)
        dstLab.setResidueLabelIndex(srcLblInd)
        dstLab.setSequenceIndex(srcSeqPos)
        dstLab.setResidueType(srcResType)
        dstIdNew = dstLab.pack()
        #
        editOpNext = int(editOpLast) + 1
        editType = "replaceid"
        #
        # source ---
        sE = self.__makeSequenceEdit(targetId=srcResidueId, editType=editType, newValueList=[newSrcVal], priorValue=srcCode1, opId=editOpNext, newId=srcIdNew)
        edList.append(sE)
        if self.__verbose:
            sE.printIt(self.__lfh)
        #
        rD = self.__makeResponseDict(sid=srcResidueId, val=newSrcVal, val3=newSrcVal, editType=editType, editOpId=editOpNext, newId=srcIdNew)
        rDDict[srcResidueId] = rD
        #
        # destination --
        sE = self.__makeSequenceEdit(targetId=dstResidueId, editType=editType, newValueList=[srcCode3], priorValue=dstVal, opId=editOpNext, newId=dstIdNew)
        edList.append(sE)
        if self.__verbose:
            sE.printIt(self.__lfh)
        #
        rD = self.__makeResponseDict(sid=dstResidueId, val=srcCode1, val3=srcCode3, editType=editType, editOpId=editOpNext, newId=dstIdNew)
        rDDict[dstResidueId] = rD
        #
        ses.storeEditList(edList)
        #
        rDTop = {}
        rDTop["editlist"] = rDDict
        return rDTop

    def __openEditStore(self):
        #
        esfn = getEditStoreFilename(self.__alignmentTag)
        if self.__verbose:
            self.__lfh.write("+AlignmentFrontEndUIEditor.__openEditStore() sessionid %s filepath %s\n" % (self.__sessionId, esfn))
        #
        ses = SequenceEditStore(sessionObj=self.__reqObj.newSessionObj(), fileName=esfn, verbose=self.__verbose, log=self.__lfh)
        return ses

    def __makeResponseDict(self, sid, val, val3, editType, editOpId=None, newId=None, errorText="Unable to save the Edit.<br />\n"):
        """Convenience method for constructing the response for the edit operation.

        A response consists of a standard dictionary of the following attributes:
        - sid       the unique element identifier of the edited residue
        - val      displayed of the edited element
        - val3     value of edited element as the residue 3-letter-code
        - edittype     the type of edit performed (insert,replace,delete,undo,replaceid,...)
        - editiopid    integer identifier for this edit operation
        - classAdd     the CSS classes to be added to the edited element
        - classRemove  the CSS classes to be removed from the edited element
        - errortext    text to describe any error related to the edit operation
        - debug        text string to aid in debugging the edit operation
        - tooltip      informational text to describe the edit operation
        - newid        if the editType is a replaceId type then the newId will be included.
        """
        rD = {}
        rD["id"] = str(sid)
        rD["val"] = str(val)
        rD["val3"] = str(val3)
        rD["edittype"] = str(editType)
        if editOpId is not None:
            rD["editopid"] = editOpId
        #
        rD["errorflag"] = False
        rD["errortext"] = "none"
        rD["debug"] = "%s = %s(%s)" % (sid, val, val3)
        #
        rD["tooltip"] = self.__setToolTipText(editType, val)
        rD["classAdd"] = self.__setCssClassAdd(editType, sid, val)
        rD["classRemove"] = self.__setCssClassRemove(editType)

        if editType == "replaceid":
            if newId is not None:
                rD["newid"] = newId
            #
        elif editType == "error":
            rD["errorflag"] = True
            rD["errortext"] = errorText
        elif editType == "details":
            rD["edittype"] = "replace"
        else:
            pass
        #
        return rD

    def __makeSequenceEdit(self, targetId, editType, newValueList, priorValue, opId, newId=None):
        """Convenience method for constructing sequence edit objects."""
        sE = SequenceEdit(self.__verbose)
        sE.setTargetElementId(targetId)
        sE.setEditType(editType)
        sE.setValueNew(newValueList)
        sE.setValuePrevious(priorValue)
        sE.setEditOpId(opId)
        sE.setNewElementId(newId)
        #
        return sE

    def __setToolTipText(self, editType, value=""):
        text = ""
        if editType == "replace":
            text = "Edited residue:  %s" % value
        elif editType == "insert":
            text = "Edited residue:  %s" % value
        elif editType == "undo":
            text = "Restored residue:  %s" % value
        elif editType == "delete":
            text = "Residue %s marked for deletion" % value
        elif editType == "replaceid":
            text = "Edited residue:  %s" % value
        elif editType == "noop":
            text = "Residue:  %s" % value
        elif editType == "details":
            text = "Residue details:  %s" % value
        elif editType == "error":
            text = "Error processing edit for residue:  %s" % value
        else:
            pass
        #
        return text

    def __setCssClassAdd(self, editType="replace", tId="", val=""):
        eventClassString = ""
        if tId:
            rLabel = ResidueLabel()
            rLabel.unpack(tId)
            if rLabel.getSequenceType() == "xyz":
                eventClassString = "draggable "
            elif (rLabel.getSequenceType() == "auth") or (rLabel.getSequenceType() == "ref" and val == self.__gapSymbol):
                eventClassString = "dblclick draggable "
            #
        #
        cssClassString = ""
        if editType == "replace":
            cssClassString = eventClassString + "bgcolreplace "
        elif editType == "insert":
            cssClassString = eventClassString + "bgcolinsert ovflow "
        elif editType == "undo":
            cssClassString = " bgcolundo"
        elif editType == "delete":
            cssClassString = eventClassString + "bgcolinsert ovflow "
        elif editType == "replaceid":
            cssClassString = eventClassString + "bgcolreplace "
        elif editType == "noop":
            cssClassString = eventClassString
        elif editType == "details":
            cssClassString = "dblclickselect bgcolreplace "
        elif editType == "error":
            cssClassString = eventClassString
        else:
            cssClassString = eventClassString
        #
        if val == self.__gapSymbol:
            cssClassString += " cf-gap-test "
        #
        return cssClassString

    def __setCssClassRemove(self, editType=""):  # pylint: disable=unused-argument
        return self.__cssClassRemove

    def __editAny(self):
        """Store individual edit operations in the the edit store -

        Distiguishes if residues or details are provided.
        """
        residueId = None
        editValue = None
        priorValue = None
        rD = {}
        try:
            residueId = self.__reqObj.getValue("id")
            editValue = str(self.__reqObj.getValue("value")).upper()
            # empty residues treated as gaps -
            if len(editValue) < 1:
                editValue = self.__gapSymbol
            #
            priorValue = self.__reqObj.getValue("priorvalue")
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
                self.__lfh.write(
                    "+AlignmentFrontEndUIEditor.__editAny()  Missing input residueId %s value %s prior value %s session id %s\n"
                    % (residueId, editValue, priorValue, self.__sessionId)
                )
            #
            rD = self.__makeResponseDict(sid=residueId, val=editValue, val3=editValue, editType="error", editOpId=None, newId=None)
            return rD
        #
        # Distinguish what we are editing using extensions on the residue id.
        #   A "_D" extension marks a details field.
        #
        if residueId.endswith("_D"):
            if self.__verbose:
                self.__lfh.write("+AlignmentFrontEndUIEditor.__editAny()  Editing details for residueId %s value %s\n" % (residueId, editValue))
            #
            return self.__editDetails()
        #
        # Handle NOOP -- edit
        #
        if editValue == priorValue:
            if self.__verbose:
                self.__lfh.write("+AlignmentFrontEndUIEditor.__editAny()  Edit(NOOP) for residueId %s value %s prior value %s\n" % (residueId, editValue, priorValue))
            #
            rD = self.__makeResponseDict(sid=residueId, val=editValue, val3=editValue, editType="noop", editOpId=None, newId=None)
            return rD
        #
        return self.__editResidue(residueId, editValue, priorValue)

    def __editResidue(self, residueId=None, editValue=None, priorValue=None):
        """Store individual edit operations for this residue value in the the edit store -"""
        #
        # Need residue type to interpret edits - using the residue type
        rLabel = ResidueLabel()
        rLabel.unpack(residueId)

        polymerTypeCode = rLabel.getResidueType()

        # Get the list of 1-letter and 3-letter edits from the input editValue
        (r1L, r3L) = self.__srd.parseSequence(editValue, polymerTypeCode)
        #
        # Is this a replacement or replacement + insertion
        #
        r1 = " ".join(r1L)
        r3 = " ".join(r3L)
        #
        #
        if len(r3L) > 1:
            eType = "insert"
        else:
            eType = "replace"
            #
            # Check if residue name replacement is allowed
            #
            if r3L[0] != "UNK":
                aU = AlignmentTools(reqObj=self.__reqObj, entityId=self.__alignmentTag, verbose=self.__verbose, log=self.__lfh)
                errorMsg = aU.checkResidueNameReplaceCompatibility(rLabel.getAlignmentIndex(), rLabel.getResidueCode3(), rLabel.getResidueLabelIndex(), r3L[0])
                if errorMsg:
                    return self.__makeResponseDict(sid=residueId, val=editValue, val3=editValue, editType="error", editOpId=None, newId=None, errorText=errorMsg)
                #
            #
        #
        ses = self.__openEditStore()
        editOpLast = ses.getLastEditOp()
        editOpNext = int(editOpLast) + 1
        sE = self.__makeSequenceEdit(targetId=residueId, editType=eType, newValueList=r3L, priorValue=priorValue, opId=editOpNext, newId=None)
        if self.__verbose:
            self.__lfh.write("+AlignmentFrontEndUIEditor.__editResidue() dumping the current residue edit\n")
            sE.printIt(self.__lfh)
        #
        ses.storeEdit(sE)
        #
        rD = self.__makeResponseDict(sid=residueId, val=r1, val3=r3, editType=eType, editOpId=editOpNext, newId=None)
        #
        return rD

    def __editDetails(self):
        """Store individual detail edit operations for this residue value in the the edit store -"""
        residueId = None
        editValue = None
        priorValue = None
        rD = {}
        try:
            residueId = self.__reqObj.getValue("id")
            editValue = self.__reqObj.getValue("value")
            priorValue = self.__reqObj.getValue("priorvalue")
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
                self.__lfh.write(
                    "+Alignment.__editDetails()  Missing input residueId %s value %s prior value %s session id %s\n" % (residueId, editValue, priorValue, self.__sessionId)
                )
            #
            rD = self.__makeResponseDict(sid=residueId, val=editValue, val3=editValue, editType="error", editOpId=None, newId=None)
            return rD
        #
        # Check for sequence details -
        #
        if not residueId.endswith("_D"):
            if self.__verbose:
                self.__lfh.write("+Alignment.__editDetails()  Editing details with bad residueId %s value %s\n" % (residueId, editValue))
            # should never get here !!
            return rD
        #
        eType = "details"

        if self.__verbose:
            self.__lfh.write("+Alignment.__editDetails() residue %s edit type  %s value (%s)  sessionid %s\n" % (residueId, eType, editValue, self.__sessionId))
        #
        ses = self.__openEditStore()

        editOpLast = ses.getLastEditOp()
        editOpNext = int(editOpLast) + 1
        #
        sE = self.__makeSequenceEdit(targetId=residueId[:-2], editType=eType, newValueList=editValue, priorValue=priorValue, opId=editOpNext, newId=None)
        if self.__verbose:
            sE.printIt(self.__lfh)
        #
        ses.storeEdit(sE)
        #
        # Setup the return values -
        #
        rD = self.__makeResponseDict(sid=residueId, val=editValue, val3=editValue, editType=eType, editOpId=editOpNext, newId=None)
        #
        return rD

    def __undo(self):
        """Undo the last edit operation..

        2010-02-15 -  Special handling of undo operations for menu edits
        """
        rD = {}
        edOp = 0
        if self.__verbose:
            self.__lfh.write("+Alignment.__undo()  Starting in session id %s\n" % (self.__sessionId))
        #
        try:
            edOp = self.__reqObj.getValue("editopid")
            if len(edOp) < 1:
                edOp = 0
            edOp = int(edOp)
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
                self.__lfh.write("+Alignment.__undo()  Missing target editOp for the undo operation in session id %s\n" % (self.__sessionId))
            #
            rD = self.__makeResponseDict(sid="", val="", val3="", editType="error", editOpId=None, newId=None)
            return rD
        #
        if self.__verbose:
            self.__lfh.write("+Alignment.__undo()  input editOp is %s for session id %s\n" % (edOp, self.__sessionId))
        #
        cssUndoClass = " bgcolundo"
        cssClassRemove = self.__cssClassRemove

        if int(edOp) != 0:
            ses = self.__openEditStore()
            editOpLast = ses.getLastEditOp()
            editOpNext = max(editOpLast - 1, 0)
            #
            if self.__verbose:
                self.__lfh.write("+Alignment.__undo() undo edit op id %d last op id %d edit list length %d\n" % (edOp, editOpLast, ses.length()))
                ses.printIt(self.__lfh)
            #
            eObjList = ses.get(editOpLast)
            for eObj in eObjList:
                tD = {}
                # undo this edit.
                eType = eObj.getEditType()
                prevVal = eObj.getValuePrevious()
                idTarget = eObj.getTargetElementId()
                prevStyle = eObj.getStylePrevious()
                newId = eObj.getNewElementId()
                #
                if self.__verbose:
                    self.__lfh.write("+Alignment.__undo()  Undo edit type %s previous value %s on target id %s editOpNext %d\n" % (eType, prevVal, idTarget, editOpNext))
                #
                if eType == "replace":
                    tD["id"] = idTarget
                    tD["val"] = prevVal
                    if (prevStyle is not None) and (len(prevStyle) > 0):
                        tD["classAdd"] = prevStyle
                        tD["classRemove"] = cssClassRemove + cssUndoClass + " " + prevStyle
                    else:
                        tD["classAdd"] = cssUndoClass
                        tD["classRemove"] = cssClassRemove
                    #
                    tD["editopid"] = editOpNext
                    tD["tooltip"] = "Restored previous value: " + prevVal

                    rD[idTarget] = tD
                elif eType == "replaceid":
                    tD["id"] = newId
                    tD["newid"] = idTarget
                    tD["val"] = prevVal
                    if (prevStyle is not None) and (len(prevStyle) > 0):
                        tD["classAdd"] = prevStyle
                        tD["classRemove"] = cssClassRemove + cssUndoClass + " " + prevStyle
                    else:
                        tD["classAdd"] = cssUndoClass
                        tD["classRemove"] = cssClassRemove
                    #
                    tD["editopid"] = editOpNext
                    tD["tooltip"] = "Restored previous value: " + prevVal

                    if tD["val"] == self.__gapSymbol:
                        tD["classAdd"] += " cf-gap-test "
                    #
                    rD[newId] = tD
                elif eType == "insert":
                    tD["id"] = idTarget
                    tD["val"] = prevVal
                    tD["classAdd"] = cssUndoClass
                    tD["classRemove"] = cssClassRemove
                    tD["editopid"] = editOpNext
                    tD["tooltip"] = "Restored previous value: " + prevVal
                    rD[idTarget] = tD
                elif eType == "details":
                    # add the "_D" suffix onto the iD to distinguish this is a details item.
                    tId = idTarget + "_D"
                    tD["id"] = tId
                    tD["val"] = prevVal
                    tD["classAdd"] = cssUndoClass
                    tD["classRemove"] = cssClassRemove
                    tD["editopid"] = editOpNext
                    tD["tooltip"] = "Restored previous details: " + prevVal
                    rD[tId] = tD
                elif eType == "delete":
                    tD["id"] = idTarget
                    tD["val"] = prevVal
                    tD["classAdd"] = cssUndoClass
                    tD["classRemove"] = cssClassRemove
                    tD["editopid"] = editOpNext
                    tD["tooltip"] = "Restored previous value: " + prevVal
                    rD[idTarget] = tD
                elif eType == "deletelist":
                    tD["id"] = idTarget
                    tD["val"] = prevVal
                    tD["classAdd"] = cssUndoClass
                    tD["classRemove"] = cssClassRemove
                    tD["editopid"] = editOpNext
                    tD["tooltip"] = "Restored previous value: " + prevVal
                    rD[idTarget] = tD
                else:
                    pass
                #
            #
            if self.__verbose:
                self.__lfh.write("+Alignment.__undo() removing edit op %d\n" % editOpLast)
            #
            ses.remove(editOpLast)
        #
        return rD

    def __delete(self):
        """Mark the selected residues for deletion and record the deletion in the edit store.

        Operation is performed using separate delete operations for each residue.
        """
        rD = {}
        try:
            dList = self.__reqObj.getDeleteList()
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
                self.__lfh.write("+Alignment.__delete()  Missing delete selection for session id %s\n" % (self.__sessionId))
            #
            rD = self.__makeResponseDict(sid="", val="", val3="", editType="error", editOpId=None, newId=None)
            return rD
        #
        cssClass = "bgcoldelete"
        #
        if len(dList) > 0:
            rD["classAdd"] = cssClass
            rD["classRemove"] = self.__cssClassRemove
            rD["errortext"] = ""
            rD["errorflag"] = False
            rD["debug"] = "Marking deletions: " + " ".join(dTup[0] for dTup in dList)
            #
            ses = self.__openEditStore()
            editOpLast = ses.getLastEditOp()
            editOpNext = int(editOpLast) + 1
            edList = []
            for dTup in dList:
                rLab = ResidueLabel()
                rLab.unpack(dTup[0])
                sE = SequenceEdit(self.__verbose)
                sE.setValueNew("")
                sE.setValuePrevious(dTup[1])
                sE.setEditType("delete")
                sE.setTargetElementId(dTup[0])
                sE.setEditOpId(editOpNext)
                edList.append(sE)

                if self.__verbose:
                    sE.printIt(self.__lfh)
                #
            #
            rD["editopid"] = editOpNext
            rD["edittype"] = "delete"
            #
            ses.storeEditList(edList)
        #
        return rD

    def __global_edit(self):
        """Perform edit operations from global edit form."""
        aU = AlignmentTools(reqObj=self.__reqObj, entityId=self.__alignmentTag, verbose=self.__verbose, log=self.__lfh)
        editList, unCompatibleList = aU.getGlobalEditList()
        if len(unCompatibleList) > 0:
            message = ""
            for unCompatibleTup in unCompatibleList:
                message += "Residue '" + unCompatibleTup[2] + " " + unCompatibleTup[3] + " " + unCompatibleTup[4] + "' can not be moved to position '" + unCompatibleTup[6]
                if (unCompatibleTup[5] != "") and (unCompatibleTup[5] != self.__gapSymbol):
                    message += "' where it is '" + unCompatibleTup[5] + "' in the auth sequence.<br/>\n"
                else:
                    message += "' where it is a gap in the auth sequence.<br/>\n"
                #
            #
            return self.__makeResponseDict(
                sid=unCompatibleList[0][0], val=unCompatibleList[0][1], val3=unCompatibleList[0][3], editType="error", editOpId=None, newId=None, errorText=message
            )
        #
        ses = self.__openEditStore()
        editOpLast = ses.getLastEditOp()
        editOpNext = int(editOpLast) + 1
        # editType = "replaceid"
        #
        rDDict = {}
        edList = []
        #
        for eTup in editList:
            sE = self.__makeSequenceEdit(targetId=eTup[0], editType="replaceid", newValueList=[eTup[2]], priorValue=eTup[3], opId=editOpNext, newId=eTup[4])
            edList.append(sE)
            #
            rD = self.__makeResponseDict(sid=eTup[0], val=eTup[1], val3=eTup[2], editType="replaceid", editOpId=editOpNext, newId=eTup[4])
            rDDict[eTup[0]] = rD
        #
        if self.__verbose:
            self.__lfh.write("+Alignment.__global_edit()  Saving edit list of length %d\n" % len(edList))
        #
        ses.storeEditList(edList)
        #
        rDTop = {}
        rDTop["editlist"] = rDDict
        return rDTop

    def __global_edit_auth_seq(self):
        """Perform to insert missing author sequence part(s) based on reference sequence."""
        aU = AlignmentTools(reqObj=self.__reqObj, entityId=self.__alignmentTag, verbose=self.__verbose, log=self.__lfh)
        errMsg, editList = aU.getGlobalEditAuthSeqList()
        if errMsg:
            rD = {}
            rD["errorflag"] = True
            rD["errortext"] = errMsg
            return rD
        #
        ses = self.__openEditStore()
        editOpLast = ses.getLastEditOp()
        editOpNext = int(editOpLast) + 1
        # editType = "replace"
        #
        rDDict = {}
        edList = []
        #
        for eTup in editList:
            sE = self.__makeSequenceEdit(targetId=eTup[0], editType="replace", newValueList=[eTup[2]], priorValue=eTup[3], opId=editOpNext, newId=eTup[4])
            edList.append(sE)
            #
            rD = self.__makeResponseDict(sid=eTup[0], val=eTup[1], val3=eTup[2], editType="replace", editOpId=editOpNext, newId=eTup[4])
            rDDict[eTup[0]] = rD
        #
        if self.__verbose:
            self.__lfh.write("+Alignment.__global_edit_auth_seq()  Saving edit list of length %d\n" % len(edList))
        #
        ses.storeEditList(edList)
        #
        rDTop = {}
        rDTop["editlist"] = rDDict
        return rDTop

    def __global_edit_menu(self):
        """Perform edit operations from the global edit menu."""
        #
        sList = self.__reqObj.getEditMenuSelectList()
        if self.__verbose:
            self.__lfh.write("+AlignmentFrontEndUIEditor.__global_edit_menu()  Select list length is %d\n" % len(sList))
        #
        if len(sList) < 1:
            return {}
        #
        ses = self.__openEditStore()
        editOpLast = ses.getLastEditOp()
        editOpNext = int(editOpLast) + 1
        #
        rDTop = {}
        edList = []
        for sTup in sList:
            tId = sTup[0]
            tPriorValue = sTup[1]
            # up case the residue value
            tNewValue = str(sTup[2]).upper()
            tPriorCss = sTup[3]
            sE, rD = self.__get_edit_util("menu", tId, tId, tPriorValue, tNewValue, tPriorCss, editOpNext, " ovflow", "Global menu edit on residue: %s")
            edList.append(sE)
            rDTop[tId] = rD
        #
        if self.__verbose:
            self.__lfh.write("+Alignment.__global_edit_menu()  Saving edit list of length %d\n" % len(edList))
        #
        ses.storeEditList(edList)
        #
        return rDTop

    def __get_edit_util(self, oType, tId, newId, tPriorValue, editValue, tPriorCss, editOpNext, overflow, tooltip):
        #
        # Need residue type to interpret edits - using the residue type
        rLabel = ResidueLabel()
        rLabel.unpack(tId)
        polymerTypeCode = rLabel.getResidueType()
        #
        eventClassString = ""
        if rLabel.getSequenceType() == "xyz":
            eventClassString = "draggable"
        elif rLabel.getSequenceType() == "auth":
            eventClassString = "dblclick"
            if oType == "menu":
                eventClassString += " draggable"
            #
        #
        # Get the list of 1-letter and 3-letter edits from the input editValue
        (r1L, r3L) = self.__srd.parseSequence(editValue, polymerTypeCode)
        #
        # Is this a replacement or replacement + insertion
        #
        r1 = " ".join(r1L)
        r3 = " ".join(r3L)
        #
        if len(r3L) > 1:
            eType = "insert"
            cssClass = eventClassString + " bgcolinsert" + overflow
        else:
            eType = "replace"
            cssClass = eventClassString + " bgcolreplace"
        #
        sE = SequenceEdit(self.__verbose)
        sE.setValueNew(r3L)
        sE.setValuePrevious(tPriorValue)
        sE.setEditType(eType)
        sE.setTargetElementId(newId)
        sE.setEditOpId(editOpNext)
        if oType == "menu":
            sE.setStylePrevious(tPriorCss)
        #
        if self.__verbose:
            sE.printIt(self.__lfh)
        #
        rD = {}
        rD["val"] = r1
        rD["val3"] = r3
        rD["id"] = newId
        rD["classAdd"] = cssClass
        rD["classRemove"] = self.__cssClassRemove
        rD["errortext"] = "none"
        rD["debug"] = "%s = %s(%s)" % (tId, r1, r3)
        rD["tooltip"] = "Global menu edit on residue: %s" % r3
        rD["tooltip"] = tooltip % r3
        rD["editopid"] = editOpNext
        rD["edittype"] = eType
        #
        return sE, rD
