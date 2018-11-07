##
# File:    AlignmentViewRcsbTests.py
# Date:    7-Feb-2010
#
# Update:
#  12-Mar-2010 jdw  Revised example handling.
#  20-Apr-2010 jdw  Ported to module seqmodule.
##
"""
Test cases for alignment view class.

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"


import sys, unittest, traceback
import time, os, os.path

from wwpdb.apps.seqmodule.webapp.WebRequest     import SequenceInputRequest
from wwpdb.apps.seqmodule.align.AlignmentView   import AlignmentView
from wwpdb.apps.seqmodule.io.SequenceDataImport import SequenceDataImportRcsb
from wwpdb.apps.seqmodule.io.SequenceDataStore  import SequenceDataStore
from wwpdb.apps.seqmodule.util.SequenceLabel    import SequenceLabel

class AlignmentViewRcsbTests(unittest.TestCase):
    def setUp(self):
        #
        self.__verbose=True
        self.__lfh = sys.stderr
        # Create a session object and session directories for test cases
        self.__reqObj=SequenceInputRequest(paramDict={},verbose=self.__verbose,log=self.__lfh)
        self.__reqObj.setValue("SessionPath","../sessions")
        self.__sobj=self.__reqObj.newSessionObj()
        self.__sessionPath=self.__sobj.getPath()
        #
        # Local example data -
        #
        self.__reqObj.setValue("RcsbDataPath","../examples/rcsb-data")
        self.__reqObj.setValue("RcsbReferenceSequencePath","../examples/rcsb-sequence")
        self.__rcsbIdExample='RCSB101544'
        #
        self.__seqIds=''
        self.__reqObj.setValue("alignids",self.__seqIds)
        self.__reqObj.printIt(self.__lfh)        
        #

        #
        #self.__rcsbId='RCSB101544'
        #self.__rcsbId='RCSB055847'
        # ribosome test
        #self.__rcsbId= 'RCSB056215'
        #self.__rcsbId='rcsb057823'
        #
        # Import data 
        self.__sdi=SequenceDataImportRcsb(reqObj=self.__reqObj,rcsbId=self.__rcsbIdExample,fileSource='rcsb-repository',
                                          verbose=self.__verbose,log=self.__lfh,fetchPdbFile=False)
        self.__sdi.doImport()
        self.__sda=SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose,log=self.__lfh)
        #
        

    def tearDown(self):
        pass
    
    def testRcsbAlign1(self): 
        """ Create alignment view for a specific rcsb test case.
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            idString="['auth_A_1_1,xyz_A_1_1,ref_A_1_1']"
            self.__reqObj.setValue("operation","load")                        
            self.__reqObj.setValue("alignids",idString)         
            self.__reqObj.printIt(sys.stderr)
            aV=AlignmentView(reqObj=self.__reqObj, verbose=self.__verbose,log=self.__lfh)
            oL=aV.loadAlign()
            aV.dump(self.__lfh)
            self.__lfh.write("\n\nRendered alignment length %d:\n" % len(oL))
            self.__lfh.write("%s\n" % "\n".join(oL))
            #
            cL=aV.renderConflicts()
            self.__lfh.write("\n\nRendered conflict report length %d:\n" % len(cL))
            self.__lfh.write("%s\n" % "\n".join(cL))            

        except:
            traceback.print_exc(file=sys.stderr)
            self.fail()

    def testAlignReference(self): 
        """ Create alignment view for specific rcsb test case (reference sequences).
        """
        self.__lfh.write("\n------------------------ ")
        self.__lfh.write("Starting test function %s" % sys._getframe().f_code.co_name)
        self.__lfh.write(" -------------------------\n")
        try:
            sLab=SequenceLabel(self.__verbose)
            seqIds=self.__sda.getIds(dataType="sequence",seqType="ref")
            idList=[]
            for seqId in seqIds:
                altIds=self.__sda.getAlternativeIds(seqId, dataType="sequence", seqType="ref")
                for altId in altIds:
                    sLab.set(seqType='ref', seqInstId=seqId,  seqAltId=altId, seqVersion=1)
                    idList.append(sLab.pack())

            alignIds=",".join(idList)
            print alignIds

            self.__reqObj.setValue("operation","load")                        
            self.__reqObj.setValue("alignids",alignIds)         
            self.__reqObj.printIt(self.__lfh)
            aV=AlignmentView(reqObj=self.__reqObj, verbose=self.__verbose,log=self.__lfh)
            oL=aV.loadAlign()
            self.__lfh.write("\n\nRendered alignment length %d:\n" % len(oL))
            self.__lfh.write("%s\n" % "\n".join(oL))

        except:
            traceback.print_exc(file=sys.stderr)
            self.fail()

def suite():
    return unittest.makeSuite(AlignmentViewRcsbTests,'test')

if __name__ == '__main__':
    unittest.main()


