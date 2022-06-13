##
# File:    SequenceDataPrepTests.py
# Date:    03-Mar-2013
#
# Updates:
#
#   21-May-2014 jdw   add default alignments to preparation process --
#
##
"""
Test cases for preparing stating data for sequence module test cases.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.07"


import sys
import unittest
import traceback
import time
import inspect
import os
import os.path

from wwpdb.apps.seqmodule.control.SequenceDataAssemble_v2 import SequenceDataAssemble
from wwpdb.apps.seqmodule.align.AlignmentStatistics import AlignmentStatistics

# from wwpdb.apps.seqmodule.align.MultiAlignPseudo import MultiAlignPseudo
from wwpdb.apps.seqmodule.control.SummaryView_v2 import SummaryView
from wwpdb.apps.seqmodule.control.SummaryViewDepiction import SummaryViewDepiction

from wwpdb.utils.config.ConfigInfo import ConfigInfo, getSiteId
from wwpdb.apps.seqmodule.webapp.SeqModWebRequest import SeqModInputRequest

from wwpdb.io.locator.PathInfo import PathInfo
from wwpdb.utils.dp.DataFileAdapter import DataFileAdapter


class SequenceDataPrepTests(unittest.TestCase):
    def setUp(self):
        #
        self.__verbose = True
        self.__lfh = sys.stderr

        self.__pathExamplesRel = "../tests"
        self.__pathExamples = os.path.abspath(self.__pathExamplesRel)

        self.__siteId = None
        self.__cI = None
        self.__topPath = None
        self.__topSessionPath = None
        self.__reqObj = None
        self.__templatePath = None
        # self.__sessionId = None
        self.__sessionObj = None
        self.__sessionPath = None
        #
        # self.__examFileList    = ['1cbs.cif','3rer.cif','rcsb056751.cif']
        # self.__exampleFileList    = ['rcsb056751.cif']
        # self.__exampleIdList    = ['rcsb056751','rcsb057750','rcsb055453','rcsb057171','rcsb052365','rcsb056215',
        #              'rcsb053659','rcsb057584','rcsb076237','rcsb056751','rcsb095269','rcsb057584',
        # 'rcsb078538,','rcsb052365','rcsb078027','rcsb051659','rcsb075964','rcsb078385']
        # self.__exampleIdList= ['rcsb075964']
        #
        # self.__exampleIdList= ['rcsb052365','rcsb078488']
        self.__exampleIdList = ["rcsb078858"]
        self.__exampleFileType = "rcsb-mmcif"
        # SS example with depositor details -
        # self.__exampleIdList= ['D_800715']
        # JY example for Demo
        self.__exampleIdList = ["D_123563"]
        self.__exampleIdList = ["D_3000346"]
        self.__exampleIdList = ["4erd"]
        self.__exampleIdList = ["D_800003"]
        self.__exampleIdList = ["D_123584"]
        self.__exampleIdList = ["D_123578"]
        self.__exampleIdList = ["D_123607"]
        # ST-47
        self.__exampleIdList = ["D_1004021"]
        # ST-37
        self.__exampleIdList = ["D_4000100"]
        # ST-51
        self.__exampleIdList = ["D_3000362"]
        # irina jcsg test case
        self.__exampleIdList = ["D_123689"]
        self.__exampleIdList = ["D_123694"]
        self.__exampleIdList = ["D_058993"]
        #
        # self.__exampleIdList= ['D_123614']
        # self.__exampleIdList= ['D_123592']
        # self.__exampleIdList= ['D_123869']
        #
        # self.__exampleIdList= ['D_123889']
        # self.__exampleFileType='pdbx-mmcif'
        # self.__exampleFileType='rcsb-cifeps'
        #
        # self.__exampleIdList= ['rcsb055453']
        #
        # case sensitive test
        # self.__exampleIdList= ['rcsb080603']
        # self.__exampleFileType='rcsb-cifeps'
        #
        # self.__exampleIdList= ['D_000002']
        # self.__exampleIdList= ['D_180262']
        # self.__exampleIdList= ['D_180252']
        # self.__exampleIdList= ['D_1100200023']
        # self.__exampleIdList= ['D_180254']
        # self.__exampleIdList= ['D_1000200017']
        # self.__exampleIdList= ['D_180233']
        # latest example with conflict table issue
        # self.__exampleIdList= ['D_1100200179']
        # self.__exampleIdList= ['D_080826']
        # self.__exampleIdList= ['D_083112']
        # self.__exampleFileType='pdbx-mmcif'
        # chimera example
        # self.__exampleIdList= ['rcsb056751']
        # multiple polymer entity example
        # self.__exampleIdList= ['rcsb082726']
        # self.__exampleIdList= ['D_1100200023']
        #
        # self.__exampleIdList= ['D_1000000000','D_1000000001']
        # self.__exampleIdList= ['D_1000000007']
        # self.__exampleFileType='rcsb-mmcif'
        self.__exampleFileType = "pdbx-mmcif"
        #
        # self.__dsList             = ['D_083014']
        # self.__dsList             = ['D_1000000004']
        # self.__dsList             = ['D_083043']
        # self.__dsList             = ['D_083325']
        # self.__dsList             = ['D_082926']
        # self.__dsList             = ['D_1000000006']
        # self.__dsList             = ['D_082773']
        # self.__dsList =['D_1100201098']
        # self.__dsList             = ['D_083747']
        # self.__dsList = ['D_083314']
        # self.__dsList = ['D_1000200095']
        # self.__dsList=['D_1000200160']
        # self.__dsList=['D_084021']
        # self.__dsList=['D_1000200142']
        # self.__dsList=['D_1000200101']
        # self.__dsList=['D_1000200353']
        # self.__dsList=['D_1000200141']
        # self.__dsList=['D_1000200502']
        # self.__dsList=['D_1000200473']
        # self.__dsList=['D_1000200641']
        # self.__dsList=['D_1000200076']
        # self.__dsList=['D_1000200573']
        #
        # self.__dsList=['D_1000200581']
        # self.__dsList=['D_1000200747']
        # self.__dsList=['D_1000200782']
        # self.__dsList=['D_1000200141']
        # self.__dsList=['D_1000200774']
        # self.__dsList=['D_1000200246']
        # self.__dsList=['D_1000000000']

        # self.__dsList= ['D_1000201233']
        # self.__dsList=['D_1000201428']
        # self.__dsList=['D_1000201443']
        # pre-defined chimera example --
        # self.__dsList             = ['D_1000000005']
        # multi-entity example -
        # self.__dsList=['D_1000000000']
        # PYL/O examples
        # self.__dsList=['D_1000000013','D_1000000014']
        # self.__dsList=['D_1000000013','D_1000000014','D_1000000005','D_1000000000']
        #  ambiguous alignment edit --
        # self.__dsList=['D_1000200868']
        # sample sequence deletion with move around gap issue
        # self.__dsList=['D_1000201045']
        #
        # Ezra's ribosome entry -
        # self.__dsList=['D_1000200788']
        #
        #  Synthetic chimera --
        # self.__dsList=['D_1000201612']
        #
        # Ezra deletion edit --
        # self.__dsList=[' D_1100201537']
        # self.__dsList=['D_1000201880']
        #
        # color test 3rd entity --
        # self.__dsList=['D_1000201632']
        #
        # ref db alignment update
        # self.__dsList=['D_1000201762']
        #
        # spurious methionine
        # self.__dsList=['D_9999999997']
        # self.__dsList=['D_1000201769']
        # self.__dsList=['D_1000201492']
        # self.__dsList=['D_1000201681']
        # D_1000202201,/D_1000201218
        # Chenghua ribosome missing taxa search length problem
        # self.__dsList=['D_1000201218']
        # self.__dsList=['D_1000201820']
        # Chenghua ribosome with GB id problem --
        # self.__dsList=['D_1000201234']
        # Ezra reverse sense sequence
        # self.__dsList=['D_1000202354']
        # brian missing annotations
        # self.__dsList=['D_1000202211']
        #  another ribsome from Ezra requiring extended search for E.coli. RNA
        # self.__dsList=['D_1000202036']
        # PDBe issue with linker annotation
        # self.__dsList=['D_1000202812']
        # yuhe
        # self.__dsList=['D_1000202921']
        # Microhet - Buvna/Irina
        # self.__dsList=['D_1000202462']
        # add N at 182 in auth sequence
        # self.__dsList=['D_1000203231']
        # c/n terminal edits
        # self.__dsList=['D_1000203026']
        #
        # self.__dsList= ['D_1000000009']
        # self.__dsList=['D_1000203724']
        #
        # anomalous '.' in one-letter code sequence
        # self.__dsList=['D_1000203111']
        # self.__dsList=['D_1000201737']
        #
        # Marina/Ezra deletion issue with Valine
        # self.__dsList=['D_1000204532']
        #
        #        self.__dsList = ['D_1000202254']
        # self.__dsList = ['D_1000206113']
        # self.__dsList = ['D_1000205510']
        # Luigi crash on save - entity 3.
        # self.__dsList = ['D_1000206034']
        # self.__dsList = ['D_1000212844']
        # chimera bug
        # self.__dsList = ['D_8000200376']
        #
        # isoform
        # self.__dsList = ['D_1000215836']
        # isoform input isoform Q9D6K9-2 entity 3
        self.__dsList = ["D_1000214529"]
        #
        self.__maxRefAlign = 100
        #

    def __setup(self, sessionId="SEQUENCE-DATA-CACHE"):
        """Simulate the web application environment for managing session storage of
        temporary data files.
        """
        self.__siteId = getSiteId(defaultSiteId="WWPDB_DEPLOY_TEST")
        self.__cI = ConfigInfo(self.__siteId)
        self.__topPath = self.__cI.get("SITE_WEB_APPS_TOP_PATH")
        self.__topSessionPath = self.__cI.get("SITE_WEB_APPS_TOP_SESSIONS_PATH")
        #
        self.__reqObj = SeqModInputRequest({}, verbose=self.__verbose, log=self.__lfh)
        self.__templatePath = os.path.join(self.__topPath, "htdocs", "seqmodule")
        self.__reqObj.setValue("TopSessionPath", self.__topSessionPath)
        self.__reqObj.setValue("TemplatePath", self.__templatePath)
        self.__reqObj.setValue("TopPath", self.__topPath)
        self.__reqObj.setValue("WWPDB_SITE_ID", self.__siteId)
        os.environ["WWPDB_SITE_ID"] = self.__siteId

        if sessionId is not None:
            self.__reqObj.setValue("sessionid", sessionId)

        # self.__sessionId = self.__reqObj.getSessionId()
        self.__sessionObj = self.__reqObj.newSessionObj()
        self.__sessionPath = self.__sessionObj.getPath()
        self.__reqObj.printIt(ofh=self.__lfh)

        # self.__cI.dump(self.__lfh)

    def testSearchAndAssembleFromUpload(self):
        """Test search each entity sequence against appropriate reference sequence database
        service storing the matching sequences.

        Using upload file source. rcsb-mmcif
        """
        startTime = time.time()
        self.__lfh.write("\n\n========================================================================================================\n")
        self.__lfh.write("Starting %s %s at %s\n" % (self.__class__.__name__, inspect.currentframe().f_code.co_name, time.strftime("%Y %m %d %H:%M:%S", time.localtime())))

        try:
            self.__setup(sessionId="SEQUENCE-DATA-CACHE")
            pI = PathInfo(siteId=self.__siteId, sessionPath=self.__sessionPath, verbose=self.__verbose, log=self.__lfh)
            for f in self.__exampleIdList:

                (idCode, _fExt) = os.path.splitext(f)
                pdbxFilePath = pI.getModelPdbxFilePath(idCode, fileSource="session")
                inpFilePath = os.path.join(self.__pathExamples, f + ".cif")
                self.__lfh.write(
                    "+testSearchAndAssembleFromUpload() Starting with id %s \n + session path %s \n + input file %s \n   + pdbxfile %s\n"
                    % (idCode, self.__sessionPath, inpFilePath, pdbxFilePath)
                )
                dfa = DataFileAdapter(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
                _ok = dfa.modelConvertToPdbx(filePath=inpFilePath, fileType=self.__exampleFileType, pdbxFilePath=pdbxFilePath)  # noqa: F841
                self.__reqObj.setValue("identifier", idCode)
                #
                if not os.access(inpFilePath, os.F_OK):
                    self.__lfh.write("+testSearchAndAssembleFromUpload() input file missing %s\n" % inpFilePath)
                    self.fail()
                    break
                if not os.access(pdbxFilePath, os.F_OK):
                    self.__lfh.write("+testSearchAndAssembleFromUpload() format conversion failed for %s\n" % idCode)
                    self.fail()
                    break
                #
                sda = SequenceDataAssemble(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
                # sda.doAssemble(fileSource="local-upload")
                sda.doAssemble()
                alstat = AlignmentStatistics(reqObj=self.__reqObj, maxRefAlign=self.__maxRefAlign, verbose=self.__verbose, log=self.__lfh)
                alstat.doUpdate()
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write(
            "\nCompleted %s %s at %s (%.2f seconds)\n"
            % (self.__class__.__name__, inspect.currentframe().f_code.co_name, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        )

    #    def testSearchAndAssembleFromArchive(self):
    #        """ Test search each entity sequence against appropriate reference sequence database
    #            service storing the matching sequences.
    #
    #            Using archive file source.
    #        """
    #        startTime = time.time()
    #        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__,
    #                                                       inspect.currentframe().f_code.co_name,
    #                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime())))
    #        try:
    #            self.__setup(sessionId=None)
    #            for dsId in self.__dsList:
    #                self.__reqObj.setValue("identifier", dsId)
    #                #
    #                sda = SequenceDataAssemble(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
    #                sda.doAssemble(fileSource='archive')
    #                alstat = AlignmentStatistics(reqObj=self.__reqObj, maxRefAlign=self.__maxRefAlign, verbose=self.__verbose, log=self.__lfh)
    #                alstat.doUpdate()
    #                map = MultiAlignPseudo(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
    #                map.runSelected(identifier=dsId)
    #                #
    #                dI = DataImporter(reqObj=self.__reqObj, fileSource='wf-archive', verbose=self.__verbose, log=self.__lfh)
    #                ok = dI.copyFilesOnClose()
    #
    #        except:
    #            traceback.print_exc(file=self.__lfh)
    #            self.fail()
    #
    #        endTime = time.time()
    #        self.__lfh.write("\nCompleted %s %s at %s (%.2f seconds)\n" % (self.__class__.__name__,
    #                                                                       inspect.currentframe().f_code.co_name,
    #                                                                       time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
    #                                                                       endTime - startTime))

    def testSummaryView(self):
        """Test construction of a summary view from a sequence data store containing pairwise alignment stats."""
        startTime = time.time()
        self.__lfh.write("\nStarting %s %s at %s\n" % (self.__class__.__name__, inspect.currentframe().f_code.co_name, time.strftime("%Y %m %d %H:%M:%S", time.localtime())))
        try:
            self.__setup(sessionId="SEQUENCE-DATA-CACHE")
            self.__reqObj.setValue("identifier", self.__dsList[0])
            op = "load"
            sV = SummaryView(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            sumObj = sV.loadSummary(operation=op)

            sVD = SummaryViewDepiction(verbose=self.__verbose, log=self.__lfh)
            oL = sVD.buildSummaryView(sumObj)
            htmlPath = "current-alignment-summary.html"
            fp = open(htmlPath, "w")
            fp.write("%s" % "".join(oL))
            fp.close()
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

        endTime = time.time()
        self.__lfh.write(
            "\nCompleted %s %s at %s (%.2f seconds)\n"
            % (self.__class__.__name__, inspect.currentframe().f_code.co_name, time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime)
        )


def suiteSearchAndAssembleTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(SequenceDataPrepTests("testSearchAndAssembleFromArchive"))
    # suiteSelect.addTest(SequenceDataPrepTests("testSearchAndAssembleFromUpload"))
    return suiteSelect


def suiteSequenceDataPrepTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(SequenceDataPrepTests("testSummaryView"))
    return suiteSelect


if __name__ == "__main__":
    if True:  # pylint: disable=using-constant-test
        mySuite = suiteSearchAndAssembleTests()
        unittest.TextTestRunner(verbosity=2).run(mySuite)

    if False:  # pylint: disable=using-constant-test
        mySuite = suiteSequenceDataPrepTests()
        unittest.TextTestRunner(verbosity=2).run(mySuite)
