##
# File:    DataImporter.py
# Date:    25-Apr-2010
#
# Update:
# 25-Apr-2010  jdw Adapted from DataImportView
# 05-May-2010  jdw Adjust  handling of restart statistics data ...
#                  Leaving this hardwired with old non-checking method.
# 08-May-2010  jdw make polymer linkage table in rcsb mode
# 18-Dec-2012  jdw pass environment on to dependent functions
# 20-Feb-2013  jdw generalize to local file import
# 03-Mar-2013  jdw handle polymer linkage downstream
# 04-Mar-2013  jdw change prototype for AlignmentStatistics() class and verify this step.
# 07-Mar-2013  jdw add method to set log handle -
# 08-Mar-2013  jdw provide support for cached repository results.
# 03-Apr-2013  jdw add updater method at wf completion
#  1-May-2013  jdw add updateDataInstance()
# 20-Apr-2014  jdw add copyFilesOnStart() and copyFilesOnClose()
# 27-Apr-2014  jdw move selective logic for transfering files on open and close caller.
# 28-Apr-2014  jdw retire updateData* methods -
#  1-May-2014  jdw pre-check output files before copy --
#  8-May-2014  jdw add model file with naming conventions for JMol/JSMol
# 22-May-2014  jdw add pull/store alignment data files
# 26-May-2014  jdw precompute alignments on upload --
# 24-Aug-2017  zf  add copyFiles(), runBlastSearch(), __processArchiveOrInstanceCase(),
#                  replace __loadData() with loadSeqDataAssemble(),
#                  replace __archiveImport() with copyModelFile()
##

"""
Controlling class for data import - includes sequence alignment statistics generation -

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.07"


import sys
import os.path
import shutil
import glob
import traceback
from wwpdb.apps.seqmodule.control.SequenceDataAssemble_v2 import SequenceDataAssemble

from wwpdb.utils.dp.DataFileAdapter import DataFileAdapter
from wwpdb.io.locator.PathInfo import PathInfo


class DataImporter(object):
    """Controlling class for data import operations

    Supported file sources:
    + local-repository -  RCSB data repository or new archive repository
    + local-upload    -  Upload of Local CIF or CIFEPS data file -
    + archive         -  WF archive storage
    + wf-instance     -  WF instance storage

    """

    def __init__(self, reqObj=None, fileSource="local-upload", maxRefAlign=100, verbose=False, log=sys.stderr):  # pylint: disable=unused-argument
        self.__reqObj = reqObj
        self.__fileSource = fileSource
        self.__verbose = verbose
        self.__lfh = log
        # self.__maxRefAlign = maxRefAlign
        #
        self.__instance = ""
        self.__setup()
        #

    def __setup(self):
        """Extract the global input environment from the input request  -"""
        self.__sessionObj = self.__reqObj.getSessionObj()
        self.__sessionPath = self.__sessionObj.getPath()
        self.__siteId = self.__reqObj.getValue("WWPDB_SITE_ID")
        self.__operation = self.__reqObj.getValue("operation")
        self.__identifier = str(self.__reqObj.getValue("identifier")).upper()
        self.__instance = str(self.__reqObj.getValue("instance")).upper()
        self.__uploadFileName = str(self.__reqObj.getValue("UploadFileName"))
        self.__uploadFileType = str(self.__reqObj.getValue("UploadFileType"))
        self.__pI = PathInfo(siteId=self.__siteId, sessionPath=self.__sessionPath, verbose=self.__verbose, log=self.__lfh)
        self.__fileTypeList_wf_instance = [
            ["model", "pdbx", "model", "1", "next", True],
            ["seq-data-stats", "pic", "seq stats", "latest", "next", True],
            ["polymer-linkage-report", "html", "link report", "latest", "next", True],
            ["mismatch-warning", "pic", "mismatch warning", "latest", "latest", False],
        ]
        self.__fileTypeList_archive = [
            ["model", "pdbx", "model", "1", "next", True],
            ["seq-data-stats", "pic", "seq stats", "latest", "next", False],
            ["polymer-linkage-report", "html", "link report", "latest", "next", False],
        ]
        self.__seqFileTypeList = [["seq-align-data", "pic", "seq align", "latest", "next", False], ["seqdb-match", "pic", "seqdb match", "latest", "latest", False]]
        self.__polyLinkType = ["polymer-linkage-distances", "pdbx", "polymer linkage", "latest", "latest", False]
        self.__seqAssignType = ["seq-assign", "pdbx", "seq assign", "latest", "next", False]
        #
        if self.__verbose:
            self.__lfh.write("+DataImporter.__setup() session path    %s\n" % self.__sessionPath)
            self.__lfh.write("+DataImporter.__setup() file source     %s\n" % self.__fileSource)
            self.__lfh.write("+DataImporter.__setup() identifier      %s\n" % self.__identifier)
            self.__lfh.write("+DataImporter.__setup() instance        %s\n" % self.__instance)
            self.__lfh.write("+DataImporter.__setup() operation       %s\n" % self.__operation)
            self.__lfh.write("+DataImporter.__setup() UploadFileType  %s\n" % self.__uploadFileType)
            self.__lfh.write("+DataImporter.__setup() UploadFileName  %s\n" % self.__uploadFileName)
        #

    def copyFilesOnStart(self, includePolyLinkFile=False, includeSeqAssignFile=False):
        """Manage moving selected data files from archive or instance storage to session storage

        fileSource         is the target destination (archive|wf-instance)

        Return: True on success or False otherwise
        """
        ok = self.copyFiles(
            inputFileSource=self.__fileSource,
            inputWfInstanceId=self.__instance,
            includeModelFile=True,
            includePolyLinkFile=includePolyLinkFile,
            includeSeqAssignFile=includeSeqAssignFile,
            checkFileAvailability=True,
            messageHead="DataImporter.copyFilesOnStart()",
        )
        #
        pdbxModelFilePath = self.__pI.getModelPdbxFilePath(self.__identifier, fileSource="session", versionId="1")
        self.__copyJMolFile(pdbxModelFilePath)
        #
        return ok

    def copyFilesOnClose(self, includePolyLinkFile=False, includeSeqAssignFile=False):
        """Manage moving selected data files from the current session to archive or instance storage.

        fileSource         is the target destination (archive|wf-instance)

        Return: True on success or False otherwise

        """
        return self.copyFiles(
            inputFileSource="session",
            outputFileSource=self.__fileSource,
            outputWfInstanceId=self.__instance,
            includeModelFile=True,
            includePolyLinkFile=includePolyLinkFile,
            includeSeqAssignFile=includeSeqAssignFile,
            versionIndex=4,
            messageHead="DataImporter.copyFilesOnClose()",
        )

    def copyFiles(
        self,
        inputFileSource="archive",
        inputWfInstanceId=None,
        outputFileSource="session",
        outputWfInstanceId=None,
        versionIndex=3,
        includeModelFile=False,
        includePolyLinkFile=False,
        includeSeqAssignFile=False,
        checkFileAvailability=False,
        entityIdList=None,
        messageHead="DataImporter.copyFiles()",
    ):
        """ """
        if entityIdList is None:
            entityIdList = []
        startIndex = 1
        if includeModelFile:
            startIndex = 0
        #
        if (inputFileSource == "wf-instance") or (outputFileSource == "wf-instance"):
            fileTypeList = self.__fileTypeList_wf_instance[startIndex:]
        else:
            fileTypeList = self.__fileTypeList_archive[startIndex:]
        #
        returnOk = True
        afpL = []
        for fileType in fileTypeList:
            inputFilePath = self.__pI.getFilePath(self.__identifier, wfInstanceId=inputWfInstanceId, contentType=fileType[0], formatType=fileType[1], fileSource=inputFileSource)
            if not os.access(inputFilePath, os.R_OK):
                if checkFileAvailability and fileType[5]:
                    returnOk = False
                    self.__lfh.write("+%s missing input %s file %s\n" % (messageHead, fileType[2], inputFilePath))
                #
                continue
            #
            outputFilePath = self.__pI.getFilePath(
                self.__identifier, wfInstanceId=outputWfInstanceId, contentType=fileType[0], formatType=fileType[1], fileSource=outputFileSource, versionId=fileType[versionIndex]
            )
            afpL.append((inputFilePath, outputFilePath, fileType[2]))
        #
        if checkFileAvailability and (not returnOk):
            return False
        #
        if includePolyLinkFile:
            inputFilePath = self.__pI.getFilePath(
                self.__identifier, wfInstanceId=inputWfInstanceId, contentType=self.__polyLinkType[0], formatType=self.__polyLinkType[1], fileSource=inputFileSource
            )
            if os.access(inputFilePath, os.R_OK):
                outputFilePath = self.__pI.getFilePath(
                    self.__identifier,
                    wfInstanceId=outputWfInstanceId,
                    contentType=self.__polyLinkType[0],
                    formatType=self.__polyLinkType[1],
                    fileSource=outputFileSource,
                    versionId=self.__polyLinkType[versionIndex],
                )
                afpL.append((inputFilePath, outputFilePath, self.__polyLinkType[2]))
            #
        #
        if includeSeqAssignFile:
            inputFilePath = self.__pI.getFilePath(
                self.__identifier, wfInstanceId=inputWfInstanceId, contentType=self.__seqAssignType[0], formatType=self.__seqAssignType[1], fileSource=inputFileSource
            )
            if os.access(inputFilePath, os.R_OK):
                outputFilePath = self.__pI.getFilePath(
                    self.__identifier,
                    wfInstanceId=outputWfInstanceId,
                    contentType=self.__seqAssignType[0],
                    formatType=self.__seqAssignType[1],
                    fileSource=outputFileSource,
                    versionId=self.__seqAssignType[versionIndex],
                )
                afpL.append((inputFilePath, outputFilePath, self.__seqAssignType[2]))
            #
        #
        if entityIdList:
            for eId in entityIdList:
                for fileType in self.__seqFileTypeList:
                    inputFilePath = self.__pI.getFilePath(
                        self.__identifier, wfInstanceId=inputWfInstanceId, contentType=fileType[0], formatType=fileType[1], fileSource=inputFileSource, partNumber=eId
                    )
                    if not os.access(inputFilePath, os.R_OK):
                        continue
                    #
                    outputFilePath = self.__pI.getFilePath(
                        self.__identifier,
                        wfInstanceId=outputWfInstanceId,
                        contentType=fileType[0],
                        formatType=fileType[1],
                        fileSource=outputFileSource,
                        versionId=fileType[versionIndex],
                        partNumber=eId,
                    )
                    afpL.append((inputFilePath, outputFilePath, fileType[2]))
                #
            #
        else:
            gId = 1
            while True:
                inputFilePath = self.__pI.getFilePath(
                    self.__identifier,
                    wfInstanceId=inputWfInstanceId,
                    contentType=self.__seqFileTypeList[1][0],
                    formatType=self.__seqFileTypeList[1][1],
                    fileSource=inputFileSource,
                    partNumber=str(gId),
                )
                if os.access(inputFilePath, os.R_OK):
                    outputFilePath = self.__pI.getFilePath(
                        self.__identifier,
                        wfInstanceId=outputWfInstanceId,
                        contentType=self.__seqFileTypeList[1][0],
                        formatType=self.__seqFileTypeList[1][1],
                        fileSource=outputFileSource,
                        versionId=self.__seqFileTypeList[1][versionIndex],
                        partNumber=str(gId),
                    )
                    afpL.append((inputFilePath, outputFilePath, self.__seqFileTypeList[1][2]))
                #
                inputFilePath = self.__pI.getFilePath(
                    self.__identifier,
                    wfInstanceId=inputWfInstanceId,
                    contentType=self.__seqFileTypeList[0][0],
                    formatType=self.__seqFileTypeList[0][1],
                    fileSource=inputFileSource,
                    partNumber=str(gId),
                )
                if os.access(inputFilePath, os.R_OK):
                    outputFilePath = self.__pI.getFilePath(
                        self.__identifier,
                        wfInstanceId=outputWfInstanceId,
                        contentType=self.__seqFileTypeList[0][0],
                        formatType=self.__seqFileTypeList[0][1],
                        fileSource=outputFileSource,
                        versionId=self.__seqFileTypeList[0][versionIndex],
                        partNumber=str(gId),
                    )
                    afpL.append((inputFilePath, outputFilePath, self.__seqFileTypeList[0][2]))
                    gId += 1
                else:
                    break
                #
            #
        #
        if not afpL:
            return False
        #
        if self.__verbose:
            for ftuple in afpL:
                self.__lfh.write("+%s input %s file %s\n" % (messageHead, ftuple[2], ftuple[0]))
            #
            for ftuple in afpL:
                self.__lfh.write("+%s export %s file %s\n" % (messageHead, ftuple[2], ftuple[1]))
            #
        #
        try:
            for ftuple in afpL:
                shutil.copyfile(ftuple[0], ftuple[1])
            #
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+%s data file update FAILED\n" % messageHead)
                traceback.print_exc(file=self.__lfh)
            #
            returnOk = False
        #
        return returnOk

    def copyModelFile(self, inputFileSource="archive", inputWfInstanceId=None, outputFileSource="session", outputWfInstanceId=None, versionIndex=3):
        """Copies target model file. Returns the imported model file Name or None"""
        try:
            fileType = self.__fileTypeList_wf_instance[0]
            inputFilePath = self.__pI.getFilePath(self.__identifier, wfInstanceId=inputWfInstanceId, contentType=fileType[0], formatType=fileType[1], fileSource=inputFileSource)
            if not os.access(inputFilePath, os.R_OK):
                return None
            #
            outputFilePath = self.__pI.getFilePath(
                self.__identifier, wfInstanceId=outputWfInstanceId, contentType=fileType[0], formatType=fileType[1], fileSource=outputFileSource, versionId=fileType[versionIndex]
            )
            if not os.access(outputFilePath, os.R_OK):
                shutil.copyfile(inputFilePath, outputFilePath)
            #
            return outputFilePath
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            return None
        #

    def setLogHandle(self, log=sys.stderr):
        """Reset the stream for logging output."""
        try:
            self.__lfh = log
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            return False

    def __testCache(self, identifier):
        """Return the path to file cache the cache directory contains any target files -"""
        cP = os.path.join(self.__reqObj.getSessionPath(), "SEQUENCE-DATA-CACHE")
        if os.access(cP, os.R_OK):
            pat = cP + "/" + identifier.upper() + "_*"
            pathList = glob.glob(pat)
            if len(pathList) > 0:
                return cP
        return None

    def loadData(self, useCache=True):
        """Using input data from the input request file source, assemble and persist data into
        the storage system used by the sequence editor components.
        """
        #
        if useCache:
            cachePath = self.__testCache(self.__identifier)
        else:
            cachePath = None
        #

        if self.__verbose:
            self.__lfh.write("+DataImporter.loadData() with file source %r : session path %s cache path %s\n" % (self.__fileSource, self.__sessionPath, cachePath))
            self.__lfh.flush()
        #
        # if ((cachePath is not None)  and (self.__fileSource == "local-repository")):
        #     entityIdList,ok=self.loadSeqDataAssemble(cachePath=cachePath)
        #     return ok
        # elif self.__fileSource == "local-repository":
        if self.__fileSource == "local-repository":
            #
            # Here we have a repository input identifier (e.g. RCSB ID or D_xxxxxx) -
            #
            idCode = self.__identifier
            if idCode.lower().startswith("rcsb"):
                modelFileName = self.__repositoryImport(idCode)  # pylint: disable=assignment-from-none
                if modelFileName is None:
                    return False
                #
                filePath = os.path.join(self.__sessionPath, modelFileName)
                pdbxFilePath = self.__pI.getModelPdbxFilePath(idCode, fileSource="session", versionId="1")
                dfa = DataFileAdapter(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
                ok = dfa.modelConvertToPdbx(filePath=filePath, fileType="rcsb-mmcif", pdbxFilePath=pdbxFilePath)
                if not ok:
                    return False
                #
                self.__copyJMolFile(pdbxFilePath)
                _entityIdList, ok = self.loadSeqDataAssemble()
                return ok
            elif idCode.upper().startswith("D_"):
                # reset file source to archive --
                self.__fileSource = "archive"
                self.__instance = ""
                return self.__processArchiveOrInstanceCase()
            else:
                return False
            #
        elif self.__fileSource == "local-upload":
            #
            # Here we have a copy of the input model file in the current session directory -
            #
            idCode = self.__identifier
            modelFileName = self.__uploadFileName
            filePath = os.path.join(self.__sessionPath, modelFileName)
            pdbxFilePath = self.__pI.getModelPdbxFilePath(idCode, fileSource="session", versionId="1")
            if self.__verbose:
                self.__lfh.write("+DataImporter.loadData() idCode %s filePath %s pdbxFilePath %s fileType %s\n" % (idCode, filePath, pdbxFilePath, self.__uploadFileType))
            #
            dfa = DataFileAdapter(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            ok = dfa.modelConvertToPdbx(filePath=filePath, fileType=self.__uploadFileType, pdbxFilePath=pdbxFilePath)
            if self.__verbose:
                self.__lfh.write("+DataImporter.loadData() pdbx convert status  %r\n" % ok)
            if not ok:
                return False
            #
            self.__copyJMolFile(pdbxFilePath)
            _entityIdList, ok = self.loadSeqDataAssemble(doRefSearch=True)
            return ok
        elif self.__fileSource in ["archive", "wf-archive", "wf-instance"]:
            return self.__processArchiveOrInstanceCase()
        elif self.__fileSource in ["sequence-database", "input-taxid"]:
            # place holder for data update.
            return True
        else:
            return False
        #

    def loadSeqDataAssemble(self, calPolyLink=True, doRefSearch=False, doAutoProcess=False):
        if self.__verbose:
            self.__lfh.write(
                "\n+DataImporter.loadSeqDataAssemble() RECALCULATING INITIAL DATA STATE calPolyLink=%r doRefSearch=%r doAutoProcess=%r\n"
                % (calPolyLink, doRefSearch, doAutoProcess)
            )
        #
        try:
            if self.__verbose:
                self.__lfh.write("+DataImporter.loadSeqDataAssemble() using file source %s identifier %s instance %s\n" % (self.__fileSource, self.__identifier, self.__instance))
            #
            self.__lfh.flush()

            sda = SequenceDataAssemble(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
            sda.doAssemble(calPolyLink=calPolyLink, doRefSearch=doRefSearch, doAutoProcess=doAutoProcess)
            #
            return sda.getEntityIdList(), True
        except:  # noqa: E722 pylint: disable=bare-except
            self.__lfh.write("+DataImporter.loadSeqDataAssemble() fails with file source %s identifier %s instance %s\n" % (self.__fileSource, self.__identifier, self.__instance))
            traceback.print_exc(file=self.__lfh)

        return [], False

    def __repositoryImport(self, idCode):  # pylint: disable=unused-argument
        """Copies target model file from repository to the current session directory.

        Returns the imported model file Name or None
        """
        # Not used anymore
        return None

    def __processArchiveOrInstanceCase(self):
        if self.__verbose:
            self.__lfh.write("\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
            self.__lfh.write("+DataImporter.loadData() loading curent data from a workflow file source %s\n\n" % self.__fileSource)

        # First test if we can restart -
        #
        includePolyLinkFile = False
        if self.__fileSource == "wf-instance":
            includePolyLinkFile = True
        #
        ok = self.copyFilesOnStart(includePolyLinkFile=includePolyLinkFile)
        if ok and self.__fileSource == "wf-instance":
            if self.__verbose:
                self.__lfh.write(
                    "\n+DataImporter.__processArchiveOrInstanceCase() SUCCESSFULLY COPIED STARTING DATA from file source %s with status code %r\n" % (self.__fileSource, ok)
                )
            return ok

        if not ok:
            modelFileName = self.copyModelFile()
            if modelFileName is None:
                return False
                #
            self.__copyJMolFile(modelFileName)
            #
        _entityIdList, ok = self.loadSeqDataAssemble()
        return ok

    def __copyJMolFile(self, pdbxModelFilePath):
        pdbxJMolFilePath = os.path.join(self.__sessionPath, self.__identifier + ".cif")
        if os.access(pdbxModelFilePath, os.R_OK):
            shutil.copyfile(pdbxModelFilePath, pdbxJMolFilePath)
        #


if __name__ == "__main__":
    di = DataImporter()
