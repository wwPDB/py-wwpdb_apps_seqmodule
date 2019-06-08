##
# File:  SequenceDataAssemble.py
# Date:  25-Feb-2013
#
# Updates:
#  26-Feb-2013 jdw  assign polymer parts to sequences and features on input
#  03-Mar-2013 jdw  create rendered table of long intra-polymer linkages
#  03-Mar-2013 jdw  upldate new part data items in sequence data store.
#  05-Mar-2013 jdw  nuke and pave
#  29-Nov-2013 jdw  save model file name in UtilDataStore() for later use by visualization app.
#   3-Dec-2013 jdw  skip PDB format data file production
#  19-Jan-2014 jdw  add 'HOST_ORG_CELL_LINE'
#  19-Jan-2014 jdw  standardize path polymer link report using getPolyLinkReportFilePath()
#  01-Feb-2015 jdw  add host org variant
#  24-Aug-2017  zf  remove __doAssembleFromArchive() & __doAssembleFromModel(), move all file copy operations into DataImporter.py
#                   modify __calcPolymerLinkages(), also modify backend C++ code to speed up file process
#                   modify __doReferenceSearch() to use previous blast results if existing & applicable
##
"""
Assemble sequence and other required data to start/restart the sequence editor tool.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.09"

import sys
import os
import shutil
import time
import traceback

from wwpdb.apps.seqmodule.io.SequenceDataStore import SequenceDataStore
from wwpdb.apps.seqmodule.util.SequenceReferenceData import SequenceReferenceData
from wwpdb.apps.seqmodule.util.SequenceLabel import SequenceLabel, SequenceFeature

from wwpdb.apps.seqmodule.io.PdbxIoUtils import PdbxFileIo, ModelFileIo, ReferenceSequenceIo, PolymerLinkageIo, PolymerInstanceIo

from wwpdb.apps.seqmodule.io.ModelSequenceUtils import ModelSequenceUtils
from wwpdb.apps.seqmodule.io.ReferenceSequenceUtils import ReferenceSequenceUtils

from wwpdb.apps.seqmodule.link.PolymerLinkageDepict import PolymerLinkageDepict

from wwpdb.io.locator.PathInfo import PathInfo
from wwpdb.utils.dp.RcsbDpUtility import RcsbDpUtility
from wwpdb.utils.session.UtilDataStore import UtilDataStore


class SequenceDataAssemble(object):

    """
     This class encapsulates the data assembly operations
     of model sequence and reference data matching data for the sequence
     editor tool.

     Storage model - imported data is loaded into the sequence data store
                     where it is managed by the SequenceDataStore() class.

     Methods in this class extract the author and coordinate sequence data
     and source information from PDBx CIF.

     Reference sequence data is extracted from the processed BLAST results
     (ie. wwPDB seqdb-match data files) for each polymer entity.


    """

    def __init__(self, reqObj=None, verbose=False, log=sys.stderr):

        self.__reqObj = reqObj
        self.__verbose = verbose
        self.__lfh = log
        self.__debug = True

        self.__sessionObj = self.__reqObj.getSessionObj()
        self.__sessionPath = self.__sessionObj.getPath()
        self.__siteId = self.__reqObj.getValue("WWPDB_SITE_ID")
        self.__pI = PathInfo(siteId=self.__siteId, sessionPath=self.__sessionPath, verbose=self.__verbose, log=self.__lfh)
        #
        # minimum sequence length to search in reference databases -
        #
        self.__minSearchSequenceLengthAA = 12
        self.__minSearchSequenceLengthNA = 50
        self.__entityIdList = []
        #

    def doAssemble(self, calPolyLink=True, selectedEntityIdList=[], doRefSearch=False):
        """  Perform all steps required to construct the session sequence data store using data from within
             the input 'fileSource'.
        """
        dataSetId = str(self.__reqObj.getValue("identifier")).upper()
        #
        if (self.__verbose):
            self.__lfh.write("\n+SequenceDataAssemble.doAssemble() STARTING with depId %r calPolyLink %r selectedEntityIdList %r doRefSearch %r\n" %
                             (dataSetId, calPolyLink, selectedEntityIdList, doRefSearch))
        #
        pdbxFilePath = self.__pI.getModelPdbxFilePath(dataSetId, fileSource='session', versionId="1")
        polyLinkPath = self.__pI.getPolyLinkFilePath(dataSetId, fileSource='session')
        #
        # Calculate polymer linkage information
        #
        if calPolyLink or (not os.access(polyLinkPath, os.R_OK)):
            self.__calcPolymerLinkages(pdbxFilePath=pdbxFilePath, polyLinkPath=polyLinkPath)
        #
        entryD, entityD, instanceD, depSeqAssign, seqAssign = self.__importModel(dataSetId, pdbxFilePath=pdbxFilePath, polyLinkPath=polyLinkPath)
        self.__entityIdList = entityD.keys()
        #
        if selectedEntityIdList:
            entityIdList = selectedEntityIdList
            forceSearch = True
        else:
            entityIdList = self.__entityIdList
            forceSearch = doRefSearch
        #
        self.__doReferenceSearch(dataSetId, entityD, entityIdList, forceSearch, minSearchSequenceLengthAA=self.__minSearchSequenceLengthAA, \
                                 minSearchSequenceLengthNA=self.__minSearchSequenceLengthNA)
        #
        # Read reference sequences -
        #
        refStatus, eRefD = self.__readReferenceSearchResults(dataSetId, entityD)
        #
        if (self.__verbose):
            self.__lfh.write("+SequenceDataAssemble.doAssemble() - Reference read status     %r\n" % refStatus)
            self.__lfh.write("+SequenceDataAssemble.doAssemble() - Reference match length    %d\n" % len(eRefD))
        #
        self.__updateUtilDataStore(dataSetId, pdbxFilePath)
        #
        self.__loadSequenceDataStore(dataSetId, entryD, entityD, instanceD, eRefD, depSeqAssign, seqAssign)

        if (self.__verbose):
            self.__lfh.write("+SequenceDataAssemble.doAssemble() COMPLETED for %r\n\n" % dataSetId)
        return True

    def getEntityIdList(self):
        return self.__entityIdList

    def __calcPolymerLinkages(self, pdbxFilePath=None, polyLinkPath=None):
        #
        # Calculate intra-polymer linkage distances -
        #
        try:
            dp = RcsbDpUtility(tmpPath=self.__sessionPath, siteId=self.__siteId,verbose=self.__verbose,log=self.__lfh)
            dp.imp(pdbxFilePath)
            dp.op("annot-poly-link-dist")
            dp.exp(polyLinkPath)
            dp.cleanup()
            if (self.__verbose):
                self.__lfh.write("+SequenceDataAssemble.__calcPolymerLinkages() - Polymer link distance file copied to: %s\n" % polyLinkPath)
            return True
        except:
            if (self.__verbose):
                self.__lfh.write("+SequenceDataAssemble.__calcPolymerLinkages() - Polymer link distance calculation failed for file %s\n" % pdbxFilePath)
                traceback.print_exc(file=self.__lfh)
            return False

    def __importModel(self, dataSetId, pdbxFilePath=None, polyLinkPath=None):
        """  Return entry, entity, instance detail dictionary objects containing content extracted
             from the input model data file.
        """
        instanceD = {}
        polymerLinkDistList = []
        #
        # Read input model file -
        #
        if os.access(polyLinkPath, os.R_OK):
            c0 = PdbxFileIo(verbose=self.__verbose, log=self.__lfh).getContainer(polyLinkPath)
            pl = PolymerLinkageIo(dataContainer=c0, verbose=self.__verbose, log=self.__lfh)
            polymerLinkDistList =pl.getPolymerLinkDistances()
            #
            pi = PolymerInstanceIo(dataContainer=c0, verbose=self.__verbose, log=self.__lfh)
            instanceD = pi.getPolymerInstances()
        else:
            c0 = PdbxFileIo(verbose=self.__verbose, log=self.__lfh).getContainer(pdbxFilePath)
        #
        # Render polymer linkage information
        #
        polyLinkHtmlPath = self.__pI.getPolyLinkReportFilePath(dataSetId, fileSource='session')
        self.__renderPolymerLinkages(polymerLinkDistList=polymerLinkDistList, polyLinkHtmlPath=polyLinkHtmlPath)
        #
        #   Model identifier -
        #
        mfIo = ModelFileIo(dataContainer=c0, verbose=self.__verbose, log=self.__lfh)
        entryD = {}
        entryD['PDB_ID'] = mfIo.getDbCode('PDB')
        entryD['PDBX_MODEL_FILE_PATH'] = pdbxFilePath
        entryD['STRUCT_TITLE'] = mfIo.getStructTitle()
        entryD['CITATION_TITLE'] = mfIo.getCitationTitle()
        #
        #  Extract entity and coordinate sequence details from the model
        #
        msu = ModelSequenceUtils(dataContainer=c0, verbose=self.__verbose, log=self.__lfh)
        entityD = msu.getEntitySequenceDetails()
        if not instanceD:
            instanceD = msu.getCoordinateSequenceDetails(linkInstD=linkInstD)
        #
        depSeqAssign, seqAssign = msu.getSequenceAssignmentDetails()
        #
        if (self.__verbose):
            self.__lfh.write("+SequenceDataAssemble.__importModel() - session path              %s\n" % self.__sessionPath)
            self.__lfh.write("+SequenceDataAssemble.__importModel() - data set id               %s\n" % dataSetId)
            self.__lfh.write("+SequenceDataAssemble.__importModel() - PDBx file                 %s\n" % pdbxFilePath)
            self.__lfh.write("+SequenceDataAssemble.__importModel() - polymer linkage file      %s\n" % polyLinkPath)
            self.__lfh.write("+SequenceDataAssemble.__importModel() - polymer linkage HTML file %s\n" % polyLinkHtmlPath)
        #
        return entryD, entityD, instanceD, depSeqAssign, seqAssign

    def __renderPolymerLinkages(self, polymerLinkDistList, polyLinkHtmlPath):
        """  Create HTML table rendering of the input polymer linkage distance list.

             HTML output is written to the 'polyLinkHtmlPath':

             Returns True for success of False otherwise.
        """
        try:
            pld = PolymerLinkageDepict(upperBound=1.75, lowerBound=1.2, verbose=self.__verbose, log=self.__lfh)
            oL = pld.buildPolymerLinkageTable(polymerLinkDistList)
            ofh = open(polyLinkHtmlPath, 'w')
            ofh.write("%s" % "".join(oL))
            ofh.close()
            return True
        except:
            self.__lfh.write("+SequenceDataAssemble.__renderPolymerLinkage() - length linkage list %d  html path %s\n" %
                             (len(oL), polyLinkHtmlPath))
            traceback.print_exc(file=self.__lfh)
        return False

    def __doReferenceSearch(self, dataSetId, entityD, entityIdList, forceSearch, minSearchSequenceLengthAA=20, minSearchSequenceLengthNA=50):
        """  Perform the reference sequence database search using the input entity dictionary.
             Store matching results in the local session directory.
        """
        try:
            startTime = time.clock()
            #
            sdr = SequenceReferenceData(self.__verbose, self.__lfh)
            rsu = ReferenceSequenceUtils(siteId=self.__siteId, verbose=self.__verbose, log=self.__lfh)
            rsio = ReferenceSequenceIo(verbose=self.__verbose, log=self.__lfh)
            for eId in entityIdList:
                if not entityD.has_key(eId):
                    continue
                #
                eD = entityD[eId]
                polyTypeCode = sdr.getPolymerTypeCode(eD['POLYMER_LINKING_TYPE'])
                seqLen = len(eD['SEQ_ENTITY_1_CAN'])
                #
                skip = ((polyTypeCode in ['SAC']) or ((polyTypeCode in ['DNA', 'RNA', 'XNA']) and (seqLen < minSearchSequenceLengthNA)) or \
                       ((polyTypeCode in ['AA']) and (seqLen < minSearchSequenceLengthAA)))
                if (self.__verbose):
                    self.__lfh.write("+SequenceDataAssemble.__doReferenceSearch() search for entity id %s type %s length %d skip status %r\n" %
                                     (eId, polyTypeCode, seqLen, skip))
                    self.__lfh.flush()
                #
                if skip:
                    continue
                #
                fn = self.__pI.getReferenceSequenceFilePath(dataSetId, entityId=eId, fileSource='session')
                if (not forceSearch) and os.access(fn, os.R_OK):
                    if (self.__verbose):
                        self.__lfh.write("+SequenceDataAssemble.__doReferenceSearch() found previous ReferenceSequenceFile %s for entity id %s\n" % (fn, eId))
                    #
                    currSeq = eD['SEQ_ENTITY_1']
                    currSeq = currSeq.replace(' ', '').replace('\n', '').replace('\t', '')
                    #
                    rc0 = PdbxFileIo(verbose=self.__verbose, log=self.__lfh).getContainer(fn)
                    if rc0.exists('info'):
                        table = rc0.getObj('info')
                        prevSeq = table.getValue('sequence', 0)
                        prevSeq = prevSeq.replace(' ', '').replace('\n', '').replace('\t', '')
                        if currSeq == prevSeq:
                            if (self.__verbose):
                                self.__lfh.write("+SequenceDataAssemble.__doReferenceSearch() using previous search result for entity id %s\n" % eId)
                            continue
                        #
                    #
                #
                mR = rsu.searchEntities(entityD=eD, saveBlast=True, filePath=self.__sessionPath, filePrefix=dataSetId)
                rsio.writeMatchResults(eD, outFilePath=fn, matchResults=mR)

            if (self.__verbose):
                endTime = time.clock()
                self.__lfh.write("+SequenceDataAssemble.__doReferenceSearch()  completed at %s (%.2f seconds)\n" %
                                 (time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime))
            #
            return True
        except:
            self.__lfh.write("+SequenceDataAssemble.__doReferenceSearch() - failing \n")
            traceback.print_exc(file=self.__lfh)
            return False
        #

    def __readReferenceSearchResults(self, dataSetId, entityD):
        """  Read the reference sequence database search results for entities input entity dictionary.
             Return (status,eRefD) completion status and a dictionary of reference matching details.
        """
        eRefD = {}
        #
        try:
            for eId, eD in entityD.items():
                fn = self.__pI.getReferenceSequenceFilePath(dataSetId, entityId=eId, fileSource='session')
                if (os.access(fn, os.R_OK)):
                    rc0 = PdbxFileIo(verbose=self.__verbose, log=self.__lfh).getContainer(fn)
                    rsin = ReferenceSequenceIo(dataContainer=rc0, verbose=self.__verbose, log=self.__lfh)
                    rL = rsin.readMatchResults()
                    eRefD[eId] = rL
                else:
                    if (self.__verbose):
                        self.__lfh.write("+SequenceDataAssemble.__readReferenceSearchResults() - no reference file for entity %s in path %s\n" % (eId, fn))
                    eRefD[eId] = []
            return True, eRefD
        except:
            self.__lfh.write("+SequenceDataAssemble.__readReferenceSearchResults() - failing \n")
            traceback.print_exc(file=self.__lfh)
            return False, eRefD
        #

    def __updateUtilDataStore(self, dataSetId, pdbxFilePath):
        uds = UtilDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        pdbxFileName = os.path.basename(pdbxFilePath)
        uds.set('modelfilename', pdbxFileName)
        uds.serialize()
        #
        # Make an extra copy of the CIF file with a standard file extension  --  (may not be needed)
        #
        dstPath = os.path.join(self.__sessionPath, dataSetId + ".cif")
        self.__copyFile(pdbxFilePath, dstPath)

    def __loadSequenceDataStore(self, dataSetId, entryD, entityD, instanceD, eRefD, depSeqAssign, seqAssign):
        """ Store sequence data in a SequenceDataStore repository in the current session directory.

            Inputs sources -- entityD[entityId]      = ModelSequenceUtils.getEntitySequenceDetails()
                              instanceD[instanceID]  = ModelSequenceUtils.getCoordinateSequenceDetails()
                              eRefD[entityID]        = ReferenceSequenceIo.readMatchResults()
                              entryD[ENTRY_ID,PDB_ID]
                              depSeqAssign[entityId] = depositor reference sequence assignments -
                              seqAssign[entity]      = archive reference sequence assignments
        """
        #
        startTime = time.clock()
        #
        if (self.__verbose):
            self.__lfh.write("+SequenceDataAssemble.__loadSequenceDataStore() for sessionId %s started at %s \n" %
                             (self.__sessionObj.getId(), time.strftime("%Y %m %d %H:%M:%S", time.localtime())))
        #
        sdr = SequenceReferenceData(self.__verbose, self.__lfh)
        groupDict = {}
        groupPartD = {}
        polymerTypeCode = {}
        for id, eD in entityD.items():
            polyType = eD['POLYMER_LINKING_TYPE']
            polyTypeCode = sdr.getPolymerTypeCode(polyType)
            groupDict[id] = eD['INSTANCE_LIST']

            groupPartD[id] = range(1, len(eD['PART_LIST']) + 1)
            if self.__verbose:
                self.__lfh.write("+SequenceDataAssemble.__loadSequenceDataStore() entity %s PART_LIST %r\n" % (id, eD['PART_LIST']))
                self.__lfh.write("+SequenceDataAssemble.__loadSequenceDataStore() entity %s numParts %d groupPartD  %r\n" %
                                 (id, len(eD['PART_LIST']), groupPartD[id]))
            for inst in eD['INSTANCE_LIST']:
                polymerTypeCode[inst] = polyTypeCode
                if (self.__verbose):
                    self.__lfh.write("+SequenceDataAssemble.__loadSequenceDataStore()  instance %s  type %s\n" % (inst, polymerTypeCode[inst]))
        #
        seqFeature = SequenceFeature(verbose=self.__verbose, log=self.__lfh)
        sds = SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
        sds.reset()
        #
        #
        sds.setDepositorReferenceAssignments(assignD=depSeqAssign)
        sds.setReferenceAssignments(assignD=seqAssign)

        #
        # Load coordinate sequences & features  (coordinate sequences have no parts i.e. partId=1)
        #
        if ((entryD['PDB_ID'] is not None) and (len(entryD['PDB_ID']))) == 4:
            pdbId = entryD['PDB_ID']
        else:
            pdbId = dataSetId

        for id, S3L in instanceD.items():
            #
            # Only store sequences if we know the polymer type details -
            if id not in polymerTypeCode:
                continue
            sds.setSequence(S3L, id, 'xyz', partId=1, altId=1, version=1)
            seqFeature.clear()
            seqFeature.setId(dbName='PDB', dbCode=pdbId, dbAccession=pdbId)
            seqFeature.setPolymerType(polymerTypeCode[id])
            #
            sds.setFeature(seqFeature.get(), id, 'xyz', partId=1, altId=1, version=1)

        #
        #
        # Create intermediate dictionary of source part details --
        #
        #   Add sequence annotations of expression tag and linker if these are
        #   defined within entity source category --
        #
        ePartD = {}
        for eId, eD in entityD.items():
            seqBegMin = 100
            seqEndMax = 0
            rList = []
            for pD in eD['PART_LIST']:
                seqNumBeg = pD['SEQ_NUM_BEG']
                seqNumEnd = pD['SEQ_NUM_END']
                seqPartId = pD['SEQ_PART_ID']
                tpT = pD['SEQ_PART_TYPE']
                pT = tpT.lower()
                spt = ''
                if pT in ['biological sequence']:
                    spt = ''
                elif pT in ['n-terminal tag', 'c-terminal tag']:
                    spt = 'expression tag'
                elif pT in ['linker']:
                    spt = 'linker'
                seqBegMin = min(seqBegMin, seqNumBeg)
                seqEndMax = max(seqEndMax, seqNumEnd)
                rList.append((seqNumBeg, seqNumEnd, spt))
            #
            ePartD[eId] = (rList, seqBegMin, seqEndMax)
        #
        # ---------------------
        # Load author/entity sequence and features  - for multipart entities, the full residue sequence
        #  is loaded for each part. Distinguishing sequence features are stored for each part.
        #

        for eId, eD in entityD.items():
            sId = eD['IDENTIFIER']
            polyType = eD['POLYMER_LINKING_TYPE']
            polyTypeCode = sdr.getPolymerTypeCode(polyType)
            seqS = eD['SEQ_ENTITY_1']
            (r1L, r3L) = sdr.parseSequence(seqS, polyTypeCode)
            sTupL = []
            ir = 1
            for r3 in r3L:
                comment = self.__getAuthSeqComment(ir, ePartD[eId])
                sTupL.append((r3, str(ir), comment, ir))
                ir += 1
            #
            pDList = eD['PART_LIST']
            for (partNo, pD) in enumerate(pDList, start=1):
                seqNumBeg = pD['SEQ_NUM_BEG']
                seqNumEnd = pD['SEQ_NUM_END']
                seqPartId = pD['SEQ_PART_ID']
                seqPartType = pD['SEQ_PART_TYPE']

                sds.setSequence(sTupL, sId, 'auth', partId=partNo, altId=1, version=1)
                if (self.__debug):
                    sds.dumpSequence(sId, 'auth', partId=partNo, altId=1, version=1)

                seqFeature.clear()
                seqFeature.setPolymerType(polyTypeCode)
                seqFeature.setPolymerLinkingType(polyType)
                seqFeature.setId(dbName='PDB', dbCode=pdbId, dbAccession=pdbId)

                seqFeature.setSource(organism=pD['SOURCE_NAME'], strain=pD['SOURCE_STRAIN'], taxid=pD['SOURCE_TAXID'], gene=pD['SOURCE_GENE_NAME'],
                                     method=eD['SOURCE_METHOD'], commonName=pD['SOURCE_COMMON_NAME'], variant=pD['SOURCE_VARIANT'])
                seqFeature.setSourceOrig(organism=pD['SOURCE_NAME'], strain=pD['SOURCE_STRAIN'], taxid=pD['SOURCE_TAXID'], gene=pD['SOURCE_GENE_NAME'],
                                         method=eD['SOURCE_METHOD'], commonName=pD['SOURCE_COMMON_NAME'], variant=pD['SOURCE_VARIANT'])

                seqFeature.setAuthPartDetails(partNo, seqNumBeg, seqNumEnd, seqPartType)
                seqFeature.setAuthPartDetailsOrig(partNo, seqNumBeg, seqNumEnd, seqPartType)

                seqFeature.setEntityDescription(description=eD['ENTITY_DESCRIPTION'])
                seqFeature.setEntityDescriptionOrig(description=eD['ENTITY_DESCRIPTION'])

                seqFeature.setEntitySynonyms(synonyms=eD['ENTITY_SYNONYMS'])
                seqFeature.setEntitySynonymsOrig(synonyms=eD['ENTITY_SYNONYMS'])
                #
                seqFeature.setEntityEnzymeClass(ec=eD['ENZYME_CLASS'])
                seqFeature.setEntityEnzymeClassOrig(ec=eD['ENZYME_CLASS'])

                seqFeature.setEntityFragmentDetails(details=eD['FRAGMENT_DETAILS'])
                seqFeature.setEntityFragmentDetailsOrig(details=eD['FRAGMENT_DETAILS'])

                seqFeature.setEntityMutationDetails(details=eD['MUTATION_DETAILS'])
                seqFeature.setEntityMutationDetailsOrig(details=eD['MUTATION_DETAILS'])

                seqFeature.setEntityDetails(details=eD['ENTITY_DETAILS'])
                seqFeature.setEntityDetailsOrig(details=eD['ENTITY_DETAILS'])

                seqFeature.setHostOrgDetails(source=pD['HOST_ORG_SOURCE'], strain=pD['HOST_ORG_STRAIN'], taxid=pD['HOST_ORG_TAXID'],
                                             vector=pD['HOST_ORG_VECTOR'], vectorType=pD['HOST_ORG_VECTOR_TYPE'], plasmid=pD['HOST_ORG_PLASMID'],
                                             commonName=pD['HOST_ORG_COMMON_NAME'], cellLine=pD['HOST_ORG_CELL_LINE'], variant=pD['HOST_ORG_VARIANT'])
                seqFeature.setHostOrgDetailsOrig(source=pD['HOST_ORG_SOURCE'], strain=pD['HOST_ORG_STRAIN'], taxid=pD['HOST_ORG_TAXID'],
                                                 vector=pD['HOST_ORG_VECTOR'], vectorType=pD['HOST_ORG_VECTOR_TYPE'], plasmid=pD['HOST_ORG_PLASMID'],
                                                 commonName=pD['HOST_ORG_COMMON_NAME'], cellLine=pD['HOST_ORG_CELL_LINE'], variant=pD['HOST_ORG_VARIANT'])
                #
                sds.setFeature(seqFeature.get(), sId, 'auth', partId=partNo, altId=1, version=1)

        # ------------------------------------------------------------------------------------------
        # Load reference sequences and features -
        #
        for eId, eD in entityD.items():
            sId = eD['IDENTIFIER']
            polyType = eD['POLYMER_LINKING_TYPE']
            polyTypeCode = sdr.getPolymerTypeCode(polyType)
            #
            rList = eRefD[eId]
            pDList = eD['PART_LIST']
            for (partNo, pD) in enumerate(pDList, start=1):
                seqNumBeg = pD['SEQ_NUM_BEG']
                seqNumEnd = pD['SEQ_NUM_END']
                seqPartId = pD['SEQ_PART_ID']
                seqPartType = pD['SEQ_PART_TYPE']
                #
                altId = 1
                for rD in rList:
                    #
                    # select out the references sequences corresponding to this entity part - mapped to attribute 'fragment_id'
                    #
                    if partNo == int(rD['fragment_id']):
                        myOrderId = int(rD['id'])
                        iBegin = int(rD['hitFrom'])
                        iEnd = int(rD['hitTo'])
                        seqS = rD['subject']
                        # Here the sequence comes from Blast and if iBegin > iEnd is already reverse sense.
                        if iBegin < iEnd:
                            sTup3L = sdr.cnv1To3ListIdx(seqS, iBegin, polyTypeCode)
                        else:
                            sTup3L = sdr.cnv1To3ListIdx(seqS, iBegin, polyTypeCode, indexStep=-1)
                        sds.setSequence(sTup3L, sId, 'ref', partId=partNo, altId=altId, version=1)
                        #
                        seqFeature.clear()
                        # disambiguate the organism and strain data -
                        if rD['db_name'] in ['SP', 'TR', 'UNP']:
                            org, strain = seqFeature.decodeUniProtSourceOrganism(rD['source_scientific'])
                            seqFeature.setSource(organism=org, strain=strain, taxid=rD['taxonomy_id'], commonName=rD['source_common'])
                        else:
                            seqFeature.setSource(organism=rD['source_scientific'], taxid=rD['taxonomy_id'], commonName=rD['source_common'])
                        #
                        seqFeature.setId(dbName=rD['db_name'], dbCode=rD['db_code'], dbAccession=rD['db_accession'], dbIsoform=rD['db_isoform'])
                        seqFeature.setDbIsoformDescription(description=rD['db_isoform_description'])

                        seqFeature.setRefSeqNames(proteinName=rD['name'], synonyms=rD['synonyms'], geneName=rD['gene'])
                        seqFeature.setRefSeqDetails(enzymeClass=rD['ec'], description=rD['db_description'], comments=rD['comments'], keywords=rD['keyword'])
                        seqFeature.setItem('REF_MATCH_BEGIN', iBegin)
                        seqFeature.setItem('REF_MATCH_END', iEnd)
                        seqFeature.setItem('ORG_ORDER_ID', myOrderId)
                        seqFeature.setItem('FULL_LENGTH', int(rD['db_length']))
                        self.__lfh.write("+SequenceDataAssemble.FULL_LENGTH=%d\n" % int(rD['db_length']))
                        seqFeature.setPolymerType(polyTypeCode)
                        seqFeature.setRefSortOrder(sortIndex=rD['sort_order'], sortMetric=rD['sort_metric'])
                        #
                        seqFeature.setItem('AUTH_REF_SEQ_SIM_BLAST', rD['seq_sim'])
                        #
                        seqFeature.setAuthPartDetails(partNo, seqNumBeg, seqNumEnd, seqPartType)
                        sds.setFeature(seqFeature.get(), sId, 'ref', partId=partNo, altId=altId, version=1)
                        altId += 1
                        #

        entryDict = {}
        entryDict['PDB_ID'] = pdbId
        entryDict['RCSB_ID'] = dataSetId
        entryDict['DEPOSITION_DATA_SET_ID'] = dataSetId
        entryDict['PDBX_MODEL_FILE_PATH'] = entryD['PDBX_MODEL_FILE_PATH']
        # JDW
        entryDict['STRUCT_TITLE'] = entryD['STRUCT_TITLE']
        entryDict['CITATION_TITLE'] = entryD['CITATION_TITLE']
        #
        for k, v in groupDict.items():
            sds.setGroup(k, v)
        for k, v in groupPartD.items():
            sds.setGroupParts(k, v)
        for k, v in entryDict.items():
            sds.setEntryDetail(k, v)

        # JDW add default selection to the starting object --
        defSelList = self.makeDefaultSelection(sds=sds)
        ok = sds.setSelectedIds(idList=defSelList)

        sds.serialize()
        self.__sequenceDataStorePath = sds.getFilePath()
        #
        if (self.__verbose):
            sds.dump(self.__lfh)

        if (self.__verbose):
            endTime = time.clock()
            self.__lfh.write("+SequenceDataAssembleWf.__loadSequenceDataStore()  completed at %s (%.2f seconds)\n" %
                             (time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
                              endTime - startTime))
        return True

    def __getAuthSeqComment(self, idx, partTup):
        """      Return initial annotation assignment for residue position 'idx' based on the
                 the input part definition OR an out-of-part assignment of a terminal expression
                 tag or an internal linker.

                 input    idx   sequence position in author/sample sequence (1-N)
                          partTup=(rangeList,seqBegMin,seqEndMax)

                 where    rangeList =[(beg,end,type),(beg,end,type),...]
        """
        # in part ?
        rangeList, seqBegMin, seqEndMax = partTup
        for rl in rangeList:
            if idx < rl[0] or idx > rl[1]:
                continue
            else:
                return rl[2]

        # not in part --
        if idx < seqBegMin or idx > seqEndMax:
            return 'expression tag'
        else:
            return 'linker'

    def __copyFile(self, srcFilePath, dstFilePath):
        """   File copy
        """
        try:
            shutil.copyfile(srcFilePath, dstFilePath)
            return dstFilePath
        except:
            traceback.print_exc(file=self.__lfh)
            return None

    def makeDefaultSelection(self, sds=None):
        """ Make the initial sequence selection and return the corresponding list of selected sequence identifiers
            for each entity/group.
        """
        if (self.__verbose):
            self.__lfh.write("+SequenceDataAssemble.makeDefaultSelection() starting in sessionId %s\n" % (self.__sessionObj.getId()))
        #
        defSelectList = []
        # Get the entity group list -
        #
        gIdList = sds.getGroupIds()

        for gId in gIdList:
            defSelectList.extend(self.makeEntityDefaultSelection(gId, sds=sds))

        if (self.__verbose):
            self.__lfh.write("+SequenceDataAssemble.makeDefaultSelection() returning complete selection list: %r\n" % defSelectList)

        return defSelectList

    def makeEntityDefaultSelection(self, entityId, sds=None):
        """ Make the initial sequence selection and return the corresponding list of selected sequence identifiers
            for the input entity group.
        """
        if (self.__verbose):
            self.__lfh.write("+SequenceDataAssemble.makeEntityDefaultSelection() starting with entityId %s\n" % entityId)
        #
        defSelectList = []
        seqLabel = SequenceLabel()
        seqFeature = SequenceFeature()

        #
        # Get seqIds in group -
        #
        seqIdList = sds.getGroup(entityId)
        if (self.__verbose):
            self.__lfh.write("+SequenceDataAssemble.makeEntityDefaultSelection() group %s contains sequence list %r\n" % (entityId, seqIdList))

        if len(seqIdList) < 1:
            if (self.__verbose):
                self.__lfh.write("+SequenceDataAssemble.makeEntityDefaultSelection() group %s is empty\n" % entityId)
            return defSelectList
        #
        # Coordinate sequences - select the latest version (altId/partId=1) for coordinate sequences.
        #
        altId = 1
        partId = 1
        for seqId in seqIdList:
            verList = sds.getVersionIds(seqId, partId=partId, altId=altId, dataType="sequence", seqType='xyz')
            if len(verList) > 0:
                ver = verList[0]
                seqLabel.set(seqType='xyz', seqInstId=seqId, seqPartId=partId, seqAltId=altId, seqVersion=ver)
                idXyzSeq = seqLabel.pack()
                defSelectList.append(idXyzSeq)
        #
        seqId0 = entityId
        partIdList = sds.getPartIds(seqId0, dataType='sequence', seqType='auth')
        for partId in partIdList:
            #
            # Author sequence - select the latest version
            #
            altId = 1
            verList = sds.getVersionIds(seqId0, altId=altId, partId=partId, dataType="sequence", seqType='auth')
            if len(verList) > 0:
                ver = verList[0]
                seqLabel.set(seqType='auth', seqInstId=seqId0, seqPartId=partId, seqAltId=altId, seqVersion=ver)
                idAuthSeq = seqLabel.pack()
                # JDW only save the selection of the first part of the auth  --
                if int(partId) == 1:
                    defSelectList.append(idAuthSeq)
            #
            # Get selected features of the author sequence -
            #
            seqAuthFD = sds.getFeature(seqId=seqId0, seqType='auth', partId=partId, altId=altId, version=ver)
            seqAuthIdx = sds.getSequence(seqId=seqId0, seqType='auth', partId=partId, altId=altId, version=ver)
            seqFeature.set(seqAuthFD)
            #
            authPolyTypeCode = seqFeature.getPolymerType()
            authTaxId = seqFeature.getSourceTaxId()
            authSeqLen = len(seqAuthIdx)

            # Special cases for default self-reference assignments -- TaxId = 0 --
            #
            #   if (len(authTaxId) < 1  or authTaxId == '0'):
            if (authTaxId == '0'):
                defSelectList.append('selfref_' + str(entityId) + '_' + str(partId))
                if self.__verbose:
                    self.__lfh.write("+SequenceDataAssemble.makeEntityDefaultSelection() Self reference (no taxid) for entity/group %s\n" % entityId)
            #
            # Short nucleic acid polymers --
            #
            # if ((authPolyTypeCode in ['RNA', 'DNA']) and (authSeqLen < 50)) :
            #    defSelectList.append('selfref_'+str(entityId)+'_'+str(partId) )
            #    if self.__verbose:
            #        self.__lfh.write("+SummaryView.__makeEntityDefaultSelection() Self reference short NA polymer for entity/group %s\n" % entityId)
            #
            # Reference sequence   -- minimum sort index --
            #
            sD = {}
            altIdList = sds.getAlternativeIds(seqId0, dataType="sequence", seqType="ref", partId=partId)
            if len(altIdList) == 0:
                defSelectList.append('selfref_' + str(entityId) + '_' + str(partId))
                if self.__verbose:
                    self.__lfh.write("+SequenceDataAssemble.makeEntityDefaultSelection() Self reference entity/group %s no reference sequences\n" % entityId)
                continue
            #
            for altId in altIdList:
                #
                verList = sds.getVersionIds(seqId=seqId0, partId=partId, altId=altId, dataType="sequence", seqType='ref')
                for ver in verList:
                    seqRefFD = sds.getFeature(seqId0, 'ref', partId=partId, altId=altId, version=ver)
                    seqFeature.set(seqRefFD)
                    idx = seqFeature.getRefSortOrderIndex()
                    sVal = seqFeature.getRefSortMetric()
                    sD[idx] = (altId, ver, sVal)
            #
            #
            kys = sD.keys()
            #
            kys.sort(reverse=True)
            altId = sD[kys[0]][0]
            ver = sD[kys[0]][1]
            seqLabel.set(seqType='ref', seqInstId=seqId0, seqPartId=partId, seqAltId=altId, seqVersion=ver)
            idRefSeq = seqLabel.pack()
            defSelectList.append(idRefSeq)

        if (self.__verbose):
            self.__lfh.write("+SequenceDataAssemble.makeEntityDefaultSelection() for entity %r returning selection list: %r\n" % (entityId, defSelectList))

        return defSelectList
