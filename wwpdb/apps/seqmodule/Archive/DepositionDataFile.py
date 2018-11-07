##
# File:  DepositionDataFile.py
# Date:  3-Feb-2010
#
# Updates:
#   20-Apr-2010 jdw Ported to module seqmodule.
##
"""
Utility methods for accessing entity, chain and sequence details from
deposited data files.

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"


import os

from wwpdb.apps.seqmodule.util.SequenceReferenceData import SequenceReferenceData
from wwpdb.api.core.wrapper.CorePyWrap               import ParseCif


class CifFile(object):
    """
    CifFile 
    """

    def __init__(self, fileName):
        self.__fileName = fileName

        self.__cifFile = ParseCif(self.__fileName)


    def getCifFile(self):
        return (self.__cifFile)


    @classmethod
    def getFileExt(cls):
        return ('cif')


    def write(self, fileName):
        self.__cifFile.Write(fileName)


    @classmethod
    def read(cls, fileName):
        return cls(fileName)


class ReferenceSequenceDataFile(object):
    """
    The () class encapsulates access to data in reference sequence files
    
    """

    def __init__(self, dataFile, blockIndex=0):
        self.__cifFile = dataFile.getCifFile()
        #
        self.__blockList = self.__cifFile.GetBlockNames()
        self.__currentBlockIndex=blockIndex
        self.__block = self.__cifFile.GetBlock(self.__blockList[self.__currentBlockIndex])
        #print self.__block.GetTableNames()
            
    def setBlock(self, blockIndex):
        """ Set the current data block to the input zero-based block index in this file object
        """
        self.__currentBlockIndex=blockIndex
        self.__block = self.__cifFile.GetBlock(self.__blockList[self.__currentBlockIndex])        

    def getBlockNames(self):
        """ Get a list of the data blocks in the current file object
        """
        self.__blockList = self.__cifFile.GetBlockNames()
        return self.__blockList

    def getCurrentBlockName(self):
        """ Get the string name for the current data block.
        """
        return (self.__blockList[self.__currentBlockIndex])

    def getReferenceDetails(self):
        """Returns a dictionary of reference sequence details.
        """
        if (not self.__block.IsTablePresent('match_entity')):
            return []
        
        rD={}
        refTable = self.__block.GetTable('match_entity')
        nRows = refTable.GetNumRows()
        #
        colNames = list(refTable.GetColumnNames())
        myList  = ['id', 'db_name', 'db_code', 'db_accession', 'match_length', 'queryFrom', 'queryTo',
                   'hitFrom', 'hitTo', 'identity', 'positive', 'gaps', 'alignLen', 'query', 'subject',
                   'midline', 'query_length', 'name', 'source_scientific', 'source_common', 'taxonomy_id',
                   'gene', 'synonyms', 'comments', 'keyword', 'ec']
        myIntList = ['match_length', 'queryFrom', 'queryTo','hitFrom', 'hitTo', 'identity', 'positive',
                     'gaps', 'alignLen', 'taxonomy_id','query_length']

        rList=[]
        for iRow in range(0,nRows):
            rD={}
            row = refTable.GetRow(iRow)
            for col in myList:
                if col in colNames:
                    if col in myIntList:
                        rD[col] = int(str(row[colNames.index(col)]))
                    elif col in ['query', 'subject', 'midline']:
                        rD[col] = str(row[colNames.index(col)]).upper()                        
                    else:
                        rD[col] = str(row[colNames.index(col)])
            try:
                seq_sim=float(rD['identity']) / float(rD['alignLen'])
            except:
                seq_sim = 0.0
            rD['seq_sim'] = seq_sim
            # placeholder 
            rD['source_string_distance']=100
            rD['taxid_string_distance']=100
            rList.append(rD)
        return rList
        
    def OLDgetReferenceDetails(self):
        """Returns a dictionary of reference sequence details.
        """
        if (not self.__block.IsTablePresent('match_entity')):
            return []
        
        rD={}
        refTable = self.__block.GetTable('match_entity')
        nRows = refTable.GetNumRows()
        #
        colNames = list(refTable.GetColumnNames())
        myList  = ['id', 'db_name', 'db_code', 'db_accession', 'match_length', 'queryFrom', 'queryTo',
                   'hitFrom', 'hitTo', 'identity', 'positive', 'gaps', 'alignLen', 'query', 'subject',
                   'midline', 'query_length', 'name', 'source_scientific', 'source_common', 'taxonomy_id',
                   'gene', 'synonyms', 'comments', 'keyword', 'ec']  
        rList=[]
        for iRow in range(0,nRows):
            rD={}
            row = refTable.GetRow(iRow)
            for col in myList:
                if col in colNames:
                    rD[col] = row[colNames.index(col)]
            try:
                seq_sim=float(rD['identity']) / float(rD['alignLen'])
            except:
                seq_sim = 0.0
            rD['seq_sim'] = seq_sim
            # placeholder 
            rD['source_string_distance']=100
            rD['taxid_string_distance']=100
            rList.append(rD)
        return rList
        

class DepositionDataFile(object):
    """
    The DepositionDataFile() class encapsulates access to data items with deposited entry files.
    
    """

    def __init__(self, dataFile, blockIndex=0):
        self.__cifFile = dataFile.getCifFile()
        #
        self.__blockList = self.__cifFile.GetBlockNames()
        self.__currentBlockIndex=blockIndex
        self.__block = self.__cifFile.GetBlock(self.__blockList[self.__currentBlockIndex])
        self.__polymerEntityChainDict = {}
        self.__chainPolymerEntityDict = {}        
        self.__buildPolymerEntityChainDict()
            
    def setBlock(blockIndex):
        """ Set the current data block to the input zero-based block index in this file object
        """
        self.__currentBlockIndex=blockIndex
        self.__block = self.__cifFile.GetBlock(self.__blockList[self.__currentBlockIndex])        

    def getBlockNames(self):
        """ Get a list of the data blocks in the current file object
        """
        self.__blockList = self.__cifFile.GetBlockNames()
        return self.__blockList

    def getCurrentBlockName(self):
        """ Get the string name for the current data block.
        """
        return (self.__blockList[self.__currentBlockIndex])
        
    def getPolyPeptideLEntityCount(self):
        """  Get the number of entities of type polypeptide(L)
        """
        return (self.getPolymerEntityCount(type="polypeptide(L)"))

    def getPolyPeptideDEntityCount(self):
        """ Get the number of entities of type polypeptide(D)
        """
        return (self.getPolymerEntityCount(type="polypeptide(D)"))    
        

    def getEntityCount(self,type):
        """Returns the integer count of entities of the input 'type'
           (polymer, non-polymer, macrolide, water)
        """
        if type in SequenceReferenceData._entityTypes:
            entityTable = self.__block.GetTable('entity')
            indices = entityTable.Search((type,), ('type',))
            return len(indices)
        else:
            return 0
        
    def getPolymerEntityCount(self,type):
        """Returns the integer count of polymer entities of the input 'type'
             ++ allowed types are reference._polymerEntityTypes - 
        """
        if type in SequenceReferenceData._polymerEntityTypes:
            entityPolyTable = self.__block.GetTable('entity_poly')
            indices = entityPolyTable.Search((type,), ('type',))
            return len(indices)
        else:
            return 0

    def getPolymerEntityType(self,entityId):
        entityPolyTable = self.__block.GetTable('entity_poly')
        indices = entityPolyTable.Search((entityId,), ('entity_id',))
        if len(indices) > 0:
            ty = (entityPolyTable(indices[0], "type"))
            if self.isEmptyValue(ty):
                return ''
            else:
                return ty
        return ''

    def isEmptyValue(self,val):
        if ((val is None) or (len(val) == 0) or (val in ['.','?'])):
            return True
        else:
            return False
        
    def getPolymerEntityList(self,type=None):
        """Returns a list of polymer entity id's  of the input 'type'
           type is an entity type (all, polymer, non-polymer,  any)  or
               one of the polymer entity types.
        """
        if type in ['any','all','polymer']:
            tType='polymer'
            entityTable = self.__block.GetTable('entity')
            indices = entityTable.Search((tType,), ('type',))
            eList=[]
            for ii in indices:
                eList.append(entityTable(ii, "id"))
            return eList
        elif type in SequenceReferenceData._polymerEntityTypes:
            entityPolyTable = self.__block.GetTable('entity_poly')
            indices = entityPolyTable.Search((type,), ('type',))
            eList=[]
            for ii in indices:
                eList.append(entityPolyTable(ii, "entity_id"))
            return eList
        else:
            return []

    def getDbCode(self, dbId):
        """ Return the database code for the input database id/name
        """
        database2 = self.__block.GetTable('database_2')
        indices = database2.Search((dbId,), ('database_id',))
        return (database2(indices[0], "database_code"))

    def getSequence(self, entityId):
        """ Return one-letter-code sequence for the input entity.
        """
        entityPolyTable = self.__block.GetTable('entity_poly')
        indices = entityPolyTable.Search((entityId,), ('entity_id',))
        return (entityPolyTable(indices[0], "ndb_seq_one_letter_code"))

    def getSourceMethod(self, entityId):
        """ Return source method for the input entity - 
        """
        if (not self.__block.IsTablePresent('entity')):
            return  ''
        entityTable = self.__block.GetTable('entity')
        indices = entityTable.Search((entityId,), ('id',))
        if len(indices) > 0:
            return (entityTable(indices[0], "src_method"))
        return ''

    def getSourceGenDetails(self, entityId):
        """ source details -  type = MAN 
        """
        if (not self.__block.IsTablePresent('entity_src_gen')):
            return  ('','','')
        
        entitySrcGenTable = self.__block.GetTable('entity_src_gen')
        indices = entitySrcGenTable.Search((entityId,), ('entity_id',))

        if (len(indices) > 0):
            scientificName=''
            taxId=''
            strain=''
            try: 
                scientificName=entitySrcGenTable(indices[0], "ndb_gene_src_scientific_name")
                strain=entitySrcGenTable(indices[0], "gene_src_strain")
                taxId=entitySrcGenTable(indices[0], "pdbx_gene_src_ncbi_taxonomy_id")
                return (scientificName,strain,taxId)
            except:
                return (scientificName,strain,taxId)                
        else:
            return ('','','')

    def getSourceNatDetails(self, entityId):
        """ source details -  type = NAT 
        """
        if (not self.__block.IsTablePresent('entity_src_nat')):
            return  ('','','')        
        entitySrcNatTable = self.__block.GetTable('entity_src_nat')
        indices = entitySrcNatTable.Search((entityId,), ('entity_id',))
        
        if (len(indices) > 0):
            scientificName=''
            taxId=''
            strain=''            
            try:
                scientificName=entitySrcNatTable(indices[0], "ndb_organism_scientific")
                strain=entitySrcNatTable(indices[0], "strain")
                taxId=entitySrcNatTable(indices[0], "pdbx_ncbi_taxonomy_id")                
                return (scientificName,strain,taxId)
            except:
                return (scientificName,strain,taxId)                
        else:
            return ('','','')

    def getSourceSynDetails(self, entityId):
        """ source details -  type = SYN
        """
        if (not self.__block.IsTablePresent('rcsb_entity_src_syn')):
            return  ('','','')
        
        entitySrcSynTable = self.__block.GetTable('rcsb_entity_src_syn')
        indices = entitySrcSynTable.Search((entityId,), ('entity_id',))
        if (len(indices) > 0):
            scientificName=''
            taxId=''
            strain=''
            try:
                scientificName=entitySrcSynTable(indices[0], "organism_scientific")
                taxId=entitySrcSynTable(indices[0], "ncbi_taxonomy_id")
                strain=entitySrcSynTable(indices[0], "strain")
                return (scientificName,strain,taxId)
            except:
                return (scientificName,strain,taxId)                
        else:
            return ('','','')
        
    def getResidueTableLength(self):
        if (self.__block.IsTablePresent('ndb_poly_seq_scheme')):
            tab = self.__block.GetTable('ndb_poly_seq_scheme')
            return (tab.GetNumRows())
        else:
            return 0
            
    def getEncapsulatedCoordinates(self):
        """  _pdbx_original_pdb_coordinates.coord_section   
        """
        try:
            if (self.__block.IsTablePresent('pdbx_original_pdb_coordinates')):
                encapCoordTable = self.__block.GetTable('pdbx_original_pdb_coordinates')
            elif (self.__block.IsTablePresent('ndb_original_pdb_coordinates')):
                encapCoordTable = self.__block.GetTable('ndb_original_pdb_coordinates')
            elif (self.__block.IsTablePresent('ndb_original_ndb_coordinates')):
                encapCoordTable = self.__block.GetTable('ndb_original_ndb_coordinates')                
            row = encapCoordTable.GetRow(0)                
        except:
            row=['']

        return row[0]

    def getEntitySequence1Auth(self, kd):
        """ Get the dictoinary of author provided one-letter-code sequences corresponding
            to the input keyword arguments:

            kd  -->   entityId=1|2...
                      chainId=A|B...

            example:  getEntitySequenceAuth1(chainId="B")
                      getEntitySequenceAuth1(entityId='1')
        """
        sq=None
        for k,v in kd.items():
            eId=None
            if (k == "entityId"):
                eId = v
            elif (k == "chainId"):
                if self.__chainPolymerEntityDict.has_key(v):
                    eId = self.__chainPolymerEntityDict[v]


            if eId is not None:
                entityPolyTable = self.__block.GetTable('entity_poly')
                indices = entityPolyTable.Search((eId,), ('entity_id',))
                sqR=entityPolyTable(indices[0], "ndb_seq_one_letter_code")
                sq=""
                for cc in sqR:
                    if cc.isspace():
                        continue
                    sq+=cc
                return sq

        return sq

    def getEntityCanSequence1Auth(self,kd):
        """ Get the dictionary of author provided cannonical one-letter-code sequences corresponding
            to the input keyword arguments:
            kd  -->   entityId=1|2...
                      chainId=A|B...

            example:  getEntityCanSequenceAuth1(chainId="B")
                      getEntityCanSequenceAuth1(entityId='1',entityId='2')
        """
        sq=None
        for k,v in kd.items():
            eId=None
            if (k == "entityId"):
                eId = v
            elif (k == "chainId"):
                if self.__chainPolymerEntityDict.has_key(v):
                    eId = self.__chainPolymerEntityDict[v]

            if eId is not None:
                entityPolyTable = self.__block.GetTable('entity_poly')
                indices = entityPolyTable.Search((eId,), ('entity_id',))
                sqR=entityPolyTable(indices[0], "ndb_seq_one_letter_code_can")
                sq=""
                for cc in sqR:
                    if cc.isspace():
                        continue
                    sq+=cc
                return sq
        return sq


    def getPdbChainIdList(self, entityId):
        entityPolyTable = self.__block.GetTable('entity_poly')
        indices = entityPolyTable.Search((entityId,), ('entity_id',))

        st=entityPolyTable(indices[0], "ndb_chain_id")
        tList=[]
        if st is not None and len(st) > 0:
            if ((len(st) > 1) and (st.count(',') > 0)):
                tList=st.split(',')
            elif ((len(st) > 1) and (st.count(' ') > 0)):
                tList=st.split()                
            else:
                tList=st.split(',')                
        rList=[]
        for ch in tList:
            if len(ch) == 0 or ch in ['.','?']:
                continue
            rList.append(str(ch).strip())
        return rList

    def getPolymerEntityChainDict(self):
        if len(self.__polymerEntityChainDict) == 0:
            self.__buildPolymerEntityChainDict()
        return self.__polymerEntityChainDict

    def getChainPolymerEntityDict(self):
        if len(self.__chainPolymerEntityDict) == 0:
            self.__buildPolymerEntityChainDict()
        return self.__chainPolymerEntityDict
    
    def __buildPolymerEntityChainDict(self):
        """
        """
        self.__polymerEntityChainDict = {}
        pEntityList=self.getPolymerEntityList('all')
        for eId in pEntityList:
            self.__polymerEntityChainDict[eId]= self.getPdbChainIdList(eId)
        self.__chainPolymerEntityDict = {}
        for eId,cList in  self.__polymerEntityChainDict.items():
            for cId in cList:
                self.__chainPolymerEntityDict[cId]=eId
    

    def getSequence3AlignList(self,chainId):
        """
        Return a list of tuples containing (PDB Monomer id (xyz), PDBx mononer(entity), Ent. Poly Seq Num,  PDB residue Num)
        for the aligned entity and coordinate sequences (expressed as 3-letter-codes).
        
        """
        tabO = self.__block.GetTable('ndb_poly_seq_scheme')
        indices = tabO.Search((chainId,), ('pdb_id',))
        aa3List=[]
        for ii in indices:
            aa3List.append((tabO(ii, "pdb_mon_id"), tabO(ii, "mon_id"), tabO(ii, "seq_id"), tabO(ii, "pdb_num")) )
            
        return aa3List

    def residueMapTableExists(self):
        if (self.__block.IsTablePresent('ndb_poly_seq_scheme')):
            return True
        else:
            return False
        
    def getCoordinateSequenceList(self,chainId):
        """
        Return a list of tuples containing (PDB Monomer 3-letter-code, PDB residue num, '', Ent. Poly Seq Num, )
        for the aligned entity and coordinate sequences (expressed as 3-letter-codes).
        
        """
        tabO = self.__block.GetTable('ndb_poly_seq_scheme')
        indices = tabO.Search((chainId,), ('pdb_id',))
        aa3List=[]
        for ii in indices:
            monId=tabO(ii, "pdb_mon_id")
            if ((len(monId) > 0) and (monId != "?") and (monId != ".")):
                aa3List.append((tabO(ii, "pdb_mon_id"), tabO(ii, "pdb_num"), '', tabO(ii, "seq_id") ) )
            
        return aa3List


    def getEntityPolyList(self):
        """Returns a list of polymer entity ids
        """
        rD={}
        polyTable = self.__block.GetTable('entity_poly')
        nRows = polyTable.GetNumRows()
        #
        colNames = list(polyTable.GetColumnNames())

        eList=[]
        for iRow in range(0,nRows):
            row = polyTable.GetRow(iRow)
            eId = row[colNames.index('entity_id')]
            eList.append(eId)
        return eList

    def getSequence1Xyz(self,chainId):
        """
        Get the one-letter-code sequence from the deposited coordinate data
        """
        tabO = self.__block.GetTable('ndb_poly_seq_scheme')
        indices = tabO.Search((chainId,), ('pdb_id',))
        mon3List=[]
        for ii in indices:
            mon3List.append(tabO(ii, "pdb_mon_id") )

        # translate to one-letter-codes --
        sq=""
        for mon in mon3List:
            if SequenceReferenceData._monDict3.has_key(mon):
                sq+=SequenceReferenceData._monDict3[mon]
            elif ((mon in ['.','?']) or (len(mon) == 0)):
                sq+="-"
            else:
                sq+="X"
        
        return sq

    def getMisalignCount(self,chainId):
        """
        Return the count of misaligned residues in maxit sequence alignment.
        
        """
        tabO = self.__block.GetTable('ndb_poly_seq_scheme')
        indices = tabO.Search((chainId,), ('pdb_id',))
        mCount = 0
        for ii in indices:
            pdb_mon_id = tabO(ii, "pdb_mon_id")
            mon_id = tabO(ii, "mon_id")
            if pdb_mon_id != mon_id:
                mCount +=1
            
        return mCount

    def getAtomSiteTableLength(self):
        if (self.__block.IsTablePresent('atom_site')):
            tab = self.__block.GetTable('atom_site')
            return (tab.GetNumRows())
        else:
            return 0

    def getSequenceFromAtomSite(self):
        """ Extract the polymer sequence from the PDBx atom_site category for the first model.

            self.__chD{} =[(a3,orglabel index, comment/details,  align index placeholder),(),..] is a dictionary with chain id key
                          containing sequences stored as a list of tuples.

            Each tuple has the value (3-letter-code, original label residue index, comment/details, alignment index )            
        """
        
        if (not self.__block.IsTablePresent('atom_site')):
            return {}
        
        atTable = self.__block.GetTable('atom_site')
        nRows = atTable.GetNumRows()
        #
        colNames = list(atTable.GetColumnNames())
        
        i1=colNames.index('auth_asym_id')
        i2=colNames.index('auth_comp_id')
        i3=colNames.index('auth_seq_id')
        #
        if 'pdbx_PDB_model_num' in colNames:
            i4=colNames.index('pdbx_PDB_model_num')
        elif 'ndb_model' in colNames:
            i4=colNames.index('ndb_model')
        #
        chD={}        
        rList=[]
        idx=1
        row = atTable.GetRow(0)
        rList.append((row[i2],row[i3],'',idx))
        idx+=1
        chId_cur=row[i1]
        compId_cur=row[i2]
        seqId_cur=row[i3]
        modelId_cur=row[i4]                
        #
        #
        for iRow in range(1,nRows):
            row = atTable.GetRow(iRow)
            chId=row[i1]
            compId=row[i2]
            seqId=row[i3]
            modelId=row[i4]
            #
            if modelId != modelId_cur:
                break
            #
            if chId != chId_cur:
                if not chD.has_key(chId_cur):
                    chD[chId_cur]=rList
                rList=[]
                idx=1
                chId_cur=chId

            
            if ((compId == compId_cur) and (seqId == seqId_cur)):
                continue
            else:
                # jdw - filter waters - 
                if (compId != 'HOH'):
                    rList.append((compId,seqId,'',idx))
                    idx+=1
                    compId_cur=compId
                    seqId_cur=seqId

                

        # assign the last chain
        if not chD.has_key(chId_cur):
            chD[chId_cur]=rList        

        #
        #for k,vL in chD.items():
        #    print "Chain ",k, len(vL)
        #   for v in vL:
        #        print v
        #
        return chD
        


if __name__ == "__main__":
    pass
