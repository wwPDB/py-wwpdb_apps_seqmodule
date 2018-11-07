##
# File:  SequenceReferenceData.py
# Date:  22-Dec-2009
#
# Updates:
# 20-Apr-2010 jdw Port to module seqmodule.
# 19-Mar-2014 jdw Add 3-letter to 1-letter conversions
# 21-Mar-2014 jdw Add 3-letter to 1-letter conversions with formatting
# 20-May-2014 jdw Add support for PYL and SEC 
#  3-Jun-2014 jdw add convertDbNameToResource(dbName)
#  7-Jul-2014 jdw add indexStep=+/-1 to cnv1To3ListIdx() and related methods -
#  7-Sep-2017 zf  add unknown nucleotide 'N'
#
##
"""
Sequence reference data and utility methods on reference data.

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"

import sys, string


class SequenceReferenceData(object):
    """ This class contains reference data for identifying standard polymer types and
        mapping residue nomenclature.
    """
    _polymerEntityTypes = ['polypeptide(D)','polypeptide(L)','polydeoxyribonucleotide','polyribonucleotide',
                                 'polysaccharide(D)','polysaccharide(L)','polydeoxyribonucleotide/polyribonucleotide hybrid',
                                 'cyclic-pseudo-peptide','other']
    
    _entityTypes = ['polymer','non-polymer','macrolide','water']
    
    _monDict3 = {
        "ALA":"A",
        "ARG":"R",
        "ASN":"N",
        "ASP":"D",
        "ASX":"B",
        "CYS":"C",
        "GLN":"Q",
        "GLU":"E",
        "GLX":"Z",
        "GLY":"G",
        "HIS":"H",
        "ILE":"I",
        "LEU":"L",
        "LYS":"K",
        "MET":"M",
        "PHE":"F",
        "PRO":"P",
        "SER":"S",
        "THR":"T",
        "TRP":"W",
        "TYR":"Y",
        "VAL":"V",
        "PYL":"O",
        "SEC":"U",
        "DA":"A",
        "DC":"C",
        "DG":"G",
        "DT":"T",
        "DU":"U",
        "DI":"I",
        "A":"A",
        "C":"C",
        "G":"G",
        "I":"I",
        "N":"N",
        "T":"T",
        "U":"U",
        "UNK":"X",
        #         "MSE":"M",
        ".":"."
        }            

    _monDictDNA1 = {
        ".":".",                
        "A":"DA",
        "C":"DC",
        "G":"DG",
        "I":"DI",
        "N":"N",
        "T":"DT",
        "U":"DU",
        "X":"UNK"
        }
    
    _monDictRNA1 = {
        ".":".",        
        "A":"A",
        "C":"C",
        "G":"G",
        "I":"I",
        "N":"N",
        "T":"T",
        "U":"U",
        "X":"UNK"
        }
        
    _monDict1 = {
        "." :   ".",
        "A" :   "ALA",
        "R" :   "ARG",
        "N" :   "ASN",
        "D" :   "ASP",
        "B" :   "ASX",
        "C" :   "CYS",
        "Q" :   "GLN",
        "E" :   "GLU",
        "Z" :   "GLX",
        "G" :   "GLY",
        "H" :   "HIS",
        "I" :   "ILE",
        "L" :   "LEU",
        "K" :   "LYS",
        "M" :   "MET",
        "F" :   "PHE",
        "P" :   "PRO",
        "S" :   "SER",
        "T" :   "THR",
        "W" :   "TRP",
        "Y" :   "TYR",
        "V" :   "VAL",
        "U" :   "SEC",
        "O" :   "PYL"   }
    
    _monDNAList = ["DA","DG","DT","DC","DI","DU"]

    _complimentRNA={'A':'U',
                    'T':'A',
                    'G':'C',
                    'C':'G',
                    'U':'A'}
    
    _complimentDNA={'A':'T',
                    'T':'A',
                    'G':'C',
                    'C':'G'}

    def __init__(self,verbose=False,log=sys.stderr):
        self.__verbose=verbose
        self.__lfh=log
        self.__gapSymbol='.'


    def convertDbNameToResource(self,dbName):
        if dbName in ['SP','TR']:
            return 'UNP'
        elif dbName in ['GB','DBJ','EMB','EMBL']:
            return 'GB'
        else:
            return dbName

    def isDNA(self,r3):
        return (r3 in SequenceReferenceData._monDNAList)
        
    def getGapSymbol(self):
        return self.__gapSymbol
    
    def getPolymerTypeCode(self,polymerType):
        # if (self.__verbose):
        #    sys.stderr.write("+SequenceDataImportRcsb.__getPolymerTypeCode()  polymerType |%s|\n" % (polymerType))                
        if (polymerType is None or len(polymerType) < 2):
            return 'AA'
        elif (polymerType.find('peptide') != -1):
            return 'AA'
        elif ( polymerType.find('hybrid') != -1):
            return 'XNA'
        elif ( polymerType.find('polyribonuc') != -1):
            return 'RNA'
        elif  ( polymerType.find('polydeoxyribonuc') != -1):
            return 'DNA'
        elif ( polymerType.find('polysac') != -1):
            return 'SAC'
        else:
            return 'AA'

    def parseSequence(self,ss,typeCode):
        if typeCode == 'AA':
            return self.__parseSequenceAA(ss)
        elif typeCode == "RNA" or typeCode == "XNA" :
            return self.__parseSequenceRNA(ss)            
        elif typeCode == "DNA" :
            return self.__parseSequenceDNA(ss)
        else:
            return self.__parseSequenceAA(ss)            

    def guessPolymerType(self,ss):
        """ Parse a string like AAA(MOD)(MOD)AAA and convert
            this to a list of 1- and 3-letter code residue lists.
        """
        if ss is None or len(ss) < 1:
            return ''

        r1L=[]
        inP=False
        for s in ss:
            if s in string.whitespace:
                continue
            if s == "(":
                inP=True
                continue
            if s == ")":
                inP=False

            if not inP:
                r1L.append(s)

        aaL=SequenceReferenceData._monDict1.keys()
        naL=SequenceReferenceData._monDictDNA1.keys()

        #if (self.__verbose):
        #    self.__lfh.write("Sequence list : %r\n" % r1L)            
        #    self.__lfh.write("Sequence list AA : %r\n" % aaL)
        #    self.__lfh.write("Sequence list NA : %r\n" % naL)
            
        notAA=0
        notNA=0
        numT=0
        for r1 in r1L:
            if r1 not in aaL:
                notAA += 1
            if r1 not in naL:
                notNA += 1
            if r1 == 'T':
                numT += 1

        if notAA < notNA:
            return 'polypeptide'
        elif numT > 0:
            return 'polydeoxyribonucleotide'
        else:
            return 'polyribonucleotide'
        
    def __parseSequenceAA(self,ss):
        """ Parse a string like AAA(MOD)(MOD)AAA and convert
            this to a list of 1- and 3-letter code residue lists.
        """
        r1L=[]
        r3L=[]
        if ss is None or len(ss) < 1:
            return r1L,r3L

        inP=False
        for s in ss:
            if s in string.whitespace:
                continue
            if s == "(":
                inP=True
                r3=""
                continue
            if s == ")":
                inP=False
                r3L.append(r3)
                if SequenceReferenceData._monDict3.has_key(r3):
                    r1L.append(SequenceReferenceData._monDict3[r3])
                else:
                    r1L.append('X')
                continue

            if inP:
                r3 += s
            else:
                if SequenceReferenceData._monDict1.has_key(s):
                    r3L.append(SequenceReferenceData._monDict1[s])
                else:
                    r3L.append('UNK')
                r1L.append(s)

        return (r1L,r3L)

    def __parseSequenceRNA(self,ss):
        """ Parse a string like AAA(MOD)(MOD)AAA and convert
            this to a list of 1- and 3-letter code residue lists.
        """

        r1L=[]
        r3L=[]
        if ss is None or len(ss) < 1:
            return r1L,r3L

        inP=False
        for s in ss:
            if s in string.whitespace:
                continue
            if s == "(":
                inP=True
                r3=""
                continue
            if s == ")":
                inP=False
                r3L.append(r3)
                if SequenceReferenceData._monDict3.has_key(r3):
                    r1L.append(SequenceReferenceData._monDict3[r3])
                else:
                    r1L.append('X')
                continue

            if inP:
                r3 += s
            else:
                if SequenceReferenceData._monDictRNA1.has_key(s):
                    r3L.append(SequenceReferenceData._monDictRNA1[s])
                else:
                    r3L.append('UNK')
                r1L.append(s)

        return (r1L,r3L)

    def __parseSequenceDNA(self,ss):
        """ Parse a string like AAA(MOD)(MOD)AAA and convert
            this to a list of 1- and 3-letter code residue lists.
        """


        r1L=[]
        r3L=[]
        if ss is None or len(ss) < 1:
            return r1L,r3L

        inP=False
        for s in ss:
            if s in string.whitespace:
                continue
            if s == "(":
                inP=True
                r3=""
                continue
            if s == ")":
                inP=False
                r3L.append(r3)
                if SequenceReferenceData._monDict3.has_key(r3):
                    r1L.append(SequenceReferenceData._monDict3[r3])
                else:
                    r1L.append('X')
                continue

            if inP:
                r3 += s
            else:
                if SequenceReferenceData._monDictDNA1.has_key(s):
                    r3L.append(SequenceReferenceData._monDictDNA1[s])
                else:
                    r3L.append('UNK')
                r1L.append(s)

        return (r1L,r3L)

    def compliment1NA(self,s1S,polyTypeCode):
        s1L=self.__toList(s1S)
        o1L=[]
        if polyTypeCode == "RNA" or polyTypeCode == "XNA":
            for r1 in s1L:
                if self._complimentRNA.has_key(r1):
                    o1L.append(self._complimentRNA[r1])
                else:
                    if (self.__verbose):
                        self.__lfh.write("+SequenceReferenceData.compliment1NA() unexpected %s nucleotide code %r\n" % (polyTypeCode,r1))
                    o1L.append(r1)
        elif polyTypeCode == "DNA":
            for r1 in s1L:
                if self._complimentDNA.has_key(r1):
                    o1L.append(self._complimentDNA[r1])
                else:
                    if (self.__verbose):
                        self.__lfh.write("+SequenceReferenceData.compliment1NA() unexpected %s nucleotide code %r\n" % (polyTypeCode,r1))
                    o1L.append(r1)

        return ''.join(o1L)

    def cnv1To3ListIdx(self,s1S,iBegin,polyTypeCode,indexStep=1):
        if polyTypeCode == "AA":
            return self.__cnv1To3ListIdxAA(s1S,iBegin,indexStep=indexStep)
        elif polyTypeCode == "RNA" or polyTypeCode == "XNA":
            return self.__cnv1To3ListIdxRNA(s1S,iBegin,indexStep=indexStep)
        elif polyTypeCode == "DNA":
            return self.__cnv1To3ListIdxDNA(s1S,iBegin,indexStep=indexStep)            
        else:
            return self.__cnv1To3ListIdxAA(s1S,iBegin,indexStep=indexStep)

    def __cnv1To3ListIdxAA(self,s1S,iBegin,indexStep=1):
        s1L=self.__toList(s1S)
        sTup3L=[]
        ir=int(iBegin)
        idx=1
        for r1 in s1L:
            if r1 in ['-']:
                continue
            if SequenceReferenceData._monDict1.has_key(r1):
                sTup3L.append((SequenceReferenceData._monDict1[r1],str(ir),'',idx))
            else:
                sTup3L.append(('UNK',str(ir),'',idx))
                if (self.__verbose):
                    self.__lfh.write("+SequenceReferenceData.__cnv1To3ListIdAA() Failed to map one-letter-code %s using UNK\n" % r1)
            ir+=indexStep
            idx+=1
        return sTup3L

    def __cnv1To3ListIdxRNA(self,s1S,iBegin,indexStep=1):
        s1L=self.__toList(s1S)
        sTup3L=[]
        ir=int(iBegin)
        idx=1
        for r1 in s1L:
            if SequenceReferenceData._monDictRNA1.has_key(r1):
                sTup3L.append((SequenceReferenceData._monDictRNA1[r1],str(ir),'',idx))
            else:
                sTup3L.append(('UNK',str(ir),'',idx))
                if (self.__verbose):
                    self.__lfh.write("+SequenceReferenceData.__cnv1To3ListIdRNA() Failed to map one-letter-code %s using UNK\n" % r1)
            ir+=indexStep
            idx+=1            
        return sTup3L

    def __cnv1To3ListIdxDNA(self,s1S,iBegin,indexStep=1):
        s1L=self.__toList(s1S)
        sTup3L=[]
        ir=int(iBegin)
        idx=1
        for r1 in s1L:
            if SequenceReferenceData._monDictRNA1.has_key(r1):
                sTup3L.append((SequenceReferenceData._monDictDNA1[r1],str(ir),'',idx))
            else:
                sTup3L.append(('UNK',str(ir),'',idx))
                if (self.__verbose):
                    self.__lfh.write("+SequenceReferenceData.__cnv1To3ListIdRNA() Failed to map one-letter-code %s using UNK\n" % r1)
            ir+=indexStep
            idx+=1
        return sTup3L


    def cnv1To3List(self,s1S,polyTypeCode):
        if polyTypeCode == "AA":
            return self.__cnv1To3ListAA(s1S)
        elif polyTypeCode == "RNA" or  polyTypeCode == "XNA":
            return self.__cnv1To3ListRNA(s1S)
        elif polyTypeCode == "DNA":
            return self.__cnv1To3ListDNA(s1S)            
        else:
            return self.__cnv1To3ListAA(s1S)


    def __cnv1To3ListAA(self,s1S):
        s1L=self.__toList(s1S)
        s3L=[]
        for r1 in s1L:
            if SequenceReferenceData._monDict1.has_key(r1):
                s3L.append(SequenceReferenceData._monDict1[r1])
            else:
                s3L.append('UNK')
        return s3L

    def __cnv1To3ListDNA(self,s1S):
        s1L=self.__toList(s1S)
        s3L=[]
        for r1 in s1L:
            if SequenceReferenceData._monDict1.has_key(r1):
                s3L.append(SequenceReferenceData._monDictDNA1[r1])
            else:
                s3L.append('UNK')
        return s3L

    def __cnv1To3ListRNA(self,s1S):
        s1L=self.__toList(s1S)
        s3L=[]
        for r1 in s1L:
            if SequenceReferenceData._monDict1.has_key(r1):
                s3L.append(SequenceReferenceData._monDictRNA1[r1])
            else:
                s3L.append('UNK')
        return s3L

    def __toList(self,strIn):
        sL=[]
        for ss in strIn:
            if ss in string.whitespace:
                continue
            sL.append(ss)
        return sL

    def cnvList3to1(self,r3List):
        oL=[]
        for r3 in r3List:
            if r3 == self.__gapSymbol:
                oL.append( self.__gapSymbol )
            elif SequenceReferenceData._monDict3.has_key(r3):
                oL.append( SequenceReferenceData._monDict3[r3] )
            else:
                oL.append('X')
        return ''.join(oL)

    def cnvList3to1WithMods(self,r3List):
        oL=[]
        for r3 in r3List:
            if r3 == self.__gapSymbol:
                oL.append( self.__gapSymbol )
            elif SequenceReferenceData._monDict3.has_key(r3):
                oL.append( SequenceReferenceData._monDict3[r3] )
            else:
                oL.append('('+r3+')')
        return ''.join(oL)

    def cnvList3to1WithModsFormatted(self,r3List,maxLine=60):
        """  Return one-letter-code including parenthetical residue modifications formatted with 
             maximum line length -
        """
        oL=[]
        myLen=0
        for r3 in r3List:
            if r3 == self.__gapSymbol:
                tS=self.__gapSymbol
            elif SequenceReferenceData._monDict3.has_key(r3):
                tS=SequenceReferenceData._monDict3[r3]
            else:
                tS='('+r3+')'

            myLen+=len(tS)
            if myLen > maxLine:
                oL.append('\n')
                myLen=1
            oL.append(tS)
            
        return ''.join(oL)

        
    def cnv3To1(self,r3):
        if r3 == self.__gapSymbol:
            return self.__gapSymbol
        elif SequenceReferenceData._monDict3.has_key(r3):
            return SequenceReferenceData._monDict3[r3]
        else:
            return 'X'        

    def isStandard3(self,r3):
        return SequenceReferenceData._monDict3.has_key(r3)


