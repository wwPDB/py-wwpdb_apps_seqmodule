##
# File:    BlastPlusReader.py
# Date:    17-Apr-2013  JDW
#
# Adapted from  seqdb_v2/ReadNcbiBlastXml.py
#
#
# Updates:
#
#  18-Apr-2013  jdw  modify to handle current conventions in UNP fasta files comment lines
#                    Extract isoform information from comment lines -- 
#  19-Apr-2013  jdw  remove any identity-based filtering of output
##
 


from xml.dom import minidom
import sys, math, getopt,  traceback
#from wwpdb.utils.seqdb_v2.FetchNcbiXml import FetchNcbiXml

class BlastPlusReader:
    """Read Blast+ result file (xml format)  and return the list of dictionaries containing the following -- 

           dict['db_name']
           dict['db_code']
           dict['db_accession']
           dict['db_isoform']
           dict['identity']
           dict['positive']
           dict['gaps']
           dict['midline']
           dict['query']
           dict['queryFrom']
           dict['queryTo']
           dict['subject']
           dict['hitFrom']
           dict['hitTo']
           dict['alignLen']
           dict['match_length']

    """

    def __init__(self,verbose=True,log=sys.stderr):
        self.__verbose=verbose
        self.__debug=True
        self.__lfh=log
        self._resultList=[]
        self.__sequenceType='polypeptide(L)'

    def setSequenceType(self,polyType):
        if polyType in ['polypeptide(L)', 'polypeptide(D)','polypeptide','polyribonucleotide']:
            self.__sequenceType=polyType
            return True
        else:
            return False

    def readFile(self,filePath):
        try:
            if (self.__verbose):
                self.__lfh.write("+BlastPlusReader.readFile() reading %s\n" % filePath)
            ifh=open(filePath,'r')
            inputText=ifh.read()
            ifh.close()
            if len(inputText) > 0:
                domObj = minidom.parseString(inputText)
                self._resultList = self._parse(domObj)
        except:
            if (self.__verbose):
                self.__lfh.write("+BlastPlusReader.readFile() failed for %s\n" % filePath)
                traceback.print_exc(file=self.__lfh)


        return self._resultList

    def GetResultList(self):
        return self._resultList

    def _parse(self, domObj):
        resultList = []
        list = domObj.getElementsByTagName('Hit')
        if not list:
            return resultList

        length = None
        for node in domObj.getElementsByTagName('BlastOutput_query-len'):
            if node.nodeType != node.ELEMENT_NODE:
                continue

            if node.tagName == 'BlastOutput_query-len':
                length = node.firstChild.data
                break


        for node in list:
            if node.nodeType != node.ELEMENT_NODE:
                continue

            alignlist = self._ProcessHitTag(node.childNodes, length)
            if alignlist:
                for align in alignlist:
                    if length:
                        align['query_length'] = length
                    resultList.append(align)

        return resultList

    def _ProcessHitTag(self, nodelist, length):
        resultList = []
        alignlist = []
        description =''
        length=''
        hitId = ''
        code = ''
        for node in nodelist:
            if node.nodeType != node.ELEMENT_NODE:
                continue

            if node.tagName == 'Hit_id':
                dbName,dbCode,accCode,isoForm=self._parseID(node.firstChild.data)
                if not accCode:
                    return resultList
            elif node.tagName == 'Hit_accession':
                code = node.firstChild.data
            elif node.tagName == 'Hit_def':
                description = node.firstChild.data
            elif node.tagName == 'Hit_len':
                length = node.firstChild.data
            elif node.tagName == 'Hit_hsps':
                list = self._ProcessHit_hspsTag(node.childNodes, length)
                if list:
                    for li in list:
                        alignlist.append(li)
    
        if not accCode or not alignlist:
            return resultList

        #dict={}

        for align in alignlist:
            align['db_name'] = dbName
            align['db_code'] = dbCode
            align['db_accession'] = accCode
            align['db_isoform'] = isoForm
            align['db_description'] = description
            align['db_length'] = length
            resultList.append(align)

        return resultList

    def _parseID(self, data):
        accCode = ''
        dbCode  = ''
        dbName  = ''
        isoForm = ''
        list = data.split('|')
        if len(list) >= 2 and list[2] != 'pdb':
            f0 = str(list[0]).upper()
            if f0 in ['TR','SP']:
                dbName = f0
                accCode = str(list[1])
                dbCode  = str(list[2])
                isoForm=''
                if (accCode is not None and (accCode.find('-') != -1)):
                    tL=accCode.split('-')
                    if len(tL)>1 and len(tL[1])>0:
                        #  Isoform is the acession + '-' + variant #
                        #isoForm=tL[1]
                        isoForm=accCode
                        accCode=tL[0]
            elif f0 in ['GI']:
                dbName = str(list[2]).upper()
                accCode = str(list[1])
                dbCode  = str(list[3])

        return dbName,dbCode,accCode,isoForm

    def _ProcessHit_hspsTag(self, nodelist, length):
        resultList = []
        for node in nodelist:
            if node.nodeType != node.ELEMENT_NODE:
                continue
            if node.tagName != 'Hsp':
                continue

            dict = self._GetMatchAlignment(node.childNodes, length)
            if dict:
                resultList.append(dict)
        return resultList

    def _GetMatchAlignment(self, nodelist, length):
        dict = {}
        hsp_num=''
        for node in nodelist:
            if node.nodeType != node.ELEMENT_NODE:
                continue
    
            if node.tagName == 'Hsp_identity':
                dict['identity'] = node.firstChild.data
            elif node.tagName == 'Hsp_positive':
                dict['positive'] = node.firstChild.data
            elif node.tagName == 'Hsp_gaps':
                dict['gaps'] = node.firstChild.data 
            elif node.tagName == 'Hsp_midline':
                if self.__sequenceType == 'polyribonucleotide':
                    dict['midline'] = node.firstChild.data.replace('T', 'U')
                else:
                    dict['midline'] = node.firstChild.data
            elif node.tagName == 'Hsp_qseq':
                if self.__sequenceType == 'polyribonucleotide':
                    dict['query'] = node.firstChild.data.replace('T', 'U')
                else:
                    dict['query'] = node.firstChild.data
            elif node.tagName == 'Hsp_query-from':
                dict['queryFrom'] = node.firstChild.data
            elif node.tagName == 'Hsp_query-to':
                dict['queryTo'] = node.firstChild.data
            elif node.tagName == 'Hsp_hseq':
                if self.__sequenceType == 'polyribonucleotide':
                    dict['subject'] = node.firstChild.data.replace('T', 'U')
                else:
                    dict['subject'] = node.firstChild.data
            elif node.tagName == 'Hsp_hit-from':
                dict['hitFrom'] = node.firstChild.data
            elif node.tagName == 'Hsp_hit-to':
                dict['hitTo'] = node.firstChild.data
            elif node.tagName == 'Hsp_align-len':
                dict['alignLen'] = node.firstChild.data
                dict['match_length'] = node.firstChild.data
            elif node.tagName == 'Hsp_num':
                hsp_num = node.firstChild.data
                  

        # only take the first alignment within the hit 
        if str(hsp_num) != "1":
            dict.clear()

        if (False):
            if dict.has_key('identity'):
                if length:
                    identity = int(dict['identity']) * 100 / int(length)
                    if identity < 60:
                        dict.clear()
                elif dict.has_key('alignLen'):
                    identity = int(dict['identity']) * 100 / int(dict['alignLen'])
                    if identity < 60:
                        dict.clear()

        
        return dict

