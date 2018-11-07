##
# File:    WebRequest.py
# Date:    18-Jan-2010
#
# Updated:
# 20-Apr-2010 Ported to seqmodule package
# 27-Jul-2010  RPS: Added support for accommodating different ordering of sequence types as per user preferences
##
"""
WebRequest provides containers and accessors for managing request parameter information.

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"


import sys
from json import loads, dumps

from wwpdb.apps.seqmodule.io.SessionManager import SessionManager

class WebRequest(object):
    """ Base container and accessors for input and output parameters and control information. 
    """
    def __init__(self,paramDict={},verbose=False):
        self.__verbose=verbose
        #
        #  Input and storage model is dictionary of lists (e.g. dict[myKey] = [,,,])
        #  Single values are stored in the leading element of the list (e.g. dict[myKey][0])
        #
        self.__dict=paramDict
        

    def printIt(self,ofh=sys.stdout):
        try:
            ofh.write("\n--------------------------------------------:\n")                         
            ofh.write("\nWebRequest.printIt() Input/Output Request Dictionary Contents:\n")             
            for k,vL in self.__dict.items():
                ofh.write("  Key: %s  value(s): " % k)
                for v in vL:
                    ofh.write(" %s " % v)
                ofh.write("\n")
            ofh.write("\n--------------------------------------------\n\n")                                         
        except:
            pass


    def dump(self,format='text'):
        oL=[]
        try:
            if (format == 'html'):
                oL.append('<pre>\n')
            oL.append("\n--------------------------------------------:\n")                         
            oL.append("\nWebRequest.printIt() Input/Output Request Dictinary Contents:\n")             
            for k,vL in self.__dict.items():
                oL.append("  Key: %s  value(s): " % k)
                for v in vL:
                    oL.append(" %r " % v)
                oL.append("\n")
            oL.append("\n--------------------------------------------\n\n")
            if (format == 'html'):
                oL.append('</pre>\n')
        except:
            pass
        
        return oL        

    def getJSON(self):
        return dumps(self.__dict)

    def setJSON(self,JSONString):
        self.__dict = loads(JSONString)
                    
    def getValue(self,myKey):
        return(self._getStringValue(myKey))

    def getValueList(self,myKey):
        return(self._getStringList(myKey))

    def getRawValue(self,myKey):
        return(self._getRawValue(myKey))
    
    #
    def setValue(self,myKey,aValue):
        self.__dict[myKey]=[aValue]

    def setValueList(self,myKey,valueList):
        self.__dict[myKey]=valueList        

    
    def exists(self,myKey):
        try:
            return self.__dict.has_key(myKey)
        except:
            return False

    #
    def _getRawValue(self,myKey):
        try:
            return self.__dict[myKey][0]
        except:
            return None
        
    def _getStringValue(self,myKey):
        try:
            return str(self.__dict[myKey][0]).strip()
        except:
            return ''


    def _getIntegerValue(self,myKey):
        try:
            return int(self.__dict[myKey][0])
        except:
            return None

    def _getDoubleValue(self,myKey):
        try:
            return double(self.__dict[myKey][0])
        except:
            return None

    def _getStringList(self,myKey):
        try:
            return self.__dict[myKey]
        except:
            return []
            
class SequenceInputRequest(WebRequest):
    def __init__(self,paramDict,verbose=False,log=sys.stderr):
        super(SequenceInputRequest,self).__init__(paramDict,verbose)
        self.__verbose = verbose
        self.__lfh=log

    def getDeleteList(self):
        """  Return a list of tuples containing [(elementId,previousValue),(),...]
        """
        oL=[]
        try:
            dString=self._getStringValue('deleteselect')
            #self.__lfh.write("\n+SequenceInputRequest.getDeleteList() dString length %d : %r\n" % (len(dString),dString))
            dList=dString.split(',')
            if (self.__verbose):
                self.__lfh.write("\n+SequenceInputRequest.getDeleteList() dList length %d : %r\n" % (len(dList),dList))            
            for d in dList:
                #self.__lfh.write("\n+SequenceInputRequest.getDeleteList() dList %s\n" % d)                
                dTup=d.split('|')
                #self.__lfh.write("\n+SequenceInputRequest.getDeleteList() dList %s %s\n" % (dTup[0],dTup[1]))
                oL.append((dTup[0].strip(),dTup[1].strip()))
        except:
            self.__lfh.write("\n+SequenceInputRequest.getDeleteList() error processing delete select list\n")
        return oL

    def getRequestPath(self):
        return (self._getStringValue('request_path'))

    def getSessionId(self):
        return (self._getStringValue('sessionid'))

    def getSessionPath(self):
        return (self._getStringValue('SessionPath'))    

    def getSemaphore(self):
        return (self._getStringValue('semaphore'))        

    def getEditMenuSelectList(self):
        """  Process the global edit menu select list --

             'targetId|priorVal|newValue|priorCss', 'targetId|priorVal|newValue|priorCss', 'targetId|priorVal|newValue|priorCss', ...,

             Return a list of tuples -
             
             [(targetId,priorValue,NewValue,priorCss),(,,),,,]
             
        """
        oL=[]
        #        try:
        selectString=self._getStringValue('selectedids')
        rowList=selectString.split(',')
        if (self.__verbose):
            self.__lfh.write("\n+SequenceInputRequest.getEditMenuSelectList() rowlist length %d : %r\n" % (len(rowList),rowList))      
        for row in rowList:
            tp=row.split('|')
            oL.append( (tp[0],tp[1],tp[2],tp[3] ))

            #        except:
            #            self.__lfh.write("\n+SequenceInputRequest.getEditMenuSelectList(self): failed processing select list\n")
        return oL
        
    def getSelectList(self):
        """  Return a list of tuples containing [(id,value,type),(id,value,type),(id,value),(),...]

             3 id,value pairs are returned per conflict record in the order ref, aligned, details.
        """
        oL=[]
        colLab=['reference','aligned','details']
        try:
            selectString=self._getStringValue('selectedids')
            rowList=selectString.split(',')
            if (self.__verbose):
                self.__lfh.write("\n+SequenceInputRequest.getSelectList() rowlist length %d : %r\n" % (len(rowList),rowList))      
            for row in rowList:
                idTupList=row.split('|')
                #
                #self.__lfh.write("row = %r\n" % row)                
                #self.__lfh.write("idtuplist = %r\n" % idTupList)            
                #if (self.__verbose):
                #    self.__lfh.write("\n+SequenceInputRequest.getSelectList() idTupList length %d : %r\n" % (len(idTupList),idTupList))
                
                ii=0
                for idTup in idTupList:
                    sTup=idTup.split("~")
                    oL.append( (sTup[0],sTup[1],colLab[ii]) )
                    ii+=1                    
                    #self.__lfh.write("idTup %r\n" % idTup)                    
                    #self.__lfh.write("sTup %r\n" % sTup)
                    #self.__lfh.write("length sTup %d\n" % len(sTup))                
                    #self.__lfh.write("\n+SequenceInputRequest.getSelectList() id %s value %s type %s\n" % (sTup[0],sTup[1],"static"))
                    #self.__lfh.write(" sTup[0] = %s\n" % sTup[0])
                    #self.__lfh.write(" sTup[1] = %s\n" % sTup[1])
                    #self.__lfh.write(" ii = %d\n" % ii)
                    #self.__lfh.write(" colLab  = %s\n" % colLab[ii])

        except:
            self.__lfh.write("\n+SequenceInputRequest.getSelectList() error processing select list\n")

        return oL


    def getAlignmentOp(self):
        tV=self._getStringValue('operation')
        if tV in ['re-align','reload','re-load']:
            return "realign"
        elif tV in ['reset']:
            return "reset"
        elif tV in ['load']:                
            return "load"
        elif tV in ['loadandstore']:                
            return "loadandstore"                        
        else:
            return "load"

    def getAlignmentOrdering(self):
        """  RPS, 20100727: Return ordering code as sourced from UI drop-down box via "viewalign_ordering" parameter.
        """
        tV=self._getStringValue('viewalign_order')
        #self.__lfh.write("\n+WebRequest.getAlignmentOrdering() got viewalign_order as: %s\n" % tV )
        
        if tV in ['auth-xyz-ref']:
            return "auth-xyz-ref"
        elif tV in ['auth-ref-xyz']:
            return "auth-ref-xyz"
        elif tV in ['xyz-auth-ref']:                
            return "xyz-auth-ref"
        elif tV in ['xyz-ref-auth']:                
            return "xyz-ref-auth" 
        elif tV in ['ref-auth-xyz']:                
            return "ref-auth-xyz"
        elif tV in ['ref-xyz-auth']:                
            return "ref-xyz-auth" 
        else:
            #self.__lfh.write("\n+WebRequest.getAlignmentOrdering() got nothing from front-end!!" )
            return "auth-xyz-ref"
            
        
        
        #return tV
        


    def getEditOp(self):
        tV=self._getStringValue('operation')
        if tV in ['edit']:
            return "edit"
        elif tV in ['reset']:
            return "reset"
        elif tV in ['undo']:
            return "undo"        
        elif tV in ['delete']:                
            return "delete"
        elif tV in ['global_edit_form']:
            return "global_edit_form"
        elif tV in ['global_edit']:
            return "global_edit"
        elif tV in ['global_edit_menu','globalmenuedit']:
            return "global_edit_menu"
        elif tV in ['move','edit_move']:
            return "move"                
        else:
            return "edit"

    def getAlignIdList(self):
        idString=self._getStringValue('alignids')
        if (self.__verbose):
            self.__lfh.write("\n+SequenceInputRequest() original idString %s\n" % idString)
        if ((len(idString) > 0) and idString.startswith("['") and idString.endswith("']")):
            tString = idString[2:-2]
        else:
            tString = idString
        if (self.__verbose):
            self.__lfh.write("\n+SequenceInputRequest() processed  idString %s\n" % tString)

        tList=tString.split(',')
        idL=[]
        for tId in tList:
            idL.append(tId.strip())
        return idL
        
    def getSummaryOp(self):
        tV=self._getStringValue('operation')
        if tV in ['load','reload']:
            return "load"
        elif tV in ['loadexample']:
            return "loadexample"
        elif tV in ['loadrcsb']:
            return "loadrcsb"        
        else:
            return "loadexample"

    def getSummaryAlignList(self):
        idString=self._getStringValue('alignids')
        if len(idString) > 0:
            idList = str(idString).split(',')
            return idList
        else:
            return []
    
    def getSummarySelectList(self):
        idString=self._getStringValue('selectids')
        if len(idString) > 0:
            idList = str(idString).split(',')
            return idList
        else:
            return []


    def getSessionObj(self):
        if (self.exists("SessionPath")):
            sObj=SessionManager(topPath=self._getStringValue("SessionPath"))
        else:
            sObj=SessionManager()
        sObj.setId(uid=self._getStringValue("sessionid"))
        return sObj

    def newSessionObj(self):
        if (self.exists("SessionPath")):
            sObj=SessionManager(topPath=self._getStringValue("SessionPath"))
        else:
            sObj=SessionManager()

        sessionId = self._getStringValue("sessionid")

        if (len(sessionId) > 0):
            sObj.setId(sessionId)
            sObj.makeSessionPath()
        else:
            sObj.assignId()
            sObj.makeSessionPath()
            self.setValue('sessionid',sObj.getId())
        
        return sObj

    def getSeqIdList(self):
        return []

#

class ResponseContent(object):

    def __init__(self, reqObj=None, verbose=False,log=sys.stderr):
        """
        Manage content items to be transfered as part of the
        the application response.
        
        """
        self.__verbose=verbose
        self.__lfh=log
        self.__reqObj=reqObj
        #
        self.__cD={}
        self.__setup()
        
    def __setup(self):
        """ Default response content is set here.
        """
        self.__cD['htmlcontent']=''
        self.__cD['errorflag']=False
        self.__cD['errortext']=''
        if self.__reqObj is not None:
            self.__cD['sessionid']=self.__reqObj.getSessionId()
            self.__cD['semaphore']=self.__reqObj.getSemaphore()
        else:
            self.__cD['sessionid']=''
            self.__cD['semaphore']=''

    def get(self):
        return self.__cD

    def setHtmlList(self,htmlList=[]):
        self.__cD['htmlcontent']=''.join(htmlList)

    def setHtmlText(self,htmlText=''):
        self.__cD['htmlcontent']=htmlText
        
    def setError(self,errMsg='',semaphore=''):
        self.__cD['errorflag']=True
        self.__cD['errortext']=errMsg
        self.__cD['semaphore']=semaphore

    def setStatusCode(self,aCode):
        self.__cD['statuscode']=aCode

    def setHtmlContentPath(self,aPath):
        self.__cD['htmlcontentpath']=aPath

  
    def setConflictReportPath(self,aPath):        
        self.__cD['conflictreportpath']=aPath

    def setConflictReportFlag(self,aFlag):        
        self.__cD['conflictreportflag']=aFlag        

    def setAlignTag(self,aTag):
        self.__cD['aligntag']   = aTag
        
    def setViewAlignOrder(self,aOrder):
        """  RPS, 20100727: Set code that determines display ordering of seq types for alignment interface.
        """
        self.__cD['viewalign_order']   = aOrder

    def setEditOp(self,aOp):
        self.__cD['editopid'] = aOp

    def setEditTimeStamp(self,tS):
        self.__cD['edittimestamp']=int(tS)

    def setAlignmentTimeStamp(self,tS):
        self.__cD['aligntimestamp']=int(tS)
        
    


if __name__ == '__main__':
    rC=ResponseContent()








