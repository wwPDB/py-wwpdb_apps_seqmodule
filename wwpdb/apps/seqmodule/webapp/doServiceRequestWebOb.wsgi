##
# File:     doServiceRequestWebOb.wsgi     
# Created:  16-Dec-2009
# Updates:
# 14-Jan-2010  Refactored to better decouple web services
#              handler and backend application functionality.
#
# 15-Feb-2010  Set web application path and verbose flag in
#                        this module for all applications.
#  2-Mar-2010  Integrate path with parameter dictionary.                        
#
# 20-Apr-2010 Ported to seqmodule package
##
"""
This top-level responder for requests to /services/.... url for the
wwPDB sequence editor application framework.

This version depends on WSGI and WebOb.

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"

import sys

import traceback
from webob import Request, Response

#  - URL mapping and application specific classes are launched from SeqModWebApp()
from wwpdb.apps.seqmodule.webapp.SeqModWebApp import SeqModWebApp

class MyRequestApp(object):
    """  Handle server interaction using WSGI and WebOb Request
         and Response objects.
    """
    def __init__(self,textString="Initialized from contructor",verbose=False,log=sys.stderr):
        """ 
        """
        self.__text=textString
        self.__verbose=verbose
        self.__debug=False
        self.__lfh=log
        self.__siteId=None
        self._myParameterDict={}        
        
    def __dumpEnv(self,request):
        outL=[]
        outL.append("\n+MyRequestApp.__dumpEnv() ------------------------------------------------\n")
        outL.append("Web server request data content:\n")                
        try:
            outL.append("Text initialization:   %s\n" % self.__text)        
            outL.append("Host:         %s\n" % request.host)
            outL.append("Path:         %s\n" % request.path)
            outL.append("Method:       %s\n" % request.method)        
            outL.append("Query string: %s\n" % request.query_string)
            outL.append("Parameter List:\n")
            for name,value in request.params.items():
                outL.append("Request parameter:    %s:  %r\n" % (name,value))
        except:
            traceback.print_exc(file=self.__lfh)
        outL.append("------------------------------------------------\n")
        return outL

    def __call__(self, environment, responseApplication):
        """          WSGI callable entry point


        """
        myRequest  = Request(environment)
        #
        self._myParameterDict={}   
        try:
            if environment.has_key('WWPDB_SITE_ID'):
                self.__siteId=environment['WWPDB_SITE_ID']
                if self.__debug:
                    self.__lfh.write("+MyRequestApp.__call__() - WWPDB_SITE_ID environ variable captured as %s\n" % self.__siteId)
            #
            if (self.__debug):
                for name,value in environment.items():
                    self.__lfh.write("+MyRequestApp.__call__() - request environment:    %s:  %r\n" % (name,value))

            for name,value in myRequest.params.items():
                if (not self._myParameterDict.has_key(name)):
                    self._myParameterDict[name]=[]
                self._myParameterDict[name].append(value)
            self._myParameterDict['request_path']=[myRequest.path.lower()]

            if environment.has_key('HTTP_HOST'):            
                self._myParameterDict['request_host']=[environment['HTTP_HOST']]
            else:
                self._myParameterDict['request_host']=['']
        except:
            self.__lfh.write("+MyRequestApp.__call__() - Exception processing in request setup\n")            
            traceback.print_exc(file=self.__lfh)
            self.__lfh.write("+MyRequestApp.__call__() - contents of request data\n")
            self.__lfh.write("%s" % ("".join(self.__dumpEnv(request=myRequest))))
            
        ###
        ### At this point we have everything needed from the request !
        ###
        myResponse = Response()
        myResponse.status       = '200 OK'
        myResponse.content_type = 'text/html'       
        ###
        ###  Application specific functionality called here --
        ###  Application receives path and parameter info only!
        ###
        seqT= SeqModWebApp(parameterDict=self._myParameterDict,verbose=self.__verbose,log=self.__lfh,siteId=self.__siteId)
        rspD=seqT.doOp()
        myResponse.content_type=rspD['CONTENT_TYPE']
        myResponse.body=rspD['RETURN_STRING']
        if rspD.has_key('ENCODING'):
            myResponse.content_encoding=rspD['ENCODING']
        if rspD.has_key('DISPOSITION'):
            myResponse.content_disposition=rspD['DISPOSITION']

        if (self.__debug):
            self.__lfh.write("+MyRequestApp.__call__() - Response content_type %r \n" % rspD['CONTENT_TYPE'])
            self.__lfh.write("+MyRequestApp.__call__() - Response body  %r \n" % rspD['RETURN_STRING'])       
        ####
        ###
        return myResponse(environment,responseApplication)
##
##  NOTE -  Path to top of the web application tree and verbose setting are set here ONLY! 
##
application = MyRequestApp(textString="doServiceRequest() - WebOb version",verbose=True,log=sys.stderr)













