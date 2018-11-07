##
# File:    SessionManager.py
# Date:    14-Dec-2009
#
# Updates:
# 20-Apr-2010 jdw Ported to module seqmodule.
##
"""
Provides containment and access for session information.  Methods
are provided to create temporary directories to preserve session files.

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"


import sys, sha, time, os.path, shutil


class SessionManager(object):
    """
        Utilities for session directory maintenance.
        
    """
    def __init__(self,topPath=".",verbose=False):

        self.__verbose = verbose
        self.__topSessionPath = topPath
        self.__uid=None

    def setId(self,uid):
        self.__uid=uid

    def getId(self):
        return self.__uid
        
    def assignId(self):
        self.__uid = sha.new(repr(time.time())).hexdigest()
        return self.__uid        

    def getPath(self):
        try:
            pth=os.path.join(self.__topSessionPath,self.__uid)
            if (self.__verbose):
                sys.stderr.write("+SessionManager.getPath() path %s\n" % pth)
            if os.access(pth,os.F_OK):
                return pth
            else:
                return None
        except:
            return None

    def getTopPath(self):
        return self.__topSessionPath
    
    def makeSessionPath(self):
        try:
            pth=os.path.join(self.__topSessionPath,self.__uid)        
            if (not os.access(pth,os.F_OK)):
                os.mkdir(pth)
            return pth
        except:
            return None

    def remakeSessionPath(self):
        try:
            pth=os.path.join(self.__topSessionPath,self.__uid)        
            if (os.access(pth,os.F_OK)):
                shutil.rmtree(pth,True)
            os.mkdir(pth)
            return pth
        except:
            return None

