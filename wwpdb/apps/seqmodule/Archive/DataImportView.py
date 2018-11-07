##
# File:    DataImportView.py
# Date:    3-March-2010
#
# Update:
#  3-Mar-2010  jdw Split off data import operations from SummaryView Class
#  5-Mar-2010  jdw Add file source 'repository|upload' option
# 20-Apr-2010  jdw Ported to module seqmodule
##

"""
Controlling class for data import view - includes sequence alignment statistics generation -

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"


import sys, os.path, shutil
from wwpdb.apps.seqmodule.io.SequenceDataImport     import SequenceDataImportRcsb,SequenceDataImportExampleRcsb
from wwpdb.apps.seqmodule.align.AlignmentStatistics import AlignmentStatistics

class DataImportView(object):
    """ Controlling class for data import view.

        Supported operations:
        
         loadrcsb:        create sequence data store with data for an rcsb entry and then load summary.
         loadexample:     create test data set from internal example class and then load summary
         loadrcsbfile:    upload file option - 
            
    """
    def __init__(self,reqObj=None,verbose=False,log=sys.stderr):
        self.__verbose=verbose
        self.__reqObj=reqObj
        self.__lfh=log
        #
        self.__operation  = None
        self.__sessionObj = None
        self.__identifier = None
        #
        if (self.__verbose):
            self.__lfh.write("+DataImportView() starting\n")
            self.__lfh.flush()
        #
        self.__setup()
        #

    def __setup(self):
        try:
            self.__sessionObj  = self.__reqObj.getSessionObj()
            self.__sessionPath = self.__sessionObj.getPath()
            self.__operation   = self.__reqObj.getValue("operation")
            self.__identifier  = str(self.__reqObj.getValue("identifier")).upper()
            #
            if (self.__verbose):
                self.__lfh.write("+SeqToolWebApp.__setup() operation  %s\n" % self.__reqObj.getValue("operation"))
                self.__lfh.write("+SeqToolWebApp.__setup() identifier %s\n" % self.__reqObj.getValue("identifier"))
                self.__lfh.write("+SeqToolWebApp.__setup() UploadFileType  %s\n" % self.__reqObj.getValue("UploadFileType"))
                self.__lfh.write("+SeqToolWebApp.__setup() UploadFileName  %s\n" % self.__reqObj.getValue("UploadFileName"))
                #
                self.__lfh.write("+DataImportView.__setup() operation  %s\n" %  str(self.__operation))
                self.__lfh.write("+DataImportView.__setup() identifier %s\n" %  str(self.__identifier))
                self.__lfh.flush()                
        except:
            if (self.__verbose):
                self.__lfh.write("+DataImportView.__setup() with operation %s : sessionId %s failed\n" %
                                 (self.__operation,self.__sessionObj.getId()))                                         
    def loadExampleData(self):
        return(self.__loadExampleData())

    def loadData(self):
        if (self.__verbose):
            self.__lfh.write("+DataImportView.loadData() with operation %s : sessionId %s\n" %
                             (self.__operation,self.__sessionObj.getId()))
            self.__lfh.flush()

        if self.__operation == "loadrcsb":
            rcsbId=str(self.__reqObj.getValue("identifier")).upper()
            self.__loadRcsbData(rcsbId,fileSource='repository')
            alstat=AlignmentStatistics(sessionObj=self.__sessionObj,verbose=self.__verbose,log=self.__lfh)
            alstat.doUpdate()
            return {}

        elif self.__operation == "loadrcsbfile":
            rcsbId=str(self.__reqObj.getValue("identifier")).upper()
            self.__loadRcsbData(rcsbId,fileSource='upload')
            alstat=AlignmentStatistics(sessionObj=self.__sessionObj,verbose=self.__verbose,log=self.__lfh)
            alstat.doUpdate()
            return {}

        elif self.__operation == "loadexample":
            self.__loadExampleData()
            alstat=AlignmentStatistics(sessionObj=self.__sessionObj,verbose=self.__verbose,log=self.__lfh)
            alstat.doUpdate()                        
            return {}

        elif self.__operation == "loadexampledata":
            return(self.__loadExampleData())

    def __loadExampleData(self):
        if (self.__verbose):
            self.__lfh.write("+DataImportView.__loadExampleData() sessionId %s\n" % self.__sessionObj.getId())
            self.__lfh.write("+DataImportView.__loadExampleData() operation %s\n" % self.__operation)
            
        sdi=SequenceDataImportExampleRcsb(reqObj=self.__reqObj,verbose=self.__verbose,log=self.__lfh)
        sdi.doImport()

    def __loadRcsbData(self,rcsbId,fileSource='repository'):
        if (self.__verbose):
            self.__lfh.write("+DataImportView.__loadRcsbData() with Op %s : rcsb id %s file source %s\n" %
                             (self.__operation,rcsbId,fileSource ))
            self.__lfh.flush()
            
        sdi=SequenceDataImportRcsb(reqObj=self.__reqObj,rcsbId=rcsbId,fileSource=fileSource,verbose=self.__verbose,log=self.__lfh,fetchPdbFile=True)
        sdi.doImport()
        
        
                
if __name__ == '__main__':
    pass
