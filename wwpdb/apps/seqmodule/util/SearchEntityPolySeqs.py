##
# File:  SearchEntityPolySeqs.py
# Date:  21-Mar-2018
# Updates:
##
"""

This software was developed as part of the World Wide Protein Data Bank
Common Deposition and Annotation System Project

Copyright (c) 2018 wwPDB

This software is provided under a Creative Commons Attribution 3.0 Unported
License described at http://creativecommons.org/licenses/by/3.0/.

"""
__docformat__ = "restructuredtext en"
__author__ = "Zukang Feng"
__email__ = "zfeng@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.07"

# try:
#     import cPickle as pickle
# except ImportError:
#     import pickle as pickle
#
# import os
import sys

from wwpdb.utils.config.ConfigInfo import ConfigInfo
from wwpdb.utils.wf.dbapi.DbApiUtil import DbApiUtil


class SearchEntityPolySeqs(object):
    #
    __schemaMap = {
        "SELECT_ALL": "select r.structure_id, r.pdb_id, r.rcsb_annotator, r.status_code, r.date_begin_processing, e.entity_id, "
        + "e.pdbx_seq_one_letter_code from rcsb_status as r, entity_poly as e where r.structure_id = e.structure_id "
        + "and r.status_code in ( 'REL', 'AUTH', 'WAIT', 'HOLD', 'HPUB' ) and r.date_begin_processing is not null "
        + "and r.date_begin_processing != '0000-00-00' order by r.date_begin_processing desc",
        "SELECT_BY_SEQ": "select r.structure_id, r.pdb_id, r.rcsb_annotator, r.status_code, r.date_begin_processing, e.entity_id from "
        + "rcsb_status as r, entity_poly as e where r.structure_id = e.structure_id and ( r.status_code in ( 'REL', 'AUTH', 'HOLD', 'HPUB', 'WAIT', 'PROC' ) "
        + "or r.structure_id = '%s' ) and r.date_begin_processing is not null and r.date_begin_processing != '0000-00-00' and "
        + "replace(replace(replace(e.pdbx_seq_one_letter_code, '\n', ''), '\t', ''), ' ' ,'') = '%s' order by r.date_begin_processing desc",
    }
    #

    def __init__(self, siteId=None, sessionPath=None, verbose=False, log=sys.stderr):
        """ """
        self.__siteId = siteId
        self.__sessionPath = sessionPath  # pylint: disable=unused-private-member
        self.__verbose = verbose
        self.__lfh = log
        #
        self.__cI = ConfigInfo(self.__siteId)
        self.__dbServer = self.__cI.get("SITE_DB_SERVER")
        self.__dbHost = self.__cI.get("SITE_DB_HOST_NAME")
        self.__dbName = "da_internal"
        self.__dbUser = self.__cI.get("SITE_DB_USER_NAME")
        self.__dbPw = self.__cI.get("SITE_DB_PASSWORD")
        self.__dbSocket = self.__cI.get("SITE_DB_SOCKET")
        self.__dbPort = int(self.__cI.get("SITE_DB_PORT_NUMBER"))
        #
        self.__dbApi = DbApiUtil(
            dbServer=self.__dbServer,
            dbHost=self.__dbHost,
            dbName=self.__dbName,
            dbUser=self.__dbUser,
            dbPw=self.__dbPw,
            dbSocket=self.__dbSocket,
            dbPort=self.__dbPort,
            verbose=self.__verbose,
            log=self.__lfh,
        )
        self.__dbApi.setSchemaMap(self.__schemaMap)
        #

    def searchSameSequenceEntity(self, entryId="D_1000000000", seq=""):
        """ """
        seq = str(seq).replace("\n", "").replace("\t", "").replace(" ", "").strip().upper()
        rows = self.__dbApi.selectData(key="SELECT_BY_SEQ", parameter=(entryId, seq))
        #
        selfList = []
        otherList = []
        for row in rows:
            hasMissingValue = False
            metaDataList = []
            for item in ("structure_id", "entity_id", "pdb_id", "rcsb_annotator", "status_code", "date_begin_processing"):
                if (item not in row) or (not row[item]):
                    hasMissingValue = True
                    break
                #
                metaDataList.append(str(row[item]).strip())
            #
            if hasMissingValue:
                continue
            #
            if metaDataList[0] == entryId:
                selfList.append(metaDataList)
            else:
                otherList.append(metaDataList)
            #
        #
        return selfList, otherList

    # Commented out as unused.  June 2022
    # def __getSeqMap(self):
    #     """ """
    #     self.__pickleFileName = "EntityPolySeqs.pickle"
    #     #
    #     pickleFile = os.path.join(self.__sessionPath, self.__pickleFileName)
    #     if os.access(pickleFile, os.F_OK):
    #         self.__seqMap = self.__loadPickle(pickleFile)
    #         return
    #     #
    #     rows = self.__dbApi.selectData(key="SELECT_ALL", parameter=())
    #     for row in rows:
    #         hasMissingValue = False
    #         for item in ("structure_id", "entity_id", "pdb_id", "rcsb_annotator", "status_code", "date_begin_processing", "pdbx_seq_one_letter_code"):
    #             if (item not in row) or (not row[item]):
    #                 hasMissingValue = True
    #                 break
    #             #
    #         #
    #         if hasMissingValue:
    #             continue
    #         #
    #         metaDataList = []
    #         for item in ("structure_id", "entity_id", "pdb_id", "rcsb_annotator", "status_code", "date_begin_processing"):
    #             val = ""
    #             if (item in row) and row[item]:
    #                 val = str(row[item]).strip()
    #             #
    #             metaDataList.append(val)
    #         #
    #         seq = str(row["pdbx_seq_one_letter_code"]).replace("\n", "").replace("\t", "").replace(" ", "").strip().upper()
    #         seqLen = len(seq)
    #         if seqLen in self.__seqMap:
    #             if seq in self.__seqMap[seqLen]:
    #                 self.__seqMap[seqLen][seq].append(metaDataList)
    #             else:
    #                 self.__seqMap[seqLen][seq] = [metaDataList]
    #             #
    #         else:
    #             self.__seqMap[seqLen] = {seq: [metaDataList]}
    #         #
    #     #
    #     self.__dumpPickle(pickleFile, self.__seqMap)

    # def __loadPickle(self, pickleFile):
    #     """ """
    #     fb = open(pickleFile, "rb")
    #     pickleData = pickle.load(fb)
    #     fb.close()
    #     return pickleData

    # def __dumpPickle(self, pickleFile, pickleData):
    #     fb = open(pickleFile, "wb")
    #     pickle.dump(pickleData, fb)
    #     fb.close()
