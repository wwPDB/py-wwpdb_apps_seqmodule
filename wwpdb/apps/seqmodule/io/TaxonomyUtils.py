##
# File:  TaxonomyIoUtils.py
# Date:  23-Feb-2013
#
# Updates:
#   24-Feb-2013  jdw add serialization
#   25-Feb-2013  jdw replace dependency request object with siteId on the request object
#   20-Mar-2013  jdw add method to lookup organism name from taxid
#   17-Oct-2022   zf add checking if pickle file exists
##
"""
Accessors for NCBI taxonomy names and organizational data hierarchy.

Python style serialization and deserialization of taxonomy data structures is provided
and coordinated within the project reference data storage model.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.07"


import os
import sys
import traceback

try:
    import cPickle as pickle
except ImportError:
    import pickle
from wwpdb.utils.config.ConfigInfoApp import ConfigInfoAppCommon


class TaxonomyUtils(object):
    """Accessors for NCBI taxonomy names and organizational data hierarchy.

    Python style serialization and deserialization of taxonomy data structures is provided
    and coordinated within the project reference data storage model.
    """

    def __init__(self, siteId="WWPDB_DEPLOY_TEST", verbose=True, log=sys.stderr):
        self.__verbose = verbose
        # self.__debug = True
        self.__lfh = log
        self.__siteId = siteId

        self.__cIAppCommon = ConfigInfoAppCommon(self.__siteId)
        self.__taxPath = self.__cIAppCommon.get_taxdump_path()

        self.__pickleProtocol = pickle.HIGHEST_PROTOCOL
        #
        self.__namesPicPath = os.path.join(self.__taxPath, "names.pic")
        self.__names = None
        #
        self.__nodesPicPath = os.path.join(self.__taxPath, "nodes.pic")
        self.__nodes = None

        self.__orgNameListPicPath = os.path.join(self.__taxPath, "org-name-list.pic")
        self.__nodes = None

    def lookUpSource(self, taxId):
        try:
            self.getNames()
            return self.__names[str(taxId)]["name"]
        except:  # noqa: E722 pylint: disable=bare-except
            return []

    def getParentDict(self):
        if not self.__nodes:
            self.__nodes = self.__deserialize(self.__nodesPicPath)
        return self.__nodes

    def getNames(self):
        if not self.__names:
            self.__names = self.__deserialize(self.__namesPicPath)
        return self.__names

    def __serialize(self, d, fn):
        try:
            fb = open(fn, "wb")
            pickle.dump(d, fb, self.__pickleProtocol)
            fb.close()
        except:  # noqa: E722 pylint: disable=bare-except
            pass
        #

    def __deserialize(self, fn):
        try:
            if os.access(fn, os.F_OK):
                fb = open(fn, "rb")
                d = pickle.load(fb)
                fb.close()
                if self.__verbose:
                    self.__lfh.write("+TaxonomyUtils.__deserialize() return %d records for file %s\n" % (len(d), fn))
                return d
            else:
                return {}
            #
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+TaxonomyUtils.__deserialize() failed for file %s\n" % fn)
            #
            return {}
        #

    def readTaxonomyNodeData(self, serialize=False):
        """Read the NCBI taxonomy data file 'nodes.dmp' and
        return a dictionary of parent taxonomy id's.

        d[tax_id]=parent_tax_id
        """
        d = {}
        try:
            fn = os.path.join(self.__taxPath, "nodes.dmp")
            if os.access(fn, os.F_OK):
                ifh = open(fn, "r")
                for line in ifh:
                    if not line:
                        continue
                    #
                    fields = line.split("\t|\t")
                    d[str(fields[0]).strip()] = str(fields[1]).strip()
                #
                ifh.close()
                if serialize:
                    self.__serialize(d, self.__nodesPicPath)
                #
            #
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
            #
            self.__lfh.write("+TaxonomyUtils.readTaxonomyNodeData() failed to read taxonomy data file %s\n" % fn)
        #
        return d

    def readTaxonomyNameData(self, serialize=False):
        """Read the NCBI taxonomy data file 'names.dmp' and
        return a dictionary of taxonomy names and synonyms.

        d[tax_id]=(name: [] , synoynm: [])

        """
        d = {}
        orgD = {}
        try:
            fn = os.path.join(self.__taxPath, "names.dmp")
            if os.access(fn, os.F_OK):
                ifh = open(fn, "r")
                for line in ifh:
                    if line is None or len(line) < 1:
                        continue
                    #
                    fields = line[:-1].split("\t|")
                    taxId = str(fields[0]).strip()
                    name = str(fields[1]).strip()
                    nType = str(fields[3]).strip()
                    #
                    if taxId not in d:
                        d[taxId] = {}
                        d[taxId]["name"] = []
                        d[taxId]["synonym"] = []
                    #
                    if nType in ["scientific name"]:
                        d[taxId]["name"].append(name)
                        orgD[name] = name
                    elif nType in ["synonym", "equivalent name"]:
                        d[taxId]["synonym"].append(name)
                    #
                #
                ifh.close()
                if serialize:
                    self.__serialize(d, self.__namesPicPath)
                    self.__serialize(sorted(orgD.keys()), self.__orgNameListPicPath)
                #
            #
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                traceback.print_exc(file=self.__lfh)
            #
            self.__lfh.write("+TaxonomyUtils.readTaxonomyNameData() failed to read taxonomy data file %s\n" % fn)
        #
        return d

    def getAncestors(self, taxId):
        """Fetch the ancestors for the input taxId.

        Return ancestors dictionary containing keys:

                id    - input TaxId
                p_id  - parent TaxId
               gp_id  - grand parent TaxId
        """
        if not self.__nodes:
            self.__nodes = self.__deserialize(self.__nodesPicPath)
        #
        anD = {}
        if taxId is None or len(taxId) < 1:
            return anD
        #
        try:
            sTaxId = str(taxId).strip()
            anD["id"] = sTaxId
            if self.__nodes:
                if sTaxId in self.__nodes:
                    parent_id = self.__nodes[sTaxId]
                    anD["p_id"] = parent_id
                    if parent_id in self.__nodes:
                        gparent_id = self.__nodes[parent_id]
                        anD["gp_id"] = gparent_id
                    #
                #
            #
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+TaxonomyUtils.getAncestos() failed for input taxId %s\n" % taxId)
                traceback.print_exc(file=self.__lfh)
            #
        #
        return anD
