##
# File:  TaxonomyDbUtils.py
# Date:  11-Oct-2021
#
##
"""
Get NCBI taxonomy names from "taxonomy" table in "status" database

"""
__docformat__ = "restructuredtext en"
__author__ = "Zukang Feng"
__email__ = "zfeng@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.07"


import sys
from wwpdb.utils.wf.dbapi.WfDbApi import WfDbApi


class TaxonomyDbUtils(object):
    """Get NCBI taxonomy names from "taxonomy" table in "status" database"""

    def __init__(self, verbose=True, log=sys.stderr):
        self.__verbose = verbose
        self.__lfh = log  # pylint:  disable=unused-private-member
        self.__wfApi = WfDbApi(verbose=self.__verbose)

    def getTaxonomyNames(self, taxId):
        """ """
        try:
            scientific_name = ""
            common_name = ""
            sql = "select tax_id, name, class from taxonomy where tax_id = '" + taxId + "' and class in ( 'scientific name', 'genbank common name' )"
            rows = self.__wfApi.runSelectSQL(sql)
            for row in rows:
                if (str(row[0]).strip() == taxId) and (row[2].strip().lower() == "scientific name"):
                    scientific_name = str(row[1].strip())
                #
                if (str(row[0]).strip() == taxId) and (row[2].strip().lower() == "genbank common name"):
                    common_name = str(row[1].strip())
                #
            #
            return scientific_name, common_name
        except:  # noqa: E722 pylint: disable=bare-except
            return "", ""
        #


def __main():
    taxUtil = TaxonomyDbUtils()
    scientific_name, common_name = taxUtil.getTaxonomyNames(sys.argv[1])
    print("scientific_name=%s" % scientific_name)
    print("common_name=%s" % common_name)


if __name__ == "__main__":
    __main()
