##
# File: ImportTests.py
# Date:  06-Oct-2018  E. Peisach
#
# Updates:
##
"""Test cases for sequence module"""

__docformat__ = "restructuredtext en"
__author__ = "Ezra Peisach"
__email__ = "peisach@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"

import os
import sys
import unittest

if __package__ is None or __package__ == "":
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from commonsetup import HERE  # noqa:  F401 pylint: disable=import-error,unused-import
else:
    from .commonsetup import HERE  # noqa: F401 pylint: disable=relative-beyond-top-level

from wwpdb.apps.seqmodule.webapp.SeqModWebApp_v2 import SeqModWebApp


class ImportTests(unittest.TestCase):
    def setUp(self):
        pass

    def testInstantiate(self):
        """Tests simple instantiation"""
        SeqModWebApp()
