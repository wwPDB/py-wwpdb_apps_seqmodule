##
# File: WorkerCountTests.py
#
"""Worker-process counts used by the sequence module searches.

Without configuration the pool size is half the CPUs the process may use (at least 1); the two searches that run
per entity use one worker per entity when the compute cluster is enabled.  The resource allocation file can set a
count for each call site.  The tests replace the process pool, the site configuration and the CPU counts, so nothing
real is started.
"""

# pylint: disable=protected-access,unused-argument,no-member,abstract-method

import json
import multiprocessing
import os
import sys
import tempfile
import unittest
from unittest import mock

try:
    import MySQLdb  # noqa: F401  # pylint: disable=unused-import
except ImportError:
    # The database driver is not needed here; a stand-in lets the modules under test import.
    _MYSQL = mock.MagicMock()
    sys.modules["MySQLdb"] = _MYSQL
    sys.modules["MySQLdb.cursors"] = _MYSQL.cursors

import wwpdb.apps.seqmodule.control.SequenceDataAssemble_v2 as assemble_module  # noqa: E402
import wwpdb.apps.seqmodule.util.LocalBlastSearchUtils as blast_module  # noqa: E402


class _FakePool(object):
    """Stands in for MultiProcUtil and remembers the requested number of processes."""

    requested = []

    def __init__(self, *args, **kwargs):
        pass

    def set(self, **kwargs):
        pass

    def setWorkingDir(self, path):
        pass

    def setOptions(self, options):
        pass

    def runMulti(self, dataList=None, numProc=None, numResults=None):
        _FakePool.requested.append(numProc)
        return True, [], [[]], []


class _FakeConfigInfo(object):
    """Site configuration with a resource allocation file and an optional compute-cluster flag."""

    path = None
    cluster = False

    def __init__(self, *args, **kwargs):
        pass

    def get(self, key, default=None):
        if key == "RESOURCE_ALLOCATION_CONFIG_FILE":
            return _FakeConfigInfo.path
        if key == "USE_COMPUTE_CLUSTER":
            return _FakeConfigInfo.cluster
        return default


class _FakeSequenceData(object):
    @staticmethod
    def toList(sequence):
        return list(sequence)


class _Base(object):
    key = None
    # Whether the call site uses one worker per entity when the compute cluster is enabled.
    per_entity_on_cluster = False

    def setUp(self):
        self._tmp = tempfile.mkdtemp()
        _FakePool.requested = []
        _FakeConfigInfo.path = None
        _FakeConfigInfo.cluster = False
        mock.patch("wwpdb.io.misc.ResourceAllocation.ConfigInfo", _FakeConfigInfo).start()
        self.addCleanup(mock.patch.stopall)

    def _configure(self, jobs=None, defaults=None):
        path = os.path.join(self._tmp, "resources.yaml")
        with open(path, "w") as out:
            json.dump({"jobs": jobs or {}, "defaults": defaults or {}}, out)
        _FakeConfigInfo.path = path

    def _cpus(self, available, node=None):
        mock.patch.object(os, "sched_getaffinity", lambda pid: set(range(available)), create=True).start()
        mock.patch.object(multiprocessing, "cpu_count", lambda: available if node is None else node).start()

    def _pool_size(self):
        raise NotImplementedError

    def testUnconfiguredUsesHalfTheCpus(self):
        self._cpus(16)
        self.assertEqual(8, self._pool_size())

    def testUnconfiguredUsesCpusAvailableToTheProcessNotTheNode(self):
        self._cpus(4, node=64)
        self.assertEqual(2, self._pool_size())

    def testOneCpuStillStartsOneWorker(self):
        self._cpus(1)
        self.assertEqual(1, self._pool_size())

    def testConfiguredCountIsUsed(self):
        self._cpus(16)
        self._configure(jobs={self.key: {"cpus": 3}})
        self.assertEqual(3, self._pool_size())

    def testConfiguredCountAboveAvailableCpusIsClamped(self):
        self._cpus(4)
        self._configure(jobs={self.key: {"cpus": 12}})
        self.assertEqual(4, self._pool_size())

    def testMalformedConfigurationFallsBackToTheDefault(self):
        self._cpus(16)
        self._configure(jobs={self.key: {"cpus": "lots"}})
        self.assertEqual(8, self._pool_size())


class _PerEntityBase(_Base):
    """Call sites that start one worker per entity on the compute cluster."""

    def testClusterUsesOneWorkerPerEntityEvenAboveTheCpus(self):
        _FakeConfigInfo.cluster = True
        self._cpus(2)
        self.assertEqual(5, self._pool_size())

    def testConfiguredCountIsUsedWithCluster(self):
        _FakeConfigInfo.cluster = True
        self._cpus(16)
        self._configure(jobs={self.key: {"cpus": 3}})
        self.assertEqual(3, self._pool_size())

    def testMalformedConfigurationWithClusterFallsBackToOnePerEntity(self):
        _FakeConfigInfo.cluster = True
        self._cpus(16)
        self._configure(jobs={self.key: {"cpus": "lots"}})
        self.assertEqual(5, self._pool_size())


class LocalBlastWorkerCountTests(_Base, unittest.TestCase):
    key = "seqmodule.blast"

    def _pool_size(self):
        search = object.__new__(blast_module.LocalBlastSearchUtils)
        search._LocalBlastSearchUtils__siteId = "TEST"
        search._LocalBlastSearchUtils__verbose = False
        search._LocalBlastSearchUtils__lfh = sys.stderr
        search._LocalBlastSearchUtils__sessionPath = self._tmp
        search._LocalBlastSearchUtils__dataSetId = "D_1"
        search._LocalBlastSearchUtils__entityId = "1"
        search._LocalBlastSearchUtils__ncbilock = None
        search._LocalBlastSearchUtils__srd = _FakeSequenceData()
        search._LocalBlastSearchUtils__entityD = {
            "SEQ_ENTITY_1_CAN": "MKVLSGIVGLLLQWERTYPE",
            "POLYMER_LINKING_TYPE": "polypeptide(L)",
            "PART_LIST": [{"SEQ_PART_TYPE": "biological sequence", "SEQ_PART_ID": 1, "SEQ_NUM_BEG": 1, "SEQ_NUM_END": 20, "SOURCE_TAXID": "9606"}],
        }
        _FakePool.requested = []
        with mock.patch.object(blast_module, "MultiProcUtil", _FakePool), mock.patch.object(
            blast_module.LocalBlastSearchUtils, "_LocalBlastSearchUtils__fetchAuthProvidedRefSequence", lambda *args, **kwargs: (False, False, [])
        ):
            search._LocalBlastSearchUtils__runBlastSearch()
        self.assertEqual(1, len(_FakePool.requested))
        return _FakePool.requested[0]


class ReferenceAlignmentWorkerCountTests(_Base, unittest.TestCase):
    key = "seqmodule.ref_alignments"

    def _pool_size(self):
        search = object.__new__(blast_module.LocalBlastSearchUtils)
        search._LocalBlastSearchUtils__siteId = "TEST"
        search._LocalBlastSearchUtils__verbose = False
        search._LocalBlastSearchUtils__lfh = sys.stderr
        search._LocalBlastSearchUtils__sessionPath = self._tmp
        _FakePool.requested = []
        with mock.patch.object(blast_module, "MultiProcUtil", _FakePool):
            search._LocalBlastSearchUtils__getRefSeqAlignments([("hit",)])
        self.assertEqual(1, len(_FakePool.requested))
        return _FakePool.requested[0]


def _assemble(tmp):
    assembler = object.__new__(assemble_module.SequenceDataAssemble)
    assembler._verbose = False
    assembler._lfh = sys.stderr
    assembler._SequenceDataAssemble__cI = _FakeConfigInfo()
    assembler._SequenceDataAssemble__sessionPath = tmp
    return assembler


class SameSequenceAnnotationWorkerCountTests(_PerEntityBase, unittest.TestCase):
    key = "seqmodule.same_seq_annotation"

    def _pool_size(self):
        _FakePool.requested = []
        with mock.patch.object(assemble_module, "MultiProcUtil", _FakePool):
            _assemble(self._tmp)._SequenceDataAssemble__runSameSeqAnnotationSearch([(str(i), "SEQ", "1") for i in range(5)])
        self.assertEqual(1, len(_FakePool.requested))
        return _FakePool.requested[0]


class SequenceBlastSearchWorkerCountTests(_PerEntityBase, unittest.TestCase):
    key = "seqmodule.blast"

    def _pool_size(self):
        _FakePool.requested = []
        with mock.patch.object(assemble_module, "MultiProcUtil", _FakePool), mock.patch.object(assemble_module, "MultiProcLimit", mock.Mock()):
            _assemble(self._tmp)._SequenceDataAssemble__runSeqBlastSearch([(str(i), "SEQ", "1") for i in range(5)])
        self.assertEqual(1, len(_FakePool.requested))
        return _FakePool.requested[0]


if __name__ == "__main__":
    unittest.main()
