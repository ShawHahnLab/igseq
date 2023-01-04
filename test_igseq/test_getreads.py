"""Tests for igseq.getreads."""

import re
import logging
from collections import defaultdict
from tempfile import TemporaryDirectory
from pathlib import Path
from unittest.mock import Mock
from igseq import getreads
from .util import TestBase, TestLive

class TestGetreads(TestBase):
    """Basic test of getreads.

    This fakes the bcl2fastq call but checks inputs and outputs and log
    messages.
    """

    def setUp(self):
        # pylint: disable=protected-access
        getreads._run_bcl2fastq_orig = getreads._run_bcl2fastq
        self.mock = self.make_mock_bcl2fastq()
        getreads._run_bcl2fastq = self.mock
        super().setUp()

    def tearDown(self):
        # pylint: disable=protected-access
        getreads._run_bcl2fastq = getreads._run_bcl2fastq
        super().tearDown()

    @classmethod
    def make_mock_bcl2fastq(cls):
        """Create a mock getreads._run_bcl2fastq object.

        When called it will make empty files with the expected filenames.
        """
        def side_effect(args, extra_args=None):
            idx = args.index("--output-dir") + 1
            path = Path(args[idx])
            path.mkdir(parents=True, exist_ok=True)
            for read in ["I1", "R1", "R2"]:
                (path / f"Undetermined_S0_L001_{read}_001.fastq.gz").touch(exist_ok=False)
            (path/"Stats").mkdir()
            with open(path/"Stats/FastqSummaryF1L1.txt", "wt") as f_out:
                f_out.write("SampleNumber\tTile\tNumberOfReadsRaw\tNumberOfReadsPF\n")
        return Mock(side_effect=side_effect)

    def test_getreads(self):
        """Test that getreads calls bcl2fastq with default options and output path."""
        with TemporaryDirectory() as tmpdir:
            with self.assertLogs(logger=getreads.LOGGER, level=logging.INFO) as logcm:
                getreads.getreads(self.path/"input/run", tmpdir)
                with open(Path(tmpdir)/"getreads.counts.csv") as f_in:
                    cts_obs = f_in.read()
                with open(self.path/"output/getreads.counts.csv") as f_in:
                    cts_exp = f_in.read()
                self.assertEqual(cts_obs, cts_exp)
                logcounts = defaultdict(int)
                for record in logcm.records:
                    logcounts[record.levelname] += 1
            args_exp = [
                '--runfolder-dir', self.path/'input/run',
                '--output-dir', Path(tmpdir),
                '--interop-dir', Path(tmpdir)/'InterOp',
                '--create-fastq-for-index-reads',
                '--barcode-mismatches', 0,
                '--sample-sheet', '/tmp/tmpqkpjltq5',
                '--loading-threads', 1,
                '--demultiplexing-threads', 1,
                '--processing-threads', 1,
                '--writing-threads', 1,
                '--min-log-level', 'ERROR']
        call_args = self.mock.call_args
        args_obs = call_args.args[0]
        def tmpmunge(txt):
            try:
                return re.sub("^/tmp/.*$", "/tmp/tmpfile", txt)
            except TypeError:
                return txt
        args_exp = [tmpmunge(txt) for txt in args_exp]
        args_obs = [tmpmunge(txt) for txt in args_obs]
        # The bcl2fastq function should have been called once with the expected
        # arguments (after putting a placeholder for any /tmp/... paths)
        self.mock.assert_called_once()
        self.assertEqual(args_exp, args_obs)
        # There should be some info message but nothing higher than that
        self.assertEqual(logcounts, {"INFO": 5})


class TestGetreadsMissingFiles(TestGetreads):
    """Test of getreads with fewer than three fastq.gz files out."""

    @classmethod
    def make_mock_bcl2fastq(cls):
        """Mock _run_bcl2fastq that doesn't make any of the expected output files."""
        def side_effect(args, extra_args=None):
            idx = args.index("--output-dir") + 1
            path = Path(args[idx])
            path.mkdir(parents=True, exist_ok=True)
        return Mock(side_effect=side_effect)

    def test_getreads(self):
        """Test that getreads complains when fastq.gz files are missing."""
        with TemporaryDirectory() as tmpdir:
            with self.assertLogs(logger=getreads.LOGGER, level=logging.INFO) as logcm:
                getreads.getreads(self.path/"input/run", tmpdir)
                logcounts = defaultdict(int)
                for record in logcm.records:
                    logcounts[record.levelname] += 1
        # The bcl2fastq function should have been called once
        self.mock.assert_called_once()
        # There should be some info messages and also errors about the missing
        # files
        self.assertEqual(logcounts, {"INFO": 5, "CRITICAL": 2})


class TestGetreadsExtraFiles(TestGetreads):
    """Test of getreads with more than three fastq.gz files out."""

    @classmethod
    def make_mock_bcl2fastq(cls):
        """Mock _run_bcl2fastq that create extra files."""
        def side_effect(args, extra_args=None):
            idx = args.index("--output-dir") + 1
            path = Path(args[idx])
            path.mkdir(parents=True, exist_ok=True)
            for read in ["I1", "R1", "R2"]:
                (path / f"Undetermined_S0_L001_{read}_001.fastq.gz").touch(exist_ok=False)
            for read in ["I1", "R1", "R2"]:
                (path / f"sample1_S1_L001_{read}_001.fastq.gz").touch(exist_ok=False)
            (path/"Stats").mkdir()
            with open(path/"Stats/FastqSummaryF1L1.txt", "wt") as f_out:
                f_out.write("SampleNumber\tTile\tNumberOfReadsRaw\tNumberOfReadsPF\n")
        return Mock(side_effect=side_effect)

    def test_getreads(self):
        """Test that getreads complains when fastq.gz files are missing."""
        with TemporaryDirectory() as tmpdir:
            with self.assertLogs(logger=getreads.LOGGER, level=logging.INFO) as logcm:
                getreads.getreads(self.path/"input/run", tmpdir)
                logcounts = defaultdict(int)
                for record in logcm.records:
                    logcounts[record.levelname] += 1
        # The bcl2fastq function should have been called once
        self.mock.assert_called_once()
        # There should be some info messages and also a warning about the extra
        # files
        self.assertEqual(logcounts, {"INFO": 5, "WARNING": 1})


class TestGetreadsLive(TestGetreads, TestLive):
    """Basic test of getreads with actual bcl2fastq."""

    def test_getreads(self):
        # this would need an actual run directory tree to use
        self.skipTest("not yet implemented")
