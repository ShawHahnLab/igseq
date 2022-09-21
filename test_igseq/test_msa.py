import sys
from abc import ABC
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import Mock
from igseq import msa
from igseq.record import RecordReader
from .util import TestBase, TestLive


def setup_mock_popen(orig, recs_in, recs_out):
    """Replace subprocess.Popen with fake version for testing."""
    mock = Mock(spec=orig)
    mock.return_value = mock
    mock.__enter__ = mock
    mock.__exit__ = Mock(return_value=False)
    mock.returncode = None
    txt_in = [f">{rec['sequence_id']}\n{rec['sequence']}" for rec in recs_in]
    mock.expected_input = "\n".join(txt_in) + "\n"
    txt_out = [f">{rec['sequence_id']}\n{rec['sequence']}" for rec in recs_out]
    mock.expected_output = "\n".join(txt_out) + "\n"
    def popen_mock_exit(*_, **__):
        # when exiting, Popen() waits for the process so it should now have
        # an exit code
        sys.stderr.write("\nmuscle 5")
        mock.returncode = 0
    def popen_mock_check_input(*_, **__):
        # This mock gives output assuming this particular input, so check
        # for that.
        if mock.stdin.getvalue() != mock.expected_input:
            raise ValueError("Input FAST mismatch for mock")
    mock.__exit__.side_effect = popen_mock_exit
    mock.stdin = StringIO()
    mock.stderr = StringIO("""
muscle 5.1.linux64 []  396Gb RAM, 56 cores
Built Feb 24 2022 03:16:15
(C) Copyright 2004-2021 Robert C. Edgar.
https://drive5.com

Input: 2 seqs, avg length 30, max 30

00:00 18Mb   CPU has 56 cores, defaulting to 20 threads
00:00 177Mb   100.0% Calc posteriors
00:00 177Mb   100.0% UPGMA5""")
    mock.stdin.close = Mock(side_effect=popen_mock_check_input)
    mock.stdout = StringIO(mock.expected_output)
    mock.popen_real = orig
    return mock


class CommonMSATests(ABC):
    """Shared test logic for actual test case classes below."""

    def setUp(self):
        super().setUp()
        with RecordReader(self.path/"input.fasta") as reader:
            self.records_in = list(reader)
        with RecordReader(self.path/"output.fasta") as reader:
            self.records_out_exp = list(reader)
        self.records_out = None

    def test_msa(self):
        """Test msa with file input and output."""
        with TemporaryDirectory() as tmpdir:
            stdout, stderr = self.redirect_streams(lambda:
                msa.msa(self.path/"input.fasta", Path(tmpdir)/"output.fasta"))
            self.assertEqual("", stdout)
            self.assertTrue(stderr.startswith("\nmuscle 5"))
            self.assertTxtsMatch(
                self.path/"output.fasta",
                Path(tmpdir)/"output.fasta")

    def test_run_muscle(self):
        """Test muscle alignment with record inputs and outputs."""
        def do_muscle(self, recs):
            self.records_out = msa.run_muscle(recs)
        stdout, stderr = self.redirect_streams(lambda: do_muscle(self, self.records_in))
        self.assertEqual("", stdout)
        self.assertTrue(stderr.startswith("\nmuscle 5"))
        self.assertEqual(self.records_out, self.records_out_exp)


class TestMSA(CommonMSATests, TestBase):
    """Basic test of msa, with a Mock MUSCLE call."""

    def setUp(self):
        super().setUp()
        msa.Popen = setup_mock_popen(msa.Popen, self.records_in, self.records_out_exp)

    def tearDown(self):
        super().tearDown()
        msa.Popen = msa.Popen.popen_real


class TestMSALive(CommonMSATests, TestBase, TestLive):
    """Basic test of msa with actual MUSCLE calls."""
