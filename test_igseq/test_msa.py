from abc import ABC
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from igseq import msa
from igseq.record import RecordReader
from .util import TestBase, TestLive, MockPopen


class MockPopenMuscle(MockPopen):

    def setup_texts(self, kwargs):
        texts = {"stdin": "", "stdout": "", "stderr": ""}
        for key in texts.keys():
            try:
                texts[key] = kwargs[key]
                del kwargs[key]
            except KeyError:
                texts[key] = None
        if texts["stderr"] is None:
            texts["stderr"] = """
muscle 5.1.linux64 []  396Gb RAM, 56 cores
Built Feb 24 2022 03:16:15
(C) Copyright 2004-2021 Robert C. Edgar.
https://drive5.com

Input: 2 seqs, avg length 30, max 30

00:00 18Mb   CPU has 56 cores, defaulting to 20 threads
00:00 177Mb   100.0% Calc posteriors
00:00 177Mb   100.0% UPGMA5"""
        for key, val in texts.items():
            if val is None:
                texts[key] = ""
        return texts


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
        txt_in = [f">{rec['sequence_id']}\n{rec['sequence']}" for rec in self.records_in]
        txt_in = "\n".join(txt_in) + "\n"
        txt_out = [f">{rec['sequence_id']}\n{rec['sequence']}" for rec in self.records_out_exp]
        txt_out = "\n".join(txt_out) + "\n"
        msa.Popen = MockPopenMuscle(
            orig=msa.Popen,
            stdin=txt_in,
            stdout=txt_out)

    def tearDown(self):
        super().tearDown()
        msa.Popen = msa.Popen.orig


class TestMSALive(CommonMSATests, TestBase, TestLive):
    """Basic test of msa with actual MUSCLE calls."""


class TestMSAEmpty(TestMSA):
    """Test of MSA with no input sequences, with Mock MUSCLE call."""

    def test_msa(self):
        with TemporaryDirectory() as tmpdir:
            with self.assertLogs(level="WARNING") as log_cm:
                stdout, stderr = self.redirect_streams(lambda:
                    msa.msa(self.path/"input.fasta", Path(tmpdir)/"output.fasta"))
                self.assertEqual(len(log_cm.output), 1,
                    "a warning should be logged for empty input")
                self.assertEqual("", stdout)
                self.assertEqual("", stderr)
                self.assertTxtsMatch(
                    self.path/"output.fasta",
                    Path(tmpdir)/"output.fasta")

    def test_run_muscle(self):
        def do_muscle(self, recs):
            self.records_out = msa.run_muscle(recs)
        with self.assertLogs(level="WARNING") as log_cm:
            stdout, stderr = self.redirect_streams(lambda: do_muscle(self, self.records_in))
            self.assertEqual(len(log_cm.output), 1,
                "a warning should be logged for empty input")
            self.assertEqual("", stdout)
            self.assertEqual("", stderr)
            self.assertEqual(self.records_out, self.records_out_exp)


class TestMSAEmptyLive(TestMSALive):
    """Test of MSA with no input sequences, with actual MUSCLE calls."""

    def test_msa(self):
        with TemporaryDirectory() as tmpdir:
            with self.assertLogs(level="WARNING") as log_cm:
                stdout, stderr = self.redirect_streams(lambda:
                    msa.msa(self.path/"input.fasta", Path(tmpdir)/"output.fasta"))
                self.assertEqual(len(log_cm.output), 1,
                    "a warning should be logged for empty input")
                self.assertEqual("", stdout)
                self.assertEqual("", stderr)
                self.assertTxtsMatch(
                    self.path/"output.fasta",
                    Path(tmpdir)/"output.fasta")

    def test_run_muscle(self):
        def do_muscle(self, recs):
            self.records_out = msa.run_muscle(recs)
        with self.assertLogs(level="WARNING") as log_cm:
            stdout, stderr = self.redirect_streams(lambda: do_muscle(self, self.records_in))
            self.assertEqual(len(log_cm.output), 1,
                "a warning should be logged for empty input")
            self.assertEqual("", stdout)
            self.assertEqual("", stderr)
            self.assertEqual(self.records_out, self.records_out_exp)
