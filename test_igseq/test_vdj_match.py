import re
from pathlib import Path
from tempfile import TemporaryDirectory
from igseq import vdj_match
from igseq.util import IgSeqError
from .util import TestBase, TestLive

class TestVDJMatchLive(TestBase, TestLive):

    def test_vdj_match(self):
        # by default, will print output table
        stdout, stderr = self.redirect_streams(
            lambda: vdj_match.vdj_match(["rhesus/imgt"], self.path/"input/query.fasta"))
        self.assertEqual(stderr, "")
        with open(self.path/"output/stdout.txt") as f_in:
            stdout_exp = f_in.read()
        self.assertEqual(stdout_exp, stdout)

    def test_vdj_match_multi(self):
        # with all matching rhesus refs selected, we should get a block of
        # output rows for each ref found.  but, the zhang2019 ref only has V
        # and J so it's not usable as a standalone IgBLAST DB.  We'll get a
        # warning about skipping that one, and output rows for imgt/
        # sonarramesh/bernat2021.
        stdout, stderr = self.redirect_streams(
            lambda: vdj_match.vdj_match(["rhesus"], self.path/"input/query.fasta"))
        self.assertEqual(stderr, 'No FASTA for D from rhesus/zhang2019; skipping\n')
        with open(self.path/"output/stdout_multi.txt") as f_in:
            stdout_exp = f_in.read()
        self.assertEqual(stdout_exp, stdout)

    def test_vdj_match_missing(self):
        # If we explicitly request a reference that's missing segments, though,
        # that warning becomes an error.
        def wrapper():
            with self.assertRaisesRegex(IgSeqError, "Missing VDJ input for database"):
                vdj_match.vdj_match(["rhesus/zhang2019"], self.path/"input/query.fasta")
        stdout, stderr = self.redirect_streams(wrapper)
        self.assertEqual(stdout, "")
        self.assertEqual(stderr, 'No FASTA for D from rhesus/zhang2019\n')

    def test_vdj_match_csv(self):
        # CSV input should be supported
        stdout, stderr = self.redirect_streams(
            lambda: vdj_match.vdj_match(["rhesus/imgt"], self.path/"input/query.csv"))
        self.assertEqual(stderr, "")
        with open(self.path/"output/stdout.txt") as f_in:
            stdout_exp = f_in.read()
        self.assertEqual(stdout_exp, stdout)

    def test_vdj_match_outfile(self):
        # with an output filename, it will write CSV there, and won't print the
        # table
        with TemporaryDirectory() as temp:
            output_csv = Path(temp)/"output.csv"
            def func():
                vdj_match.vdj_match(
                    ["rhesus/imgt"], self.path/"input/query.fasta", output_csv)
            stdout, stderr = self.redirect_streams(func)
            self.assertEqual(stdout, "")
            self.assertEqual(stderr, "")
            self.assertTxtsMatch(output_csv, self.path/"output/output.csv")

    def test_vdj_match_outfile_show(self):
        # we can optionally both print and save the table
        with TemporaryDirectory() as temp:
            output_csv = Path(temp)/"output.csv"
            def func():
                vdj_match.vdj_match(
                    ["rhesus/imgt"], self.path/"input/query.fasta", output_csv, showtxt=True)
            stdout, stderr = self.redirect_streams(func)
            self.assertEqual(stderr, "")
            with open(self.path/"output/stdout.txt") as f_in:
                stdout_exp = f_in.read()
            self.assertEqual(stdout_exp, stdout)
            self.assertTxtsMatch(output_csv, self.path/"output/output.csv")
