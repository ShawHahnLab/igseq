import re
from pathlib import Path
from tempfile import TemporaryDirectory
from igseq import summarize
from igseq.util import IgSeqError
from .util import TestBase, TestLive

class TestSummarizeLive(TestBase, TestLive):

    def test_summarize(self):
        # by default, will print output table
        stdout, stderr = self.redirect_streams(
            lambda: summarize.summarize(["rhesus/imgt"], self.path/"input/query.fasta"))
        self.assertEqual(stderr, "")
        with open(self.path/"output/stdout.txt") as f_in:
            stdout_exp = f_in.read()
        self.assertEqual(stdout_exp, stdout)

    def test_summarize_outfile(self):
        # with an output filename, it will write CSV there, and won't print the
        # table
        with TemporaryDirectory() as temp:
            output_csv = Path(temp)/"output.csv"
            def func():
                summarize.summarize(
                    ["rhesus/imgt"], self.path/"input/query.fasta", output_csv)
            stdout, stderr = self.redirect_streams(func)
            self.assertEqual(stdout, "")
            self.assertEqual(stderr, "")
            self.assertTxtsMatch(output_csv, self.path/"output/output.csv")

    def test_summarize_multi(self):
        # more than one germline reference included
        stdout, stderr = self.redirect_streams(
            lambda: summarize.summarize(["rhesus"], self.path/"input/query.fasta"))
        self.assertEqual(stderr, "")
        with open(self.path/"output/stdout_ref_rhesus.txt") as f_in:
            stdout_exp = f_in.read()
        self.assertEqual(stdout_exp, stdout)
