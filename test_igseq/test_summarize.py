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

    def test_summarize_csv(self):
        # CSV input should be supported
        stdout, stderr = self.redirect_streams(
            lambda: summarize.summarize(["rhesus/imgt"], self.path/"input/query.csv"))
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

    def test_detect_ref(self):
        """Test that appropriate reference is selected based on species when needed."""
        # (this is all equivalent to the corresponding test for igblast.)
        # Should be able to infer ref from exact species name or synonym
        with self.assertLogs(level="INFO") as log_cm:
            summarize.summarize(None, "-", species="rhesus", dry_run=True)
            self.assertTrue(any("inferred ref path: rhesus" in msg for msg in log_cm.output))
        with self.assertLogs(level="INFO") as log_cm:
            summarize.summarize(None, "-", species="rhesus_monkey", dry_run=True)
            self.assertTrue(any("inferred ref path: rhesus" in msg for msg in log_cm.output))
        # Shouldn't infer ref when explicitly given
        with self.assertLogs(level="INFO") as log_cm:
            summarize.summarize(["rhesus"], "-", species="rhesus", dry_run=True)
            self.assertFalse(any("inferred ref path: rhesus" in msg for msg in log_cm.output))
        # Should still catch unknown names
        with self.assertRaises(IgSeqError) as err_cm:
            summarize.summarize(None, "-", species="unknown", dry_run=True)
            self.assertIn("species not recognized", err_cm.exception.message)
