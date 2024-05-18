from tempfile import TemporaryDirectory
from pathlib import Path
from igseq.cut import cut, parse_regions
from .util import TestBase, TestLive

class TestCutLive(TestBase, TestLive):
    """Basic tests of cut feature (with real igblastn)."""

    def test_cut(self):
        """Basic test of cut function"""
        with TemporaryDirectory() as tmpdir:
            cut(
                "rhesus",
                self.path/"input_query.fasta",
                Path(tmpdir)/"output.fasta",
                "fwr1-fwr3")
            self.assertTxtsMatch(self.path/"output_fwr1-fwr3.fasta", Path(tmpdir)/"output.fasta")
