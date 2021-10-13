from tempfile import TemporaryDirectory
from pathlib import Path
from igseq.trim import trim
from .util import TestBase

class TestTrim(TestBase):
    """Basic tests of trim.

    Here we have a simple case with two sample with perfect adapters in R1 and
    R2.  We should see the adapters get removed in the output, cutadapt's JSON
    report written, and our counts.csv files written.  This should work with
    either a directory or individual R1/R2 pairs as input.

    Each input sequence has 200 random nucleotides and then the adapter that
    should be trimmed, so the output should be a 200 nt sequence each time.
    """

    def test_trim_dir_input(self):
        """Test that adapters are trimmed from R1 and R2 pairs with dir input."""
        with TemporaryDirectory() as temp:
            trim([self.path/"input/run"], self.path/"samples.csv", dir_out=temp)
            files = sorted([p.name for p in Path(temp).glob("*")])
            #import shutil
            #for path in files:
            #    shutil.copy(Path(temp)/path, ".")
            files_expected = sorted([p.name for p in (self.path/"output").glob("*")])
            self.assertEqual(files_expected, files)
            for path in files_expected:
                if path.endswith(".gz"):
                    self.assertGzipsMatch(Path(temp)/path, self.path/"output"/path)
                if path.endswith(".counts.csv"):
                    self.assertTxtsMatch(Path(temp)/path, self.path/"output"/path)

    def test_trim_file_input(self):
        """Test that adapters are trimmed from R1 and R2 pairs with file input."""
        with TemporaryDirectory() as temp:
            trim([
                self.path/"input/run/sample1.R1.fastq.gz",
                self.path/"input/run/sample1.R2.fastq.gz"],
                self.path/"samples.csv", dir_out=temp)
            trim([
                self.path/"input/run/sample2.R1.fastq.gz",
                self.path/"input/run/sample2.R2.fastq.gz"],
                self.path/"samples.csv", dir_out=temp)
            files = sorted([p.name for p in Path(temp).glob("*")])
            files_expected = sorted([p.name for p in (self.path/"output").glob("*")])
            self.assertEqual(files_expected, files)
            for path in files_expected:
                if path.endswith(".gz"):
                    self.assertGzipsMatch(Path(temp)/path,
                            self.path/"output"/path)
                if path.endswith(".counts.csv"):
                    self.assertTxtsMatch(Path(temp)/path, self.path/"output"/path)
