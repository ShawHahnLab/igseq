from tempfile import TemporaryDirectory
from pathlib import Path
from igseq.trim import trim, get_adapters_fwd, get_adapter_rev
from igseq.util import IgSeqError
from .util import TestBase, TestLive

class TestTrim(TestBase):
    """Basic tests of trim-related functions"""

    def setUp(self):
        super().setUp()
        self.sample = {
            "Sample": "sample1",
            "BarcodeFwdSeq": "NNNNAACCACTA",
            "BarcodeRevSeq": "TAGTGGTT",
            "Type": "gamma"}

    def test_get_adapters_fwd(self):
        """Test getting the adapter sequences to trim from the end of R1"""
        with self.subTest(case="basic"):
            adapter_fwd = get_adapters_fwd(self.sample, "rhesus")
            self.assertEqual(
                adapter_fwd,
                {"rhesus_gamma": "TCCACCAAGGGCCCATCGGTCTTCCCCCTGGC"})
        with self.subTest(case="bad species"):
            with self.assertRaises(IgSeqError):
                adapter_fwd = get_adapters_fwd(self.sample, "rhino")
        with self.subTest(case="bad chain type"):
            samp = self.sample.copy()
            samp["Type"] = "heavy"
            with self.assertRaises(IgSeqError):
                adapter_fwd = get_adapters_fwd(samp)
        with self.subTest(case="no species"):
            adapter_fwd = get_adapters_fwd(self.sample)
            self.assertEqual(
                adapter_fwd, {
                "rhesus_gamma": "TCCACCAAGGGCCCATCGGTCTTCCCCCTGGC",
                "human_gamma": "TCCACCAAGGGCCCATCGGTCTTCCCCCTGGC"})
        with self.subTest(case="no chain type"):
            adapter_fwd = get_adapters_fwd({}, "rhesus")
            self.assertEqual(
                adapter_fwd, {
                "rhesus_gamma": "TCCACCAAGGGCCCATCGGTCTTCCCCCTGGC",
                "rhesus_alpha": "CCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTC",
                "rhesus_mu": "TGCATCCGCCCCAACCCTTTTCCCCCTCGTCTC",
                "rhesus_epsilon": "CACACAGAGCCCATCCGTCTTCCCCTTGACCCG",
                "rhesus_delta": "CCAAGGCTCCGGATGTGTTCCCCATCATATCAG",
                "rhesus_kappa": "CTGTGGCTGCACCATCTGTCTTCATCTTCCCGC",
                "rhesus_lambda": "CCCTCCTCTGAGGAGCTTCAAGCCAACAAGGCC"})
        with self.subTest(case="blank chain type"):
            adapter_fwd = get_adapters_fwd({"Type": ""}, "rhesus")
            self.assertEqual(
                adapter_fwd, {
                "rhesus_gamma": "TCCACCAAGGGCCCATCGGTCTTCCCCCTGGC",
                "rhesus_alpha": "CCGACCAGCCCCAAGGTCTTCCCGCTGAGCCTC",
                "rhesus_mu": "TGCATCCGCCCCAACCCTTTTCCCCCTCGTCTC",
                "rhesus_epsilon": "CACACAGAGCCCATCCGTCTTCCCCTTGACCCG",
                "rhesus_delta": "CCAAGGCTCCGGATGTGTTCCCCATCATATCAG",
                "rhesus_kappa": "CTGTGGCTGCACCATCTGTCTTCATCTTCCCGC",
                "rhesus_lambda": "CCCTCCTCTGAGGAGCTTCAAGCCAACAAGGCC"})

    def test_get_adapter_rev(self):
        """Test getting the adapter sequence to trim from the end of R2"""
        adapter_rev = get_adapter_rev(self.sample)
        self.assertEqual(adapter_rev, "TAGTGGTTNNNNAGATCGGAAGAGCGTCGTGTAGGGAAAGA")

    def test_trim(self):
        """Test the high-level trim() function"""
        self.skipTest("not yet implemented")

    def test_trim_pair(self):
        """Test the lower-level single-sample trim_pair() function"""
        self.skipTest("not yet implemented")

class TestTrimLive(TestBase, TestLive):
    """Basic tests of trim with actual cutadapt.

    Here we have a simple case with two sample with perfect adapters in R1 and
    R2.  We should see the adapters get removed in the output, cutadapt's JSON
    report written, and our counts.csv files written.  This should work with
    either a directory or individual R1/R2 pairs as input.

    Each input sequence has 200 random nucleotides and then the adapter that
    should be trimmed, so the output should be a 200 nt sequence each time.

    samples 1 and 2 are a basic test with different barcode pairs.  Sample 3
    has a valid barcode pair but is missing the expected 5' RACE anchor so the
    read gets filtered out.  sample 4 has the anchor but not the barcode pair
    but those aren't required so it's kept.
    """

    def test_trim_dir_input(self):
        """Test that adapters are trimmed from R1 and R2 pairs with dir input."""
        with TemporaryDirectory() as temp:
            trim([self.path/"input/run"], self.path/"samples.csv", dir_out=temp)
            files = sorted([p.name for p in Path(temp).glob("*")])
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
            for idx in [1, 2, 3, 4]:
                trim([
                    self.path/f"input/run/sample{idx}.R1.fastq.gz",
                    self.path/f"input/run/sample{idx}.R2.fastq.gz"],
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

    def test_trim_custom_args(self):
        """Test giving custom cutadapt arguments"""
        with TemporaryDirectory() as temp:
            temp = Path(temp)
            trim([
                self.path/"input/run/sample1.R1.fastq.gz",
                self.path/"input/run/sample1.R2.fastq.gz"],
                self.path/"samples.csv", dir_out=temp,
                extra_cutadapt_args=["--too-short-output", temp/"short.fastq.gz"])
            # there actually should be no too-short sequences, but, the file
            # should still be created (even though empty)
            self.assertEmpty(temp/"short.fastq.gz")

class TestTrimLiveNoChainType(TestBase, TestLive):
    """Test trimming without specifying chain type (for the constant region primer).

    In this case cutadapt should fall back on using all applicable primer
    sequences (per-species, if specified, or across species, if not).  The test
    files are the same as TestTrimLive, where sample 4 R1 happens to end with
    CCG which, here, ends up matching the alpha primer and getting slightly
    trimmed.
    """

    def test_trim_dir_input(self):
        """Test that adapters are trimmed from R1 and R2 pairs with dir input."""
        with TemporaryDirectory() as temp:
            trim([self.path/"input/run"], self.path/"samples.csv", dir_out=temp)
            files = sorted([p.name for p in Path(temp).glob("*")])
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
            for idx in [1, 2, 3, 4]:
                trim([
                    self.path/f"input/run/sample{idx}.R1.fastq.gz",
                    self.path/f"input/run/sample{idx}.R2.fastq.gz"],
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
