import subprocess
from unittest import expectedFailure
from gzip import open as gzopen
from io import BytesIO
from pathlib import Path
from tempfile import TemporaryDirectory
from igseq.convert import convert
from igseq.util import IgSeqError
from .util import TestBase

def gunzip(path):
    subprocess.run(["gunzip", path], check=True)

def gzip(path):
    subprocess.run(["gzip", path], check=True)

class TestConvert(TestBase):
    """Basic tests of convert.

    Tests every feasible input/output pair, including with replacement, for:

      * fasta/fasta.gz
      * fastq/fastq.gz
      * csv/csv.gz
      * tsv/tsv.gz
    """

    # From FASTA

    def test_convert_fa_fa(self):
        """Test converting fasta to fasta.

        This should just reformat by unwrapping.
        """
        with TemporaryDirectory() as tmpdir:
            # fasta -> fasta
            convert(self.path/"wrapped.fasta", Path(tmpdir)/"unwrapped.fasta")
            self.assertTxtsMatch(
                self.path/"unwrapped.fasta",
                Path(tmpdir)/"unwrapped.fasta")
            # fasta.gz -> fasta
            gzip(Path(tmpdir)/"unwrapped.fasta")
            convert(Path(tmpdir)/"unwrapped.fasta.gz", Path(tmpdir)/"unwrapped.fasta")
            self.assertTxtsMatch(
                self.path/"unwrapped.fasta",
                Path(tmpdir)/"unwrapped.fasta")
        with TemporaryDirectory() as tmpdir:
            # fasta -> fasta.gz
            convert(self.path/"wrapped.fasta", Path(tmpdir)/"unwrapped.fasta.gz")
            gunzip(Path(tmpdir)/"unwrapped.fasta.gz")
            self.assertTxtsMatch(
                self.path/"unwrapped.fasta",
                Path(tmpdir)/"unwrapped.fasta")
            # fasta.gz -> fasta.gz
            gzip(Path(tmpdir)/"unwrapped.fasta")
            convert(Path(tmpdir)/"unwrapped.fasta.gz", Path(tmpdir)/"unwrapped2.fasta.gz")
            self.assertGzipsMatch(
                Path(tmpdir)/"unwrapped.fasta.gz",
                Path(tmpdir)/"unwrapped2.fasta.gz")
        # to stdout, but we need the format
        with self.assertRaises(IgSeqError):
            convert(self.path/"wrapped.fasta", "-")
        # now, with the format
        stdout, stderr = self.redirect_streams(
            lambda: convert(self.path/"wrapped.fasta", "-", fmt_out="fa"))
        with open(self.path/"unwrapped.fasta") as f_in:
            stdout_exp = f_in.read()
        self.assertEqual(stdout, stdout_exp)
        self.assertEqual(stderr, "")
        # or as gz
        stdout, stderr = self.redirect_streams(
            lambda: convert(self.path/"wrapped.fasta", "-", fmt_out="fagz"))
        with open(self.path/"unwrapped.fasta") as f_in:
            stdout_exp = f_in.read()
        with gzopen(BytesIO(stdout), "rt", "ascii") as f_in:
            stdout_txt = f_in.read()
        self.assertEqual(stdout_txt, stdout_exp)
        self.assertEqual(stderr, "")

    def test_convert_fa_fq(self):
        """Test converting fasta to fastq.

        This should fake the quality scores, with a warning if the dummy score
        isn't explicitly set.
        """
        with TemporaryDirectory() as tmpdir:
            with self.assertLogs(level="WARNING"):
                convert(self.path/"wrapped.fasta", Path(tmpdir)/"unwrapped.fastq")
            self.assertTxtsMatch(
                self.path/"unwrapped.fastq",
                Path(tmpdir)/"unwrapped.fastq")
        with TemporaryDirectory() as tmpdir:
            convert(self.path/"wrapped.fasta", Path(tmpdir)/"unwrapped.fastq", dummyqual="I")
            self.assertTxtsMatch(
                self.path/"unwrapped.fastq",
                Path(tmpdir)/"unwrapped.fastq")
        with TemporaryDirectory() as tmpdir:
            convert(self.path/"wrapped.fasta", Path(tmpdir)/"unwrapped.fastq.gz", dummyqual="I")
            gunzip(Path(tmpdir)/"unwrapped.fastq.gz")
            self.assertTxtsMatch(
                self.path/"unwrapped.fastq",
                Path(tmpdir)/"unwrapped.fastq")

    def test_convert_fa_csv(self):
        """Test converting fasta to csv.

        There should be two columns of output, one for sequence IDs and one for
        sequences.
        """
        with TemporaryDirectory() as tmpdir:
            convert(self.path/"wrapped.fasta", Path(tmpdir)/"unwrapped.csv")
            self.assertTxtsMatch(
                self.path/"unwrapped.csv",
                Path(tmpdir)/"unwrapped.csv")
        with TemporaryDirectory() as tmpdir:
            convert(self.path/"wrapped.fasta", Path(tmpdir)/"unwrapped.csv.gz")
            gunzip(Path(tmpdir)/"unwrapped.csv.gz")
            self.assertTxtsMatch(
                self.path/"unwrapped.csv",
                Path(tmpdir)/"unwrapped.csv")

    def test_convert_fa_tsv(self):
        """Test converting fasta to tsv.

        Like fa_csv, but tab-separated.
        """
        with TemporaryDirectory() as tmpdir:
            convert(self.path/"wrapped.fasta", Path(tmpdir)/"unwrapped.tsv")
            self.assertTxtsMatch(
                self.path/"unwrapped.tsv",
                Path(tmpdir)/"unwrapped.tsv")
        with TemporaryDirectory() as tmpdir:
            convert(self.path/"wrapped.fasta", Path(tmpdir)/"unwrapped.tsv.gz")
            gunzip(Path(tmpdir)/"unwrapped.tsv.gz")
            self.assertTxtsMatch(
                self.path/"unwrapped.tsv",
                Path(tmpdir)/"unwrapped.tsv")

    # From FASTQ

    def test_convert_fq_fa(self):
        """Test converting fastq to fasta."""
        with TemporaryDirectory() as tmpdir:
            convert(self.path/"unwrapped.fastq", Path(tmpdir)/"unwrapped.fasta")
            self.assertTxtsMatch(
                self.path/"unwrapped.fasta",
                Path(tmpdir)/"unwrapped.fasta")
        with TemporaryDirectory() as tmpdir:
            convert(self.path/"unwrapped.fastq", Path(tmpdir)/"unwrapped.fasta.gz")
            gunzip(Path(tmpdir)/"unwrapped.fasta.gz")
            self.assertTxtsMatch(
                self.path/"unwrapped.fasta",
                Path(tmpdir)/"unwrapped.fasta")

    def test_convert_fq_fq(self):
        """Test converting fastq to fastq.

        This shouldn't change anything.
        """
        with TemporaryDirectory() as tmpdir:
            convert(self.path/"unwrapped.fastq", Path(tmpdir)/"unwrapped.fastq")
            self.assertTxtsMatch(
                self.path/"unwrapped.fastq",
                Path(tmpdir)/"unwrapped.fastq")
        with TemporaryDirectory() as tmpdir:
            convert(self.path/"unwrapped.fastq", Path(tmpdir)/"unwrapped.fastq.gz")
            gunzip(Path(tmpdir)/"unwrapped.fastq.gz")
            self.assertTxtsMatch(
                self.path/"unwrapped.fastq",
                Path(tmpdir)/"unwrapped.fastq")

    def test_convert_fq_csv(self):
        """Test converting fastq to csv.

        There should be three columns of output, one for sequence IDs, one for
        sequences, and one for quality scores.
        """
        with TemporaryDirectory() as tmpdir:
            convert(self.path/"unwrapped.fastq", Path(tmpdir)/"unwrapped.csv")
            self.assertTxtsMatch(
                self.path/"unwrapped_quals.csv",
                Path(tmpdir)/"unwrapped.csv")
        with TemporaryDirectory() as tmpdir:
            convert(self.path/"unwrapped.fastq", Path(tmpdir)/"unwrapped.csv.gz")
            gunzip(Path(tmpdir)/"unwrapped.csv.gz")
            self.assertTxtsMatch(
                self.path/"unwrapped_quals.csv",
                Path(tmpdir)/"unwrapped.csv")

    def test_convert_fq_tsv(self):
        """Test converting fastq to tsv.

        Like fq_csv, but tab-separated.
        """
        with TemporaryDirectory() as tmpdir:
            convert(self.path/"unwrapped.fastq", Path(tmpdir)/"unwrapped.tsv")
            self.assertTxtsMatch(
                self.path/"unwrapped_quals.tsv",
                Path(tmpdir)/"unwrapped.tsv")
        with TemporaryDirectory() as tmpdir:
            convert(self.path/"unwrapped.fastq", Path(tmpdir)/"unwrapped.tsv.gz")
            gunzip(Path(tmpdir)/"unwrapped.tsv.gz")
            self.assertTxtsMatch(
                self.path/"unwrapped_quals.tsv",
                Path(tmpdir)/"unwrapped.tsv")

    # From CSV/TSV

    def test_convert_csv_fa(self):
        """Test converting csv to fasta."""
        with TemporaryDirectory() as tmpdir:
            convert(self.path/"unwrapped.csv", Path(tmpdir)/"unwrapped.fasta")
            self.assertTxtsMatch(
                self.path/"unwrapped.fasta",
                Path(tmpdir)/"unwrapped.fasta")

    def test_convert_csv_fq(self):
        """Test converting csv to fastq."""
        with TemporaryDirectory() as tmpdir:
            convert(self.path/"unwrapped_quals.csv", Path(tmpdir)/"unwrapped.fastq")
            self.assertTxtsMatch(
                self.path/"unwrapped.fastq",
                Path(tmpdir)/"unwrapped.fastq")
        with TemporaryDirectory() as tmpdir:
            with self.assertLogs(level="WARNING"):
                convert(self.path/"unwrapped.csv", Path(tmpdir)/"unwrapped.fastq")
            self.assertTxtsMatch(
                self.path/"unwrapped.fastq",
                Path(tmpdir)/"unwrapped.fastq")
        with TemporaryDirectory() as tmpdir:
            convert(self.path/"unwrapped.csv", Path(tmpdir)/"unwrapped.fastq", dummyqual="I")
            self.assertTxtsMatch(
                self.path/"unwrapped.fastq",
                Path(tmpdir)/"unwrapped.fastq")

    def test_convert_csv_csv(self):
        """Test converting csv to csv."""
        with TemporaryDirectory() as tmpdir:
            convert(self.path/"unwrapped_quals.csv", Path(tmpdir)/"unwrapped_quals.csv")
            self.assertTxtsMatch(
                self.path/"unwrapped_quals.csv",
                Path(tmpdir)/"unwrapped_quals.csv")
        with TemporaryDirectory() as tmpdir:
            convert(self.path/"unwrapped.csv", Path(tmpdir)/"unwrapped_quals.csv", dummyqual="I")
            self.assertTxtsMatch(
                self.path/"unwrapped_quals.csv",
                Path(tmpdir)/"unwrapped_quals.csv")

    def test_convert_csv_tsv(self):
        """Test converting csv to tsv."""
        with TemporaryDirectory() as tmpdir:
            convert(self.path/"unwrapped_quals.csv", Path(tmpdir)/"unwrapped_quals.tsv")
            self.assertTxtsMatch(
                self.path/"unwrapped_quals.tsv",
                Path(tmpdir)/"unwrapped_quals.tsv")

    def test_convert_tsv_fa(self):
        """Test converting tsv to fasta."""
        with TemporaryDirectory() as tmpdir:
            convert(self.path/"unwrapped.tsv", Path(tmpdir)/"unwrapped.fasta")
            self.assertTxtsMatch(
                self.path/"unwrapped.fasta",
                Path(tmpdir)/"unwrapped.fasta")

    def test_convert_tsv_fq(self):
        """Test converting tsv to fastq."""
        with TemporaryDirectory() as tmpdir:
            convert(self.path/"unwrapped_quals.tsv", Path(tmpdir)/"unwrapped.fastq")
            self.assertTxtsMatch(
                self.path/"unwrapped.fastq",
                Path(tmpdir)/"unwrapped.fastq")

    def test_convert_tsv_csv(self):
        """Test converting tsv to csv."""
        with TemporaryDirectory() as tmpdir:
            convert(self.path/"unwrapped_quals.tsv", Path(tmpdir)/"unwrapped_quals.csv")
            self.assertTxtsMatch(
                self.path/"unwrapped_quals.csv",
                Path(tmpdir)/"unwrapped_quals.csv")

    def test_convert_tsv_tsv(self):
        """Test converting tsv to tsv."""
        with TemporaryDirectory() as tmpdir:
            convert(self.path/"unwrapped_quals.tsv", Path(tmpdir)/"unwrapped_quals.tsv")
            self.assertTxtsMatch(
                self.path/"unwrapped_quals.tsv",
                Path(tmpdir)/"unwrapped_quals.tsv")


class TestConvertDesc(TestConvert):
    """Test with sequence descriptions (not just ID)"""


class TestConvertPathological(TestBase):
    """Test with strange sequence description lines"""

    @expectedFailure
    def test_convert_fa_fa(self):
        """Test converting fasta to fasta.

        These *should* make output identical to the input, but Biopython subtly
        munges certain strings when it parsed out seq ID vs description.
        """
        with TemporaryDirectory() as tmpdir:
            convert(self.path/"seqs.fasta", Path(tmpdir)/"seqs.fasta")
            self.assertTxtsMatch(
                self.path/"seqs.fasta",
                Path(tmpdir)/"seqs.fasta")


class TestConvertDryRun(TestBase):
    """Test that convert(..., dry_run=True) doesn't touch output files."""

    def test_convert_fa_fa(self):
        """Test converting fasta to fasta with dry_run=True.

        Nothing should actually be written.
        """
        with TemporaryDirectory() as tmpdir:
            convert(self.path/"wrapped.fasta", Path(tmpdir)/"unwrapped.fasta", dry_run=True)
            self.assertTrue(not (Path(tmpdir)/"unwrapped.fasta").exists())


class TestConvertCustomCols(TestBase):
    """Test that convert(..., colmap=...) can customize column names used."""

    def test_convert_fa_csv(self):
        """Test converting fasta to csv.

        There should be two columns of output, one for sequence IDs and one for
        sequences, with custom column names.
        """
        with self.subTest("both custom columns"), TemporaryDirectory() as tmpdir:
            convert(
                self.path/"wrapped.fasta", Path(tmpdir)/"unwrapped.csv",
                colmap = {"sequence_id": "SeqID", "sequence": "Seq"})
            self.assertTxtsMatch(
                self.path/"unwrapped.csv",
                Path(tmpdir)/"unwrapped.csv")
        with self.subTest("one custom column"), TemporaryDirectory() as tmpdir:
            convert(
                self.path/"wrapped.fasta", Path(tmpdir)/"unwrapped.csv",
                colmap = {"sequence_id": "SeqID"})
            self.assertTxtsMatch(
                self.path/"unwrapped.alt1.csv",
                Path(tmpdir)/"unwrapped.csv")
        with self.subTest("no custom columns"), TemporaryDirectory() as tmpdir:
            convert(
                self.path/"wrapped.fasta", Path(tmpdir)/"unwrapped.csv",
                colmap = {})
            self.assertTxtsMatch(
                self.path/"unwrapped.alt2.csv",
                Path(tmpdir)/"unwrapped.csv")

    def test_convert_fa_csv_desc(self):
        """Test converting fasta to csv, with descriptions."""
        colmap = {"sequence_id": "SeqID", "sequence": "Seq", "sequence_description": "SeqDesc"}
        # See #59
        with self.subTest("all descs"), TemporaryDirectory() as tmpdir:
            convert(
                self.path/"wrapped_desc.fasta", Path(tmpdir)/"unwrapped_desc.csv",
                colmap=colmap)
            self.assertTxtsMatch(
                self.path/"unwrapped_desc.csv",
                Path(tmpdir)/"unwrapped_desc.csv")
        # just to double-check the edge case of seq in, tabular out, and no
        # desc on first record, since that's given me trouble (see #53)
        with self.subTest("later desc"), TemporaryDirectory() as tmpdir:
            convert(
                self.path/"wrapped_desc_alt.fasta", Path(tmpdir)/"unwrapped_desc_alt.csv",
                colmap=colmap)
            self.assertTxtsMatch(
                self.path/"unwrapped_desc_alt.csv",
                Path(tmpdir)/"unwrapped_desc_alt.csv")

    def test_convert_fq_csv(self):
        """Test converting fastq to csv.

        There should be three columns of output, one for sequence IDs, one for
        sequences, and one for quality scores.
        """
        with TemporaryDirectory() as tmpdir:
            convert(
                self.path/"unwrapped.fastq", Path(tmpdir)/"unwrapped.csv",
                colmap = {"sequence_id": "SeqID", "sequence": "Seq", "sequence_quality": "SeqQual"})
            self.assertTxtsMatch(
                self.path/"unwrapped_quals.csv",
                Path(tmpdir)/"unwrapped.csv")

    def test_convert_csv_fa(self):
        """Test converting csv to fasta.

        It should understand how to pull from custom-named columns.
        """
        with TemporaryDirectory() as tmpdir:
            convert(
                self.path/"unwrapped.csv", Path(tmpdir)/"unwrapped.fasta",
                colmap = {"sequence_id": "SeqID", "sequence": "Seq"})
            self.assertTxtsMatch(
                self.path/"unwrapped.fasta",
                Path(tmpdir)/"unwrapped.fasta")

    def test_convert_csv_fq(self):
        """Test converting csv to fastq.

        Like the csv_fa case, but it should get quality scores from custom col too.
        """
        with TemporaryDirectory() as tmpdir:
            convert(
                self.path/"unwrapped_quals.csv", Path(tmpdir)/"unwrapped.fastq",
                colmap = {"sequence_id": "SeqID", "sequence": "Seq", "sequence_quality": "SeqQual"})
            self.assertTxtsMatch(
                self.path/"unwrapped.fastq",
                Path(tmpdir)/"unwrapped.fastq")
