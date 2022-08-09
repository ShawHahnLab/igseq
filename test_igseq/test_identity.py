import unittest
from tempfile import TemporaryDirectory
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from igseq.identity import identity, align, calc_identity, calc_coverage
from .util import TestBase
from .test_convert import gunzip

class TestIdentity(TestBase):
    """Basic tests of identity subcommand."""

    def test_identity(self):
        """Basic test of identity with FASTA inputs."""
        # Here a FASTA query and FASTA ref can give CSV output
        with TemporaryDirectory() as tmpdir:
            identity(
                self.path/"input_query.fasta",
                Path(tmpdir)/"output.csv",
                self.path/"input_ref.fasta")
            self.assertTxtsMatch(self.path/"output.csv", Path(tmpdir)/"output.csv")
            # should warn if colmap given but not using tabular inputs
            with self.assertLogs(level="WARNING"):
                identity(
                    self.path/"input_query.fasta",
                    Path(tmpdir)/"output.csv",
                    self.path/"input_ref.fasta",
                    colmap={"sequence": "sequence2"})
                self.assertTxtsMatch(self.path/"output.csv", Path(tmpdir)/"output.csv")

    def test_identity_aa(self):
        """Test using amino acid sequences."""
        with TemporaryDirectory() as tmpdir:
            identity(
                self.path/"input_query_aa.fasta",
                Path(tmpdir)/"output_aa.csv",
                self.path/"input_ref_aa.fasta")

    def test_identity_stdout(self):
        """Test writing output to stdout."""
        with open(self.path/"output.csv") as f_in:
            stdout_expected = f_in.read()
        stdout, stderr = self.redirect_streams(
            lambda: identity(
                    self.path/"input_query.fasta",
                    "-",
                    self.path/"input_ref.fasta"))
        self.assertEqual(stdout, stdout_expected)
        self.assertEqual(stderr, "")

    def test_identity_dryrun(self):
        """Test that dry run won't write an output file."""
        with TemporaryDirectory() as tmpdir:
            identity(
                self.path/"input_query.fasta",
                Path(tmpdir)/"output.csv",
                self.path/"input_ref.fasta",
                dry_run=True)
            self.assertFalse((Path(tmpdir)/"output.scv").exists())

    def test_identity_single(self):
        """Identity with implicit ref via query."""
        # In this case the first record in the query will be used as the
        # reference.
        with TemporaryDirectory() as tmpdir:
            identity(
                self.path/"input_query.fasta",
                Path(tmpdir)/"output.csv")
            self.assertTxtsMatch(self.path/"output_single.csv", Path(tmpdir)/"output.csv")
            # should warn if ref format is specified but no ref path given
            with self.assertLogs(level="WARNING"):
                identity(
                    self.path/"input_query.fasta",
                    Path(tmpdir)/"output.csv",
                    fmt_in_ref="fa")
            self.assertTxtsMatch(self.path/"output_single.csv", Path(tmpdir)/"output.csv")

    def test_identity_csvgz_out(self):
        """Identity with csv.gz output."""
        with TemporaryDirectory() as tmpdir:
            identity(
                self.path/"input_query.fasta",
                Path(tmpdir)/"output.csv.gz",
                self.path/"input_ref.fasta")
            gunzip(Path(tmpdir)/"output.csv.gz")
            self.assertTxtsMatch(self.path/"output.csv", Path(tmpdir)/"output.csv")


class TestIdentityTabular(TestBase):
    """Test tabular inputs for identity subcommand."""

    def test_identity(self):
        """Basic test of identity with CSV query."""
        # Here a CSV query and CSV ref can give CSV output
        # the defaults are the same as for convert() so sequence_id and
        # sequence columns will be used.
        with TemporaryDirectory() as tmpdir:
            identity(
                self.path/"input_query.csv",
                Path(tmpdir)/"output.csv",
                self.path/"input_ref.csv")
            self.assertTxtsMatch(self.path/"output.csv", Path(tmpdir)/"output.csv")

    def test_identity_columns(self):
        """Test using different columns from input."""
        # Here a CSV query and CSV ref can give CSV output
        # the defaults are the same as for convert() so sequence_id and
        # sequence columns will be used.
        with TemporaryDirectory() as tmpdir:
            identity(
                self.path/"input_query.csv",
                Path(tmpdir)/"output.csv",
                self.path/"input_ref.csv",
                colmap={"sequence": "sequence2"})
            self.assertTxtsMatch(self.path/"output_col2.csv", Path(tmpdir)/"output.csv")

    def test_identity_single(self):
        """Identity with implicit ref via query."""
        # In this case the first record in the query will be used as the
        # reference.
        with TemporaryDirectory() as tmpdir:
            identity(
                self.path/"input_query.csv",
                Path(tmpdir)/"output.csv")
            self.assertTxtsMatch(self.path/"output_single.csv", Path(tmpdir)/"output.csv")


class TestAlign(unittest.TestCase):
    """Basic tests of sequence alignment."""

    def test_align(self):
        """Test that Seq objects can be supplied but not SeqRecords."""
        # Seq objects should work
        seq1 = Seq("ACTG")
        seq2 = Seq("ACTG")
        self.assertEqual(align(seq1, seq2).score, 4)
        # can't apply SeqReords though
        rec1 = SeqRecord(seq1, id="seq1")
        rec2 = SeqRecord(seq2, id="seq2")
        with self.assertRaises(ValueError):
            align(rec1, rec2)


class TestCalc(unittest.TestCase):
    """Basic tests of the identity and coverage calculations."""

    def setUp(self):
        # Bio.Align.PairwiseAligner as used in the implementation will match
        # characters (such as those that might represent IUPAC codes) in a
        # simple way, just literally.  So that means N matches only N rather
        # than any nucleotide, and so on.  (Makse sense, considering that's the
        # only sane way it could handle arbitrary sequence types (NT or AA)
        # without also needing a lot of extra metadata tracking.)
        #
        # Each entry in this list is the target sequence, query sequence,
        # expected identity, and expected coverage.
        self.cases = [
            ("ACTG", "ACTG", 1.00, 1.00), # identical
            ("AC-G", "A-CG", 1.00, 1.00), # identical (gaps disregarded)
            ("ACTG", "actg", 1.00, 1.00), # identical (case disregarded)
            ("ACTG",  "ACT", 0.75, 0.75), # one mismatch, 75% coverage
            ( "ACT", "ACTG", 0.75, 1.00), # one mismatch, 100% coverage
            ("ACTG",     "", 0.00, 0.00), # one blank -> 0 by definition
            (    "", "ACTG", 0.00, 0.00), # other blank -> 0 by definition
            (    "",     "", 0.00, 0.00), # both blank -> 0 by definition
            ("ACTG", "ACTA", 0.75, 1.00), # one mismatch
            ("ACTG", "ACTN", 0.75, 1.00), # one mismatch, IUPAC if NT
            ("ACTN", "ACTN", 1.00, 1.00), # identical with IUPAC if NT
            ("ACTR", "ACTR", 1.00, 1.00), # identical with IUPAC if NT
            ("ACDE", "ACDE", 1.00, 1.00), # or are these AA?  that works too.
            ("ACDE", "ACDP", 0.75, 1.00), # one mismatch, AA
            ]

    def test_calc_identity(self):
        """Test that pairs of input sequences each produce the appropriate identity score."""
        for case in self.cases:
            with self.subTest(pair=case[0:2]):
                self.assertEqual(
                    calc_identity(align(case[0], case[1])),
                    case[2])

    def test_calc_coverage(self):
        """Test that pairs of input sequences each produce the appropriate coverage score.

        Coverage is for the second sequence relative to the first.
        """
        for case in self.cases:
            with self.subTest(pair=case[0:2]):
                self.assertEqual(
                    calc_coverage(align(case[0], case[1])),
                    case[3])
