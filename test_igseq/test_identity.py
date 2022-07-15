import unittest
from tempfile import TemporaryDirectory
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from igseq.identity import identity, score_identity
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

    def test_identity_stdin(self):
        """Test reading input to stdout."""
        # assume FASTA for this case?  Or require format?
        self.skipTest("not yet implemented")

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


class TestScoreIdentity(unittest.TestCase):
    """Basic tests of the identitiy calculation."""

    def test_score_identity(self):
        """Test that pairs of input sequences each produce the appropriate score."""
        # Bio.Align.PairwiseAligner as used in the implementation will match
        # IUPAC codes in a simple way, just literally.  So that means N matches
        # only N rather than any nucleotide, and so on.
        cases = [
            ("ACTG", "ACTG", 1.00), # identical
            ("ACTG",     "", 0.00), # one blank -> 0 by definition
            (    "", "ACTG", 0.00), # other blank -> 0 by definition
            (    "",     "", 0.00), # both blank -> 0 by definition
            ("ACTG", "ACTA", 0.75), # one mismatch
            ("ACTG", "ACTN", 0.75), # one mismatch, IUPAC
            ("ACTN", "ACTN", 1.00), # identical with IUPAC
            ("ACTR", "ACTR", 1.00), # identical with IUPAC
            ]
        for case in cases:
            with self.subTest(pair=case[0:2]):
                self.assertEqual(
                    score_identity(case[0], case[1]),
                    case[2])

    def test_score_identity_objs(self):
        """Test that Seq objects can be supplied but not SeqRecords."""
        # Seq objects should work
        seq1 = Seq("ACTG")
        seq2 = Seq("ACTG")
        self.assertEqual(score_identity(seq1, seq2), 1)
        # can't apply SeqReords though
        rec1 = SeqRecord(seq1, id="seq1")
        rec2 = SeqRecord(seq2, id="seq2")
        with self.assertRaises(ValueError):
            self.assertEqual(score_identity(rec1, rec2), 1)
