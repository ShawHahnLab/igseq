"""Tests for igseq.identity."""

import unittest
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
        identity(
            self.path/"input_query.fasta",
            self.tmp/"output.csv",
            self.path/"input_ref.fasta")
        self.assertTxtsMatch(self.path/"output.csv", self.tmp/"output.csv")

    def test_identity_colmap(self):
        """Test identity with FASTA inputs but column names specified."""
        # should warn if colmap given but not using tabular inputs
        with self.assertLogs(level="WARNING"):
            identity(
                self.path/"input_query.fasta",
                self.tmp/"output.csv",
                self.path/"input_ref.fasta",
                colmap={"sequence": "sequence2"})
            self.assertTxtsMatch(self.path/"output.csv", self.tmp/"output.csv")

    def test_identity_aa(self):
        """Test using amino acid sequences."""
        identity(
            self.path/"input_query_aa.fasta",
            self.tmp/"output_aa.csv",
            self.path/"input_ref_aa.fasta")

    def test_identity_stdout(self):
        """Test writing output to stdout."""
        with open(self.path/"output.csv", encoding="ASCII") as f_in:
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
        identity(
            self.path/"input_query.fasta",
            self.tmp/"output.csv",
            self.path/"input_ref.fasta",
            dry_run=True)
        self.assertFalse((self.tmp/"output.csv").exists())

    def test_identity_single(self):
        """Identity with implicit ref via query."""
        # In this case the first record in the query will be used as the
        # reference.
        identity(
            self.path/"input_query.fasta",
            self.tmp/"output.csv")
        self.assertTxtsMatch(self.path/"output_single.csv", self.tmp/"output.csv")
        # should warn if ref format is specified but no ref path given
        with self.assertLogs(level="WARNING"):
            identity(
                self.path/"input_query.fasta",
                self.tmp/"output.csv",
                fmt_in_ref="fa")
        self.assertTxtsMatch(self.path/"output_single.csv", self.tmp/"output.csv")

    def test_identity_csvgz_out(self):
        """Identity with csv.gz output."""
        identity(
            self.path/"input_query.fasta",
            self.tmp/"output.csv.gz",
            self.path/"input_ref.fasta")
        gunzip(self.tmp/"output.csv.gz")
        self.assertTxtsMatch(self.path/"output.csv", self.tmp/"output.csv")


class TestIdentityTabular(TestBase):
    """Test tabular inputs for identity subcommand."""

    def test_identity(self):
        """Basic test of identity with CSV query."""
        # Here a CSV query and CSV ref can give CSV output
        # the defaults are the same as for convert() so sequence_id and
        # sequence columns will be used.
        identity(
            self.path/"input_query.csv",
            self.tmp/"output.csv",
            self.path/"input_ref.csv")
        self.assertTxtsMatch(self.path/"output.csv", self.tmp/"output.csv")

    def test_identity_columns_seq(self):
        """Test using different columns from input (seq)."""
        # Here a CSV query and CSV ref can give CSV output.
        # The defaults are the same as for convert() so sequence_id and
        # sequence columns will be used unless overridden.
        identity(
            self.path/"input_query.csv",
            self.tmp/"output.csv",
            self.path/"input_ref.csv",
            colmap={"sequence": "sequence2"})
        self.assertTxtsMatch(self.path/"output_col2.csv", self.tmp/"output.csv")

    def test_identity_columns_seqid(self):
        """Test using different columns from input (seq ID)."""
        identity(
            self.path/"input_query.csv",
            self.tmp/"output.csv",
            self.path/"input_ref.csv",
            colmap={"sequence_id": "sequence_id2"})
        self.assertTxtsMatch(self.path/"output_col3.csv", self.tmp/"output.csv")

    def test_identity_single(self):
        """Identity with implicit ref via query."""
        # In this case the first record in the query will be used as the
        # reference.
        identity(
            self.path/"input_query.csv",
            self.tmp/"output.csv")
        self.assertTxtsMatch(self.path/"output_single.csv", self.tmp/"output.csv")


class TestScoreIdentity(unittest.TestCase):
    """Basic tests of the identitiy calculation."""

    def test_score_identity(self):
        """Test that pairs of input sequences each produce the appropriate score."""
        # Bio.Align.PairwiseAligner as used in the implementation will match
        # characters (such as those that might represent IUPAC codes) in a
        # simple way, just literally.  So that means N matches only N rather
        # than any nucleotide, and so on.  (Makse sense, considering that's the
        # only sane way it could handle arbitrary sequence types (NT or AA)
        # without also needing a lot of extra metadata tracking.)
        cases = [
            ("ACTG", "ACTG", 1.00), # identical
            ("AC-G", "A-CG", 1.00), # identical (gaps disregarded)
            ("ACTG", "actg", 1.00), # identical (case disregarded)
            ("ACTG",     "", 0.00), # one blank -> 0 by definition
            (    "", "ACTG", 0.00), # other blank -> 0 by definition
            (    "",     "", 0.00), # both blank -> 0 by definition
            ("ACTG", "ACTA", 0.75), # one mismatch
            ("ACTG", "ACTN", 0.75), # one mismatch, IUPAC if NT
            ("ACTN", "ACTN", 1.00), # identical with IUPAC if NT
            ("ACTR", "ACTR", 1.00), # identical with IUPAC if NT
            ("ACDE", "ACDE", 1.00), # or are these AA?  that works too.
            ("ACDE", "ACDP", 0.75), # one mismatch, AA
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
