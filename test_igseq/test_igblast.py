import re
from igseq import igblast
from igseq.util import IgSeqError
from .util import TestBase, TestLive

class TestIgblast(TestBase):
    """Basic test of igblast."""

    def test_igblast(self):
        self.skipTest("not yet implemented")

    def test_detect_organism(self):
        """Test that detect_organism can infer correct IgBLAST organism name."""
        # first argument is set of possible inferred species names, second is
        # explicitly given species name, if any.  The matching is
        # case-insensitive and recognizes some synonyms.
        self.assertEqual(igblast.detect_organism(["human"]), "human")
        self.assertEqual(igblast.detect_organism(["homo sapiens"]), "human")
        self.assertEqual(igblast.detect_organism(["rhesus"]), "rhesus_monkey")
        self.assertEqual(igblast.detect_organism(["RHESUS"]), "rhesus_monkey")
        self.assertEqual(igblast.detect_organism(["rhesus-monkey"]), "rhesus_monkey")
        # There has to be at least one, or the second argument must be given
        with self.assertRaises(IgSeqError) as err_cm:
            igblast.detect_organism(set())
        self.assertIn("species not detected from input", err_cm.exception.message)
        self.assertEqual(
            igblast.detect_organism(set(), "rhesus"),
            "rhesus_monkey")
        # Multiples in the first argument also require the second to be given
        with self.assertRaises(IgSeqError) as err_cm:
            igblast.detect_organism(["species1", "species2"])
        self.assertIn("multiple species detected from input", err_cm.exception.message)
        self.assertEqual(
            igblast.detect_organism(["species1", "species2"], "rhesus"),
            "rhesus_monkey")


class TestIgblastInternal(TestBase):
    """Test igblast with internal db."""

    def test_igblast(self):
        self.skipTest("not yet implemented")


class TestIgblastMissing(TestBase):
    """Test igblast with missing db."""

    def test_igblast(self):
        self.skipTest("not yet implemented")


class TestIgblastLive(TestBase, TestLive):
    """Basic test of igblast with actual makeblastdb+igblastn."""

    def test_igblast(self):
        stdout, stderr = self.redirect_streams(
            lambda: igblast.igblast(["rhesus/imgt"], self.path/"input/query.fasta"))
        with open(self.path/"output/stdout.txt", encoding="ascii") as f_in:
            stdout_exp = f_in.read()
        btwn = lambda txt, idx1, idx2: "".join(txt.splitlines(True)[idx1:(idx1+idx2)])
        self.assertEqual(
                btwn(stdout, 2, 97),
                btwn(stdout_exp, 2, 97))
        self.assertEqual(stderr, "")


class TestIgblastLiveCrash(TestBase, TestLive):
    """Test that an igblastn crash is handled"""

    def test_igblast(self):
        """Test that an igblastn crash is caught and re-raised as an IgSeqError.

        Standard error and output from the igblastn process should still show
        up on stderr and stdout as usual.
        """
        def catch_expected_error():
            with self.assertRaises(IgSeqError):
                igblast.igblast(
                    ["rhesus/imgt"], self.path/"input/query.fasta", extra_args=["-bad-arg"])
        stdout, stderr = self.redirect_streams(catch_expected_error)
        self.assertEqual(stdout, "")
        self.assertIn('Error: Unknown argument: "bad-arg"', stderr)
