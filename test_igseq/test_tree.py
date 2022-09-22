from abc import ABC
from pathlib import Path
from tempfile import TemporaryDirectory
from igseq import msa
from igseq import tree
from igseq.record import RecordReader
from .util import TestBase, TestLive, MockPopen
from .test_msa import MockPopenMuscle


class CommonTreeTests(ABC):
    """Shared test logic for actual test case classes below."""

    def setUp(self):
        super().setUp()
        with RecordReader(self.path/"input.fasta") as reader:
            self.records_in = list(reader)
        with RecordReader(self.path/"input.aln.fasta") as reader:
            self.records_in_aln = list(reader)
        with open(self.path/"output.tree") as f_in:
            self.tree_newick_text_exp = f_in.read()
        with open(self.path/"output.nex") as f_in:
            self.tree_nexus_text_exp = f_in.read()
        self.tree_newick_text = None
        self.tree_nexus_text = None

    def test_tree(self):
        """Test tree with unaligned file input and newick output."""
        with TemporaryDirectory() as tmpdir:
            stdout, stderr = self.redirect_streams(lambda:
                tree.tree(self.path/"input.fasta", Path(tmpdir)/"output.tree"))
            self.assertEqual("", stdout)
            self.assertTrue(stderr.startswith("\nmuscle 5"))
            self.assertTxtsMatch(
                self.path/"output.tree",
                Path(tmpdir)/"output.tree")

    def test_tree_aln(self):
        """Test tree with aligned file input and newick output."""
        with TemporaryDirectory() as tmpdir:
            stdout, stderr = self.redirect_streams(lambda:
                tree.tree(self.path/"input.aln.fasta", Path(tmpdir)/"output.tree"))
            self.assertEqual("", stdout)
            self.assertEqual("", stderr)
            self.assertTxtsMatch(
                self.path/"output.tree",
                Path(tmpdir)/"output.tree")

    def test_tree_aln_nexus(self):
        """Test tree with aligned file input and NEXUS output."""
        with TemporaryDirectory() as tmpdir:
            stdout, stderr = self.redirect_streams(lambda:
                tree.tree(self.path/"input.aln.fasta", Path(tmpdir)/"output.nex"))
            self.assertEqual("", stdout)
            self.assertEqual("", stderr)
            self.assertTxtsMatch(
                self.path/"output.nex",
                Path(tmpdir)/"output.nex")

    def test_run_fasttree_aln(self):
        """Test run_fasttree with aligned record inputs and newick output."""
        # This doesn't handle alignments so use already-aligned seqs
        def do_fasttree(self, recs):
            self.tree_newick_text = tree.run_fasttree(recs)
        stdout, stderr = self.redirect_streams(lambda: do_fasttree(self, self.records_in_aln))
        self.assertEqual("", stdout)
        self.assertEqual("", stderr)
        self.assertEqual(self.tree_newick_text, self.tree_newick_text_exp)


class TestTree(CommonTreeTests, TestBase):
    """Basic test of tree, with Mock fasttree and MUSCLE calls."""

    def setUp(self):
        super().setUp()
        txt_in = [f">{rec['sequence_id']}\n{rec['sequence']}" for rec in self.records_in]
        txt_in = "\n".join(txt_in) + "\n"
        txt_aln = [f">{rec['sequence_id']}\n{rec['sequence']}" for rec in self.records_in_aln]
        txt_aln = "\n".join(txt_aln) + "\n"
        msa.Popen = MockPopenMuscle(orig=msa.Popen, stdin=txt_in, stdout=txt_aln)
        tree.Popen = MockPopen(orig=tree.Popen, stdin=txt_aln, stdout=self.tree_newick_text_exp)

    def tearDown(self):
        super().tearDown()
        msa.Popen = msa.Popen.orig
        tree.Popen = tree.Popen.orig


class TestTreeLive(CommonTreeTests, TestBase, TestLive):
    """Basic test of msa with actual MUSCLE calls."""
