import re
from abc import ABC
from pathlib import Path
from tempfile import TemporaryDirectory
from itertools import product
import json
from igseq import msa
from igseq import tree
from igseq.util import IgSeqError
from igseq.record import RecordReader
from .util import TestBase, TestLive, MockPopen
from .test_msa import MockPopenMuscle


class CommonTreeTests(ABC):
    """Shared test logic for actual test case classes below."""

    def setUp(self):
        super().setUp()
        with RecordReader(self.path/"seqs.fasta") as reader:
            self.records_in = list(reader)
        with RecordReader(self.path/"seqs.aln.fasta") as reader:
            self.records_in_aln = list(reader)
        with open(self.path/"tree.tree") as f_in:
            self.tree_newick_text_exp = f_in.read()
        with open(self.path/"tree.nex") as f_in:
            self.tree_nexus_text_exp = f_in.read()
        self.tree_newick_text = None
        self.tree_nexus_text = None

    def test_tree_conversions(self):
        """Test all tree input/output combinations"""
        paths_from = ["seqs.fasta", "seqs.aln.fasta", "tree.tree", "tree.nex"]
        paths_to = ["tree.tree", "tree.nex"]
        # all supported combinations
        for path_from, path_to in product(paths_from, paths_to):
            with self.subTest(path_from=path_from, path_to=path_to):
                with TemporaryDirectory() as tmpdir:
                    stdout, stderr = self.redirect_streams(lambda:
                        tree.tree(self.path/path_from, Path(tmpdir)/path_to))
                    self.assertEqual("", stdout)
                    if path_from == "seqs.fasta":
                        # special case, we see muscle run for unaligned sequences
                        self.assertTrue(stderr.startswith("\nmuscle 5"))
                    else:
                        self.assertEqual("", stderr)
                    self.assertTxtsMatch(self.path/path_to, Path(tmpdir)/path_to)

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


class TestTreeEmpty(TestBase):
    """Test with empty inputs.

    These should just throw exceptions.
    """

    def setUp(self):
        super().setUp()
        with RecordReader(self.path/"seqs.fasta") as reader:
            self.records_in = list(reader)

    def test_tree(self):
        """Test tree with empty file input and newick output."""
        with self.assertRaises(IgSeqError):
            with TemporaryDirectory() as tmpdir:
                tree.tree(self.path/"seqs.fasta", Path(tmpdir)/"tree.tree")

    def test_run_fasttree(self):
        """Test run_fasttree with zero records."""
        with self.assertRaises(IgSeqError):
            tree.run_fasttree([])


class TestTreeMulti(TestBase):

    def setUp(self):
        super().setUp()
        self.records_in = []
        with RecordReader(self.path/"seqs.combo.fasta") as reader:
            self.records_in = list(reader)
        with RecordReader(self.path/"seqs.aln.fasta") as reader:
            self.records_in_aln = list(reader)
        with open(self.path/"tree.tree") as f_in:
            self.tree_newick_text_exp = f_in.read()
        with open(self.path/"tree.nex") as f_in:
            self.tree_nexus_text_exp = f_in.read()
        self.tree_newick_text = None
        self.tree_nexus_text = None
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

    def test_tree_multi_input(self):
        # a special case: multiple FASTA inputs will allow implicit set
        # membership based on which files have which sequences
        with TemporaryDirectory() as tmpdir:
            stdout, stderr = self.redirect_streams(lambda:
                tree.tree(
                    [self.path/"seqs.fasta", self.path/"seqs2.fasta"],
                    Path(tmpdir)/"tree.tree"))
            self.assertTxtsMatch(self.path/"tree.tree", Path(tmpdir)/"tree.tree")
        with TemporaryDirectory() as tmpdir:
            stdout, stderr = self.redirect_streams(lambda:
                tree.tree(
                    [self.path/"seqs.fasta", self.path/"seqs2.fasta"],
                    Path(tmpdir)/"tree.nex"))
            self.assertTxtsMatch(self.path/"tree.nex", Path(tmpdir)/"tree.nex")


class TestTreeMultiGaps(TestTreeMulti):
    """Same as TestTreeMulti, but with gaps in the input sequences

    Gaps should always be disregarded with multiple input FASTAs.
    """


class TestTreeDups(TestBase):

    def test_tree_duplicates(self):
        # completely duplicated sequences will be quietly skipped,
        # but it should refuse to allow different sequences with the same
        # sequence ID
        with TemporaryDirectory() as tmpdir:
            with self.assertRaises(IgSeqError):
                tree.tree(self.path/"seqs.fasta", Path(tmpdir)/"tree.tree")


class TestParseLists(TestBase):

    def setUp(self):
        super().setUp()
        self.lists = sorted(self.path.glob("*.txt"))
        with open(self.path/"output.json") as f_in:
            self.output_exp = json.load(f_in)
            for key in self.output_exp:
                self.output_exp[key] = set(self.output_exp[key])

    def test_parse_lists(self):
        """Test that parse_lists loads filename lists into seq sets."""
        output = tree.parse_lists(self.lists)
        self.assertEqual(output, self.output_exp)
        # also try explicitly naming the sets
        lists_named = [f"{idx+1}={path}" for idx, path in enumerate(self.lists)]
        output_named = {re.sub("set", "", k): v for k, v in self.output_exp.items()}
        output = tree.parse_lists(lists_named)
        self.assertEqual(output, output_named)
        # or just one.  it should still use the second's name
        lists_named[0] = re.sub("^1=", "", lists_named[0])
        output_named["set1"] = output_named["1"]
        del output_named["1"]
        output = tree.parse_lists(lists_named)
        self.assertEqual(output, output_named)
        # None should be equivalent to no lists given
        self.assertEqual(tree.parse_lists(None), {})


class TestParseColors(TestBase):

    def test_parse_colors(self):
        # implicit naming
        color_texts = ["#ff0000", "#0000ff"]
        colors_exp = {"set1": [255, 0, 0], "set2": [0, 0, 255]}
        colors = tree.parse_colors(color_texts)
        self.assertEqual(colors, colors_exp)
        # explicit naming
        color_texts = ["1=#ff0000", "2=#0000ff"]
        colors_exp = {"1": [255, 0, 0], "2": [0, 0, 255]}
        colors = tree.parse_colors(color_texts)
        self.assertEqual(colors, colors_exp)
        # just one name
        color_texts = ["1=#ff0000", "#0000ff"]
        colors_exp = {"1": [255, 0, 0], "set2": [0, 0, 255]}
        colors = tree.parse_colors(color_texts)
        self.assertEqual(colors, colors_exp)
        # None should be equivalent to no colors given
        self.assertEqual(tree.parse_colors(None), {})


class TestLooksAligned(TestBase):

    def test_looks_aligned(self):
        self.assertFalse(tree.looks_aligned([
            {"sequence_id": "seq1", "sequence": "ACG"},
            {"sequence_id": "seq2", "sequence": "ACTG"},
            ]))
        self.assertTrue(tree.looks_aligned([
            {"sequence_id": "seq1", "sequence": "ACTG"},
            {"sequence_id": "seq2", "sequence": "ACTG"},
            ]))
        self.assertTrue(tree.looks_aligned([
            {"sequence_id": "seq1", "sequence": "A-TG"},
            {"sequence_id": "seq2", "sequence": "ACTG"},
            ]))
        # one sequence "looks aligned" I guess
        self.assertTrue(tree.looks_aligned([
            {"sequence_id": "seq1", "sequence": "ATG"},
            ]))
        # no sequences doesn't
        self.assertFalse(tree.looks_aligned([]))
