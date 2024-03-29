"""Tests for igseq.vdj."""

from pathlib import Path
from igseq import vdj
from igseq.util import IgSeqError, DATA
from .util import TestBase

class TestParseVDJPaths(TestBase):
    """Basic test of parse_vdj_paths."""

    def test_parse_vdj_paths(self):
        """parse_vdj_paths should give a list of dicts with info parsed from the given paths"""
        shared = {"fasta": True, "input": str(self.path/"input"), "type": "dir"}
        attrs_list_exp = [
                {"path": self.path/"input/D.fasta", "segment": "D"},
                {"path": self.path/"input/J.fasta", "segment": "J"},
                {"path": self.path/"input/V.fasta", "segment": "V"}]
        for attrs in attrs_list_exp:
            attrs.update(shared)
        # It should be able to take Path objects, strings, or lists
        ref_path_inputs = [
            self.path / "input",
            str(self.path / "input"),
            [self.path / "input"]]
        for ref_paths in ref_path_inputs:
            attrs_list = vdj.parse_vdj_paths(ref_paths)
            self.assertEqual(attrs_list, attrs_list_exp)

    def test_parse_vdj_paths_duplicates(self):
        """parse_vdj_paths shouldn't repeat paths that come from multiple input names"""
        # Instead it'll just squish both input names into the same entries.
        shared = {"fasta": True, "input": str(self.path/"input"), "type": "dir"}
        attrs_list_exp = [
                {"path": self.path/"input/D.fasta", "segment": "D"},
                {"path": self.path/"input/J.fasta", "segment": "J"},
                {"path": self.path/"input/V.fasta", "segment": "V"}]
        for attrs in attrs_list_exp:
            attrs.update(shared)
        ref_paths = [self.path / "input", self.path/"input/D.fasta"]
        attrs_list_exp[0]["type"] = "file"
        attrs_list_exp[0]["input"] = "; ".join([str(p) for p in ref_paths])
        attrs_list = vdj.parse_vdj_paths(ref_paths)
        self.assertEqual(attrs_list, attrs_list_exp)

    def test_parse_vdj_paths_with_files(self):
        """parse_vdj_paths should work with filenames as inputs"""
        # In this case they're sorted as they're given, since the sorting is
        # by-ref and then by-file.
        mkdict = lambda s: {
            "path": self.path/f"input/{s}.fasta",
            "input": str(self.path/f"input/{s}.fasta"),
            "segment": s,
            "fasta": True,
            "type": "file"}
        attrs_list_exp = [mkdict(s) for s in ["V", "D", "J"]]
        paths = [self.path/f"input/{segment}.fasta" for segment in ["V", "D", "J"]]
        attrs_list = vdj.parse_vdj_paths(paths)
        self.assertEqual(attrs_list, attrs_list_exp)

    def test_parse_vdj_paths_with_ref(self):
        """parse_vdj_paths should work with builtin filename fragments"""
        shared = {"fasta": True, "input": str(self.path/"input"), "type": "dir"}
        attrs_list_exp = [
                {"path": self.path/"input/D.fasta", "segment": "D"},
                {"path": self.path/"input/J.fasta", "segment": "J"},
                {"path": self.path/"input/V.fasta", "segment": "V"}]
        for attrs in attrs_list_exp:
            attrs.update(shared)
        for locus_segment in sorted(["IGHV", "IGHD", "IGHJ", "IGKV", "IGKJ", "IGLV", "IGLJ"]):
            attrs_list_exp.append({
                'fasta': True,
                'type': 'internal',
                'input': 'rhesus/imgt',
                'locus': locus_segment[:3],
                'path': DATA/f"germ/rhesus/imgt/{locus_segment[:3]}/{locus_segment}.fasta",
                'ref': 'imgt',
                'segment': locus_segment[3],
                'species': 'rhesus'})
        attrs_list = vdj.parse_vdj_paths([self.path/"input", "rhesus/imgt"])
        self.assertEqual(attrs_list, attrs_list_exp)


class TestParseVDJId(TestBase):
    """Basic test of parse_vdj_id."""

    def test_parse_vdj_id(self):
        """parse_vdj_id should give dictionary of parsed attributes"""
        # Regular IMGT
        self.assertEqual(
            vdj.parse_vdj_id("IGHV1-1*01"), {
            "seqid":   "IGHV1-1*01",
            "prefix":  "",
            "suffix":  "",
            "allele":  "IGHV1-1*01",
            "gene":    "IGHV1-1",
            "family":  "IGHV1",
            "segment": "IGHV",
            "locus":   "IGH"})
        # Bernat 2021
        self.assertEqual(
            vdj.parse_vdj_id("IGHV1-NL_1*01_S2052"), {
            "seqid":   "IGHV1-NL_1*01_S2052",
            "prefix":  "",
            "suffix":  "",
            "allele":  "IGHV1-NL_1*01_S2052",
            "gene":    "IGHV1-NL_1",
            "family":  "IGHV1",
            "segment": "IGHV",
            "locus":   "IGH"})
        # Ramesh 2017 as shown in SONAR, with ORF prefix and unknown scaffold
        # placement label (-X)
        self.assertEqual(
            vdj.parse_vdj_id("ORF_IGHV3-AHA-X*01"), {
            "seqid":   "ORF_IGHV3-AHA-X*01",
            "prefix":  "ORF_",
            "suffix":  "",
            "allele":  "IGHV3-AHA-X*01",
            "gene":    "IGHV3-AHA-X",
            "family":  "IGHV3",
            "segment": "IGHV",
            "locus":   "IGH"})
        # Zhang 2019 typical syntax
        self.assertEqual(
            vdj.parse_vdj_id("VH3.9B*01c"), {
            "seqid":   "VH3.9B*01c",
            "prefix":  "",
            "suffix":  "",
            "allele":  "VH3.9B*01c",
            "gene":    "VH3.9B",
            "family":  "IGHV3",
            "segment": "IGHV",
            "locus":   "IGH"})
        # Zhang 2019 other syntax
        self.assertEqual(
            vdj.parse_vdj_id("IGHV3-9A*02c"), {
            "seqid":   "IGHV3-9A*02c",
            "prefix":  "",
            "suffix":  "",
            "allele":  "IGHV3-9A*02c",
            "gene":    "IGHV3-9A",
            "family":  "IGHV3",
            "segment": "IGHV",
            "locus":   "IGH"})


class TestParseVDJPathsMissing(TestBase):
    """Test parse_vdj_paths with missing input."""

    def test_parse_vdj_paths(self):
        with self.assertRaises(IgSeqError):
            vdj.parse_vdj_paths(self.path / "missing")


class TestGetInternalVDJ(TestBase):

    def test_get_internal_vdj(self):
        self.assertEqual(
            vdj.get_internal_vdj("rhesus/imgt/IGH"),
            [
                DATA/"germ/rhesus/imgt/IGH/IGHD.fasta",
                DATA/"germ/rhesus/imgt/IGH/IGHJ.fasta",
                DATA/"germ/rhesus/imgt/IGH/IGHV.fasta"])
        self.assertEqual(
            vdj.get_internal_vdj("notfound"),
            [])


class TestParseVDJFilename(TestBase):

    def test_parse_vdj_filename(self):
        equivalents = [
            "IGHV.fasta",
            Path("IGHV.fasta"),
            "IGHV.fa",
            "ighv.FA",
            "IGH/V.fasta"]
        for filename in equivalents:
            attrs = vdj.parse_vdj_filename(filename)
            attrs_exp = {
                "path": Path(filename),
                "locus": "IGH",
                "fasta": True,
                "segment": "V"}
            self.assertEqual(attrs, attrs_exp)


class TestCombineVDJ(TestBase):

    def test_combine_vdj(self):
        # at the simplest this will just take one input FASTA and re-write it
        # to the output.
        input_v = self.path/"input/V.fasta"
        attrs_list = [{
            "path": input_v, "segment": "V", "fasta": True, "type": "file"}]
        fasta = self.tmp/"output/V.fasta"
        vdj.combine_vdj_by_attrs(attrs_list, fasta)
        self.assertTxtsMatch(fasta, self.path/"output/V.fasta")


class TestCombineVDJWithInternal(TestBase):

    def test_combine_vdj(self):
        # a list of per-segment FASTAs should be combined into one, with a
        # suffix appended to each seq ID if needed to differentiate them.
        # I thought this seemed like a pragmatic way to do it but IgBLAST wants
        # sequence IDs shorter than 50 characters, so jamming long filesystem
        # paths into the IDs breaks that.  Maybe autogenerated sequence-based
        # suffixes and a lookup table?
        input_v = self.path/"input/V.fasta"
        attrs_list = [{
            "path": input_v, "segment": "V", "fasta": True, "type": "file"}]
        attrs_list.append({
            'fasta': True,
            'type': 'internal',
            'input': 'rhesus/imgt',
            'locus': "IGH",
            'path': DATA/"germ/rhesus/imgt/IGH/IGHV.fasta",
            'ref': 'imgt',
            'segment': "V",
            'species': 'rhesus'})
        fasta = self.tmp/"output/V.fasta"
        vdj.combine_vdj_by_attrs(attrs_list, fasta)
        self.assertEqual(
            fasta.read_text(),
            (self.path/"output/V.fasta").read_text().replace("INPUTV", str(input_v)))


class TestGroup(TestBase):

    def test_group(self):
        attrs_in = [
            {"path": "path1", "segment": "V"},
            {"path": "path2", "segment": "J"}]
        groups_exp = {
            "V": [{"path": "path1", "segment": "V"}],
            "D": [],
            "J": [{"path": "path2", "segment": "J"}] }
        groups = vdj.group(attrs_in)
        self.assertEqual(groups, groups_exp)

    def test_group_with_keyfunc(self):
        attrs_in = [
            {"path": "path1", "segment": "V", "group": "group1"},
            {"path": "path2", "segment": "J", "group": "group1"},
            {"path": "path3", "segment": "V", "group": "group2"},
            {"path": "path4", "segment": "J", "group": "group2"}]
        groups_exp = {
            "group1": {
                "V": [{"path": "path1", "segment": "V", "group": "group1"}],
                "D": [],
                "J": [{"path": "path2", "segment": "J", "group": "group1"}]},
            "group2": {
                "V": [{"path": "path3", "segment": "V", "group": "group2"}],
                "D": [],
                "J": [{"path": "path4", "segment": "J", "group": "group2"}]}}
        groups = vdj.group(attrs_in, keyfunc=lambda entry: entry["group"])
        self.assertEqual(groups, groups_exp)
