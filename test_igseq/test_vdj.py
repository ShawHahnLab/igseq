"""Tests for igseq.vdj."""

from pathlib import Path
from igseq import vdj
from igseq.util import IgSeqError, DATA
from .util import TestBase

class TestParseVDJPaths(TestBase):
    """Basic test of parse_vdj_paths."""

    def test_parse_vdj_paths(self):
        # It should give a list of dictionaries with info parsed from the given
        # paths.
        shared = {"fasta": True, "input": str(self.path/"input"), "type": "dir"}
        attrs_list_exp = [
                {"path": self.path/"input/V.fasta", "segment": "V"},
                {"path": self.path/"input/D.fasta", "segment": "D"},
                {"path": self.path/"input/J.fasta", "segment": "J"}]
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

    def test_parse_vdj_paths_with_ref(self):
        # If we also ask for a fragment of the filenames of builtin FASTA
        # files, it should find those too
        shared = {"fasta": True, "input": str(self.path/"input"), "type": "dir"}
        attrs_list_exp = [
                {"path": self.path/"input/V.fasta", "segment": "V"},
                {"path": self.path/"input/D.fasta", "segment": "D"},
                {"path": self.path/"input/J.fasta", "segment": "J"}]
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
