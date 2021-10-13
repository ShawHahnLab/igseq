from pathlib import Path
from igseq import util as igseq_util
from .util import TestBase

class TestCommonParent(TestBase):

    def test_common_parent(self):
        parent = igseq_util.common_parent([
            "path/to/something/dir1/subdir1",
            "path/to/something/dir2/subdir2"])
        self.assertEqual(parent, Path("path/to/something"))

    def test_common_parent_paths(self):
        parent = igseq_util.common_parent([
            Path("path/to/something/dir1/subdir1"),
            Path("path/to/something/dir2/subdir2")])
        self.assertEqual(parent, Path("path/to/something"))

    def test_common_parent_dict(self):
        parent = igseq_util.common_parent({
            "R1": Path("path/to/something/dir1/R1.fastq.gz"),
            "R2": Path("path/to/something/dir2/R2.fastq.gz")})
        self.assertEqual(parent, Path("path/to/something"))

    def test_common_parent_none(self):
        parent = igseq_util.common_parent([
            "dir1/subdir1",
            "dir2/subdir2"])
        self.assertEqual(parent, Path("."))

    def test_common_parent_none_absolute(self):
        parent = igseq_util.common_parent([
            "/dir1/subdir1",
            "/dir2/subdir2"])
        self.assertEqual(parent, Path("/"))
