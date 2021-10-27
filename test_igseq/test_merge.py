from igseq.merge import merge
from .util import TestBase

class TestMerge(TestBase):
    """Basic tests of merge."""

    def test_merge_dir_input(self):
        raise self.skipTest("not yet implemented")

    def test_merge_file_input(self):
        raise self.skipTest("not yet implemented")


class TestMergeEmpty(TestBase):
    """Test merge behavior with empty files.

    PEAR crashes with empty inputs, so the merge function should detect this
    and skip actually calling the pear executable.
    """

    def test_merge_dir_input(self):
        raise self.skipTest("not yet implemented")

    def test_merge_file_input(self):
        raise self.skipTest("not yet implemented")
