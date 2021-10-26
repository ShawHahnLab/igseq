from unittest.mock import Mock
from igseq import igblast
from .util import TestBase, TestLive

class TestIgblast(TestBase):
    """Basic test of igblast."""

    def test_igblast(self):
        self.skipTest("not yet implemented")


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
        self.skipTest("not yet implemented")
