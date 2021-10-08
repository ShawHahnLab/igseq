import os
from tempfile import TemporaryDirectory
from pathlib import Path
from shutil import copytree
from igseq.phix import phix
from .util import TestBase


class TestPhix(TestBase):
    """Basic test of phix read mapping.

    Here we have one read pair intended to map to PhiX and one not.
    """

    def test_phix(self):
        with TemporaryDirectory() as temp:
            path = Path(temp)
            copytree(self.path/"input", path/"input")
            phix([path/"input/run"])
            # Check that the output filenames all match up
            files = sorted([p.name for p in (path/"phix/run").glob("*")])
            files_expected = sorted(["phix.bam", "phix.counts.csv"])
            self.assertEqual(files_expected, files)
            # Check the read counts.  If this looks right we'll assume the bam
            # does too.
            with open(path/"phix/run/phix.counts.csv", encoding="ASCII") as f_in:
                self.assertEqual(f_in.read(), "Category,Item,NumSeqs\nphix,phix,1\n")
