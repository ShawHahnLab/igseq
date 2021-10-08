from tempfile import TemporaryDirectory
from pathlib import Path
from igseq.demux import demux
from .util import TestBase


class TestDemux(TestBase):
    """Basic test of demux.

    Here we have one read each intended for sample1, sample2, and unassigned.

    For the assigned samples, the barcode (including random prefix) should be
    removed from R1, while R2 and I1 should be copied as-is.  The unassigned
    read should be copied as-is for all three (I1/R1/R2).
    """

    def test_demux(self):
        """Test that I1/R1/R2 trio is correctly split into samples."""
        with TemporaryDirectory() as temp:
            demux([self.path/"input/run"], self.path/"samples.csv", dir_out=temp)
            # Check that the output filenames all match up
            files = sorted([p.name for p in Path(temp).glob("*")])
            files_expected = sorted([p.name for p in (self.path/"output").glob("*")])
            self.assertEqual(files_expected, files)
            # Check the contents of each file
            for path in files_expected:
                if path.endswith(".gz"):
                    self.assertGzipsMatch(Path(temp)/path, self.path/"output"/path)

    def test_demux_run_mismatch(self):
        """Test handling of missing sample metadata for this run."""
        with TemporaryDirectory() as temp:
            with self.assertRaises(ValueError):
                demux([self.path/"input/run"], self.path/"samples.csv", run_id="run2", dir_out=temp)
