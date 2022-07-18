from abc import ABC
from igseq import show
from igseq.util import FILES
from .util import TestBase

def pathlist(paths):
    return "\n".join([str(p) for p in paths]) + "\n"

class TestListFiles(TestBase):
    """Basic tests of list_files"""

    def test_list_files_empty_list(self):
        """an empty list in should list all files"""
        stdout, stderr = self.redirect_streams(
            lambda: show.list_files([]))
        self.assertEqual(stdout, pathlist(FILES))
        self.assertEqual(stderr, "")

    def test_list_files_empty_string(self):
        """an empty string in should list all files"""
        stdout, stderr = self.redirect_streams(
            lambda: show.list_files([""]))
        self.assertEqual(stdout, pathlist(FILES))
        self.assertEqual(stderr, "")

    def test_list_files_one_item(self):
        """one item should list just matching files"""
        stdout, stderr = self.redirect_streams(
            lambda: show.list_files(["primers"]))
        stdout_exp = pathlist([p for p in FILES if "primers" in str(p)])
        self.assertEqual(stdout, stdout_exp)
        self.assertEqual(stderr, "")

    def test_list_files_two_items(self):
        """multiple items should list all files matching all items"""
        stdout, stderr = self.redirect_streams(
            lambda: show.list_files(["rhesus/imgt", "V"]))
        stdout_exp = pathlist([p for p in FILES if "rhesus/imgt" in str(p) and "V" in str(p)])
        self.assertEqual(stdout, stdout_exp)
        self.assertEqual(stderr, "")


class TestShowFiles(TestBase):
    """Basic tests of show_files"""

    def test_show_files_empty_list(self):
        """An empty list in should show nothing"""
        with self.assertLogs(level="WARNING") as log_cm:
            stdout, stderr = self.redirect_streams(
                lambda: show.show_files([]))
            self.assertEqual(len(log_cm.output), 1)
            self.assertIn("No items to show", log_cm.output[0])
        self.assertEqual(stdout, "")
        self.assertEqual(stderr, "")

    def test_show_files_empty_string(self):
        """An empty string in should show nothing"""
        with self.assertLogs(level="WARNING") as log_cm:
            stdout, stderr = self.redirect_streams(
                lambda: show.show_files([""]))
            self.assertEqual(len(log_cm.output), 1)
            self.assertIn("No items to show", log_cm.output[0])
        self.assertEqual(stdout, "")
        self.assertEqual(stderr, "")

    def test_show_files_one_item(self):
        """One item should show all matching files"""
        stdout, stderr = self.redirect_streams(
            lambda: show.show_files(["primers"]))
        stdout_exp = """ Species    Type                               Seq
   human   gamma  GCCAGGGGGAAGACCGATGGGCCCTTGGTGGA
   human   alpha GAGGCTCAGCGGGAAGACCTTGGGGCTGGTCGG
   human      mu AGGAGACGAGGGGGAAAAGGGTTGGGGCGGATG
   human epsilon GCGGGTCAAGGGGAAGACGGATGGGCTCTGTGT
   human   delta CTGATATGATGGGGAACACATCCGGAGCCTTGG
   human   kappa GCGGGAAGATGAAGACAGATGGTGCAGCCACAG
   human  lambda GGCCTTGTTGGCTTGAAGCTCCTCAGAGGAGGG
  rhesus   gamma  GCCAGGGGGAAGACCGATGGGCCCTTGGTGGA
  rhesus   alpha GAGGCTCAGCGGGAAGACCTTGGGGCTGGTCGG
  rhesus      mu GAGACGAGGGGGAAAAGGGTTGGGGCGGATGCA
  rhesus epsilon CGGGTCAAGGGGAAGACGGATGGGCTCTGTGTG
  rhesus   delta CTGATATGATGGGGAACACATCCGGAGCCTTGG
  rhesus   kappa GCGGGAAGATGAAGACAGATGGTGCAGCCACAG
  rhesus  lambda GGCCTTGTTGGCTTGAAGCTCCTCAGAGGAGGG
"""
        self.assertEqual(stdout, stdout_exp)
        self.assertEqual(stderr, "")

    def test_show_files_two_items(self):
        """Multiple items should show all files matching all items"""
        stdout, stderr = self.redirect_streams(
            lambda: show.show_files(["rhesus/imgt", "V"]))
        stdout_lines = stdout.splitlines()
        self.assertEqual(stdout_lines[0], ">IGHV1-111*01")
        self.assertEqual(len(stdout_lines), 914)
        self.assertEqual(stderr, "")

class TestShowOneFile(ABC):
    """Abstract class for below specific tests"""

    def test_show(self):
        for path in (self.path / "input").glob("*"):
            stdout, stderr = self.redirect_streams(
                lambda: show.show_files([path]))
            try:
                with open(self.path/"output/stdout.txt") as f_in:
                    stdout_exp = f_in.read()
            except FileNotFoundError:
                stdout_exp = ""
            self.assertEqual(stdout, stdout_exp)
            self.assertEqual(stderr, "")

class TestShowCSV(TestBase, TestShowOneFile):
    """Test show_files with external CSV file"""

class TestShowStubCSV(TestBase, TestShowOneFile):
    """Test show_files with header-only CSV file"""

class TestShowTree(TestBase, TestShowOneFile):
    """Test show_files with external tree file"""

class TestShowMissingCSV(TestBase):
    """Test show_files with non-existent path"""

    def test_show(self):
        with self.assertLogs(level="WARNING") as log_cm:
            stdout, stderr = self.redirect_streams(
                lambda: show.show_files([self.path / "input/missing.csv"]))
            self.assertEqual(len(log_cm.output), 1)
            self.assertIn("No files found", log_cm.output[0])
        self.assertEqual(stdout, "")
        self.assertEqual(stderr, "")
