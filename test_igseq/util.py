"""Utils for the tests, not tests for igseq.util."""

import sys
import unittest
from unittest.mock import Mock
import shutil
import random
from io import StringIO, BytesIO
from contextlib import redirect_stdout, redirect_stderr
from tempfile import TemporaryDirectory
from pathlib import Path
import os
import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from igseq.demux import BARCODES_FWD, BARCODES_REV

def __boolish(txt):
    if txt.upper() in ["TRUE", "YES", "1"]:
        return True
    if txt.upper() in ["FALSE", "NO", "0", ""]:
        return False
    raise ValueError

# run "live" (real bcl2fastq/cutadapt/pear/etc.) tests?
LIVE = __boolish(os.getenv("TEST_IGSEQ_LIVE", "True"))


class MockPopen(Mock):
    # Mock for Popen used in a particular way:
    # checks for expected text written to its stdin
    # produces canned stderr and stdout, if so
    def __init__(self, *args, **kwargs):
        texts = self.setup_texts(kwargs)
        orig = kwargs["orig"]
        del kwargs["orig"]
        kwargs["spec"] = orig
        # for some reason touching attributes before calling Mock's __init__
        # leads to AttributeError: _mock_methods so we call this first (but
        # only after some changes to kwargs)
        # https://stackoverflow.com/a/62550261/
        super().__init__(*args, **kwargs)
        # NOW we can actually set attributes
        self.orig = orig
        self.texts = texts
        # return itself, to pretend we've called Popen()
        self.return_value = self
        self.__enter__ = self
        self.__exit__ = Mock(return_value=False)
        # returncode is None until a process ends
        self.returncode = None
        self.__exit__.side_effect = self.exit
        self.reset_streams()

    def reset_streams(self):
        self.stdin = StringIO()
        self.stderr = StringIO(self.texts["stderr"])
        self.stdin.close = Mock(side_effect=self.stdin_close)
        self.stdout = StringIO(self.texts["stdout"])

    def setup_texts(self, kwargs):
        texts = {"stdin": "", "stdout": "", "stderr": ""}
        for key in texts.keys():
            try:
                texts[key] = kwargs[key]
                del kwargs[key]
            except KeyError:
                texts[key] = ""
        return texts

    def exit(self, *_, **__):
        sys.stderr.write(self.texts["stderr"])
        # when exiting, Popen() waits for the process so it should now have an
        # exit code
        self.returncode = 0
        # resetting stdin/stdout/stderr means we can call it again, though
        # still only with the same expected use
        self.reset_streams()

    def stdin_close(self, *_, **__):
        # This mock gives output assuming this particular input, so check for
        # that.
        if self.stdin.getvalue() != self.texts["stdin"]:
            raise ValueError("standard input mismatch for mock")


def random_nt(length):
    return "".join(random.choices(["A", "T", "C", "G"], k = length))

def randomize_n(txt):
    return "".join([random.choice(["A", "T", "C", "G"]) if x == "N" else x for x in txt])

def simulate(i1_path, r1_path, r2_path, num=50000, random_fraction=0.5):
    """Randomize I1/R1/R2 reads

    The only thing this does other than shuffle A/C/T/G is select barcodes at
    random for the beginning of R1 and for I1.
    """
    def writerec(seq, seqid, desc, qual, hndl):
        SeqIO.write(
            SeqRecord(
                Seq(seq),
                id=seqid,
                description=desc,
                letter_annotations={"phred_quality": qual}),
            hndl,
            "fastq")
    with gzip.open(i1_path, "wt") as f_i1, \
        gzip.open(r1_path, "wt") as f_r1, \
        gzip.open(r2_path, "wt") as f_r2:
        for idx in range(1, num+1):
            bcfwd = None
            bcrev = None
            if random.random() > random_fraction:
                # sample read
                bcpair = [random.choice(BARCODES_FWD), random.choice(BARCODES_REV)]
                bcfwd = bcpair[0]
                bcrev = bcpair[1]
                bcpair[0] = randomize_n(bcpair[0])
                seq_i1 = Seq(bcpair[1]).reverse_complement()
                seq_r1 = bcpair[0] + random_nt(309 - len(bcpair[0]))
                seq_r2 = random_nt(309)
            else:
                # random read
                seq_i1 = random_nt(8)
                seq_r1 = random_nt(309)
                seq_r2 = random_nt(309)
            qual_i1 = [37 for _ in seq_i1]
            qual_r1 = [37 for _ in seq_r1]
            qual_r2 = [37 for _ in seq_r2]
            writerec(seq_i1, f"read{idx}", f"BarcodeFwd={bcfwd} BarcodeRev={bcrev}", qual_i1, f_i1)
            writerec(seq_r1, f"read{idx}", f"BarcodeFwd={bcfwd} BarcodeRev={bcrev}", qual_r1, f_r1)
            writerec(seq_r2, f"read{idx}", f"BarcodeFwd={bcfwd} BarcodeRev={bcrev}", qual_r2, f_r2)


class TestBase(unittest.TestCase):
    def setUp(self):
        if isinstance(self, TestLive) and not LIVE:
            self.skipTest("skipping tests with external commands")
        self.path = self.__setup_path()
        self.__tmp = TemporaryDirectory()
        self.tmp = Path(self.__tmp.name)
        self.__startdir = os.getcwd()
        self.maxDiff = None

    def tearDown(self):
        # adapted from https://stackoverflow.com/a/39606065
        result = self.defaultTestResult()
        self._feedErrorsToResult(result, self._outcome.errors)
        if result.errors or result.failures:
            shutil.copytree(self.tmp, Path("/tmp/igseq-testdirs")/str(self.tmp).lstrip("/"))
        self.__tmp.cleanup()
        os.chdir(self.__startdir)

    def __setup_path(self):
        """Path for supporting files for each class."""
        path = self.__class__.__module__.split(".") + [self.__class__.__name__]
        path.insert(1, "data")
        path = Path("/".join(path)).resolve()
        return path

    def assertGzipsMatch(self, path1, path2):
        """Assert that the contents of the two .gz files are identical."""
        self.__compare_files(path1, path2, gzip.open)

    def assertTxtsMatch(self, path1, path2):
        """Assert that the contents of the two text files are identical."""
        self.__compare_files(path1, path2, open)

    def __compare_files(self, path1, path2, opener):
        with opener(path1, "rt") as f1_in, opener(path2, "rt") as f2_in:
            contents1 = f1_in.read()
            contents2 = f2_in.read()
            if contents1 != contents2:
                raise AssertionError(f"mismatch between {path1} and {path2}")

    def redirect_streams(self, func):
        """Capture and return stdout and stderr as strings when calling func."""
        # (the underlying buffer is used by convert() to write binary gzipped
        # data to stdout.)
        stdout = StringIO()
        stdout.buffer = BytesIO()
        stderr = StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            func()
        # kludgy way to get whichever "version" of stdout has something in it,
        # preferring the regular one.
        stdout_buffer = stdout.buffer.getvalue()
        stdout = stdout.getvalue()
        if not stdout and stdout_buffer:
            stdout = stdout_buffer
        stderr = stderr.getvalue()
        return stdout, stderr


class TestLive:
    """Parent class for tests that call external commands like igblastn.

    This lets us handle those tests differently if needed, like skipping tests
    that use external commands but running all others.
    """
