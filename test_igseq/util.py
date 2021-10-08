import unittest
import random
from pathlib import Path
import os
import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from igseq.demux import BARCODES_FWD, BARCODES_REV

def random_nt(length):
    return "".join(random.choices(["A", "T", "C", "G"], k = length))

def randomize_n(txt):
    return "".join([random.choice(["A", "T", "C", "G"]) if x == "N" else x for x in txt])

def simulate(i1_path, r1_path, r2_path, num=50000, random_fraction=0.5, barcodes=None):
    with gzip.open(i1_path, "wt") as f_i1, gzip.open(r1_path, "wt") as f_r1, gzip.open(r2_path, "wt") as f_r2:
        for idx in range(1, num+1):
            if random.random() > random_fraction:
                # sample read
                # if barcodes given, select from those.  otherwise select at random
                if barcodes:
                    pair = random.choice(barcodes)
                else:
                    pair = (random.choice(BARCODES_FWD), random.choice(BARCODES_REV))
                pair = list(pair)
                pair[0] = randomize_n(pair[0])
                seq_i1 = Seq(pair[1]).reverse_complement()
                seq_r1 = pair[0] + random_nt(309 - len(pair[0]))
                seq_r2 = random_nt(309)
            else:
                # random read
                seq_i1 = random_nt(8)
                seq_r1 = random_nt(309)
                seq_r2 = random_nt(309)
            qual_i1 = [37 for _ in seq_i1]
            qual_r1 = [37 for _ in seq_r1]
            qual_r2 = [37 for _ in seq_r2]
            def writerec(seq, seqid, qual, hndl):
                SeqIO.write(
                    SeqRecord(
                        Seq(seq),
                        id=seqid,
                        description="",
                        letter_annotations={"phred_quality": qual}),
                    hndl,
                    "fastq")
            writerec(seq_i1, f"read{idx}", qual_i1, f_i1)
            writerec(seq_r1, f"read{idx}", qual_r1, f_r1)
            writerec(seq_r2, f"read{idx}", qual_r2, f_r2)
            #SeqIO.write(SeqRecord(Seq(seq_i1), id=f"read{idx}", letter_annotations={"phred_quality": qual_i1}), f_i1, "fastq")
            #SeqIO.write(SeqRecord(Seq(seq_r1), id=f"read{idx}", letter_annotations={"phred_quality": qual_r1}), f_r1, "fastq")
            #SeqIO.write(SeqRecord(Seq(seq_r2), id=f"read{idx}", letter_annotations={"phred_quality": qual_r2}), f_r2, "fastq")

class TestBase(unittest.TestCase):
    def setUp(self):
        self.path = self.__setup_path()
        self.__startdir = os.getcwd()

    def tearDown(self):
        os.chdir(self.__startdir)

    def __setup_path(self):
        """Path for supporting files for each class."""
        path = self.__class__.__module__.split(".") + [self.__class__.__name__]
        path.insert(1, "data")
        path = Path("/".join(path)).resolve()
        return path

    def assertGzipsMatch(self, path1, path2):
        """Assert that the contents of the two .gz files are identical."""
        with gzip.open(path1, "rt") as f1_in, gzip.open(path2, "rt") as f2_in:
            contents1 = f1_in.read()
            contents2 = f2_in.read()
            if contents1 != contents2:
                raise AssertionError(f"mismatch between {path1} and {path2}")
