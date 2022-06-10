"""
Helper classes for records stored in FASTA/FASTQ/CSV/TSV with/without gzip.
"""

import sys
import csv
import gzip
import logging
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from igseq import util

LOGGER = logging.getLogger(__name__)

DEFAULT_DUMMY_QUAL = "I"
DEFAULT_COLUMNS = {k: k for k in ["sequence_id", "sequence", "sequence_quality", "sequence_description"]}

# Mapping from file extensions to file format names used here
FMT_EXT_MAP = {
    ".csv": "csv",
    ".tsv": "tsv",
    ".tab": "tsv",
    ".fa": "fa",
    ".fasta": "fa",
    ".afa": "fa",
    ".fq": "fq",
    ".fastq": "fq"}

# Mapping from file format names to functions that take file handles and give
# format-specific record iterators.
READERS = {
    "csv": lambda hndl: csv.DictReader(hndl),
    "tsv": lambda hndl: csv.DictReader(hndl, delimiter="\t"),
    "fa": lambda hndl: SeqIO.parse(hndl, "fasta"),
    "fq": lambda hndl: SeqIO.parse(hndl, "fastq")}

class RecordHandler:
    """Abstract class for shared generic record reading and writing
    functionality.  This behaves as a context manager to handle file opening
    and closing and provides converters between SeqRecord objects and
    dictionaries, Phred quality score integers and Illumina text encoded
    quality scores.

    See RecordReader/RecordWriter for the main stuff.
    """

    def __init__(self, pathlike, fmt=None, colmap=None, dummyqual=None, dry_run=False):
        self.pathlike = pathlike
        self.fmt = self.infer_fmt(fmt)
        self.colmap = {}
        self.colmap.update(DEFAULT_COLUMNS)
        if colmap is not None:
            self.colmap.update(colmap)
        LOGGER.info("colmap: %s", self.colmap)
        self.dummyqual = dummyqual
        self.dry_run = dry_run
        self.handle = None
        self.autoclose = True

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.close()

    def open(self):
        """Open the path as needed for reading or writing."""
        raise NotImplementedError(self)

    def close(self):
        """Close the file handle used here, if it's not stdin/stdout."""
        std_streams = [sys.stdin, sys.stdin.buffer, sys.stdout, sys.stdout.buffer]
        # Flush to handle the case of an already-opened handle that won't be
        # closed here
        if self.handle and not self.handle.closed:
            self.handle.flush()
            if self.handle not in std_streams and self.autoclose:
                self.handle.close()

    def infer_fmt(self, fmt=None):
        """Guess file format from our path, with optional default."""
        fmt_inferred = self._infer_fmt(self.pathlike)
        LOGGER.info("given format: %s", fmt)
        LOGGER.info("inferred format: %s", fmt_inferred)
        if not fmt:
            if not fmt_inferred:
                raise util.IgSeqError(
                    f"format not detected from filename ({self.pathlike}).  "
                    "specify a format manually.")
            fmt = fmt_inferred
        return fmt

    @staticmethod
    def _infer_fmt(path):
        try:
            path = Path(path)
        except TypeError:
            path = Path(str(path.name))
        ext = path.suffix.lower()
        ext2 = Path(path.stem).suffix.lower()
        fmt2 = FMT_EXT_MAP.get(ext2)
        if ext == ".gz" and fmt2:
            return fmt2 + "gz"
        else:
            return FMT_EXT_MAP.get(ext)

    def encode_record(self, record):
        """Convert record dictionary into a SeqRecord object."""
        seq = record[self.colmap["sequence"]]
        seqid = record[self.colmap["sequence_id"]]
        desc = record.get(self.colmap["sequence_description"], "")
        seqrecord = SeqRecord(Seq(seq), id=seqid, description=desc)
        quals = record.get(self.colmap["sequence_quality"])
        if quals:
            nums = self.decode_phred(quals)
            seqrecord.letter_annotations["phred_quality"] = nums
        elif self.dummyqual:
            nums = self.decode_phred(self.dummyqual * len(seqrecord))
            seqrecord.letter_annotations["phred_quality"] = nums
        return seqrecord

    def decode_record(self, obj):
        """Convert object (SeqRecord or dictionary) to dictionary."""
        if isinstance(obj, SeqRecord):
            record = {
                self.colmap["sequence_id"]: obj.id,
                self.colmap["sequence"]: str(obj.seq)}
            quals = obj.letter_annotations.get("phred_quality")
            if quals:
                record[self.colmap["sequence_quality"]] = self.encode_phred(quals)
        else:
            record = obj
        return record

    @staticmethod
    def encode_phred(qualnums):
        """Encode quality scores from list of integers as Illumina text string."""
        return "".join([chr(qual + 33) for qual in qualnums])

    @staticmethod
    def decode_phred(qualtxt):
        """Decode quality scores from Illumina text string to list of integers."""
        return [ord(letter) - 33 for letter in qualtxt]


class RecordReader(RecordHandler):
    """Generic record reader class.

    This will work as an iterator to produce record dictionaries for each
    record from the file (FASTA/FASTQ sequence or CSV/TSV row).
    """

    def __init__(self, pathlike, fmt=None, colmap=None, dry_run=False):
        super().__init__(pathlike, fmt, colmap=colmap, dry_run=dry_run)
        self.reader = None

    def open(self):
        if self.pathlike == "-":
            if "gz" in self.fmt:
                LOGGER.info("reading gzip from stdin")
                self.handle = gzip.open(sys.stdin.buffer, "rt", encoding="ascii")
            else:
                LOGGER.info("reading text from stdin")
                self.handle = sys.stdin
        else:
            if "gz" in self.fmt:
                LOGGER.info("reading from gzip")
                self.handle = gzip.open(self.pathlike, "rt", encoding="ascii")
            else:
                LOGGER.info("reading from text")
                self.handle = open(self.pathlike, "rt", encoding="ascii")
        self.reader = READERS[self.fmt.removesuffix("gz")](self.handle)

    def  __iter__(self):
        return self

    def __next__(self):
        obj = next(self.reader)
        record = self.decode_record(obj)
        return record


class RecordWriter(RecordHandler):
    """Generic record writer class.

    This will work as a context manager for opening and closing the underlying
    file and provides a write() method to write record dictionaries (as
    sequences to FASTA/FASTQ or rows to CSV/TSV).
    """

    def __init__(self, pathlike, fmt=None, colmap=None, dummyqual=None, dry_run=False):
        super().__init__(pathlike, fmt, colmap, dummyqual, dry_run)
        self.writer = None

    def open(self):
        if hasattr(self.pathlike, "fileno"):
            self.handle = self.pathlike
            self.autoclose = False
        elif self.pathlike == "-":
            if "gz" in self.fmt:
                LOGGER.info("writing gzip to stdout")
                if not self.dry_run:
                    self.handle = gzip.open(sys.stdout.buffer, "wt", encoding="ascii")
            else:
                LOGGER.info("writing text to stdout")
                if not self.dry_run:
                    self.handle = sys.stdout
        else:
            if "gz" in self.fmt:
                LOGGER.info("writing to gzip")
                if not self.dry_run:
                    self.handle = gzip.open(self.pathlike, "wt", encoding="ascii")
            else:
                LOGGER.info("writing to text")
                if not self.dry_run:
                    self.handle = open(self.pathlike, "wt", encoding="ascii")

    def write(self, record):
        """Write a record dictionary to the output stream."""
        if self.dummyqual and not self.colmap["sequence_quality"] in record:
            quals = self.dummyqual * len(record[self.colmap["sequence"]])
            record[self.colmap["sequence_quality"]] = quals
        if not self.dry_run and not self.writer:
            if self.fmt in ["csv", "csvgz"]:
                self.writer = csv.DictWriter(
                    self.handle, fieldnames=record.keys(), lineterminator="\n")
                self.writer.writeheader()
            elif self.fmt in ["tsv", "tsvgz"]:
                self.writer = csv.DictWriter(
                    self.handle, fieldnames=record.keys(), lineterminator="\n", delimiter="\t")
                self.writer.writeheader()
        if self.fmt in ["csv", "tsv", "csvgz", "tsvgz"]:
            if not self.dry_run:
                self.writer.writerow(record)
        elif self.fmt in ["fa", "fagz"]:
            seqrecord = self.encode_record(record)
            if not self.dry_run:
                SeqIO.write(seqrecord, self.handle, "fasta-2line")
        elif self.fmt in ["fq", "fqgz"]:
            seqrecord = self.encode_record(record)
            if not "phred_quality" in seqrecord.letter_annotations:
                LOGGER.warning(
                    "No quality scores available, using default dummy value: %s",
                    DEFAULT_DUMMY_QUAL)
                self.dummyqual = DEFAULT_DUMMY_QUAL
            seqrecord = self.encode_record(record)
            if not self.dry_run:
                SeqIO.write(seqrecord, self.handle, "fastq")
        else:
            raise ValueError(f"Unknown format {self.fmt}")
