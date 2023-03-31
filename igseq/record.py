"""
Helper classes for records stored in FASTA/FASTQ/CSV/TSV with/without gzip.
"""

import re
import sys
import csv
import gzip
import logging
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from igseq import util

LOGGER = logging.getLogger(__name__)

DEFAULT_DUMMY_QUAL = "I"
DEFAULT_COLUMNS = {
    k: k for k in ["sequence_id", "sequence", "sequence_quality", "sequence_description"]}

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
    "csv": csv.DictReader,
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

    @property
    def tabular(self):
        return self.fmt in ["csv", "tsv", "csvgz", "tsvgz"]

    @staticmethod
    def _infer_fmt(path):
        try:
            try:
                # Ordinary file paths
                path = Path(path)
            except TypeError:
                # for example, already-open file handles
                path = Path(str(path.name))
        except AttributeError:
            # for example, StringIO objects
            return None
        ext = path.suffix.lower()
        ext2 = Path(path.stem).suffix.lower()
        fmt2 = FMT_EXT_MAP.get(ext2)
        if ext == ".gz" and fmt2:
            return fmt2 + "gz"
        return FMT_EXT_MAP.get(ext)

    def decode_record(self, obj):
        """Convert object (SeqRecord or dictionary) to dictionary."""
        if isinstance(obj, SeqRecord):
            # (Biopython has some weirdly asymmetric behavior with sequence IDs
            # and descriptions, so I'm trying to handle that part myself.)
            if obj.description != "<unknown description>":
                if obj.description.startswith(obj.id):
                    # If there's a space, use that as a separator and record
                    # the description after that space.  If not, just ID.
                    match = re.match("([^ ]+)( ?)(.*)", obj.description)
                    seq_id, spacer, seq_desc = match.groups()
                    if not spacer:
                        seq_desc = None
                else:
                    # If the description *doesn't* have the ID at the start,
                    # just take both as-is.
                    seq_id = obj.id
                    seq_desc = obj.description
            else:
                seq_id = obj.id
                seq_desc = None
            record = {
                self.colmap["sequence_id"]: seq_id,
                self.colmap["sequence"]: str(obj.seq)}
            quals = obj.letter_annotations.get("phred_quality")
            if quals:
                record[self.colmap["sequence_quality"]] = self.encode_phred(quals)
            if seq_desc is not None:
                record[self.colmap["sequence_description"]] = seq_desc
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
        if hasattr(self.pathlike, "fileno"):
            self.handle = self.pathlike
        elif self.pathlike == "-":
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
            # order columns like in DEFAULT_COLUMNS with any custom ones on the
            # end
            def colsort(val):
                for idx, key in enumerate(self.colmap):
                    if self.colmap[key] == val:
                        return idx
                return len(DEFAULT_COLUMNS)
            fieldnames = sorted(record.keys(), key=colsort)
            if self.fmt in ["csv", "csvgz"]:
                self.writer = csv.DictWriter(
                    self.handle, fieldnames=fieldnames, lineterminator="\n")
                self.writer.writeheader()
            elif self.fmt in ["tsv", "tsvgz"]:
                self.writer = csv.DictWriter(
                    self.handle, fieldnames=fieldnames, lineterminator="\n", delimiter="\t")
                self.writer.writeheader()
        if self.fmt in ["csv", "tsv", "csvgz", "tsvgz"]:
            if not self.dry_run:
                self.writer.writerow(record)
        elif self.fmt in ["fa", "fagz"]:
            self._write_fa(record)
        elif self.fmt in ["fq", "fqgz"]:
            self._write_fq(record)
        else:
            raise ValueError(f"Unknown format {self.fmt}")

    def _write_fa(self, record):
        seq = record[self.colmap["sequence"]]
        defline = record[self.colmap["sequence_id"]]
        desc = record.get(self.colmap["sequence_description"])
        if desc is not None:
            defline += f" {desc}"
        if not self.dry_run:
            self.handle.write(f">{defline}\n{seq}\n")

    def _write_fq(self, record):
        seq = record[self.colmap["sequence"]]
        defline = record[self.colmap["sequence_id"]]
        desc = record.get(self.colmap["sequence_description"])
        if self.colmap["sequence_quality"] in record:
            quals = record[self.colmap["sequence_quality"]]
        else:
            LOGGER.warning(
                "No quality scores available, using default dummy value: %s",
                DEFAULT_DUMMY_QUAL)
            quals = "".join(DEFAULT_DUMMY_QUAL * len(seq))
        if desc is not None:
            defline += f" {desc}"
        if not self.dry_run:
            self.handle.write(f"@{defline}\n{seq}\n+\n{quals}\n")
