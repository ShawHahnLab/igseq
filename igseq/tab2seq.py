"""
Convert TSV/CSV to FASTA/FASTQ.

If converting to FASTQ the quality column must be stored as PHRED+33 (the text
shown in modern Illumina FASTQ files).
"""

import sys
import logging
from csv import DictReader
from pathlib import Path
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from . import util
LOGGER = logging.getLogger(__name__)

TAB_EXTS = {
    ".csv": "csv",
    ".tsv": "tsv",
    ".tab": "tsv"}

SEQ_EXTS = {
    ".fasta": "fasta",
    ".fa": "fasta",
    ".fastq": "fastq",
    ".fq": ".fastq"}

DASH = Path("-")

def tab2seq(tab_path_in, seq_path_out, seq_col, seq_id_col,
    seq_desc_col=None, qual_col=None, tab_fmt=None, seq_fmt=None):
    tab_path_in = Path(tab_path_in)
    seq_path_out = Path(seq_path_out)
    tab_ext = tab_path_in.suffix.lower()
    seq_ext = seq_path_out.suffix.lower()
    tab_fmt_det = TAB_EXTS.get(tab_ext)
    seq_fmt_det = SEQ_EXTS.get(seq_ext)
    if not tab_fmt_det and not tab_fmt:
        raise util.IgSeqError("could not detect input format automatically")
    if not seq_fmt_det and not seq_fmt:
        raise util.IgSeqError("could not detect output format automatically")
    tab_fmt = tab_fmt or tab_fmt_det
    seq_fmt = seq_fmt or seq_fmt_det
    delim = {"csv": ",", "tsv": "\t"}[tab_fmt]
    fmt = {"fasta": "fasta-2line", "fastq": "fastq"}[seq_fmt]
    if seq_fmt == "fastq" and not qual_col:
        raise util.IgSeqError("quality column name required for FASTQ output")
    phred_decode = lambda txt: [ord(letter)-33 for letter in txt]
    try:
        f_in = sys.stdin if tab_path_in == DASH else open(tab_path_in)
        f_out = sys.stdout if seq_path_out == DASH else open(seq_path_out, "wt")
        reader = DictReader(f_in, delimiter=delim)
        for row in reader:
            record = SeqRecord(
                id=row[seq_id_col],
                description=row[seq_desc_col] if seq_desc_col else "",
                seq=Seq(row[seq_col]))
            if fmt == "fastq":
                record.letter_annotations["phred_quality"] = phred_decode(row[qual_col])
            SeqIO.write(record, f_out, fmt)
    finally:
        if tab_path_in != DASH:
            f_in.close()
        if seq_path_out != DASH:
            f_out.close()
