"""
Convert between various sequence and tabular file formats.

Input and output formats are by default inferred from filenames but can be
given explicitly if needed.  The formats are:

    fa:    FASTA
    fagz:  gzipped FASTA
    fq:    FASTQ
    fqgz:  gzipped FASTQ
    csv:   comma-separated values
    csvgz: gzipped comma-separated values
    tsv:   tab-separated values
    tsvgz: gzipped tab-separated values
"""

import logging
from .record import RecordReader, RecordWriter

def convert(path_in, path_out, fmt_in=None, fmt_out=None, dummyqual=None):
    with RecordReader(path_in, fmt_in) as reader, \
        RecordWriter(path_out, fmt_out, dummyqual=dummyqual) as writer:
        for record in reader:
            writer.write(record)
