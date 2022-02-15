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

from .record import RecordReader, RecordWriter

def convert(path_in, path_out, fmt_in=None, fmt_out=None, colmap=None, dummyqual=None, dry_run=False):
    with RecordReader(path_in, fmt_in, colmap, dry_run=dry_run) as reader, \
        RecordWriter(path_out, fmt_out, colmap, dummyqual=dummyqual, dry_run=dry_run) as writer:
        for record in reader:
            writer.write(record)
