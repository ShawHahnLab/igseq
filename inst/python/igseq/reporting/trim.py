"""
Summarizing and reporting helper functions - adapter trimming.
"""

import csv
import gzip
import logging
from Bio import SeqIO
from cutadapt import qualtrim
LOGGER = logging.getLogger(__name__)

def make_qualtrim_grid(fqgzs_in, qual_breaks=None, len_breaks=None):
    """Tally how many reads would be trimmed to what lengths at what quality cutoffs.

    This runs per-sample but takes a list of fastq.gz paths in to support
    chunked data.  Just give a list of one to disregard that.

    This produces a grid (in the form of a per-quality-beak dictionary with
    each value being a per-length dictionary containing read counts).  Each
    per-quality-value dictionary contains numbers summing to the total number
    of reads, since all reads are considered for each quality value.

    https://cutadapt.readthedocs.io/en/stable/algorithms.html#quality-trimming-algorithm
    https://github.com/marcelm/cutadapt/blob/master/src/cutadapt/qualtrim.pyx#L6
    """

    # Go from highest qual cutoff to lowest
    if not qual_breaks:
        qual_breaks = list(range(0, 42))
    qual_breaks = sorted(qual_breaks)[::-1]
    # Go from lowest to highets length
    if not len_breaks:
        len_breaks = list(range(0, 401))
    len_breaks = sorted(len_breaks)

    tally = {b: {lenb: 0 for lenb in len_breaks} for b in qual_breaks}
    for fqgz_in in fqgzs_in:
        with gzip.open(fqgz_in, "rt") as f_in:
            for record in SeqIO.parse(f_in, "fastq"):
                qual = [chr(val+33) for val in record.letter_annotations["phred_quality"]]
                qual = ''.join(qual)
                for cutoff in qual_breaks:
                    _, trim = qualtrim.quality_trim_index(qual, 0, cutoff)
                    diffs = [abs(len_break - trim) for len_break in len_breaks]
                    len_break = len_breaks[diffs.index(min(diffs))]
                    tally[cutoff][len_break] += 1
    return tally

def make_qualtrim_csv(fqgzs_in, fp_csv, qual_breaks=None, len_breaks=None):
    """CSV writer for make_qualtrim_grid.

    This writes a CSV file with quality cutoffs on rows, sequence lengths on
    columns, and counts of occurrences in each cell.
    """
    grid = make_qualtrim_grid(fqgzs_in, qual_breaks, len_breaks)
    keys_qual = sorted(grid.keys())
    keys_len = sorted(grid[keys_qual[0]].keys())
    with open(fp_csv, "wt") as f_out:
        writer = csv.writer(f_out)
        header = ["Q"] + keys_len
        writer.writerow(header)
        for key in keys_qual:
            row = [key] + [grid[key][key2] for key2 in keys_len]
            writer.writerow(row)
