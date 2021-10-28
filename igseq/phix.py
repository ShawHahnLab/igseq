"""
Map reads unassigned after demultiplexing to the PhiX genome.
"""

import logging
import tempfile
from subprocess import Popen, PIPE, DEVNULL
from pathlib import Path
from . import util
from .util import DATA, default_path, save_counts, parse_fqgz_paths

LOGGER = logging.getLogger(__name__)

REF = DATA / "phix/phix.fasta"

BWA = "bwa"
SAMTOOLS = "samtools"

def phix(paths_input, bam_out="", counts_out="", dry_run=False, threads=1):
    """Map reads to the PhiX genome.

    paths_input: list of inputs, either one directory or two files, for R1 and
                 R2, in order.
    bam_out: path to bam file for output.  By default a path is inferred from
             the input filenames.
    path_counts: path to csv to write per-sample read counts to.  If an empty
                 string, the bam_out filename is used but with .counts.csv
                 instead of .bam.  If None, this file isn't written.
    threads: number of threads to use for parallel processing.
    """
    try:
        paths_input = parse_fqgz_paths(paths_input, ["R1", "R2"], "unassigned")
    except ValueError as err:
        raise ValueError(
            "input path should be one directory or two (R1/R2) fastq.gz files.") from err
    if not bam_out:
        bam_out = default_path(paths_input.values(), "phix", "phix.bam")
    else:
        bam_out = Path(bam_out)
    if counts_out == "":
        counts_out = bam_out.with_suffix(".counts.csv")
    LOGGER.info("input R1: %s", paths_input["R1"])
    LOGGER.info("input R2: %s", paths_input["R2"])
    LOGGER.info("output bam: %s", bam_out)
    LOGGER.info("output counts: %s", counts_out)
    if not dry_run:
        bam_out.parent.mkdir(parents=True, exist_ok=True)
        map_reads(REF, paths_input["R1"], paths_input["R2"], bam_out, threads)
        if counts_out:
            counts_out.parent.mkdir(parents=True, exist_ok=True)
            if counts_out:
                num_reads = _count_bam_reads(bam_out)
                save_counts(
                    counts_out,
                    [{"Category": "phix", "Item": "mapped", "NumSeqs": num_reads}])

def map_reads(ref_path, r1_path, r2_path, bam_out, threads):
    """Map paired reads to reference (e.g. PhiX genome).

    ref_path: Path to FASTA of reference to map to.  The associated BWA index
              files need to exist already.
    r1_path: Path to forward reads fastq.gz
    r1_path: Path to reverse reads fastq.gz
    bam_out: Path to write sorted BAM file containing mapped reads.
    threads: number of threads to use for parallel processing.
    """

    with tempfile.TemporaryDirectory() as samtmp:
        cmd_bwa = [BWA, "mem", "-t", threads, ref_path, r1_path, r2_path]
        cmd_samtools_view = [SAMTOOLS, "view", "-b", "-F", "0x4"]
        cmd_samtools_sort = [SAMTOOLS, "sort", "-T", samtmp, "-@", threads]
        cmd_bwa = [str(obj) for obj in cmd_bwa]
        cmd_samtools_view = [str(obj) for obj in cmd_samtools_view]
        cmd_samtools_sort = [str(obj) for obj in cmd_samtools_sort]
        with open(bam_out, "wb") as f_out, \
            Popen(cmd_bwa, stdout=PIPE, stderr=DEVNULL) as bwa, \
            Popen(cmd_samtools_view, stdin=bwa.stdout, stdout=PIPE) as sam_view, \
            Popen(cmd_samtools_sort, stdin=sam_view.stdout, stdout=f_out) as sam_sort:
            sam_sort.wait()
            for proc in [bwa, sam_view, sam_sort]:
                if proc.returncode:
                    LOGGER.critical("%s exited with code %d", proc.args[0], proc.returncode)
                    raise util.IgSeqError(proc.args[0] + " crashed")

def _count_bam_reads(bam_path):
    # https://qnot.org/2012/04/14/counting-the-number-of-reads-in-a-bam-file/
    # Except that counts R1 and R2 separately.  Is there a better way than the
    # below to count read pairs with either/both R1/R2 mapping?
    with Popen([SAMTOOLS, "fasta", "-n", bam_path], stdout=PIPE, stderr=DEVNULL, text=True) as sam:
        desclines = set()
        for line in sam.stdout:
            if line.startswith(">"):
                desclines.add(line)
        return len(desclines)
