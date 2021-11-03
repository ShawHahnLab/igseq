"""
Get raw reads from an Illumina run directory.

This is a wrapper around Illumina's bcl2fastq program, just with some automatic
handling of the weirdness of the IgSeq protocol.  A dummy sample sheet is used
to try to ensure that all reads are placed in the Undetermined files (since
demultiplexing must be performed separately anyway) and index reads are enabled
so they can be used later during demultiplexing.

This step can be skipped if you already have all of your reads in single trio
of I1/R1/R2 fastq.gz files.
"""

import logging
import subprocess
from tempfile import NamedTemporaryFile
from pathlib import Path
from csv import DictReader
from . import util

LOGGER = logging.getLogger(__name__)

BCL2FASTQ = "bcl2fastq"

def getreads(path_input, dir_out, path_counts="", threads_load=1, threads_proc=1, dry_run=False):
    """Get reads directly from Illumina run directory.

    path_input: path to an Illumina run directory
    dir_out: path to write fastq.gz files to
    path_counts: path to csv to write read counts to.  If an empty
                 string, dir_out/getreads.counts.csv is used.  If None, this
                 file isn't written.
    threads_load: number of threads for parallel loading
    threads_proc: number of threads for parallel processing
    dry_run: If True, don't actually call any commands or write any files.

    This is a wrapper around Illumina's bcl2fastq with its built-in
    demultiplexing features (mostly) disabled and the I1 output enabled.  All
    reads should end up in the I1/R1/R2 Undetermined files in the output
    directory.

    There's a chance reads will be written to a dummy sample's output files due
    to bcl2fastq's insistence of having at least one sample defined in order to
    write an I1 file.  We use an implausible N-heavy sequence for that dummy
    barcode, so if that happens it most likely indicates a quality issue.  If
    any reads are written for that sample, a message will be logged complaining
    about it.
    """
    path_input = Path(path_input)
    if not dir_out:
        name = path_input.name
        dir_out = Path("analysis/reads") / name
    else:
        dir_out = Path(dir_out)
    LOGGER.info("input: %s", path_input)
    LOGGER.info("output: %s", dir_out)

    if path_counts == "":
        path_counts = dir_out / "getreads.counts.csv"
        LOGGER.info("path_counts inferred: %s", path_counts)
    else:
        LOGGER.info("path_counts: %s", path_counts)

    # Choose a bcl2fastq log level based on our own curent level.  By default
    # we'll be one increment less verbose, more or less.
    log_levels = ["TRACE", "DEBUG", "INFO", "WARNING", "ERROR", "FATAL", "NONE"]
    log_level = log_levels[
        min(len(log_levels)-1, max(0, round(logging.getLogger().getEffectiveLevel()/10 + 1)))]
    LOGGER.info("bcl2fastq log level: %s", log_level)
    if not dry_run:
        # We don't want any samples defined since we're just using the program
        # to convert .bcl to .fastq.gz (no demultiplexing).  We also don't want
        # any weirdness coming from the sample sheet like the ReverseComplement
        # setting.  Apparently the sample sheet can be completely empty as far
        # as bcl2fastq is concerned, but if there's isn't at least one sample
        # defined it won't create the Undetermined I1 file.
        # So we'll make a stub sample sheet just with that minimum info.
        with NamedTemporaryFile("wt") as sample_sheet:
            sample_sheet.write("[Data]\nSample_ID,index\nsample1,NANANANA\n")
            sample_sheet.flush()
            args = [
                "--runfolder-dir", path_input,
                "--output-dir", dir_out,
                # By default it'll try to write some files to the input run directory.
                # we'll send those to the output directory instead.
                "--interop-dir", dir_out / "InterOp",
                # We want the I1 file, so we enable that here.
                "--create-fastq-for-index-reads",
                # don't bother doing a fuzzy match on any supposed barcodes
                "--barcode-mismatches", 0,
                "--sample-sheet", sample_sheet.name,
                # parallel processing during loading can help a bit in my tests
                "--loading-threads", threads_load,
                # parallel processing does *not* help during the bcl2fastq
                # demultiplexing step, go figure, when we don't have any
                # demultiplexing to perform here
                "--demultiplexing-threads", 1,
                # parallel processing in the processing step helps quite a bit
                "--processing-threads", threads_proc,
                # help text says "this must not be higher than number of
                # samples"
                "--writing-threads", 1,
                "--min-log-level", log_level]
            try:
                _run_bcl2fastq(args)
            except subprocess.CalledProcessError as err:
                msg = f"bcl2fastq exited with code {err.returncode}"
                LOGGER.critical(msg)
                raise util.IgSeqError(msg) from err
        # We should have three (I1/R1/R2 Undetermined) fastq.gz files.  If
        # there are more, maybe some reads got dumped to other "samples."  If
        # there are fewer...?
        fqgzs = list(dir_out.glob("*.fastq.gz"))
        if len(fqgzs) > 3:
            LOGGER.warning("more .fastq.gz outputs than expected.")
        elif len(fqgzs) < 3:
            LOGGER.critical("Missing .fastq.gz outputs")
        if not dry_run and path_counts:
            try:
                cts = count_bcl2fastq_reads(dir_out/"Stats/FastqSummaryF1L1.txt")
            except FileNotFoundError as err:
                LOGGER.critical("Missing FastqSummaryF1L1.txt output")
            else:
                cts = [{"Category": "getreads", "Item": k, "NumSeqs": v} for k, v in cts.items()]
                util.save_counts(path_counts, cts)

def count_bcl2fastq_reads(summary_txt):
    """Summarize read counts from bcl2fastq's Stats/FastqSummaryF1L1.txt.

    This gives a simple, flat dictionary output with read counts by category.

    unassigned-raw: Number of reads not assigned to a sample, before filtering
    unassigned-pf: Number of reads not assigned to a sample that pass filters
    extra-raw: Number of reads assigned to a sample, before filtering
    extra-pf: Number of reads assigned to a sample that pass filters
    """

    cts = {"unassigned-raw": 0, "unassigned-pf": 0, "extra-raw": 0, "extra-pf": 0}
    with open(summary_txt) as f_in:
        reader = DictReader(f_in, delimiter="\t")
        for row in reader:
            if row["SampleNumber"] == "0":
                cts["unassigned-raw"] += int(row["NumberOfReadsRaw"])
                cts["unassigned-pf"] += int(row["NumberOfReadsPF"])
            else:
                cts["extra-raw"] += int(row["NumberOfReadsRaw"])
                cts["extra-pf"] += int(row["NumberOfReadsPF"])
    return cts

def _run_bcl2fastq(args):
    args = [BCL2FASTQ] + [str(arg) for arg in args]
    LOGGER.info("bcl2fastq command: %s", args)
    subprocess.run(args, check=True)
