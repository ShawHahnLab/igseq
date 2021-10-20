"""
Merge read pairs.
"""

import re
import sys
import logging
from pathlib import Path
import subprocess
from subprocess import Popen, PIPE, DEVNULL
from . import util

LOGGER = logging.getLogger(__name__)
PEAR = "pear"

def merge(paths_input, dir_out="", path_counts="", dry_run=False, threads=1):
    """Merge reads from one or more file pairs.

    paths_input: list of paths to fastq.gz files (one directory or a
                 R1/R2 pair)
    dir_out: path to write merged fastq.gz file(s) to
    path_counts: path to csv to write per-sample read counts to.  If an empty
                 string, the R1 filename is used but with .counts.csv
                 instead of .R1.fastq.gz.  If None, this file isn't written.
    dry_run: If True, don't actually call any commands or write any files.
    threads: number of threads for parallel processing with PEAR.
    """
    pairs = util.parse_multi_fqgz_paths(paths_input)
    if len(pairs) == 0:
        raise ValueError(f"No input files found for: {paths_input}")
    if len(pairs) == 1 and path_counts:
        pairs[0]["path_counts"] = path_counts
    elif path_counts:
        LOGGER.warning("multiple R1/R2 pairs found; ignoring counts path %s", str(path_counts))
    # Use a given output directory, or By default use an output directory based
    # on the input directory
    if dir_out:
        dir_out = Path(dir_out)
    else:
        dir_out = util.default_path(pairs[0], "merge")
    # Assign sample name and file paths for each R1/R2 pair
    for pair in pairs:
        if not pair.get("sample_name"):
            pair["sample_name"] = re.match(r"(.*)\.R1\.fastq\.gz", pair["R1"].name).group(1)
        samp = pair["sample_name"]
        if not pair.get("path_counts"):
            if path_counts is not None:
                pair["path_counts"] = dir_out / f"{samp}.merge.counts.csv"
            else:
                pair["path_counts"] = None
        pair["out"] = dir_out / f"{samp}.fastq.gz"
        pair["log"] = dir_out / f"{samp}.pear.log"
    # Loop over each pair and call pear with the appropriate sequences to
    # trim
    LOGGER.info("output dir: %s", dir_out)
    if not dry_run:
        dir_out.mkdir(parents=True, exist_ok=True)

    for pair in pairs:
        samp = pair["sample_name"]
        LOGGER.info("sample %s: R1 in: %s", samp, pair["R1"])
        LOGGER.info("sample %s: R2 in: %s", samp, pair["R2"])
        LOGGER.info("sample %s: out: %s", samp, pair["out"])
        LOGGER.info("sample %s: log: %s", samp, pair["log"])
        if pair["path_counts"]:
            LOGGER.info("sample %s: counts: %s", samp, pair["path_counts"])
        if not dry_run:
            quiet = logging.getLogger().getEffectiveLevel() >= logging.INFO
            logtxt = pear(pair["R1"], pair["R2"], pair["out"], threads, pair["log"], quiet)
            if pair["path_counts"]:
                cts= _count_pear_reads(logtxt)
                cts = [{"Category": "merge", "Item": k, "NumSeqs": v} for k, v in cts.items()]
                util.save_counts(pair["path_counts"], cts)

def pear(r1_in, r2_in, outfile, threads, log_path=None, quiet=True):
    """Call PEAR to merge paired-end data.

    This will post-process pear's four uncompresed output files by compressing
    the assembled reads as the specified outfile and removing the others.
    pear's status text is returned and optionally written to a log file (if
    log_path is not None) and/or stdout here (if quiet is False).
    """
    Path(outfile).parent.mkdir(parents=True, exist_ok=True)
    args = [PEAR, "-f", r1_in, "-r", r2_in, "-o", outfile, "-j", threads]
    args = [str(arg) for arg in args]
    log_path = log_path or "/dev/null"
    logtxt = ""
    # text=False because pear uses carriage returns (0x0D) to update a status
    # line in-place and newlines (0x0A) as line breaks.
    with Popen(args, stdout=PIPE, text=False, bufsize=0) as pear_proc, open(log_path, "wt") as log:
        while True:
            txt = pear_proc.stdout.read(1).decode("ascii")
            if pear_proc.poll() is not None and not txt:
                break
            if not quiet:
                sys.stderr.write(txt)
            if log:
                log.write(txt)
            logtxt += txt
        pear_proc.wait()
        if pear_proc.returncode:
            LOGGER.error("pear exited with code %s", pear_proc.returncode)

    # compress the assembled output with the given filename and remove the
    # uncompressed version
    out_merged = f"{outfile}.assembled.fastq"
    with open(out_merged, "rb") as f_in, open(outfile, "wb") as f_out:
        subprocess.run("gzip", stdin=f_in, stdout=f_out, check=True)
    Path(out_merged).unlink()
    # remove the other output files
    for name in ["discarded", "unassembled.forward", "unassembled.reverse"]:
        Path(f"{outfile}.{name}.fastq").unlink()

    return logtxt

def _count_pear_reads(logtxt):
    """Parse read counts from stdout from PEAR run."""
    # We're looking for something like this:
    # Assembled reads ...................: 2,874,079 / 3,179,377 (90.398%)
    # Discarded reads ...................: 0 / 3,179,377 (0.000%)
    # Not assembled reads ...............: 305,298 / 3,179,377 (9.602%)
    cts = {}
    pat = r"^(Assembled|Discarded|Not assembled) reads [\.]*: ([0-9,]*) / ([0-9,]*).*$"
    for line in logtxt.split("\n"):
        match = re.match(pat, line)
        if match:
            item = match.group(1)
            val = re.sub(",", "", match.group(2))
            total = re.sub(",", "", match.group(3))
            if "input" not in cts:
                cts["input"] = total
            if item == "Assembled":
                item = "output"
            cts[item] = val
    return cts
