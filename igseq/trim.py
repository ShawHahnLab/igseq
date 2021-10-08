"""
Trim adapter/low-quality parts from sequences.
"""

import re
import logging
import subprocess
from pathlib import Path
from . import util
LOGGER = logging.getLogger(__name__)

CUTADAPT = "cutadapt"

# https://cutadapt.readthedocs.io/en/stable/guide.html#quality-trimming
# https://cutadapt.readthedocs.io/en/stable/algorithms.html#quality-trimming-algorithm
def trim(paths_input, path_samples, dir_out="", path_counts="", species="rhesus", sample_name=None, min_length=50, quality_cutoff=15, dry_run=False, threads=1):
    """

    paths_input: list of paths to demultiplexed samples (one directory or a
                 R1/R2 pair)
    path_samples: path to samples CSV file
    dir_out: path to write trimmed R1/R2 fastq.gz files to
    path_counts: path to write read counts CSV
    """

    # paths_input can be:
    #    dir/to/pairs
    #    single_r1, single_r2
    # sample_name:
    #    always inferred, if paths_input is dir
    #    can be given (must be?) if paths_input is single pair

    if not dir_out:
        dir_out = util.default_path(paths_input, "trim")
    samples = util.load_samples(path_samples)
    samples = util.assign_barcode_seqs(samples)

    #### parse paths_input
    if len(paths_input) == 1 and Path(paths_input[0]).is_dir():
        # detect R1/R2 pairs
        path = Path(paths_input[0])
        pairs = []
        for fpr1, fpr2 in zip(path.glob("*.R1.fastq.gz"), path.glob("*.R2.fastq.gz")):
            prefix1 = re.sub(r"\.R1\.fastq\.gz", "", fpr1.name)
            prefix2 = re.sub(r"\.R2\.fastq\.gz", "", fpr2.name)
            if prefix1 == "unassigned":
                continue
            if prefix1 != prefix2:
                raise ValueError
            pairs.append({"R1": fpr1, "R2": fpr2})
    elif len(paths_input) == 2:
        # take R1/R2 as given
        pairs = [{"R1": Path(paths_input[0]), "R2": Path(paths_input[1])}]
    else:
        raise ValueError

    LOGGER.info("input samples: %s", path_samples)
    LOGGER.info("output dir: %s", dir_out)
    LOGGER.info("output counts: %s", path_counts)
    if not dry_run:
        dir_out.mkdir(parents=True, exist_ok=True)
    for pair in pairs:
        # what sample attributes go with this file pair?
        if sample_name:
            # use name if one given
            samp_name = sample_name
            sample = [samp for samp in samples.values() if samp["Sample"] == samp_name][0]
        else:
            # otherwise infer from paths
            samp_name = re.match(r"(.*)\.R1\.fastq\.gz", pair["R1"].name).group(1)
            sample = [samp for samp in samples.values() if samp["Sample"] == samp_name][0]
        adapter_fwd = get_adapter_fwd(sample, species)
        adapter_rev = get_adapter_rev(sample)
        output_r1 = dir_out / f"{samp_name}.R1.fastq.gz"
        output_r2 = dir_out / f"{samp_name}.R2.fastq.gz"
        LOGGER.info("sample %s: Fwd Adapter: %s", samp_name, adapter_fwd)
        LOGGER.info("sample %s: Rev Adapter: %s", samp_name, adapter_rev)
        LOGGER.info("sample %s: R1 in: %s", samp_name, pair["R1"])
        LOGGER.info("sample %s: R2 in: %s", samp_name, pair["R2"])
        LOGGER.info("sample %s: R1 in: %s", samp_name, output_r1)
        LOGGER.info("sample %s: R2 in: %s", samp_name, output_r2)
        if not dry_run:
            cutadapt(
                pair["R1"], pair["R2"],
                output_r1, output_r2,
                adapter_fwd, adapter_rev,
                min_length, quality_cutoff, threads)

# cutadapt on a single pair
def cutadapt(r1_in, r2_in, r1_out, r2_out, adapter_fwd, adapter_rev, min_length,
        quality_cutoff, threads):
    args = [
        CUTADAPT, "--cores", threads,
        "--quality-cutoff", quality_cutoff,
        "--minimum-length", min_length,
        "-a", adapter_fwd,
        "-A", adapter_rev,
        "-o", r1_out,
        "-p", r2_out,
        r1_in,
        r2_in]
    args = [str(arg) for arg in args]
    subprocess.run(args, check=True)

def get_adapter_fwd(sample, species):
    """Get the adapter sequence to trim off the end of R1.
    sample: dictionary of sample attributes
    """
    # The PCR primer specific to the antibody type occurs just *after* the
    # beginning of R2 (in the 3' direction, that is), so we'll trim that off
    # the end of R1.
    for row in util.PRIMERS:
        if row["Species"] == species and row["Type"] == sample["Type"]:
            adapter = row["Seq"]
            break
    else:
        raise ValueError("Unknown antibody type %s" % sample["Type"])
    return util.revcmp(adapter)

def get_adapter_rev(sample):
    """Get the adapter sequence to trim off the end of R2.
    sample: dictionary of sample attributes
    """
    # Reverse is trickier than forward.  The 2nd round PCR Forward Primer 1 /
    # P5 Sequencing Primer Site occurs just *before* the beginning of R1, so we
    # could cut on that, but we don't want to leave a dangling barcode segment
    # to get paired back in when combining R1 and R2 later (it's preferable if
    # our paired sequences start exactly after the forward barcode).  So we'll
    # include the barcode for each specific sample as part of the adapter here.
    # Note that the "N" characters in the barcode sequences should be handled
    # just fine by cutadapt.
    bcfwd = sample["BarcodeFwdSeq"]
    # Reverse complement the barcode and prepend to the constant region.
    return util.revcmp(bcfwd) + util.revcmp(util.P5SEQ)
