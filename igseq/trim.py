"""
Trim adapter/low-quality parts from sequences.

By default this will remove any instances of the R2 adapter found toward the
end of R1 and the R1 adapter found toward the end of R2.  It will also insist
that the 5' RACE Anchor be found at the start of R1, discarding read pairs that
are mising the anchor.  The adapter sequences will be determined from the
barcodes used for each sample and the selected species.
"""

import re
import logging
from subprocess import Popen, PIPE, DEVNULL
import json
from pathlib import Path
from . import util
LOGGER = logging.getLogger(__name__)

CUTADAPT = "cutadapt"

DEFAULTS = {
    # For selecting constant region primer to trim from 3' end of R1
    "species": "rhesus",
    # cutadapt settings
    "min_length": 50,
    "quality_cutoff": 15
    }

# https://cutadapt.readthedocs.io/en/stable/guide.html#quality-trimming
# https://cutadapt.readthedocs.io/en/stable/algorithms.html#quality-trimming-algorithm
def trim(
    paths_input, path_samples, dir_out="", path_counts="", species=DEFAULTS["species"],
    sample_name=None, dry_run=False, **kwargs):
    """Trim sample-specific adapter sequences from one or more file pairs.

    paths_input: list of paths to demultiplexed samples (one directory or a
                 R1/R2 pair)
    path_samples: path to samples CSV file
    dir_out: path to write trimmed R1/R2 fastq.gz files to
    path_counts: path to csv to write per-sample read counts to.  If an empty
                 string, the R1 filename is used but with .counts.csv
                 instead of .R1.fastq.gz.  If None, this file isn't written.
    species: species name, for choosing appropriate constant region primer
             sequence
    sample_name: explicit sample name (default: infer from filenames)
    dry_run: If True, don't actually call any commands or write any files.
    kwargs: additional keyword arguments for cutadapt()
    """

    samples = util.load_samples(path_samples)
    # filter to just samples with both forward and reverse barocdes IDs
    # specified, and then match with the sequences
    # TODO should filter more sensibly to just relevant samples, and sanity
    # check for unique barcode combos and such
    samples = {k: v for k, v in samples.items() if v["BarcodeFwd"] and v["BarcodeRev"]}
    samples = util.assign_barcode_seqs(samples)

    # Initial setup for one or more pairs of R1/R2 input
    pairs = util.parse_multi_fqgz_paths(paths_input)
    if len(pairs) == 0:
        raise ValueError(f"No input files found for: {paths_input}")
    if len(pairs) == 1 and path_counts:
        pairs[0]["path_counts"] = path_counts
    elif path_counts:
        LOGGER.warning("multiple R1/R2 pairs found; ignoring counts path %s", str(path_counts))
    if len(pairs) == 1 and sample_name:
        pairs[0]["sample_name"] = sample_name
    elif sample_name:
        LOGGER.warning("multiple R1/R2 pairs found; ignoring sample name %s", str(sample_name))

    # Use a given output directory, or By default use an output directory based
    # on the input directory
    if dir_out:
        dir_out = Path(dir_out)
    else:
        dir_out = util.default_path(pairs[0], "trim")

    # Assign sample name and file paths for each R1/R2 pair
    for pair in pairs:
        if not pair.get("sample_name"):
            pair["sample_name"] = re.match(r"(.*)\.R1\.fastq\.gz", pair["R1"].name).group(1)
        samp = pair["sample_name"]
        if not pair.get("path_counts"):
            if path_counts is not None:
                pair["path_counts"] = dir_out / f"{samp}.trim.counts.csv"
            else:
                pair["path_counts"] = None
        pair["R1_out"] = dir_out / f"{samp}.R1.fastq.gz"
        pair["R2_out"] = dir_out / f"{samp}.R2.fastq.gz"
        pair["JSON_out_1"] = dir_out / f"{samp}.cutadapt1.json"
        pair["JSON_out_2"] = dir_out / f"{samp}.cutadapt2.json"

    # Loop over each pair and call cutadapt with the appropriate sequences to
    # trim
    LOGGER.info("input samples: %s", path_samples)
    LOGGER.info("output dir: %s", dir_out)
    LOGGER.info("output counts: %s", path_counts)
    LOGGER.info("5PIIA seq: %s", util.ANCHOR5P)
    if not dry_run:
        dir_out.mkdir(parents=True, exist_ok=True)
    for pair in pairs:
        # what sample attributes go with this file pair?
        sample = [samp for samp in samples.values() if samp["Sample"] == pair["sample_name"]][0]
        adapter_fwd = get_adapter_fwd(sample, species)
        adapter_rev = get_adapter_rev(sample)
        samp = pair["sample_name"]
        LOGGER.info("sample %s: Fwd Adapter: %s", samp, adapter_fwd)
        LOGGER.info("sample %s: Rev Adapter: %s", samp, adapter_rev)
        LOGGER.info("sample %s: R1 in: %s", samp, pair["R1"])
        LOGGER.info("sample %s: R2 in: %s", samp, pair["R2"])
        LOGGER.info("sample %s: R1 out: %s", samp, pair["R1_out"])
        LOGGER.info("sample %s: R2 out: %s", samp, pair["R2_out"])
        LOGGER.info("sample %s: JSON out 1: %s", samp, pair["JSON_out_1"])
        LOGGER.info("sample %s: JSON out 2: %s", samp, pair["JSON_out_2"])
        if pair["path_counts"]:
            LOGGER.info("sample %s: counts: %s", samp, pair["path_counts"])
        # tell cutadapt to be quiet if we're at a less-verbose log level, but
        # not quiet if we're at a more verbose log level (in effect this means
        # we'd have to be at DEBUG to get quiet=False)
        quiet = logging.getLogger().getEffectiveLevel() >= logging.INFO
        # combine 5PIIA and forward adapter to get the sequence cutadapt will
        # expect at the start and end of R1, respectively.  5PIIA will be
        # anchored so it is implicitly requred.  The other sequence may or may
        # not be found.
        adapter_fwd = f"^{util.ANCHOR5P}...{adapter_fwd}"
        if not dry_run:
            trim_pair(
                pair["R1"], pair["R2"], pair["R1_out"], pair["R2_out"],
                pair["JSON_out_1"], pair["JSON_out_2"],
                adapter_fwd, adapter_rev,
                discard_untrimmed = True,
                quiet=quiet,
                **kwargs)
            if pair["path_counts"]:
                cts= _count_cutadapt_reads(pair["JSON_out_1"], pair["JSON_out_2"])
                cts = [{"Category": "trim", "Sample": samp, "Item": k, "NumSeqs": v} for k, v in cts.items()]
                util.save_counts(pair["path_counts"], cts)

def trim_pair(r1_in, r2_in, r1_out, r2_out, json1_out, json2_out, adapter_fwd, adapter_rev,
    discard_untrimmed=True,
    min_length=DEFAULTS["min_length"],
    quality_cutoff=DEFAULTS["quality_cutoff"],
    threads=1, quiet=True):
    """Trim adapters on paired-end fastq.gz files with cutadapt.

    This chains two cutadapt calls together to get the trimming and filtering behavior.

    r1_in: Path to input R1 fastq.gz file
    r2_in: Path to input R2 fastq.gz file
    r1_out: Path to output R1 fastq.gz file
    r2_out: Path to output R2 fastq.gz file
    json1_out: Path to output JSON-format report file for first cutadapt
               command.  If empty or None the report is not written.
    json2_out: Path to output JSON-format report file for second cutadapt
               command.  If empty or None the report is not written.
    adapter_fwd: Sequence to trim from 3' end of R1 (for -a argument)
    adapter_rev: Sequence to trim from 3' end of R2 (for -A argument)
    discard_untrimmed: should reads without required adapters found be
                       discarded?
    """
    # arguments shared by both cutadapt commands
    args_common = [CUTADAPT, "--interleaved", "--cores", threads]
    if quiet:
        args_common.append("--quiet")
    # first command: trim R2 adapter and interleave
    args1 = args_common + ["-A", adapter_rev, r1_in, r2_in]
    if json1_out:
        args1.extend(["--json", json1_out])
    args1 = [str(arg) for arg in args1]
    # second command: trim linked R1 adapter and filter those that don't start
    # with the expected sequence.  also apply all other filtering criteria.
    args2 = args_common + [
        "-a", adapter_fwd,
        "--quality-cutoff", quality_cutoff,
        "--minimum-length", min_length,
        "-o", r1_out, "-p", r2_out, "-"]
    if discard_untrimmed:
        args2.append("--discard-untrimmed")
    if json2_out:
        args2.extend(["--json", json2_out])
    args2 = [str(arg) for arg in args2]
    _run_cutadapt_pair(args1, args2)

def _run_cutadapt_pair(args1, args2):
    # Pipe the output of the first command into the second, and wait for the
    # second to finish
    with \
        Popen(args1, stdout=PIPE, stderr=DEVNULL) as proc1, \
        Popen(args2, stdin=proc1.stdout) as proc2:
        proc2.wait()

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

def _count_cutadapt_reads(json_path_1, json_path_2):
    """Pull out some read counts of interest from cutadapt's JSON reports.

    This gives a simple, flat dictionary output with a subset of the available
    read counts.
    """

    with open(json_path_1) as f_in:
        report1 = json.load(f_in)
    with open(json_path_2) as f_in:
        report2 = json.load(f_in)
    cts1 = report1["read_counts"]
    cts2 = report2["read_counts"]
    output = {
        "input": cts1["input"],
        "output": cts2["output"],
        "too_short": cts2["filtered"]["too_short"],
        "discard_untrimmed": cts2["filtered"]["discard_untrimmed"],
        "read1_with_adapter": cts2["read1_with_adapter"],
        "read2_with_adapter": cts1["read2_with_adapter"]}
    return output
