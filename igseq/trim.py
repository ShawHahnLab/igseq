"""
Trim adapter/low-quality parts from sequences.

By default this will remove any instances of the R2 adapter found toward the
end of R1 and the R1 adapter found toward the end of R2.  It will also insist
that the 5' RACE Anchor be found at the start of R1, discarding read pairs that
are missing the anchor.  The adapter sequences will be determined from the
barcodes used for each sample and the selected species.

The R1 adapter to trim from R2 is the forward barcode sequence at the start of
R1 and the constant P5 sequence just upstream of that.

The R2 adapter to trim from R1 is whatever constant region primers are
applicable based on the supplied sample metadata and specified species.  If no
species is specified and/or no chain type is specified via the sample metadata,
it will recognize all primers that could be applicable.

Any command-line arguments not recognized here are passed as-is to the
cutadapt command, like the igblast command allows.  See cutadapt --help for
those options.
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
    paths_input, path_samples, dir_out="", path_counts="", *,
    species=DEFAULTS["species"], extra_cutadapt_args=None, sample_name=None, dry_run=False,
    **kwargs):
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
    extra_cutadapt_args: optional list of command-line arguments to include in
                         the second cutadapt call
    dry_run: If True, don't actually call any commands or write any files.
    kwargs: additional keyword arguments for trim_pair()
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
    LOGGER.info("extra cutadapt arguments: %s", extra_cutadapt_args)
    if not dry_run:
        dir_out.mkdir(parents=True, exist_ok=True)
    for pair in pairs:
        # what sample attributes go with this file pair?
        sample = [samp for samp in samples.values() if samp["Sample"] == pair["sample_name"]][0]
        adapters_fwd = get_adapters_fwd(sample, species)
        adapter_rev = get_adapter_rev(sample)
        samp = pair["sample_name"]
        LOGGER.info("sample %s: Fwd Adapter: %s", samp, adapters_fwd)
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
        # combine 5PIIA and forward adapter(s) to get the sequence cutadapt
        # will expect at the start and end of R1, respectively.  5PIIA will be
        # anchored so it is implicitly requred.  The other sequence, on the 3'
        # end, may or may not be found.
        adapters_fwd_lnk = {
            f"race_anchor_and_{k}": f"^{util.ANCHOR5P}...{v}" for k, v in adapters_fwd.items()}
        if not dry_run:
            trim_pair(
                pair["R1"], pair["R2"], pair["R1_out"], pair["R2_out"],
                pair["JSON_out_1"], pair["JSON_out_2"],
                adapters_fwd_lnk, adapter_rev,
                discard_untrimmed = True,
                quiet=quiet,
                extra_cutadapt_args=extra_cutadapt_args,
                **kwargs)
            if pair["path_counts"]:
                cts= _count_cutadapt_reads(pair["JSON_out_1"], pair["JSON_out_2"])
                cts = [{"Category": "trim", "Sample": samp, "Item": k, "NumSeqs": v} \
                    for k, v in cts.items()]
                util.save_counts(pair["path_counts"], cts)

def trim_pair(r1_in, r2_in, r1_out, r2_out, json1_out, json2_out, adapters_fwd, adapter_rev,
    *, extra_cutadapt_args=None, discard_untrimmed=True,
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
    adapters_fwd: Dictionary of sequences to trim from 3' end of R1 (for -a
                  arguments)
    adapter_rev: Sequence to trim from 3' end of R2 (for -A argument)
    extra_cutadapt_args: list of additional arguments to pass to the second
                         cutadapt call
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
    # second command: trim linked R1 adapter(s) and filter those that don't
    # start with the expected sequence.  also apply all other filtering
    # criteria.
    args_adapt_fwd = []
    for name, seq in adapters_fwd.items():
        args_adapt_fwd += ["-a", f"{name}={seq}"]
    args2 = args_common + args_adapt_fwd + [
        "--quality-cutoff", quality_cutoff,
        "--minimum-length", min_length,
        "-o", r1_out, "-p", r2_out, "-"]
    if discard_untrimmed:
        args2.append("--discard-untrimmed")
    if json2_out:
        args2.extend(["--json", json2_out])
    if extra_cutadapt_args:
        args2.extend(extra_cutadapt_args)
    args2 = [str(arg) for arg in args2]
    _run_cutadapt_pair(args1, args2)

def _run_cutadapt_pair(args1, args2):
    # Pipe the output of the first command into the second, and wait for the
    # second to finish
    with \
        Popen(args1, stdout=PIPE, stderr=DEVNULL) as proc1, \
        Popen(args2, stdin=proc1.stdout) as proc2:
        proc2.wait()
        if proc1.returncode:
            LOGGER.critical("cutadapt proc 1 exited with code %s", proc1.returncode)
            raise util.IgSeqError("cutadapt crashed")
        if proc2.returncode:
            LOGGER.critical("cutadapt proc 2 exited with code %s", proc2.returncode)
            raise util.IgSeqError("cutadapt crashed")

def get_adapters_fwd(sample, species=None):
    """Get the adapter sequence to trim off the end of R1.

    Returns a dictionary of all applicable options, with adapter names as keys
    as sequences as values.  If species is specified and chain type is given
    via the Type key of the sample dictionary, this will just be a single name
    for a single adapter sequence, like {"rhesus_gamma": "..."}.

    sample: dictionary of sample attributes
    species: species name from the 
    """
    # The PCR primer specific to the antibody type occurs just *after* the
    # beginning of R2 (in the 3' direction, that is), so we'll trim that off
    # the end of R1.
    chain_type = sample.get("Type")
    options = [row for row in util.PRIMERS if not chain_type or row["Type"] == chain_type]
    if not options:
        raise util.IgSeqError(f"Unknown antibody chain type {chain_type}")
    matches = {}
    for row in options:
        if species is None or row["Species"] == species:
            key = row["Species"] + "_" + row["Type"]
            matches[key] = util.revcmp(row["Seq"])
    if not matches:
        raise util.IgSeqError(f"Unknown species {species}")
    return matches

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

    with open(json_path_1, encoding="UTF8") as f_in:
        report1 = json.load(f_in)
    with open(json_path_2, encoding="UTF8") as f_in:
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
