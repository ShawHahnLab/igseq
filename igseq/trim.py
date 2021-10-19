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
import subprocess
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
    pairs = __parse_multi_fqgz_paths(paths_input)
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
        pair["JSON_out"] = dir_out / f"{samp}.cutadapt.json"

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
        LOGGER.info("sample %s: JSON out: %s", samp, pair["JSON_out"])
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
            cutadapt(
                pair["R1"], pair["R2"],
                pair["R1_out"], pair["R2_out"], pair["JSON_out"],
                adapter_fwd, adapter_rev,
                discard_untrimmed = True,
                quiet=quiet,
                **kwargs)
            if pair["path_counts"]:
                LOGGER.info("sample %s: counts: %s", samp, pair["path_counts"])
            if not dry_run:
                cts= _count_cutadapt_reads(pair["JSON_out"])
                cts = [{"Category": "trim", "Item": k, "NumSeqs": v} for k, v in cts.items()]
                util.save_counts(pair["path_counts"], cts)

def __parse_multi_fqgz_paths(paths_input):
    if len(paths_input) == 1 and Path(paths_input[0]).is_dir():
        # detect R1/R2 pairs
        path = Path(paths_input[0])
        pairs = []
        for fpr1, fpr2 in zip(
            sorted(path.glob("*.R1.fastq.gz")),
            sorted(path.glob("*.R2.fastq.gz"))):
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
    return pairs

def cutadapt(r1_in, r2_in, r1_out, r2_out, json_out, adapter_fwd, adapter_rev,
        discard_untrimmed=False,
        min_length=DEFAULTS["min_length"],
        quality_cutoff=DEFAULTS["quality_cutoff"],
        threads=1, quiet=True):
    """Call cutadapt on paired-end fastq.gz files.

    r1_in: Path to input R1 fastq.gz file
    r2_in: Path to input R2 fastq.gz file
    r1_out: Path to output R1 fastq.gz file
    r2_out: Path to output R2 fastq.gz file
    json_out: Path to output JSON-format report file.  If empty or None the
              report is not written.
    adapter_fwd: Sequence to trim from 3' end of R1 (for -a argument)
    adapter_rev: Sequence to trim from 3' end of R2 (for -A argument)
    discard_untrimmed: should reads without required adapters found be
                       discarded?
    """
    args = [
        CUTADAPT, "--cores", threads,
        "--quality-cutoff", quality_cutoff,
        "--minimum-length", min_length,
        "-a", adapter_fwd,
        "-A", adapter_rev] + \
        (["--discard-untrimmed"] if discard_untrimmed else []) + \
        ["-o", r1_out,
        "-p", r2_out] + \
        (["--json", json_out] if json_out else []) + \
        (["--quiet"] if quiet else []) + \
        [r1_in, r2_in]
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

def _count_cutadapt_reads(json_path):
    """Pull out some read counts of interest from cutadapt's JSON report.

    This gives a simple, flat dictionary output with a subset of the available
    read counts.
    """

    with open(json_path) as f_in:
        report = json.load(f_in)
    cts = report["read_counts"]
    output = {
        "input": cts["input"],
        "output": cts["output"],
        "too_short": cts["filtered"]["too_short"],
        "discard_untrimmed": cts["filtered"]["discard_untrimmed"],
        "read1_with_adapter": cts["read1_with_adapter"],
        "read2_with_adapter": cts["read2_with_adapter"]}
    return output
