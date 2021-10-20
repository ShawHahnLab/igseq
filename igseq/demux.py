"""
Demultiplex reads into per-sample files.
"""

import logging
from collections import defaultdict
from pathlib import Path
from csv import DictWriter
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from . import util
from .util import READS, BARCODES, assign_barcode_seqs
from .fastgzip import fast_gzip_open

LOGGER = logging.getLogger(__name__)

# The built-in gzip/zlib modules actually perform far slower than just piping
# data in and out of the gzip/gunzip binaries, and for hundreds of MBs or above
# it can make a big difference.  It looks like people have been grumbling about
# this for a while:
# https://www.reddit.com/r/Python/comments/2olhrf/fast_gzip_in_python/
# https://codebright.wordpress.com/2011/03/25/139/
# This wrapper tries to act as a drop-in replacement for gzip.open for our
# purposees here (no promises for other uses), just faster.  Any problems, just
# set this to gzip.open instead.
GZIP_OPEN = fast_gzip_open

# lists of all forward and reverse barcodes
BARCODES_FWD = [attrs["Seq"] for attrs in BARCODES if attrs["Direction"] == "F"]
BARCODES_REV = [attrs["Seq"] for attrs in BARCODES if attrs["Direction"] == "R"]

def demux(paths_input, path_samples, run_id=None, dir_out="", path_counts="",
    path_details=None, dry_run=False):
    """Demultiplex reads from combined I1/R1/R2 files into per-sample files.

    This is the high-level demux function based on file inputs and outputs.
    See demux_by_barcode for the low-level function.

    paths_input: list of inputs, either one directory or three files, for I1,
                 R1, and R2, in order.
    path_samples: path to CSV file with sample barcode information.  The Run
                  column will be used to identify just the relevant samples
                  here.  Give a run_id argument explicitly if needed.
    run_id: Illumina Run ID for these reads.  Assumed to be the last common
            parent directory of paths_input if not given.
    dir_out: output directory for all files.
    path_counts: path to csv to write per-sample read counts to.  If an empty
                 string, <dir_out>/demux.counts.csv is used.  If None, this
                 file isn't written.
    path_details: csv.gz file to write detailed per-read information.  If
                 an empty string, <dir_out>/demux.info.csv.gz is used.  If None,
                 this file isn't written.
    dry_run: If True, don't actually call any commands or write any files.
    """
    try:
        paths_input_trio = util.parse_fqgz_paths(paths_input)
    except ValueError as err:
        raise ValueError(
            "input path should be one directory or three (I1/R1/R2) fastq.gz files.") from err
    samples = util.load_samples(path_samples)
    if not run_id:
        run_id = util.common_parent(paths_input_trio.values()).name
    if not dir_out:
        dir_out = util.default_path(paths_input_trio.values(), "demux")
    else:
        dir_out = Path(dir_out)
    if path_counts == "":
        path_counts = Path(dir_out) / "demux.counts.csv"
    if path_details == "":
        path_details = Path(dir_out) / "demux.info.csv.gz"
    # filter to just the relevant run
    samples = {key: val for key, val in samples.items() if val["Run"] == run_id}
    if not samples:
        raise ValueError(f"No matching samples for run {run_id} found in {path_samples}")
    samples = assign_barcode_seqs(samples)
    for rdid in READS:
        LOGGER.info("input %s: %s", rdid, paths_input_trio[rdid])
    LOGGER.info("input samples: %s (%d samples for run \"%s\")", path_samples, len(samples), run_id)
    LOGGER.info("output dir: %s", dir_out)
    LOGGER.info("output counts: %s", path_counts)
    LOGGER.info("output details: %s", path_details)
    demux_by_barcode(samples, paths_input_trio, dir_out, path_counts, path_details, dry_run)

def demux_by_barcode(samples, fps, dir_out, path_counts, path_details, dry_run=False):
    """Demultiplex one trio of files from dictionaries of sample and barcode data.

    samples: dictionary of sample attributes for this run, with BarcodeFwdSeq
             and BarcodeRevSeq defined
    fps: dict of "I1", "R1", and "R2" keys pointing to file paths to
         fastq.gz inputs
    dir_out: output directory to write demultiplexed fastq.gz files to
    path_counts: path to csv to write per-sample read counts to.  If empty this
                 file isn't written.
    path_details: csv.gz file to write detailed per-read information.  If
                  empty, this file isn't written.
    dry_run: If True, don't actually call any commands or write any files.
    """
    counts = defaultdict(int)
    # nested dictionary of sample name -> trios of I1/R1/R2 paths
    # NOTE
    # with too many samples at once, this will cause an OS error due to too
    # many open files. In that case we'd have to open/close as needed.  It's
    # easy here to just open a bunch and store handles in a dictionary, though.
    fp_outs = {s: {rdid: Path(dir_out) / f"{s}.{rdid}.fastq.gz" for rdid in READS} for s in samples}
    fp_outs["None"] = {rdid: Path(dir_out) / f"unassigned.{rdid}.fastq.gz" for rdid in READS}
    for samp in fp_outs:
        LOGGER.info("output I1 for %s: %s", samp, fp_outs[samp]["I1"])
    # lookup table between pairs of barcodes and sample names
    bc_map = {(v["BarcodeFwdSeq"], v["BarcodeRevSeq"]): k for k, v in samples.items()}
    if not dry_run:
        Path(dir_out).mkdir(parents=True, exist_ok=True)
        try:
            f_outs = {s: {rdid: GZIP_OPEN(fp_outs[s][rdid], "wt") for rdid in READS} for s in fp_outs}
            details_writer = None
            if path_details:
                Path(path_details).parent.mkdir(parents=True, exist_ok=True)
                f_details = GZIP_OPEN(path_details, "wt")
                details_writer = DictWriter(
                    f_details,
                    fieldnames=[
                        "SeqID", "BarcodeFwdSeq", "BarcodeRevSeq",
                        "BarcodeFwdQualMin", "BarcodeRevQualMin"],
                    lineterminator="\n")
                details_writer.writeheader()
            with GZIP_OPEN(fps["I1"], "rt") as f_i1, \
                    GZIP_OPEN(fps["R1"], "rt") as f_r1, \
                    GZIP_OPEN(fps["R2"], "rt") as f_r2:
                for trio in zip(
                    # each of these is a tuple of (seqid, seq, qual) text
                    FastqGeneralIterator(f_i1),
                    FastqGeneralIterator(f_r1),
                    FastqGeneralIterator(f_r2)):
                    trio = list(trio)
                    trio.extend([
                        assign_barcode_fwd(trio[1][1], BARCODES_FWD),
                        assign_barcode_rev(trio[0][1], BARCODES_REV)])
                    _write_chunk([trio], bc_map, f_outs, counts, details_writer)
        finally:
            for trio in f_outs.values():
                for f_rd in trio.values():
                    f_rd.close()
            if path_details:
                f_details.close()

        if path_counts:
            _write_counts(path_counts, counts)

def _write_chunk(chunk, bc_map, f_outs, counts, details_writer):
    for entry in chunk:
        assigned = (entry[3], entry[4])
        recs = [list(entry[0]), list(entry[1]), list(entry[2])]
        sample_name = bc_map.get(assigned)
        counts[sample_name] += 1
        # trim barcode from forward read (both seq and qual)
        if sample_name:
            recs[1][1] = recs[1][1][len(assigned[0]):]
            recs[1][2] = recs[1][2][len(assigned[0]):]
        fmt = lambda rec: f"@{rec[0]}\n{rec[1]}\n+\n{rec[2]}\n"
        for idx, rdid in enumerate(READS):
            f_outs[str(sample_name)][rdid].write(fmt(recs[idx]))
        if details_writer:
            details_writer.writerow({
                "SeqID": entry[0][0].split(" ", 1)[0],
                "BarcodeFwdSeq": entry[3],
                "BarcodeRevSeq": entry[4],
                "BarcodeFwdQualMin": min(util.parse_quals(entry[1][2][:16])),
                "BarcodeRevQualMin": min(util.parse_quals(entry[0][2])),
                })

def _write_counts(path_counts, samp_counts):
    counts = []
    for key, val in samp_counts.items():
        row = {"Category": "demux", "NumSeqs": val}
        if key:
            row["Sample"] = key
            row["Item"] = "output"
        else:
            row["Item"] = "unassigned"
    util.save_counts(path_counts, counts)

def assign_barcode_fwd(seq, barcodes, max_mismatch=1):
    """Assign a forward barcode to a single R1 record.

    The forward barcodes come with a varying (known per barcode but not before
    we've assigned barcodes) prefix of at least four random nucleotides.

    seq: text of one R1 sequence
    barcodes: list of known forward barcodes, N included
    max_mismatch: how many mismatches to a barcode should be allowed?
    """
    mismatches = []
    for barcode in barcodes:
        prefix = len(barcode) - len(barcode.lstrip("N"))
        ref = barcode[prefix:]
        query = seq[prefix:len(barcode)]
        mismatch = sum([txt1 != txt2 for txt1, txt2 in zip(query, ref)])
        mismatches.append(mismatch)
    min_mismatch = min(mismatches)
    if min_mismatch <= max_mismatch:
        idx = mismatches.index(min_mismatch)
        return barcodes[idx]
    return None

def assign_barcode_rev(seq, barcodes, max_mismatch=1):
    """Assign a reverse barcode to a single I1 record.

    The reverse barcode is stored in the I1 read as the reverse-complement of
    what's listed in the protocol tables.

    seq: text of one I1 sequence
    barcodes: list of known reverse barcodes
    max_mismatch: how many mismatches to a barcode should be allowed?
    """
    # Reverse is simpler than forward, since we just need to take the
    # reverse-complement of the entire I1 read.
    query = Seq(seq).reverse_complement()
    mismatches = []
    for barcode in barcodes:
        mismatch = sum([txt1 != txt2 for txt1, txt2 in zip(query, barcode)])
        mismatches.append(mismatch)
    min_mismatch = min(mismatches)
    if min_mismatch <= max_mismatch:
        idx = mismatches.index(min_mismatch)
        return barcodes[idx]
    return None
