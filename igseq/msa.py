"""
Create multiple sequence alignments.

Supported formats are the same as for the convert command, with some of the
same column options for tabular input (though the same options will be applied
to both query and reference inputs for whichever are tabular).

Currently the aligner used is MUSCLE.
"""

import sys
import logging
from subprocess import Popen, PIPE
from io import StringIO
from .util import IgSeqError
from .record import RecordReader, RecordWriter

LOGGER = logging.getLogger(__name__)

def msa(path_in, path_out, fmt_in=None, fmt_out=None, colmap=None, dry_run=False):
    """Align sequences into a multiple sequence alignment."""
    LOGGER.info("given input path: %s", path_in)
    LOGGER.info("given output path: %s", path_out)
    LOGGER.info("given input format: %s", fmt_in)
    LOGGER.info("given output format: %s", fmt_out)
    LOGGER.info("given colmap: %s", colmap)
    with RecordReader(path_in, fmt_in, colmap, dry_run=dry_run) as reader:
        records = list(reader)
    if not dry_run:
        LOGGER.info("aligning sequences")
        records_aln = run_muscle(records)
        LOGGER.info("writing output")
        with RecordWriter(path_out, fmt_out, colmap) as writer:
            for record in records_aln:
                writer.write(record)

def run_muscle(records):
    """Align a set of records with MUSCLE."""
    args = ["muscle", "-align", "/dev/stdin", "-output", "/dev/stdout"]
    with Popen(args, stdin=PIPE, stdout=PIPE, stderr=PIPE, text=True) as proc:
        with RecordWriter(proc.stdin, "fa") as writer:
            for rec in records:
                writer.write(rec)
        proc.stdin.close()
        msa_txt = proc.stdout.read()
        stderr = proc.stderr.read()
    sys.stderr.write(stderr)
    if proc.returncode:
        LOGGER.critical("%s exited with code %d", proc.args[0], proc.returncode)
        raise IgSeqError("MUSCLE crashed")
    with RecordReader(StringIO(msa_txt), "fa") as reader:
        records_out = list(reader)
    # sort according to input records
    record_lut = {rec["sequence_id"]: rec for rec in records_out}
    seqids = [rec["sequence_id"] for rec in records]
    records_out = [record_lut[seqid] for seqid in seqids]
    return records_out
