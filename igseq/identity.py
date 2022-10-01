"""
Calculate fraction of matching nucleotides between sequences.

Output is a long-format CSV with each row being one query/ref comparison.  If
no reference is given, the first sequence in the query is used as the
reference.  Supported formats are the same as for the convert command, with
some of the same column options for tabular input (though the same options will
be applied to both query and reference inputs for whichever are tabular).  This
allows, for example, junction sequences from an AIRR tsv as queries and a FASTA
of known junctions as references.

The scoring is based on a simple global pairwise alignment, with matches scored
as 1, mismatches and gaps 0.  Any existing gaps are removed before comparing
sequences, and differeces in case (lower/upper) are disregarded.
"""

import re
import logging
from Bio import Align
from .record import RecordReader, RecordWriter

LOGGER = logging.getLogger(__name__)
ALIGNER = Align.PairwiseAligner()

def identity(path_in, path_out, path_ref=None, fmt_in=None, fmt_in_ref=None, colmap=None, dry_run=False):
    """Create CSV of fraction of matching nucleotides between query and reference sequences."""
    LOGGER.info("given ref path: %s", path_ref)
    LOGGER.info("given query path: %s", path_in)
    LOGGER.info("given output path: %s", path_out)
    LOGGER.info("given input format: %s", fmt_in)
    LOGGER.info("given colmap: %s", colmap)
    refs = None
    fmt_ref = ""
    if path_ref:
        with RecordReader(path_ref, fmt_in_ref, colmap, dry_run=dry_run) as reader:
            refs = list(reader)
            fmt_ref = reader.fmt
    else:
        if fmt_in_ref:
            LOGGER.warning("given ref format ignored if ref not given")
    # Trying to leverage the same (otherwise-sequence-focused) Record code here
    # for a more generic table.  Will implicitly also allow csv.gz/tsv/tsv.gz
    fmt_out = None
    if path_out == "-":
        fmt_out = "csv"
    colmap_out = {field: field for field in ["query", "ref", "identity"]}
    tabular = lambda txt: any(ext in txt for ext in ["csv", "tsv"])
    with RecordReader(path_in, fmt_in, colmap, dry_run=dry_run) as reader, \
        RecordWriter(path_out, fmt_out, colmap_out, dry_run=dry_run) as writer:
        if not fmt_ref:
            fmt_ref = reader.fmt
        if not tabular(reader.fmt) and not tabular(fmt_ref) and colmap:
            LOGGER.warning("given column names ignored for non-tabular inputs")
        for record in reader:
            if refs is None:
                refs = [record]
            for ref in refs:
                score = calc_identity(align(
                    record[reader.colmap["sequence"]],
                    ref[reader.colmap["sequence"]]))
                writer.write({
                    "query": record["sequence_id"],
                    "ref": ref["sequence_id"],
                    "identity": score})

def align(seq1, seq2):
    """Simple pairwise alignment, ignoring case and gaps in input.

    If either or both inputs are empty, the result is None.
    """
    if seq1 and seq2:
        seq1 = _munge_seq(seq1)
        seq2 = _munge_seq(seq2)
        return ALIGNER.align(seq1, seq2)[0]
    return None

def calc_identity(aln):
    """Calculate identity fraction from an alignment object."""
    if not aln:
        return 0
    return aln.score/max(len(aln.target), len(aln.query))

def calc_coverage(aln):
    """Calculate coverage fraction from an alignment object."""
    # See Biopython #3991
    if not aln:
        return 0
    parts = aln.format().split("\n")
    ref, query = [parts[0], parts[2]]
    # trimming, based on SONAR
    left_gap = re.match( "-+", ref) or re.match( "-+", query)
    if left_gap:
        query = query[left_gap.end():]
        ref = ref[left_gap.end():]
    # NOTE: not quite the same as SONAR which has re.match("-+", query).  I
    # don't see how that would work though.
    right_gap = re.search("-+$", ref) or re.search("-+$", query)
    if right_gap:
        query = query[:right_gap.start()]
        ref = ref[:right_gap.start()]
    coverage = len(re.sub("-", "", ref))/len(aln.target)
    return coverage

def _munge_seq(seq):
    # We'll disregard gaps and case.  If whatever these objects are can't
    # handle those methods we'll throw a ValueError.
    try:
        return seq.replace("-", "").upper()
    except AttributeError as err:
        raise ValueError("argument should be Seq object or string") from err
