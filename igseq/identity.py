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
sequences, and differences in case (lower/upper) are disregarded.
"""

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
                score = score_identity(
                    record[reader.colmap["sequence"]],
                    ref[reader.colmap["sequence"]])
                writer.write({
                    "query": record[reader.colmap["sequence_id"]],
                    "ref": ref["sequence_id"],
                    "identity": score})

def score_identity(seq1, seq2):
    """Get identity fraction between two sequences.

    If one or both are zero-length the result will be 0.
    """
    # These can be strings or Seq objects, but not SeqRecords.
    if seq1 and seq2:
        # We'll disregard gaps and case.  If whatever these objects are can't
        # handle those methods we'll throw a ValueError.
        try:
            seq1 = seq1.replace("-", "").upper()
            seq2 = seq2.replace("-", "").upper()
        except AttributeError as err:
            raise ValueError("score_identity arguments should be Seq objects or strings") from err
        return ALIGNER.score(seq1, seq2)/max(len(seq1), len(seq2))
    return 0
