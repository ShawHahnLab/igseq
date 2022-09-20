"""
Create and format phylogenetic trees.

Takes input as aligned or unaligned sequences (same formats as the convert
command supports) or an existing tree file and creates tree files in various
formats.

Tree formats:

    newick:  Plain Newick text, like (:0.1,:0.2,(...
    nex:     NEXUS file with Newick+styling info
    png/pdf: Rendered tree image (for output)

Nodes and branches can be color-coded according to each node's membership in
one or more sets.
"""

import logging
from pathlib import Path
from subprocess import Popen, PIPE
from . import util
from . import record

LOGGER = logging.getLogger(__name__)

FMT_EXT_MAP = {
    ".nex": "nex",
    ".nxs": "nex",
    ".nexus": "nex",
    ".tree": "newick",
    ".newick": "newick",
    ".pdf": "pdf",
    ".png": "png",
    }

# supported for input and output
FMT_IN = {"nex", "newick"}
FMT_OUT = {"nex", "newick", "pdf", "png"}

def tree(path_in, path_out, fmt_in=None, fmt_out=None, aligned=None, pattern=None, lists=None, colmap=None, dry_run=False):
    LOGGER.info("given input path: %s", path_in)
    LOGGER.info("given output path: %s", path_out)
    LOGGER.info("given input format: %s", fmt_in)
    LOGGER.info("given output format: %s", fmt_out)
    LOGGER.info("given aligned attribute: %s", aligned)
    LOGGER.info("given colmap: %s", colmap)
    LOGGER.info("given set pattern: %s", pattern)
    LOGGER.info("given set lists: %s", lists)

    # Infer input format
    fmt_in_tree = infer_tree_fmt(path_in)
    fmt_in_seqs = record.RecordHandler._infer_fmt(path_in)
    fmt_in_inf = fmt_in or fmt_in_tree or fmt_in_seqs
    if not fmt_in_inf:
        raise util.IgSeqError(
            f"input format not detected from filename ({path_in}).  "
            "specify a format manually.")
    if not fmt_in:
        LOGGER.info("inferred input format: %s", fmt_in_inf)
    if fmt_in_inf not in FMT_IN and fmt_in_inf not in record.FMT_EXT_MAP.values():
        raise util.IgSeqError(f"Can't use format {fmt_in_inf} for input")

    # Infer output format
    fmt_out_inf = fmt_out or infer_tree_fmt(path_out)
    if not fmt_out:
        LOGGER.info("inferred output format: %s", fmt_out_inf)
    if fmt_out_inf not in FMT_OUT:
        raise util.IgSeqError(f"Can't use format {fmt_out_inf} for output")

    # Handle input
    if fmt_in_inf in record.FMT_EXT_MAP.values():
        # Sequence input
        with record.RecordReader(path_in, fmt_in_inf, colmap, dry_run=dry_run) as reader:
            records = list(reader)
        # TODO handle empty case
        if aligned is None:
            aligned = looks_aligned(records)
            LOGGER.info("inferred aligned attribute: %s", aligned)
        if not aligned:
            # TODO: align!
            raise NotImplementedError("MSA not yet implemented")
        # TODO make tree
        newick_data = run_fasttree(records)
    else:
        # tree input
        if fmt_in_inf == "newick":
            newick_data = load_newick(path_in)
        else:
            raise NotImplementedError("NEXUS input not yet implemented")

    # Handle output
    if fmt_out_inf == "newick":
        with open(path_out, "wt") as f_out:
            f_out.write(newick_data)
    elif fmt_out_inf == "nex":
        raise NotImplementedError("NEXUS output not yet implemented")
    else:
        raise NotImplementedError("image output not yet implemented")

def looks_aligned(records):
    # all the same length looks aligned, multiple lengths doesn't.
    lengths = {len(rec["sequence"]) for rec in records}
    return len(lengths) == 1

def infer_tree_fmt(path):
    try:
        path = Path(path)
    except TypeError:
        path = Path(str(path.name))
    ext = path.suffix.lower()
    return FMT_EXT_MAP.get(ext)

def run_fasttree(records):
    # TODO -quiet unless -v was given via CLI
    args = ["fasttree", "-nt", "-quiet"]
    with Popen(args, stdin=PIPE, stdout=PIPE, text=True) as proc:
        with record.RecordWriter(proc.stdin, "fa") as writer:
            for rec in records:
                writer.write(rec)
            writer.close()
        proc.stdin.close()
        treedata = proc.stdout.read()
    return treedata

def load_newick(path):
    with open(path) as f_in:
        return f_in.read()
