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
one or more sets, as defined by a regular expression matching sequence IDs
(-P) or explicit lists of sequence IDs via text files (-L).  Note that set
membership is combined from all sources; seq1 matching set1 via the
pattern-matching can also be listed via -L as a member of set2.
"""

import re
import logging
import random
from pathlib import Path
from subprocess import Popen, PIPE
from collections import defaultdict
from . import util
from . import record
from . import msa

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

def tree(path_in, path_out, fmt_in=None, fmt_out=None, aligned=None, pattern=None, lists=None, colors=None, colmap=None, dry_run=False):
    LOGGER.info("given input path: %s", path_in)
    LOGGER.info("given output path: %s", path_out)
    LOGGER.info("given input format: %s", fmt_in)
    LOGGER.info("given output format: %s", fmt_out)
    LOGGER.info("given aligned attribute: %s", aligned)
    LOGGER.info("given set pattern: %s", pattern)
    LOGGER.info("given set lists: %s", lists)
    LOGGER.info("given set colors: %s", colors)
    LOGGER.info("given colmap: %s", colmap)

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

    # Parse lists
    lists = parse_lists(lists)
    if lists:
        LOGGER.info("parsed lists: %s", lists)

    # Parse colors
    colors = parse_colors(colors)
    if colors:
        LOGGER.info("parsed colors: %s", colors)

    # Handle input
    if fmt_in_inf in record.FMT_EXT_MAP.values():
        # Sequence input
        with record.RecordReader(path_in, fmt_in_inf, colmap, dry_run=dry_run) as reader:
            records = list(reader)
        # I've written run_muscle so it'll just pass through an empty set of
        # records with a warning, but even if we did that, fasttree wouldn't
        # know what to do with it.  better to halt and catch fire.
        if not records:
            raise util.IgSeqError(f"No sequences provided from {path_in}")
        if aligned is None:
            aligned = looks_aligned(records)
            LOGGER.info("inferred aligned attribute: %s", aligned)
        if not aligned:
            LOGGER.info("aligning records")
            records = msa.run_muscle(records)
        newick_text = run_fasttree(records)
    else:
        # tree input
        if fmt_in_inf == "newick":
            raise NotImplementedError("Newick input not yet implemented")
            newick_text = load_newick_text(path_in)
        else:
            raise NotImplementedError("NEXUS input not yet implemented")

    # Handle output
    if fmt_out_inf == "newick":
        with open(path_out, "wt") as f_out:
            f_out.write(newick_text)
    elif fmt_out_inf == "nex":
        seq_ids = [rec["sequence_id"] for rec in records]
        seq_sets = build_seq_sets(seq_ids, pattern, lists)
        seq_colors = color_seqs(seq_ids, seq_sets, colors)
        save_nexus(seq_colors, newick_text, path_out)
    else:
        raise NotImplementedError("image output not yet implemented")

def parse_lists(lists):
    lists = lists or []
    # TODO allow missing name
    lists = [seq_list.split("=", 1) for seq_list in lists]
    lists = {p[0]: load_seq_list(p[1]) for p in lists}
    return lists

def parse_colors(color_texts):
    """Parse list of "setname=colorcode" strings into dictionary of RGB tuples."""
    colors = color_texts or []
    # TODO allow missing name
    colors = [color.split("=", 1) for color in colors]
    colors = {p[0]: color_str_to_trio(p[1]) for p in colors}
    return colors

def load_seq_list(path):
    seq_list = []
    with open(path) as f_in:
        for line in f_in:
            seq_list.append(line.strip())
    return seq_list

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
    if not records:
        raise util.IgSeqError("No seq records provided")
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

def load_newick_text(path):
    with open(path) as f_in:
        return f_in.read()

def build_seq_sets(seq_ids, pattern=None, lists=None):
    # mapping of set names to sets of seq IDs
    seq_sets = defaultdict(set)
    if pattern:
        for seq_id in seq_ids:
            match = re.search(pattern, seq_id)
            if match:
                try:
                    set_name = match.group(1)
                except IndexError:
                    set_name = match.group(0)
                seq_sets[set_name].add(seq_id)
    if lists:
        for set_name, seq_set_ids in lists.items():
            seq_sets[set_name].update(seq_set_ids)
    return seq_sets

def color_seqs(seq_ids, seq_sets, seq_set_colors=None):
    """Assign colors to sequences based on set membership."""
    seq_set_colors_combo = make_seq_set_colors(seq_sets)
    if seq_set_colors:
        seq_set_colors_combo.update(seq_set_colors)
    seq_colors = {}
    for seq_id in seq_ids:
        # sets that contain this sequence, and the colors for those sets
        sets_here = {set_name for set_name in seq_sets if seq_id in seq_sets[set_name]}
        colors_here = [seq_set_colors_combo[set_name] for set_name in sets_here]
        combo_color = merge_colors(colors_here, len(seq_set_colors_combo))
        seq_colors[seq_id] = combo_color
    return seq_colors

def make_seq_set_colors(seq_sets):
    # TODO do smarter than just randomizing!
    seq_set_colors = {}
    for set_name in seq_sets:
        seq_set_colors[set_name] = [random.randint(0, 255) for _ in range(3)]
    return seq_set_colors

def merge_colors(colors, scale):
    """Take an average of a list of colors and shift toward black."""
    result = [0, 0, 0]
    if not colors:
        return result
    if len(colors) == 1:
        return colors[0]
    for color in colors:
        for idx in range(3):
            result[idx] += color[idx]
    # not quite right, should rotate, really, not move directly toward the
    # middle... but it'll do for now
    scaling = ((scale - len(colors))/scale)**0.3
    for idx in range(3):
        result[idx] = result[idx] / len(colors)
        result[idx] = int(result[idx] * scaling)
    return result

def color_str_to_trio(color_txt):
    """Convert hex color string to trio of 0:255 ints."""
    color_txt = color_txt.removeprefix("#")
    color = [int(color_txt[idx:(idx+2)], 16) for idx in range(0, 6, 2)]
    return color

def color_trio_to_str(color):
    """Convert trio of 0:255 ints to hex color string."""
    return "#" + "".join([f"{c:02x}" for c in color])

def save_nexus(seq_colors, newick_text, path):
    with open(path, "wt") as f_out:
        f_out.write("#NEXUS\n")
        f_out.write("begin taxa;\n")
        f_out.write(f"dimensions ntax={len(seq_colors)};\n")
        f_out.write("taxlabels\n")
        for seq_id, seq_color in seq_colors.items():
            color_txt = color_trio_to_str(seq_color)
            f_out.write(f"'{seq_id}'[&!color={color_txt}]\n")
        f_out.write(";\n")
        f_out.write("end;\n")
        f_out.write("\n")
        # TODO mark rooted/unrooted?
        f_out.write("begin trees;\n")
        f_out.write(f"tree tree_1 = {newick_text}\n")
        f_out.write("end;\n")
        f_out.write("\n")
