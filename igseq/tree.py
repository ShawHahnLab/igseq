"""
Create and format phylogenetic trees.

Takes input as aligned or unaligned sequences (same formats as the convert
command supports) or an existing tree file and creates tree files in various
formats.

Tree formats:

    newick:  Plain Newick text, like (:0.1,:0.2,(...
    nex:     NEXUS file with Newick+styling info

Multiple input files are allowed if all are FASTA.  In this case sequences will
implicitly be grouped into sets based on which sequences are present in which
input files.  Sequences with the same ID must share identical sequence content,
and gaps are removed and the sequences are then aligned.

Nodes and branches can be color-coded according to each node's membership in
one or more sets, as defined by a regular expression matching sequence IDs
(-P) or explicit lists of sequence IDs via text files (-L).  Note that set
membership is combined from all sources; seq1 matching set1 via the
pattern-matching can also be listed via -L as a member of set2.

Options defined per-set can optionally specify a set name, for example -L
set1=path/list.txt or -C set1=#ff0000.  Unnamed sets will be named "set#" in
the order they are encountered.
"""

import re
import logging
import random
from pathlib import Path, PurePath
from subprocess import Popen, PIPE
from collections import defaultdict
import newick
from Bio.Nexus.Nexus import Nexus
from . import util
from . import record
from . import msa
from . import colors as colortools

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

def tree(paths_in, path_out,
        fmt_in=None, fmt_out=None, aligned=None, pattern=None, lists=None,
        colors=None, merge_colors=False, figtree_opts=None, colmap=None, dry_run=False):
    LOGGER.info("given input path(s): %s", paths_in)
    LOGGER.info("given output path: %s", path_out)
    LOGGER.info("given input format: %s", fmt_in)
    LOGGER.info("given output format: %s", fmt_out)
    LOGGER.info("given aligned attribute: %s", aligned)
    LOGGER.info("given set pattern: %s", pattern)
    LOGGER.info("given set lists: %s", lists)
    LOGGER.info("given set colors: %s", colors)
    LOGGER.info("given merge colors setting: %s", merge_colors)
    LOGGER.info("given figtree options: %s", figtree_opts)
    LOGGER.info("given colmap: %s", colmap)

    paths_in_parsed = parse_paths_in(paths_in)
    fmt_in_inf = _infer_input_format(paths_in_parsed.values(), fmt_in)
    fmt_out_inf = _infer_output_format(path_out, fmt_out)
    seq_sets = parse_lists(lists)
    colors = parse_colors(colors)

    if fmt_out_inf != "nex" and figtree_opts:
        LOGGER.warning(
            "Given FigTree options are ignored for non-nex output format %s", fmt_out_inf)

    # Handle input
    if fmt_in_inf in record.FMT_EXT_MAP.values():
        # Sequence input.  Only include one sequence per ID because fasttree
        # doesn't allow duplicates.
        seen = {}
        records = []
        if len(paths_in_parsed) > 1 and aligned is not False:
            # If the input is multiple FASTA files, require alignment, and
            # disregard any existing gaps while reading in seqs below.
            LOGGER.info("requiring alignment for multi-FASTA input")
            aligned = False
        for set_name, path_in in paths_in_parsed.items():
            if set_name not in seq_sets and len(paths_in_parsed) > 1:
                seq_sets[set_name] = set()
            with record.RecordReader(path_in, fmt_in_inf, colmap, dry_run=dry_run) as reader:
                for rec in reader:
                    seq = rec["sequence"].upper()
                    if len(paths_in_parsed) > 1:
                        seq = re.sub("-", "", seq)
                        seq_sets[set_name].add(rec["sequence_id"])
                    if rec["sequence_id"] in seen:
                        # if we've seen this seq ID before, is the sequence
                        # the same?  If not complain loudly
                        if seq != seen[rec["sequence_id"]]:
                            raise util.IgSeqError(
                                "Duplicate sequence IDs encountered " \
                                f"for different sequences from {paths_in}")
                    else:
                        seen[rec["sequence_id"]] = seq
                        rec["sequence"] = seq
                        records.append(rec)
        # I've written run_muscle so it'll just pass through an empty set of
        # records with a warning, but even if we did that, fasttree wouldn't
        # know what to do with it.  better to halt and catch fire.
        if not records:
            raise util.IgSeqError(f"No sequences provided from {paths_in}")
        if aligned is None:
            aligned = looks_aligned(records)
            LOGGER.info("inferred aligned attribute: %s", aligned)
        if not aligned:
            LOGGER.info("aligning records")
            records = msa.run_muscle(records)
        newick_text = run_fasttree(records)
        newick_obj = newick.loads(newick_text)[0]
    else:
        path_in = list(paths_in_parsed.values())[0]
        # tree input
        if fmt_in_inf == "newick":
            newick_obj = newick.read(path_in)[0]
        else:
            newick_obj = load_newick_from_nexus(path_in)
            if not newick_obj:
                raise util.IgSeqError("No tree found in NEXUS file {path_in}")

    # Handle output
    if fmt_out_inf == "newick":
        with open(path_out, "wt", encoding="ASCII") as f_out:
            newick.dump(newick_obj, f_out)
            f_out.write("\n")
    elif fmt_out_inf == "nex":
        seq_ids = all_leaf_ids(newick_obj)
        seq_sets_combo = build_seq_sets(seq_ids, pattern, seq_sets)
        seq_colors = color_seqs(seq_ids, seq_sets_combo, merge_colors, colors)
        custom_blocks = None
        if figtree_opts:
            custom_blocks = {"figtree": [f"set {opt}" for opt in figtree_opts]}
        save_nexus(newick_obj, path_out, seq_colors, seq_sets_combo, custom_blocks)
    else:
        raise NotImplementedError("image output not yet implemented")

def _infer_input_format(paths_in, fmt_in):
    # multiple paths only supported for FASTA
    fmt_in_inf_set = set()
    for path_in in paths_in:
        fmt_in_tree = infer_tree_fmt(path_in)
        fmt_in_seqs = record.RecordHandler._infer_fmt(path_in)
        fmt_in_inf = fmt_in or fmt_in_tree or fmt_in_seqs
        if not fmt_in_inf:
            raise util.IgSeqError(
                f"input format not detected from filename ({path_in}).  "
                "specify a format manually.")
        fmt_in_inf_set.add(fmt_in_inf)
    if len(fmt_in_inf_set) > 1:
        raise util.IgSeqError("Multiple inputs only supported with FASTA")
    fmt_in_inf = fmt_in_inf_set.pop()
    if not fmt_in:
        LOGGER.info("inferred input format: %s", fmt_in_inf)
    if fmt_in_inf not in FMT_IN and fmt_in_inf not in record.FMT_EXT_MAP.values():
        raise util.IgSeqError(f"Can't use format {fmt_in_inf} for input")
    return fmt_in_inf

def _infer_output_format(path_out, fmt_out):
    fmt_out_inf = fmt_out or infer_tree_fmt(path_out)
    if not fmt_out:
        LOGGER.info("inferred output format: %s", fmt_out_inf)
    if fmt_out_inf not in FMT_OUT:
        if not fmt_out_inf:
            raise util.IgSeqError(f"Unrecognized output format for {path_out}")
        raise util.IgSeqError(f"Can't use format {fmt_out_inf} for output")
    return fmt_out_inf

def all_leaf_ids(tree_obj):
    """Get a set of all non-empty leaf node IDs in a tree object."""
    node_ids = []
    for child in tree_obj.descendants:
        for leaf_id in all_leaf_ids(child):
            if leaf_id not in node_ids:
                node_ids.append(leaf_id)
    if tree_obj.name and tree_obj.is_leaf:
        node_ids.append(tree_obj.name)
    return node_ids

def parse_paths_in(paths_in):
    if isinstance(paths_in, (str, PurePath)):
        paths_in = [paths_in]
    paths = [str(named_path).split("=", 1) for named_path in paths_in]
    def parse_named_path(idx, pair):
        try:
            key, path = pair
        except ValueError:
            key = f"set{idx+1}"
            path = pair[0]
        return key, path
    named_paths = dict(parse_named_path(idx, p) for idx, p in enumerate(paths))
    return named_paths

def parse_lists(lists):
    """Parse a list of filenames into a dictionary of sequence sets."""
    lists = lists or []
    lists = [str(seq_list).split("=", 1) for seq_list in lists]
    def load_seq_list(idx, pair):
        try:
            key, path = pair
        except ValueError:
            key = f"set{idx+1}"
            path = pair[0]
        seq_set = set()
        with open(path) as f_in:
            for line in f_in:
                seq_set.add(line.strip())
        return key, seq_set
    seq_sets = dict(load_seq_list(idx, p) for idx, p in enumerate(lists))
    if seq_sets:
        LOGGER.info("parsed lists: %s", seq_sets)
    return seq_sets

def parse_colors(color_texts):
    """Parse list of "setname=colorcode" strings into dictionary of RGB tuples."""
    color_texts = color_texts or []
    color_texts = [str(txt).split("=", 1) for txt in color_texts]
    def parse_color(idx, pair):
        try:
            key, colortxt = pair
        except ValueError:
            key = f"set{idx+1}"
            colortxt = pair[0]
        color = colortools.color_str_to_trio(colortxt)
        return key, color
    colors_out = dict(parse_color(idx, p) for idx, p in enumerate(color_texts))
    if colors_out:
        LOGGER.info("parsed colors: %s", colors_out)
    return colors_out

def looks_aligned(records):
    """True if one or more records are all the same length, False otherwise."""
    lengths = {len(rec["sequence"]) for rec in records}
    return len(lengths) == 1

def infer_tree_fmt(path):
    try:
        path = Path(path)
    except TypeError:
        path = Path(str(path.name))
    ext = path.suffix.lower()
    return FMT_EXT_MAP.get(ext)

def load_newick_from_nexus(path_in):
    nexus_obj = Nexus()
    nexus_obj.read(path_in)
    for obj in nexus_obj.structured:
        if obj.title == "trees":
            for cmd in obj.commandlines:
                if cmd.command == "tree":
                    newick_text = re.sub("^[^=]+= *", "", cmd.options)
                    newick_obj = newick.loads(newick_text)[0]
                    return newick_obj
    return None

def run_fasttree(records):
    """Run fasttree on a list of sequence records and return newick text"""
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

def build_seq_sets(seq_ids, pattern=None, lists=None):
    """Given a list of seq IDs, place each in one or more sets."""
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
    # Sort by key so we get a consistent order (yes, dicts are ordered these
    # days)
    seq_sets = {key: seq_sets[key] for key in sorted(seq_sets.keys())}
    return seq_sets

def color_seqs(seq_ids, seq_sets, merge_colors=False, seq_set_colors=None):
    """Assign colors to sequences based on set membership."""
    seq_set_colors_combo = make_seq_set_colors(seq_sets)
    if seq_set_colors:
        seq_set_colors_combo.update(seq_set_colors)
    seq_colors = {}
    for seq_id in seq_ids:
        # sets that contain this sequence, and the colors for those sets
        sets_here = [set_name for set_name in seq_sets if seq_id in seq_sets[set_name]]
        colors_here = [seq_set_colors_combo[set_name] for set_name in sets_here]
        if merge_colors:
            combo_color = colortools.merge_colors(colors_here, len(seq_set_colors_combo))
        elif colors_here:
            combo_color = colors_here[-1]
        else:
            combo_color = [0, 0, 0]
        LOGGER.debug(
            "color_seqs: %s belongs to %s with colors %s -> %s",
            seq_id, sets_here, colors_here, combo_color)
        seq_colors[seq_id] = combo_color
    return seq_colors

def make_seq_set_colors(seq_sets):
    seq_set_colors = {}
    for idx, set_name in enumerate(seq_sets):
        # adapted from SONAR
        # this stretches across COLORS in even increments for as many as we need here
        num = len(colortools.COLORS)
        subset = [int( a * (num-1) / max(1, (len(seq_sets)-1)) ) for a in range(num)]
        try:
            seq_set_colors[set_name] = colortools.color_str_to_trio(colortools.COLORS[subset[idx]])
        except IndexError:
            seq_set_colors[set_name] = [random.randint(0, 255) for _ in range(3)]
    return seq_set_colors

def save_nexus(newick_obj, path, seq_colors, seq_sets=None, custom_blocks=None):
    newick_text = newick.dumps(newick_obj)
    with open(path, "wt", encoding="ASCII") as f_out:
        f_out.write("#NEXUS\n")
        if seq_sets:
            # https://plewis.github.io/nexus/#sets-block
            f_out.write("begin sets;\n")
            for set_key, set_seqids in seq_sets.items():
                set_seqids_txt = " ".join(sorted(set_seqids))
                f_out.write(f"taxset {set_key} = {set_seqids_txt}\n")
            f_out.write("end;\n")
        f_out.write("begin taxa;\n")
        f_out.write(f"dimensions ntax={len(seq_colors)};\n")
        f_out.write("taxlabels\n")
        seq_color_pairs = sorted([[k, v] for k, v in seq_colors.items()])
        for seq_id, seq_color in seq_color_pairs:
            color_txt = colortools.color_trio_to_str(seq_color)
            f_out.write(f"'{seq_id}'[&!color={color_txt}]\n")
        f_out.write(";\n")
        f_out.write("end;\n")
        f_out.write("\n")
        # TODO mark rooted/unrooted?
        f_out.write("begin trees;\n")
        f_out.write(f"tree tree_1 = {newick_text}\n\n")
        f_out.write("end;\n")
        f_out.write("\n")
        if custom_blocks:
            for block_name, block_lines in custom_blocks.items():
                f_out.write(f"begin {block_name};\n")
                for line in block_lines:
                    f_out.write(f"\t{line};\n")
                f_out.write("end;\n")
