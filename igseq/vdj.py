"""
VDJ-handling helpers used elsewhere.

There are no user-facing commands in here.
"""

import re
import logging
from os import PathLike
from pathlib import Path
from Bio import SeqIO
from . import util

LOGGER = logging.getLogger(__name__)

LOCI = ["IGH", "IGK", "IGL"]
SEGMENTS = ["V", "D", "J"]
LOCUS_SEGMENTS = ["IGHV", "IGHD", "IGHJ", "IGKV", "IGKJ", "IGLV", "IGLJ"]

def parse_vdj_paths(ref_paths):
    """Take a list of file/directory/builtin names and return dicts with more info.

    For example, "some/dir/V.fasta" would be recognized as a single FASTA
    input, "some/dir" split into separate V/D/J FASTA found inside, and
    "rhesus" split into the built-in reference FASTAs provided by the package.
    """

    if isinstance(ref_paths, (str, PathLike)):
        ref_paths = [ref_paths]
    attrs_list = []
    for entry in ref_paths:
        internal_matches = get_internal_vdj(entry)
        path = Path(entry)
        # first priority: actual file name
        if path.is_file():
            attrs = parse_vdj_filename(path)
            attrs["input"] = str(entry)
            attrs["type"] = "file"
            if "segment" not in attrs:
                raise ValueError("couldn't determine segment for file: %s" % path)
            attrs_list.append(attrs)
        # second priority: actual directory name
        if path.is_dir():
            for path2 in sorted(path.glob("*")):
                attrs = parse_vdj_filename(path2)
                attrs["input"] = str(entry)
                if attrs["fasta"]:
                    if "segment" not in attrs:
                        raise ValueError("couldn't determine segment for file: %s" % path)
                    attrs["type"] = "dir"
                    attrs_list.append(attrs)
        # third priority: internal reference
        elif internal_matches:
            for fasta in internal_matches:
                fasta_rel = fasta.relative_to(util.DATA / "germ")
                attrs = parse_vdj_filename(fasta_rel)
                attrs["input"] = str(entry)
                attrs["path"] = fasta
                attrs["type"] = "internal"
                parents = [parent.name for parent in fasta_rel.parents if parent.name]
                attrs["species"] = parents[-1]
                attrs["ref"] = parents[-2]
                attrs_list.append(attrs)
        else:
            raise util.IgSeqError("ref path not recognized: %s" % entry)
    return attrs_list

def get_internal_vdj(name):
    """Get list of builtin germline FASTA files matching a path fragment."""
    name = str(name)
    germ = util.DATA / "germ"
    candidates = list(germ.glob("**/*.fasta"))
    output = []
    for candidate in candidates:
        if name in str(candidate.relative_to(germ)):
            output.append(candidate)
    return sorted(output)

def parse_vdj_filename(txt):
    """Parse recognized VDJ fields in a filename into a dictionary.

    For example, "IGH/V.fasta" ->
    {"locus": "IGH", "segment": "V", "fasta": True}
    """
    # "IGHV.fasta" -> {"path": "IGHV.fasta", "locus": "IGH", "segment": "V"}
    txt = str(txt)
    attrs = {"path": Path(txt)}
    prefix = attrs["path"].stem
    suffix = attrs["path"].suffix.lower()
    attrs["fasta"] = suffix in [".fasta", ".fa", ".fna"]
    parents = [parent.name for parent in Path(txt).parents]
    attrs.update(_parse_vdj_tokens(parents))
    name_fields = re.split("[^A-Za-z]+", prefix)
    attrs.update(_parse_vdj_tokens(name_fields))
    return attrs

def combine_vdj(attrs_list, fasta):
    """Combine FASTAs from a list from parse_vdj_paths into one file."""
    fasta = Path(fasta)
    fasta.parent.mkdir(parents=True, exist_ok=True)
    with open(fasta, "wt") as f_out:
        by_locus = {"IGH": [], "IGK": [], "IGL": [], None: []}
        if len(attrs_list) > 1:
            for attrs in attrs_list:
                by_locus[attrs.get("locus")].append(attrs)
            if by_locus[None]:
                add_suffix = {key: True for key in by_locus}
            else:
                add_suffix = {key: len(val) > 1 for key, val in by_locus.items()}
        else:
            add_suffix = None
        for attrs in attrs_list:
            add_suffix_here = add_suffix and add_suffix[attrs.get("locus")]
            with open(attrs["path"]) as f_in:
                for record in SeqIO.parse(f_in, "fasta"):
                    if add_suffix_here:
                        if attrs["type"] == "internal":
                            suffix = attrs["species"] + "_" + attrs["ref"]
                        else:
                            suffix = str(attrs["path"])
                        record.id += "_" + suffix
                        record.description=""
                    SeqIO.write(record, f_out, "fasta-2line")

def group(attrs_list, keyfunc=None):
    """Group list items from parse_vdj_paths in a dictionary.

    The keys are the segments from each item, and optionally an outer grouping
    using keys from the given keyfunc.
    """
    groups = {}
    for entry in attrs_list:
        if keyfunc:
            key = keyfunc(entry)
        else:
            key = "all"
        if key not in groups:
            groups[key] = {seg: [] for seg in SEGMENTS}
        LOGGER.debug(f"grouping {entry['path']} as \"{key}\"")
        groups[key][entry["segment"]].append(entry)
    if not keyfunc:
        return groups["all"]
    return groups

def _parse_vdj_tokens(tokens):
    # ["IGH", "blah", "V"] -> {"locus": "IGH", "segment": "V"}
    attrs = {}
    for token in tokens:
        token = token.upper()
        if token in LOCUS_SEGMENTS:
            attrs["locus"] = token[0:3]
            attrs["segment"] = token[3]
            break
        if not attrs.get("locus") and token in LOCI:
            attrs["locus"] = token
        if not attrs.get("segment") and token in SEGMENTS:
            attrs["segment"] = token
    return attrs
