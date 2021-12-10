"""
VDJ-handling helpers used elsewhere.

There are no user-facing commands in here.
"""

import re
import csv
import hashlib
import logging
from os import PathLike
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
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
        elif path.is_dir():
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
    # handle duplicates produced from multiple ref_paths leading to the same
    # files
    groups = {}
    for attrs in attrs_list:
        if attrs["path"] not in groups:
            groups[attrs["path"]] = attrs
        else:
            new_input = groups[attrs["path"]]["input"] + "; " + str(attrs["input"])
            groups[attrs["path"]].update(attrs)
            groups[attrs["path"]]["input"] = new_input
    attrs_list = list(groups.values())
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

def combine_vdj_by_attrs(attrs_list, fasta):
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

def combine_vdj(fastas_in, dir_out, csv_lookup_table=None, dry_run=False):
    """Combine FASTAs from a list into a trio of V/D/J files.

    fastas_in: list of paths to individual FASTA files
    csv_lookup_table: path to output a per-sequence lookup table.  When
                      multiple files are given with entries for the same locus
                      and segment this maps modified sequence IDs used here to
                      the original ID and input file.
    dry_run: If True, don't actually write any files
    """
    # totally different approach from the other functions here: parse the
    # details from each sequence ID.
    dir_out = Path(dir_out)
    if not csv_lookup_table:
        csv_lookup_table = dir_out / "seqids.csv"
    csv_lookup_table = Path(csv_lookup_table)
    fastas_out = {segment: dir_out/f"{segment}.fasta" for segment in SEGMENTS}

    attrs_list = []
    for fasta_in in fastas_in:
        for record in SeqIO.parse(fasta_in, "fasta"):
            attrs = parse_vdj_id(record.id)
            attrs["path"] = fasta_in
            attrs["seq"] = str(record.seq)
            attrs_list.append(attrs)

    # modify IDs where needed, using a lookup table from segments to FASTA
    # paths.  segments with multiple paths will get a suffix appended for each
    # ID, if needed, to differentiate.
    path_suffix_lut = _make_seqid_suffix_lut(attrs_list)
    for attrs in attrs_list:
        seqid_here = attrs["allele"]
        suffix = path_suffix_lut.get(attrs["path"])
        if suffix:
            seqid_here += "_" + suffix
        attrs["seqid_here"] = seqid_here

    # write output: both the FASTA trio and the lookup table as CSV
    if not dry_run:
        csv_lookup_table.parent.mkdir(parents=True, exist_ok=True)
        dir_out.mkdir(parents=True, exist_ok=True)
        with open(csv_lookup_table, "wt") as csv_out:
            writer = csv.DictWriter(
                csv_out,
                fieldnames=["sequence_id", "sequence_id_original", "input_path"],
                lineterminator="\n")
            writer.writeheader()
            for segment, attrs_group in group(attrs_list).items():
                attrs_group = sorted(
                    attrs_group,
                    key=lambda r: (r["path"], r["locus"], r["segment"], r["family"], r["seqid"]))
                with open(fastas_out[segment], "wt") as f_out:
                    for attrs in attrs_group:
                        SeqIO.write(
                            SeqRecord(Seq(attrs["seq"]), id=attrs["seqid_here"], description=""),
                            f_out, "fasta-2line")
                        writer.writerow({
                            "sequence_id": attrs["seqid_here"],
                            "sequence_id_original": attrs["seqid"],
                            "input_path": attrs["path"]})
    return attrs_list

def _make_seqid_suffix_lut(attrs_list):
    # mapping from path -> suffix to append to sequence IDs from that path
    segment_path_lut = {}
    for attrs in attrs_list:
        if attrs["segment"] not in segment_path_lut:
            segment_path_lut[attrs["segment"]] = set()
        segment_path_lut[attrs["segment"]].add(attrs["path"])
    path_suffix_lut = {}
    for paths in segment_path_lut.values():
        if len(paths) == 1:
            # we only need to worry about any of this if there's more than one
            # path
            continue
        internal = set()
        external = set()
        for path in paths:
            if Path(path).is_relative_to(util.DATA/"germ"):
                internal.add(path)
            else:
                external.add(path)
        for path in paths:
            if path in external and len(external) == 1:
                # special case: with just one external ref we don't need a suffix
                pass
            else:
                if path in internal:
                    fasta_rel = Path(path).relative_to(util.DATA/"germ")
                    parents = [parent.name for parent in fasta_rel.parents if parent.name]
                    species = parents[-1]
                    ref = parents[-2]
                    suffix = f"{species}/{ref}"
                else:
                    suffix = hashlib.sha224(path.encode("ascii")).hexdigest()[:6]
                path_suffix_lut[path] = suffix
    return path_suffix_lut

def parse_vdj_id(seqid):
    """Parse a "standard" V/D/J segment gene/allele ID into a dictionary."""
    seqid = str(seqid)
    # IGHV1-NL_1*01_S2052
    match = re.search(
        r"(IG[HKL][VDJ][0-9])(-?[A-Za-z0-9_]*-?[A-Za-z0-9]*)\*?([0-9]*_?S?[0-9]*)", seqid)
    if not match:
        raise ValueError("Seq ID not recognized: %s" % seqid)
    prefix = seqid[0:match.start()]
    suffix = seqid[match.end():]
    attrs = {
        "seqid": seqid,
        "prefix": prefix,
        "suffix": suffix,
        "allele": match.group(0) if match.group(3) else "",
        "gene": match.group(1) + match.group(2) if match.group(2) else "",
        "family": match.group(1),
        "segment": match.group(1)[:4],
        "locus": match.group(1)[:3]}
    return attrs

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
        segment = re.sub("^IG[HKL]", "", entry["segment"])
        groups[key][segment].append(entry)
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
