"""
Gather VDJ sequences into one directory.

This is handy for starting off an IgDiscover run from one of the built-in
references, or manually prepping an IgBLAST database.
"""

import logging
from pathlib import Path
from Bio import SeqIO
from . import util

LOGGER = logging.getLogger(__name__)

def vdj_gather(db_paths, dir_path_out, dry_run=False):
    """Gather V/D/J FASTA into one directory as V.fasta, D.fasta, J.fasta."""
    # front-end for _vdj_gather with support for internal paths
    LOGGER.info("input db path(s): %s", db_paths)
    LOGGER.info("input dir_path_out: %s", dir_path_out)
    paths = parse_vdj_paths(db_paths)
    LOGGER.info("detected db type: %s", paths["type"])
    LOGGER.info("detected db details: %s", {k: v for k, v in paths.items() if k != "type"})
    if paths["type"] == "filetrio":
        raise NotImplementedError
    if paths["type"] == "dir":
        if not dry_run:
            _vdj_gather(paths["path"], dir_path_out)
    elif paths["type"] == "internal":
        if not dry_run:
            _vdj_gather(paths["path"], dir_path_out)

def _vdj_gather(dir_path_in, dir_path_out):
    # take VDJ found and aggregate by segment in one place.
    dir_path_in = Path(dir_path_in)
    dir_path_out  = Path(dir_path_out)
    # gather up any recognized file paths
    parts = {"V": [], "D": [], "J": []}
    for path in dir_path_in.glob("**/*"):
        if path.suffix.lower() in [".fasta", ".fa", ".fna"]:
            segment = parse_segment(path.name)
            if segment:
                parts[segment].append(path)
            else:
                LOGGER.warning("Found FASTA but doesn't look like V/D/J: %s", path)
    if not all(parts.values()):
        LOGGER.error(
            "Missing input for: %s",
            [key for key, val in parts.items() if not val])
    dir_path_out.mkdir(parents=True, exist_ok=True)
    # combine by segment
    for segment, paths in parts.items():
        with open(dir_path_out/f"{segment}.fasta", "wt") as f_out:
            for path in paths:
                for record in SeqIO.parse(path, "fasta"):
                    SeqIO.write(record, f_out, "fasta-2line")

def parse_segment(txt):
    # e.g. "IGHV.fasta" -> "V"
    segments = ["V", "D", "J"]
    match = [letter for letter in segments if letter in txt.upper()]
    if len(match) != 1:
        return None
    return match[0]

def parse_vdj_paths(db_paths):
    if len(db_paths) == 3 and all([Path(path).is_file() for path in db_paths]):
        return {
            "type": "filetrio",
            "V": Path(db_paths[0]),
            "D": Path(db_paths[1]),
            "J": Path(db_paths[2])}
    if len(db_paths) == 1 and Path(db_paths[0]).is_dir():
        return {
            "type": "dir",
            "path": Path(db_paths[0])}
    if len(db_paths) == 1:
        for filepath in util.FILES:
            for parent in filepath.parents:
                if str(parent).endswith(db_paths[0]):
                    return {
                        "type": "internal",
                        "path": parent}
    raise ValueError("Input not recognized: %s" % str(db_paths))
