"""
Run IgBLAST.
"""

import logging
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
from Bio import SeqIO
from . import util

LOGGER = logging.getLogger(__name__)

IGBLASTN = "igblastn"
MAKEBLASTDB = "makeblastdb"

# from our generic names to the IgBLAST names
SPECIESMAP = {
    "rhesus": "rhesus_monkey",
    "human": "human"}

def gather_germline(dir_path_in, dir_path_out):
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
            "Missing germline input for: %s",
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

def igblast(db_paths, query_path, species=None, dry_run=False, threads=1):
    LOGGER.info("input db path(s): %s", db_paths)
    LOGGER.info("input query path: %s", query_path)
    LOGGER.info("input species: %s", species)
    LOGGER.info("input threads: %s", threads)
    # if db_paths is a single dir, find all fastas
    # if three files, use in order
    # if species given, use that, if not, try to parse from db_paths
    paths = _parse_db_paths(db_paths)
    LOGGER.info("detected db type: %s", paths["type"])
    LOGGER.info("detected db details: %s", {k: v for k, v in paths.items() if k != "type"})

    try:
        if not dry_run:
            tmpdb = TemporaryDirectory()
            tmp = Path(tmpdb.name)

        if paths["type"] == "filetrio":
            for segment in ["V", "D", "J"]:
                if not dry_run:
                    with open(tmp/f"{segment}.fasta", "wt") as f_out, open(paths[segment]) as f_in:
                        f_out.write(f_in.read())
        elif paths["type"] == "dir":
            if not dry_run:
                gather_germline(paths["path"], tmp)
        elif paths["type"] == "internal":
            if not species:
                for key in SPECIESMAP:
                    if key in str(paths["path"]):
                        species = key
                        LOGGER.info("detected species from db: %s", species)
                        break
            if not species:
                LOGGER.error("Could not determine species from database input %s", paths["path"])
            if not dry_run:
                gather_germline(paths["path"], tmp)
        if not species:
            LOGGER.error("Species is required")
            raise ValueError
        # TODO allow synonyms
        if species not in SPECIESMAP:
            LOGGER.error("Species not recognized: %s", species)
            raise ValueError
        species_igblast = SPECIESMAP[species]
        LOGGER.info("species name for igblastn: %s", species_igblast)
        if not dry_run:
            _run_makeblastdb(f"{tmp}/V.fasta")
            _run_makeblastdb(f"{tmp}/D.fasta")
            _run_makeblastdb(f"{tmp}/J.fasta")
            args = [
                "-germline_db_V", f"{tmp}/V",
                "-germline_db_D", f"{tmp}/D",
                "-germline_db_J", f"{tmp}/J",
                "-query", query_path,
                "-auxiliary_data", f"optional_file/{species_igblast}_gl.aux",
                "-organism", species_igblast,
                "-ig_seqtype", "Ig",
                "-num_threads", threads]
            _run_igblastn(args)
    finally:
        if not dry_run:
            tmpdb.cleanup()

def _run_makeblastdb(path_fasta):
    path_fasta = Path(path_fasta)
    args = [
        MAKEBLASTDB,
        "-dbtype", "nucl",
        "-parse_seqids",
        "-in", path_fasta,
        "-out", path_fasta.parent/path_fasta.stem]
    args = [str(arg) for arg in args]
    subprocess.run(args, check=True, stdout=subprocess.DEVNULL)

def _run_igblastn(args):
    args = [IGBLASTN] + [str(arg) for arg in args]
    LOGGER.info("igblastn command: %s", args)
    subprocess.run(args, check=True)

def _parse_db_paths(db_paths):
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
    raise ValueError("IgBLAST database input not recogized: %s" % str(db_paths))
