"""
Run IgBLAST.
"""

import logging
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
from . import vdj_gather

LOGGER = logging.getLogger(__name__)

IGBLASTN = "igblastn"
MAKEBLASTDB = "makeblastdb"

# from our generic names to the IgBLAST names
SPECIESMAP = {
    "rhesus": "rhesus_monkey",
    "human": "human"}

def igblast(db_paths, query_path, species=None, dry_run=False, threads=1):
    LOGGER.info("input db path(s): %s", db_paths)
    LOGGER.info("input query path: %s", query_path)
    LOGGER.info("input species: %s", species)
    LOGGER.info("input threads: %s", threads)
    # if db_paths is a single dir, find all fastas
    # if three files, use in order
    # if species given, use that, if not, try to parse from db_paths
    paths = vdj_gather.parse_vdj_paths(db_paths)
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
                vdj_gather._vdj_gather(paths["path"], tmp)
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
                vdj_gather._vdj_gather(paths["path"], tmp)
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

def _run_igblastn(args, **kwargs):
    args = [IGBLASTN] + [str(arg) for arg in args]
    LOGGER.info("igblastn command: %s", args)
    return subprocess.run(args, check=True, **kwargs)
