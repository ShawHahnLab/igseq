"""
Find closest-matching VDJ sequences from one or more references.

This uses IgBLAST to assign germline genes, but with a separate query for each
database.
"""

import logging
import subprocess
from csv import DictReader, DictWriter
from io import StringIO
from tempfile import TemporaryDirectory
from pathlib import Path
from . import util
from . import vdj_gather
from . import igblast
from . import show

LOGGER = logging.getLogger(__name__)

def vdj_match(db_paths, query,  output=None, showtxt=None, species=None, dry_run=False, threads=1):
    LOGGER.info("given DB path(s): %s", db_paths)
    LOGGER.info("given query path: %s", query)
    LOGGER.info("given output: %s", output)
    LOGGER.info("given showtxt: %s", showtxt)
    LOGGER.info("given species: %s", species)
    dbs = parse_vdj_paths(db_paths)
    LOGGER.info("detected DBs: %d", len(dbs))
    # if not specified, show text when not saving output
    if showtxt is None:
        showtxt = not output
        LOGGER.info("detected showtxt: %s", showtxt)
    species_det = set()
    for attrs in dbs:
        LOGGER.info("detected db type: %s", attrs["type"])
        LOGGER.info("detected db path: %s", attrs["path"])
        if attrs["type"] == "internal":
            spec = attrs["path"].parent.name
            LOGGER.info("detected db species: %s", spec)
            species_det.add(spec)
    if not species and not species_det:
        raise util.IgSeqError(
            "species not detected from input.  specify a species manually.")
    if not species and len(species_det) > 1:
        raise util.IgSeqError(
            "multiple species detected from input.  specify a species manually.")
    if not species:
        species = species_det.pop()
        LOGGER.info("detected species: %s", species)
    try:
        species_igblast = igblast.SPECIESMAP[species]
    except KeyError as err:
        keys = str(igblast.SPECIESMAP.keys())
        raise util.IgSeqError("species not recognized.  should be one of: %s" % keys) from err
    LOGGER.info("detected IgBLAST organism: %s", species_igblast)
    if not dry_run:
        results = []
        for attrs in dbs:
            with TemporaryDirectory() as tmp:
                vdj_gather._vdj_gather(attrs["path"], tmp)
                igblast._run_makeblastdb(f"{tmp}/V.fasta")
                igblast._run_makeblastdb(f"{tmp}/D.fasta")
                igblast._run_makeblastdb(f"{tmp}/J.fasta")
                args = [
                    "-germline_db_V", f"{tmp}/V",
                    "-germline_db_D", f"{tmp}/D",
                    "-germline_db_J", f"{tmp}/J",
                    "-query", query,
                    "-auxiliary_data", f"optional_file/{species_igblast}_gl.aux",
                    "-organism", species_igblast,
                    "-ig_seqtype", "Ig",
                    "-outfmt", 19,
                    "-num_threads", threads]
                proc = igblast._run_igblastn(args, stdout=subprocess.PIPE, text=True)
                reader = DictReader(StringIO(proc.stdout), delimiter="\t")
                for row in reader:
                    if attrs["type"] == "internal":
                        name = str(attrs["path"].relative_to(util.DATA / "germ"))
                    else:
                        name = str(attrs["path"])
                    for segment in ["v", "d", "j"]:
                        results.append({
                            "query": row["sequence_id"],
                            "database": name,
                            "segment": segment.upper(),
                            "call": row[f"{segment}_call"],
                            "identity": row[f"{segment}_identity"]})
        results = sorted(results, key=lambda r: (r["query"], r["database"]))
        if showtxt:
            show.show_grid(results)
        if output:
            output = Path(output)
            output.parent.mkdir(parents=True, exist_ok=True)
            with open(output, "wt") as f_out:
                writer = DictWriter(f_out, fieldnames=results[0].keys(), lineterminator="\n")
                writer.writerows(results)

def parse_vdj_paths(db_paths):
    # TODO unify with vdj_gather's version
    parsed = []
    for entry in db_paths:
        internal_matches = get_internal_vdj(entry)
        path = Path(entry)
        # first priority: actual directory name
        if path.is_dir():
            parsed.append({
                "type": "dir",
                "path": path})
        # second priority: internal reference
        elif internal_matches:
            for ref in internal_matches:
                parsed.append({
                    "type": "internal",
                    "path": ref})
        else:
            raise util.IgSeqError(
                "db path isn't a directory or an internal germline reference: %s" % entry)
    return parsed

def get_internal_vdj(name):
    candidates = []
    germ = util.DATA / "germ"
    for path in germ.glob("*"):
        if path.is_dir():
            for subpath in path.glob("*"):
                if subpath.is_dir():
                    candidates.append(subpath)
    output = []
    for candidate in candidates:
        if name in str(candidate.relative_to(germ)):
            output.append(candidate)
    return output
