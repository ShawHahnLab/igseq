"""
Find closest-matching VDJ sequences from one or more references.

This uses IgBLAST to assign germline genes, but with a separate query for each
reference specified as a separate database.
"""

import logging
import subprocess
from csv import DictReader, DictWriter
from io import StringIO
from tempfile import TemporaryDirectory
from pathlib import Path
from . import util
from . import igblast
from . import show
from . import vdj

LOGGER = logging.getLogger(__name__)

def vdj_match(ref_paths, query,  output=None, showtxt=None, species=None, dry_run=False, threads=1):
    LOGGER.info("given ref path(s): %s", ref_paths)
    LOGGER.info("given query path: %s", query)
    LOGGER.info("given output: %s", output)
    LOGGER.info("given showtxt: %s", showtxt)
    LOGGER.info("given species: %s", species)
    # if not specified, show text when not saving output
    if showtxt is None:
        showtxt = not output
        LOGGER.info("detected showtxt: %s", showtxt)

    vdj_files = vdj.parse_vdj_paths(ref_paths)
    vdj_files_grouped = _group_by_input_by_segment(vdj_files)
    for key, trio in vdj_files_grouped.items():
        LOGGER.info("detected V FASTA from %s: %d", key, len(trio["V"]))
        LOGGER.info("detected D FASTA from %s: %d", key, len(trio["D"]))
        LOGGER.info("detected J FASTA from %s: %d", key, len(trio["J"]))
        for segment, attrs in trio.items():
            if not attrs:
                LOGGER.critical("No FASTA for %s from %s", segment, key)
                raise util.IgSeqError("Missing VDJ input for database")

    species_det = set()
    for attrs in vdj_files:
        LOGGER.info("detected ref path: %s", attrs["path"])
        LOGGER.info("detected ref type: %s", attrs["type"])
        if attrs["type"] == "internal":
            LOGGER.info("detected db species: %s", attrs["species"])
            species_det.add(attrs["species"])
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
        for key, trio in vdj_files_grouped.items():
            with TemporaryDirectory() as tmp:
                aux = f"{tmp}/{species_igblast}_gl.aux"
                for segment, attrs_list in trio.items():
                    fasta = f"{tmp}/{segment}.fasta"
                    vdj.vdj_gather(attrs_list, fasta)
                    if segment == "J":
                        igblast.make_aux_file(fasta, aux)
                igblast.makeblastdbs(tmp)
                args = [
                    "-germline_db_V", f"{tmp}/V",
                    "-germline_db_D", f"{tmp}/D",
                    "-germline_db_J", f"{tmp}/J",
                    "-query", query,
                    "-auxiliary_data", aux,
                    "-organism", species_igblast,
                    "-ig_seqtype", "Ig",
                    "-outfmt", 19,
                    "-num_threads", threads]
                proc = igblast._run_igblastn(args, stdout=subprocess.PIPE, text=True)
                reader = DictReader(StringIO(proc.stdout), delimiter="\t")
                for row in reader:
                    for segment in ["v", "d", "j"]:
                        start = int(row[f"{segment}_sequence_start"])
                        stop = int(row[f"{segment}_sequence_end"])
                        results.append({
                            "query": row["sequence_id"],
                            "reference": key,
                            "segment": segment.upper(),
                            "call": row[f"{segment}_call"],
                            "length": stop - start + 1,
                            "identity": row[f"{segment}_identity"]})
        results = sorted(results, key=lambda r: (r["query"], r["reference"]))
        if showtxt:
            show.show_grid(results)
        if output:
            output = Path(output)
            output.parent.mkdir(parents=True, exist_ok=True)
            with open(output, "wt") as f_out:
                writer = DictWriter(f_out, fieldnames=results[0].keys(), lineterminator="\n")
                writer.writerows(results)

def _group_by_input_by_segment(vdj_path_attrs):
    # for prepping multiple separate databases (rather than one big one
    # combining by segment)
    refs = {}
    for attrs in vdj_path_attrs:
        # is this enough?
        key = attrs["input"]
        if key not in refs:
            refs[key] = {s: [] for s in vdj.SEGMENTS}
        refs[key][attrs["segment"]].append(attrs)
    return refs
