"""
Find closest-matching VDJ sequences from one or more references.

This uses IgBLAST to assign germline genes, but with a separate query for each
reference specified as a separate database.
"""

import logging
import subprocess
from csv import DictReader, DictWriter
from io import StringIO
from pathlib import Path
from . import util
from . import igblast
from . import show
from . import vdj

LOGGER = logging.getLogger(__name__)

def vdj_match(ref_paths, query, output=None, showtxt=None, species=None, dry_run=False, threads=1):
    LOGGER.info("given ref path(s): %s", ref_paths)
    LOGGER.info("given query path: %s", query)
    LOGGER.info("given output: %s", output)
    LOGGER.info("given showtxt: %s", showtxt)
    LOGGER.info("given species: %s", species)
    # if not specified, show text when not saving output
    if showtxt is None:
        showtxt = not output
        LOGGER.info("detected showtxt: %s", showtxt)
    attrs_list = vdj.parse_vdj_paths(ref_paths)

    species_det = {attrs.get("species") for attrs in attrs_list}
    species_det = {s for s in species_det if s}
    organism = igblast.detect_organism(species_det, species)

    vdj_files_grouped = vdj.group(
        attrs_list,
        lambda x: f"{x['species']}/{x['ref']}" if x["type"] == "internal" else x["input"])
    for key, trio in vdj_files_grouped.items():
        LOGGER.info("detected V FASTA from %s: %d", key, len(trio["V"]))
        LOGGER.info("detected D FASTA from %s: %d", key, len(trio["D"]))
        LOGGER.info("detected J FASTA from %s: %d", key, len(trio["J"]))
        for segment, attrs_group in trio.items():
            if not attrs_group:
                LOGGER.critical("No FASTA for %s from %s", segment, key)
                raise util.IgSeqError("Missing VDJ input for database")

    if not dry_run:
        results = []
        for key, trio in vdj_files_grouped.items():
            paths = [attrs["path"] for attrs in trio["V"] + trio["D"] + trio["J"]]
            proc = igblast.setup_db_dir_and_igblast(
                paths, organism, query, threads=threads, extra_args=["-outfmt", "19"],
                stdout=subprocess.PIPE, text=True, check=True)
            reader = DictReader(StringIO(proc.stdout), delimiter="\t")
            for row in reader:
                for segment in ["v", "d", "j"]:
                    try:
                        start = int(row[f"{segment}_sequence_start"])
                        stop = int(row[f"{segment}_sequence_end"])
                        length = stop - start + 1
                    except ValueError:
                        length = ""
                    results.append({
                        "query": row["sequence_id"],
                        "reference": key,
                        "segment": segment.upper(),
                        "call": row[f"{segment}_call"],
                        "length": length,
                        "identity": row[f"{segment}_identity"]})
        results = sorted(results, key=lambda r: (r["query"], r["reference"]))
        if showtxt:
            show.show_grid(results)
        if output:
            output = Path(output)
            output.parent.mkdir(parents=True, exist_ok=True)
            with open(output, "wt") as f_out:
                writer = DictWriter(f_out, fieldnames=results[0].keys(), lineterminator="\n")
                writer.writeheader()
                writer.writerows(results)
