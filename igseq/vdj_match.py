"""
Find closest-matching VDJ sequences from one or more references.

This uses IgBLAST to assign germline genes, but with a separate query for each
reference specified as a separate database.  The query can be in any file
format supported by the convert command and can be given as "-" for standard
input.
"""

import logging
from csv import DictReader, DictWriter
from pathlib import Path
from . import util
from . import igblast
from . import show
from . import vdj

LOGGER = logging.getLogger(__name__)

def vdj_match(ref_paths, query, output=None, showtxt=None, species=None, fmt_in=None, colmap=None, dry_run=False, threads=1):
    LOGGER.info("given ref path(s): %s", ref_paths)
    LOGGER.info("given query path: %s", query)
    LOGGER.info("given output: %s", output)
    LOGGER.info("given showtxt: %s", showtxt)
    LOGGER.info("given species: %s", species)
    LOGGER.info("given input format: %s", fmt_in)
    LOGGER.info("given colmap: %s", colmap)
    LOGGER.info("given threads: %s", threads)
    # if not specified, show text when not saving output
    if showtxt is None:
        showtxt = not output
        LOGGER.info("detected showtxt: %s", showtxt)
    attrs_list = vdj.parse_vdj_paths(ref_paths)

    species_det = {attrs.get("species") for attrs in attrs_list}
    species_det = {s for s in species_det if s}
    organism = igblast.detect_organism(species_det, species)

    # Group references into sets with V/D/J.  If a segment is missing (note
    # that D is always needed by igblast even for light chain) we'll remove the
    # whole reference with a warning *if* it was found internally from a
    # partial text match, but will fail outright if it was selected explicitly
    # or supplied from external files.
    def groupkey(row):
        return(f"{row['species']}/{row['ref']}" if row["type"] == "internal" else row["input"])
    vdj_files_grouped = vdj.group(attrs_list, groupkey)
    # This is sloppy because it'll just take whichever is the last matching row
    # (with its segment-specific stuff included) but i'll do for now
    orig_attrs = {groupkey(row): row for row in attrs_list}
    skips = set()
    for key, trio in vdj_files_grouped.items():
        LOGGER.info("detected V FASTA from %s: %d", key, len(trio["V"]))
        LOGGER.info("detected D FASTA from %s: %d", key, len(trio["D"]))
        LOGGER.info("detected J FASTA from %s: %d", key, len(trio["J"]))
        for segment, attrs_group in trio.items():
            if key not in skips and not attrs_group:
                # Only case where we'll just warn and continue is: it's an
                # internal ref, and it was found implicitly (not by name).
                if orig_attrs[key]["type"] == "internal" and orig_attrs[key]["input"] != key:
                    LOGGER.warning("No FASTA for %s from %s; skipping", segment, key)
                    skips.add(key)
                else:
                    LOGGER.critical("No FASTA for %s from %s", segment, key)
                    raise util.IgSeqError("Missing VDJ input for database")
    for key in skips:
        del vdj_files_grouped[key]

    if not dry_run:
        results = []
        for key, trio in vdj_files_grouped.items():
            paths = [attrs["path"] for attrs in trio["V"] + trio["D"] + trio["J"]]
            with igblast.setup_db_dir(paths) as (db_dir, _):
                with igblast.run_igblast(db_dir, organism, query, threads,
                        fmt_in, colmap, extra_args=["-outfmt", "19"]) as proc:
                    reader = DictReader(proc.stdout, delimiter="\t")
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
