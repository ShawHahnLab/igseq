"""
Extract subsequences based on AIRR region via IgBLAST.

Input and output can be in any file format supported by the convert command and
can be given as "-" for standard input/output.

    igseq cut -S rhesus -R cdr3 seqs.fa cdr3.fa
    igseq cut -S rhesus -R fwr1-fwr3 V.fa cysTruncated.fa
"""

import re
import logging
from csv import DictReader
from . import util
from . import igblast
from . import vdj
from .record import RecordWriter

LOGGER = logging.getLogger(__name__)

# These things have _start and _end columns defined that specify start and end
# positions relative to a query sequence
# https://docs.airr-community.org/en/stable/datarep/rearrangements.html
AIRR_REGIONS = [
    "v_sequence",
    "d_sequence",
    "j_sequence",
    "c_sequence",
    "cdr1",
    "cdr2",
    "cdr3",
    "fwr1",
    "fwr2",
    "fwr3",
    "fwr4",
    ]

def cut(
        ref_paths, query, output, regions,
        species=None, fmt_in=None, fmt_out=None, colmap=None, dry_run=False, threads=1):
    """Extract subsequences based on AIRR region via IgBLAST"""
    LOGGER.info("given ref path(s): %s", ref_paths)
    LOGGER.info("given query path: %s", query)
    LOGGER.info("given output: %s", output)
    LOGGER.info("given regions: %s", regions)
    LOGGER.info("given species: %s", species)
    LOGGER.info("given input format: %s", fmt_in)
    LOGGER.info("given colmap: %s", colmap)
    LOGGER.info("given threads: %s", threads)
    if species and not ref_paths:
        # If only species is given, default to using all available reference
        # sets for that species
        ref_paths = [igblast.fuzzy_species_match(species)]
        LOGGER.info("inferred ref path: %s", ref_paths[0])
    attrs_list = vdj.parse_vdj_paths(ref_paths)
    species_det = {attrs.get("species") for attrs in attrs_list}
    species_det = {s for s in species_det if s}
    organism = igblast.detect_organism(species_det, species)
    attrs_list_grouped = vdj.group(attrs_list)
    for key, attrs_group in attrs_list_grouped.items():
        LOGGER.info("detected %s references: %d", key, len(attrs_group))
        if len(attrs_group) == 0:
            raise util.IgSeqError(f"No references for segment {key}")

    parts = parse_regions(regions)
    LOGGER.info("parsed regions: {parts}")

    if not dry_run:
        with igblast.setup_db_dir(
            [attrs["path"] for attrs in attrs_list]) as (db_dir, _):
            with igblast.run_igblast(
                    db_dir, organism, query, threads,
                    fmt_in, colmap, extra_args=["-outfmt", "19"]) as proc:
                reader = DictReader(proc.stdout, delimiter="\t")
                with RecordWriter(output, fmt_out, colmap, dry_run=dry_run) as writer:
                    for row in reader:
                        # TODO what should happen if a region isn't present at
                        # all?  Fall back on nearest position that does exist?
                        # TODO should this also handle the non-positioned based
                        # entries like np1 and junction?  (But then, what to do
                        # with something like "np1-J"?)
                        pos_start = int(row[parts[0] + "_start"])
                        pos_end = int(row[parts[1] + "_end"])
                        seq = row["sequence"][pos_start-1:pos_end]
                        record = {"sequence_id": row["sequence_id"], "sequence": seq}
                        writer.write(record)

def parse_regions(regions):
    parts = [p.lower() for p in re.split(r"[-:]", regions, 1)]
    # support implicit ..._sequence
    for idx, part in enumerate(parts):
        if f"{part}_sequence" in AIRR_REGIONS:
            parts[idx] = f"{part}_sequence"
    if len(parts) == 1:
        parts *= 2
    for part in parts:
        if part not in AIRR_REGIONS:
            raise util.IgSeqError(f"region \"{part}\" not recognized")
    if AIRR_REGIONS.index(parts[0]) > AIRR_REGIONS.index(parts[1]):
        LOGGER.warning("regions out of order; will reverse automatically: %s", parts)
        parts.reverse()
    return parts
