"""
Create identity/divergence and coverage percentage tables like SONAR does.

Output is mainly a pair of wide-format CSV with each row being one query and
comparisons to germline V and mAbs on columns (see SONAR's
2.1-calculate_id-div.py).  Rather than multiple sequence alignments via MUSCLE
or CLUSTAL as SONAR uses, the identity and coverage values here are calculated
with pairwise alignments, and germline divergence is determined with IgBLAST
and a full set of gene segment references (rather than trimmed V).
"""

import re
import logging
from pathlib import Path
from csv import DictReader
from tempfile import NamedTemporaryFile
from .record import RecordReader, RecordWriter
from .identity import align, calc_identity, calc_coverage
from .igblast import detect_organism, setup_db_dir, run_igblast
from .vdj import parse_vdj_paths
from .parallel import parallelize

LOGGER = logging.getLogger(__name__)

def id_div(path_in, path_abs, path_id_div, path_cov, path_airr_update=None, paths_ref=None, species=None, dry_run=False, threads=1):
    """Create SONAR-style ID/DIV CSV file."""
    LOGGER.info("given query path: %s", path_in)
    LOGGER.info("given antibodies path: %s", path_abs)
    LOGGER.info("given ID/DIV path: %s", path_id_div)
    LOGGER.info("given coverage path: %s", path_cov)
    LOGGER.info("given ref path(s): %s", paths_ref)
    LOGGER.info("given AIRR for update: %s", path_airr_update)
    LOGGER.info("given species: %s", species)
    LOGGER.info("given threads: %s", threads)

    if not dry_run:
        with RecordReader(path_abs) as reader:
            abs_records = sorted(list(reader), key=lambda row: row["sequence_id"])
        fields_abs = [ab["sequence_id"] for ab in abs_records]
        fields_id_div = ["sequence_id", "v_gene", "germ_div"] + fields_abs
        fields_cov = ["sequence_id", "germ_cov"] + fields_abs
        with RecordWriter(path_id_div, colmap={f: f for f in fields_id_div}) as writer_id_div, \
            RecordWriter(path_cov, colmap={f: f for f in fields_cov}) as writer_cov:
            # If the input itself looks like AIRR, skip IgBLAST (assuming the TSV is from IgBLAST)
            suffixes = Path(path_in.lower()).suffixes
            if (len(suffixes) > 1 and suffixes[-2:] == [".tsv", ".gz"]) or \
                (len(suffixes) > 0 and suffixes[-1] == ".tsv"):
                LOGGER.info("Reading germline info directly for input AIRR")
                with RecordReader(path_in) as reader:
                    v_identity_map = _process_rows(reader, writer_id_div, writer_cov, abs_records, threads)
            else:
                attrs_list = parse_vdj_paths(paths_ref)
                for attrs in attrs_list:
                    LOGGER.info("detected ref path: %s", attrs["path"])
                    LOGGER.info("detected ref type: %s", attrs["type"])
                species_det = {attrs.get("species") for attrs in attrs_list}
                species_det = {s for s in species_det if s}
                organism = detect_organism(species_det, species)
                vdj_ref_paths = [attrs["path"] for attrs in attrs_list]
                LOGGER.info("Running IgBLAST")
                with setup_db_dir(vdj_ref_paths) as (db_dir, _):
                    with run_igblast(
                        db_dir, organism, path_in, threads, extra_args=["-outfmt", "19"]) as proc:
                        v_identity_map = _process_rows(
                            DictReader(proc.stdout, delimiter="\t"),
                            writer_id_div, writer_cov, abs_records, threads)
                        LOGGER.debug("process rows done")
        # If we're updating a previous AIRR table, add a column in a temporary
        # copy, and then move it into place.
        if path_airr_update:
            LOGGER.info("Updating SONAR AIRR Table with v_identity")
            _update_sonar_airr(path_airr_update, v_identity_map)

def _update_sonar_airr(path_airr_update, v_identity_map):
    with NamedTemporaryFile("wt") as airr_tmp:
        with RecordReader(path_airr_update) as reader, \
            RecordWriter(airr_tmp.name, fmt="tsv") as writer:
            for row in reader:
                # The AIRR rearrangements definition states that the
                # _identity values should be fractional, not
                # percentages, but IgBLAST (as of 1.17.1 at least) gives
                # percentages.  SONAR follows the spec and uses fractional
                # values.
                # We'll try to handle either here and follow the spec.
                # This assumes we won't encounter a value below 1% but
                # in that case we've got bigger problems anyway.
                v_ident = v_identity_map.get(row["sequence_id"])
                if v_ident is not None:
                    v_ident = float(v_ident)
                    if v_ident > 1:
                        v_ident = v_ident / 100
                    row["v_identity"] = f"{v_ident:0.3}"
                writer.write(row)
        # Path().rename() doesn't work cross-filesystem so some
        # equivalent of cp+rm is needed instead.
        with open(airr_tmp.name) as f_in, open(path_airr_update, "wt") as f_out:
            for line in f_in:
                f_out.write(line)

def _worker(rows, abs_records):
    id_div_rows = []
    cov_rows = []
    v_identity_map = {}
    LOGGER.debug("worker func started")
    for row in rows:
        alns = {ab["sequence_id"]: align(row["sequence"], ab["sequence"]) for ab in abs_records}
        v_identity_map[row["sequence_id"]] = row["v_identity"]
        id_div_rows.append(_prep_id_div_row(row, alns))
        cov_rows.append(_prep_cov_row(row, alns))
    return id_div_rows, cov_rows, v_identity_map

def _receiver(results, id_div_w, cov_w):
    for id_div_row in results[0]:
        id_div_w.write(id_div_row)
    for cov_row in results[0]:
        cov_w.write(cov_row)
    #v_identity_map.update(results[2])

def _process_rows(airr_iter, id_div_w, cov_w, abs_records, threads):
    v_identity_map = {}
    parallelize(airr_iter, threads, _worker, _receiver, 1000, w_args=[abs_records], r_args=[id_div_w, cov_w])
    return v_identity_map

def _prep_id_div_row(airr_row, alns):
    score_div = 100 - float(airr_row["v_identity"])
    row_out = {
        "sequence_id": airr_row["sequence_id"],
        "v_gene": re.sub(r"\*.*$", "", airr_row["v_call"]),
        "germ_div": f"{score_div:.1f}"}
    for ab_id, aln in alns.items():
        score = 100 * calc_identity(aln)
        row_out[ab_id] = f"{score:.1f}"
    return row_out

def _prep_cov_row(airr_row, alns):
    len_seq = int(airr_row["v_sequence_end"]) - int(airr_row["v_sequence_start"])
    len_germ = int(airr_row["v_germline_end"]) - int(airr_row["v_germline_start"])
    score_cov = 100 * len_seq / len_germ
    row_out = {
        "sequence_id": airr_row["sequence_id"],
        "germ_cov": f"{score_cov:.1f}"}
    for ab_id, aln in alns.items():
        score = 100 * calc_coverage(aln)
        row_out[ab_id] = f"{score:.1f}"
    return row_out
