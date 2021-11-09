"""
Make a summary table of antibody attributes, via IgBLAST.
"""

import logging
import subprocess
from csv import DictReader, DictWriter
from io import StringIO
from pathlib import Path
from Bio import SeqIO
from . import util
from . import igblast
from . import show
from . import vdj

LOGGER = logging.getLogger(__name__)

def summarize(ref_paths, query, output=None, showtxt=None, species=None, dry_run=False, threads=1):
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
    species_igblast = igblast.detect_species(attrs_list, species)
    attrs_list_grouped = vdj.group(attrs_list)
    for key, attrs_group in attrs_list_grouped.items():
        LOGGER.info("detected %s references: %d", key, len(attrs_group))
        if len(attrs_group) == 0:
            raise util.IgSeqError("No references for segment %s" % key)

    if not dry_run:
        results = []
        proc = igblast.setup_db_and_igblast(
            attrs_list_grouped, species_igblast, query, threads=threads,
            extra_args=["-outfmt", "19"],
            stdout=subprocess.PIPE, text=True)
        reader = DictReader(StringIO(proc.stdout), delimiter="\t")
        lengthmap = {}
        for attrs_group in attrs_list_grouped.values():
            for attrs in attrs_group:
                with open(attrs["path"]) as f_in:
                    for record in SeqIO.parse(f_in, "fasta"):
                        lengthmap[record.id] = len(record)
        for row in reader:
            cdrlens = []
            for cdr in ["cdr1_aa", "cdr2_aa", "cdr3_aa"]:
                cdrlens.append(str(len(row[cdr])) if row[cdr] else "?")
            cdrlens = ".".join(cdrlens)
            idents_pct = []
            idents_nt = []
            for segment in ["v", "d", "j"]:
                try:
                    start = int(row[f"{segment}_sequence_start"])
                    stop = int(row[f"{segment}_sequence_end"])
                    length1 = stop - start + 1
                    pct = row[f"{segment}_identity"]
                    num = round(float(pct)/100.0*length1)
                    idents_pct.append(f"{pct}%")
                    calls = row[f"{segment}_call"]
                    if calls:
                        length2 = lengthmap[calls.split(",")[0]]
                        idents_nt.append(f"{num}/{length1} of {length2}")
                    else:
                        idents_nt.append("")
                except ValueError:
                    idents_pct.append("")
                    idents_nt.append("")
            results.append({
                "query": row["sequence_id"],
                "cdr_lens": cdrlens,
                "v_call": row["v_call"],
                "v_ident_pct": idents_pct[0],
                "v_ident_nt": idents_nt[0],
                "d_call": row["d_call"],
                "d_ident_pct": idents_pct[1],
                "d_ident_nt": idents_nt[1],
                "j_call": row["j_call"],
                "j_ident_pct": idents_pct[2],
                "j_ident_nt": idents_nt[2]})
        results = sorted(results, key=lambda r: (r["query"]))
        if showtxt:
            show.show_grid(results)
        if output:
            output = Path(output)
            output.parent.mkdir(parents=True, exist_ok=True)
            with open(output, "wt") as f_out:
                writer = DictWriter(f_out, fieldnames=results[0].keys(), lineterminator="\n")
                writer.writeheader()
                writer.writerows(results)
