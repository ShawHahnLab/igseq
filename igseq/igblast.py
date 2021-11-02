"""
Run IgBLAST, automatically building databases and the auxiliary data file.

This will gather up whatever combination of reference sequences are given,
build one set of V/D/J database files and one J gene auxiliary data file, and
run igblastn with a query FASTA.
"""

import re
import logging
import csv
from pathlib import Path
import subprocess
import shlex
from tempfile import TemporaryDirectory
from Bio import SeqIO
from . import util
from . import vdj

LOGGER = logging.getLogger(__name__)

IGBLASTN = "igblastn"
MAKEBLASTDB = "makeblastdb"

# from our generic names to the IgBLAST names
SPECIESMAP = {
    "rhesus": "rhesus_monkey",
    "human": "human"}

# lifted from https://github.com/scharch/SONAR/blob/master/annotate/1.3-finalize_assignments.py
J_MOTIF_REGEX = {
    "JH": "TGGGG",
    "JK": "TT[C|T][G|A]G",
    "JL": "TT[C|T][G|A]G"}

J_MOTIFS = {
    "JH": ["TGGGG"],
    "JK": ["TTCGG", "TTTGG"],
    "JL": ["TTCGG", "TTCTG"]}


# TODO: handle collisions between the entries.  makeblastdb will crash if there are duplicate IDs.
# (duplicate seqs are fine, though; it'll just show multiple matches.)
# maybe:
#   for internal: add _species_ref suffix if there are multiple internal
#   for external: add _path_to_each suffix
# but wait, then the auxiliary data names don't match.

def igblast(ref_paths, query_path, db_path=None, species=None, extra_args=None, dry_run=False, threads=1):
    LOGGER.info("given ref path(s): %s", ref_paths)
    LOGGER.info("given query path: %s", query_path)
    LOGGER.info("given species: %s", species)
    LOGGER.info("given extra args: %s", extra_args)
    LOGGER.info("given threads: %s", threads)
    LOGGER.info("given db_path: %s", db_path)
    vdj_files = vdj.parse_vdj_paths(ref_paths)

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
    if species not in SPECIESMAP:
        raise util.IgSeqError("Species not recognized: %s" % species)
    species_igblast = SPECIESMAP[species]
    LOGGER.info("species name for igblastn: %s", species_igblast)

    vdj_files_grouped = _group_by_segment(vdj_files)
    for key, attrs in vdj_files_grouped.items():
        LOGGER.info("detected %s references: %d", key, len(attrs))
        if len(attrs) == 0:
            raise util.IgSeqError("No references for segment %s" % key)

    if not dry_run:
        with TemporaryDirectory() as tmp:
            if db_path:
                db_dir = Path(db_path)
                db_dir.mkdir(parents=True, exist_ok=True)
            else:
                db_dir = Path(tmp)
            LOGGER.info("inferred DB directory: %s", db_dir)
            for segment, attrs_list in vdj_files_grouped.items():
                fasta = db_dir/f"{segment}.fasta"
                vdj.vdj_gather(attrs_list, fasta)
                if segment == "J":
                    make_aux_file(fasta, db_dir/f"{species_igblast}_gl.aux")
            makeblastdbs(db_dir)
            args = [
                "-germline_db_V", f"{db_dir}/V",
                "-germline_db_D", f"{db_dir}/D",
                "-germline_db_J", f"{db_dir}/J",
                "-query", query_path,
                "-auxiliary_data", f"{db_dir}/{species_igblast}_gl.aux",
                "-organism", species_igblast,
                "-ig_seqtype", "Ig",
                "-num_threads", threads]
            if extra_args:
                args += shlex.split(extra_args)
            _run_igblastn(args)

def makeblastdbs(dir_path):
    for segment in ["V", "D", "J"]:
        _run_makeblastdb(f"{dir_path}/{segment}.fasta")

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

def make_aux_file(j_fasta_in, aux_txt_out):
    """Autogenerate an IgBLAST auxiliary data file from a J FASTA."""
    rows = []
    with open(j_fasta_in) as f_in, open(aux_txt_out, "wt") as f_out:
        writer = csv.writer(f_out, delimiter="\t", lineterminator="\n")
        for record in SeqIO.parse(f_in, "fasta"):
            # NOTE: positions are 0-based!
            # fields are:
            # gene/allele name
            # first coding frame start position
            # chain type
            # CDR3 stop
            # extra bps beyond J coding end
            match = re.match("IG([HKL])J", record.id)
            if not match:
                match = re.match("J([HKL])", record.id)
                if not match:
                    LOGGER.warning("Sequence ID not recognized: %s", record.id)
                    continue
            chain_type = "J" + match.group(1)
            match = _best_motif_match(record.seq, J_MOTIFS[chain_type])
            cdr3_stop = match[0] - 1
            frame = (cdr3_stop + 1) % 3
            extra_bps = (len(record.seq) - frame) % 3
            writer.writerow([
                record.id,
                frame,
                chain_type,
                cdr3_stop,
                extra_bps])

def _best_motif_match(jgene, motifs):
    matches = [_align_motif(jgene, motif) + [motif] for motif in motifs]
    dists = [m[1] for m in matches]
    idx = dists.index(min(dists))
    return matches[idx]

def _align_motif(jgene, motif):
    matches = []
    for idx in range(len(jgene) - len(motif) + 1):
        fragment = jgene[idx:(idx+len(motif))]
        dist = sum([f != m for f, m in zip(fragment, motif)])
        matches.append([idx, dist])
    scores = [m[1] for m in matches]
    idx = scores.index(min(scores))
    return matches[idx]

def _run_igblastn(args, **kwargs):
    """Call igblastn with the given list of arguments.

    Any extra keyword arguments are passed to subprocess.run.
    """
    args = [IGBLASTN] + [str(arg) for arg in args]
    LOGGER.info("igblastn command: %s", args)
    return subprocess.run(args, check=True, **kwargs)

def _group_by_segment(vdj_path_attrs):
    # for prepping one big database from various inputs
    refs = {"V": [], "D": [], "J": []}
    for attrs in vdj_path_attrs:
        refs[attrs["segment"]].append(attrs)
    return refs
