"""
Run IgBLAST, automatically building databases and the auxiliary data file.

This will gather up whatever combination of reference sequences are given,
build one set of V/D/J database files and one J gene auxiliary data file, and
run igblastn with a query FASTA.

Any command-line arguments not recognized here are passed as-is to the igblastn
command, so you can configure things like the output format and file path.  See
igblastn -help for those options.  Any igblastn argument can be given with two
dashes if needed to force igseq to handle it correctly (for example,
-num_alignments_V will be interprted as -n um_alignments_V, but
--num_alignments_V will work).
"""

import re
import sys
import logging
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
from . import util
from . import aux
from . import vdj

LOGGER = logging.getLogger(__name__)

IGBLASTN = "igblastn"
MAKEBLASTDB = "makeblastdb"

# from our generic names to the IgBLAST names
SPECIESMAP = {
    "rhesus": "rhesus_monkey",
    "human": "human"}

# from synonyms to our generic names
SPECIESOTHER = {
    "human": "human",
    "homosapiens": "human",
    "rhesus": "rhesus",
    "rhesusmonkey": "rhesus"}

def igblast(
    ref_paths, query_path, db_path=None, species=None, extra_args=None, dry_run=False, threads=1):
    """Make temporary IgBLAST DB files and run a query with them.

    ref_paths: list of FASTA files/directories/built-in reference names to use
               for the databases
    query_path: path to FASTA with sequences to check
    db_path: If given, store database files in this directory name and don't
             remove them after running
    species: species name ("human" or "rhesus")
    extra_args: list of arguments to pass to igblastn command.  Must not
                overlap with the arguments set here.
    dry_run: If True, don't actually call any commands or write any files
    threads: number of threads for parallel processing
    """
    LOGGER.info("given ref path(s): %s", ref_paths)
    LOGGER.info("given query path: %s", query_path)
    LOGGER.info("given species: %s", species)
    LOGGER.info("given extra args: %s", extra_args)
    LOGGER.info("given threads: %s", threads)
    LOGGER.info("given db_path: %s", db_path)
    attrs_list = vdj.parse_vdj_paths(ref_paths)
    for attrs in attrs_list:
        LOGGER.info("detected ref path: %s", attrs["path"])
        LOGGER.info("detected ref type: %s", attrs["type"])
    species_det = {attrs.get("species") for attrs in attrs_list}
    species_det = {s for s in species_det if s}
    organism = detect_organism(species_det, species)
    if not dry_run:
        try:
            proc, _ = setup_db_dir_and_igblast(
                [attrs["path"] for attrs in attrs_list],
                organism, query_path, db_path, threads, extra_args,
                capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as err:
            sys.stdout.write(err.stdout)
            sys.stderr.write(err.stderr)
            raise util.IgSeqError("IgBLAST crashed") from err
        sys.stdout.write(proc.stdout)
        sys.stderr.write(proc.stderr)

def detect_organism(species_det, species=None):
    """Determine IgBLAST organism name from multiple species name inputs

    species_det: set of possible species names to use
    species: optional overriding species name

    This includes some fuzzy matching so things like "rhesus_monkey", "RHESUS",
    "rhesus" will all map to "rhesus_monkey" for IgBLAST.
    """
    if not species and not species_det:
        raise util.IgSeqError(
            "species not detected from input.  specify a species manually.")
    if not species and len(species_det) > 1:
        raise util.IgSeqError(
            "multiple species detected from input.  specify a species manually.")
    if not species:
        species = species_det.pop()
        LOGGER.info("detected species: %s", species)
    # match species names if needed
    species_key = re.sub("[^a-z]", "", species.lower())
    if species not in SPECIESMAP and species_key in SPECIESOTHER:
        species_new = SPECIESOTHER[species_key]
        LOGGER.info(
            "detected species as synonym: %s -> %s -> %s", species, species_key, species_new)
        species = species_new
    try:
        organism = SPECIESMAP[species]
    except KeyError as err:
        keys = str(SPECIESMAP.keys())
        raise util.IgSeqError(f"species not recognized.  should be one of: {keys}") from err
    LOGGER.info("detected IgBLAST organism: %s", organism)
    return organism

def setup_db_dir_and_igblast(vdj_ref_paths, organism, query_path,
        db_path=None, threads=1, extra_args=None, **runargs):
    """Run igblastn with automatic database setup"""
    with TemporaryDirectory() as tmp:
        if db_path:
            db_dir = Path(db_path)
            db_dir.mkdir(parents=True, exist_ok=True)
        else:
            db_dir = Path(tmp)
            LOGGER.info("inferred DB directory: %s", db_dir)

        attrs_list = vdj.combine_vdj(vdj_ref_paths, db_dir)
        for segment, attrs in vdj.group(attrs_list).items():
            LOGGER.info("detected %s references: %d", segment, len(attrs))
            if len(attrs) == 0:
                raise util.IgSeqError(f"No references for segment {segment}")
        aux.make_aux_file(db_dir/"J.fasta", db_dir/f"{organism}_gl.aux")
        makeblastdbs(db_dir)
        args = [
            "-germline_db_V", f"{db_dir}/V",
            "-germline_db_D", f"{db_dir}/D",
            "-germline_db_J", f"{db_dir}/J",
            "-query", query_path,
            "-auxiliary_data", f"{db_dir}/{organism}_gl.aux",
            "-organism", organism,
            "-ig_seqtype", "Ig",
            "-num_threads", threads]
        if extra_args:
            # remove any extra - at the start
            extra_args = [re.sub("^--", "-", arg) for arg in extra_args]
            # make sure none of the extra arguments, if there are any, clash
            # with the ones we've defined above.
            args_dashes = {arg for arg in args if str(arg).startswith("-")}
            extra_dashes = {arg for arg in extra_args if str(arg).startswith("-")}
            shared = args_dashes & extra_dashes
            if shared:
                raise util.IgSeqError(f"igblastn arg collision from extra arguments: {shared}")
            args += extra_args
        proc = _run_igblastn(args, **runargs)
        return proc, attrs_list

def makeblastdbs(dir_path):
    """Run makeblastdb for existing V.fasta, D.fasta, J.fasta in a directory."""
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

def _run_igblastn(args, **runargs):
    """Call igblastn with the given list of arguments.

    Any extra keyword arguments are passed to subprocess.run.
    """
    args = [IGBLASTN] + [str(arg) for arg in args]
    LOGGER.info("igblastn command: %s", args)
    LOGGER.info("igblastn subprocess.run args: %s", runargs)
    return subprocess.run(args, **runargs) # pylint: disable=subprocess-run-check
