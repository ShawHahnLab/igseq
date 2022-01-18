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

import os
import re
import sys
import logging
from contextlib import contextmanager
from pathlib import Path
from subprocess import run, Popen, PIPE, DEVNULL
from tempfile import TemporaryDirectory
from . import util
from . import aux
from . import vdj
from .record import RecordReader, RecordWriter

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
        with setup_db_dir([attrs["path"] for attrs in attrs_list], db_path) as (db_dir, _):
            with run_igblast(db_dir, organism, "-", threads, extra_args, stdin=PIPE) as proc:
                # https://stackoverflow.com/a/66410605
                os.set_blocking(proc.stdout.fileno(), False)
                os.set_blocking(proc.stderr.fileno(), False)
                # read whatever format from the query file, write FASTA to the
                # igblastn proc's stdin
                with RecordReader(query_path, None) as reader, RecordWriter(proc.stdin, "fa") as writer:
                    while proc.poll() is None:
                        try:
                            rec = next(reader)
                        except StopIteration:
                            proc.stdin.close()
                        else:
                            writer.write(rec)
                        for stdout_line in proc.stdout:
                            sys.stdout.write(stdout_line)
                        for stderr_line in proc.stderr:
                            sys.stderr.write(stderr_line)
                if proc.returncode:
                    LOGGER.critical("%s exited with code %d", proc.args[0], proc.returncode)
                    raise util.IgSeqError("IgBLAST crashed")

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

@contextmanager
def run_igblast(
    db_dir, organism, query_path, threads=1, extra_args=None, stdin=None, stdout=PIPE, stderr=PIPE):
    """Start an igblastn process with an already-prepped database directory.

    Run as a context manager; this will yield a Popen object for the running
    process and automatically handle cleanup.

    db_dir: path to directory with V, D, J, and auxiliary files. (See
            setup_db_dir.)
    organism: name of organism as expected by igblastn (e.g. "rhesus_monkey")
    query_path: path to FASTA for query
    threads: number of threads for parallel processing with IgBLAST
    extra_args: list of arguments to pass to igblastn command.  Must not
                overlap with the arguments set here.
    stdin: argument to Popen; default does nothing with stdin
    stdout: argument to Popen; default captures stdout text as a stream
    stderr: argument to Popen; default captures stderr text as a stream
    """
    args = [
        "-germline_db_V", f"{db_dir}/V",
        "-germline_db_D", f"{db_dir}/D",
        "-germline_db_J", f"{db_dir}/J",
        "-auxiliary_data", f"{db_dir}/gl.aux",
        "-organism", organism,
        "-query", query_path,
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
    args = [IGBLASTN] + [str(arg) for arg in args]
    LOGGER.info("igblastn command: %s", args)
    with Popen(args, stdin=stdin, stdout=stdout, stderr=stderr, text=True) as proc:
        yield proc

@contextmanager
def setup_db_dir(vdj_ref_paths, db_path=None):
    """Set up a directory with database files for igblast.

    Run as a context manager; this will yield the path to the prepared
    directory with V, D, and J databases and the auxiliary data file, and will
    automatically handle cleanup if a temporary directory is used (the
    default).  If a db_path is given the files will not be removed afterward.

    vdj_ref_paths: list of V/D/J FASTA files to gather into the db dir.
    db_path: custom path for database output (default: temporary directory)
    """
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
        aux.make_aux_file(db_dir/"J.fasta", db_dir/"gl.aux")
        for segment in ["V", "D", "J"]:
            path_fasta = Path(f"{db_dir}/{segment}.fasta")
            args = [
                MAKEBLASTDB,
                "-dbtype", "nucl",
                "-parse_seqids",
                "-in", path_fasta,
                "-out", path_fasta.parent/path_fasta.stem]
            args = [str(arg) for arg in args]
            run(args, check=True, stdout=DEVNULL)
        yield db_dir, attrs_list
