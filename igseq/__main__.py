"""
Utilities for common IgSeq tasks.
"""

import sys
import argparse
import logging
import shutil
import textwrap
from . import getreads
from . import demux
from . import phix
from . import trim
from . import merge
from . import igblast
from . import vdj_gather
from . import vdj_match
from . import tab2seq
from .show import show_files, list_files
from .util import IgSeqError

LOGGER = logging.getLogger()

def rewrap(txt):
    """Re-wrap text at 80 columns or less, preserving paragraphs."""
    # inspired roughly by what https://github.com/marcelm/cutadapt does
    width = min(80, shutil.get_terminal_size().columns)
    wrap = lambda txt: "\n".join(textwrap.wrap(txt, width=width))
    chunks = txt.strip().split("\n\n")
    return "\n\n".join([wrap(chunk) for chunk in chunks])

def main(arglist=None):
    """Command-line interface.

    args will be parsed from sys.argv by deafult, or a given list.
    """
    parser = __setup_arg_parser()
    if arglist is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(arglist)
    prefix = ""
    if args.dry_run:
        prefix = "[DRYRUN] "
    _setup_log(args.verbose, args.quiet, prefix)
    try:
        args.func(args)
    except IgSeqError as err:
        sys.stderr.write(
            f"\nigseq failed because: {err.message}\n"
            "Considering adding -v or -vv to the command if the problem isn't clear.\n")
        sys.exit(1)

def _main_getreads(args):
    if args.no_counts:
        args.countsfile = None
    getreads.getreads(
        path_input=args.input,
        dir_out=args.outdir,
        path_counts=args.countsfile,
        threads_load=args.threads_load,
        threads_proc=args.threads,
        dry_run=args.dry_run)

def _main_demux(args):
    if args.no_counts:
        args.countsfile = None
    demux.demux(
        paths_input=args.input,
        path_samples=args.samples,
        run_id=args.run,
        dir_out=args.outdir,
        path_counts=args.countsfile,
        path_details=args.details,
        dry_run=args.dry_run)

def _main_phix(args):
    if args.no_counts:
        args.countsfile = None
    phix.phix(
        paths_input=args.input,
        bam_out=args.outfile,
        counts_out=args.countsfile,
        dry_run=args.dry_run,
        threads=args.threads)

def _main_trim(args):
    if args.no_counts:
        args.countsfile = None
    trim.trim(
        paths_input=args.input,
        path_samples=args.samples,
        dir_out=args.outdir,
        path_counts=args.countsfile,
        species=args.species,
        sample_name=args.sample_name,
        min_length=args.min_length,
        quality_cutoff=args.quality_cutoff,
        dry_run=args.dry_run,
        threads=args.threads)

def _main_merge(args):
    if args.no_counts:
        args.countsfile = None
    merge.merge(
        paths_input=args.input,
        dir_out=args.outdir,
        path_counts=args.countsfile,
        dry_run=args.dry_run,
        threads=args.threads)

def _main_show(args):
    show_files(text_items=args.text, force=args.force)

def _main_list(args):
    list_files(text_items=args.text)

def _main_igblast(args):
    igblast.igblast(
        query_path=args.query,
        ref_paths=args.reference,
        db_path=args.database,
        species=args.species,
        extra_args=args.extra_args,
        dry_run=args.dry_run,
        threads=args.threads)

def _main_vdj_gather(args):
    vdj_gather.vdj_gather(
        db_paths=args.input,
        dir_path_out=args.outdir,
        dry_run=args.dry_run)

def _main_vdj_match(args):
    vdj_match.vdj_match(
        ref_paths=args.reference,
        query=args.query,
        output=args.output,
        showtxt=args.show,
        species=args.species,
        dry_run=args.dry_run)

def _main_tab2seq(args):
    tab2seq.tab2seq(
        tab_path_in=args.input,
        seq_path_out=args.output,
        seq_col=args.seq_col,
        seq_id_col=args.seq_id_col,
        seq_desc_col=args.seq_desc_col,
        qual_col=args.seq_qual_col,
        tab_fmt=args.tab_fmt,
        seq_fmt=args.seq_fmt)

def _setup_log(verbose, quiet, prefix):
    # Handle warnings via logging
    logging.captureWarnings(True)
    # Configure the root logger
    # each -v or -q decreases or increases the log level by 10, starting from
    # WARNING by default.
    lvl_current = LOGGER.getEffectiveLevel()
    lvl_subtract = (verbose - quiet) * 10
    verbosity = max(0, lvl_current - lvl_subtract)
    logging.basicConfig(
        format=prefix+"%(levelname)s: %(message)s",
        stream=sys.stderr,
        level=verbosity)

def __setup_arg_parser():
    parser = argparse.ArgumentParser(
        description=rewrap(__doc__),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    #__add_common_args(parser)
    subps = parser.add_subparsers(metavar="", description="igseq "
            "features are split up into these sub-commands.  Call igseq "
            "subcommand --help to get more detailed information on each one.")
    p_get = subps.add_parser("getreads",
        help="get raw read data with Illumina bcl2fastq",
        description=rewrap(getreads.__doc__),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p_demux = subps.add_parser("demux",
        help="demultiplex raw read data into separate samples",
        description=rewrap(demux.__doc__),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p_phix = subps.add_parser("phix",
        help="align unassigned reads post-demux to PhiX genome",
        description=rewrap(phix.__doc__),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p_trim = subps.add_parser("trim",
        help="trim off low-quality and adapter sequences at the ends of the reads with cutadapt",
        description=rewrap(trim.__doc__),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p_merge = subps.add_parser("merge",
        help="merge read pairs with pear",
        description=rewrap(merge.__doc__),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p_igblast = subps.add_parser("igblast",
        help="Run IgBLAST on a set of sequences",
        description=rewrap(igblast.__doc__),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p_vdj_gather = subps.add_parser("vdj-gather",
        help="Gather VDJ sequences into one directory",
        description=rewrap(vdj_gather.__doc__),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p_vdj_match = subps.add_parser("vdj-match",
        help="Find closest-matching germline VDJ sequences",
        description=rewrap(vdj_match.__doc__),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p_tab2seq = subps.add_parser("tab2seq",
        help="Convert CSV/TSV to FASTA/FASTQ",
        description=rewrap(tab2seq.__doc__),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p_show = subps.add_parser("show", help="show builtin reference data")
    p_list = subps.add_parser("list", help="list builtin reference data files")

    __add_common_args(p_get)
    p_get.add_argument("input", help="one Illumina run directory")
    p_get.add_argument("-o", "--outdir", default="", help="Output directory")
    p_get.add_argument("-c", "--countsfile", default="",
        help="file to write read counts to (default: <outdir>/getreads.counts.csv)")
    p_get.add_argument("--no-counts", action="store_true",
        help="don't write a counts file")
    p_get.add_argument("-t", "--threads", type=int, default=28,
        help="number of threads for parallel processing (default: 28)")
    p_get.add_argument("--threads-load", type=int, default=4,
        help="number of threads for parallel loading (default: 4)")
    p_get.set_defaults(func=_main_getreads)

    __add_common_args(p_demux)
    p_demux.add_argument("-s", "--samples", default="metadata/samples.csv",
        help="CSV of sample attributes")
    p_demux.add_argument("-r", "--run",
        help="Run ID (default: parsed from input paths)")
    p_demux.add_argument("-o", "--outdir", default="",
        help="Output directory")
    p_demux.add_argument("-c", "--countsfile", default="",
        help="file to write read counts to (default: <outdir>/demux.counts.csv)")
    p_demux.add_argument("--no-counts", action="store_true",
        help="don't write a counts file")
    p_demux.add_argument("-d", "--details",
        help=".csv.gz file to write with per-read attributes such as assigned barcodes")
    p_demux.add_argument("input", nargs="+",
        help="one directory or individual I1/R1/R2 files in order")
    p_demux.set_defaults(func=_main_demux)

    __add_common_args(p_phix)
    p_phix.add_argument("-o", "--outfile", default="",
        help="Output filename")
    p_phix.add_argument("-c", "--countsfile", default="",
        help="file to write read counts to")
    p_phix.add_argument("--no-counts", action="store_true",
        help="don't write a counts file")
    p_phix.add_argument("-t", "--threads", type=int, default=1,
        help="number of threads for parallel processing (default: 1)")
    p_phix.add_argument("input", nargs="+",
        help="one directory or individual R1/R2 files in order")
    p_phix.set_defaults(func=_main_phix)

    __add_common_args(p_trim)
    p_trim.add_argument("-s", "--samples", default="metadata/samples.csv",
        help="CSV of sample attributes")
    p_trim.add_argument("-o", "--outdir", default="",
        help="Output directory")
    p_trim.add_argument("-c", "--countsfile", default="",
        help="file to write read counts to")
    p_trim.add_argument("--no-counts", action="store_true",
        help="don't write a counts file")
    p_trim.add_argument("-S", "--species", default="rhesus",
        help="species to use for selecting appropriate primer sequences (human or rhesus)")
    p_trim.add_argument("--sample-name",
        help="use this sample name rather than inferring from filenames")
    p_trim.add_argument("--min-length", type=int, default=trim.DEFAULTS["min_length"],
        help="minimum length setting passed to cutadapt "
        f"(default: {trim.DEFAULTS['min_length']})")
    p_trim.add_argument("--quality-cutoff", type=int, default=trim.DEFAULTS["quality_cutoff"],
        help="quality cutoff setting passed to cutadapt "
        f"(default: {trim.DEFAULTS['quality_cutoff']})")
    p_trim.add_argument("-t", "--threads", type=int, default=1,
        help="number of threads for parallel processing (default: 1)")
    p_trim.add_argument("input", nargs="+",
        help="one directory or individual R1/R2 files in order")
    p_trim.set_defaults(func=_main_trim)

    __add_common_args(p_merge)
    p_merge.add_argument("-o", "--outdir", default="",
        help="Output directory")
    p_merge.add_argument("-c", "--countsfile", default="",
        help="file to write read counts to")
    p_merge.add_argument("--no-counts", action="store_true",
        help="don't write a counts file")
    p_merge.add_argument("-t", "--threads", type=int, default=1,
        help="number of threads for parallel processing (default: 1)")
    p_merge.add_argument("input", nargs="+",
        help="one directory or individual R1/R2 files in order")
    p_merge.set_defaults(func=_main_merge)

    __add_common_args(p_show)
    p_show.add_argument("text", nargs="+",
        help="partial filename to show")
    p_show.add_argument("-f", "--force", action="store_true",
        help="force display of possibly non-text files")
    p_show.set_defaults(func=_main_show)

    __add_common_args(p_list)
    p_list.add_argument("text", nargs="*", help="partial filename to list")
    p_list.set_defaults(func=_main_list)

    __add_common_args(p_igblast)
    p_igblast.add_argument("-Q", "--query", required=True,
        help="query FASTA")
    p_igblast.add_argument("-r", "--reference", nargs="+",
            help="one or more FASTA/directory/builtin names pointing to V/D/J FASTA files")
    p_igblast.add_argument("-d", "--database",
            help="optional persistent database directory name (default: use temp directory)")
    p_igblast.add_argument("-S", "--species",
            help="species to use (human or rhesus).  Default: infer from database if possible")
    p_igblast.add_argument("-e", "--extra-args",
            help="extra igblastn arguments as a single string")
    p_igblast.add_argument("-t", "--threads", type=int, default=1,
        help="number of threads for parallel processing (default: 1)")
    p_igblast.set_defaults(func=_main_igblast)

    __add_common_args(p_vdj_gather)
    p_vdj_gather.add_argument("input", nargs="+",
        help="one directory with one or more each of V, D, J FASTA files.")
    p_vdj_gather.add_argument("-o", "--outdir",
        help="Output directory")
    p_vdj_gather.set_defaults(func=_main_vdj_gather)

    __add_common_args(p_vdj_match)
    p_vdj_match.add_argument("-r", "--reference", nargs="+",
        help="one or more FASTA/directory/builtin names pointing to V/D/J FASTA files")
    p_vdj_match.add_argument("-Q", "--query", required=True,
        help="query FASTA")
    p_vdj_match.add_argument("-S", "--species",
            help="species to use (human or rhesus).  Default: infer from database if possible")
    p_vdj_match.add_argument("-o", "--output",
        help="Output filename")
    p_vdj_match.add_argument("--show", action=argparse.BooleanOptionalAction,
        help="Explicitly enable/disable showing the results directly on standard output "
        "(default: disabled if using file output, enabled otherwise)")
    p_vdj_match.set_defaults(func=_main_vdj_match)

    __add_common_args(p_tab2seq)
    p_tab2seq.add_argument("input",
        help="one CSV or TSV file path, or a literal '-' for standard input")
    p_tab2seq.add_argument("output",
        help="one FASTA or FASTQ file path, or a literal '-' for standard output")
    p_tab2seq.add_argument("--seq-col", required=True,
        help="name of table column containing sequences")
    p_tab2seq.add_argument("--seq-id-col", required=True,
        help="name of table column containing sequence IDs")
    p_tab2seq.add_argument("--seq-desc-col",
        help="name of table column containing sequence descriptions (optional)")
    p_tab2seq.add_argument("--seq-qual-col",
        help="name of table column containing sequence quality "
        "scores as PHRED+33 (for FASTQ output only)")
    p_tab2seq.add_argument("--tab-fmt",
            help="Format of input: tsv or csv.  "
            "default is detected from input filename if possible")
    p_tab2seq.add_argument("--seq-fmt",
            help="Format of output: fasta or fastq.  "
            "default is detected from output filename if possible")
    p_tab2seq.set_defaults(func=_main_tab2seq)

    return parser

def __add_common_args(obj):
    obj.add_argument("-v", "--verbose", action="count", default=0,
            help="Increment verbosity.  This can be specified multiple times.")
    obj.add_argument("-q", "--quiet", action="count", default=0,
            help="Decrement verbosity.  This can be specified multiple times.")
    obj.add_argument("-n", "--dry-run", action="store_true", help="Don't "
            "actually write any files.  This is useful with --verbose to see "
            "what would be happening without actually doing anything.")

if __name__ == "__main__":
    main()
