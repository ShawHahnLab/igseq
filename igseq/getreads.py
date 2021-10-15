"""
Get raw reads from an Illumina run directory.

This is a wrapper around Illumina's bcl2fastq program, just with some automatic
handling of the weirdness of the IgSeq protocol.  A dummy sample sheet is used
to try to ensure that all reads are placed in the Undetermined files (since
demultiplexing must be performed separately anyway) and index reads are enabled
so they can be used later during demultiplexing.

This step can be skipped if you already have all of your reads in single trio
of I1/R1/R2 fastq.gz files.
"""

import gzip
import logging
from tempfile import NamedTemporaryFile
from pathlib import Path
from subprocess import Popen

LOGGER = logging.getLogger(__name__)

BCL2FASTQ = "bcl2fastq"

def getreads(path_input, dir_out, threads_load=1, threads_proc=1, dry_run=False):
    """Get reads directly from Illumina run directory.

    This is a wrapper around Illumina's bcl2fastq with its built-in
    demultiplexing features (mostly) disabled and the I1 output enabled.  All
    reads should end up in the I1/R1/R2 Undetermined files in the output
    directory.

    There's a chance reads will be written to a dummy sample's output files due
    to bcl2fastq's insistence of having at least one sample defined in order to
    write an I1 file.  We use all Ns for that dummy barcode, so if that happens
    it most likely indicates a quality issue.  If any reads are written for
    that sample, a message will be logged complaining about it.
    """
    path_input = Path(path_input)
    if not dir_out:
        name = path_input.name
        dir_out = Path("analysis/reads" / name)
    else:
        dir_out = Path(dir_out)
    LOGGER.info("input: %s", path_input)
    LOGGER.info("output: %s", dir_out)
    if not dry_run:
        # We don't want any samples defined since we're just using the program
        # to convert .bcl to .fastq.gz (no demultiplexing).  We also don't want
        # any weirdness coming from the sample sheet like the ReverseComplement
        # setting.  Apparently the sample sheet can be completely empty as far
        # as bcl2fastq is concerned, but if there's isn't at least one sample
        # defined it won't create the Undetermined I1 file.
        # So we'll make a stub sample sheet just with that minimum info.
        with NamedTemporaryFile("wt") as sample_sheet:
            sample_sheet.write("[Data]\nSample_ID,index\nsample1,NANANANA\n")
            sample_sheet.flush()
            args = [
                BCL2FASTQ,
                "--runfolder-dir", str(path_input),
                "--output-dir", str(dir_out),
                # By default it'll try to write some files to the input run directory.
                # we'll send those to the output directory instead.
                "--interop-dir", str(dir_out / "InterOp"),
                # We want the I1 file, so we enable that here.
                "--create-fastq-for-index-reads",
                # don't bother doing a fuzzy match on any supposed barcodes
                "--barcode-mismatches", "0",
                "--sample-sheet", sample_sheet.name,
                # parallel processing during loading can help a bit in my tests
                "--loading-threads", str(threads_load),
                # parallel processing does *not* help during the bcl2fastq
                # demultiplexing step, go figure, when we don't have any
                # demultiplexing to perform here
                "--demultiplexing-threads", "1",
                # parallel processing in the processing step helps quite a bit
                "--processing-threads", str(threads_proc),
                # help text says "this must not be higher than number of
                # samples"
                "--writing-threads", "1",
                # it's pretty verbose by default so we'll set a higher log level
                "--min-log-level", "WARNING"
                ]
            LOGGER.info("bcl2fastq command: %s", args)
            with Popen(args) as proc:
                proc.wait()
        notempties = []
        # peek into each fastq.gz and see if it has anything
        for path in dir_out.glob("*.fastq.gz"):
            with gzip.open(path, "rt") as f_in:
                if f_in.read(1):
                    notempties.append(path)
        # We should have three (I1/R1/R2 Undetermined) non-empty fastq.gz
        # files.  If there are more, maybe some reads got dumped to other
        # "samples."  If there are fewer...?
        if len(notempties) > 3:
            LOGGER.warning("more .fastq.gz outputs than expected.")
        elif len(notempties) < 3:
            LOGGER.error("Missing .fastq.gz outputs")

#def _aggregate_read_trios(entries):
#    extras = []
#    for key, val in entries.items():
#        if key == "Undetermined":
#            outputs = [val[read]["path"] for read in ["I1", "R1", "R2"]]
#        else:
#            extras.append([val[read]["path"] for read in ["I1", "R1", "R2"]])
#    for idx, output in enumerate(outputs):
#        with open(output, "wb") as f_out:
#            for trio in extras:
#                with open(trio[idx], "rb") as f_in:
#                    f_out.write(f_in.read())
#                trio[idx].unlink()
#
#def parse_fqgz_paths(
#        dirpath, pattern=r"([-A-za-z0-9]+)_S([0-9]+)_L([0-9]+)_([IR][12])_001\.fastq\.gz"):
#    dirpath = Path(dirpath)
#    result = []
#    for path in dirpath.glob("*.fastq.gz"):
#        match = re.fullmatch(pattern, path.name)
#        if match:
#            result.append({
#                "sample": match.group(1),
#                "sample_num": int(match.group(2)),
#                "lane": int(match.group(3)),
#                "read": match.group(4),
#                "path": path})
#    result_groups = {}
#    for entry in result:
#        if not entry["sample"] in result_groups:
#            result_groups[entry["sample"]] = {}
#        result_groups[entry["sample"]][entry["read"]] = result
#    return result_groups
