"""
Various helper functions.  Not much to see here.
"""

import re
from pathlib import Path
from csv import DictReader, DictWriter
from Bio.Seq import Seq

DATA = Path(__file__).parent / "data"
FILES = [path for path in DATA.glob("**/*") if path.is_file()]
READS = ["I1", "R1", "R2"]

def __load_primers(path):
    with open(path, encoding="ascii") as f_in:
        reader = DictReader(f_in)
        return list(reader)

def __load_barcodes(path):
    with open(path, encoding="ascii") as f_in:
        reader = DictReader(f_in)
        return list(reader)

BARCODES = __load_barcodes(DATA / "barcodes.csv")
PRIMERS = __load_primers(DATA / "primers.csv")
# P5 Sequencing Primer Site
# comes just before forward barcode
P5SEQ = "TCTTTCCCTACACGACGCTCTTCCGATCT"
# 5' RACE Anchor
# comes just after forward barcode
ANCHOR5P = "AAGCAGTGGTATCAACGCAGAGTACATGGG"

def load_samples(path_samples):
    """"Load a CSV of sample info into a per-sample dictionary."""
    samples = {}
    with open(path_samples, newline="", encoding="ASCII") as f_in:
        reader = DictReader(f_in)
        for row in reader:
            samples[row["Sample"]] = row
    return samples

def save_counts(path_counts, counts):
    """Save seq count info (lists of Category/Item/NumSeqs dicts) to CSV file."""
    path_counts = Path(path_counts)
    path_counts.parent.mkdir(parents=True, exist_ok=True)
    with open(path_counts, "wt", encoding="ASCII") as f_out:
        writer = DictWriter(f_out, fieldnames=["Category", "Item", "NumSeqs"], lineterminator="\n")
        writer.writeheader()
        for row in counts:
            writer.writerow(row)

def parse_fqgz_paths(paths, keys=None, txt=""):
    """Take list of paths and match up with I1/R1/R2 files as a dictionary.

    If the list has one item and it's a directory, try to match by filename.
    If the list contains one or more files, pair them up with the keys in
    order.

    keys: what items to match.  If None, will be I1/R1/R2.
    txt: substring to require in fastq.gz filenames.  If empty, all fastq.gz
         files will be considered.
    """
    if not keys:
        keys = READS.copy()
    result = {}
    if len(paths) == 1:
        path = Path(paths[0])
        if path.is_dir():
            files = list(path.glob("*.fastq.gz"))
            for key in keys:
                for fpath in files:
                    if txt in fpath.name and key in fpath.name:
                        result[key] = Path(fpath)
        else:
            raise ValueError
    elif len(paths) == len(keys):
        for key, path in zip(keys, paths):
            result[key] = Path(path)
    else:
        raise ValueError
    if len(set(result.values())) != len(keys):
        raise ValueError
    return result

def default_path(paths, name1, name2=""):
    parent = common_parent(paths)
    return parent.parent.parent / name1 / parent.name / name2

def common_parent(paths):
    """Get the parent directory in common for all the given paths.

    For example:

    path/to/something/dir1/subdir1
    path/to/something/dir2/subdir2

    would give:

    path/to/something
    """
    # work with either dict of paths or list of paths
    try:
        paths = [Path(path) for path in paths.values()]
    except AttributeError:
        paths = [Path(path) for path in paths]
    parents = [reversed(path.parents) for path in paths]
    parent = None
    for level in zip(*parents):
        if len(set(level)) > 1:
            break
        parent = level[0]
    return parent

def parse_quals(txt, shift=33):
    """Parse FASTQ-encoded Phred scores from ASCII text."""
    return [ord(letter)-shift for letter in txt]

def revcmp(record):
    """Reverse complement a SeqRcord, keeping ALL metadata.
    If a string is given, reverse-complement that.
    """
    try:
        return record.reverse_complement(
            id=True, name=True, description=True, features=True,
            annotations=True, letter_annotations=True, dbxrefs=True)
    except AttributeError:
        return Seq(record).reverse_complement()

def assign_barcode_seqs(samples):
    # copy sample attributes so we can safely add barcode sequences to our copy
    # TODO handle unexpected keys
    samples = {key: samples[key].copy() for key in samples}
    # lookup table of direction+barcode ID, with ID as integer.
    bc_lut = {}
    for attrs in BARCODES:
        bc_lut[(attrs["Direction"], int(attrs["BC"]))] = attrs["Seq"]
    # remove a prefix on the text, if there is one, and cast to int.  Or if
    # that doesn't work assume it already is.
    def bcmunge(txt):
        try:
            return int(re.sub("^.*_", "", txt))
        except TypeError:
            return txt
    for attrs in samples.values():
        attrs["BarcodeFwdSeq"] = bc_lut[("F", bcmunge(attrs["BarcodeFwd"]))]
        attrs["BarcodeRevSeq"] = bc_lut[("R", bcmunge(attrs["BarcodeRev"]))]
    return samples
