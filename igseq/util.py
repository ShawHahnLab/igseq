"""
Various helper functions.  Not much to see here.
"""

from pathlib import Path
from csv import DictReader, DictWriter

DATA = Path(__file__).parent / "data"
FILES = [path for path in DATA.glob("**/*") if path.is_file()]
READS = ["I1", "R1", "R2"]

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
