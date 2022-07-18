"""
List or show contents of files.

The list and show commands can display files built into the igseq package, plus
some recognized file types from outside the package.  Give one or more filename
fragments (for builtin files) or complete filenames (for all others).  Builtin
files will be filtered to those matching ALL fragments, and external files will
include exact matches only.

    igseq list                list all data files built into the package
    igseq list germ           list germline reference files
    igseq list germ V         list germline ref files with V in their names
    igseq show primers        show igseq's primers.csv with nice formatting
    igseq show IGHV.fasta     show all builtin IGHV FASTAs concatenated
    igseq show some/file.csv  show an arbitrary CSV file
"""

import os
import sys
import logging
from pathlib import Path
from csv import DictReader
import newick
from .util import FILES

LOGGER = logging.getLogger(__name__)

def show_files(text_items, force=False):
    """Format and print matching files to stdout."""
    if not text_items or not [item for item in text_items if item]:
        LOGGER.warning("No items to show for %s", text_items)
        return
    text_items = [str(item) for item in text_items]
    other_paths = []
    for item in text_items:
        if Path(item).is_file():
            other_paths.append(item)
            show_file(item, force)
    text_items = [item for item in text_items if item not in other_paths]
    internal_paths = []
    if text_items:
        for path in FILES:
            if all(text in str(path) for text in text_items):
                show_file(path, force)
                internal_paths.append(path)
    if not other_paths and not internal_paths:
        LOGGER.warning("No files found for %s", text_items)

def list_files(text_items):
    """List all matching builtin data files on stdout."""
    text_items = [str(item) for item in text_items]
    for path in FILES:
        if all(text in str(path) for text in text_items):
            print(path)

def show_file(path, force=False):
    """Infer filetype from extention and pretty-print to stdout."""
    path = Path(path)
    if path.suffix in [".csv"]:
        show_csv(path)
    elif path.suffix in [".tsv", ".tab"]:
        show_csv(path, delimiter="\t")
    elif path.suffix in [".txt"]:
        show_text(path)
    elif path.suffix in [".yml", ".yaml"]:
        show_text(path)
    elif path.suffix in [".fasta", ".fa", ".fastq", ".fq"]:
        show_text(path)
    elif path.suffix in [".tree", ".newick"]:
        show_tree(path)
    else:
        if force:
            show_raw(path)
        else:
            sys.stderr.write(
                f"{path} might not display well.  "
                "Use --force if you're sure.\n")

def show_csv(path, **kwargs):
    """Pretty-print CSV file to stdout."""
    with open(path, encoding="UTF8") as f_in:
        reader = DictReader(f_in, **kwargs)
        show_grid(list(reader))

def show_grid(grid):
    """Pretty-print list of dictionaries to stdout."""
    if not grid:
        return
    fieldnames = grid[0].keys()
    widths = {}
    for key in fieldnames:
        widths[key] = max(len(key), widths.get(key, 0))
    for row in grid:
        for key in row:
            widths[key] = max(len(str(row[key])), widths.get(key, 0))
    for key in fieldnames:
        print(key.rjust(widths[key]+1), end="")
    print("")
    for row in grid:
        for key in fieldnames:
            print(str(row[key]).rjust(widths[key]+1), end="")
        print("")

def show_tree(path):
    """Print a newick tree file to stdout."""
    tree = newick.read(path)[0]
    sys.stdout.write(tree.ascii_art())
    sys.stdout.write("\n")

def show_text(path):
    """Print a plaintext file to stdout."""
    with open(path, encoding="UTF8") as f_in:
        sys.stdout.write(f_in.read())

def show_raw(path):
    """Print a binary file to stdout."""
    # https://stackoverflow.com/a/54073813/4499968
    with os.fdopen(sys.stdout.fileno(), "wb", closefd=False) as stdout, open(path, "rb") as f_in:
        stdout.write(f_in.read())
        stdout.flush()
