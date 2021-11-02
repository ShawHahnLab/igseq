"""
Show contents of files built into the igseq package.
"""

import os
import sys
from pathlib import Path
from csv import DictReader
from .util import FILES

def show_files(text_items, force=False):
    other_paths = []
    for item in text_items:
        if Path(item).is_file():
            other_paths.append(item)
            show_file(item, force)
    text_items = [item for item in text_items if item not in other_paths]
    if text_items:
        for path in FILES:
            if all([text in str(path) for text in text_items]):
                show_file(path, force)

def list_files(text_items):
    for path in FILES:
        if all([text in str(path) for text in text_items]):
            print(path)

def show_file(path, force=False):
    path = Path(path)
    if path.suffix in [".csv"]:
        show_csv(path)
    elif path.suffix in [".txt"]:
        show_text(path)
    elif path.suffix in [".yml", ".yaml"]:
        show_text(path)
    elif path.suffix in [".fasta", ".fa", ".fastq", ".fq"]:
        show_text(path)
    else:
        if force:
            show_raw(path)
        else:
            sys.stderr.write(
                f"{path} might not display well.  "
                "Use --force if you're sure.\n")

def show_csv(path):
    with open(path, encoding="UTF8") as f_in:
        reader = DictReader(f_in)
        show_grid(list(reader))

def show_grid(grid):
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
            print(row[key].rjust(widths[key]+1), end="")
        print("")

def show_text(path):
    with open(path, encoding="UTF8") as f_in:
        sys.stdout.write(f_in.read())

def show_raw(path):
    # https://stackoverflow.com/a/54073813/4499968
    with os.fdopen(sys.stdout.fileno(), "wb", closefd=False) as stdout, open(path, "rb") as f_in:
        stdout.write(f_in.read())
        stdout.flush()
