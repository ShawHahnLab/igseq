"""
Gather VDJ sequences into one directory.

This is handy for starting off an IgDiscover run from one of the built-in
references, or manually prepping an IgBLAST database.  If multiple FASTA files
look like they match the same locus and segment (for example two IGHV files
from different references) a suffix is added to each sequence ID for those to
differentiate them.
"""

import logging
from pathlib import Path
from . import vdj

LOGGER = logging.getLogger(__name__)

def vdj_gather(ref_paths, dir_path_out, dry_run=False):
    """Gather V/D/J FASTA into one directory as V.fasta, D.fasta, J.fasta."""
    LOGGER.info("given ref path(s): %s", ref_paths)
    LOGGER.info("given output: %s", dir_path_out)
    vdj_files = vdj.parse_vdj_paths(ref_paths)
    vdj_files_grouped = vdj.group(vdj_files)
    for segment, attrs_list in vdj_files_grouped.items():
        LOGGER.info("detected %s FASTA: %d", segment, len(attrs_list))
        if not dry_run:
            vdj.combine_vdj(attrs_list, Path(dir_path_out)/f"{segment}.fasta")
