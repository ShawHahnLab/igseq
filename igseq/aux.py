"""
Helper functions for IgBLAST auxiliary data files

These optional files provide IgBLAST with details like CDR3 3' boundaries.
https://ncbi.github.io/igblast/cook/How-to-set-up.html
"""

import re
import csv
import logging
from Bio import SeqIO

LOGGER = logging.getLogger(__name__)

# Adapter from SONAR approach but using alignments instead of regex
# See: https://github.com/scharch/SONAR/blob/master/annotate/1.3-finalize_assignments.py
J_MOTIFS = {
    "JH": ["TGGGG"],
    "JK": ["TTCGG", "TTTGG"],
    "JL": ["TTCGG", "TTCTG"]}

def make_aux_file(j_fasta_in, aux_txt_out):
    """Autogenerate an IgBLAST auxiliary data file from a J FASTA."""
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
