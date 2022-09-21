#!/usr/bin/env bash

[ -v EXAMPLES ] || EXAMPLES=$(python -c 'import igseq.util; print(igseq.util.DATA)')/examples

# A simple MSA case
igseq msa $EXAMPLES/inputs/msa/seqs.fasta seqs.aln.fasta
