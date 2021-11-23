#!/usr/bin/env bash

[ -v EXAMPLES ] || EXAMPLES=$(python -c 'import igseq.util; print(igseq.util.DATA)')/examples

QUERY=$EXAMPLES/inputs/igblast/query.fasta

# regular text output
igseq igblast -r rhesus/imgt -Q $QUERY

# AIRR output, saved to a file
igseq igblast -r rhesus/imgt -Q $QUERY -outfmt 19 -out igblast.tsv

# attributes are stored in columns with one row for each query.  For example,
# the V call:
cut -f 10 igblast.tsv
