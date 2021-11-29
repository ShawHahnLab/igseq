#!/usr/bin/env bash

[ -v EXAMPLES ] || EXAMPLES=$(python -c 'import igseq.util; print(igseq.util.DATA)')/examples

# An arbitrary antibody sequence pulled from one of our datasets that looks
# complete and in-frame
QUERY=$EXAMPLES/inputs/igblast/query.fasta

# using the built-in Rhesus germline reference from IMGT and using the default
# text output
igseq igblast -r rhesus/imgt -Q $QUERY

# Trying with AIRR output, saved to a file
igseq igblast -r rhesus/imgt -Q $QUERY -outfmt 19 -out igblast.tsv

# With the AIRR output, attributes are stored in columns with one row for each
# query.  For example, the V call and % identity:
cut -f 10,62 igblast.tsv

# or we can include all rhesus references, but igseq will append suffixes to
# the names for any overlapping loci/segments from the reference FASTA files.
# (IgBLAST can't handle repeated sequence IDs for the V/D/J databases.)
# Here we have two equivalent possible V genes, one from IMGT and one from
# Bernat et al.  2021 (https://doi.org/10.1016/j.immuni.2020.12.018)
igseq igblast -r rhesus -Q $QUERY -outfmt 19 | cut -f 10,62
