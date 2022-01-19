#!/usr/bin/env bash

[ -v EXAMPLES ] || EXAMPLES=$(python -c 'import igseq.util; print(igseq.util.DATA)')/examples

# An arbitrary antibody sequence pulled from one of our datasets that looks
# complete and in-frame
QUERY=$EXAMPLES/inputs/igblast/query.fasta
# A .fastq.gz version, to show off flexibility in query formats
QUERY_FQGZ=$EXAMPLES/inputs/igblast/query.fastq.gz
QUERY_CSV=$EXAMPLES/inputs/igblast/query.csv

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

# The AIRR format is one row per query sequence, but there's also IgBLAST's own
# tabular output format that can give more information on separate lines, like
# additional gene assignments.
# The -num_alignments_V argument clashes with iseq's -n, so we need to use --
# to clarify.  igseq will remove the extra - when calling igblastn.
igseq igblast -r rhesus -Q $QUERY -outfmt 7 --num_alignments_V 5

# like the first example, except giving fastq.gz as the query.  It'll
# automatically be converted to FASTA while being passed to the igblastn
# command.
igseq igblast -r rhesus/imgt -Q $QUERY_FQGZ

# or using tabular (CSV/TSV) input, and specifying which columns have the IDs
# and sequences
igseq igblast -r rhesus/imgt -Q $QUERY_CSV --col-seq-id SeqID --col-seq Seq -outfmt 19 -out igblast2.tsv
