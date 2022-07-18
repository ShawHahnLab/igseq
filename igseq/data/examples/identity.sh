#!/usr/bin/env bash

[ -v EXAMPLES ] || EXAMPLES=$(python -c 'import igseq.util; print(igseq.util.DATA)')/examples

# A few example query sequences derived from the igblast example
QUERY=$EXAMPLES/inputs/identity/query.fasta
# AIRR TSV table for those sequences
QUERY_AIRR=$EXAMPLES/inputs/identity/query_airr.tsv
# A contrived reference (slightly altered version of the first query)
REF=$EXAMPLES/inputs/identity/ref.fasta
# Just the junction region from that sequence, as NT and AA
REF_JUNCTION=$EXAMPLES/inputs/identity/ref_junction.fasta
REF_JUNCTION_AA=$EXAMPLES/inputs/identity/ref_junction_aa.fasta
# AIRR TSV for the reference
REF_AIRR=$EXAMPLES/inputs/identity/ref_airr.tsv

# Check identity of each query to the one reference, with each row in the CSV
# output giving query ID, reference ID, and fraction of matching bases
igseq identity -r $REF $QUERY idents.csv
# use "igseq show" command for a human-readable view of the CSV file
igseq show idents.csv

# With no reference specified the first query sequence will be used as the reference
igseq identity $QUERY idents_first_query.csv

# To compare all-to-all, explicitly give the query path as the reference too
igseq identity -r $QUERY $QUERY idents_all.csv

# We can use a "-" to read from/to standard input and/or output, but we need to
# specify the input format for that case
igseq identity -r $REF --input-format fa - - < $QUERY > idents_stdout.csv

# We can work with the formats supported by the convert command, like AIRR TSV
# for example.  As with that command, the defaults are to use the sequence_id
# and sequence columns.
igseq identity -r $REF $QUERY_AIRR idents_from_airr.csv

# We can specify other columns to use instead
igseq identity -r $REF_JUNCTION --col-seq junction $QUERY_AIRR idents_from_airr_junctions.csv

# AA works too
igseq identity -r $REF_JUNCTION_AA --col-seq junction_aa $QUERY_AIRR idents_from_airr_junctions_aa.csv

# We can also give the reference in a format other than FASTA, though the same
# column settings are used for query and reference
igseq identity -r $REF_AIRR --col-seq junction_aa $QUERY_AIRR idents_from_airr_junctions_aa_2.csv
