#!/usr/bin/env bash

[ -v EXAMPLES ] || EXAMPLES=$(python -c 'import igseq.util; print(igseq.util.DATA)')/examples

# converting FASTA to FASTA just unwraps it
igseq convert $EXAMPLES/inputs/convert/wrapped.fasta unwrapped.fasta

# or, convert to CSV/TSV
igseq convert $EXAMPLES/inputs/convert/wrapped.fasta unwrapped.csv

# or .fastq.gz to .fasta
igseq convert $EXAMPLES/inputs/convert/unwrapped.fastq.gz unwrapped2.fasta

# a - can be used for stdin/stdout, but the format has to be given explicitly:
igseq convert --input-format fa --output-format fa - - < $EXAMPLES/inputs/convert/wrapped.fasta > unwrapped3.fasta

# other table formats can be converted to FASTA or FASTQ if the column names to
# use are specified.  the default would find sequence_id and sequence columns
# from AIRR:
igseq convert $EXAMPLES/inputs/convert/airr.tsv from_airr.fasta
# or maybe we want the junctions instead:
igseq convert --col-seq junction $EXAMPLES/inputs/convert/airr.tsv from_airr_junctions.fasta
