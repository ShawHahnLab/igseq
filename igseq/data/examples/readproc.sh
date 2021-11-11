#!/usr/bin/env bash

EXAMPLES=$(python -c 'import igseq.util; print(igseq.util.DATA)')/examples

# A whole start-to-finish example for one run.
#
# Each step generally creates a set of output files in one directory with one
# or more counts.csv files summarizing read counts.  For each command the
# output directory is left implicit, but given as the input for the next
# command.

# create R1/R2/I1 fastq.gz files using bcl2fastq and save to
# the default output, analysis/reads/211105_M05588_0469_000000000-JWV49
igseq getreads $EXAMPLES/analysis/runs/211105_M05588_0469_000000000-JWV49

# separate reads from that run per-sample using the per-run per-sample barcode
# information from a samples.csv file.
igseq demux -s $EXAMPLES/samples.csv analysis/reads/211105_M05588_0469_000000000-JWV49

# some reads don't have a recognized barcode combination so they aren't
# assigned to any sample.  These *should* be the PhiX phage reads that were
# added to the run, but we can double-check that.  This command will find the
# unasigned.R1.fastq.gz and unassigned.R2.fastq.gz files and map them to the
# PhiX genome.
igseq phix analysis/demux/211105_M05588_0469_000000000-JWV49

# Take the reads for each sample and trim off primer/adapter sequences and
# lower-quality 3' regions.  This needs to know the antibody chain type,
# species, and barcodes for each sample via the samples.csv file to figure out
# what sequences to trim.
igseq trim -s $EXAMPLES/samples.csv analysis/demux/211105_M05588_0469_000000000-JWV49

# Merge the trimmed R1+R2 pairs for each sample
igseq merge analysis/trim/211105_M05588_0469_000000000-JWV49

# From here we can take the merged reads and give them to IgDiscover, SONAR,
# or anything else.

# As one example, use the IgBLAST wraper command to automatically set up IgBLAST database
# files for all builtin rhesus gene sequences and make an AIRR-format TSV table
# (format 19 in IgBLAST jargon)
igseq igblast -r rhesus -Q analysis/merge/211105_M05588_0469_000000000-JWV49/wk12.fastq.gz -outfmt 19 -out wk14_airr.tsv
