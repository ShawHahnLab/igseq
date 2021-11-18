#!/usr/bin/env bash

# A whole start-to-finish example for one run.
#
# Each step generally creates a set of output files in one directory with one
# or more counts.csv files summarizing read counts.  For each command the
# output directory is left implicit, but given as the input for the next
# command.

EXAMPLES=$(python -c 'import igseq.util; print(igseq.util.DATA)')/examples

# Extract an example Illumina run directory with just a handful of reads
tar xzf $EXAMPLES/runs/YYMMDD_M05588_0232_000000000-CLL8M.tgz

# create R1/R2/I1 fastq.gz files using bcl2fastq and save to
# the default output, analysis/reads/YYMMDD_M05588_0232_000000000-CLL8M
igseq getreads YYMMDD_M05588_0232_000000000-CLL8M

# separate reads from that run per-sample using the per-run per-sample barcode
# information from a samples.csv file.
igseq demux -s $EXAMPLES/samples.csv analysis/reads/YYMMDD_M05588_0232_000000000-CLL8M

# some reads don't have a recognized barcode combination so they aren't
# assigned to any sample.  These *should* be the PhiX phage reads that were
# added to the run, but we can double-check that.  This command will find the
# unasigned.R1.fastq.gz and unassigned.R2.fastq.gz files and map them to the
# PhiX genome.
igseq phix analysis/demux/YYMMDD_M05588_0232_000000000-CLL8M

# Take the reads for each sample and trim off primer/adapter sequences and
# lower-quality 3' regions.  This needs to know the antibody chain type,
# species, and barcodes for each sample via the samples.csv file to figure out
# what sequences to trim.
igseq trim -s $EXAMPLES/samples.csv -S rhesus analysis/demux/YYMMDD_M05588_0232_000000000-CLL8M

# Merge the trimmed R1+R2 pairs for each sample
# I've left out reads for the wk24H sample so there's a warning here about
# skipping read merging via pear (which crashes with empty input files) and
# instead empty outputs are created for that one.
igseq merge analysis/trim/YYMMDD_M05588_0232_000000000-CLL8M

# From here we can take the merged reads for sample "wk48H" and give them to
# IgDiscover, SONAR, or anything else.
