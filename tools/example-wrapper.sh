#!/usr/bin/env bash

# To be called by tools/check-examples.sh to intercept getreads command if
# needed

export PATH=$PATH_ORIG
if [[ $1 == getreads && -v SKIPBCL2FASTQ ]]; then
	echo "SKIPBCL2FASTQ set; Skipping igseq getreads"
	rundir=$(basename "$2")
	dir_out=analysis/reads/$rundir
	mkdir -p $dir_out
	for read in I1 R1 R2; do
		cp $EXAMPLES_HERE/runs/$rundir.${read}.fastq.gz $dir_out
	done
	cp $EXAMPLES_HERE/outputs/readproc/$dir_out/getreads.counts.csv $dir_out
else
	igseq "$@"
fi
