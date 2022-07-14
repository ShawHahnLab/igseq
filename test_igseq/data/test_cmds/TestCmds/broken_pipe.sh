#!/usr/bin/env bash
# For #20
# Run a big igblast query but only read one line from stdout
set -e -o pipefail
python -m igseq igblast -r rhesus/imgt -Q <(python -m igseq show IGHV fasta 2> /dev/null) --input-format fa 2> /dev/null | head -n 1 > /dev/null
