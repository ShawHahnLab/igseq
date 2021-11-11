#!/usr/bin/env bash

# Take trimmed R1+R2 pairs from a directory, and merge each pair, saving
# outputs in new directory "merged"
igseq merge $EXAMPLES/analysis/trim/211105_M05588_0469_000000000-JWV49 merged
