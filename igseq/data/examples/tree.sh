#!/usr/bin/env bash

[ -v EXAMPLES ] || EXAMPLES=$(python -c 'import igseq.util; print(igseq.util.DATA)')/examples

FASTA=$EXAMPLES/inputs/tree/seqs.fasta
FASTA2=$EXAMPLES/inputs/tree/seqs2.fasta
WK16=$EXAMPLES/inputs/tree/wk16.txt
WK20=$EXAMPLES/inputs/tree/wk20.txt

# Very basic example: Take FASTA, align with muscle if needed, and create
# newick tree with fasttree.
igseq tree $FASTA example.tree

# NEXUS is also supported.  By default it just wraps the same newick data.
igseq tree $FASTA example.nex

# But NEXUS can also support formatting info, which we can use to add
# per-sequence color codes with a set-based system.  For example, grouping
# sequences into sets based on a pattern in the sequence IDs.  This internally
# tracks set names (wk16, wk20) even though nothing is done with those
# explicitly.
igseq tree $FASTA -P '^wk[0-9]+' example_pattern.nex

# For NEXUS output FigTree's formatting options can be specified as key=value
# pairs and they'll be included in a figtree block at the end of the file.
igseq tree $FASTA -P '^wk[0-9]+' -F 'branchLabels.fontSize=8' example_figtree.nex

# We can also use a newick tree file (like created above) as input and create a
# color-coded NEXUS file from that.
igseq tree example.tree -P '^wk[0-9]+' example_pattern_from_newick.nex

# Or, define sets by lists of sequence IDs.
igseq tree $FASTA -L $WK16 -L $WK20 example_lists.nex

# We can override colors using the set names found in the IDs based on that
# pattern.
igseq tree $FASTA -P '^wk[0-9]+' -C 'wk16=#BB0000' -C 'wk20=#6600EE' example_colors.nex

# Or, we can give multiple FASTA as input, in which case sequences will be
# automatically placed in sets based on which are in which input file.  The IDs
# and sequences have to be consistent between files (ignoring gaps and
# uppercase/lowercase).
igseq tree $FASTA $FASTA2 example_multi.nex
