#!/usr/bin/env bash

set -e

EXAMPLES=$(readlink -f igseq/data/examples)
for script in "$EXAMPLES"/*.sh; do
	example=$(basename $script)
	example=${example%.sh}
	echo $example...
	# should trap this for safer cleanup
	mkdir -p examples-tmp
	tmpdir=$(mktemp -d -p examples-tmp ${example}.XXXXXXXX)
	pushd $tmpdir
	$script
	rsyncout=$(rsync -crin "$EXAMPLES/outputs/$example/" .)
	if [[ "$rsyncout" != "" ]]; then
		echo "$rsyncout"
		exit 1
	fi
	popd
done
