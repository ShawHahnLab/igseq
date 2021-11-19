#!/usr/bin/env bash

set -e

export EXAMPLES_HERE=$(readlink -f igseq/data/examples)
export HERE=$PWD
for script in "$EXAMPLES_HERE"/*.sh; do
	example=$(basename $script)
	example=${example%.sh}
	echo $example...
	# should trap this for safer cleanup
	mkdir -p examples-tmp
	tmpdir=$(mktemp -d -p examples-tmp ${example}.XXXXXXXX)
	pushd $tmpdir
	ln -s $HERE/tools/example-wrapper.sh igseq
	export PATH_ORIG=$PATH
	export PATH=".:$PATH"
	$script
	rsyncout=$(rsync -crin "$EXAMPLES_HERE/outputs/$example/" .)
	if [[ "$rsyncout" != "" ]]; then
		echo "$rsyncout"
		exit 1
	fi
	popd
done
