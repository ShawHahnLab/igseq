#!/usr/bin/env bash

set -e

export EXAMPLES=$(readlink -f igseq/data/examples)
export HERE=$PWD
export PATH_ORIG=$PATH
export PATH=".:$PATH"
for script in "$EXAMPLES"/*.sh; do
	example=$(basename $script)
	example=${example%.sh}
	echo $example...
	# should trap this for safer cleanup
	mkdir -p examples-tmp
	tmpdir=$(mktemp -d -p examples-tmp ${example}.XXXXXXXX)
	pushd $tmpdir
	ln -s $HERE/tools/example-wrapper.sh igseq
	(set -e; source $script)
	rsyncout=$(rsync -crin "$EXAMPLES/outputs/$example/" .)
	if [[ "$rsyncout" != "" ]]; then
		echo "$rsyncout"
		exit 1
	fi
	popd
done
