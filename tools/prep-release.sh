#!/usr/bin/env bash

version=$(python -c 'from igseq.version import __version__; print(__version__)')
date=$(date +%Y-%m-%d)
echo Version: $version
echo
sed -i 's/{% set version = ".*" %}/{% set version = "'$version'" %}'/ conda/meta.yaml
sed -i "s/# dev$/# $version - $date/" CHANGELOG.md
echo conda/meta.yaml vs igseq/data/environment.yml:
diff -W60 -y \
	<(grep '^  run:$' conda/meta.yaml -A 1000 | tail -n +2 | grep ^$ -B 1000 | head -n -1 | sed 's/^ *//') \
	<(grep dependencies igseq/data/environment.yml -A 1000 | tail -n +2 | sed 's/^ *//')
[ $? -eq 0 ] && echo ...OK! || echo "...fix that!"
