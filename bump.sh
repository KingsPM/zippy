#!/bin/sh

# this script bumps version and sets status for all python scripts in folder
ARGS=("$@")

if [ $# == 2 ]; then
	echo "Substituting ${ARGS[0]} with ${ARGS[1]} in all python scripts"
	for f in `find . -name '*.py'`; do
		echo 's/^(__'${ARGS[0]}'__\s*=\s*)[^\n]+/${1}'${ARGS[1]}'/x' $f
		perl -i -pe 's/^(__'${ARGS[0]}'__\s*=\s*)[^\n]+/${1}'${ARGS[1]}'/x' $f
	done
else
	echo "Substituting give field and new value. Example: bump.sh version 1.2"
fi

