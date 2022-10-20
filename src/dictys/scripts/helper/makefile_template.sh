#!/bin/bash
# Lingfei Wang, 2022. All rights reserved.

function usage()
{
	echo "Usage: $(basename "$0") [-h] [makefile1.mk ...]" >&2
	echo "Generate network inference pipeline makefiles in current working folder from template" >&2
	fmt='  %-20s%s\n'
	printf "$fmt" 'makefile1.mk ...' 'Name of each makefile to generate from template.' >&2
	printf "$fmt" '' 'If omitted, all available makefiles will be generated.' >&2
	printf "$fmt" '-h' 'Display this help' >&2
	exit 1
}

#Parse arguments
while getopts ':h' o; do case "$o" in
	:)	echo "Error: -${OPTARG} requires an argument." >&2;echo >&2;usage;;
	*)	usage;;
	esac
done
shift $((OPTIND-1))

set -eo pipefail
dmake="$(dirname "$(realpath "$0")")/../makefiles"

if [ "a$1" != "a" ]; then
	#Generate specific files
	while [ "a$1" != "a" ]; do
		cp "$dmake/$1" ./
		shift
	done
else
	#Generate all available files
	ls -1 "$dmake" | grep '[.]mk$' | while read f; do
		cp "$dmake/$f" ./
	done
fi
























#
