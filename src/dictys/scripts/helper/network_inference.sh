#!/bin/bash
# Lingfei Wang, 2022. All rights reserved.

#Default values
clean=0
ncpu=4
ngpu=1

function usage()
{
	echo "Usage: $(basename "$0") [-h] [-j ncpu] [-J ngpu] [-c] mode" >&2
	echo "Runs network inference pipeline for input folder 'data' and output folder 'output'" >&2
	fmt='  %-12s%s\n'
	printf "$fmt" 'mode' 'Mode to run pipeline. Accepts:' >&2
	printf "$fmt" '' 'static: for context specific networks.' >&2
	printf "$fmt" '' 'dynamic: for dynamic networks.' >&2	 
	printf "$fmt" '-h' 'Display this help' >&2
	printf "$fmt" '-j ncpu' "Number of parallel jobs in CPU phase. Defaults to $ncpu." >&2
	printf "$fmt" '-J ngpu' "Number of parallel jobs in GPU phase. Ignored if no GPU is used. Defaults to $ngpu." >&2
	printf "$fmt" '-c' 'Clean up intermediate files after inference.' >&2
	exit 1
}

#Parse arguments
while getopts ':j:J:hc' o; do case "$o" in
	c)	clean="1";;
	j)	ncpu="$OPTARG";;
	J)	ngpu="$OPTARG";;
	:)	echo "Error: -${OPTARG} requires an argument." >&2;echo >&2;usage;;
	*)	usage;;
	esac
done
shift $((OPTIND-1))

if [ "a$1" == "a" ] || [ "a$2" != "a" ]; then
	usage
elif [ "a$1" != "astatic" ] && [ "a$1" != "adynamic" ]; then
	usage
fi

set -eo pipefail

mkfile="makefiles/$1.mk"
if [ "a$1" == "adynamic" ]; then
	make -f $mkfile subset
fi
make -f $mkfile -j $ncpu -k cpu 2> >( grep -v '^make: [*][*][*] No rule to make target' >&2 ) || true
#Run again in case of segfault from homer due to preparsed genome being overwritten in the first runs.
make -f $mkfile -j $ncpu -k cpu 2> >( grep -v '^make: [*][*][*] No rule to make target' >&2 ) || true
make -f $mkfile -j $ngpu -k gpu 2> >( grep -v '^make: [*][*][*] No rule to make target' >&2 ) || true
make -f $mkfile -j $ncpu -k cpu 2> >( grep -v '^make: [*][*][*] No rule to make target' >&2 ) || true
make -f $mkfile combine
if [ "a$clean" == "a1" ]; then
	make -f $mkfile clean
fi
























#
