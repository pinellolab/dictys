#!/bin/bash

set -eo pipefail


function usage()
{
	echo "Usage: $(basename "$0") [-h] repo file [arguments ...]" >&2
	echo "Downloads and decompresses split tar.xz file from repository to current folder." >&2
	echo "Overwrites destination and intermediate files." >&2
	fmt='  %-12s%s\n'
	printf "$fmt" 'repo' 'Repository name to download from. Determines which helper script to use. Accepts: dataverse.' >&2
	printf "$fmt" 'file' 'File name to download and decompress' >&2
	printf "$fmt" 'arguments' 'Arguments passed to repository based script' >&2
	printf "$fmt" '-h' 'Display this help' >&2
	exit 1
}

set -eo pipefail

#Parse arguments
while getopts ':h' o; do case "$o" in
	:)	echo "Error: -${OPTARG} requires an argument." >&2;echo >&2;usage;;
	*)	usage;;
	esac
done
shift $((OPTIND-1))
if [ "a$2" == "a" ] ; then usage; fi

repo="$1"
fname="$2"
shift 2

#repo dependent commands
if [ "a$repo" == "adataverse" ]; then
	cmd_ls='dictys_helper repo_dataverse.py ls '"$@"
	cmd_down1='dictys_helper repo_dataverse.py down '"$@"
	cmd_down2=''
else
	usage
fi

#Download files
rm -f "$fname"
fs="$($cmd_ls | grep '^'"$fname"'[.][0-9]*$' | sort )"
for f in $fs; do
	$cmd_down1 $f $cmd_down2
	if ! [ -e "$f" ]; then 
		echo "Cannot find file $f"
		exit 1
	fi
	cat "$f" >> "$fname"
	rm "$f"
done
#Decompress
tar xf "$fname"
rm -f "$fname"

