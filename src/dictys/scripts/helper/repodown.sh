#!/bin/bash
# Lingfei Wang, 2022. All rights reserved.

set -eo pipefail


function usage()
{
	echo "Usage: $(basename "$0") [-h] [-n] repo file [arguments ...]" >&2
	echo "Downloads and decompresses split tar.xz file from repository to current folder." >&2
	echo "Overwrites destination and intermediate files." >&2
	fmt='  %-12s%s\n'
	printf "$fmt" 'repo' 'Repository name to download from. Determines which helper script to use. Accepts: zenodo, dataverse.' >&2
	printf "$fmt" 'file' 'File name to download and decompress' >&2
	printf "$fmt" 'arguments' 'Arguments passed to repository based script' >&2
	printf "$fmt" '-n' 'No decompression with tar xf' >&2
	printf "$fmt" '-h' 'Display this help' >&2
	exit 1
}

decompress='1'
#Parse arguments
while getopts ':hn' o; do case "$o" in
	n)	decompress='0';;
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
elif [ "a$repo" == "azenodo" ]; then
	cmd_ls='dictys_helper repo_zenodo.py ls '"$@"
	cmd_down1='dictys_helper repo_zenodo.py down '"$@"
	cmd_down2=''
else
	usage
fi

#Download files
rm -f "$fname"
f1="$($cmd_ls | grep '^'"$fname"'$' | sort || true )"
if [ "a$f1" != "a" ]; then
	#Download single file
	echo "Downloading $fname"
	$cmd_down1 $fname $cmd_down2
	if ! [ -e "$fname" ]; then 
		echo "Downloading failed for $fname. Now exit."
		exit 1
	fi
else
	#Download split files
	fs="$($cmd_ls | grep '^'"$fname"'[.][0-9]*$' | sort || true )"
	if [ "a$fs" == "a" ]; then
		echo "File not found: $f. Now exit."
		exit 1
	fi
	for f in $fs; do
		echo "Downloading $f"
		$cmd_down1 $f $cmd_down2
		if ! [ -e "$f" ]; then 
			echo "Downloading failed for $f. Now exit."
			exit 1
		fi
		cat "$f" >> "$fname"
		rm "$f"
	done
fi

if [ "a$decompress" == "a1" ]; then
	#Decompress
	tar xf "$fname"
	rm -f "$fname"
fi
