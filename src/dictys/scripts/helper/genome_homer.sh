#!/bin/bash
# Lingfei Wang, 2022. All rights reserved.

function detect_basedir()
{
	b="$(which homer)"
	if [ "a$b" == "a" ]; then printf '%s' 'Failed: check if you have homer installed'; return; fi
	echo "$(dirname "$(dirname "$(realpath "$b")")")"
}


function usage()
{
	echo "Usage: $(basename "$0") [-b basedir] [-h] refgenome output_dir" >&2
	echo "Extracts reference genome from homer installation to output directory" >&2
	fmt='  %-14s%s\n'
	printf "$fmt" 'refgenome' 'Name of reference genome in homer format, e.g. hg38.' >&2
	printf "$fmt" '' 'You can get reference genomes available in homer with $basedir/configureHomer.pl -list' >&2
	printf "$fmt" 'output_dir' 'Output directory to export reference genome as' >&2
	printf "$fmt" '-b basedir' 'Base directory of homer installation' >&2
	printf "$fmt" '' 'Default: autodetect ('"$(detect_basedir)"')' >&2
	printf "$fmt" '-h' 'Display this help' >&2
	exit 1
}

#Parse arguments
while getopts ':b:h' o; do case "$o" in
	b)	basedir="$(realpath "$OPTARG")";;
	:)	echo "Error: -${OPTARG} requires an argument." >&2;echo >&2;usage;;
	*)	usage;;
	esac
done
shift $((OPTIND-1))
if [ "a$2" == "a" ] || [ "a$3" != "a" ]; then usage; fi
refgenome="$1"
output_dir="$2"

set -eo pipefail
#Remove possible trailing /
output_dir="$(realpath "$output_dir")"
output_dir="$(dirname "$output_dir")/$(basename "$output_dir")"

if [ "a$basedir" == "a" ]; then basedir="$(detect_basedir)"; fi
if ! [ -e "$basedir" ]; then echo "Error: cannot find homer base directory $basedir." >&2; exit 1; fi

cd "$basedir"
chmod u+x configureHomer.pl
if [ "a$(./configureHomer.pl -list 2> /dev/null | awk '$1=="-" && $2=="'"$refgenome"'"' | wc -l)" == "a1" ]; then
	#Download unavailable refgenome
	echo "Downloading reference genome $refgenome in homer"
	err="$(./configureHomer.pl -install $refgenome 2>&1 || true)"
	chmod u+x configureHomer.pl
fi
if [ "a$(./configureHomer.pl -list 2> /dev/null | awk '$1=="+" && $2=="'"$refgenome"'"' | wc -l)" != "a1" ]; then
	echo "Error: reference genome $refgenome not found in homer"
	if [ "a$err" != "a" ]; then
		echo "Possible error during reference genome installation:"
		echo "$err"
	fi
	exit 1
fi

#Copy reference genome to output_dir
cp -R "data/genomes/$refgenome" "$output_dir"
#Remove preparsed data
rm -Rf "$output_dir"/preparsed






















#
