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
	echo "Usage: $(basename "$0") [-b basedir] [ (-m mapfile) | (-o organism) ] [-c capitalization] [-h]" >&2
	echo "Extracts motif file from homer installation to stdout" >&2
	fmt='  %-20s%s\n'
	printf "$fmt" '-b basedir' 'Base directory of homer installation' >&2
	printf "$fmt" '' 'Default: autodetect ('"$(detect_basedir)"')' >&2
	printf "$fmt" '-m mapfile' 'Mapfile mode: use motifs in $basedir/motifs/ by specifying a motif to gene mapping file.' >&2
	printf "$fmt" '' 'The mapping file is a two-column headered text file mapping motif file (column 0) to gene name (column 1).' >&2
	printf "$fmt" '' 'Default: $basedir/motifs/extras/motifs2symbol.txt' >&2
	printf "$fmt" '-o organism' 'Organism mode: use motifs in $basedir/data/knownTFs/organism/known.motifs by specifying an organism.' >&2
	printf "$fmt" '' 'Each motifs is directly mapped to the gene in the front of its name (separated by :).' >&2
	printf "$fmt" '' 'Those without gene names are kept but will be disgarded in the inference pipeline.' >&2
	printf "$fmt" '' 'If option is unspecified, uses -m with its default setting.' >&2
	printf "$fmt" '-c capitalization' 'Capitalization conversion for gene name. Accepts:' >&2
	printf "$fmt" '' '0: no conversion (default)' >&2
	printf "$fmt" '' '1: uncapitalized' >&2
	printf "$fmt" '' '2: CAPITALIZED' >&2
	printf "$fmt" '' '3: First Character Capitalized' >&2
	printf "$fmt" '-h' 'Display this help' >&2
	exit 1
}

set -eo pipefail

#Parse arguments
capitalization='0'
while getopts ':b:m:o:c:h' o; do case "$o" in
	b)	basedir="$(realpath "$OPTARG")";;
	m)	mapfile="$(realpath "$OPTARG")";;
	o)	organism="$OPTARG";;
	c)	capitalization="$OPTARG"
		case "$capitalization" in
			[0123])	;;
			*)	echo "Error: -c has an unrecognized argument value." >&2;echo >&2;usage;;
		esac;;
	:)	echo "Error: -${OPTARG} requires an argument." >&2;echo >&2;usage;;
	*)	usage;;
	esac
done
shift $((OPTIND-1)) || true
if [ "a$1" != "a" ]; then usage; fi


if [ "a$basedir" == "a" ]; then basedir="$(detect_basedir)"; fi
if ! [ -e "$basedir" ]; then echo "Error: cannot find homer base directory $basedir." >&2; exit 1; fi
if [ "a$organism" != "a" ]; then
	#Organism mode
	awkscript1=$(cat <<-END
	BEGIN {RS=">";ORS=">"}
	{
		t1=\$2
		gsub(/[\(\)\/]/,"_",t1)
		ind=match(t1,/_/)-1
		suf=substr(t1,ind+1)
		pref=substr(\$2,1,ind)
		split(pref,prefs,":")
		for (pref1 in prefs)
		{
			pref21=substr(prefs[pref1],1,1)
			pref22=substr(prefs[pref1],2)
	END
	)
	awkscript2=$(cat <<-END
			\$2=pref31 pref32 suf
			print
		}
	}
	END
	)
	fname="$basedir/data/knownTFs/$organism/known.motifs"
	if ! [ -e "$fname" ]; then echo "Error: cannot find motif file $fname." >&2; exit 1; fi
	case "$capitalization" in
		0) awkscriptc='pref31=pref21; pref32=pref22';;
		1) awkscriptc='pref31=tolower(pref21); pref32=tolower(pref22)';;
		2) awkscriptc='pref31=toupper(pref21); pref32=toupper(pref22)';;
		3) awkscriptc='pref31=toupper(pref21); pref32=tolower(pref22)';;
	esac
	printf ">"
	awk -F "\t" "$awkscript1"";$awkscriptc;""$awkscript2" "$fname" | head -n -2
else
	#Mapfile mode
	if [ "a$mapfile" == "a" ]; then mapfile="$basedir/motifs/extras/motifs2symbol.txt"; fi
	if ! [ -e "$mapfile" ]; then echo "Error: cannot find gene name mapping file $mapfile." >&2; exit 1; fi

	cd "$basedir"/motifs

	tail -n +2 "$mapfile" | while read f sym; do
		if ! [ -f "$f" ]; then echo "Warning: skipped motif file because not found: $f" >&2; continue; fi
		if [ "a$(head -n 1 "$f" | grep "^>" | wc -l)" != "a1" ]; then echo "Warning: skipped motif file because not starting from first line: $f" >&2; continue; fi
		if [ "a$(grep "^>" "$f" | wc -l)" != "a1" ]; then echo "Warning: skipped motif file because not containing exactly one motif: $f" >&2; continue; fi
		#Process gene name
		case "$capitalization" in
			0) symname="$sym";;
			1) symname="$(echo "$sym" | tr '[:upper:]' '[:lower:]')";;
			2) symname="$(echo "$sym" | tr '[:lower:]' '[:upper:]')";;
			3) symname="$(echo "$sym" | head -c 1 | tr '[:lower:]' '[:upper:]'; echo "$sym" | tail -c +2 | tr '[:upper:]' '[:lower:]')";;
		esac
		#Process file path
		symname="$symname""_$(echo "$f" | tr '/\\.' '_')"
		awk 'NR==1{$2="'"$symname"'"} {print}' "$f"
	done
fi
























#
