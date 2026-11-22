#!/bin/bash
# Nikolaos Trasanidis, Lingfei Wang, 2022. All rights reserved.

# TF Footprinting with wellington
# Parameters:
# $1 : Path of input bam file of all reads
# $2 : Path of input bai file of all reads
# $3 : Path of input bed file of peaks
# $4 : Path of output bed file of footprints
# $5 : Cutoff for wellington score
# $6 : Number of threads
# $7 : Maximum number of footprints to retain, ordered by wellington score. Use 0 for no limit.
# $8 : Path of input bed file of blacklisted genome regions to be removed. Use None to disable.

set -eo pipefail

bamfile="$1"
baifile="$2"
bedfile="$3"
bedfile_out="$4"
Wellington_cutoff=$5
nodes="$6"
threshold="$7"
blacklist="$8"

if [ "a$blacklist" != "a" ] && [ "a$blacklist" != "aNone" ]; then
	#Remove blacklisted regions
	bedtools subtract -a "$bedfile" -b "$blacklist" -A > 11-filtered.bed
	bedfile="11-filtered.bed"
fi
#Refine bed file
bedtools merge -i "$bedfile" > 11-merged.bed
bedtools sort -i 11-merged.bed > 12-sorted.bed
rm -f 11-merged.bed

#Wellington
mkdir -p 13-wellington
# pyDNase's own wellington_footprints.py also runs the FDR null and writes files this script never reads.
# Use the dictys driver, which skips them, in place of it; any other program of that name is called as is.
wellington_cmd=(wellington_footprints.py)
wellington_path="$(command -v wellington_footprints.py || true)"
if [ -n "$wellington_path" ] && grep -q 'import pyDNase' "$wellington_path"; then
	read -r -a wellington_shebang <<< "$(head -n 1 "$wellington_path" | sed -e 's/^#![[:space:]]*//')"
	wellington_python="${wellington_shebang[0]:-}"
	if [ "${wellington_python##*/}" = env ]; then
		# '#!/usr/bin/env python3': the interpreter is the first word of env's command, not env itself
		wellington_python=
		for wellington_word in "${wellington_shebang[@]:1}"; do
			case "$wellington_word" in -*|*=*) ;; *) wellington_python="$(command -v "$wellington_word" || true)"; break;; esac
		done
	fi
	[ -x "$wellington_python" ] || wellington_python=python3
	wellington_cmd=("$wellington_python" "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/wellington_footprints.py")
fi
"${wellington_cmd[@]}" 12-sorted.bed "$bamfile" 13-wellington -fdr 0.999 -fdrlimit -5 -pv -$Wellington_cutoff -A -p $(( nodes>1 ? $nodes - 1 : 1 ))
rm -f 12-sorted.bed

#Reformat output bed file
( awk '$5< -'$Wellington_cutoff'' 13-wellington/'p value cutoffs'/*-$Wellington_cutoff.bed | awk '{printf("%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4":"$1":"$2":"$3":"$5,$5)}' > 14-reform-full.bed ) || true

if [ "a$threshold" != "a0" ]; then
	# Limit footprint count
	( sort -k5,5n -T . 14-reform-full.bed | head -n $threshold > "$bedfile_out" ) || true
else
	mv 14-reform-full.bed "$bedfile_out"
fi






















#
