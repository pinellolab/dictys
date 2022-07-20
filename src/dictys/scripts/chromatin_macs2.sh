#!/bin/bash
# Nikolaos Trasanidis, Lingfei Wang, 2022. All rights reserved.

# Peak calling with macs2
# Parameters:
# $1 : Path of input file containing one sample/cell name per line for macs2 peak calling
# $2 : Path of input folder that contains each cell's bam file by name in *f_names*
# $3 : Path of output bam file for select samples/cells
# $4 : Path of output bai file for select samples/cells
# $5 : Path of output bed file of peaks
# $6 : Genome size input of macs2. Use shortcuts hs or mm for human or mouse. See macs2 parameter -g for other custom numbers.
# $7 : Qvalue cutoff for macs2
# $8 : Number of threads

set -eo pipefail

#arguements
cells_list="$1"
cells_dir="$2"
output_bam="$3"
output_bai="$4"
output_bed="$5"
genome_size=$6
cutoff=$7
nodes=$8

#Create bam file for custom cells list, filter out chrM/chrUn/chrRandom, sort and index
awk '{printf("%s\n","'$cells_dir/'"$1)}' "$cells_list" > "00-cells.txt"
( samtools view -h -@ "$nodes" "$(head -n 1 "00-cells.txt" )" | grep -v '^@HD' | grep -v '^@PG' ; tail -n +2 "00-cells.txt" | while read l; do samtools view -@ "$nodes" "$l"; done ) | awk '$3!="chrM"' |  grep -v chrUn_ | grep -v GL00 | grep -v -e "random" | samtools view -1 -@ "$nodes" -o "02-filtered.bam" -

#filter, sort and index bam file.
samtools sort -o "$output_bam" -@ "$nodes" -l 1 02-filtered.bam
rm 02-filtered.bam
samtools index -@ "$nodes" "$output_bam" "$output_bai"

#Step3A. Peak calling on aggregate population [Keep only significant peaks]
OMP_NUM_THREADS=$nodes MKL_NUM_THREADS=$nodes NUMEXPR_NUM_THREADS=$nodes OPENBLAS_NUM_THREADS=$nodes OMP_MAX_THREADS=$nodes MKL_MAX_THREADS=$nodes NUMEXPR_MAX_THREADS=$nodes OPENBLAS_MAX_THREADS=$nodes VECLIB_MAXIMUM_THREADS=$nodes macs2 callpeak -t "$output_bam" -n 04 -g $genome_size --nomodel --shift -75 --extsize 150 --keep-dup all --verbose 4 --call-summits -q $cutoff
mv 04_peaks.narrowPeak "$output_bed"
