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

if [ "a$nodes" == "a1" ]; then
  nodes2=""
else
  nodes2="-@ $nodes"
fi

#Create bam file for custom cells list, filter out chrM/chrUn/chrRandom, sort and index
awk '{printf("%s\n","'$cells_dir/'"$1)}' "$cells_list" > "00-cells.txt"
# Each contiguous chunk of the list is decoded and filtered to SAM text in its own pipeline.
# The chunks are concatenated in list order, so the header (first file's) and the record order match the serial pipeline.
ncells=$(wc -l < "00-cells.txt")
# nodes=0 (as in the serial script) runs one chunk; samtools still gets -@ 0
nchunks=$nodes
[ "$nchunks" -ge 1 ] || nchunks=1
chunk_size=$(( (ncells + nchunks - 1) / nchunks ))
[ "$chunk_size" -ge 1 ] || chunk_size=1
split -l "$chunk_size" "00-cells.txt" "00-chunk_"
# grep -v exits 1 when a chunk's reads are all removed; that is not an error
filter() { grep -v "$@" || [ $? -eq 1 ]; }
first="$(head -n 1 "00-cells.txt")"
chunks=(00-chunk_*)
pids=()
for chunk in "${chunks[@]}"; do
  (
    if [ "$chunk" == "${chunks[0]}" ]; then
      samtools view -h $nodes2 "$first" | filter '^@HD' | filter '^@PG' ; tail -n +2 "$chunk" | while read l; do samtools view $nodes2 "$l"; done
    else
      while read l; do samtools view $nodes2 "$l"; done < "$chunk"
    fi | awk '$3!="chrM"' | filter chrUn_ | filter GL00 | filter -e "random" > "$chunk.sam"
  ) &
  pids+=($!)
done
for pid in "${pids[@]}"; do wait "$pid"; done
cat 00-chunk_*.sam | samtools view -1 $nodes2 -o "02-filtered.bam" -
rm 00-chunk_*

#filter, sort and index bam file.
samtools sort -o "$output_bam" $nodes2 -l 1 02-filtered.bam
rm 02-filtered.bam
samtools index $nodes2 "$output_bam" "$output_bai"

#Step3A. Peak calling on aggregate population [Keep only significant peaks]
OMP_NUM_THREADS=$nodes MKL_NUM_THREADS=$nodes NUMEXPR_NUM_THREADS=$nodes OPENBLAS_NUM_THREADS=$nodes OMP_MAX_THREADS=$nodes MKL_MAX_THREADS=$nodes NUMEXPR_MAX_THREADS=$nodes OPENBLAS_MAX_THREADS=$nodes VECLIB_MAXIMUM_THREADS=$nodes macs2 callpeak -t "$output_bam" -n 04 -g $genome_size --nomodel --shift -75 --extsize 150 --keep-dup all --verbose 4 --call-summits -q $cutoff
mv 04_peaks.narrowPeak "$output_bed"
