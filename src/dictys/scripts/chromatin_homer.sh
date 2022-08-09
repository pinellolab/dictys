#!/bin/bash
# Nikolaos Trasanidis, Lingfei Wang, 2022. All rights reserved.

# Motif scan with homer
# Parameters:
# $1 : Path of input bed file of regions
# $2 : Path of input motif list file in homer format
# $3 : Path of reference genome folder for homer
# $4 : Python script path
# $5 : number of threads

set -e -o pipefail

#arguements
bedfile="$1"
motifdb="$2"
genome="$3"
fi_r="$4"
nodes="$5"

mkdir 15-motifscan 14-reform-split 15-tfs
#Split bed file
nline="$(( ( $( cat "$bedfile" | wc -l ) / nodes ) + 1 ))"
split -l $nline -a 5 "$bedfile" 14-reform-split/
#Run homer in parallel
ls -1 14-reform-split | while read l; do
{
	mkdir 15-motifscan/$l
	#Remove error message: ERROR: outstripped buffer
	findMotifsGenome.pl 14-reform-split/$l $genome 15-motifscan/$l -size given -mask -find "$motifdb" 2>&1 > 15-tfs/$l | grep -v 'ERROR: outstripped buffer' >&2 || true
	touch 15-tfs/$l.done
} &
done
#Get first file with header
f1="$(ls -1 14-reform-split | head -n 1)"
while ! [ -e 15-tfs/"$f1".done ]; do
	sleep 1
done
mv 15-tfs/"$f1" 15-tf.bed
#Get other files sequentially
ls -1 14-reform-split | tail -n +2 | while read l; do
	while ! [ -e 15-tfs/"$l".done ]; do
		sleep 1
	done
	tail -n +2 15-tfs/"$l" >> 15-tf.bed
	rm -f 15-tfs/"$l"
done
#Cleanup
rm -Rf 14-reform-split 15-motifscan 15-tfs

#Step3. Reformat and return files
python3 "$fi_r"























#
