#!/bin/bash
# Lingfei Wang, 2022. All rights reserved.

set -eo pipefail


function usage()
{
	echo "Usage: $(basename "$0") [-h] whole.bam output_folder [arguments ...]" >&2
	echo "Splits input whole.bam file by cell barcode and per-barcode bam files to output folder" >&2
	fmt='  %-18s%s\n'
	printf "$fmt" 'whole.bam' 'Input whole bam file containing reads with different barcodes' >&2
	printf "$fmt" 'output_folder' 'Output folder with one text file per barcode' >&2
	printf "$fmt" 'arguments' 'Arguments passed to split_bam_text.py' >&2
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
if [ "a$2" == "a" ] ; then usage; fi

fbam="$1"
dout="$(realpath "$2")"
#Temporary folder for text format reads
dtext="$(dirname "$dout")/$(basename "$dout")_text"
#Temporary file for bam header
fheader="$(dirname "$dout")/$(basename "$dout")_header"
shift 2

if [ -e "$dout" ] || [ -e "$dout"_text ] || [ -e "$dout"_header ]; then
	echo "Error: file/folder already exists: $dout, $dtext, or $fheader. Please remove these files/folders first."
fi

# Split bam file to per-cell text files
samtools view "$fbam" | dictys_helper split_bam_text.py "$@" "$dtext"
# Extract bam header
samtools view -H "$fbam" > "$fheader"
# Convert per-cell text to bam files with header
pushd . &> /dev/null
mkdir -p "$dout"
cd "$dtext"
fs="$(ls -1)"
echo "$fs" | while read f; do
	cat "$fheader" "$f" | samtools view -b -1 -o "$dout/$f".bam
done
popd &> /dev/null
# Clean up
rm -Rf "$dtext" "$fheader"
