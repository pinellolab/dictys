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
	echo "Usage: $(basename "$0") [-f field] [-h] gtf_file bed_file" >&2
	echo "Extracts gene region from GTF file into bed file" >&2
	fmt='  %-14s%s\n'
	printf "$fmt" 'gtf_file' 'Path of input GTF file' >&2
	printf "$fmt" 'bed_file' 'Path of output BED file' >&2
	printf "$fmt" '-f field' 'Field name to extract. Default: gene_name' >&2
	printf "$fmt" '-h' 'Display this help' >&2
	exit 1
}

#Parse arguments
field="gene_name"
while getopts ':f:h' o; do case "$o" in
	f)	field="$OPTARG";;
	:)	echo "Error: -${OPTARG} requires an argument." >&2;echo >&2;usage;;
	*)	usage;;
	esac
done
shift $((OPTIND-1))
if [ "a$2" == "a" ] || [ "a$3" != "a" ]; then usage; fi
fgtf="$1"
fbed="$2"

set -eo pipefail

#Gene entries only with given field
grep -v '^#' "$fgtf" | grep "$field" | awk -F "\t" '$3=="gene" && $9!="" {print $1"\t"$4"\t"$5"\t.\t"$7"\t"$9}' > "$fbed".step1
#Columns without designated field
awk -F "\t" '{print $6}' "$fbed".step1 | awk -F "; " '{v=""; for(i=1;i<=NF;i++) if($i ~ /^'"$field"' ".*"$/)v=v $i " ";print v}' | tr -d '"' | awk '{for(i=0;i<=NF;i++) if($i=="'"$field"'")$i="";print}' | sed 's/^ *//g' > "$fbed".step2
#Designated field
awk -F "\t" '{if(!($1 ~ /^chr/))$1="chr"$1;print $1"\t"$2"\t"$3"\t"$4"\t"$5}' "$fbed".step1 > "$fbed".step3
#Merge to single file
paste "$fbed".step3 "$fbed".step2 | awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5}'> "$fbed".step4
#Remove entries with duplicates in designated field
sort "$fbed".step2 | uniq -c | awk '$1!="1"{print "\t"$2"\t"}' > "$fbed".step5
grep -vf "$fbed".step5 -F "$fbed".step4 > "$fbed".step6
#Then add back the first entry for each duplicate value in designated field
grep -f "$fbed".step5 -F "$fbed".step4 > "$fbed".step7
cat "$fbed".step5 | while read l; do grep -F -m 1 "$l" "$fbed".step7 >> "$fbed".step6; done
#Remove duplicate genes
t1="$(awk '{print $4}' "$fbed".step6 | uniq -c | grep -v '^[[:blank:]]*1[[:blank:]]' | awk '{print $2}' || true )"
if [ "a$(printf "$t1" | wc -c )" != "a0" ]; then
	echo "Warning. Duplicate TSS found. Only one is used for genes: $(printf "$t1" | tr "\n" ",")" > /dev/stderr
	sort -k4 -u "$fbed".step6 > "$fbed".step8
else
	cp "$fbed".step6 "$fbed".step8
fi
#Sort
grep '^chr[1-9]' "$fbed".step8 | sort -k2g,3g | sort -s -k1.4g > "$fbed"
grep -v '^chr[1-9]' "$fbed".step8 | sort -k2g,3g | sort -s -k1 >> "$fbed"
rm -f "$fbed".step[1-8]
























#
