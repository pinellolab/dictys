#!/usr/bin/env python3
# Lingfei Wang, 2022. All rights reserved.
import sys
from os import linesep,mkdir
from os.path import join as pjoin
from collections import defaultdict
from shutil import rmtree
import logging
import argparse
import pandas as pd


def cleanup(diro,funknown):
	import os
	rmtree(diro,ignore_errors=True)
	if funknown is not None:
		try:
			os.remove(funknown)
		except FileNotFoundError:
			pass


parser = argparse.ArgumentParser(
	description="Splits input bam file (stdin from samtools view) by cell barcode and outputs headerless individual text file per barcode to output folder.",
	usage='samtools view whole.bam | python3 split_bam_text.py [-h] [--output_unknown OUTPUT_UNKNOWN] [--section SECTION] [--buffer_size BUFFER_SIZE] [--ref_expression REF_EXPRESSION] [--namemap NAMEMAP] output_folder')
parser.add_argument('output_folder',type=str,help='Output folder with one text file per barcode')
parser.add_argument('--output_unknown',type=str,default=None,help='Output text file for reads without barcodes or with unknown barcodes (see --ref_expression)')
parser.add_argument('--section',type=str,default='RG:',help='Section header that contains cell barcode. Must be the same list of cell barcodes/names as other places in the pipeline, e.g. `subsets/*/names_atac.txt` and `coord_atac.tsv.gz`. Default: "RG:".')
parser.add_argument('--buffer_size',type=int,default=10000,help='Buffer read counts for output of each barcode. Larger value consumes more memory for faster speed. Reduce if you see a MemoryError. Default: 10000.')
parser.add_argument('--ref_expression',type=str,default=None,help='Cell RNA barcode reference file as expression.tsv.gz. If specified, cell barcodes not contained in the reference file are also regarded as unknown.')
parser.add_argument('--namemap',type=str,default=None,help='Cell barcode map from RNA read barcodes to ATAC read barcodes in format file_path,RNA_column_ID,ATAC_column_ID. File should be in tsv format. If unset, will regard RNA and ATAC barcdoes identical (identity map).')

args=parser.parse_args()
diro=args.output_folder
funknown=args.output_unknown
section_cell=args.section
buffline=args.buffer_size
ref=args.ref_expression
namemap=args.namemap

if ref is not None:
	ref=pd.read_csv(ref,header=0,index_col=0,sep='\t',nrows=1)
	ref=set(ref.columns)
if namemap is not None:
	namemap=namemap.split(',')
	assert len(namemap)==3
	t1=pd.read_csv(namemap[0],index_col=None,header=None,sep='\t')
	t1=t1[[int(x) for x in namemap[1:]]].values[:,::-1]
	namemap=dict(t1)
	if ref is not None:
		ref=set(x[0] for x in namemap.items() if x[1] in ref)

ans=defaultdict(list)
err=False
cleanup(diro,funknown)
try:
	mkdir(diro)
except FileNotFoundError:
	raise
except:
	pass
try:
	while True:
		ls=sys.stdin.readlines(buffline)
		if len(ls)==0:
			break
		#Go through individual lines
		t1=[list(filter(lambda x:x.startswith(section_cell),x.split('\t'))) for x in ls]
		if any(len(x)>1 for x in t1):
			raise RuntimeError('Input line has multiple fields satisfying grouping criteria: '+section_cell)
		t1=[x[0][len(section_cell):] if len(x)>0 else None for x in t1]
		if ref is not None:
			t1=[x if x in ref else None for x in t1]
		for xi in range(len(t1)):
			ans[t1[xi]].append(ls[xi].strip())
		if not err and len(ans)>100000:
			logging.warning('This script is designed for <100000 cells/barcodes. Raise an issue to request a more advance script if your dataset is larger.')
			err=True
		for xi in list(filter(lambda x:len(ans[x])>=buffline,set(t1))):
			#Append to output file with buffer
			if xi is None:
				fo=funknown
			elif namemap is not None:
				if xi not in namemap:
					fo=funknown
				else:
					fo=pjoin(diro,namemap[xi])
			else:
				fo=pjoin(diro,xi)
			if fo is not None:
				with open(fo,'a') as f:
					f.write(linesep.join(ans[xi])+linesep)
			ans[xi]=[]
except:
	#Cleanup output folder
	cleanup(diro,funknown)
	raise

#Output remaining lines
for xi in ans:
	if len(ans[xi])==0:
		continue
	if xi is None:
		fo=funknown
	elif namemap is not None:
		if xi not in namemap:
			fo=funknown
		else:
			fo=pjoin(diro,namemap[xi])
	else:
		fo=pjoin(diro,xi)
	if fo is not None:
		with open(fo,'a') as f:
			f.write(linesep.join(ans[xi])+linesep)
