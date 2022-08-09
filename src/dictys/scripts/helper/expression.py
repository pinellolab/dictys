#!/usr/bin/env python3
from os.path import join as pjoin
import logging
import argparse
import pandas as pd
from scipy.io import mmread

parser = argparse.ArgumentParser(description="Converts mtx.gz format expression file to tsv.gz format.")
parser.add_argument('input_folder',type=str,help='Input folder that contains matrix.mtx.gz, features.tsv.gz, and barcodes.tsv.gz.')
parser.add_argument('output_file',type=str,help='Output file in tsv.gz format')
parser.add_argument('--column',default=1,type=str,help='Column ID in features.tsv.gz for gene name. Starts with 0. Defaults: 1.')

args=parser.parse_args()
diri=args.input_folder
fo=args.output_file
colid=args.column

#Read files
d=mmread(pjoin(diri,'matrix.mtx.gz'))
d=d.toarray()
names=[pd.read_csv(pjoin(diri,x+'.tsv.gz'),header=None,index_col=None,sep='\t') for x in ['features','barcodes']]
assert names[1].shape[1]==1
names[0]=names[0][colid].values
names[1]=names[1][0].values
assert d.shape==tuple(len(x) for x in names)

#Select unique rows & columns
inames=['gene','cell']
for xi in range(2):
	td=set()
	ids=[]
	for xj in range(len(names[xi])):
		if names[xi][xj] not in td:
			td.add(names[xi][xj])
			ids.append(xj)
		else:
			logging.warning('Skipped duplicate occurence of {} name: {}'.format(inames[xi],names[xi][xj]))
	d=d.swapaxes(0,xi)[ids].swapaxes(0,xi)
	names[xi]=names[xi][ids]
	assert len(names[xi])==len(set(names[xi]))
d=pd.DataFrame(d,index=names[0],columns=names[1])

#Remove unneeded genes
t1=[':' not in x and '.' not in x for x in d.index]
#Output
d=d.loc[t1]
d=d.loc[sorted(d.index)]
d=d[sorted(d.columns)]
d.to_csv(fo,header=True,index=True,sep='\t')
