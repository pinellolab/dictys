#!/usr/bin/env python3
import json
import argparse
import sys
from os.path import join as pjoin
from os.path import exists as pexists
from os.path import isdir
from os import listdir
import logging
import itertools
import numpy as np
import pandas as pd
from dictys.traj import trajectory,point

parser = argparse.ArgumentParser(description="Validates makefile and input data of network inference pipeline.",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--dir_data',type=str,default='data',help='Directory of input data folder.')
parser.add_argument('--dir_makefiles',type=str,default='makefiles',help='Directory of makefiles folder.')

args=parser.parse_args()
dirdata=args.dir_data
dirmakefiles=args.dir_makefiles

if not isdir(dirdata) or not isdir(dirmakefiles):
	raise FileNotFoundError('Cannot find input data or makefiles folder. See python {} -h for help.'.format(sys.argv[0]))

#Detect whether joint profiles from makefile
with open(pjoin(dirmakefiles,'config.mk'),'r') as f:
	s=f.readlines()
s=[x.strip() for x in s if x.startswith('JOINT=')][-1]
s=s[len('JOINT='):]
if s not in {'0','1'}:
	raise ValueError('Invalid JOINT variable in '+pjoin(dirmakefiles,'config.mk'))
isjoint=s=='1'
print(f'Joint profile: {isjoint}')

#Gene & cell names from RNA
namec_rna=pd.read_csv(pjoin(dirdata,'expression.tsv.gz'),header=0,index_col=0,sep='\t',nrows=1)
namec_rna=np.array(namec_rna.columns)
snamec_rna=set(namec_rna)
print(f'Found {len(namec_rna)} cells with RNA profile')
if len(namec_rna)<100:
	raise ValueError('<100 cells found with RNA profile')
nameg=pd.read_csv(pjoin(dirdata,'expression.tsv.gz'),header=0,index_col=0,sep='\t',usecols=[0,1])
nameg=np.array(nameg.index)
snameg=set(nameg)
print(f'Found {len(nameg)} genes with RNA profile')
if len(nameg)<100:
	raise ValueError('<100 genes found with RNA profile')

#Cell names from ATAC
namec_atac=listdir(pjoin(dirdata,'bams'))
namec_atac=np.array(['.'.join(x.split('.')[:-1]) for x in namec_atac if x.endswith('.bam')])
snamec_atac=set(namec_atac)
print(f'Found {len(namec_atac)} cells with ATAC profile')
if isjoint and len(snamec_rna-snamec_atac)>0:
	raise FileNotFoundError('Not all cells with RNA profile has a bam file in bams folder for the joint profiling dataset.')

#TF names from motifs
with open(pjoin(dirdata,'motifs.motif'),'r') as f:
	nameg_motif=f.readlines()
nameg_motif=[x.split('\t')[1] for x in nameg_motif if x.startswith('>')]
if len(nameg_motif)!=len(set(nameg_motif)):
	raise ValueError('Duplicate motifs found')
print(f'Found {len(nameg_motif)} motifs')
nameg_motif=set([x.split('_')[0] for x in nameg_motif])
print(f'Found {len(nameg_motif)} TFs')
nameg_motif_notfound=sorted(nameg_motif-set(nameg))
nameg_motif=sorted(nameg_motif&set(nameg))
snameg_motif=set(nameg_motif)
print(f'Found {len(nameg_motif)} TFs in current dataset')
print(f'Missing {len(nameg_motif_notfound)} TFs in current dataset: '+','.join(nameg_motif_notfound))
if len(nameg_motif)<10:
	raise ValueError('<10 TFs found')

#Chromosomes and gene names from bed file
with open(pjoin(dirdata,'gene.bed'),'r') as f:
	nameg_bed=f.readlines()
nameg_bed=[x.strip() for x in nameg_bed]
nameg_bed=[x.split('\t') for x in nameg_bed if len(x)>0]
if len(set([len(x) for x in nameg_bed]))!=1:
	raise ValueError('Unequal number of columns in gene.bed')
nameg_bed=list(filter(lambda x:x[3] in snameg,nameg_bed))
if len(nameg_bed[0])<6:
	logging.warn('No strand information in gene.bed')
namechr_bed=sorted(set([x[0] for x in nameg_bed]))
snamechr_bed=set(namechr_bed)
nameg_bed=[x[3] for x in nameg_bed]
snameg_bed=set(nameg_bed)
if len(nameg_bed)!=len(snameg_bed):
	from collections import Counter
	raise ValueError('Duplicate gene found in gene.bed')
nameg_bed=sorted(nameg_bed)
print(f'Found {len(nameg_bed)} genes with TSS information')
if len(nameg_bed)<100:
	raise ValueError('<100 genes with TSS information')

#Ref genome chromosomes
namechr=listdir(pjoin(dirdata,'genome'))
namechr=list(filter(lambda x:x.endswith('.fa'),namechr))
t1=[]
for xi in namechr:
	with open(pjoin(dirdata,'genome',xi),'r') as f:
		t1.append(list(filter(lambda x:x.startswith('>'),f.readlines())))
namechr=[x.strip().strip('>') for x in itertools.chain.from_iterable(t1)]
snamechr=set(namechr)
if len(snamechr_bed-snamechr)>0:
	logging.warning('Chromosomes not found in reference genome: '+','.join(sorted(snamechr_bed-snamechr)))

#Blacklist chromosomes
if pexists(pjoin(dirdata,'blacklist.bed')):
	with open(pjoin(dirdata,'blacklist.bed'),'r') as f:
		namechr_bl=f.readlines()
	namechr_bl=[x.strip() for x in namechr_bl]
	namechr_bl=[x.split('\t')[0] for x in namechr_bl if len(x)>0]
	snamechr_bl=set(namechr_bl)
	if len(snamechr_bl-namechr_bl)>0:
		logging.warn('Blacklist chromosomes not found in genome: '+','.join(sorted(snamechr_bl-namechr_bl)))

#############################################
# Context specific GRN inference checks
#############################################

if not pexists(pjoin(dirmakefiles,'static.mk')) or not pexists(pjoin(dirdata,'subsets.txt')):
	logging.warn('Cannot find static.mk or subsets.txt. Skipping static network inference checks.')
else:
	#Subsets
	with open(pjoin(dirdata,'subsets.txt'),'r') as f:
		names=f.readlines()
	names=[x.strip() for x in names]
	names=list(filter(lambda x:len(x)>0,names))

	#Cell names in each subset
	#Loading cell names
	namec_srna=[]
	namec_satac=[]
	for xi in names:
		with open(pjoin(dirdata,'subsets',xi,'names_rna.txt'),'r') as f:
			t1=f.readlines()
		t1=[x.strip() for x in t1]
		t1=list(filter(lambda x:len(x)>0,t1))
		namec_srna.append(t1)
		with open(pjoin(dirdata,'subsets',xi,'names_atac.txt'),'r') as f:
			t1=f.readlines()
		t1=[x.strip() for x in t1]
		t1=list(filter(lambda x:len(x)>0,t1))
		namec_satac.append(t1)
	snamec_srna,snamec_satac=[[set(y) for y in x] for x in [namec_srna,namec_satac]]

	#Checking cell names
	t1=list(itertools.chain.from_iterable(namec_srna))
	if len(t1)!=len(set(t1)):
		print(len(t1),len(set(t1)))
		raise ValueError('Found RNA cells assigned to multiple subsets')
	t1=set(t1)
	if len(t1-snamec_rna)>0:
		raise ValueError('Subset RNA cells not found in population')
	if len(snamec_rna-t1)>0:
		logging.warn('Found RNA cells not assigned to any subset')

	if isjoint:
		if any([frozenset(x)!=frozenset(y) for x,y in zip(snamec_srna,snamec_satac)]):
			raise ValueError('Subset assignments must be identical for joint profiles.')
	else:
		t1=list(itertools.chain.from_iterable(namec_satac))
		if len(t1)!=len(set(t1)):
			raise ValueError('Found ATAC cells assigned to multiple subsets')
		t1=set(t1)
		if len(t1-snamec_atac)>0:
			raise ValueError('Subset ATAC cells not found in population')
		if len(snamec_atac-t1)>0:
			logging.warn('Found ATAC cells not assigned to any subset')

#############################################
# Dynamic GRN inference checks
#############################################

if not pexists(pjoin(dirmakefiles,'dynamic.mk')) or not pexists(pjoin(dirdata,'traj_node.h5')):
	logging.warn('Cannot find dynamic.mk or traj_node.h5. Skipping dynamic network inference checks.')
else:
	#For RNA
	traj=trajectory.from_file(pjoin(dirdata,'traj_node.h5'))
	pt_rna=point.from_file(traj,pjoin(dirdata,'traj_cell_rna.h5'))
	coord_rna=pd.read_csv(pjoin(dirdata,'coord_rna.tsv.gz'),header=0,index_col=0,sep='\t')
	if len(set(coord_rna.index)-snamec_rna)>0:
		raise ValueError('Low dimensional RNA cells not found in population')
	if len(coord_rna)!=len(pt_rna):
		raise ValueError('Inconsistent RNA cell count between coord_rna.tsv.gz and traj_cell_rna.h5')
	#For ATAC
	if not isjoint:
		pt_atac=point.from_file(traj,pjoin(dirdata,'traj_cell_atac.h5'))
		coord_atac=pd.read_csv(pjoin(dirdata,'coord_atac.tsv.gz'),header=0,index_col=0,sep='\t')
		if len(set(coord_atac.index)-snamec_atac)>0:
			raise ValueError('Low dimensional ATAC cells not found in population')
		if len(coord_atac)!=len(pt_atac):
			raise ValueError('Inconsistent ATAC cell count between coord_atac.tsv.gz and traj_cell_atac.h5')
