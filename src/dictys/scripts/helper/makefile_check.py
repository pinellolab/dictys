#!/usr/bin/env python3
# Lingfei Wang, 2022, 2023. All rights reserved.
import json
import argparse
import sys
from os.path import join as pjoin
from os.path import exists as pexists
from os.path import isdir
from os import listdir
from functools import reduce
from operator import add
from collections import Counter
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

if not isdir(dirmakefiles):
	raise FileNotFoundError('Cannot find makefiles folder at {}. Please use option --dir_makefiles to specify makefiles folder.'.format(dirmakefiles))
if not isdir(dirdata):
	raise FileNotFoundError('Cannot find input data folder at {}. Please use option --dir_data to specify input data folder.'.format(dirdata))

#Detect whether joint profiles from makefile
with open(pjoin(dirmakefiles,'config.mk'),'r') as f:
	s=f.readlines()
s=[x.strip() for x in s if x.startswith('JOINT=')][-1]
s=s[len('JOINT='):]
if s not in {'0','1'}:
	raise ValueError('Invalid JOINT variable in {}. Only accepts: 0 for separate quantification of single-cell transcriptome and chromatin accessibility, 1 for joint quantification of single-cell transcriptome and chromatin accessibility.'.format(pjoin(dirmakefiles,'config.mk')))
isjoint=s=='1'
print(f'Joint profile: {isjoint}')

#Gene & cell names from RNA
namec_rna=pd.read_csv(pjoin(dirdata,'expression.tsv.gz'),header=0,index_col=0,sep='\t',nrows=1)
namec_rna=np.array(namec_rna.columns)
snamec_rna=set(namec_rna)
print(f'Found {len(namec_rna)} cells with RNA profile')
if len(namec_rna)<100:
	raise ValueError('<100 cells found with RNA profile in expression.tsv.gz')
nameg=pd.read_csv(pjoin(dirdata,'expression.tsv.gz'),header=0,index_col=0,sep='\t',usecols=[0,1])
nameg=np.array(nameg.index)
snameg=set(nameg)
print(f'Found {len(nameg)} genes with RNA profile')
if len(nameg)<100:
	raise ValueError('<100 genes found with RNA profile in expression.tsv.gz')

#Cell names from ATAC
namec_atac=listdir(pjoin(dirdata,'bams'))
namec_atac=np.array(['.'.join(x.split('.')[:-1]) for x in namec_atac if x.endswith('.bam')])
snamec_atac=set(namec_atac)
print(f'Found {len(namec_atac)} cells with ATAC profile')
t1=snamec_rna-snamec_atac
if isjoint and len(t1)>0:
	raise FileNotFoundError('Not all cells with RNA profile in expression.tsv.gz has a bam file in bams folder for the joint profiling dataset. First three cells missing: '+', '.join(list(t1)[:3]))

#TF names from motifs
with open(pjoin(dirdata,'motifs.motif'),'r') as f:
	nameg_motif=f.readlines()
nameg_motif=[x.split('\t')[1] for x in nameg_motif if x.startswith('>')]
if len(nameg_motif)!=len(set(nameg_motif)):
	raise ValueError('Duplicate motif names found. Please check your motifs.motif file.')
print(f'Found {len(nameg_motif)} motifs')
nameg_motif=set([x.split('_')[0] for x in nameg_motif])
print(f'Found {len(nameg_motif)} TFs')
nameg_motif_notfound=sorted(nameg_motif-set(nameg))
nameg_motif=sorted(nameg_motif&set(nameg))
snameg_motif=set(nameg_motif)
print(f'Found {len(nameg_motif)} TFs in current dataset')
print(f'Missing {len(nameg_motif_notfound)} TFs in current dataset: '+','.join(nameg_motif_notfound))
if len(nameg_motif)<10:
	raise ValueError('<10 TFs found in motifs.motif file.')

#Chromosomes and gene names from bed file
with open(pjoin(dirdata,'gene.bed'),'r') as f:
	nameg_bed=f.readlines()
nameg_bed=[x.strip() for x in nameg_bed]
nameg_bed=[x.split('\t') for x in nameg_bed if len(x)>0]
if len(set([len(x) for x in nameg_bed]))!=1:
	raise ValueError('Unequal number of columns in gene.bed')
nameg_bed=list(filter(lambda x:x[3] in snameg,nameg_bed))
if len(nameg_bed[0])<6:
	logging.warning('No strand information in gene.bed')
namechr_bed=sorted(set([x[0] for x in nameg_bed]))
snamechr_bed=set(namechr_bed)
nameg_bed=[x[3] for x in nameg_bed]
t1=[x[0] for x in dict(Counter(nameg_bed)).items() if x[1]>1][:3]
if len(t1)>0:
	raise ValueError('Duplicate gene names found in gene.bed. First three duplicate gene names: '+', '.join(t1))
nameg_bed=sorted(nameg_bed)
print(f'Found {len(nameg_bed)} genes with TSS information')
if len(nameg_bed)<100:
	raise ValueError('<100 overlapping genes with TSS information between gene.bed and expression.tsv.gz')

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
	if len(snamechr_bl-set(namechr_bl))>0:
		logging.warning('Blacklist chromosomes in blacklist.bed not found in reference genome: '+', '.join(sorted(snamechr_bl-namechr_bl)))

#Covariates
if pexists(pjoin(dirdata,'covariate.tsv.gz')):
	dcov=pd.read_csv(pjoin(dirdata,'covariate.tsv.gz'),header=0,index_col=0,sep='\t')
	ncov=dcov.shape[1]
	print('Found {} covariates'.format(ncov))
	t1=[x[0] for x in dict(Counter(dcov.columns)).items() if x[1]>1][:3]
	if len(t1)>0:
		raise ValueError('Duplicate covariate names found in covariate.tsv.gz. First three duplicate covariate names: '+', '.join(t1))
	t1=[x[0] for x in dict(Counter(dcov.index)).items() if x[1]>1][:3]
	if len(t1)>0:
		raise ValueError('Duplicate cell names found in covariate.tsv.gz. First three duplicate cell names: '+', '.join(t1))
	t1=list(snamec_rna-set(dcov.index))[:3]
	if len(t1)>0:
		raise ValueError('Covariate information cannot be found in covariate.tsv.gz for cells from expression.tsv.gz. First three missing cells: '+', '.join(t1))

#############################################
# Context specific GRN inference checks
#############################################

if not pexists(pjoin(dirmakefiles,'static.mk')) or not pexists(pjoin(dirdata,'subsets.txt')):
	logging.warning('Cannot find static.mk or subsets.txt. Skipping static network inference checks.')
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
	snamec_srna,snamec_satac=[[Counter(y) for y in x] for x in [namec_srna,namec_satac]]

	#Checking cell names
	t1=reduce(add,snamec_srna)
	t2=[x[0] for x in dict(t1).items() if x[1]>1][:3]
	if len(t2)>0:
		t2=[(x,list(itertools.chain.from_iterable([[names[y]]*snamec_srna[y][x] for y in range(len(names))]))) for x in t2]
		t2='; '.join(['{}: {}'.format(x[0],','.join(x[1])) for x in t2])
		raise ValueError('Found RNA cells in appearing multiple times in subsets at subsets/*/names_rna.txt. First three cells and their assigned subsets: '+t2)
	t2=list(set(t1)-snamec_rna)[:3]
	if len(t2)>0:
		raise ValueError('Subset RNA cells in subsets/*/names_rna.txt could not be found in expression.tsv.gz. First three missing cells: '+', '.join(t2))
	t2=list(snamec_rna-set(t1))[:3]
	if len(t2)>0:
		logging.warning('Found RNA cells in expression.tsv.gz not assigned to any subset according to subsets/*/names_rna.txt. First three missing cells: '+', '.join(t2))

	if isjoint:
		if any([frozenset(x)!=frozenset(y) for x,y in zip(snamec_srna,snamec_satac)]):
			raise ValueError('Subset assignments files subsets/*/names_rna.txt and subsets/*/names_atac.txt must be identical in each subset for joint profiles.')
	else:
		t1=reduce(add,snamec_satac)
		t2=[x[0] for x in dict(t1).items() if x[1]>1][:3]
		if len(t2)>0:
			t2=[(x,list(itertools.chain.from_iterable([[names[y]]*snamec_satac[y][x] for y in range(len(names))]))) for x in t2]
			t2='; '.join(['{}: {}'.format(x[0],','.join(x[1])) for x in t2])
			raise ValueError('Found ATAC cells appearing multiple times in subsets at subsets/*/names_atac.txt. First three cells and their assigned subsets: '+t2)
		t2=list(set(t1)-snamec_atac)[:3]
		if len(t2)>0:
			raise ValueError('Subset ATAC cells in subsets/*/names_atac.txt could not be found in bams folder. First three missing cells: '+', '.join(t2))
		t2=list(snamec_atac-set(t1))[:3]
		if len(t2)>0:
			logging.warning('Found ATAC cells in bams folder not assigned to any subset according to subsets/*/names_atac.txt. First three missing cells: '+', '.join(t2))

#############################################
# Dynamic GRN inference checks
#############################################

if not pexists(pjoin(dirmakefiles,'dynamic.mk')) or not pexists(pjoin(dirdata,'traj_node.h5')):
	logging.warning('Cannot find dynamic.mk or traj_node.h5. Skipping dynamic network inference checks.')
else:
	#For RNA
	traj=trajectory.from_file(pjoin(dirdata,'traj_node.h5'))
	pt_rna=point.from_file(traj,pjoin(dirdata,'traj_cell_rna.h5'))
	coord_rna=pd.read_csv(pjoin(dirdata,'coord_rna.tsv.gz'),header=0,index_col=0,sep='\t')
	t1=list(set(coord_rna.index)-snamec_rna)[:3]
	if len(t1)>0:
		raise ValueError('Low dimensional RNA cells in coord_rna.tsv.gz cannot be found in expression.tsv.gz. First three missing cells: '+', '.join(t1))
	if len(coord_rna)!=len(pt_rna):
		raise ValueError('Inconsistent RNA cell count between coord_rna.tsv.gz and traj_cell_rna.h5')
	#For ATAC
	if not isjoint:
		pt_atac=point.from_file(traj,pjoin(dirdata,'traj_cell_atac.h5'))
		coord_atac=pd.read_csv(pjoin(dirdata,'coord_atac.tsv.gz'),header=0,index_col=0,sep='\t')
		t1=list(set(coord_atac.index)-snamec_atac)[:3]
		if len(t1)>0:
			raise ValueError('Low dimensional ATAC cells in coord_atac.tsv.gz cannot be found in bams folder. First three missing cells: '+', '.join(t1))
		if len(coord_atac)!=len(pt_atac):
			raise ValueError('Inconsistent ATAC cell count between coord_atac.tsv.gz and traj_cell_atac.h5')
