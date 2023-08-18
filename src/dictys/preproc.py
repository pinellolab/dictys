#!/usr/bin/python3
# Lingfei Wang, 2020-2023. All rights reserved.

"""Preprocessing
"""

def selects_rna(fi_reads:str,fi_names:str,fo_reads:str)->None:
	"""
	Select samples/cells based on external table for RNA data.
	
	Parameters
	----------
	fi_reads:
		Path of input tsv file of full expression matrix
	fi_names:
		Path of input text file of sample/cell names to select
	fo_reads:
		Path of output tsv file of expression matrix of selected samples/cells
	"""
	import pandas as pd
	import logging
	from dictys.utils.file import read_txt
	#Loading data
	logging.info(f'Reading file {fi_reads}.')
	d=pd.read_csv(fi_reads,index_col=0,header=0,sep='\t')
	names=set(read_txt(fi_names,unique=True))

	#Filtering 
	t1=d.columns.isin(names)
	d=d[d.columns[t1]]
	if len(d)==0:
		raise RuntimeError('No samples/cells found.')
	t1=(d.values!=0).any(axis=1)
	d=d[t1]
	logging.info(f'Writing file {fo_reads}.')
	d.to_csv(fo_reads,index=True,header=True,sep='\t')

def selects_atac(fi_exp:str,fi_list:str,fo_list:str)->None:
	"""
	Select chromatin accessibility samples/cells based on presence in expression matrix.
	
	Parameters
	------------
	fi_exp:
		Path of input tsv file of expression matrix. Column must be sample/cell name.
	fi_list:
		Path of input text file of selected cell names, one per line
	fo_list:
		Path of output text file of selected cell names, one per line

	"""
	from os import linesep
	from dictys.utils.file import read_columns,read_txt

	ind=read_txt(fi_list,unique=True)
	dt=set(read_columns(fi_exp,unique=True))
	ind=list(filter(lambda x:len(x)>0 and x in dt,ind))
	if len(ind)==0:
		raise RuntimeError('No cell selected.')
	ind=linesep.join(ind)+linesep
	with open(fo_list,'w') as f:
		f.write(ind)

def qc_reads(fi_reads:str,fo_reads:str, n_gene:int, nc_gene:int, ncp_gene:float, n_cell:int, nt_cell:int, ntp_cell:float)->None:		# noqa: C901
	"""
	Quality control by bounding read counts.

	Quality control is perform separately on genes based on their cell statisics and on cells based on their gene statistics, iteratively until dataset remains unchanged. A gene or cell is removed if any of the QC criteria is violated at any time in the iteration. All QC parameters can be set to 0 to disable QC filtering for that criterion.

	Parameters
	-----------
	fi_reads:
		Path of input tsv file of read count matrix. Rows are genes and columns are cells.
	fo_reads:
		Path of output tsv file of read count matrix after QC
	n_gene:
		Lower bound on total read counts for gene QC
	nc_gene:
		Lower bound on number of expressed cells for gene QC
	ncp_gene:
		Lower bound on proportion of expressed cells for gene QC
	n_cell:
		Lower bound on total read counts for cell QC
	nt_cell:
		Lower bound on number of expressed genes for cell QC
	ntp_cell:
		Lower bound on proportion of expressed genes for cell QC

	"""
	import numpy as np
	import pandas as pd
	import logging
	
	logging.info(f'Reading file {fi_reads}.')
	reads0=pd.read_csv(fi_reads,header=0,index_col=0,sep='\t')
	reads=reads0.values
	if reads.ndim != 2:
		raise ValueError('reads must have 2 dimensions.')
	if not np.all([
		x >= 0 for x in [n_gene, nc_gene, ncp_gene, n_cell, nt_cell, ntp_cell]]):
		raise ValueError('All parameters must be non-negative.')
	if not np.all([x <= 1 for x in [ncp_gene, ntp_cell]]):
		raise ValueError('Proportional parameters must be no greater than 1.')

	dt = reads
	nt, ns = dt.shape
	nt0 = ns0 = 0
	st = np.arange(nt)
	ss = np.arange(ns)
	while nt0 != nt or ns0 != ns:
		nt0 = nt
		ns0 = ns
		st1 = np.ones(len(st), dtype=bool)
		ss1 = np.ones(len(ss), dtype=bool)
		# Filter genes
		if n_gene > 0:
			st1 &= dt.sum(axis=1) >= n_gene
		if nc_gene > 0 or ncp_gene > 0:
			t1 = (dt > 0).sum(axis=1)
			if nc_gene > 0:
				st1 &= t1 >= nc_gene
			if ncp_gene > 0:
				st1 &= t1 >= ncp_gene * ns
		# Filter cells
		if n_cell > 0:
			ss1 &= dt.sum(axis=0) >= n_cell
		if nt_cell > 0 or ntp_cell > 0:
			t1 = (dt > 0).sum(axis=0)
			if nt_cell > 0:
				ss1 &= t1 >= nt_cell
			if ntp_cell > 0:
				ss1 &= t1 >= ntp_cell * nt
		# Removals
		st = st[st1]
		ss = ss[ss1]
		dt = dt[st1][:, ss1]
		nt = len(st)
		ns = len(ss)
		if nt == 0:
			raise RuntimeError('All genes removed in QC.')
		if ns == 0:
			raise RuntimeError('All cells removed in QC.')
	logging.info('Removed {}/{} genes and {}/{} cells in QC.'.format(
		reads.shape[0] - len(st), reads.shape[0], reads.shape[1] - len(ss),
		reads.shape[1]))
	reads0=reads0.iloc[st,ss]
	logging.info(f'Writing file {fo_reads}.')
	reads0.to_csv(fo_reads,header=True,index=True,sep='\t')
	























































#
