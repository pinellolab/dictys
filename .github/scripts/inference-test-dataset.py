# Lingfei Wang, 2022. All rights reserved.
# Usage: python3 $0 expected_folder output_folder expected_h5 output_h5 [exclusion_filename_1 ...]
# Performs comparison between expected and actual outputs in network inference
import sys
import logging
logging.basicConfig(level=logging.INFO)
import numpy as np
import pandas as pd
from os.path import join as pjoin
from os.path import isdir
import os
from os import linesep, listdir
from collections import defaultdict

class TestError(RuntimeError):
	pass

#######################################################
# Comparison functions
#######################################################

def test(statement,*a,**ka):
	if not statement:
		raise TestError(*a,**ka)

def cmp_txt(f1,f2):
	d=[]
	for xi in [f1,f2]:
		with open(f1,'r') as f:
			d.append(f.readlines())
	d=[set([x.strip() for x in y])-{''} for y in d]
	test(d[0]==d[1])
def cmp_tsv_gz(f1,f2,postproc=lambda u,v:[u,v],unordered_row=False,header=0,index_col=0,sep='\t',**ka):
	import pandas as pd
	d=[]
	for xi in [f1,f2]:
		d.append(pd.read_csv(xi,header=header,index_col=index_col,sep=sep,**ka))
	d=postproc(*d)
	#Test dimensions
	t1=[set(x.index) for x in d]
	test(len(t1[0])==len(t1[1]),'Different row counts')
	test(set(t1[0])==set(t1[1]),'Different rows')
	t1=[set(x.columns) for x in d]
	test(len(t1[0])==len(t1[1]),'Different column counts')
	test(set(t1[0])==set(t1[1]),'Different columns')
	#Align dimensions
	d[1]=d[1][d[0].columns]
	if not unordered_row:
		d[1]=d[1].loc[d[0].index]
	else:
		[x.sort_values(list(x.columns),axis=0,inplace=True,ignore_index=True) for x in d]
	for xi in d[0].columns:
		#Comparison depends on column dtype
		if any(np.issubdtype(x[xi].dtype,object) for x in d):
			#Object comparison
			t1=d[0][xi]==d[1][xi]
			t1|=pd.isna(d[0][xi])&pd.isna(d[1][xi])
			test(t1.all(),f'Different values in object column {xi}')
		elif any(np.issubdtype(x[xi].dtype,np.floating) for x in d):
			#Float comparison
			t1=np.isclose(d[0][xi],d[1][xi])
			test(t1.all(),f'Different values in float column {xi}')
		else:
			test(np.issubdtype(d[0][xi].dtype,d[1][xi].dtype) or np.issubdtype(d[1][xi].dtype,d[0][xi].dtype),f'Different dtypes in column {xi}')
			test((d[0][xi]==d[1][xi]).all(),'Different values in {} column {}'.format(d[0][xi].dtype,xi))
def cmp_binary(f1,f2):
	d=[]
	for xi in [f1,f2]:
		with open(f1,'rb') as f:
			d.append(f.read())
	test(len(d[0])==len(d[1]),'Different sizes')
	test(d[0]==d[1],'Different contents')
def cmp_bed(f1,f2,**ka):
	return cmp_tsv_gz(f1,f2,header=None,index_col=None,sep='\t',**ka)

#######################################################
# File postprocessing and exclusion
#######################################################

#Keyword arguments for each file tested
ka_file=defaultdict(dict)
ka_file['binding.tsv.gz']={
	'index_col':None,
	'unordered_row':True,
}
ka_file['tssdist.tsv.gz']={
	'index_col':None,
	'unordered_row':True,
}
ka_file['footprints.bed']={
	'postproc':lambda u,v:[u[[0,1,2,4]],v[[0,1,2,4]]],
}

#######################################################
# Testing
#######################################################

assert len(sys.argv)>=5

diris=sys.argv[1:3]
h5s=sys.argv[3:5]
exclusions=set(x.strip() for x in sys.argv[5:])

err=False
celltypes=[list(filter(lambda y:isdir(pjoin(x,y)),listdir(x))) for x in diris]
assert all(set(celltypes[0])==set(x) for x in celltypes[1:])
celltypes=celltypes[0]
for celltype in celltypes:
	files=[listdir(pjoin(x,celltype)) for x in diris]
	t1=set(files[0])-set(files[1])-exclusions
	if len(t1)>0:
		err=True
		logging.warning('{}: Expected output files missing: '.format(celltype)+','.join(sorted(list(t1))))
	testlist=list(filter(lambda x:x not in t1,files[0]))
	for file in testlist:
		fname=pjoin(celltype,file)
		if file in exclusions:
			logging.info(f'Excluded: {fname}')
			continue
		logging.info(f'Started: {fname}')
		fs=[pjoin(x,fname) for x in diris]
		try:
			if file.endswith('.bai'):
				cmp_binary(*fs,**ka_file[file])
			elif file.endswith('.bam'):
				cmp_binary(*fs,**ka_file[file])
			elif file.endswith('.bed'):
				cmp_bed(*fs,**ka_file[file])
			elif file.endswith('.tsv.gz'):
				cmp_tsv_gz(*fs,**ka_file[file])
			elif file.endswith('.txt'):
				cmp_txt(*fs,**ka_file[file])
			else:
				raise TestError(f'Not implemented: {fname}')
			logging.info(f'Passed: {fname}')
		except TestError as e:
			err=True
			logging.warning(f'{fname}: '+' '.join(e.args))

if all(len(x)>0 for x in h5s):
	logging.info(f'Started: {h5s[0]} vs {h5s[1]}')
	try:
		cmp_binary(*h5s)
		logging.info(f'Passed: {h5s[0]} vs {h5s[1]}')
	except TestError as e:
		err=True
		logging.warning(f'{h5s[0]} vs {h5s[1]}: '+' '.join(e.args))
else:
	logging.info('Excluded: final h5 network file')

if err:
	raise TestError('At least one test failed')











#
