#!/usr/bin/python3
# Lingfei Wang, Nikolaos Trasanidis, 2022, 2023. All rights reserved.
"""Chromatin accessibility analyses
"""

from typing import Optional,Union
from dictys.utils.numpy import NDArray

################################################################
# Peak calling
################################################################

def macs2(fi_names:str,fi_bam:str,fo_bam:str,fo_bai:str,fo_bed:str,genome_size:str,qcut:float=0.05,nth:int=1,nmax:int=500000)->None:
	"""
	Peak calling using macs2.

	Needs bam files for each cell in a given folder.
	
	Parameters
	------------
	fi_names:
		Path of input text file containing one sample/cell name per line for macs2 peak calling
	fi_bam:
		Path of input folder that contains each cell's bam file by name in *fi_names*
	fo_bam:
		Path of output bam file for select samples/cells
	fo_bai:
		Path of output bai file for select samples/cells
	fo_bed:
		Path of output bed file of peaks
	genome_size:
		Genome size input of macs2. Use shortcuts hs or mm for human or mouse.
	qcut:
		Qvalue cutoff for macs2
	nth:
		Number of threads
	nmax:
		Maximum number of peaks to retain, ordered by macs2 score. Use 0 for no limit.
		
	"""
	from os.path import dirname,basename,abspath,isdir,isfile
	from os.path import join as pjoin
	from os import linesep
	import numpy as np
	import pandas as pd
	import logging
	from dictys.utils import shell
	from dictys.utils.file import read_txt
	if qcut<=0 or qcut>=1:
		raise ValueError('qcut must be between 0 and 1.')
	if not isfile(fi_names):
		raise FileNotFoundError(fi_names)
	if not isdir(fi_bam):
		raise FileNotFoundError(fi_bam)
	
	scriptpath=pjoin(abspath(dirname(__file__)),'scripts',basename(__file__)[:-3]+'_macs2.sh')
	fi_names,fi_bam,fo_bam,fo_bai,fo_bed=[abspath(x) for x in [fi_names,fi_bam,fo_bam,fo_bai,fo_bed]]
	
	#Load sample names
	names=read_txt(fi_names,unique=True)
	if len(names)==0:
		raise ValueError('No sample name found in '+fi_names)
	namestxt=linesep.join([x+'.bam' for x in names])+linesep
	
	#Run script for macs2
	cmd = scriptpath+f" cellnames.txt {fi_bam} {fo_bam} {fo_bai} {fo_bed} {genome_size} {qcut} {nth}"
	d2 = shell.cmdfile(cmd,[],infiles={'cellnames.txt': namestxt},quiet=False,cd=True)
	if d2 is None or len(d2)>0 or not all(isfile(x) for x in [fo_bam,fo_bai,fo_bed]):
		raise RuntimeError('Macs2 failed.')

	if nmax==0:
		return
	#Reduce size of peak bed file
	logging.info(f'Reading file {fo_bed}')
	d3=pd.read_csv(fo_bed,header=None,index_col=None,sep='\t')
	if len(d3)<=nmax:
		return
	d4=np.partition(d3[8].values,-nmax-1)[-nmax-1]
	d3=d3[d3[8]>d4]
	logging.info(f'Writing file {fo_bed}')
	d3.to_csv(fo_bed,header=False,index=False,sep='\t')

################################################################
# TF footprinting
################################################################

def wellington(fi_bam:str,fi_bai:str,fi_bed:str,fo_bed:str,fi_blacklist:Optional[str]=None,cut:float=10,nth:int=1,nmax:int=100000)->None:
	"""
	TF footprinting with wellington.
	
	Parameters
	------------
	fi_bam:
		Path of input bam file of all reads
	fi_bai:
		Path of input bai file of all reads
	fi_bed:
		Path of input bed file of peaks
	fo_bed:
		Path of output bed file of footprints
	fi_blacklist:
		Path of input bed file of blacklisted genome regions to be removed
	cut:
		Cutoff for wellington score
	nth:
		Number of threads
	nmax:	
		Maximum number of footprints to retain, ordered by wellington score. Use 0 for no limit.
		
	"""
	from os.path import dirname, basename, abspath, isfile
	from os.path import join as pjoin
	from .utils import shell
	if fi_blacklist is None:
		fi_blacklist='None'
	else:
		if not isfile(fi_blacklist):
			raise FileNotFoundError(fi_blacklist)
		fi_blacklist=abspath(fi_blacklist)
	if cut<=0:
		raise ValueError('cut must be positive.')
	for xi in [fi_bam,fi_bai,fi_bed]:
		if not isfile(xi):
			raise FileNotFoundError(xi)
	
	scriptpath=pjoin(abspath(dirname(__file__)),'scripts',basename(__file__)[:-3]+'_wellington.sh')
	fi_bam,fi_bai,fi_bed,fo_bed=[abspath(x) for x in [fi_bam,fi_bai,fi_bed,fo_bed]]
	cmd = scriptpath+f" {fi_bam} {fi_bai} {fi_bed} {fo_bed} {cut} {nth} {nmax} {fi_blacklist}"
	d2 = shell.cmdfile(cmd,[],quiet=False,cd=True)
	if d2 is None or len(d2)>0 or not all(isfile(x) for x in [fo_bed]):
		raise RuntimeError('Wellington failed.')

################################################################
# Motif scan
################################################################

def _motif_postproc(dret,fi_exp:str,fo_bed:str,fo_wellington:str,fo_homer:str)->None:
	"""
	Postprocess motif discovery results from HOMER.
	dret:	Return tuple of HOMER call
	fi_exp:	Expression matrix file path.
	Return:
	Processed and copied new data object with motif results
	"""
	import numpy as np
	import pandas as pd
	import logging
	from io import StringIO
	from dictys.utils.file import read_index
	if dret is None or len(dret)!=3:
		raise RuntimeError('Homer failed.')
	with StringIO(dret[0].decode()) as f:
		dw = pd.read_csv(f, header=0, index_col=0, sep='\t')
	with StringIO(dret[1].decode()) as f:
		dh = pd.read_csv(f, header=0, index_col=0, sep='\t')
	with StringIO(dret[2].decode()) as f:
		try:
			dmotif = pd.read_csv(f, header=None, index_col=None, sep='\t')
		except pd.errors.EmptyDataError:
			dmotif=None
	if dmotif is None:
		raise RuntimeError('No motif found.')
	assert dw.shape==dh.shape
	assert (dw.columns==dh.columns).all()
	assert (dw.index==dh.index).all()
	assert dmotif.shape[1]==7
	if dmotif.shape[0]==0:
		raise RuntimeError('No motif found.')

	namet=np.array(read_index(fi_exp,unique=True))
	#Set na as or function
	t1=(pd.isna(dw)|pd.isna(dh)).values
	dw.fillna(0,inplace=True)
	dh.fillna(0,inplace=True)
	#Extract dimensions
	namep=np.array([str(x) for x in dw.index])
	namem=np.array(list(dw.columns))
	#Set data values
	dw,dh=[x.values.astype(float) for x in [dw,dh]]
	dw[t1]=0
	dh[t1]=0
	#Remove TFs in motifs not found in current dataset
	t1=set(namet)
	namem=[x.split('_') for x in namem]
	namem=['_'.join([','.join(list(filter(lambda y:y in t1,x[0].split(','))))]+x[1:]) for x in namem]
	namem=np.array(namem)
	#Remove motifs having no TF in current dataset
	t1=[not x.startswith('_') for x in namem]
	dw,dh=[x[:,t1] for x in [dw,dh]]
	namem=namem[t1]
	if len(namem)!=len(set(namem)):
		from collections import Counter
		t1=[x[0] for x in Counter(namem).items() if x[1]>1][:3]
		raise ValueError('Found non-unique motif name suffices. Each motif name is recommended to contain a unique suffix. First three non-unique motif names: {}'.format(', '.join(t1)))
	assert dw.shape==(len(namep),len(namem)) and dh.shape==(len(namep),len(namem))
	assert np.isfinite(dw).all() and np.isfinite(dh).all()
	assert (dw>=0).all() and (dh>=0).all()
	#Remove regions not mapped to motif
	t1=(dw>0).any(axis=1)|(dh>0).any(axis=1)
	dw,dh,namep=[x[t1] for x in [dw,dh,namep]]
	#Output
	logging.info(f'Writing file {fo_bed}')
	dmotif.to_csv(fo_bed,header=False,index=False,sep='\t')
	dw=pd.DataFrame(dw,index=namep,columns=namem)
	logging.info(f'Writing file {fo_wellington}')
	dw.to_csv(fo_wellington,header=True,index=True,sep='\t')
	dh=pd.DataFrame(dh,index=namep,columns=namem)
	logging.info(f'Writing file {fo_homer}')
	dh.to_csv(fo_homer,header=True,index=True,sep='\t')

def homer(fi_bed:str,fi_motif:str,dirio_genome:str,fi_exp:str,fo_bed:str,fo_wellington:str,fo_homer:str,nth:int=1)->None:
	"""
	Motif scan with homer.
	
	Parameters
	------------
	fi_bed:	
		Path of input bed file of regions
	fi_motif:
		Path of input motif PWM file in homer format. Motifs must be named in format 'gene_...' where gene matches gene names in fi_exp. Should not contain duplicates.
	dirio_genome:
		Path of input & output folder or file for reference genome for homer. A separate hard copy is recommended because homer may write into the folder to preparse genome.
	fi_exp:
		Path of input expression matrix file in tsv format. Used for mapping motifs to genes.
	fo_bed:	
		Path of output bed file of detected motifs
	fo_wellington:
		Path of output tsv file of wellington scores in shape (region,motif)
	fo_homer:
		Path of output tsv file of homer scores in shape (region,motif)
	nth:	
		Number of threads
	"""
	from os.path import dirname, basename, abspath, isfile, isdir
	from os.path import join as pjoin
	from .utils import shell
	for xi in [fi_bed,fi_motif]:
		if not isfile(xi):
			raise FileNotFoundError(xi)
	if not (isfile(dirio_genome) or isdir(dirio_genome)):
		raise FileNotFoundError(dirio_genome)
	
	scriptpath=pjoin(abspath(dirname(__file__)),'scripts',basename(__file__)[:-3]+'_homer.sh')
	spath=pjoin(abspath(dirname(__file__)),'scripts',basename(__file__)[:-3]+'_homer.py')
	fi_bed,fi_motif,dirio_genome,fo_bed,fo_wellington,fo_homer=[abspath(x) for x in [fi_bed,fi_motif,dirio_genome,fo_bed,fo_wellington,fo_homer]]
	cmd = scriptpath+f" {fi_bed} {fi_motif} {dirio_genome} {spath} {nth}"
	d2 = shell.cmdfile(cmd,
		['19-w.tsv','19-h.tsv','16-long.bed'],
		quiet=False,cd=True,sizelimit=None)
	return _motif_postproc(d2,fi_exp,fo_bed,fo_wellington,fo_homer)

################################################################
# Linking TFs to target genes
################################################################

def _linking_score(vw:Union[float,NDArray[float]],vh:Union[float,NDArray[float]],dist:Union[float,NDArray[float]],mode:int=7):
	"""Footprinting I score
	mode:
		Mode to compute final score. Accepts binary flags:
		* 1:	Add log(wellington score)
		* 2:	Add log(homer score)
		* 4:	Subtract log(10)*(distance_to_tss)/1E6
	"""
	import numpy as np
	ans=0
	if mode&1:
		ans=ans+np.log(vw)
	if mode&2:
		ans=ans+np.log(vh)
	if mode&4:
		ans=ans-np.log(10)*np.abs(dist)/1E6
	return ans

def binding(fi_wellington:str,fi_homer:str,fo_bind:str,cuth:float=0,cutw:float=0,cut:Optional[float]=None,combine:str='max',mode:int=3)->None:
	"""
	Finding TF binding events.

	Combines wellington and homer outputs to infer TF binding events by merging motifs to TFs.

	Parameters
	----------
	fi_wellington:
		Path of input tsv file of wellington output
	fi_homer:
		Path of input tsv file of homer output
	fo_bind:
		Path of output tsv file of binding events
	cuth:
		Homer score cutoff
	cutw:
		Wellington score cutoff
	cut:
		Final score (integrating homer & wellington) cutoff
	combine:
		Method to combine scores of motifs of the same TF. Accepts: max, mean, sum.
	mode:
		Mode to compute final score. Accepts binary flags:

		* 1:	Add log(wellington score)

		* 2:	Add log(homer score)

		* 4:	Subtract log(10)*(distance_to_tss)/1E6
	"""
	import numpy as np
	import pandas as pd
	from dictys.utils.numpy import groupby
	import logging
	if combine=='max':
		combine=np.max
	elif combine=='sum':
		combine=np.sum
	elif combine=='mean':
		combine=np.mean
	else:
		raise ValueError(f'Unknown combine method: {combine}')

	# Read in files
	logging.info(f'Reading file {fi_wellington}')
	wellington_df=pd.read_csv(fi_wellington,header=0,index_col=0,sep='\t')
	logging.info(f'Reading file {fi_homer}')
	homer_df=pd.read_csv(fi_homer,header=0,index_col=0,sep='\t')
	assert (wellington_df.index==homer_df.index).all()
	assert (wellington_df.columns==homer_df.columns).all()
	vw=wellington_df.values
	vh=homer_df.values
	namem=np.array(list(wellington_df.columns))
	#Get motif to gene map
	dmt=groupby([x.split('_')[0] for x in namem])

	mask=(vw>cutw)&(vh>cuth)
	binds=[]
	for xi in dmt:
		t1=np.nonzero(mask[:,dmt[xi]])
		if len(t1)==0:
			continue
		t1=(t1[0],dmt[xi][t1[1]])
		ds=_linking_score(vw[t1],vh[t1],None,mode=mode)
		if cut is not None:
			t2=ds>cut
			t1=tuple(x[t2] for x in t1)
			if len(t1)==0:
				continue
			ds=ds[t2]
		t3=groupby(t1[0])
		for xj in t3:
			binds.append([xi,xj,combine(ds[t3[xj]])])
	if len(binds)==0:
		raise RuntimeError('No TF binding relation remains.')
	binds=[list(x) for x in zip(*binds)]
	ans=pd.DataFrame([])
	ans['TF']=binds[0]
	ans['loc']=homer_df.index[binds[1]]
	ans['score']=binds[2]
	logging.info(f'Writing file {fo_bind}')
	ans.to_csv(fo_bind,index=False,header=True,sep='\t')

def tssdist(fi_exp:str,fi_wellington:str,fi_tss:str,fo_dist:str,cut:int=500000,nmin:int=1,nmax:int=10000000)->None:
	"""
	Annotating TF bond regions to target genes based on distance to TSS.

	Parameters
	----------
	fi_exp:
		Path of input expression matrix file in tsv format to obtain gene names
	fi_wellington:
		Path of input tsv file of wellington scores to obtain DNA regions
	fi_tss:
		Path of input bed file for gene region and strand
	fo_dist:
		Path of output tsv file of distance from TF-bond regions to TSS
	cut:
		Distance cutoff between DNA region and target gene TSS
	nmin:
		Minimal total number of links to recover
	nmax:
		Maximal total number of links to recover
	"""
	import numpy as np
	import pandas as pd
	import logging
	from pybedtools import BedTool
	from dictys.utils.file import read_index
	assert nmin>0
	if nmin>nmax:
		raise ValueError('nmax must be greater than nmin.')

	namet=np.array(read_index(fi_exp,unique=True))
	logging.info(f'Reading file {fi_wellington}')
	namep=np.array(read_index(fi_wellington,unique=True))
	
	#Produce peak bed file
	peaks=pd.DataFrame([x.split(':') for x in namep],index=None,columns=None)
	peaks=BedTool(peaks.to_csv(None,header=False,index=False,sep='\t'),from_string=True)
	#Produce target gene bed file
	logging.info(f'Reading file {fi_tss}')
	genes=pd.read_csv(fi_tss,header=None,index_col=None,sep='\t')
	t1=set(namet)
	#Shrink to genes showed up
	genes=genes[genes[3].isin(t1)]
	if genes.shape[1]<5:
		genes[4]='.'
	if genes.shape[1]<6:
		genes[5]='+'
		signed=False
	else:
		signed=True
	#Simplify to bed file
	assert genes[5].isin({'+','-'}).all()
	genes=genes.copy()
	t1=genes.apply(lambda x:x[1] if x[5]=='+' else x[2],axis=1)
	genes[1]=t1
	genes[2]=t1
	genes.drop_duplicates(inplace=True)
	assert len(genes[3].unique())==len(genes)
	genes=genes.to_csv(None,sep='\t',header=False,index=False)
	genes=BedTool(genes,from_string=True)

	#Filter pairs
	ans=genes.window(peaks,w=cut)		# pylint: disable=E1123
	try:
		ans2=ans.to_dataframe(header=None,disable_auto_names=True)
	except pd.errors.EmptyDataError:
		ans2=[]
	if len(ans2)<nmin:
		raise RuntimeError('Too few region-target relations found: {}'.format(len(ans2)))
	if len(ans2)>nmax:
		raise RuntimeError('Too many region-target relations found: {}'.format(len(ans2)))
	assert ans2[5].isin({'+','-'}).all()
	if signed and ans2[5].size>=100 and np.abs((ans2[5]=='+').mean()-0.5)>0.4:
		logging.warning('Heavily unbalanced strand found.')
	#Simplify output
	ans2['region']=ans2.apply(lambda x:'{}:{}:{}'.format(x[6],x[7],x[8]),axis=1)
	ans2['target']=ans2[3]
	ans2['dist1']=ans2.apply(lambda x:x[7]-x[1] if x[5]=='+' else x[2]-x[7],axis=1)
	ans2['dist2']=ans2.apply(lambda x:x[8]-x[1] if x[5]=='+' else x[2]-x[8],axis=1)
	t1=np.array([ans2['dist1'].values,ans2['dist2'].values])
	t2=t1[np.argmin(np.abs(t1),axis=0),np.arange(t1.shape[1])]
	t2[np.sign(t1).prod(axis=0)<=0]=0
	ans2['dist']=-t2
	ans2=ans2[['region','target','dist']].copy()
	logging.info(f'Writing file {fo_dist}')
	ans2.to_csv(fo_dist,header=True,index=False,sep='\t')

def linking(fi_binding:str,fi_dist:str,fo_linking:str,fi_whitelist:Optional[str]=None,whitelist_mode:str='intersect',combine:str='max',mode:int=4)->None:
	"""
	Linking regulators and targets with scores.

	Parameters
	----------
	fi_binding:
		Path of input tsv file of binding events
	fi_dist:
		Path of input tsv file of distance from TF-bond regions to TSS
	fo_linking:
		Path of output matrix file of TF to potential target gene link scores
	fi_whitelist:
		Path of input bed file of potential regulatory target genes of each region.
		The fourth column (or its first item after split by _) should be taget gene name.
		Can be used to filter regulatory regions based on co-accessibility or association with target gene expression.
	whitelist_mode:
		Criterion of whitelist. Accepts:
			* intersect: the region must intersect with a whitelisted region (default)
			* within: the region must be inside a whitelisted region
	combine:
		Method to combine scores of motifs of the same TF. Accepts: max, mean, sum.
	mode:
		Mode to compute final score. Accepts binary flags:

		* 4:	Subtract log(10)*(distance_to_tss)/1E6
	"""
	import numpy as np
	import pandas as pd
	import itertools
	from dictys.utils.numpy import groupby
	import logging
	if combine=='max':
		combine=np.max
	elif combine=='sum':
		combine=np.sum
	elif combine=='mean':
		combine=np.mean
	else:
		raise ValueError(f'Unknown combine method: {combine}')

	logging.info(f'Reading file {fi_binding}')
	db=pd.read_csv(fi_binding,header=0,index_col=None,sep='\t')
	logging.info(f'Reading file {fi_dist}')
	dd=pd.read_csv(fi_dist,header=0,index_col=None,sep='\t')
	
	if fi_whitelist is not None:
		#Filter distance to TSS table
		dw=pd.read_csv(fi_whitelist,header=None,index_col=None,sep='\t')
		dw[3]=dw[3].apply(lambda x:x.split('_')[0])
		dd['chr']=dd['region'].apply(lambda x:x.split(':')[0])
		dd['start']=dd['region'].apply(lambda x:x.split(':')[1]).astype(int)
		dd['stop']=dd['region'].apply(lambda x:x.split(':')[2]).astype(int)
		groupd=dd.groupby(['chr','target']).groups
		groupw=dw.groupby([0,3]).groups
		ans=[]
		for grp,target in set(groupd)&set(groupw):
			t1=dd.loc[groupd[grp,target],['start','stop']].values.T
			t2=dw.loc[groupw[grp,target],[1,2]].values.T
			t1=np.sort(t1,axis=0)
			t2=np.sort(t2,axis=0)
			t2=t2[:,np.argsort(t2[0])]
			#Cumulative max
			t2m=np.r_[-1,np.maximum.accumulate(t2[1])]
			if whitelist_mode=='within':
				#Find regions in whitelist: start and stop inside any whitelist
				t3=(t2m[np.searchsorted(t2[0],t1.ravel(),side='right')]>=t1.ravel()).reshape(2,-1).all(axis=0)
			elif whitelist_mode=='intersect':
				#Find regions in whitelist: intersect with any whitelist
				t3=t2m[np.searchsorted(t2[0],t1[1])]>t1[0]
			else:
				raise ValueError(f'Unknown whitelist_mode: {whitelist_mode}')
			ans+=list(groupd[grp,target][t3])
				
		if len(ans)==0:
			raise RuntimeError('No potential regulation remains after whitelist selection.')
		logging.info('{}/{} potential regulations remain after whitelist selection.'.format(len(ans),len(dd)))
		dd=dd.loc[ans]

	#Chain linking
	#TF-TFBS links: processing multi-TF motifs such as in GATA1,GATA2,GATA3 format
	links=list(zip(db['TF'].values,db['loc'].values,db['score'].values))
	links=list(itertools.chain.from_iterable([[(y,x[1],x[2]) for y in x[0].split(',')] for x in links]))
	links=[np.array(x) for x in zip(*links)]
	#TFBS-target gene links by distance
	links=[links,[dd['region'].values,dd['target'].values,_linking_score(None,None,dd['dist'].values,mode=mode&4) if mode&4 else np.zeros(len(dd))]]
	#Compute regulator mask matrix using links
	t1=[groupby(x[1]) for x in links]
	#{target:[(source,strength),...]}
	links=[{y:list(zip(x[0][0][x[1][y]],x[0][2][x[1][y]])) for y in x[1]} for x in zip(links,t1)]
	while len(links)>1:
		#Connect link steps sequentially
		t2=links.pop()
		t1={x:list(itertools.chain.from_iterable([[[z[0],z[1]+y[1]]for z in links[-1][y[0]]] for y in t2[x] if y[0] in links[-1]])) for x in t2}
		t1={x[0]:[np.array(z) for z in zip(*x[1])] for x in t1.items() if len(x[1])>0}
		if len(t1)==0:
			raise RuntimeError('No link found.')
		t2={x:groupby(y[0]) for x,y in t1.items()}
		links[-1]={x:[[y,combine(t1[x][1][t2[x][y]])] for y in t2[x]] for x in t2}
		assert np.isfinite(np.concatenate([[y[1] for y in x] for x in links[-1].values()])).all()
	links=links[0]
	#Remove self links
	links={x[0]:list(filter(lambda y:y[0]!=x[0],x[1])) for x in links.items()}		# pylint: disable=W0640
	links=[[[x[0]]+y for y in x[1]] for x in links.items()]
	
	#Construct matrix
	links=list(itertools.chain.from_iterable(links))
	links=[list(x) for x in zip(*links)]
	links=[links[1],links[0]]+links[2:]
	#TF and target names
	namelink=[]
	namelink.append(sorted(list(set(links[0]))))
	namelink.append(sorted(list(set(links[0])|set(links[1]))))
	nlink=[len(x) for x in namelink]
	dictlink=[dict(zip(x,range(y))) for x,y in zip(namelink,nlink)]
	links=[[dictlink[0][x] for x in links[0]],[dictlink[1][x] for x in links[1]],links[2]]
	ans=-np.ones(nlink,dtype=float)*np.inf
	ans[links[0],links[1]]=links[2]
	assert np.isfinite(ans).any() and not np.isnan(ans).any()
	#Slim links
	t1=[(ans>-np.inf).any(axis=1-x) for x in [0,1]]
	ans=ans[t1[0]][:,t1[1]]
	namelink=[np.array(x)[y] for x,y in zip(namelink,t1)]
	ans=pd.DataFrame(ans,index=namelink[0],columns=namelink[1])
	#Adds TFs to potential targets if missing
	t1=set(ans.columns)
	t1=list(filter(lambda x:x not in t1,ans.index))
	if len(t1)>0:
		ans=pd.DataFrame(np.concatenate([ans.values,np.ones((ans.shape[0],len(t1)),dtype=ans.values.dtype)*(-np.inf)],axis=1),index=ans.index,columns=list(ans.columns)+t1)
	#Output	
	logging.info(f'Writing file {fo_linking}')
	ans.to_csv(fo_linking,index=True,header=True,sep='\t')

def _binlinking_func(data:NDArray,n:int,inf:str='never')->NDArray:
	"""
	Actual function to convert regulator-target link scores to binary.

	Parameters
	----------
	data:	numpy.ndarray(ndim=1,dtype=float)
		Float vector to convert to binary vector
	n:		int
		Number of Trues to obtain in binary vector
	inf:
		Whether to select links with -inf score. Accepts: never.
	
	Returns
	-------
	numpy.ndarray(shape=data.shape,dtype=bool)
		Converted binary vector
	"""
	import numpy as np
	assert n>0 or n==-1
	assert data.ndim==1 and data.size>0
	if n==-1 or n>=data.size:
		assert inf=='never'
		return data>-np.inf
	t1=np.argpartition(data,-n)
	if data[t1[-n]]==-np.inf:
		if inf=='never':
			t1=np.nonzero(data>-np.inf)[0]
		else:
			raise NotImplementedError(f'Unknown inf {inf}.')
	else:
		t1=t1[-n:]
	assert len(t1)<=n
	ans=np.zeros(data.size,dtype=bool)
	ans[t1]=True
	return ans

def binlinking(fi_linking:str,fo_binlinking:str,n:int,axis:Optional[int]=1,selfreg:str='error',inf:str='never')->None:
	"""
	Converting regulator-target link score matrix to binary.

	Chooses the top regulator-target links, separately for each target gene by default.

	Parameters
	----------
	fi_linking:
		Path of input matrix file of TF to potential target gene link scores
	fo_binlinking:
		Path of output matrix file of TF to potential target gene links
	n:
		Number of regulator-target links. `n` strongest links (with highest scores) are selected along axis `axis`. If greater than the maximum links available, all links will be selected subject to `inf` parameter constraint. Value -1 selects all non-inf links.
	axis:
		Axis to choose the top links. If None, n links in total are selected among all regulator-target links. Defaults to 1, indicating n strongest links/regulators are selected for each target.
	selfreg:
		How to handle self regulation. Accepts:

		* error: Raise ValueError if seen in `fi_linking`

	inf:
		Whether to select links with -inf score. Accepts: never.
	"""
	import numpy as np
	import pandas as pd
	import logging
	#Loading input
	logging.info(f'Reading file {fi_linking}')
	dl=pd.read_csv(fi_linking,header=0,index_col=0,sep='\t')
	if selfreg=='error':
		tdict=dict(zip(dl.columns,range(len(dl.columns))))
		t1=dl.values[np.arange(len(dl)),[tdict[x] for x in dl.index]]
		assert t1.ndim==1
		if (t1!=-np.inf).any():
			raise ValueError(f'Input file {fi_linking} contains self-regulation.')
	else:
		raise ValueError(f'Unknown option {selfreg} for parameter selfreg')

	#Selecting top regulators
	if axis is None:
		links=_binlinking_func(dl.values.ravel(),n,inf=inf).reshape(*dl.shape)
	else:
		t1=np.swapaxes(dl.values,0,axis)
		t2=t1.shape
		t1=t1.reshape(t2[0],-1)
		links=np.array([_binlinking_func(x,n,inf=inf) for x in t1])
		links=np.swapaxes(links.reshape(*t2),0,axis)
	assert links.shape==dl.shape
	links=pd.DataFrame(links,index=dl.index,columns=dl.columns)
	if selfreg=='error':
		tdict=dict(zip(links.columns,range(len(links.columns))))
		t1=links.values[np.arange(len(links)),[tdict[x] for x in links.index]]
		assert t1.ndim==1
		assert not t1.any()

	#Reduce matrix size
	namelink=[links.index[links.values.any(axis=1)],links.columns[links.values.any(axis=0)]]
	if len(namelink[0])==0:
		raise RuntimeError('No link found.')
	namelink[1]=sorted(list(set(namelink[0])|set(namelink[1])))
	links=links.loc[namelink[0]][namelink[1]]
	logging.info(f'Writing file {fo_binlinking}')
	links.to_csv(fo_binlinking,header=True,index=True,sep='\t')
	


















































#
