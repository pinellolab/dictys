#!/usr/bin/python3
# Lingfei Wang, 2018-2021. All rights reserved.
"""Chromatin accessibility analyses
"""

from typing import Union

################################################################
# Peak calling
################################################################

def macs2(fi_names:str,fi_bam:str,fo_bam:str,fo_bai:str,fo_bed:str,genome_size:str,qcut:float=0.05,nth:int=1,nmax:int=0)->None:
	"""
	Peak calling using macs2, based on bam files for each cell in a given folder.
	
	Parameters
	------------
	fi_names:
		Path of input file containing one sample/cell name per line for macs2 peak calling
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
	from .utils import shell
	if qcut<=0 or qcut>=1:
		raise ValueError('qcut must be between 0 and 1.')
	if not isfile(fi_names):
		raise FileNotFoundError(fi_names)
	if not isdir(fi_bam):
		raise FileNotFoundError(fi_bam)
	
	scriptpath=pjoin(abspath(dirname(__file__)),'scripts',basename(__file__)[:-3],'chromatin_macs2.sh')
	fi_names,fi_bam,fo_bam,fo_bai,fo_bed=[abspath(x) for x in [fi_names,fi_bam,fo_bam,fo_bai,fo_bed]]
	
	#Load sample names
	with open(fi_names,'r') as f:
		names=f.readlines()
	names=[x.strip() for x in names]
	names=list(filter(lambda x:len(x)>0,names))
	if len(names)==0:
		raise ValueError('No sample name found in '+fi_names)
	namestxt=linesep.join([x+'.bam' for x in names])+linesep
	
	#Run script for macs2
	cmd = scriptpath+f" cellnames.txt {fi_bam} {fo_bam} {fo_bai} {fo_bed} {genome_size} {qcut} {nth}"
	d2 = shell.cmdfile(cmd,[],infiles={'cellnames.txt': namestxt},quiet=False,cd=True)
	if d2 is None or len(d2)>0 or not all(isfile(x) for x in [fo_bam,fo_bai,fo_bed]):
		raise RuntimeError('Macs2 failed.')

	if nmax>0:
		#Reduce size of peak bed file
		d3=pd.read_csv(fo_bed,header=None,index_col=None,sep='\t')
		if len(d3)>nmax:
			d4=np.partition(d3[8].values,-nmax-1)[-nmax-1]
			d3=d3[d3[8]>d4]
			d3.to_csv(fo_bed,header=False,index=False,sep='\t')

################################################################
# TF footprinting
################################################################

def wellington(fi_bam:str,fi_bai:str,fi_bed:str,fo_bed:str,fi_blacklist:Union[str,None]=None,cut:float=10,nth:int=1,npeakmax:int=1000000000)->None:
	"""
	TF Footprinting with wellington.
	
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
	npeakmax:	
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
	
	scriptpath=pjoin(abspath(dirname(__file__)),'scripts',basename(__file__)[:-3],'chromatin_wellington.sh')
	fi_bam,fi_bai,fi_bed,fo_bed=[abspath(x) for x in [fi_bam,fi_bai,fi_bed,fo_bed]]
	cmd = scriptpath+f" {fi_bam} {fi_bai} {fi_bed} {fo_bed} {cut} {nth} {npeakmax} {fi_blacklist}"
	d2 = shell.cmdfile(cmd,[],quiet=False,cd=True)
	if d2 is None or len(d2)>0 or not all(isfile(x) for x in [fo_bed]):
		raise RuntimeError('Wellington failed.')

################################################################
# Motif scan
################################################################

def _motif_postproc(dret,fo_bed,fo_wellington,fo_homer):
	"""Postprocess motif discovery results from HOMER.
	d:		Original dataset object
	dret:	Return tuple of HOMER call
	Return:
	Processed and copied new data object with motif results"""
	import numpy as np
	import pandas as pd
	from io import StringIO
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

	namet=d.dim['nt']
	#Set na as or function
	t1=(pd.isna(dw)|pd.isna(dh)).values
	dw.fillna(0,inplace=True)
	dh.fillna(0,inplace=True)
	#Extract dimensions
	namep=np.array([str(x) for x in dw.index])
	namem=np.array(['_'.join([x.split('_')[0]]+x.split('.')[-2:]) for x in dw.columns])
	#Set data values
	dw,dh=[x.values.astype(ftype0) for x in [dw,dh]]
	dw[t1]=0
	dh[t1]=0
	#Get motif to gene map
	dmt=[x.split('_')[0] for x in namem]
	t1=dict(zip([x.upper() for x in namet],range(len(namet))))
	dmt=np.array([t1[x.upper()] if x.upper() in t1 else -1 for x in dmt],dtype='i8')
	#Remove empty motifs
	t1=dmt>=0
	dw,dh=[x[:,t1] for x in [dw,dh]]
	namem,dmt=[x[t1] for x in [namem,dmt]]
	assert dw.shape==(len(namep),len(namem)) and dh.shape==(len(namep),len(namem))
	assert np.isfinite(dw).all() and np.isfinite(dh).all()
	assert (dw>=0).all() and (dh>=0).all()
	#Output
	dmotif.to_csv(fo_bed,header=False,index=False)
	dw=pd.DataFrame(dw,index=namep,header=namem)
	dw.to_csv(fo_wellington,header=True,index=True)
	dh=pd.DataFrame(dh,index=namep,header=namem)
	dh.to_csv(fo_homer,header=True,index=True)

def homer(fi_bed:str,fi_motif:str,fo_bed:str,fo_wellington:str,fo_homer:str,genome:str,nth:int=1)->None:
	"""
	Motif scan with homer.
	
	Parameters
	------------
	fi_bed:	
		Path of input bed file of regions
	fi_motif:
		Path of input motif PWM file in homer format
	fo_bed:	
		Path of output bed file of detected motifs
	fo_wellington:
		Path of output tsv file of wellington scores in shape (region,motif)
	fo_homer:
		Path of output tsv file of homer scores in shape (region,motif)
	genome:	
		Reference genome for homer (e.g. hg19, mm10)
	nth:	
		Number of threads

	"""
	from os.path import dirname, basename, abspath, isfile
	from os.path import join as pjoin
	from .utils import shell
	for xi in [fi_bed,fi_motif]:
		if not isfile(xi):
			raise FileNotFoundError(xi)
	
	scriptpath=pjoin(abspath(dirname(__file__)),'scripts',basename(__file__)[:-3],'chromatin_homer.sh')
	spath=pjoin(abspath(dirname(__file__)),'scripts',basename(__file__)[:-3],'chromatin_homer.py')
	fi_bed,fi_motif,fo_bed,fo_wellington,fo_homer=[abspath(x) for x in [fi_bed,fi_motif,fo_bed,fo_wellington,fo_homer]]
	cmd = scriptpath+f" {fi_bed} {fi_motif} {genome} {spath} {nth}"
	d2 = shell.cmdfile(cmd,
		['19-w.tsv','19-h.tsv','16-long.bed'],
		quiet=False,cd=True,sizelimit=None)
	return _motif_postproc(d2,fo_bed,fo_wellington,fo_homer)
	































































#
