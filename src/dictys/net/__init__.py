#!/usr/bin/python3
# Lingfei Wang, 2022. All rights reserved.

"""
Module for network class and stats.

"""

from __future__ import annotations
from typing import Set,Optional,Tuple,Callable,Union
import pandas as pd

__all__=['stat']

import dictys
from dictys.utils.importing import matplotlib
from . import *

class network:
	"""
	Class for context-specific network and parent class for dynamic networks

	Required variables
	----------------------------------
	?n:		int
		Number of cells/states/nodes (for $1=c/s/n respectively)
	?name:	numpy.ndarray(shape=($1n,),dtype=str)
		Cell/state/node names
	?dict:	dict
		Map from cell/node/state name to id
	nns:	numpy.ndarray(shape=(2,),dtype=int)
		Number of source & target nodes
	nids:	[numpy.ndarray(shape=(nns[0],),dtype=int),numpy.ndarray(shape=(nns[1],),dtype=int)]
		Lists of indices of nodes that are regarded as source & target nodes in `prop`

	Optional variables
	----------------------------------
	prop:	{type:{name:value}}
		Properties of cell/state/node/edge.
		type:	str
			Indicating the owner (or owner pair) of property. Takes one of the values
			c:		Cell (with shape factor cn)
			s:		State (with shape factor sn)
			n:		Network node (with shape factor nn)
			e:		Network edge (with shape factor (nns[0],nns[1]))
			Other:	Tuple of multiple owners. Each character represents a owner. Shape factors are concatenated.
					For example, ss means state-state properties with shape factor (sn,sn).
		value:	numpy.ndarray
			Value of each property. The trailing axes in its shape must exactly match the shape factor from `type`.
	traj:	dictys.traj.trajectory
		Trajectory object. Only needed for dynamic network
	point:	{x:dictys.traj.point for x in ['c','s']}
		Point objects for cells and states on the `traj` object. Only needed for dynamic network.
	"""
	_name_='network'
	_version_=[0,1,0]
	_shape_={'e':lambda x:tuple(x.nns)}
	def _get_prop_shape_(self,propname:str)->Tuple:
		"""
		Get shape of property from _shape_, use default (prefix with n) if not available.
		"""
		if propname.endswith('prop'):
			propname=propname[:-4]
		return self._shape_[propname](self) if propname in self._shape_ else (getattr(self,propname+'n'),)
	#Constructor
	def __init__(self,**ka):
		"""
		Network class constructor that should not be called directly. Please use class functions create_dynamic or create_group for creating networks from data, or from_file from file.

		Parameters
		-------------
		ka:			Keyword arguments saved to class
		"""
		import numpy as np
		from collections import defaultdict
		#Storing values
		params_allowed=set('cname,sname,nname,nids,traj,prop,point'.split(','))
		self.prop=defaultdict(dict)
		for xi in ka:
			if xi.endswith('prop') and len(xi)>4:
				self.prop[xi[:-4]]=ka[xi]
			else:
				if xi not in params_allowed:
					raise ValueError(f'Unknown parameter {xi}.')
				setattr(self,xi,ka[xi])
		#Post processing
		for xi in [['c','cell'],['s','state'],['n','node']]:
			vnames=[xi[0]+'n',xi[0]+'name',xi[0]+'dict']
			#Fill count/name with another if needed
			setattr(self,vnames[0],len(getattr(self,vnames[1])))
			#Prepare dictionary
			setattr(self,vnames[2],dict(zip(getattr(self,vnames[1]),range(getattr(self,vnames[0])))))
		#Fill count if needed
		self.nns=np.array([len(x) for x in self.nids])
		if not hasattr(self,'traj'):
			self.traj=None
			self.point={}
		#Validation
		self.check()
	def check(self)->None:
		"""
		Validate self.
		"""
		import numpy as np
		from collections import Counter
		from dictys.traj import trajectory
		for xi in [[self.cn,self.cname,self.cdict,'cell','c'],[self.sn,self.sname,self.sdict,'state','s'],[self.nn,self.nname,self.ndict,'gene','n']]:
			if xi[0]==0:
				raise ValueError(f'No {xi[3]} included.')
			if xi[1].shape!=(xi[0],):
				raise ValueError('Expecting {} {} names, actual count: {}'.format(xi[0],xi[3],xi[1].shape[0]))
			if len(xi[2])!=xi[0]:
				t1=[x[0] for x in dict(Counter(xi[1])).items() if x[1]>1][:3]
				assert len(t1)>0
				raise ValueError('Found duplicate {} names. First three: {}'.format(xi[3],', '.join(t1)))
			assert all(x in xi[2] for x in xi[1])
		for xi in self.prop:
			t1=tuple(np.concatenate([self._get_prop_shape_(x) for x in xi]))
			t2=list(filter(lambda x:self.prop[xi][x].shape[-len(t1):]!=t1,self.prop[xi]))		# pylint: disable=W0640
			if len(t2)>0:
				raise ValueError('Properties of type {} have incorrect size. Expected trailing dimensions: {}. Actual dimensions: '.format(xi,t1)+', '.join(['{} ({})'.format(tuple(self.prop[xi][x].shape),x) for x in t2]))
		assert len(self.nids)==2 and all(len(np.unique(x))==len(x) for x in self.nids)
		assert all((x>=0).all() and (x<self.nn).all() for x in self.nids)
		assert (self.nns==[len(x) for x in self.nids]).all() and (self.nns>0).all()
		assert self.traj is None or isinstance(self.traj,trajectory)
		assert (self.traj is None and len(self.point)==0) or (self.traj is not None and frozenset(self.point)==frozenset({'c','s'}))
		assert all(len(self.point[x])==getattr(self,x+'n') for x in self.point)
	#Rename dimensions
	def rename(self,dim:str,vals:Union[list[str],dict[str,str]])->None:
		"""
		Rename any dimension
		dim:
			Dimension to rename: c/s/n for cells/states/nodes
		vals:
			New names. If as a dict, names not included will not be renamed.
		"""
		import numpy as np
		if dim not in 'csn':
			raise ValueError("Can only renames cells/states/nodes with dim='c', 's', or 'n'")
		names=getattr(self,dim+'name')
			
		if isinstance(vals,list):
			if len(vals)!=len(names):
				raise TypeError('Unequal size between new ({}) and old ({}) names when using vals as a list'.format(len(vals),len(names)))
			vals=dict(zip(names,vals))
		elif not isinstance(vals,dict):
			raise TypeError('vals should be list or dict.')
		names=[vals[x] if x in vals else x for x in names]
		if len(names)!=len(set(names)):
			raise ValueError('Duplicate names found after renaming.')
		setattr(self,dim+'name',np.array(names))
		setattr(self,dim+'dict',dict(zip(names,range(len(names)))))
	#I/O
	@classmethod
	def _from_file_0_1_0(cls,path:str)->network:
		"""
		Load class instance from file version 0.1.0.

		Parameters
		-------------
		path:	str
			Path of file to load from.

		Returns
		---------
		net:	dictys.net.network
			network class instance
		"""
		from collections import defaultdict
		from dictys.traj import trajectory,point
		import numpy as np
		import h5py
		import logging
		params={}
		logging.info(f'Reading file {path}.')
		with h5py.File(path,'r') as f:
			for xi in filter(lambda x:x not in {'prop','traj'},f):
				params[xi]=np.array(f[xi])
				if params[xi].dtype.char=='S':
					params[xi]=params[xi].astype(str)
			params['prop']=defaultdict(dict)
			for xi in f['prop']:
				for xj in f['prop'][xi]:
					params['prop'][xi][xj]=np.array(f['prop'][xi][xj])
					if params['prop'][xi][xj].dtype.char=='S':
						params['prop'][xi][xj]=params['prop'][xi][xj].astype(str)
			if 'traj' in f:
				params['traj']=trajectory.from_fileobj(f['traj'])
			if 'point' in f:
				assert 'traj' in params
				params['point']={}
				for xi in f['point']:
					params['point'][xi]=point.from_fileobj(params['traj'],f['point'][xi])
		params['nids']=[params['nids1'],params['nids2']]
		del params['nids1'],params['nids2']
		if '_version_' in params:
			del params['_version_']
		return cls(**params)
	@classmethod
	def from_file(cls,path:str)->network:
		"""
		Load class instance from file.

		Parameters
		-------------
		path:	str
			Path of file to load from.

		Returns
		---------
		net:	dictys.net.network
			network class instance
		"""
		import numpy as np
		import h5py
		import logging
		logging.info(f'Reading file {path}.')
		with h5py.File(path,'r') as f:
			if '_version_' not in f:
				v=(0,1,0)
			else:
				v=tuple(np.array(f['_version_']))
		if v==(0,1,0):
			return cls._from_file_0_1_0(path)
		raise TypeError('Unknown network version {}'.format('.'.join([str(x) for x in v])))
	def to_file(self,path:str,compression:str="gzip",links:bool=True,**ka)->None:
		"""
		Save class instance to file.
		
		Parameters
		-------------
		path:		str
			Path of file to save to.
		compression: str
			Compression used.
		links:		bool
			Whether to use HDF5 links to save space.
		ka:	dict
			Keyword arguments passed to h5py.File.create_dataset.
		"""
		import logging
		import numpy as np
		import h5py
		params={x:getattr(self,x) for x in 'cname,sname,nname,nids'.split(',')}
		params['nids1']=params['nids'][0]
		params['nids2']=params['nids'][1]
		del params['nids']
		logging.info(f'Writing file {path}.')
		past=[]
		with h5py.File(path,'w') as f:
			f.create_dataset('_version_',len(self._version_),dtype='i')[:]=np.array(self._version_)
			for xi in params:
				if links:
					t1=list(filter(lambda x:x[1] is params[xi],past))		# pylint: disable=W0640
					if len(t1)>0:
						f[xi]=t1[0][2]
						continue
				p=dict(ka)
				p['compression']=compression
				p['data']=params[xi]
				if params[xi].dtype.char=='U':
					p['data']=p['data'].astype('S')
				entry=f.create_dataset(xi,**p)
				past.append((xi,params[xi],entry))
			f.create_group('prop')
			for xi in self.prop:
				f['prop'].create_group(xi)
				for xj in self.prop[xi]:
					data=self.prop[xi][xj]
					if links:
						t1=list(filter(lambda x:x[1] is data,past))		# pylint: disable=W0640
						if len(t1)>0:
							f['prop'][xi][xj]=t1[0][2]
							continue
					if data.dtype.char=='U':
						data=data.astype('S')
					entry=f['prop'][xi].create_dataset(xj,data=data,compression=compression,**ka)
					past.append((xi+'.'+xj,data,entry))
			if self.traj is not None:
				t=f.create_group('traj')
				self.traj.to_fileobj(t)
				t=f.create_group('point')
				for xi in self.point:
					self.point[xi].to_fileobj(t.create_group(xi))
	@staticmethod
	def _load_subsets(fi_subsets:str,path_cname:Optional[str]=None)->Optional[Tuple[list[str],list[list[str]]]]:
		"""
		Load cell subsets from folder.
		
		Parameters
		----------
		fi_subsets:
			Path of input txt file of cell subset names
		path_cname:
			Path of folder containing subfolders of cell subsets. Defaults to basedir of fi_subsets.
			
		Returns
		-------
		sname:
			List of cell subset names
		cnames:
			List of list of cell names in each cell subset		
		"""
		from os.path import join as pjoin
		from os.path import dirname
		from dictys.utils.file import read_txt
		if path_cname is None:
			path_cname=pjoin(dirname(fi_subsets),'subsets')
		#Cells & states/subsets
		try:
			sname=read_txt(fi_subsets,unique=True)
		except FileNotFoundError:
			sname=[]
		n=len(sname)
		if n==0:
			return None
		cnames=[]
		for xi in range(n):
			cnames.append(read_txt(pjoin(path_cname,sname[xi],'names_rna.txt')))
		return (sname,cnames)
	@classmethod
	def from_folders(cls,path_data:str,path_work:str,fi_subsets:str,dynamic:bool=False,nettype:str='n',optional:Set[str]={'readcount'},fi_c:Optional[str]=None)->network:
		"""
		Create class instance from pipeline working folder.

		Parameters
		-------------
		path_data:
			Path of data folder to load from
		path_work:
			Path of working folder to load from
		fi_subsets:
			Path of input txt file for cell subset names
		dynamic:
			Whether to load a dynamic network instead of a set of static networks
		nettype:
			Type of network. Accepts:

			* '':	Unnormalized direct & steady-state network.

			* 'n':	Normalized direct & steady-state networks.

		optional:
			Optional data to include. Accepts:

			* readcount:	RNA read count for each cell.

		fi_c:
			Path of input tsv file for extra property columns for each cell

		Returns
		---------
		dictys.net.network
			network class instance
		"""
		import logging
		import itertools
		from os.path import join as pjoin
		from dictys.traj import trajectory,point
		from dictys.utils.numpy import dtype_min
		import numpy as np
		nettypes=[nettype,'i'+nettype]
		ntype=len(nettypes)

		optional_allowed={'readcount'}
		if len(optional-optional_allowed)>0:
			raise ValueError('Unknown optional data types: '+','.join(optional-optional_allowed))

		params={}
		#Cells & states/subsets
		ans=cls._load_subsets(fi_subsets,path_cname=path_work)
		if ans is None:
			raise RuntimeError('Could not find cell subsets.')
		sname,cnames=ans
		n=len(sname)
		cname=list(itertools.chain.from_iterable(cnames))
		if not dynamic:
			assert len(cname)==len(set(cname))
		else:
			cname=sorted(list(set(cname)))
		cname_dict=dict(zip(cname,range(len(cname))))

		#State-cell properties: cell membership in subset
		msc=np.zeros((n,len(cname)),dtype=bool)
		for xi in range(n):
			msc[xi,[cname_dict[x] for x in cnames[xi]]]=True
		msc=msc.astype(float)
		scprop={'w':msc}

		#Network
		nets=[[] for _ in range(ntype)]
		success=np.ones((ntype,n),dtype=bool)
		for xi in range(n):
			for xj in range(ntype):
				nettype=nettypes[xj]
				try:
					fi=pjoin(path_work,sname[xi],f'net_{nettype}weight.tsv.gz')
					logging.info(f'Reading file {fi}')
					nets[xj].append(pd.read_csv(fi,header=0,index_col=0,sep='\t'))
					if nets[xj][xi].shape[0]==0:
						success[xj,xi]=False
				except FileNotFoundError as e:
					logging.warning('Skipping cell subset {} due to error: {}'.format(sname[xi],repr(e)))
					success[xj,xi]=False
					nets[xj].append(pd.DataFrame([],index=[],columns=[]))
		if not success.any():
			raise RuntimeError('No network found.')

		#Network nodes: source & target
		nidss=[[x.index,x.columns] for x in itertools.chain.from_iterable(nets)]
		assert all((x[1][:len(x[0])]==x[0]).all() for x in nidss)
		nids=[sorted(list(set(itertools.chain.from_iterable(x)))) for x in zip(*nidss)]
		niddicts=[dict(zip(x,range(len(x)))) for x in nids]
		nname=nids[1]
		ndict=dict(zip(nname,range(len(nname))))
		nids=[np.array([ndict[y] for y in x]) for x in nids]
		cname=np.array(cname)
		sname=np.array(sname)
		nname=np.array(nname)
		#Counts
		ns=len(sname)
		nns=np.array([len(x) for x in nids])

		#Cell properties
		cprop={}
		if fi_c is not None:
			#External file
			logging.info(f'Reading file {fi_c}')
			t1=pd.read_csv(fi_c,header=0,index_col=0,sep='\t')
			if len(set(cname)-set(t1.index))>0:
				raise ValueError(f'Not all known cells included in {fi_c}.')
			t1=t1.loc[cname]
			for xi in t1.columns:
				cprop[xi]=t1[xi].values
		if dynamic:
			#UMAP coordindates
			fi=pjoin(path_data,'coord_rna.tsv.gz')
			logging.info(f'Reading file {fi}')
			t1=pd.read_csv(fi,header=0,index_col=0,sep='\t')
			cprop['coord']=t1.loc[cname].values.T
			#Cell type/context
			ans=cls._load_subsets(pjoin(path_data,'subsets.txt'))
			if ans is None:
				logging.warning('Could not find cell contexts.')
			else:
				t1=dict(itertools.chain.from_iterable([[(y,x[0]) for y in x[1]] for x in zip(*ans)]))
				cprop['type']=np.array([t1[x] for x in cname])

		#Node-cell properties
		fi=pjoin(path_data,'expression.tsv.gz')
		logging.info(f'Reading file {fi}')
		readcount=pd.read_csv(fi,header=0,index_col=0,sep='\t')
		readcount=readcount.loc[nname][cname]
		ncprop={}
		if 'readcount' in optional:
			#Read counts
			ncprop['readcount']=readcount.values.astype(dtype_min(readcount.values.astype(int)))

		#Node-state properties
		nsprop={}
		#Node (pseudobulk) CPM
		t1=readcount.values@msc.T
		t1=1E6*t1/t1.sum(axis=0)
		assert (t1>=0).all() and (t1<=1E6).all()
		nsprop['cpm']=t1
		
		#Edge-state properties
		esprop={}
		for xi in nettypes:
			#Edge weight
			esprop['w_'+xi]=np.zeros(np.r_[nns,ns],dtype=float)
			#Edge mask
			esprop['mask_'+xi]=np.zeros(np.r_[nns,ns],dtype=bool)
		for xi,xj in np.array(np.nonzero(success)).T:
			vsuf='_'+nettypes[xi]
			t1=[[niddicts[y][x] for x in nidss[xi*n+xj][y]] for y in range(2)]
			esprop['w'+vsuf][np.ix_(t1[0],t1[1],[xj])]=nets[xi][xj].values.reshape(*nets[xi][xj].shape,1)
			esprop['mask'+vsuf][np.ix_(t1[0],t1[1],[xj])]=True
		esprop['w']=esprop['w_'+nettypes[0]]
		esprop['mask']=esprop['mask_'+nettypes[0]]

		if dynamic:
			#State-state properties
			ssprop={}
			#Whether every state point pair on the trajectory is connected
			fi=pjoin(path_work,'subset_edges.tsv.gz')
			logging.info(f'Reading file {fi}')
			t1=pd.read_csv(fi,header=0,index_col=0,sep='\t')
			ssprop['traj-neighbor']=t1.values

			#Trajectory
			traj=trajectory.from_file(pjoin(path_data,'traj_node.h5'))
			#Points
			points={}
			points['c']=point.from_file(traj,pjoin(path_data,'traj_cell_rna.h5'))
			points['s']=point.from_file(traj,pjoin(path_work,'subset_locs.h5'))

		#Construction
		params={'cname':cname,'sname':sname,'nname':nname,'nids':nids,'cprop':cprop,'scprop':scprop,'ncprop':ncprop,'esprop':esprop,'nsprop':nsprop}
		if dynamic:
			params.update({'ssprop':ssprop,'traj':traj,'point':points})

		success=success.any(axis=0)
		if success.all():
			return cls(**params)
		if not success.any():
			raise RuntimeError('All subsets failed.')
		
		#Remove failed states
		logging.warning('Removing {}/{} failed subsets.'.format((~success).sum(),len(success)))
		success=np.nonzero(success)[0]
		params['sname']=params['sname'][success]
		if 'point' in params:
			params['point']['s']=params['point']['s'][success]
		for xi in filter(lambda x:x.endswith('prop'),params):
			t1=np.nonzero([x=='s' for x in xi[:-4]])[0]
			t2=np.cumsum(np.r_[0,[2 if x=='e' else 1 for x in xi[:-4]]])
			for xj in t1:
				params[xi]={x:y.swapaxes(t2[xj],0)[success].swapaxes(t2[xj],0) for x,y in params[xi].items()}

		return cls(**params)
	def export(self,output_folder:str,sparsities:Union[None,float,list[Union[None,float]]]=None)->None:
		"""
		Export context specific networks to tsv files. Each context specific network is exported to an accordingly named file. If certain contexts do not have exported networks, the sparsity level is too large or too small.

		Parameters
		----------
		output_folder:
			Path of output folder. Must be absent.
		sparsities:
			Network sparsity to use to convert to binary networks. Accepts:

			* None: Output continuous network into folder "Full"

			* float: Output binary network with the corresponding sparsity into folder "sparsity=..."

			* List of None and/or float: Output networks of different types or sparsities into their corresponding folders

		"""
		from os.path import join as pjoin
		from os import makedirs
		import numpy as np
		from dictys.net import stat
		if sparsities is None or isinstance(sparsities,float):
			sparsities=[sparsities]
		
		makedirs(output_folder)
		stat1_lcpm=stat.lcpm(self,cut=0)
		stat1_net0=stat.net(self)
		pts=np.arange(self.sn)
		#Output CPM
		fo=pjoin(output_folder,'cpm.tsv.gz')
		d1=stat1_lcpm.compute(pts)
		t1=(d1>0).any(axis=1)
		d1=pd.DataFrame(d1[t1],index=stat1_lcpm.names,columns=self.sname)
		d1.to_csv(fo,header=True,index=True,sep='\t')
		#Output network
		for xi in sparsities:
			output_folder1=pjoin(output_folder,'Full' if xi is None else f'sparsity={xi}')
			makedirs(output_folder1)
			stat1_net=stat.fbinarize(stat1_net0,sparsity=xi) if xi is not None else stat1_net0
			d1=stat1_net.compute(pts)
			for xj in range(self.sn):
				fo=pjoin(output_folder1,self.sname[xj]+'.tsv.gz')
				d2=d1[:,:,xj]
				#Remove empty regulators & targets
				t1=[(d2!=0).any(axis=1-x) for x in range(2)]
				if xi is not None and (any(not x.any() for x in t1) or d2[t1[0]][:,t1[1]].all()):
					#Skip sparsity for cell type because too large/small (all/no edges are positive)
					continue
				d2=pd.DataFrame(d2[t1[0]][:,t1[1]],index=stat1_net.names[0][t1[0]],columns=stat1_net.names[1][t1[1]])
				#Output to tsv.gz file
				d2.to_csv(fo,header=True,index=True,sep='\t')

class dynamic_network(network):
	"""
	Class for dynamic networks with extra functions
	"""
	_name_='dynamic_network'
	def linspace(self,start:int,stop:int,num:int,dist:float,check:bool=True)->Tuple[dictys.traj.point,Callable[[dictys.net.stat.base],dictys.net.stat.base]]:
		"""
		Returns evenly spaced points between two nodes for Gaussian kernel smoothing. Removes branches not on the path between these nodes.
		
		Parameters
		----------
		start:
			Starting node ID
		stop:
			Ending node ID
		num:
			Number of points to generate. Must be non-negative.
		dist:
			Gaussian kernel smoothing distance
		check:
			Whether to perform checks, e.g. starting and ending nodes are branching or terminal nodes.
		
		Returns
		-------
		pts:		dictys.traj.point
			Evenly spaced points with pts.p as its trajectory after removing unrelated branches
		fsmooth:	functools.partial
			Partial function that produces Gaussian kernel smoothened statistic on an original statistic as its parameter
		"""
		import numpy as np
		from functools import partial
		#Gaussian kernel smoothing parameters
		weightfunc=['conv',[dist],dict({'criterion_path':'loose'})]
		if check and ((self.traj.edges==start).sum()==2 or (self.traj.edges==stop).sum()==2):
			raise ValueError('Start or stop node is not a branching or terminal node.')
		
		#Compute node location & order based by averaging cell locations on traj
		cell=self.point['c']
		point=cell.nearest(w=self.prop['sc']['w'])
		#Find start/end nodes on points
		point0=self.point['s']
		t1=(dictys.traj.point.fromnodes(self.traj)[[start,stop]]-point0)==0
		t1=[np.nonzero(x)[0] for x in t1]
		assert all(len(x)==1 for x in t1)
		start2,stop2=[x[0] for x in t1]

		#Subtrajectory
		traj2=point.subtraj(edges=np.array(np.nonzero(np.triu(self.prop['ss']['traj-neighbor']))).T)
		pts=traj2.linspace(start2,stop2,num)
		weightfunc[2]['nodes_path']=[start2,stop2]
		fsmooth=partial(dictys.net.stat.fsmooth,pts=traj2,smoothen_func=weightfunc)
		return (pts,fsmooth)
	def compute_chars(self,start:int,stop:int,num:int=100,dist:float=1.5,mode:str='regulation',sparsity:float=0.01)->pd.DataFrame:
		"""
		Compute curve characteristics for one branch.

		Parameters
		----------
		start:
			Branch starting node ID 
		stop:
			Branch ending node ID 
		num:
			Number of points from starting to ending nodes to draw
		dist:
			Gaussian kernel smoothing distance/length scale between cells. Larger value means stronger smoothing.
		mode:
			Mode or measure to compute characteristics. Accepts:

			* 'regulation': based on target count

			* 'expression': based on CPM

		sparsity:
			The number of edges to regard as positive when binarizing network. Only relevant when mode=='regulation'

		Returns
		-------
		List of returns from dictys.plot.dynamic.draw_dynamic1 for each of the patterns drawn
		"""
		from dictys.net import stat
		from dictys.plot.dynamic import _compute_chars_
		pts,fsmooth=self.linspace(start,stop,num,dist)
		if mode=='regulation':
			#Log number of targets
			stat1_net=fsmooth(stat.net(self))
			stat1_netbin=stat.fbinarize(stat1_net,sparsity=sparsity)
			stat1_y=stat.flnneighbor(stat1_netbin)
		elif mode=='expression':
			stat1_y=fsmooth(stat.lcpm(self,cut=0))
		else:
			raise ValueError(f'Unknown mode {mode}.')
		#Pseudo time
		stat1_x=stat.pseudotime(self,pts)
		dy=pd.DataFrame(stat1_y.compute(pts),index=stat1_y.names[0])
		dx=pd.Series(stat1_x.compute(pts)[0])
		return _compute_chars_(dy,dx)
	def draw_discover(self,start:int,stop:int,num:int=100,dist:float=1.5,mode:str='regulation',sparsity:float=0.01,**ka)->list[Tuple[matplotlib.figure.Figure,list,dict[str,matplotlib.cm.ScalarMappable]]]:
		"""
		Draws TF discovery plots for one branch.

		Parameters
		----------
		start:
			Branch starting node ID 
		stop:
			Branch ending node ID 
		num:
			Number of points from starting to ending nodes to draw
		dist:
			Gaussian kernel smoothing distance/length scale between cells. Larger value means stronger smoothing.
		mode:
			Mode or measure to discover TFs. Accepts:

			* 'regulation': based on target count

			* 'weighted_regulation': based on weighted outdegree without the need to binarize network

			* 'expression': based on CPM

		sparsity:
			The number of edges to regard as positive when binarizing network. Function depends on mode:
			
			* For mode='regulation': effective

			* For mode='weighted_regulation': determines the overall scale of outdegree to be comparable with a binarized network of the specified sparsity
			
			* for mode='expression': no effect

		ka:
			Keyword arguments passed to dictys.plot.dynamic.fig_discover

		Returns
		-------
		List of returns from dictys.plot.dynamic.draw_dynamic1 for each of the patterns drawn
		"""
		from dictys.net import stat
		from dictys.plot.dynamic import fig_discover
		pts,fsmooth=self.linspace(start,stop,num,dist)
		if mode=='regulation':
			#Log number of targets
			stat1_net=fsmooth(stat.net(self))
			stat1_netbin=stat.fbinarize(stat1_net,sparsity=sparsity)
			stat1_y=stat.flnneighbor(stat1_netbin)
		elif mode=='weighted_regulation':
			#Log weighted outdegree
			stat1_net=fsmooth(stat.net(self))
			stat1_y=stat.flnneighbor(stat1_net,weighted_sparsity=sparsity)
		elif mode=='expression':
			stat1_y=fsmooth(stat.lcpm(self,cut=0))
		else:
			raise ValueError(f'Unknown mode {mode}.')
		#Pseudo time
		stat1_x=stat.pseudotime(self,pts)
		dy=pd.DataFrame(stat1_y.compute(pts),index=stat1_y.names[0])
		dx=pd.Series(stat1_x.compute(pts)[0])
		return fig_discover(dy,dx,**ka)
	def draw_regulation_heatmap(self,start:int,stop:int,regulations:list[Tuple[str,str]],num:int=100,dist:float=1.5,ax:Optional[matplotlib.axes.Axes]=None,cmap:Union[str,matplotlib.cm.ScalarMappable]='coolwarm',figsize:Tuple[float,float]=(2,0.22),vmax:Optional[float]=None)->Tuple[matplotlib.pyplot.Figure,matplotlib.axes.Axes,matplotlib.cm.ScalarMappable]:
		"""Draws pseudo-time dependent heatmap of individual regulation strengths.
	
		Parameters
		----------
		start:
			Branch starting node ID 
		stop:
			Branch ending node ID 
		regulations:
			List of regulations in (Regulator name, target name) format to draw strength
		num:
			Number of points from starting to ending nodes to draw
		dist:
			Gaussian kernel smoothing distance/length scale between cells. Larger value means stronger smoothing.
		ax:
			Axes to draw on.
		figsize:
			Figure size for each regulation as a row. Should remain unassigned when ax is assigned.
		cmap:
			Colormap in matplotlib string or matplotlib.cm.ScalarMappable
		vmax:
			Maximum value in colormap. Should remain unassigned when cmap is a matplotlib.cm.ScalarMappable instance.
		
		Returns
		-------
		fig:
			Heatmap figure
		ax:
			Heatmap axes
		cmap:
			Heatmap colormap
		"""
		from dictys.plot.dynamic import fig_regulation_heatmap
		return fig_regulation_heatmap(self,start,stop,regulations,num=num,dist=dist,ax=ax,cmap=cmap,figsize=figsize,vmax=vmax)
	def export(self,output_folder:str,start:int,stop:int,num:int,dist:float,sparsities:Union[None,float,list[Union[None,float]]]=None)->None:		# pylint: disable=W0221
		"""
		Export dynamic network to tsv files. Each time point network is exported to an accordingly named file. If certain contexts do not have exported networks, the sparsity level is too large or too small.

		Parameters
		----------
		output_folder:
			Path of output folder. Must be absent.
		start:
			Starting node ID to export network
		stop:
			Stopping node ID to export network
		num:
			Number of intermediate time points to export network
		dist:
			Gaussian kernel smoothing distance
		sparsities:
			Network sparsity to use to convert to binary networks. Accepts:

			* None: Output continuous network into folder "Full"

			* float: Output binary network with the corresponding sparsity into folder "sparsity=..."

			* List of None and/or float: Output networks of different types or sparsities into their corresponding folders

		"""
		from os.path import join as pjoin
		from os import makedirs
		import numpy as np
		from dictys.net import stat
		if sparsities is None or isinstance(sparsities,float):
			sparsities=[sparsities]
		makedirs(output_folder)
		
		pts,fsmooth=self.linspace(start,stop,num,dist)
		stat1_time=stat.pseudotime(self,pts)
		stat1_lcpm=fsmooth(stat.lcpm(self,cut=0))
		stat1_net0=fsmooth(stat.net(self))
		#Output pseudotime
		d1=stat1_time.compute(pts)[0]
		d1=pd.DataFrame(zip(np.arange(len(pts))+1,d1),columns=['State_index','Pseudo-time'])
		d1.to_csv(pjoin(output_folder,'pseudotime.tsv.gz'),header=True,index=False,sep='\t')
		#Output CPM
		fo=pjoin(output_folder,'cpm.tsv.gz')
		d1=stat1_lcpm.compute(pts)
		t1=(d1>0).any(axis=1)
		d1=pd.DataFrame(d1[t1],index=stat1_lcpm.names,columns=np.arange(len(pts))+1)
		d1.to_csv(fo,header=True,index=True,sep='\t')
		#Output network
		for xi in sparsities:
			output_folder1=pjoin(output_folder,'Full' if xi is None else f'sparsity={xi}')
			makedirs(output_folder1)
			stat1_net=stat.fbinarize(stat1_net0,sparsity=xi) if xi is not None else stat1_net0
			d1=stat1_net.compute(pts)
			for xj in range(len(pts)):
				fo=pjoin(output_folder1,str(xj+1)+'.tsv.gz')
				d2=d1[:,:,xj]
				#Remove empty regulators & targets
				t1=[(d2!=0).any(axis=1-x) for x in range(2)]
				if xi is not None and (any(not x.any() for x in t1) or d2[t1[0]][:,t1[1]].all()):
					#Skip sparsity for cell type because too large/small (all/no edges are positive)
					continue
				d2=pd.DataFrame(d2[t1[0]][:,t1[1]],index=stat1_net.names[0][t1[0]],columns=stat1_net.names[1][t1[1]])
				#Output to tsv.gz file
				d2.to_csv(fo,header=True,index=True,sep='\t')



assert __name__ != "__main__"
