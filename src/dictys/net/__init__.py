#!/usr/bin/python3

"""
Module for network class and stats.

"""

from __future__ import annotations
from typing import Set,Optional
__all__=['stat']

from . import *


class network:
	"""
	Class for context-specific and dynamic networks

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
			Vluae of each property. The trailing axes in its shape must exactly match the shape factor from `type`.
	traj:	dictys.traj.trajectory
		Trajectory object. Only needed for dynamic network
	point:	{x:dictys.traj.point for x in ['c','s']}
		Point objects for cells and states on the `traj` object. Only needed for dynamic network.
	"""
	_name_='network'
	_shape_={'e':lambda x:tuple(x.nns)}
	def _get_prop_shape_(self,propname:str):
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
		ka:		Keyword arguments saved to class
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
		from dictys.traj import trajectory
		if not all(x[0]>0 and x[1].shape==(x[0],) and len(np.unique(x[1]))==x[0] and len(x[2])==x[0] and all(y in x[2] for y in x[1]) for x in [[self.cn,self.cname,self.cdict],[self.sn,self.sname,self.sdict],[self.nn,self.nname,self.ndict]]):
			raise ValueError('Cell/state/node names must be non-empty, unique, and matching their counts exactly.')
		for xi in self.prop:
			t1=tuple(np.concatenate([self._get_prop_shape_(x) for x in xi]))
			t2=list(filter(lambda x:self.prop[xi][x].shape[-len(t1):]!=t1,self.prop[xi]))
			if len(t2)>0:
				raise ValueError('Properties of type {} have incorrect size. Expected trailing dimensions: {}. Actual dimensions: '.format(xi,t1)+', '.join(['{} ({})'.format(tuple(self.prop[xi][x].shape),x) for x in t2]))
		assert len(self.nids)==2 and all(len(np.unique(x))==len(x) for x in self.nids)
		assert all((x>=0).all() and (x<self.nn).all() for x in self.nids)
		assert (self.nns==[len(x) for x in self.nids]).all() and (self.nns>0).all()
		assert self.traj is None or isinstance(self.traj,trajectory)
		assert (self.traj is None and len(self.point)==0) or (self.traj is not None and frozenset(self.point)==frozenset({'c','s'}))
		assert all(len(self.point[x])==getattr(self,x+'n') for x in self.point)
	#I/O
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
		return cls(**params)
	def to_file(self,path:str,compression:str="gzip",**ka)->None:
		"""
		Save class instance to file.
		
		Parameters
		-------------
		path:	str
			Path of file to save to.
		compression: str
			Compression used.
		ka:	dict
			Keyword arguments passed to h5py.File.create_dataset.
		"""
		import h5py
		import logging
		params={x:getattr(self,x) for x in 'cname,sname,nname,nids'.split(',')}
		params['nids1']=params['nids'][0]
		params['nids2']=params['nids'][1]
		del params['nids']
		logging.info(f'Writing file {path}.')
		with h5py.File(path,'w') as f:
			for xi in params:
				p=dict(ka)
				p['compression']=compression
				p['data']=params[xi]
				if params[xi].dtype.char=='U':
					p['data']=p['data'].astype('S')
				f.create_dataset(xi,**p)
			f.create_group('prop')
			for xi in self.prop:
				f['prop'].create_group(xi)
				for xj in self.prop[xi]:
					data=self.prop[xi][xj]
					if data.dtype.char=='U':
						data=data.astype('S')
					f['prop'][xi].create_dataset(xj,data=data,compression=compression,**ka)
			if self.traj is not None:
				t=f.create_group('traj')
				self.traj.to_fileobj(t)
				t=f.create_group('point')
				for xi in self.point:
					self.point[xi].to_fileobj(t.create_group(xi))
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

			* '':	Unnormalized direct network.

			* 'n':	Normalized direct network.

			* 'i':	Unnormalized steady-state network.

			* 'in':	Normalized steady-state network.

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
		from dictys.utils.file import read_txt
		from dictys.utils.numpy import dtype_min
		import numpy as np
		import pandas as pd

		optional_allowed={'readcount'}
		if len(optional-optional_allowed)>0:
			raise ValueError('Unknown optional data types: '+','.join(optional-optional_allowed))

		params={}
		#Cells & states/subsets
		sname=read_txt(fi_subsets,unique=True)
		n=len(sname)
		if n==0:
			raise RuntimeError('Could not find cell subsets.')
		success=np.ones(n,dtype=bool)
		cnames=[]
		for xi in range(n):
			cnames.append(read_txt(pjoin(path_work,sname[xi],'names_rna.txt')))
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
		nets=[None]*n
		for xi in range(n):
			try:
				fi=pjoin(path_work,sname[xi],f'net_{nettype}weight.tsv.gz')
				logging.info(f'Reading file {fi}')
				nets[xi]=pd.read_csv(fi,header=0,index_col=0,sep='\t')
				if nets[xi].shape[0]==0:
					success[xi]=False
			except FileNotFoundError as e:
				logging.warning('Skipping cell subset {} due to error: {}'.format(sname[xi],repr(e)))
				success[xi]=False
				nets[xi]=pd.DataFrame([],index=[],columns=[])

		#Network nodes: source & target
		nidss=[[x.index,x.columns] for x in nets]
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
		#Edge weight
		esprop['w']=np.zeros(np.r_[nns,ns],dtype=float)
		#Edge mask
		esprop['mask']=np.zeros(np.r_[nns,ns],dtype=bool)
		for xi in np.nonzero(success)[0]:
			t1=[[niddicts[y][x] for x in nidss[xi][y]] for y in range(2)]
			esprop['w'][np.ix_(t1[0],t1[1],[xi])]=nets[xi].values.reshape(*nets[xi].shape,1)
			esprop['mask'][np.ix_(t1[0],t1[1],[xi])]=True

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

		if success.all():
			return cls(**params)
		if not success.any():
			raise RuntimeError('All subsets failed.')
		
		#Remove failed states
		# logging.warning('Removing {}/{} failed subsets.'.format((~success).sum(),len(success)))
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




assert __name__ != "__main__"
