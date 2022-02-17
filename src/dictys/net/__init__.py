#!/usr/bin/python3

"""
Module for network class and stats.

"""

from __future__ import annotations
__all__=['stat']

_docstring2argparse_ignore_=['stat','network','binarize']

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
			if not all(x.shape[-len(t1):]==t1 for x in self.prop[xi].values()):
				raise ValueError('Cell/node/edge properties and their cross properties must have the correct dimensionality.')
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
		params={}
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
		params={x:getattr(self,x) for x in 'cname,sname,nname,nids'.split(',')}
		params['nids1']=params['nids'][0]
		params['nids2']=params['nids'][1]
		del params['nids']
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
	# def binarize(self,ntop,data='w',signed=False):
	# 	"""
	# 	Convert continuous to binary directed network.

	# 	Parameters
	# 	-------------
	# 	data:		str or numpy.ndarray(shape=(n1,n2,...),dtype=float)
	# 		Continuous directed network. Larger values are stronger and indicate higher ranking in binarization. Use str to represent entry in self.esprop.
	# 	ntop:		float
	# 		Number (for >=1)/proportion (for <1) of strongest edges to keep.
	# 	signed:		bool
	# 		Whether to regard continuous network as signed so the absolute values will be used for strength comparison.

	# 	Returns
	# 	---------
	# 	net_bin:	numpy.ndarray(shape=(..,n1,n2),dtype=bool)
	# 		Binary directed network
	# 	"""
	# 	import numpy as np
	# 	if isinstance(data,str):
	# 		data=self.esprop[data]
	# 	if not signed:
	# 		dr=np.abs(dr)
	# 	if ntop<1:
	# 		ntop=int(ntop*dr.size)
	# 	assert ntop<dr.size
	# 	t1=np.partition(dr.reshape(dr.shape[0]*dr.shape[1],*dr.shape[2:]),-ntop,axis=0)[-ntop]
	# 	assert (t1!=0).all()
	# 	return dr>=t1







assert __name__ != "__main__"
