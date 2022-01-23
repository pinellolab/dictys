#!/usr/bin/python3

__all__=['stat']

_docstring2argparse_ignore_=['stat','network','binarize']

from . import *

class network:
	"""
	Class for networks

	Initialization required variables
	----------------------------------
	cname:	numpy.ndarray(shape=(cn,),dtype=str)
		Cell names
	sname:	numpy.ndarray(shape=(sn,),dtype=str)
		State names
	msc:	numpy.ndarray(shape=(sn,cn),dtype=float)
		Cell-state weight
	nname:	numpy.ndarray(shape=(nn,),dtype=str)
		Node names
	nids:	[numpy.ndarray(shape=(nns[0],),dtype=int),numpy.ndarray(shape=(nns[1],),dtype=int)]
		Indices of source & target nodes

	Initialization optional variables
	----------------------------------
	cprop:	{name:numpy.ndarray(shape=(...,cn))}
		Properties of each cell
	nprop:	{name:numpy.ndarray(shape=(...,nn))}
		Constant node properties
	nsprop:	{name:numpy.ndarray(shape=(...,nn,ns))}
		State-dependent node properties
	eprop:	{name:numpy.ndarray(shape=(...,nns[0],nns[1]))}
		Constant edge properties
	esprop:	{name:numpy.ndarray(shape=(...,nns[0],nns[1],ns))}
		State-dependent edgeproperties	
	traj:	trajc
		Trajectory object

	Other variables
	-----------------
	cn:		int
		Number of cells
	cdict:	dict
		Map from cell name to id
	sn:		int
		Number of states
	sdict:	dict
		Map from state name to id
	nn:		int
		Number of nodes
	ndict:	dict
		Map from node name to id
	nns:	numpy.ndarray(shape=(2,),dtype=int)
		Number of source & target nodes
	"""
	_name_='network'
	_shape_={'e':lambda x:tuple(x.nns)}
	def _get_prop_shape_(self,propname):
		"""Get shape of property from _shape_, use default (prefix with n) if not available."""
		if propname.endswith('prop'):
			propname=propname[:-4]
		return self._shape_[propname](self) if propname in self._shape_ else (getattr(self,propname+'n'),)
	#Constructor
	def __init__(self,**ka):
		"""
		Network class constructor that should not be called directly. Please use class functions create_dynamic or create_group for creating networks from data, or from_file from file.
		ka:		Keyword arguments saved to class
		"""
		import numpy as np
		from collections import defaultdict
		#Storing values
		params_allowed=set('cname,sname,nname,nids,traj,prop'.split(','))
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
		#Validation
		self.check()
	def check(self):
		import numpy as np
		from dictys.traj import trajectory
		if not all([x[0]>0 and x[1].shape==(x[0],) and len(np.unique(x[1]))==x[0] and len(x[2])==x[0] and all([y in x[2] for y in x[1]]) for x in [[self.cn,self.cname,self.cdict],[self.sn,self.sname,self.sdict],[self.nn,self.nname,self.ndict]]]):
			raise ValueError('Cell/state/node names must be non-empty, unique, and matching their counts exactly.')
		for xi in self.prop:
			t1=tuple(np.concatenate([self._get_prop_shape_(x) for x in xi]))
			if not all([x.shape[-len(t1):]==t1 for x in self.prop[xi].values()]):
				raise ValueError('Cell/node/edge properties and their cross properties must have the correct dimensionality.')
		assert len(self.nids)==2 and all([len(np.unique(x))==len(x) for x in self.nids])
		assert all([(x>=0).all() and (x<self.nn).all() for x in self.nids])
		assert (self.nns==[len(x) for x in self.nids]).all() and (self.nns>0).all()
		assert self.traj is None or isinstance(self.traj,trajectory)
	#I/O
	@classmethod
	def from_file(cls,path):
		from collections import defaultdict
		from dictys.traj import trajectory
		import numpy as np
		import h5py
		props='cprop,nprop,nsprop,eprop,esprop'.split(',')
		params=dict()
		with h5py.File(path,'r') as f:
			for xi in filter(lambda x:x not in {'prop','traj'},f):
				params[xi]=np.array(f[xi])
				if params[xi].dtype.char=='S':
					params[xi]=params[xi].astype(str)
			params['prop']=defaultdict(dict)
			for xi in f['prop']:
				for xj in f['prop'][xi]:
					params['prop'][xi][xj]=np.array(f['prop'][xi][xj])
			if 'traj' in f:
				params['traj']=trajectory.from_fileobj(f['traj'])
		params['nids']=[params['nids1'],params['nids2']]
		del params['nids1'],params['nids2']
		return cls(**params)
	def to_file(self,path,compression="gzip",**ka):
		import h5py
		# props='cprop,nprop,nsprop,ecprop,esprop'.split(',')
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
					f['prop'][xi].create_dataset(xj,data=self.prop[xi][xj],compression=compression,**ka)
			if self.traj is not None:
				t=f.create_group('traj')
				self.traj.to_fileobj(t)
	def binarize(self,ntop,data='w',signed=False):
		"""
		Convert continuous to binary directed network.

		Parameters
		-------------
		data:		str or numpy.ndarray(shape=(n1,n2,...),dtype=float)
			Continuous directed network. Larger values are stronger and indicate higher ranking in binarization. Use str to represent entry in self.esprop.
		ntop:		float
			Number (for >=1)/proportion (for <1) of strongest edges to keep.
		signed:		bool
			Whether to regard continuous network as signed so the absolute values will be used for strength comparison.

		Returns
		---------
		net_bin:	numpy.ndarray(shape=(..,n1,n2),dtype=bool)
			Binary directed network
		"""
		import numpy as np
		if isinstance(data,str):
			data=self.esprop[data]
		if not signed:
			dr=np.abs(dr)
		if ntop<1:
			ntop=int(ntop*dr.size)
		assert ntop<dr.size
		t1=np.partition(dr.reshape(dr.shape[0]*dr.shape[1],*dr.shape[2:]),-ntop,axis=0)[-ntop]
		assert (t1!=0).all()
		return dr>=t1







assert __name__ != "__main__"
