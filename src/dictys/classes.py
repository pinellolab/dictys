
class trajc:
	def __init__(self,edges,lens):
		"""Create trajectory object. Only supports fully connected tree trajectory.
		edges:	Start and end node IDs for each edge as np.array(shape=(n_edge,2))
		lens:	Length of each edge as np.array(shape=(n_edge,))
		"""
		import networkx as nx
		import numpy as np
		from collections import Counter
		self.ne=len(lens)
		assert self.ne>0
		assert edges.shape==(self.ne,2)
		assert (lens>=0).all()
		self.nn=edges.max()+1
		assert self.nn==self.ne+1
		assert not (edges[:,0]==edges[:,1]).any()
		assert frozenset(list(edges.ravel()))==frozenset(list(range(self.nn)))

		self.lens=lens
		self.edges=edges
		#Dictionary that converts node pairs to [edge ID,direction]
		self.edgedict={tuple(edges[x]):[x,1] for x in range(self.ne)}
		self.edgedict.update({tuple(edges[x][::-1]):[x,-1] for x in range(self.ne)})
		#Networkx graph
		self.g=nx.Graph([[x[0],x[1],{'l':y}] for x,y in zip(edges,lens)])
		try:
			nx.find_cycle(self.g)
			raise ValueError("Cyclic trajectory is not supported.")
		except nx.NetworkXNoCycle:
			pass
		if not nx.is_connected(self.g):
			raise ValueError("Only connected trajectory is supported.")
		#Distance between nodes
		self.dist=np.ones((self.nn,self.nn),dtype=lens.dtype)*np.inf
		t1=nx.shortest_path_length(self.g,weight='l')
		for xi in t1:
			t2=list(xi[1])
			self.dist[xi[0],t2]=[xi[1][y] for y in t2]
		assert (self.dist>=0).all()
		assert np.isfinite(self.dist).all()
		#Node degrees
		self.deg=np.zeros(self.nn,dtype=int)
		t1=np.array(list(Counter(self.edges.ravel()).items()))
		self.deg[t1[:,0]]=t1[:,1]
	@staticmethod
	def len_edge(dists,edges):
		"""Computes the length of each edge using sample distances to each node
		dists:	Distance matrix between each sample point on the graph and each node as np.array(shape=(n_sample,n_node))	
		edges:	Start and end node IDs for each edge as np.array(shape=(n_edge,2))
		Return:
		Length of each edge as np.array(shape=(n_edge,))
		"""
		import numpy as np
		ne=len(edges)
		ans=np.zeros(ne,dtype=float)
		for xi in range(ne):
			t1=np.abs(dists[:,edges[xi,0]]-dists[:,edges[xi,1]])
			ans[xi]=t1[np.isfinite(t1)].max()
		return ans
	@classmethod
	def fromdist(cls,edges,dists):
		"""Create trajectory object from distance matrix betwen sample points and nodes. Only supports fully connected tree trajectory.
		edges:	Start and end node IDs for each edge as np.array(shape=(n_edge,2))
		dists:	Distance matrix between each sample point on the graph and each node as np.array(shape=(n_sample,n_node))	
		"""
		return cls(edges,cls.len_edge(dists,edges))
	def topoint(self):
		return pointc.fromnodes(self)
	def conform_locs(self,locs,edges,rel_err=1E-7):
		"""Conform locs by clipping small deviations under float precision.
		locs:	np.array(shape=n) of locs to conform
		edges:	np.array(shape=edges) of respective edge ID each loc is on
		Return:
		Conformed locs as np.array(shape=n)
		"""
		import numpy as np
		errbound=rel_err*self.lens[edges]
		if (locs<-errbound).any() or (locs-self.lens[edges]>errbound).any():
			raise ValueError(f'Some locs are beyond relative error {rel_err}.')
		locs=np.clip(locs,0,self.lens[edges])
		return locs
	def points(self):
		"""Create list of points of nodes of the trajectory"""
		return pointc.fromnodes(self)
	def linspace(self,start,end,n):
		"""Find evenly spaced points on a path np.linspace
		start,
		end:	Start/end nodes' IDs to indicate the path.
		n:      Number of points including terminal nodes
		Return:
		Instance of pointc class with points go from start to end nodes.
		"""
		import numpy as np
		assert n>=2
		path=self.path(start,end)
		locs=np.linspace(0,self.lens[[self.edgedict[path[x],path[x+1]][0] for x in range(len(path)-1)]].sum(),n)
		return self.path_points(start,end,locs)
	def path_points(self,start,end,lengths):
		"""Find points at specific lengths on a path.
		start,
		end:	Start/end nodes' IDs to indicate the path.
		lengths:Lengths of movement from the starting node towards the ending node as np.array. Each length correspond to an output point.
				For lengths greater than total length of path, point at the end node will be returned.
		Return:
		Instance of pointc class with points go from start to end nodes.
		"""
		import numpy as np
		n=len(lengths)
		lengths0=lengths
		aorder=np.argsort(lengths)
		lengths=lengths[aorder]
		path=self.path(start,end)
		# lens=np.array([g.edges[path[x],path[x+1]][weight] for x in range(len(path)-1)]) if weight is not None else np.ones(len(path)-1,dtype=float)
		steps=np.zeros(n,dtype=int)
		nstart=0
		for xi in range(len(path)-2):
			l=self.lens[self.edgedict[path[xi],path[xi+1]][0]]
			t1=np.searchsorted(lengths[nstart:],l)
			if t1>=n:
				break
			lengths[nstart+t1:]-=l
			steps[nstart+t1:]=xi+1
			nstart+=t1
		lengths=np.array([lengths[x]*self.edgedict[path[steps[x]],path[steps[x]+1]][1]+self.lens[self.edgedict[path[steps[x]],path[steps[x]+1]][0]]*(1-self.edgedict[path[steps[x]],path[steps[x]+1]][1])/2 for x in range(n)])
		steps=np.array([self.edgedict[path[x],path[x+1]][0] for x in steps])
		lengths=self.conform_locs(lengths,steps)
		ans_steps=steps.copy()
		ans_lengths=lengths.copy()
		ans_steps[aorder]=steps
		ans_lengths[aorder]=lengths
		return pointc(self,ans_steps,ans_lengths)
	def path(self,start,end):
		"""Find path from start to end node as list of node IDs
		start,
		end:	Start/end nodes' IDs to indicate the path.
		Return:
		np.array of edge IDs to go from start to end nodes (inclusive).
		"""
		import networkx as nx
		import numpy as np
		assert start>=0 and start<self.nn
		assert end>=0 and end<self.nn
		assert start!=end
		return np.array(nx.shortest_path(self.g,start,end))		
	def smoothen(self,data,axis,func_name,a,ka,points):
		"""Function to compute smoothened/interpolated values on given data at provided points. Used by self.smoothened.
		data:		np.array(shape=(...,n_node,...) of base data to smoothen/interpolate from.
		axis:		Axis of data to smoothen. data must have length n_node at this axis.
		func_name:	Name of smoothening function of class pointc. Using weight_conv should provide func_name='conv',
		a:			Arguments of smoothening function.
		ka:			Keyword arguments smoothening function.
		points:		pointc instance to interpolate data for
		Return:
		Smoothened values that have the same dimensionality as data in all axes other than axis. For axis, the dimension is n_pts.
		"""
		func=getattr(points,'weight_'+func_name)
		w=func(*a,**ka).T
		assert w.ndim==2
		 # and w.shape[0]==w.shape[1]
		w/=w.sum(axis=0)
		ans=(data.swapaxes(axis,-1)@w).swapaxes(axis,-1)
		assert ans.shape[:axis]==data.shape[:axis]
		assert ans.shape[axis+1]==data.shape[axis+1]
		return ans
	def smoothened(self,data,func_name,*a,axis=-1,**ka):
		"""Create a smoothened/interpolated function on given data that computes values at provided points
		data:		np.array(shape=(...,n_node,...) of base data to smoothen/interpolate from.
		func_name:	Name of smoothening function of class pointc. Using weight_conv should provide func_name='conv',
		a:			Arguments of smoothening function.
		axis:		Axis of data to smoothen. data must have length n_node at this axis.
		ka:			Keyword arguments smoothening function.
		Return:
		function(points)=smoothened values.
		Smoothened values have the same dimensionality as data in all axes other than axis. For axis, the dimension is n_pts.
		"""
		from functools import partial
		assert data.shape[axis]==self.nn
		return partial(self.smoothen,data,axis,func_name,a,ka)
	def terminal_nodes(self):
		"""Return: IDs of terminal nodes (degree=1) as np.array."""
		from collections import Counter
		import numpy as np
		ans=Counter(self.edges.ravel()).items()
		ans=np.array([x[0] for x in ans if x[1]==1])
		return ans

class networkc:
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
	cprop:	{name:numpy.ndarray(shape=(cn,...))}
		Properties of each cell
	npropc:	{name:numpy.ndarray(shape=(nn,...))}
		Constant node properties
	nprops:	{name:numpy.ndarray(shape=(ns,nn,...))}
		State-dependent node properties
	epropc:	{name:numpy.ndarray(shape=(nns[0],nns[1],...))}
		Constant edge properties
	eprops:	{name:numpy.ndarray(shape=(ns,nns[0],nns[1],...))}
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
	def __init__(self,**ka):
		"""
		Network class constructor that should not be called directly. Please use class functions create_dynamic or create_group for creating networks from data, or from_file from file.
		ka:		Keyword arguments saved to class
		"""
		import numpy as np
		#Storing values
		params_allowed=set('cname,cprop,sname,msc,nname,npropc,nprops,nids,epropc,eprops,traj',split(','))
		for xi in ka:
			if xi not in params_allowed:
				raise ValueError(f'Unknown parameter {xi}.')
			setattr(self,xi,ka[xi])
		#Post processing
		for xi in [['c','cell'],['s','state'],['n','node']]:
			vnames=[xi[0]+'n',xi[0]+'name',xi[0]+'dict']
			#Fill count/name with another if needed
			setattr(self,vnames[0],len(getattr(self,vnames[1]))
			#Prepare dictionary
			setattr(self,vnames[2],dict(zip(getattr(self,vnames[1]),range(len(getattr(self,vnames[0]))))))
		#Fill count if needed
		self.nns=np.array([len(x) for x in self.nids])
		for xi in 'cprop,npropc,nprops,epropc,eprops',split(','):
			#Initialize with empty property list
			if not hasattr(self,xi):
				setattr(self,xi,dict())
		#Validation
		self.check()
	def check(self):
		if not all([x[0]>0 and x[1].shape==(x[0],) and len(np.unique(x[1]))==x[0] and len(x[2])==x[0] and all([y in x[2] for y in x[1]]) for x in [[self.cn,self.cname,self.cdict],[self.sn,self.sname,self.sdict],[self.nn,self.nname,self.ndict]]):
			raise ValueError('Cell/state/node names must be non-empty, unique, and matching their counts exactly.')
		if not all([all([x.shape[:len(y[1])]==y[1] for x in y[0].values()]) for y in [[self.cprop,(self.cn,)],[self.npropc,(self,nn,)],[self.nprops,(self.ns,self.nn)],[self.epropc,(self.nns[0],self.nns[1])],[self.eprops,(self.ns,self.nns[0],self.nns[1])]]):
			raise ValueError('Cell/node/edge properties must have the correct dimensionality.')
		assert self.msc.shape==(self.sn,self.cn)
		assert self.nids.shape==(2,) and all([len(np.unique(x))==len(x) for x in self.nids])
		assert (self.nids>=0).all() and (self.nids<self.nn).all()
		assert self.nns==[len(x) for x in self.nids] and (self.nns>0).all()
		assert self.traj is None or isinstance(self.traj,trajc)
	@classmethod
	def from_file(cls,file_net,file_traj=None):
		import h5py
		props='cprop,npropc,nprops,epropc,eprops',split(',')
		params=dict()
		with h5py.File(file_net,'r') as f:
			for xi in filter(lambda x:x not in props,f):
				params[xi]=f[xi]
		for namepref in props:
			params[namepref]=dict()
			for xi in f[namepref]:
				params[namepref][xi]=f[namepref][xi]
		params['nids']=[params['nids1'],params['nids2']]
		del params['nids1'],params['nids2']
		if file_traj is not None:
			raise NotImplementedError
		return cls(**params)
	def to_file(self,path):
		import h5py
		props='cprop,npropc,nprops,epropc,eprops',split(',')
		params={x:getattr(self,x) for x in 'cname,sname,msc,nname,nids'.split(',')}
		params['nids1']=params['nids'][0]
		params['nids2']=params['nids'][1]
		del params['nids']
		params_grp={x:getattr(self,x) for x in props}
		with h5py.File(path,'w') as f:
			for xi in params:
				f.create_dataset(xi,data=params[xi])
			for xi in params_grp:
				f.create_group(xi)
				for xj in params_grp[xi]:
					f[xi].create_dataset(xj,data=params_grp[xi][xj])












































#
