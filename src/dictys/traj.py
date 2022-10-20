#!/usr/bin/python3
# Lingfei Wang, 2022. All rights reserved.


"""
Classes for trajectory and points on trajectory
"""

from __future__ import annotations
from typing import Union,Tuple,Callable,Optional
import h5py
from dictys.utils.numpy import NDArray,ArrayLike

_docstring2argparse_ignore_=['argpartition','trajectory','point']

def argpartition(a:NDArray,kth:int,axis:int=-1,draw_order:str='undefined')->NDArray:
	"""
	Performs numpy.argpartition with options to specify how to handle ties.

	Parameters
	----------
	a:			numpy.ndarray
		Array to perform argpartition
	kth:		int
		Index to partition
	axis:		int
		Axis to partition on
	draw_order:	str
		How to handle ties. Accepts:
		* undefined:	Default behavior of numpy.argpartition
		* random:		Random draw from ties
		* error:		Raises a RuntimeError

	Returns
	-------
	numpy.ndarray
		Result of numpy.argpartition with custom handling of ties.
	"""
	import numpy as np
	assert draw_order in {'undefined','random','error'}
	shape=a.shape
	a2=np.swapaxes(a,axis,0)
	shape2=a2.shape
	a2=a2.reshape(shape2[0],np.prod(shape2[1:],dtype=type(shape[0])))
	ans=np.argpartition(a2,kth,axis=0)
	if draw_order=='undefined':
		ans=np.swapaxes(ans.reshape(*shape2),axis,0)
		return ans
	if draw_order=='random':
		t1=a2==a2[ans[kth],np.arange(ans.shape[1])]
		t1=[np.nonzero(x)[0] for x in t1.T]
		for xi in filter(lambda x:len(t1[x])>1,range(ans.shape[1])):
			t2=t1[xi].copy()
			np.random.shuffle(t2)
			ans[t1[xi],xi]=ans[t2,xi]
	elif draw_order=='error':
		t1=a2==a2[ans[kth],np.arange(ans.shape[1])]
		t1=[np.nonzero(x)[0] for x in t1.T]
		t1=list(filter(lambda x:len(t1[x])>1,range(ans.shape[1])))
		if len(t1)>0:
			raise RuntimeError('Found duplicate values at the partition point.')
	else:
		raise NotImplementedError
	ans=np.swapaxes(ans.reshape(*shape2),axis,0)
	return ans

class trajectory:
	"""
	Class for trajectory of cell state transitions on low dimension
	"""
	def __init__(self,edges:NDArray[int],lens:NDArray[float])->None:
		"""
		Create state trajectory object. Only supports fully connected tree trajectory.

		Parameters
		----------
		edges:	numpy.ndarray(shape=(n_edge,2))
			Start and end state node IDs for each edge
		lens:	numpy.ndarray(shape=(n_edge,))
			Length of each edge

		Attributes
		----------
		nn:			int
			Node count
		ne:			int
			Edge count
		edgedict:	dict
			Dictionary that converts node pairs to [edge ID,direction]
		g:			networkx.Graph
			Networkx graph
		dist:		numpy.ndarray(shape=(nn,nn))
			Distance between nodes
		deg:		numpy.ndarray(shape=(nn,))
			Node degrees
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
	def len_edge(dists:NDArray[float],edges:NDArray[int])->NDArray[float]:
		"""
		Computes the length of each edge using points' distances to each state node

		Parameters
		----------
		dists:	numpy.ndarray(shape=(n_sample,n_node))
			Distance matrix between each sample point on the graph and state each node
		edges:	numpy.ndarray(shape=(n_edge,2))
			Start and end node IDs for each edge

		Returns
		----------
		lengths:	numpy.ndarray(shape=(n_edge,))
			Length of each edge
		"""
		import numpy as np
		ne=len(edges)
		ans=np.zeros(ne,dtype=float)
		for xi in range(ne):
			t1=np.abs(dists[:,edges[xi,0]]-dists[:,edges[xi,1]])
			ans[xi]=t1[np.isfinite(t1)].max()
		return ans
	@classmethod
	def fromdist(cls,edges:NDArray[int],dists:NDArray[float])->trajectory:
		"""
		Create trajectory object from distance matrix betwen sample points and state nodes. Only supports fully connected tree trajectory.

		Parameters
		----------
		edges:	numpy.ndarray(shape=(n_edge,2))
			Start and end node IDs for each edge
		dists:	numpy.ndarray(shape=(n_sample,n_node))
			Distance matrix between each sample point on the graph and each node
		"""
		return cls(edges,cls.len_edge(dists,edges))
	def topoint(self)->point:
		"""
		Convert trajectory object to point object on this trajectory.

		Returns
		----------
		dictys.traj.point
			Point object converted.
		"""
		return point.fromnodes(self)
	def conform_locs(self,locs:NDArray[float],edges:NDArray[int],abs_err:float=1E-7,rel_err:float=1E-7)->NDArray:
		"""
		Conform point locations by clipping small deviations under float precision. An error will be trigger if the relative and absolute errors both exceed the given bounds.

		Parameters
		----------
		locs:	numpy.ndarray(shape=n)
			Points' locations on each edge to conform
		edges:	numpy.ndarray(shape=n)
			Points' edges to conform
		abs_err:
			Bound for absolute errors
		rel_err:
			Bound for relative errors

		Returns
		----------
		numpy.ndarray(shape=n)
			Conformed locs
		"""
		import numpy as np
		t1=self.lens[edges]+1E-300
		eabs=np.max([-locs,locs-self.lens[edges]],axis=0)
		erel=eabs/t1
		t1=np.nonzero((eabs>abs_err)&(erel>rel_err))[0]
		if len(t1)>0:
			raise ValueError(f'Absolute error {eabs[t1[0]]}>{abs_err} and relative error {erel[t1[0]]}>{rel_err} both exceed the given bounds.')
		locs=np.clip(locs,0,self.lens[edges])
		return locs
	def linspace(self,start:int,end:int,n:int)->point:
		"""
		Find evenly spaced points on a path like np.linspace

		Parameters
		----------
		start:	int
			Start nodes' IDs to indicate the path
		end:	int
			End nodes' IDs to indicate the path
		n:		int
			Number of points including terminal nodes

		Returns
		----------
		dictys.traj.point
			Instance of point class with points go from start to end nodes.
		"""
		import numpy as np
		assert n>=2
		path=self.path(start,end)
		locs=np.linspace(0,self.lens[[self.edgedict[path[x],path[x+1]][0] for x in range(len(path)-1)]].sum(),n)
		return self.path_points(start,end,locs)
	def path_points(self,start:int,end:int,lengths:ArrayLike)->point:
		"""
		Find points at specific lengths on a path.

		Parameters
		----------
		start:	int
			Start nodes' IDs to indicate the path
		end:	int
			End nodes' IDs to indicate the path
		lengths:	numpy.ndarray(dtype=float)
			Lengths of movement from the starting node towards the ending node as numpy.ndarray. Each length correspond to an output point.
			For lengths with negative values, the start node will be returned.
			For lengths greater than total length of path, the end node will be returned.
		Returns
		----------
		dictys.traj.point
			Instance of point class with points go from start to end nodes.
		"""
		import numpy as np
		n=len(lengths)
		# lengths0=lengths
		aorder=np.argsort(lengths)
		lengths=lengths[aorder]
		lengths[lengths<0]=0
		#Nodes on the path
		path=self.path(start,end)
		# lens=np.array([g.edges[path[x],path[x+1]][weight] for x in range(len(path)-1)]) if weight is not None else np.ones(len(path)-1,dtype=float)
		steps=np.zeros(n,dtype=int)
		nstart=0
		#Convert length from starting node to intermediate node and length from that node
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
		return point(self,ans_steps,ans_lengths)
	def path(self,start:int,end:int)->NDArray[int]:
		"""
		Find path from start to end node as list of node IDs

		Parameters
		----------
		start:	int
			Path sttart node's index
		end:	int
			Path end node's index

		Returns
		----------
		numpy.ndarray
			Node indices to go from start to end nodes (inclusive).
		"""
		import networkx as nx
		import numpy as np
		assert 0<=start<self.nn
		assert 0<=end<self.nn
		assert start!=end
		return np.array(nx.shortest_path(self.g,start,end))
	def smoothened(self,data:NDArray,*a,axis:int=-1,nodes:Optional[list]=None,nodes_path:Optional[Tuple[int,int]]=None,criterion_path:str='strict',**ka)->NDArray:
		"""
		Create a smoothened/interpolated function of given data on trajectory nodes that computes values at provided points.

		Parameters
		----------
		data:			numpy.ndarray(shape=(...,n_node,...))
			Base data to smoothen/interpolate from.
		a:				list
			Arguments of smoothening function.
		axis:			int
			Axis of data to smoothen. data must have length n_node at this axis.
		nodes:			list or None
			List of node indices to specify the subset of nodes used.
		nodes_path:		(int,int) or None 
			Start and end node indices to specify a path for automaticaly selecting `nodes` parameter. See `criterion_path` on selection criterion.
		criterion_path:	str
			How to use `nodes_path` to filter nodes for `nodes` parameter. Ignored if `nodes_path` is None. Accepts two options:
			* strict:	Only uses nodes on the path.

			* loose:	Excludes nodes on any other branch of the path.

		ka:			dict
			Keyword arguments smoothening function.

		Returns
		----------
		function(point)->numpy.ndarray(shape=(...,len(point),...))
			Function to compute smoothened data matrix at any points. It has the same dimensionality as data in all axes other than axis.
		"""
		import itertools
		import numpy as np
		import networkx as nx
		pts=self.topoint()
		assert data.shape[axis]==len(pts)
		if nodes_path is not None:
			assert nodes is None
			assert len(nodes_path)==2
			if criterion_path=='strict':
				nodes=self.path(*nodes_path)
			elif criterion_path=='loose':
				#Find all nodes in other branches of the path
				nodes=self.path(*nodes_path)
				#Get connectivity of modified graph
				g=self.g.copy()
				g.remove_edges_from(zip(nodes[:-1],nodes[1:]))
				dist=dict(nx.all_pairs_bellman_ford_path_length(g))
				#Remove nodes connected to intermediate nodes on the path
				nodes_exclude=set(itertools.chain.from_iterable([dist[x].keys() for x in nodes[1:-1]]))-set(nodes)
				nodes=list(filter(lambda x:x not in nodes_exclude,range(self.nn)))
			else:
				raise ValueError(f'Unknown value {criterion_path} for criterion_path.')
		if nodes is not None:
			assert np.min(nodes)>=0 and np.max(nodes)<self.nn
			nodes=set(nodes)
			nodes=np.array([x in nodes for x in range(self.nn)])
			assert data.shape[axis]==len(pts)
			pts=pts[nodes]
			data=data.swapaxes(axis,0)[nodes].swapaxes(axis,0)
		assert data.shape[axis]==len(pts)
		return pts.smoothened(data,*a,axis=axis,**ka)
	def terminal_nodes(self)->NDArray[int]:
		"""
		Finds the terminal nodes (degree=1).

		Returns
		----------
		numpy.ndarray
			Terminal node indices
		"""
		from collections import Counter
		import numpy as np
		ans=Counter(self.edges.ravel()).items()
		ans=np.array([x[0] for x in ans if x[1]==1])
		return ans
	#I/O
	@classmethod
	def from_fileobj(cls,f:Union[h5py.File,h5py.Group])->trajectory:
		"""
		Load object from file object

		Parameters
		----------
		f:	fileobj
			File object to load from

		Returns
		----------
		dictys.traj.trajectory
			Loaded class object
		"""
		import numpy as np
		params=[np.array(f['edges']),np.array(f['lens'])]
		return cls(*params)
	@classmethod
	def from_file(cls,path:str)->trajectory:
		"""
		Load object from file

		Parameters
		----------
		path:	str
			File path to load from

		Returns
		----------
		dictys.traj.trajectory
			Loaded class object
		"""
		import logging
		logging.info(f'Reading file {path}.')
		with h5py.File(path,'r') as f:
			return cls.from_fileobj(f)
	def to_fileobj(self,f:Union[h5py.File,h5py.Group],compression:str="gzip",**ka)->None:
		"""
		Save object to file object

		Parameters
		----------
		f:	fileobj
			File object to save to
		compression: str
			Type of compression. See h5py.File.create_dataset.
		ka:	dict
			Keyword arguments passed to h5py.File.create_dataset
		"""
		p=dict(ka)
		p['compression']=compression
		p['data']=self.lens
		f.create_dataset('lens',**p)
		p=dict(ka)
		p['compression']=compression
		p['data']=self.edges
		f.create_dataset('edges',**p)
	def to_file(self,path:str,**ka)->None:
		"""
		Save object to file

		Parameters
		----------
		path:		str
			File path to save to.
		ka:			dict
			Keyword arguments passed to self.to_fileobj
		"""
		import logging
		logging.info(f'Writing file {path}.')
		with h5py.File(path,'w') as f:
			return self.to_fileobj(f,**ka)

class point:
	"""
	Class for points on a trajectory on low dimension
	"""
	def __init__(self,traj:trajectory,edges:NDArray[int],locs:NDArray[float],dist:Optional[NDArray[float]]=None):
		"""
		Point list on trajectory.

		Parameters
		----------
		traj:	dictys.traj.trajectory
			Trajectory the points are on
		edges:	numpy.ndarray(shape=(n_point))
			Edge ID each point is on
		locs:	numpy.ndarray(shape=(n_point))
			Distance from the starting node on its edge for each point
		dist:	numpy.ndarray(shape=(n_point,n_node))
			Distance matrix between each point and each node. Automatically computed if not provided.

		Attributes
		----------
		npt:	int
			Number of points
		p:		dictys.traj.trajectory
			Same as `traj`
		dist:	numpy.ndarray(shape=(npt,p.nn))
			Distance between each point and each node of the trajectory
		"""
		if locs is None and dist is None:
			raise TypeError('At least one of locs and dist must be specified')
		assert edges.ndim==1
		self.npt=len(edges)
		assert isinstance(traj,trajectory)
		assert (edges>=0).all() and (edges<traj.ne).all()
		self.p=traj
		self.edges=edges
		assert locs.shape==(self.npt,)
		locs=traj.conform_locs(locs,edges)
		self.locs=locs
		if dist is None:
			dist=self.compute_dist()
		assert dist.shape==(self.npt,self.p.nn)
		assert (dist>=0).all()
		self.dist=dist
	@classmethod
	def fromdist(cls,traj:trajectory,edges:NDArray[int],dist:NDArray[float])->point:
		"""
		Class conctructor from edges and distance to all nodes of each point.

		Parameters
		----------
		traj:	dictys.traj.trajectory
			Trajectory the points are on
		edges:	numpy.ndarray(shape=(n_point) or (n_point,2))
			Edge each point is on, as edge ID (for shape=(n_point)), or as two node IDs (for shape=(n_point,2))
		dist:	numpy.ndarray(shape=(n_point,n_node))
			Distance matrix between each point and each node. Automatically computed if not provided.

		Returns
		----------
		dictys.traj.point
			Constructed point object
		"""
		import numpy as np
		assert edges.ndim in {1,2}
		assert dist.shape==(len(edges),traj.nn)
		if edges.ndim==2:
			#Convert to edge ID
			if any(tuple(x) not in traj.edgedict for x in edges):
				raise ValueError('Found edges absent in trajectory.')
			edges=np.array([traj.edgedict[tuple(x)][0] for x in edges])
		locs=dist[np.arange(len(edges)),traj.edges[edges,0]]
		return cls(traj,edges,locs,dist=dist)
	@classmethod
	def fromnodes(cls,traj:trajectory)->point:
		"""
		Convert trajectory nodes to points on the trajectory

		Parameters
		----------
		traj:	dictys.traj.trajectory
			Trajectory to convert

		Returns
		----------
		dictys.traj.point
			Constructed point object
		"""
		import numpy as np
		nodes,indices=np.unique(traj.edges,return_index=True)
		indices=[indices//2,indices%2]
		edges=np.zeros(traj.nn,dtype=int)
		edges[nodes]=indices[0]
		locs=np.zeros(traj.nn,dtype=traj.lens.dtype)
		locs[nodes]=indices[1]*traj.lens[edges]
		return cls(traj,edges,locs,dist=traj.dist)
	def copy(self)->point:
		"""
		Copy of self

		Returns
		-------
		dictys.traj.point
			Copied object
		"""
		return self.__class__(self.p,self.edges,self.locs,dist=self.dist)
	@classmethod
	def concat(cls,pts:list[point])->point:
		"""
		Concatenate point lists to a single point list

		Parameters
		----------
		pts:	List of point lists

		Returns
		-------
		Concatenated point list

		"""
		import numpy as np
		if not all(isinstance(x,cls) for x in pts):
			raise TypeError('All elements of pts should have point class.')
		if not all(x.p==pts[0].p for x in pts[1:]):
			raise TypeError('All elements of pts should be from the same trajectory.')
		if len(pts)==1:
			return pts[0].copy()
		return cls(pts[0].p,np.concatenate([x.edges for x in pts]),np.concatenate([x.locs for x in pts]),dist=np.concatenate([x.dist for x in pts],axis=0))
	def __add__(self,other:point)->point:
		return self.concat([self,other])
	def __sub__(self,other:point)->NDArray[float]:
		"""
		Subtraction computes the distance matrix between all point pairs in two point lists.

		Parameters
		----------
		other:	dictys.traj.point
			Destination point list.

		Returns
		----------
		numpy.ndarray(shape=(self.npt,other.npt))
			Distance matrix
		"""
		import numpy as np
		from dictys.utils.numpy import groupby
		assert isinstance(other,self.__class__)
		p=self.p
		assert other.p==p

		ans=-np.ones((self.npt,other.npt),dtype=float)
		#Cells on the same edge
		grp=groupby(self.edges)
		for xi in grp:
			branchsame=other.edges==xi
			if branchsame.any():
				#Cells on the same edge
				t1=np.nonzero(branchsame)[0]
				ans[np.ix_(grp[xi],t1)]=np.abs(np.repeat(self.locs[grp[xi]].reshape(-1,1),len(t1),axis=1)-other.locs[t1])
			if ~branchsame.all():
				#Cells on a different edge
				t1=np.nonzero(~branchsame)[0]
				#Four possible paths from self to other
				t2=[
					# self-start-start-other
					(np.repeat(p.dist[p.edges[xi,0],p.edges[other.edges[t1],0]].reshape(1,-1),len(grp[xi]),axis=0).T+self.locs[grp[xi]]).T+other.locs[t1],
					# self-start-end-other
					(np.repeat(p.dist[p.edges[xi,0],p.edges[other.edges[t1],1]].reshape(1,-1),len(grp[xi]),axis=0).T+self.locs[grp[xi]]).T+(p.lens[other.edges[t1]]-other.locs[t1]),
					# self-end-start-other
					(np.repeat(p.dist[p.edges[xi,1],p.edges[other.edges[t1],0]].reshape(1,-1),len(grp[xi]),axis=0).T+(p.lens[xi]-self.locs[grp[xi]])).T+other.locs[t1],
					# self-end-end-other
					(np.repeat(p.dist[p.edges[xi,1],p.edges[other.edges[t1],1]].reshape(1,-1),len(grp[xi]),axis=0).T+(p.lens[xi]-self.locs[grp[xi]])).T+(p.lens[other.edges[t1]]-other.locs[t1]),
				]
				ans[np.ix_(grp[xi],t1)]=np.min(t2,axis=0)
		assert ans.shape==(self.npt,other.npt)
		assert (ans>=0).all()
		return ans
	def compute_dist(self)->NDArray[float]:
		"""
		Compute distance to every node on trajectory given locs

		Returns
		----------
		numpy.ndarray(shape=(n_point,n_node))
			Distance matrix
		"""
		return self-self.fromnodes(self.p)
	def perturb(self,scale:float=1)->None:
		"""
		Perturb locations of points without changing order of points and nodes in the list. Overlapping points will have random order.

		Parameters
		----------
		scale:	Scale of perturbation relative to distance between points/nodes. Must be <=1 because >1 will change the order.
		"""
		import numpy as np
		from collections import defaultdict
		import itertools
		assert scale<=1

		dist=self.dist
		#Pre-filter duplicate locations with hash
		t1=[hash(tuple(x)) for x in dist]
		t2=defaultdict(list)
		for xi in range(self.npt):
			t2[t1[xi]].append(xi)
		grp=[list(x) for x in t2.values() if len(x)>1]
		#Specifically filter duplicate locations with exact comparison
		t1=[]
		for xi in grp:
			t2=[x for x in filter(lambda y:y[1]>y[0],itertools.product(xi,xi)) if (dist[x[0]]==dist[x[1]]).all()]
			t3=defaultdict(list)
			for xj in t2:
				t3[xj[0]].append(xj[1])
				t3[xj[1]].append(xj[0])
			t3=list({frozenset(x[1]+[x[0]]) for x in t3.items()})
			t1+=t3
		grp=[list(x) for x in t1]
		t1=list(itertools.chain.from_iterable(grp))
		assert len(t1)==len(set(t1))
		if len(grp)==0:
			return
		#Forcing identical branch assignment (only useful for distance=0 groups)
		self.edges=self.edges.copy()
		self.locs=self.locs.copy()
		for xi in grp:
			self.edges[xi[1:]]=self.edges[xi[0]]
			self.locs[xi[1:]]=self.locs[xi[0]]

		#Apply perturbation
		perturb_amount=np.zeros(self.npt,dtype=self.locs.dtype)
		for xi in grp:
			#Get min/max bound of perturbation
			edge=self.edges[xi[0]]
			nodes=self.p.edges[edge]
			t1=np.nonzero(self.edges==edge)[0]
			t2=np.sort(np.unique(dist[t1,nodes[0]]))
			t3=np.searchsorted(t2,dist[xi[0],nodes[0]])
			bound=[]
			if t3==0:
				bound.append(-dist[xi[0],nodes[0]])
			else:
				bound.append((t2[t3-1]-t2[t3])/2)
			if t3==len(t2)-1:
				bound.append(dist[xi[0],nodes[1]])
			else:
				bound.append((t2[t3+1]-t2[t3])/2)
			#Produce perturbations
			perturb_amount[xi]=np.random.rand(len(xi))*(bound[1]-bound[0])*scale+bound[0]
		self.locs=self.p.conform_locs(self.locs+perturb_amount,self.edges)
		self.dist=self.compute_dist()
	def path(self,start:int,end:int)->NDArray[int]:
		"""
		Find path from start to end points as list of node IDs

		Parameters
		----------
		start:	int
			Path start point's index
		end:	int
			Path end point's index

		Returns
		----------
		numpy.ndarray
			Node indices to go from start to end points.
		"""
		import networkx as nx
		import numpy as np
		assert 0<=start<len(self)
		assert 0<=end<len(self)
		assert start!=end
		if self.edges[start]==self.edges[end]:
			return np.array([])
		#Insert points with edges to graph
		g=self.p.g.copy()
		g.add_edges_from([(self.p.edges[self.edges[start]][0],self.p.nn),(self.p.edges[self.edges[start]][1],self.p.nn),
			(self.p.edges[self.edges[end]][0],self.p.nn+1),(self.p.edges[self.edges[end]][1],self.p.nn+1)])
		#Remove past edges
		g.remove_edges_from([tuple(self.p.edges[self.edges[start]]),tuple(self.p.edges[self.edges[end]])])
		#Find path
		return np.array(nx.shortest_path(g,self.p.nn,self.p.nn+1))[1:-1]
	def path_loc(self,nstart:int,nend:int,distpath:Optional[float]=None)->NDArray:
		"""
		Computes points' locations on a given path. Locations are computed after mapping each point to the path. Distances in branches away from the path are ignored (set to 0).

		Parameters
		----------
		nstart:		Starting node index to define the path
		nend:		Ending node index to define the path
		distpath:	Distance from nstart to nend node. Computed if missing.

		Returns
		----------
		numpy.ndarray(shape=(len(self)))
			Locations of each point on path as the distance from node nstart. They can be the input of self.p.path_points.
		"""
		if distpath is None:
			distpath=(self.p.topoint()[[nstart]]-self.p.topoint()[[nend]]).ravel()[0]
		loc=self.dist[:,[nstart,nend]].T
		loc=loc[0]-(loc.sum(axis=0)-distpath)/2
		return loc
	def filter_path(self,nstart:int,nend:int)->NDArray:
		"""
		Filters points to retain only those on the given path.

		Parameters
		----------
		nstart:		Starting node index to define the path
		nend:		Ending node index to define the path

		Returns
		----------
		numpy.ndarray
			Indices of points on the path.
		"""
		import numpy as np
		p=self.p.path(nstart,nend)
		edges={self.p.edgedict[p[x],p[x+1]][0] for x in range(len(p)-1)}
		#Point on the edge
		t1=np.array([x in edges for x in self.edges])
		#Point at the terminal node of edge
		t1|=(self.dist[:,p]==0).any(axis=1)
		return np.nonzero(t1)[0]
	def subsets(self,ncell:int,noverlap:int,dmax:float)->Tuple[NDArray,NDArray,NDArray,NDArray,NDArray]:
		"""Constructs overlapping cell subsets as moving window for network reconstruction on each subset. Points in self should be cells.

		Parameters
		----------
		ncell:		int
			Number of cells in each subset
		noverlap:	int
			Number of cell overlap required between neighboring subsets.
		dmax:		float
			Upper bound of distance between neighboring subsets. This bound can be violated if not achievable with data.
		
		Returns
		----------
		edges:		numpy.ndarray(shape=(n_subset))
			Edge of initial center of each cell subset. Center is the anchor point from which the nearest cells are selected.
		locs:		numpy.ndarray(shape=(n_subset))
			Location of initial center of cell subset
		subsets:	numpy.ndarray(shape=(len(self),n_subset),dtype=bool)
			Subset assignment of each cell
		nodegraph:	numpy.ndarray(shape=(n_subset,n_subset),dtype=bool)
			Neighborhood graph between cell subset centers
		subsets:	numpy.ndarray(shape=(n_subset,))
			Subset names

		Method
		-------
		A. Initialize all subsets as node subsets (nearest ncell cells for each node)
		B. For each edge:
			1. Use node subsets as nearest ncell cells for each terminal node of the edge
			2. Initialize current subsets as nearest ncell cells for each cell on the edge. Remove duplicates already in 1.
			3. At a random order, remove each subset in the current subsets if after removal, its nearest neighbor on both side would still satisfy the relation overlap count>noverlap and average distance between cells in the subsets<dmax.
			4. Add remaining current subsets to all subsets

		"""
		import numpy as np
		ns=self.npt
		nn=self.p.nn
		ne=self.p.ne
		assert 0<=noverlap<ncell
		assert dmax>=0

		s=self.copy()
		#Perturb cell locations to avoid identical locations
		subsets=argpartition(s.dist,ncell,axis=0,draw_order='random')[ncell]
		if s.dist[subsets,np.arange(nn)].min()==0:
			raise ValueError('Found trajectory nodes having more cells annotated than parameter ncell. Input finer trajectory or larger ncell.')
		s.perturb()
		subsets=argpartition(s.dist,ncell,axis=0,draw_order='random')[:ncell+1]
		if not np.isfinite(np.max([s.dist[subsets[ncell,x],x] for x in range(nn)])):
			raise ValueError('Cannot find sufficient cells for some trajectory nodes. Input smaller ncell or remove small disconnected subgraphs.')

		#Distances and orders
		distn=s.dist.T
		dist=s-s
		disto,distno=[np.argsort(x,axis=1) for x in [dist,distn]]

		ans_subsets=[]
		ans_neighbors=[]
		ans_edges=[]
		ans_locs=[]
		#A
		ans_subsets+=list(distno[:,:ncell])
		t1=s.p.topoint()
		ans_edges+=list(t1.edges)
		ans_locs+=list(t1.locs)

		#B
		for edge in range(ne):
			#B1,B2
			ids=np.nonzero(s.edges==edge)[0]
			ids=ids[np.argsort(s.locs[ids])]
			subsets=disto[ids][:,:ncell]
			subsets=[ans_subsets[s.p.edges[edge,0]]]+list(subsets)+[ans_subsets[s.p.edges[edge,1]]]
			#Remove identical subsets
			t1=np.nonzero([(subsets[x]!=subsets[x+1]).any() for x in range(1,len(subsets)-1)])[0]
			subsets=[subsets[x] for x in [0]+list(t1)+[len(subsets)-1]]
			ids=ids[t1]
			#Bidirectional linked list
			rot=np.array([np.arange(len(subsets))-1,np.arange(len(subsets))+1]).T
			rot[0,0]=len(subsets)
			#B3
			t0=np.arange(1,len(subsets)-1)
			np.random.shuffle(t0)
			for xi in t0:
				t1=[subsets[x] for x in rot[xi]]
				#Compute number of overlap
				if len(np.intersect1d(*t1,assume_unique=True,return_indices=False))<noverlap:
					continue
				#Compute average distance
				if dist[t1[0]][:,t1[1]].mean()-max(dist[t1[0]][:,t1[0]].mean(),dist[t1[1]][:,t1[1]].mean())>=dmax:
					continue
				#Delete subset
				rot[rot[xi,0],1]=rot[xi,1]
				rot[rot[xi,1],0]=rot[xi,0]
			#B4
			ids2=[0]
			while ids2[-1]<len(subsets)-1:
				ids2.append(rot[ids2[-1],1])
			ids2=np.array(ids2)
			assert ids2[-1]==len(subsets)-1
			ids2=ids2[1:-1]
			ans_edges+=list(s.edges[ids[ids2-1]])
			ans_locs+=list(s.locs[ids[ids2-1]])
			subsets=[subsets[x] for x in ids2]
			if len(subsets)>0:
				ans_neighbors+=[[s.p.edges[edge,0],len(ans_subsets)]]
				ans_neighbors+=[[x,x+1] for x in range(len(ans_subsets),len(ans_subsets)+len(subsets)-1)]
				ans_subsets+=subsets
				ans_neighbors+=[[len(ans_subsets)-1,s.p.edges[edge,1]]]
			else:
				ans_neighbors+=[[s.p.edges[edge,0],s.p.edges[edge,1]]]
		#Validation & output
		n=len(ans_subsets)
		assert len(ans_neighbors)==n-1 and len(ans_edges)==n and len(ans_locs)==n
		ans_edges,ans_locs,ans_subsets,ans_neighbors=[np.array(x) for x in [ans_edges,ans_locs,ans_subsets,ans_neighbors]]
		assert (ans_neighbors>=0).all() and (ans_neighbors<n).all()
		ans_edges=ans_edges.astype('u2')
		t1=np.zeros((n,n),dtype=bool)
		t1[ans_neighbors[:,0],ans_neighbors[:,1]]=True
		t1|=t1.T
		ans_neighbors=t1
		t1=np.zeros((ns,n),dtype=bool)
		for xi in range(n):
			t1[ans_subsets[xi],xi]=True
		ans_subsets=t1
		ansn_subsets=np.array([f'Subset{x+1}' for x in range(n)])
		return (ans_edges,ans_locs,ans_subsets,ans_neighbors,ansn_subsets)
	def nearest(self,w:Optional[NDArray]=None)->point:
		"""
		Find the point on the trajectory with shortest average distance to all current points.
		Parameters
		----------
		w:	Weight for each of current point. Can have at most 2 dimensions. The last dimension must match current points. Defaults to uniform (all 1).

		Returns
		-------
		Points with shortest average distance. One for each weight vector.

		"""
		import numpy as np
		from dictys.utils.numpy import median
		if w is None:
			w=np.ones(len(self))
		assert w.shape[-1]==len(self)
		if w.ndim==1:
			w=w.reshape(1,-1)
		elif w.ndim>2:
			raise ValueError('w must have at most 2 dimensions.')

		dist=self.dist.T@w.T
		#Find nearest node
		nodes=dist.argmin(axis=0)
		ans=[]
		for xi in range(len(nodes)):
			#Check all edges of the nearest node
			edges=np.nonzero((self.p.edges==nodes[xi]).any(axis=1))[0]
			ans1=[]
			for edge in edges:
				node=self.p.edges[edge]
				#Points on edge
				p1=self.edges==edge
				#Points on either side of edge
				p2=[((dist[node[0],xi]-dist[node[1],xi])*x<0)&~p1 for x in [1,-1]]
				#Convert to weighted median problem
				v=np.r_[0,self.locs[p1],self.p.lens[edge]]
				w2=np.r_[w[xi,p2[0]].sum(),w[xi,p1],w[xi,p2[1]].sum()]
				v=median(v,w=w2)
				ans1.append(self.__class__(self.p,np.array([edge]),np.array([v])))
			#Find minimum across all edges
			t1=((self.concat(ans1)-self)@w[xi])
			t1=t1.argmin()
			ans.append(ans1[t1])
		ans=self.concat(ans)
		return ans
	def subtraj(self,edges:Optional[NDArray]=None)->trajectory:
		"""
		Convert list of points to a finer trajectory by treating each point as a trajectory node.
		Each terminal and branching node of trajectory must have one corresponding point.

		Parameters
		----------
		edges:	numpy.ndarray(shape=(len(self),2))
			The edge each point is on, specified by the start and end nodes. Must be specified.

		Returns
		----------
		dictys.traj.trajectory
			Trajectory whose each node is each point in self.
		"""
		if edges is None:
			raise NotImplementedError
		assert edges.ndim==2 and edges.shape[1]==2
		#Check existence of points for terminal & branching nodes
		# assert (self.dist[:,self.p.deg!=2]==0).any(axis=0).all()
		dist=self-self
		return trajectory(edges,dist[edges[:,0],edges[:,1]])
	def __getitem__(self,key):
		"""
		Choose a subset of points.
		"""
		if not isinstance(key,slice) and not hasattr(key,'__len__'):
			raise TypeError('Key must be iterable with __len__.')
		return self.__class__(self.p,self.edges[key],self.locs[key],dist=self.dist[key])
	def __len__(self)->int:
		return self.npt
	def drop_duplicates(self,selection:bool=False)->Union[point,Tuple[point,NDArray]]:
		"""
		Remove duplicate points.
		
		Parameters
		----------
		selection:
			Whether to return which points are selected

		Returns
		----------
		pts_new:	point
			New object with random point order
		selection:	numpy,ndarray(dtype=int)
			Indices of points that remained. Only included if selection=True.
		"""
		import numpy as np
		from dictys.utils.numpy import groupby
		t1=groupby(self.edges).values()
		t1=np.concatenate([x[np.unique(self.locs[x],return_index=True)[1]] for x in t1])
		if selection:
			return (self[t1],t1)
		return self[t1]
	def weight_linear(self,other:point)->NDArray:
		"""
		Smoothing function to compute data on other points with data on current (self's) points with linear interpolation between points.

		Parameters
		----------
		other: dictys.traj.point
			Points to compute data for

		Returns
		----------
		numpy.ndarray(shape=[len(self),len(other)])
			Weight of self's points on other's points.
		"""
		import numpy as np
		n=len(other)
		d=self-other
		d2=self-self
		#Nearest point
		t1=d.argmin(axis=0)
		#Nearest point but on other side of nearest point
		#Other side (neccesary condition)
		t2=d2[t1]>d.T
		t2=[np.nonzero(x)[0] for x in t2]
		#Nearest point
		t2=np.array([t2[x][d[t2[x],x].argmin()] if len(t2[x])>0 else -1 for x in range(n)])
		assert t2.shape==t1.shape
		t1=np.array([t1,t2])
		assert (t1[0]!=t1[1]).all()
		#Linear formula by distance
		t3=(t1>=0).all(axis=0)
		t2=np.array([d[t1[1],np.arange(n)],d[t1[0],np.arange(n)]])
		t2=(t2[0]+1E-300)/(t2.sum(axis=0)+2E-300)
		w=np.zeros((len(self),n),dtype=float)
		w[t1[0],np.arange(n)]=t2*t3
		w[t1[1],np.arange(n)]=(1-t2)*t3
		#Single value assignment if on nearest point
		w[t1[0],np.arange(n)]+=t1[1]<0
		return w
	def weight_conv(self,other:point,radius:float,cut:float=0)->NDArray:
		"""
		Smoothing function to compute data on other points with data on current (self's) points with Gaussian kernel smoothing.

		Parameters
		----------
		other:	dictys.traj.point
			Other points to compute weight for
		radius:	float
			Radius or sigma of gaussian filter as distance
		cut:	float
			Set weight to 0 if below this threshold

		Returns
		----------
		numpy.ndarray(shape=[len(self),len(other)])
			Weight of each node on each point.
		"""
		import numpy as np
		d=self-other
		w=np.exp(-((d/radius)**2)/2)
		w/=w.sum(axis=0)
		if cut>0:
			raise NotImplementedError('cut>0')
			#TODO: Need to remove only the weakest if all weights removed for each point
			# t1=(w>0)&(w<cut)
			# while t1.any():
			# 	w[t1]=0
			# 	w/=w.sum(axis=0)
			# 	t1=(w>0)&(w<cut)
		assert (w>=0).all() and (w<=1).all()
		return w
	def smoothen(self,data:NDArray,*a,points:point=None,axis:int=-1,func_name:str='linear',nan:str='ignore',**ka)->NDArray:
		"""
		Function to compute smoothened/interpolated values of other points based on data at current points. Used by self.smoothened.

		Parameters
		----------
		data:		numpy.ndarray(shape=(...,len(self),...)
			Data to smoothen/interpolate from.
		a:			list
			Arguments of smoothening function.
		points:		dictys.traj.point
			Points to interpolate data for. Required.
		axis:		int
			Axis of data to smoothen. data must have length len(self) at this axis.
		func_name:	str
			Name of smoothening function of class point. Using self.weight_conv should provide func_name='conv'. Accepts:
			linear:	For linear interpolation
			conv: For Gaussian kernel smoothing
		nan:		str
			How to handle nan in data. Accepts
			ignore:		Weights becomes zero, so output is nan only if all data with nonzero weights are nan.
			propagate:	Weights are kept, so nan propagates
		ka:			dict
			Keyword arguments smoothening function.

		Returns
		----------
		numpy.ndarray(shape=(...,len(points),...)
			Smoothened values that have the same dimensionality as data in all axes other than axis.
		"""
		import numpy as np
		assert points is not None
		assert data.shape[axis]==self.npt
		func=getattr(self,'weight_'+func_name)
		w=func(points,*a,**ka)
		assert w.ndim==2
		if nan=='ignore':
			w2=~np.isnan(data)
			data=np.nan_to_num(data)
			ans=data.swapaxes(axis,-1)@w
			w2=w2.swapaxes(axis,-1)@w
			ans=ans/(w2+1E-300)
			ans[w2==0]=np.nan
			ans=ans.swapaxes(axis,-1)
		elif nan=='propagate':
			ans=(data.swapaxes(axis,-1)@w).swapaxes(axis,-1)
		else:
			raise ValueError('Unknown nan parameter value.')
		assert ans.shape[:axis]==data.shape[:axis]
		assert ans.shape[axis+1]==data.shape[axis+1]
		return ans
	def smoothened(self,data:NDArray,*a,**ka)->Callable[point,NDArray]:
		"""
		Creates a smoothened/interpolated function on data of self's points that computes values at other points

		Parameters
		----------
		data:	numpy.ndarray(shape=(...,len(self),...)
			Base data to smoothen/interpolate from.
		a:		list
			Other arguments for self.smoothen.
		ka:		dict
			Keyword arguments for self.smoothen, except `points`.

		Returns
		----------
		function(point)->numpy.ndarray(shape=(...,len(point),...))
			Function to compute smoothened data matrix at any points. It has the same dimensionality as data in all axes other than axis.
		"""
		from functools import partial
		return partial(self.smoothen,data,*a,**ka)
	#I/O
	@classmethod
	def from_fileobj(cls,traj:trajectory,f:Union[h5py.File,h5py.Group])->point:
		"""
		Load object from file object

		Parameters
		----------
		traj:	dictys.traj.trajectory
			Trajectory object for points
		f:	fileobj
			File object to load from

		Returns
		----------
		dictys.traj.point
			Loaded class object
		"""
		import numpy as np
		params=[np.array(f['edges']),np.array(f['locs'])]
		return cls(traj,*params)
	@classmethod
	def from_file(cls,traj:Optional[trajectory],path:str)->point:
		"""
		Load object from file

		Parameters
		----------
		traj:	dictys.traj.trajectory
			Trajectory object for points. Set to None to load and use trajectory within the same file.
		path:	str
			File path to load from

		Returns
		----------
		dictys.traj.point
			Loaded class object
		"""
		import logging
		logging.info(f'Reading file {path}.')
		with h5py.File(path,'r') as f:
			if traj is None:
				traj=trajectory.from_fileobj(f['traj'])
			return cls.from_fileobj(traj,f)
	def to_fileobj(self,f:Union[h5py.File,h5py.Group],compression:str="gzip",**ka)->None:
		"""
		Save object to file object

		Parameters
		----------
		f:	fileobj
			File object to save to
		compression: str
			Type of compression. See h5py.File.create_dataset.
		ka:	dict
			Keyword arguments passed to h5py.File.create_dataset
		"""
		p=dict(ka)
		p['compression']=compression
		p['data']=self.edges
		f.create_dataset('edges',**p)
		p=dict(ka)
		p['compression']=compression
		p['data']=self.locs
		f.create_dataset('locs',**p)
	def to_file(self,path:str,traj:bool=False,**ka)->None:
		"""
		Save object to file

		Parameters
		----------
		path:		str
			File path to save to.
		traj:		bool
			Whether to save trajectory in the same file
		ka:			dict
			Keyword arguments passed to self.to_fileobj
		"""
		import logging
		logging.info(f'Writing file {path}.')
		with h5py.File(path,'w') as f:
			if traj:
				self.p.to_fileobj(f.create_group('traj'))
			return self.to_fileobj(f,**ka)












































#
