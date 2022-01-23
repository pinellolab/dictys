
#!/usr/bin/python3

_docstring2argparse_ignore_=['trajectory','point']

class trajectory:
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
		return point.fromnodes(self)
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
		return point.fromnodes(self)
	def linspace(self,start,end,n):
		"""Find evenly spaced points on a path np.linspace
		start,
		end:	Start/end nodes' IDs to indicate the path.
		n:		Number of points including terminal nodes
		Return:
		Instance of point class with points go from start to end nodes.
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
		Instance of point class with points go from start to end nodes.
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
		return point(self,ans_steps,ans_lengths)
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
		func_name:	Name of smoothening function of class point. Using weight_conv should provide func_name='conv',
		a:			Arguments of smoothening function.
		ka:			Keyword arguments smoothening function.
		points:		point instance to interpolate data for
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
		func_name:	Name of smoothening function of class point. Using weight_conv should provide func_name='conv',
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
	#I/O
	@classmethod
	def from_fileobj(cls,f):
		import numpy as np
		params=[np.array(f['edges']),np.array(f['lens'])]
		return cls(*params)
	@classmethod
	def from_file(cls,path):
		import h5py
		if isinstance(path,h5py.File) or isinstance(path,h5py.Group):
			return cls.from_fileobj(path)
		with h5py.File(path,'r') as f:
			return cls.from_fileobj(f)
	def to_fileobj(self,f,compression="gzip",**ka):
		p=dict(ka)
		p['compression']=compression
		p['data']=self.lens
		f.create_dataset('lens',**p)
		p=dict(ka)
		p['compression']=compression
		p['data']=self.edges
		f.create_dataset('edges',**p)
	def to_file(self,path,**ka):
		import h5py
		if isinstance(path,h5py.File) or isinstance(path,h5py.Group):
			return self.to_fileobj(path,**ka)
		with h5py.File(path,'w') as f:
			return self.to_fileobj(f,**ka)

class point:
	def __init__(self,traj,edges,locs,dist=None):
		"""Point list on trajectory
		traj:	Trajectory the points are on. Instance of trajectory
		edges:	Edge ID each point is on as np.array(shape=(n_point))
		losc:	Distance from the starting node on its edge for each point as np.array(shape=(n_point))
		dist:	Distance matrix between each point and each node as np.array(shape=(n_point,n_node)).
				Automatically computed if not provided
		"""
		import numpy as np
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
	def from_dist(cls,traj,edges,dist):
		"""Compute edges and locs given dist"""
		import numpy as np
		assert edges.ndim==1
		assert dist.shape==(len(edges),traj.nn)
		locs=dist[np.arange(len(edges)),traj.edges[edges,0]]
		return cls(traj,edges,locs,dist=dist)
	@classmethod
	def fromnodes(cls,traj):
		"""Create list of points for nodes of the trajectory"""
		import numpy as np
		nodes,indices=np.unique(traj.edges,return_index=True)
		indices=[indices//2,indices%2]
		edges=np.zeros(traj.nn,dtype=int)
		edges[nodes]=indices[0]
		locs=np.zeros(traj.nn,dtype=traj.lens.dtype)
		locs[nodes]=indices[1]*traj.lens[edges]
		return cls(traj,edges,locs,dist=traj.dist)
	def copy(self):
		return self.__class__(self.p,self.edges,self.locs,dist=self.dist)
	def __sub__(self,other):
		"""Subtraction computes the distance matrix between all point pairs in two point lists.
		other:	Destination point list
		Return:
		Distance matrix as np.array(shape=(self.npt,other.npt))
		"""
		import numpy as np
		from lwang.utils.numpy import groupby
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
	def compute_dist(self):
		"""Compute distance to every node on trajectory given locs
		Return:
		Distance matrix as np.array(shape=(n_point,n_node))
		"""
		return self-self.fromnodes(self.p)
	def perturb(self,scale=1):
		"""Perturb locations of points without changing order with other nodes in the list. Overlapping points will have random order.
		"""
		import numpy as np
		from collections import defaultdict
		import itertools
		import networkx as nx
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
			t3=list(set([frozenset(x[1]+[x[0]]) for x in t3.items()]))
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
	def path_loc(self,nstart,nend,distpath=None):
		"""Computes locations on the path from node nstart to nend.
		Distances due to deviation from the main path are removed.
		nstart,
		nend:		Starting/ending node ID of path
		distpath:	Distance from nstart to nend node. Computed if missing
		Return:
		Locations on path as np.array that can be the input of self.p.path_points.
		"""
		if distpath is None:
			distpath=(self.p.topoint()[[nstart]]-self.p.topoint()[[nend]]).ravel()[0]
		loc=self.dist[:,[nstart,nend]].T
		loc=loc[0]-(loc.sum(axis=0)-distpath)/2
		return loc
	def filter_path(self,nstart,nend):
		"""Filters points to retain only those on the given path from node nstart to node nend.
		nstart,
		nend:	IDs of starting and ending nodes that define the path
		Return:	np.array of index of points on the path."""
		import numpy as np
		p=self.p.path(nstart,nend)
		edges=set([self.p.edgedict[p[x],p[x+1]][0] for x in range(len(p)-1)])
		#Point on the edge
		t1=np.array([x in edges for x in self.edges])
		#Point at the terminal node of edge
		t1|=(self.dist[:,p]==0).any(axis=1)
		return np.nonzero(t1)[0]
	def subsets(self,ncell,noverlap):
		"""Construct overlapping cell subsets for network reconstruction on each subset.
		The subset assignment is based on the following hard constraints:
			a. Number of cells in each subset from ncell
			b. Upper bound of cell overlap count between neighboring subsets from noverlap 
			c. Every cell is assigned to at least one subset.
		There are exceptions where it's impossible to satisfy all constraints simultaneously so b will be relaxed and a warning will appear. Such scenarios include:
			a. The subsets of two neighboring branching/terminal nodes have overlap cell count>noverlap
			b. Between two neighboring branching/terminal nodes the number of unassigned cells<ncell-2*noverlap
			c. Very small noverlap.
		For cells at identical (pseudotime) locations, to choose a subset random draws are performed whenever necessary.
		Error will be triggered if there are more distance=0 cells than ncell at any node.
		
		Parameters:
		dist:     Distance matrix between each cell and each node as np.array(shape=(n_cell,n_node))
		edges:    Start and end node IDs for each edge as np.array(shape=(n_edge,2))
		branch:   Cell edge membership by ID as np.array(shape=n_cell)
		ncell:    Number of cells in each subset
		noverlap: Soft upper bound on the number of overlap cells between neighboring subsets.
		
		Steps:
		A. Construct subset at each node
			1. Create a subset containing cells with distance=0. Raise error if there are more distance=0 cells than ncell.
			2. Add the nearest cells to subset to reach ncell cells.
		B. Construct subsets for each edge between nodes 
			1. Determine the number of intermediate subsets as floor(( the number of cells not in any subset + noverlap ) / (ncell-noverlap))
			2. If the number of intermediate subsets=0 and the number of unassigned cells>0, use one intermediate subset instead.
				a. If one intermediate subset, set subset central location to those of unassigned cells. Proceed to 4d.
			3. Determine noverlap_edge as noverlap for this edge with equation in 1.
			4. Assign subsets with equal overlap between each neighboring subset pair. Takes these steps:
				a. Order unassigned cells and noverlap_edge/2 cells from either side of edge on 1 dimension
				b. Split cells evenly and continuously into 2*(number of intermediate subsets) groups with linspace
				c. Middle point between indicies [0,2], [2,4], etc are the central locations of subset [0,1,...]
				d. Assign the nearest ncell cells for each central location with function distance.
		
		Return:
		Data:
		edges:   	Edge of center of each cell subset np.array(shape=(n_subset)). Center is mean of min and max.
		locs:    	Distance of center of each cell subset from starting node as np.array(shape=(n_subset))
		radius:  	Distance to center of cell subset to be included as np.array(shape=(n_subset))
		subsets: 	Subset assignment of each cell as np.array(shape=(n_cell,n_subset))
		nodegraph:	Neighborhood graph between nodes as np.array(shape=(n_subset,n_subset),dtype=bool)
		Dimenions:
		subsets: 	Subset names as np.array(shape=(n_subset,))
		"""
		import numpy as np
		import logging
		from collections import Counter
		ns=self.npt
		nn=self.p.nn
		ne=self.p.ne
		assert noverlap<ncell and noverlap>=0
		
		ans_edges=[]
		ans_locs=[]
		ans_radius=[]
		test_mask=np.zeros(ns,dtype=bool)
		s=self.copy()

		#Step A
		subsets=argpartition(s.dist,ncell,axis=0,draw_order='random')[ncell]
		if s.dist[subsets,np.arange(nn)].min()==0:
			raise ValueError('Found trajectory nodes having more cells annotated than parameter ncell. Input finer trajectory or larger ncell.')
		s.perturb()
		subsets=argpartition(s.dist,ncell,axis=0,draw_order='random')[:ncell+1]
		if not np.isfinite(np.max([s.dist[subsets[ncell,x],x] for x in range(nn)])):
			raise ValueError('Cannot find sufficient cells for some trajectory nodes. Input smaller ncell or remove small disconnected subgraphs.')
		subsets=subsets[:ncell].T
		t1=np.zeros((ns,nn),dtype=bool)
		for xi in range(nn):
			t1[subsets[xi],xi]=True
			t2=[np.nonzero(x)[0] for x in s.p.edges.T==xi]
			if len(t2[0])>0:
				ans_edges.append(t2[0][0])
				ans_locs.append(0)
			else:
				assert len(t2[1])>0
				t2=t2[1][0]
				ans_edges.append(t2)
				ans_locs.append(s.p.lens[t2])
			ans_radius.append(distance(ans_edges[-1],ans_locs[-1],s.dist[subsets[xi]],s.p.edges,s.edges[subsets[xi]]).max())
		subsets=t1
		subsets_extra=[]
		test_mask[subsets.any(axis=1)]=True
		nodeg=[]
		
		for xi in range(ne):
			#Step B1,B2
			nodes=s.p.edges[xi]
			nodestart,nodeend=nodes
			freecell=np.nonzero(s.edges==xi)[0]
			freecell=freecell[~subsets[freecell].any(axis=1)]
			if len(freecell)==0:
				nodeg.append(nodes)
				continue
			nsubset=int(np.floor((len(freecell)+noverlap)/(ncell-noverlap)))
			if ncell*nsubset<len(freecell):
				nsubset+=1
			if nsubset==1:
				#Step B2a
				cells=freecell[s.dist[freecell,nodestart].argsort()]
				centrals=np.array([s.dist[cells[[0,-1]],nodestart].mean()])
				test_mask[cells]=True
			else:
				#Step B3
				noverlap_edge_half=int(np.floor((ncell*nsubset-len(freecell))/(2*(nsubset+1))))
				assert noverlap_edge_half>=0
				#Step B4a
				#Order free cells based on distance from starting node
				cells=freecell[s.dist[freecell,nodestart].argsort()]
				t1=s.dist[cells,nodestart]
				t1=t1*(2*((s.edges[cells]==xi)|(s.dist[cells,nodeend]<s.dist[cells,nodestart]))-1)
				assert (t1[1:]>=t1[:-1]).all()
				
				t0=[]
				if len(t0)<noverlap_edge_half:
					#Include half overlapping cells on starting node side on the same edge
					t1=np.nonzero(s.edges==xi)[0]
					t1=t1[s.dist[t1,nodestart]<=s.dist[cells[0],nodestart]]
					t2=set(cells)
					t0=t1[[x not in t2 for x in t1]]
					if len(t0)>0:
						if len(t0)>noverlap_edge_half:
							t0=t0[argpartition(s.dist[t0,nodestart],-noverlap_edge_half,draw_order='error')[-noverlap_edge_half:]]
						t0=t0[s.dist[t0,nodestart].argsort()]
				if len(t0)<noverlap_edge_half:
					#Include half overlapping cells on starting node side on different edges
					t1=np.nonzero(s.edges!=xi)[0]
					t1=t1[s.dist[t1,nodestart]<s.dist[t1,nodeend]]
					if len(t1)+len(t0)>noverlap_edge_half:
						t1=t1[argpartition(s.dist[t1,nodestart],noverlap_edge_half-len(t0),draw_order='error')[:noverlap_edge_half-len(t0)]]
					t0=np.r_[t1[s.dist[t1,nodestart].argsort()[::-1]],t0]
				if len(t0)<noverlap_edge_half:
					logging.warning('Cannot reach desired overlap rate at some terminal nodes. Please check validity of cell subsetting output or reduce parameter noverlap.')
				else:
					assert len(t0)==noverlap_edge_half
				cells=np.r_[t0,cells].astype(int)
				t1=s.dist[cells,nodestart]
				t1=t1*(2*((s.edges[cells]==xi)|(s.dist[cells,nodeend]<s.dist[cells,nodestart]))-1)
				assert (t1[1:]-t1[:-1]>=0).all()

				t0=[]
				if len(t0)<noverlap_edge_half:
					#Include half overlapping cells on ending node side on the same edge
					t1=np.nonzero(s.edges==xi)[0]
					t1=t1[s.dist[t1,nodeend]<=s.dist[cells[-1],nodeend]]
					t2=set(cells)
					t0=t1[[x not in t2 for x in t1]]
					if len(t0)>0:
						if len(t0)>noverlap_edge_half:
							t0=t0[argpartition(s.dist[t0,nodeend],-noverlap_edge_half,draw_order='error')[-noverlap_edge_half:]]
						t0=t0[s.dist[t0,nodeend].argsort()[::-1]]
				if len(t0)<noverlap_edge_half:
					#Include half overlapping cells on ending node side on different edges
					t1=np.nonzero(s.edges!=xi)[0]
					t1=t1[s.dist[t1,nodeend]<s.dist[t1,nodestart]]
					t1=t1[[x not in t2 for x in t1]]
					if len(t1)+len(t0)>noverlap_edge_half:
						t1=t1[argpartition(s.dist[t1,nodeend],noverlap_edge_half-len(t0),draw_order='error')[:noverlap_edge_half-len(t0)]]
					t0=np.r_[t0,t1[s.dist[t1,nodeend].argsort()]]
				if len(t0)<noverlap_edge_half:
					logging.warning('Cannot reach desired overlap rate at some terminal nodes. Please check validity of cell subsetting output or reduce parameter noverlap.')
				else:
					assert len(t0)==noverlap_edge_half
				cells=np.r_[cells,t0].astype(int)
				t1=s.dist[cells,nodestart]
				t1=t1*(2*((s.edges[cells]==xi)|(s.dist[cells,nodeend]<s.dist[cells,nodestart]))-1)
				assert (t1[1:]>=t1[:-1]).all()
				test_mask[cells]=True
				
				#Step B4bc
				centrals=(np.linspace(0,1,2*nsubset+1)[::2]*len(cells)).astype(int)
				centrals=np.concatenate([[centrals[x],centrals[x+1]-1] for x in range(len(centrals)-1)])
				centrals=cells[centrals]
				centrals=[
					s.dist[centrals,nodestart]*(2*((s.edges[centrals]==xi)|(s.dist[centrals,nodestart]>s.dist[centrals,nodeend]))-1),
					s.dist[centrals,nodeend]*(2*((s.edges[centrals]==xi)|(s.dist[centrals,nodestart]<s.dist[centrals,nodeend]))-1),
				]
				centrals=np.array(centrals).reshape(2,len(centrals[0])//2,2)
				centrals=centrals.mean(axis=2).T
				t1=(centrals>0).all(axis=1)
				if (~t1).any():
					centrals=centrals[t1]
				assert (centrals>0).all()
				centrals=centrals[:,0]
			#Step B4d
			t3=distance(xi,centrals,s.dist,s.p.edges,s.edges)
			t2=argpartition(t3,ncell,axis=0,draw_order='error')[:ncell].T
			t1=np.zeros((ns,len(centrals)),dtype=bool)
			for xj in range(len(centrals)):
				t1[t2[xj],xj]=True
			subsets_extra.append(t1)
			nodeg+=[[nodestart,len(ans_radius)]]+[[len(ans_radius)+x,len(ans_radius)+x+1] for x in range(len(centrals)-1)]+[[len(ans_radius)+len(centrals)-1,nodeend]]
			ans_edges+=[xi]*len(centrals)
			ans_locs+=list(centrals)
			ans_radius+=list(t3[t2[:,ncell-1],np.arange(len(centrals))])

		ans_subsets=np.concatenate([subsets]+subsets_extra,axis=1)
		ans_edges,ans_locs,ans_radius=[np.array(x) for x in [ans_edges,ans_locs,ans_radius]]
		ans_edges=ans_edges.astype('u2')
		ans_locs=self.p.conform_locs(ans_locs,ans_edges)
		ansname_subsets=np.array([f'subset{x}' for x in range(ans_subsets.shape[1])])
		nodeg=np.array(nodeg).T
		ans_nodegraph=np.zeros((ans_subsets.shape[1],ans_subsets.shape[1]),dtype=bool)
		ans_nodegraph[nodeg[0],nodeg[1]]=True
		ans_nodegraph[nodeg[1],nodeg[0]]=True
		assert ans_subsets.ndim==2
		assert all([x.shape==(ans_subsets.shape[1],) for x in [ans_edges,ans_locs,ans_radius]])
		assert ans_subsets.any(axis=1).all()
		assert (ans_subsets.sum(axis=0)==ncell).all()
		assert (ans_nodegraph.T==ans_nodegraph).all()
		assert (ans_nodegraph.sum(axis=0)>=1).all()
		assert not np.diag(ans_nodegraph).any()
		return ans_edges,ans_locs,ans_radius,ans_subsets,ans_nodegraph,ansname_subsets
	def subsets2(self,ncell,noverlap,dmax):
		"""Construct overlapping cell subsets for network reconstruction on each subset.
		Parameters
		ncell:		Number of cells in each subset
		noverlap:	Average number of cell overlap between neighboring subsets.
		dmax:		Upper bound of distance between neighboring subsets. Violation exception: not possible with data.

		A. Initialize all subsets as node subsets (nearest ncell cells for each node)
		B. For each edge:
			1. Use node subsets as nearest ncell cells for each terminal node of the edge
			2. Initial current subsets as nearest ncell cells for each cell on the edge. Remove duplicates already in 1.
			3. At a random order, remove each subset in the current subsets if after removal, its nearest neighbor on both side would still satisfy the relation overlap count>noverlap and average distance between cells in the subsets<dmax.
			4. Add remaining current subsets to all subsets

		Return:
		Data:
		edges:   	Edge of initial center of each cell subset np.array(shape=(n_subset)). Center is the anchor point from which the nearest cells are selected.
		locs:    	Location of initial center of cell subset as np.array(shape=(n_subset))
		subsets: 	Subset assignment of each cell as np.array(shape=(n_cell,n_subset))
		nodegraph:	Neighborhood graph between nodes as np.array(shape=(n_subset,n_subset),dtype=bool)
		Dimenions:
		subsets: 	Subset names as np.array(shape=(n_subset,))
		"""
		import numpy as np
		import logging
		from collections import Counter
		import itertools
		ns=self.npt
		nn=self.p.nn
		ne=self.p.ne
		assert noverlap<ncell and noverlap>=0
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
		return [ans_edges,ans_locs,ans_subsets,ans_neighbors,ansn_subsets]
	def subtraj(self,edges=None):
		"""Convert list of points to a finer trajectory by treating each point as a node.
		Each terminal and branching node of trajectory must have one corresponding point.
		edges:	Optional parameter to specify edges between points
		Return:
		Trajectory as instance of trajectory where each node is each point in self.
		"""
		if edges is None:
			raise NotImplementedError
		assert edges.ndim==2 and edges.shape[1]==2
		import numpy as np
		#Check existence of points for terminal & branching nodes
		# assert (self.dist[:,self.p.deg!=2]==0).any(axis=0).all()
		dist=self-self
		return trajectory(edges,dist[edges[:,0],edges[:,1]])
	def __getitem__(self,key):
		if not isinstance(key,slice) and not hasattr(key,'__len__'):
			raise TypeError('Key must be iterable with __len__.')
		return self.__class__(self.p,self.edges[key],self.locs[key],dist=self.dist[key])
	def __len__(self):
		return self.npt
	def weight_linear(self):
		"""Treat properties of each point as a weighted sum of nodes.
		Returns weight matrix between points and nodes.
		Uses linear interpolation between nodes.
		Return:
		Weight of each node on each point as np.array(shape=[npt,nn]).
		Sums to 1 over axis 1
		"""
		import numpy as np
		w=self.locs/self.p.lens[self.edges]
		assert (w>=0).all() and (w<=1).all()
		ans_w=np.zeros((self.npt,self.p.nn),dtype=float)
		ans_w[np.arange(self.npt),self.p.edges[self.edges,0]]=1-w
		ans_w[np.arange(self.npt),self.p.edges[self.edges,1]]=w
		return ans_w
	def weight_conv(self,radius,cut=0,nodes=None,nodes_path=None):
		"""Treat properties of each point as a weighted sum of nodes.
		Returns weight matrix between points and nodes.
		Uses convolution based on distance between each point and each node.
		radius:	Radius or sigma of gaussian filter in terms of distance
		cut:	Set node weight to 0 if below this threshold
		nodes:	Set of nodes to use for convolution. Defaults to all.
		nodes_path: [start,end] node IDs whose shortest path defines the nodes to use for convolution. Disabled by default.
		Return:
		Weight of each node on each point as np.array(shape=[npt,nn]).
		Sums to 1 over axis 1
		"""
		import numpy as np
		w=np.exp(-((self.dist/radius)**2)/2)
		if nodes_path is not None:
			assert nodes is None
			assert len(nodes_path)==2
			nodes=self.p.path(*nodes_path)
		if nodes is not None:
			assert np.min(nodes)>=0 and np.max(nodes)<self.p.nn
			nodes=set(nodes)
			nodes=np.array([x in nodes for x in range(self.p.nn)])
			w[:,~nodes]=0
		w=(w.T/w.sum(axis=1)).T
		if cut>0:
			t1=(w>0)&(w<cut)
			while t1.any():
				w[t1]=0
				w=(w.T/w.sum(axis=1)).T
				t1=(w>0)&(w<cut)
		assert (w>=0).all() and (w<=1).all()
		return w












































#
