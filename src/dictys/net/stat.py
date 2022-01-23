#!/usr/bin/python3
# Lingfei Wang, 2020-2022. All rights reserved.

"""
Statistics of networks
"""
import dictys.traj

class stat_base:
	def __init__(self,names=None,label=None):
		"""Base class for statistics for each gene
		label:	Label of stat that is shown as coordindate label.
		names:	List of names of each axis of output stat, except last axis which is always time. Default is obtained from default_names function.
		"""
		if label is None:
			label=self.default_label()
		self.label=label
		if names is None:
			names=self.default_names()
		assert len(names)>0
		self.names=names
		self.ndict=[dict(zip(x,range(len(x)))) for x in names]
	def compute(self,pts):
		"""Use this function to compute stat values at each state or point
		pts:	Point list instance of lwang.pipeline.dynet2.traj.pointc, or state list either by name or by id as list of str or list of int
		Return:
		Stat value for each gene as np.array(shape=(...,len(pts))) . Use nan to hide value or set as invalid.
		"""
		raise NotImplementedError
	def default_names(self):
		"""Use this function to determine the default names for each axis.
		Only names shared with other stats will be visualized
		Return:
		List of list of names for each axis.
		"""
		raise NotImplementedError
	def default_lims(self,pts=None,names=None,expansion=0.02):
		"""Use this function to determine the default limits of the stat.
		This implementation uses min/max of stat values.
		pts:	Point list instance to compute min/max
		namet:	Gene names used to compute limits
		expansion:	Expand limits by this relative amount on each side.
		Return:
		Limits as [min,max]
		"""
		if pts is not None:
			raise NotImplementedError
		ans=self.compute(pts)
		if names is not None:
			assert len(names)==len(self.names)
			assert all([all([y in self.ndict[x] for y in names[x]]) for x in range(len(self.names))])
			for xi in range(len(names)):
				ans=ans.swapaxes(0,xi)[[self.ndict[xi][x] for x in names[xi]]].swapaxes(0,xi)
		ans=[ans.min(),ans.max()]
		t1=(ans[1]-ans[0])*expansion
		return [ans[0]-t1,ans[1]+t1]
	def default_label(self):
		"""Use this function to determine the label of this stat
		Return:
		Label as str
		"""	
		raise NotImplementedError
	def __add__(self,other):
		from operator import add
		if isinstance(other,stat_base):
			return statf_function(add,[self,other],label='({})+({})'.format(self.label,other.label))
		else:
			raise NotImplementedError
	def __sub__(self,other):
		from operator import sub
		if isinstance(other,stat_base):
			return statf_function(sub,[self,other],label='({})-({})'.format(self.label,other.label))
		else:
			raise NotImplementedError
	def __mul__(self,other):
		from operator import truediv
		if isinstance(other,stat_base):
			return statf_function(mul,[self,other],label='({})*({})'.format(self.label,other.label))
		else:
			raise NotImplementedError
	def __truediv__(self,other):
		from operator import truediv
		if isinstance(other,stat_base):
			return statf_function(truediv,[self,other],label='({})/({})'.format(self.label,other.label))
		else:
			raise NotImplementedError

class statf_function(stat_base):
	def __init__(self,func,stats,**ka):
		"""Stat that is any function of any other stat(s)
		func:	Function to combine other stats. Should have self.compute_points=func(*[x.compute_points(...) for x in stats]).
		stats:	List of stats whose final outputs (compute_points) will be operated on by func.
		label:	Label of stat that is shown as coordindate label.
		"""
		self.n=len(stats)
		assert self.n>0
		self.func=func
		self.stats=stats
		super().__init__(**ka)
	def default_names(self):
		"""Use names shared by all stats by default. Only suitable when all stats have the same number of dimensions
		"""
		from functools import reduce
		from operator import and_
		if any([len(x.names)!=len(self.stats[0].names) for x in self.stats[1:]]):
			raise NotImplementedError
		names=[[set(z) for z in y] for y in zip(*[x.names for x in self.stats])]
		names=[sorted(list(reduce(and_,x))) for x in names]
		return names
	def default_label(self):
		return 'Function'
	def compute(self,pts):
		"""Computes stat for states or points
		ans:	Result of computation by given stats
		Return:
		stat as np.array(shape=[len(x) for x in self.names]+[n])
		"""
		import numpy as np
		n=len(pts)
		ans=[x.compute(pts) for x in self.stats]
		assert all([ans[x].ndim==len(self.stats[x].names)+1 for x in range(self.n)])
		assert all([ans[x].shape==tuple([len(y) for y in self.stats[x].names]+[n]) for x in range(self.n)])
		t1=[[[x.ndict[y][z] for z in self.names[y]] for y in range(len(self.names))] for x in self.stats]
		for xi in range(len(self.names)):
			ans=[ans[x].swapaxes(0,xi)[t1[x][xi]].swapaxes(0,xi) for x in range(self.n)]
		assert all([x.shape==ans[0].shape for x in ans[1:]])
		ans2=self.func(*ans)
		assert ans2.shape==tuple([len(x) for x in self.names]+[n])
		ans2[np.isnan(ans).any(axis=0)]=np.nan
		return ans2

class statf_singlestat(stat_base):
	def __init__(self,func_stat,stat,pts,**ka):
		"""Show constant value for stat
		func_stat:	Function to combine different points to one for the stat.
		stat:	Stats to be combined to single point
		pts:	Points to combine for the stat
		"""
		self.stat=stat
		super().__init__(**ka)
		self.val=func_stat(stat.compute(pts))
		if self.val.ndim==1:
			self.val=self.val.reshape(-1,1)
		assert self.val.shape==tuple([len(x) for x in self.names]+[len(pts)])
	def default_names(self):
		"""Use genes of original stat by default
		"""
		return self.stat.names
	def default_label(self):
		return 'Single value of '+self.stat.label
	def compute(self,pts):
		"""Computes stat
		pts:	Point list instance to compute stat
		Return:
		stat as np.array(shape=[len(x) for x in self.names]+[len(pts)])
		"""
		import numpy as np
		return np.repeat(self.val,len(pts),axis=-1)

def statf_initial(stat,pts,**ka):
	return statf_singlestat(lambda x,y:x,lambda x,y:y,stat,pts[[0]] if isinstance(pts,dictys.traj.point) else [pts[0]],**ka)

def statf_mean(stat,pts,**ka):
	"""Use mean value for stat
	Return:
	stat
	"""
	return statf_singlestat(lambda x,y:((x*y).sum(axis=-1)/y.sum(axis=-1)),stat,pts,**ka)
def statf_diff(stat,stat_base,label=None,**ka):
	"""Use value difference for stat
	Return:
	stat
	"""
	s1=stat-stat_base
	if label is None:
		s1.label='Difference in '+stat.label
	else:
		s1.label=label
	return s1

class statf_smooth(stat_base):
	def __init__(self,stat,traj,smoothen_func,**ka):
		"""Base class for statistics obtained from smoothening of their values at nodes
		stat:	Stat to smoothen
		traj:	Trajectory instance of lwang.pipeline.dynet2.traj.trajectory
		smoothen_func:	[name,args,keyword args] as in lwang.pipeline.dynet2.traj.trajc.smoothened
		namet:	Names of genes to compute stat for. Default is obtained from default_namet function.
		"""
		import numpy as np
		assert len(smoothen_func)==3
		self.traj=traj
		self.stat=stat
		super().__init__(**ka)
		self.func_smooth=self.traj.smoothened(stat.compute(np.arange(traj.nn)),smoothen_func[0],*smoothen_func[1],**smoothen_func[2])
	def default_names(self):
		"""Use genes of original stat by default
		"""
		return self.stat.names
	def default_label(self):
		return self.stat.label
	def compute(self,pts):
		"""Use designated smoothening function to compute stat values at points from nodes.
		pts:	Point list instance of lwang.pipeline.dynet2.traj.pointc
		Return:
		Stat value for each gene as np.array(shape=[len(x) for x in self.names]+[len(pts)])
		"""
		if not isinstance(pts,dictys.traj.point):
			raise TypeError('Smooth function only available at points not states.')
		ans=self.func_smooth(pts)
		assert ans.shape==tuple([len(x) for x in self.names]+[len(pts)])
		return ans

class statf_binarize(stat_base):
	def __init__(self,stat,statmask,*a,posrate=0.01,signed=True,**ka):
		"""Compute all edge weights of whole network with binarization if specified.
		varname:	Variable name used for network
		posrate:	Proportion of top edges when converting to binary network. Set to None to retain continuous network.
		"""
		self.posrate=posrate
		self.signed=signed
		self.stat=stat
		self.mask=statmask
		super().__init__(*a,**ka)
	def default_names(self):
		"""Use genes of original stat by default
		"""
		return self.stat.names
	def default_label(self):
		return self.stat.label
	def compute(self,pts):
		"""Computes network strength at each node 
		Return:
		Networks as np.array(shape=(regulators,targets,n_node))
		"""
		#Load individual networks for each node
		import numpy as np
		ans=self.stat.compute(pts)
		mask=self.mask.compute(pts)
		if self.signed:
			ans=np.abs(ans)
		cut=(self.posrate*mask.sum(axis=0).sum(axis=0)).astype(int)
		cut=[np.partition(x.ravel(),-y)[-y] for x,y in zip(ans.transpose(2,0,1),cut)]
		ans=ans>=cut
		assert ans.shape==tuple([len(y) for y in self.names]+[len(pts)])
		return ans

class stat_smooth_old(stat_base):
	def __init__(self,d,traj,smoothen_func,raw_states=True,**ka):
		"""Base class for statistics obtained from smoothening of their values at nodes
		d:		network object
		traj:	Trajectory instance of lwang.pipeline.dynet2.traj.trajc
		label:	Label of stat that is shown as coordindate label.
		smoothen_func:	[name,args,keyword args] as in lwang.pipeline.dynet2.traj.trajc.smoothened
		err_states:	Whether to throw error when smoothening for state computations.
		namet:	Names of genes to compute stat for. Default is obtained from default_namet function.
		"""
		assert len(smoothen_func)==3
		self.d=d
		self.traj=traj
		super().__init__(**ka)
		self.err_states=err_states
		self.func_smooth=self.traj.smoothened(self.compute_allstates(),smoothen_func[0],*smoothen_func[1],**smoothen_func[2])
	def compute(self,pts):
		"""Use designated smoothening function to compute stat values at points from nodes.
		pts:	Point list instance of lwang.pipeline.dynet2.traj.pointc
		Return:
		Stat value for each gene as np.array(shape=[len(x) for x in self.names]+[len(pts)])
		"""
		if not isinstance(pts,dictys.traj.point):
			raise TypeError('Values at states should not be smoothened.')
		ans=self.func_smooth(pts)
		assert ans.shape==tuple([len(x) for x in self.names]+[len(pts)])
		return ans

class stat_pseudotime(stat_base):
	def __init__(self,d,traj,pts,*a,**ka):
		"""Pseudotime statistic
		d:		Dataset object
		traj:	Trajectory instance of lwang.pipeline.dynet2.traj.trajc
		pts:	Actual point list instance to use for visualization
		"""
		self.d=d
		self.traj
		self.pts=pts
		super().__init__(*a,**ka)
	def default_names(self):
		"""Use all genes by default.
		"""
		return [self.d.nname]
	def default_label(self):
		return 'Pseudo-time'
	def compute(self,pts):
		"""Computes pseudotime at each point. All genes have the same value.
		pts:	Point list instance to compute pseudotime
		Return:
		pseudotime at each point as np.array(shape=[len(x) for x in self.names]+[len(pts)])
		"""
		import numpy as np
		if not isinstance(pts,dictys.traj.point):
			pts=self.traj.topoint()[pts]
		ans=(pts-self.pts[[0]]).T
		ans=np.repeat(ans,len(self.names[0]),axis=0)
		assert ans.shape==tuple([len(x) for x in self.names]+[len(pts)])
		return ans

class stat_lcpm(stat_base):
	def __init__(self,d,*a,cut=0.01,const=1,**ka):
		"""LogCPM stat. Specifically: log2 (CPM+const)
		d:		Dataset object
		const:	Constant to add to CPM before log.
		cut:	CPM above cut will be shown
		"""
		self.d=d
		self.cut=cut
		self.const=const
		super().__init__(*a,**ka)
	def default_names(self):
		"""Use all genes by default.
		"""
		return [self.d.nname]
	def default_label(self):
		return 'Log2 CPM'
	def compute(self,pts):
		"""Computes logCPM at each node
		Return:
		log2(CPM+const) as np.array(shape=(len(namet),n_node))
		"""
		import numpy as np
		if isinstance(pts,dictys.traj.point):
			raise ValueError('stat_lcpm should not be computed at any point. Use existing states or wrap with stat_smooth instead.')

		t1=set(self.d.nname)
		t1=np.nonzero([x not in t1 for x in self.names[0]])[0]
		if len(t1)>0:
			raise ValueError('Genes not found: {}'.format(','.join(self.d.nname[t1])))
		#Gene IDs for each TF
		namet=self.d.nname
		tdict=dict(zip(namet,range(len(namet))))
		t1=np.nonzero([x not in tdict for x in self.names[0]])[0]
		if len(t1)>0:
			raise ValueError('Genes not found: {}'.format(','.join(namet[t1])))
		#Compute logCPM
		if 'cpm' not in self.d.prop['ns']:
			raise ValueError('CPM results not found. Please recompute.')
		cpm=self.d.prop['ns']['cpm'][[tdict[x] for x in namet]][:,pts]
		ans=np.log2(cpm+self.const)
		ans[cpm<=self.cut]=np.nan
		return ans

class stat_net(stat_base):
	def __init__(self,d,*a,varname='w',**ka):
		"""Compute all edge weights of whole network with binarization if specified.
		varname:	Variable name used for network
		posrate:	Proportion of top edges when converting to binary network. Set to None to retain continuous network.
		"""
		self.d=d
		self.varname=varname
		super().__init__(*a,**ka)
	def default_label(self):
		return 'Network edge strength'
	def default_names(self):
		"""Use all genes of the specified role by default.
		"""
		# dynet=self.d.eprops['w']
		# assert dynet.ndim==3
		namet=[self.d.nname[x] for x in self.d.nids]
		return namet
	def compute(self,pts):
		"""Computes network strength at each node 
		Return:
		Networks as np.array(shape=(regulators,targets,n_node))
		"""
		import numpy as np
		if isinstance(pts,dictys.traj.point):
			raise ValueError('stat_net should not be computed at any point. Use existing states or wrap with stat_smooth instead.')
		#Load individual networks for each node
		dynet=np.take(self.d.prop['es'][self.varname],pts,axis=-1)
		#Gene IDs for each TF
		t1=[np.nonzero([x not in self.ndict[y] for x in self.names[y]])[0] for y in range(len(self.names))]
		for xi in range(len(self.names)):
			if len(t1[xi])>0:
				raise ValueError('Genes not found in given role: {}'.format(','.join(self.names[xi][t1[xi]])))
		dynet=dynet[[self.ndict[0][x] for x in self.names[0]]][:,[self.ndict[1][x] for x in self.names[1]]]
		return dynet

class statf_centrality_base(stat_base):
	def __init__(self,statnet,func,directed=False,roleaxis=0,label=None):
		self.stat=statnet
		self.func_centrality=func
		self.directed=directed
		self.roleaxis=roleaxis
		super().__init__(label=label)
	def default_names(self):
		return [self.stat.names[self.roleaxis]]
	def default_label(self):
		return 'Centrality'
	def compute(self,pts):
		"""Computes log number of neighbors from networks at each point
		pts:	Point list instance to compute stat
		Return:
		log2(number of targets+1) as np.array(shape=(len(namet),len(pts)))
		"""
		import numpy as np
		import networkx as nx
		dynet=self.stat.compute(pts)
		ans=[]
		for xi in range(len(pts)):
			g=list(zip(*np.nonzero(dynet[:,:,xi])))
			g=(nx.DiGraph if self.directed else nx.Graph)(g)
			t1=self.func_centrality(g)
			t2=np.array(list(t1.keys()))
			t1=np.array([t1[x] for x in t2])
			t3=np.ones([len(self.names[0])],dtype=float)*np.nan
			t4=t2<len(self.names[0])
			t3[t2[t4]]=t1[t4]
			assert t3.shape==(len(self.names[0]),)
			ans.append(t3)
		ans=np.array(ans).T
		assert ans.shape==(len(self.names[0]),len(pts))
		return ans

class statf_centrality_degree(stat_base):
	"""Use roleaxis=0 for outdegree and =1 for indegree."""
	def __init__(self,statnet,statmask=None,roleaxis=0):
		"""
		Compute degree centrality for network.
		statnet:	stat for network
		statmask:	stat for mask. When specified, computes the degree centrality rate, i.e. degree/node_count, instead of degree
		roleaxis:	Axis to compute degree for. 0 for outdegree and 1 for in degree.
		"""
		self.stat=statnet
		self.mask=statmask
		self.roleaxis=roleaxis
		super().__init__()
	def default_names(self):
		return [self.stat.names[self.roleaxis]]
	def default_label(self):
		return 'Outdegree centrality' if self.roleaxis==0 else 'Indegree centrality'
	def compute(self,pts):
		"""Computes log number of neighbors from networks at each point
		pts:	Point list instance to compute stat
		Return:
		log2(number of targets+1) as np.array(shape=(len(namet),len(pts)))
		"""
		import numpy as np
		dynet=self.stat.compute(pts)
		if self.mask is None:
			dynet=dynet.sum(axis=1-self.roleaxis)
		else:
			mask=self.mask.compute(pts)
			dynet=(dynet*mask).sum(axis=1-self.roleaxis)/(mask.sum(axis=1-self.roleaxis)+1E-300)
		assert dynet.shape==(len(self.names[0]),len(pts))
		return dynet

def statf_centrality_eigenvector(statnet,label='Eigenvalue centrality',**ka):
	import networkx as nx
	return statf_centrality_base(statnet,nx.eigenvector_centrality,label=label,**ka)
def statf_centrality_betweenness(statnet,label='Betweenness centrality',**ka):
	import networkx as nx
	print('Betweenness centrality: very slow!')
	return statf_centrality_base(statnet,nx.betweenness_centrality,label=label,**ka)
def statf_centrality_closeness(statnet,label='Closeness centrality',**ka):
	import networkx as nx
	print('Closeness centrality: very slow!')
	return statf_centrality_base(statnet,nx.closeness_centrality,label=label,**ka)
def statf_lnneighbor(statnet,label='Log2 (Outdegree + 1)',const=1,**ka):
	import numpy as np
	return statf_function(lambda x:np.log2(x+const),[statf_centrality_degree(statnet,**ka)],label=label)












































#
