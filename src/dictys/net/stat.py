#!/usr/bin/python3
# Lingfei Wang, 2022. All rights reserved.

"""
Statistics of networks for data visualization and export.
"""

from __future__ import annotations
import abc
from typing import Union,Callable,Tuple,Optional
import numpy as np
import dictys.traj
import dictys.net
from dictys.utils.numpy import NDArray,ArrayLike

def _getitem(key,v:NDArray)->NDArray:
	"""
	Get items from numpy.array
	key:	iterable of keys. Each key is a iterable or individual value.
	v:		numpy.array to get items from
	Return:	numpy.array
	"""
	sid=0
	for xi in key:
		if hasattr(xi,'__len__'):
			v=v.swapaxes(sid,0)[xi].swapaxes(sid,0)
			sid+=1
		else:
			v=np.take(v,xi,axis=sid)
	return v

class base(metaclass=abc.ABCMeta):
	"""
	Abstract base class for stat of network.
	"""
	def __init__(self,names:Optional[list[ArrayLike[str]]]=None,label:Optional[str]=None):
		"""
		Base class for statistics for each gene
		names:	List of names of each axis of output stat, except last axis which is always time. Default is obtained from default_names function.
		label:	Label of stat that is shown as coordindate label unless overidden.
		"""
		if label is None:
			label=self.default_label()
		self.label=label
		if names is None:
			names=self.default_names()
		assert len(names)>0
		self.names=[np.array(x) for x in names]
		self.ndict=[dict(zip(x,range(len(x)))) for x in names]
	@abc.abstractmethod
	def compute(self,pts:Union[dictys.traj.point,NDArray[np.int_]])->NDArray:
		"""
		Use this function to compute stat values at each state or point
		pts:	Point list instance of dictys.traj.point, or state list as list of int
		Return:
		Stat values as numpy.array(shape=(...,len(pts))) . Use nan to hide value or set as invalid.
		"""
	@abc.abstractmethod
	def default_names(self)->list[NDArray[str]]:
		"""
		Use this function to determine the default names for each axis.
		Note that only names shared with other stats will be visualized.
		Return:
		List of list of names for each axis.
		"""
	def default_lims(self,pts:Union[dictys.traj.point,NDArray[np.int_],None]=None,names:Optional[list[ArrayLike[str]]]=None,expansion:float=0.02)->ArrayLike:
		"""
		Use this function to determine the default limits of the stat.
		This implementation uses min/max of stat values.
		pts:	Point list instance to compute min/max
		namet:	Names used to compute limits. Defaults to all.
		expansion:	Expand limits by this relative amount on each side.
		Return:
		Limits as [min,max]
		"""
		if pts is None:
			raise NotImplementedError
		ans=self.compute(pts)
		if names is not None:
			assert len(names)==len(self.names)
			assert all(all(y in self.ndict[x] for y in names[x]) for x in range(len(self.names)))
			for xi in range(len(names)):
				ans=ans.swapaxes(0,xi)[[self.ndict[xi][x] for x in names[xi]]].swapaxes(0,xi)
		t1=np.isfinite(ans)
		#Ignore NAN except all are NAN
		if not t1.any():
			return [np.nan,np.nan]
		ans=[ans[t1].min(),ans[t1].max()]
		t1=(ans[1]-ans[0])*expansion
		return [ans[0]-t1,ans[1]+t1]
	@abc.abstractmethod
	def default_label(self)->str:
		"""
		Use this function to determine the label of this stat
		Return:
		Label as str
		"""
	#Arithmetic operations between stats
	def _operator_(self,op:Callable,other,symbol:str)->base:
		if isinstance(other,base):
			return function(op,[self,other],label=f'({self.label}){symbol}({other.label})')
		return function(op,[self],other,label='({}){}{}'.format(self.label,symbol,str(base)))
	def __add__(self,other:base)->base:
		from operator import add as op
		return self._operator_(op,other,'+')
	def __sub__(self,other:base)->base:
		from operator import sub as op
		return self._operator_(op,other,'-')
	def __mul__(self,other:base)->base:
		from operator import mul as op
		return self._operator_(op,other,'*')
	def __truediv__(self,other:base)->base:
		from operator import truediv as op
		return self._operator_(op,other,'/')
	def __lt__(self,other:base)->base:
		from operator import lt as op
		return self._operator_(op,other,'<')
	def __le__(self,other:base)->base:
		from operator import le as op
		return self._operator_(op,other,'<=')
	def __eq__(self,other:base)->base:
		from operator import eq as op
		return self._operator_(op,other,'==')
	def __ne__(self,other:base)->base:
		from operator import ne as op
		return self._operator_(op,other,'!=')
	def __gt__(self,other:base)->base:
		from operator import gt as op
		return self._operator_(op,other,'>')
	def __ge__(self,other:base)->base:
		from operator import ge as op
		return self._operator_(op,other,'>=')
	def __and__(self,other:base)->base:
		from operator import and_ as op
		return self._operator_(op,other,'&')
	def __or__(self,other:base)->base:
		from operator import or_ as op
		return self._operator_(op,other,'|')
	def __getitem__(self,key)->base:
		"""Subset stat as a substat"""
		from functools import partial
		if not isinstance(key, tuple):
			key=(key,)
		if len(key)>len(self.names):
			raise KeyError('Takes at most {} dimensions'.format(len(self.names)))
		names=[]
		keys=[]
		for xi in range(len(key)):
			key1=key[xi]
			if isinstance(key1,slice):
				key1=list(range(*key1.indices(len(self.names[xi]))))
			if hasattr(key1,'__len__') and not isinstance(key1,str):
				if isinstance(key1[0],str):
					key1=[self.ndict[xi][x] for x in key1]
				elif not isinstance(key1[0],int):
					raise TypeError('Key must be type int or str.')
				names.append(self.names[xi][key1])
			keys.append(key1)
		names+=self.names[len(key):]
		return function(partial(_getitem,keys),[self],names=names,label=self.label)
	# def match_names(self,*others):
	# 	raise NotImplementedError

class const(base):
	"""
	Show constant(/state-invariant) value for stat.
	"""
	isconst=True
	def __init__(self,val:NDArray,names:list[ArrayLike[str]],label:str='Constant',**ka):
		"""
		Show constant(/state-invariant) value for stat
		val:	Value to show
		names:	Names of each value
		label:	Label of stat
		ka:		Keyword arguments passed to parent class
		"""
		assert val.shape==tuple(len(x) for x in names)
		self.val=val
		self.default_names_=names
		self.default_label_=label
		super().__init__(**ka)
	def default_names(self)->list[ArrayLike[str]]:
		return self.default_names_
	def default_label(self)->str:
		return self.default_label_
	def compute(self,pts:Union[dictys.traj.point,NDArray[np.int_]])->NDArray:
		return np.repeat(self.val.reshape(*self.val.shape,1),len(pts),axis=-1)

class function(base):
	"""
	Stat that is a function of other stat(s)
	"""
	def __init__(self,func:Callable[Tuple[NDArray,...],NDArray],stats:list[base],*a,isconst:Optional[bool]=None,**ka):
		"""
		Stat that is a function of other stat(s)
		func:	Function to combine other stats. Should have self.compute=func(*[x.compute(...) for x in stats],*a).
		stats:	List of stats whose final outputs (compute) will be operated on by func.
		a:		Other arguments passed to func
		isconst:Overide whether the function stat is a constant. By default, it is constant only if all dependent stats are constant.
		label:	Label of stat that is shown as coordindate label.
		"""
		self.n=len(stats)
		assert self.n>0
		self.func=func
		self.a=a
		self.stats=stats
		if isconst is None:
			isconst=all(hasattr(x,'isconst') and x.isconst for x in stats)
		self.isconst=isconst
		super().__init__(**ka)
	def default_names(self):
		"""
		Use intersection names of all stats by default. Only suitable when all stats have the same number of dimensions
		"""
		from functools import reduce
		from operator import and_
		if any(len(x.names)!=len(self.stats[0].names) for x in self.stats[1:]):
			raise NotImplementedError
		names=[[set(z) for z in y] for y in zip(*[x.names for x in self.stats])]
		names=[sorted(list(reduce(and_,x))) for x in names]
		return names
	def default_label(self):
		return 'Function'
	def compute(self,pts):
		"""
		Computes stat for states or points
		pts:	Point list instance of dictys.traj.point, or state list as list of int
		Return:
		Stat values as numpy.array(shape=(...,len(pts))) . Use nan to hide value or set as invalid.
		"""
		n=len(pts)
		ans=[x.compute(pts) for x in self.stats]
		assert all(ans[x].ndim==len(self.stats[x].names)+1 for x in range(self.n))
		assert all(ans[x].shape==tuple([len(y) for y in self.stats[x].names]+[n]) for x in range(self.n))
		ans2=self.func(*ans,*self.a)
		assert ans2.shape==tuple([len(x) for x in self.names]+[n])
		return ans2

class fsinglestat(const):
	"""
	Show constant value for stat by combining different points.
	"""	
	def __init__(self,func_stat:Callable[Tuple[NDArray,...],NDArray],stat:list[base],pts:Union[dictys.traj.point,NDArray[np.int_]],**ka):
		"""
		Show constant value for stat by combining different points.
		func_stat:	Function to combine different points to one for the stat.
		stat:	Stats to be combined to single point
		pts:	Points to combine for the stat
		"""
		val=func_stat(stat.compute(pts))
		ka1={'label':'Single value of '+stat.label}
		ka1.update(ka)
		super().__init__(val,stat.names,**ka1)

def finitial(stat:base,pts:Union[dictys.traj.point,NDArray[np.int_]],**ka)->base:
	"""
	Using initial value as a constant stat.
	"""
	return fsinglestat(lambda x:x.ravel(),stat,pts[[0]],**ka)

def fmean(stat:base,pts:Union[dictys.traj.point,NDArray[np.int_]],**ka)->base:
	"""
	Use mean value for stat.
	"""
	return fsinglestat(lambda x,y:((x*y).sum(axis=-1)/y.sum(axis=-1)),stat,pts,**ka)

def fdiff(stat:base,stat_base:base,label:str=None)->base:
	"""
	Use value difference for stat: stat-stat_base
	"""
	s1=stat-stat_base
	if label is None:
		s1.label='Difference in '+stat.label
	else:
		s1.label=label
	return s1

def fmasked(stat:base,stat_mask:base,**ka)->base:
	"""
	Masking entries to NAN with a boolean mask matrix. False items will be masked.
	"""
	def masking(s,m):
		assert s.shape==m.shape
		s=s.copy()
		s[~m]=np.nan
		return s
	ans=function(masking,[stat,stat_mask],**ka)
	return ans

class fsmooth(base):
	"""
	Base class for statistics obtained from smoothing of their values at points/states
	"""
	def __init__(self,stat:base,pts:Union[dictys.traj.trajectory,dictys.traj.point],smoothen_func:Tuple[str,Tuple,dict],**ka):
		"""
		Base class for statistics obtained from smoothing of their values at points/states
		stat:	Stat to smoothen
		pts:	Trajectory or point instance of dictys.traj.trajectory or dictys.traj.point
		smoothen_func:	[name,args,keyword args] as in dictys.point.smoothened
		"""
		assert len(smoothen_func)==3
		if isinstance(pts,dictys.traj.point):
			n=len(pts)
			pts2=pts
		elif isinstance(pts,dictys.traj.trajectory):
			n=pts.nn
			pts2=np.arange(n)
		else:
			raise TypeError('pts must be dictys.traj.trajectory or distys.traj.point')
		self.stat=stat
		super().__init__(**ka)
		self.func_smooth=pts.smoothened(stat.compute(pts2),*smoothen_func[1],func_name=smoothen_func[0],**smoothen_func[2])
	def default_names(self):
		return self.stat.names
	def default_label(self):
		return self.stat.label
	def compute(self,pts):
		"""
		Use designated smoothening function to compute stat values at points from nodes.
		pts:	Point list instance of dictys.traj.point, or state list as list of int
		Return:
		Stat values as numpy.array(shape=(...,len(pts))) . Use nan to hide value or set as invalid.
		"""
		if not isinstance(pts,dictys.traj.point):
			raise TypeError('Smooth function only available at points not states.')
		ans=self.func_smooth(points=pts)
		assert ans.shape==tuple([len(x) for x in self.names]+[len(pts)])
		return ans

class fbinarize(base):
	"""
	Convert continuous network stat to binary network stat.
	"""
	def __init__(self,stat:base,*a,statmask:Optional[base]=None,sparsity:Optional[float]=0.01,signed:bool=True,**ka):
		"""
		Convert continuous network stat to binary network stat.
		stat:		Stat for continuous network as numpy.ndarray(shape=(n_reg,n_target),dtype=float)
		statmask:	Stat for network mask indicating which edges are tested in the continuous network,
					as numpy.ndarray(shape=(n_reg,n_target),dtype=bool)
		sparsity:	Proportion of significant edges when converting to binary network.
					If None, sets all non-zero (if signed) or positive (if not signed) edges as significant.
		signed:		Whether continuous network is signed.
					If so, larger absolute values indicate stronger edge.
					If not, larger values indicate stronger edge.
		a,
		ka:			Arguments and keyword arguments passed to parent class
		"""
		self.sparsity=sparsity
		self.signed=signed
		self.stat=stat
		self.mask=statmask
		super().__init__(*a,**ka)
	def default_names(self):
		return self.stat.names
	def default_label(self):
		return self.stat.label
	def compute(self,pts):
		"""
		Computes binary network strength at each point
		pts:	Point list instance of dictys.traj.point, or state list as list of int
		Return:
		Networks as numpy.array(shape=(n_reg,n_target,n_pts))
		"""
		#Load individual networks for each node
		ans=self.stat.compute(pts)
		if self.signed:
			ans=np.abs(ans)
		#Determine cutoff for each point
		if self.sparsity is None:
			cut=np.array([ans[...,x][ans[...,x]>0].min() for x in range(ans.shape[-1])])
		else:
			if self.mask is None:
				cut=np.repeat(int(self.sparsity*np.prod(ans.shape[:-1])),ans.shape[-1])
			else:
				mask=self.mask.compute(pts)
				cut=(self.sparsity*mask.sum(axis=0).sum(axis=0)).astype(int)
			assert cut.shape==(ans.shape[-1],)
			cut=np.array([np.partition(x.ravel(),-y)[-y] for x,y in zip(ans.transpose(2,0,1),cut)])
		assert cut.shape==(ans.shape[-1],)
		ans=ans>=cut
		assert ans.shape==tuple([len(y) for y in self.names]+[len(pts)])
		return ans

class pseudotime(base):
	"""
	Statistic to output pseudotime
	"""
	def __init__(self,d:dictys.net.network,pts:Union[dictys.traj.point,NDArray[np.int_]],*a,traj:Optional[dictys.traj.trajectory]=None,**ka):
		"""
		Statistic to output pseudotime
		d:		Dataset object
		pts:	Actual point list instance to use for visualization
		traj:	Trajectory instance of dictys.traj.trajectory. Defaults to pts.p.
		"""
		self.d=d
		if traj is None:
			traj=pts.p
		self.traj=traj
		self.pts=pts
		super().__init__(*a,**ka)
	def default_names(self):
		"""
		Use all genes by default.
		"""
		return [self.d.nname]
	def default_label(self):
		return 'Pseudo-time'
	def compute(self,pts):
		"""
		Computes pseudotime at each point. All genes have the same value.
		pts:	Point list instance to compute pseudotime
		Return:
		pseudotime at each point as np.array(shape=[len(x) for x in self.names]+[len(pts)])
		"""
		if not isinstance(pts,dictys.traj.point):
			pts=self.traj.topoint()[pts]
		ans=(pts-self.pts[[0]]).T
		ans=np.repeat(ans,len(self.names[0]),axis=0)
		assert ans.shape==tuple([len(x) for x in self.names]+[len(pts)])
		return ans

class lcpm(base):
	"""
	LogCPM stat. Specifically: log2 (CPM+const)
	"""
	def __init__(self,d:dictys.net.network,*a,cut:float=0.01,constant:float=1,**ka):
		"""
		LogCPM stat. Specifically: log2 (CPM+const)
		d:			Dataset object
		constant:	Constant to add to CPM before log.
		cut:		CPM below cut will be hidden
		"""
		self.d=d
		self.cut=cut
		self.const=constant
		super().__init__(*a,**ka)
	def default_names(self):
		return [self.d.nname]
	def default_label(self):
		return 'Log2 CPM'
	def compute(self,pts):
		"""Computes logCPM at each node
		Return:
		log2(CPM+const) as np.array(shape=(n_gene,len(pts)))
		"""
		if isinstance(pts,dictys.traj.point):
			raise TypeError('lcpm should not be computed at any point. Use existing states or wrap with smooth instead.')

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
		ans[cpm<self.cut]=np.nan
		return ans

class sprop(base):
	"""
	Base class of state dependent properties directly read from dataset object
	"""
	def __init__(self,d:dictys.net.network,ptype:str,pname:str,*a,names_pref:list[str]=[],**ka):
		"""
		Base class of state dependent properties directly read from dataset object
		d:		Dataset object
		ptype:	Property type (key in d.prop)
		pname:	Property name (key in d.prop[ptype])
		names_pref:	Prefixes added to the names of each dimension
		"""
		import itertools
		assert ptype in d.prop and pname in d.prop[ptype]

		#Name collection
		sid=np.array(list(itertools.chain.from_iterable([[False,False] if x=='e' else [x=='s'] for x in ptype])))
		names=list(itertools.chain.from_iterable([[getattr(d,x+'name')] if x!='e' else [d.nname[y] for y in d.nids] for x in ptype]))
		if not sid.any():
			raise ValueError('Stat sprop is only for state-dependent properties.')
		names=[names[x] for x in np.nonzero(~sid)[0]]
		t1=[d.prop[ptype][pname].ndim,len(names_pref)+len(names)+sid.sum()]
		if t1[0]!=t1[1]:
			raise ValueError("Incorrect prefix dimension for property {}. Final dimension: {}. Correct dimension: {}.".format(pname,t1[1],t1[0]))
		t1=[d.prop[ptype][pname].shape[:len(names_pref)],tuple(len(x) for x in names_pref)]
		if t1[0]!=t1[1]:
			raise ValueError("Incorrect prefix shape for property {}. Final shape: {}. Correct shape: {}.".format(pname,t1[1],t1[0]))
		names=names_pref+names
		self.d=d
		self.ptype=ptype
		self.pname=pname
		self.sid=np.nonzero(sid)[0]
		self.default_names_=names
		super().__init__(*a,**ka)
	def default_label(self):
		return self.pname
	def default_names(self):
		return self.default_names_
	def compute(self,pts):
		"""
		Obtains property at each state.
		Return:
		Networks as np.array(shape=(...,len(pts)))
		"""
		if isinstance(pts,dictys.traj.point):
			raise ValueError('Property should not be computed at any point. Use existing states or wrap with smooth instead.')
		#Load individual networks for each node
		ans=[]
		for xi in pts:
			ans1=self.d.prop[self.ptype][self.pname]
			for xj in self.sid[::-1]:
				ans1=np.take(ans1,xi,axis=xj)
			ans.append(ans1)
		assert all(x.shape==ans[0].shape for x in ans[1:])
		ans=np.array(ans)
		ans=ans.transpose(*list(range(1,ans.ndim)),0)
		return ans

class net(sprop):
	"""
	Compute all edge weights of whole network.
	"""
	def __init__(self,d:dictys.net.network,*a,varname:str='w',**ka):
		super().__init__(d,'es',varname,*a,label='Network edge strength',**ka)

class netmask(sprop):
	"""
	Compute all edge masks of whole network.
	"""
	def __init__(self,d:dictys.net.network,*a,varname:str='mask',**ka):
		super().__init__(d,'es',varname,*a,label='Network edge mask',**ka)

class fcentrality_base(base):
	"""
	Base class for centrality measure from network stat
	"""
	def __init__(self,statnet:base,func:Callable,directed:bool=False,roleaxis:int=0,label:Optional[str]=None):
		"""
		Base class for centrality measure from network stat
		statnet:	Stat for network
		func:		networkx function to compute the centrality from a networkx.Graph or networkx.DiGraph
		directed:	Whether to created a networkx.DiGraph instead of networkx.Graph
		roleaxis:	Axis to compute centrality for. Use 0 for regulators' centrality and 1 for targets' centrality.
		label:		Custom label for the statistic
		"""
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
		"""
		Computes centrality for states or points
		pts:	Point list instance of dictys.traj.point, or state list as list of int
		Return:
		Centrality as numpy.array(shape=(n,len(pts))) . Use nan to hide value or set as invalid.
		"""
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

class fcentrality_degree(base):
	"""
	Degree centrality stat for network. 
	"""
	def __init__(self,statnet:base,roleaxis:int=0,statmask:Optional[base]=None,weighted_sparsity:Optional[float]=None,constant:float=0):
		"""
		Degree centrality stat for network.
		statnet:	stat for network. Both binarized or continuous networks are accepted.
		roleaxis:	Axis to compute degree for. 0 for outdegree and 1 for indegree.
		statmask:	stat for mask. When specified, computes the degree centrality rate, i.e. (degree+constant)/(node_count+constant), instead of degree
		weighted_sparsity:	If set, treat statnet as a continuous network and compute weighted degree centrality.
							Value is used to normalize centrality so the sum equals to that of unweighted degree centrality of the same network binarized under the specified sparsity.
							Ignored if statmask is set.
		constant:	Constant to add to degree centrality. Particularly useful when statmask is set. Always used.
		"""
		self.stat=statnet
		self.mask=statmask
		self.roleaxis=roleaxis
		self.constant=constant
		self.weighted_sparsity=weighted_sparsity
		super().__init__()
	def default_names(self):
		return [self.stat.names[self.roleaxis]]
	def default_label(self):
		return 'Outdegree centrality' if self.roleaxis==0 else 'Indegree centrality'
	def compute(self,pts):
		"""
		Computes degree centrality for states or points
		pts:	Point list instance of dictys.traj.point, or state list as list of int
		Return:
		Centrality as numpy.array(shape=(n,len(pts))). Use nan to hide value or set as invalid.
		"""
		dynet=np.abs(self.stat.compute(pts))
		if self.mask is None:
			if self.weighted_sparsity is not None:
				#Current sparsity
				t1=dynet.sum(axis=0).sum(axis=0)/(dynet.shape[0]*dynet.shape[1])
				#Normalized to specified sparsity
				dynet=dynet*(self.weighted_sparsity/(t1+1E-300))
			dynet=dynet.sum(axis=1-self.roleaxis)+self.constant
		else:
			mask=self.mask.compute(pts)
			dynet=((dynet*mask).sum(axis=1-self.roleaxis)+self.constant)/((mask!=0).any(axis=self.roleaxis).sum(axis=0)+self.constant)
		assert dynet.shape==(len(self.names[0]),len(pts))
		return dynet

def fcentrality_eigenvector(statnet:base,label:str='Eigenvalue centrality',**ka)->base:
	"""
	Eigenvalue centrality stat
	statnet:	stat for network
	"""
	import networkx as nx
	return fcentrality_base(statnet,nx.eigenvector_centrality,label=label,**ka)
def fcentrality_betweenness(statnet:base,label:str='Betweenness centrality',**ka)->base:
	"""
	Betweenness centrality stat
	statnet:	stat for network
	"""
	import networkx as nx
	print('Betweenness centrality: very slow!')
	return fcentrality_base(statnet,nx.betweenness_centrality,label=label,**ka)
def fcentrality_closeness(statnet:base,label:str='Closeness centrality',**ka)->base:
	"""
	Closeness centrality stat
	statnet:	stat for network
	"""
	import networkx as nx
	print('Closeness centrality: very slow!')
	return fcentrality_base(statnet,nx.closeness_centrality,label=label,**ka)

def flnneighbor(statnet:base,label:str='Log2 (Outdegree + 1)',constant:float=1,statmask:Optional[base]=None,**ka)->base:
	"""
	Log2 outdegree centrality stat. Specifically, log2 (Outdegree + const)
	statnet:	stat for network
	constant:	Constant to add before log2.
	statmask:	If set, computes relative log2 outdegree: log2 (outdegree + const) - log2 (max possible outdegree allowd by statmask + const)
	ka:			Keyword arguments passed to fcentrality_degree
	"""
	return function(np.log2,[fcentrality_degree(statnet,statmask=statmask,constant=constant,**ka)],label=label)

class flayout_base(base):
	"""
	Compute layout coordinates of nodes from network.
	"""
	def __init__(self,statnet:base,layout_func:Callable[NDArray,NDArray],ndim:int=2,pts:Union[dictys.traj.point,NDArray,None]=None,netscale:float=1,**ka):
		"""
		Compute layout coordinates of nodes from network.
		statnet:	stat of network edge weight
		layout_func:Function to compute layout of network in networkx format, e.g. dictys.net.layout._fruchterman_reingold.
		ndim:		Number of dimensions of coordinates
		pts:		Points on the trajectory to compute layout for
		netscale:	Normalization scale for network edge strength.
		ka:			Keyword arguments passed to flayout_base.compute_all, including
			nodrop_reg:	Whether to keep regulators even if they do not have any edge
			abs:		Whether to use the absolute value of edge strength. If not, uses raw value where negative values are weaker than positive.
			scale:		How to scale network layout coordindates. Accepts:
				none:	No scaling
				size:	Use fixed average distance from figure center
				length:	Use fixed average length of each edge
			rand_expand:	Extra distance to put newly added nodes from the existing network
		"""
		assert len(statnet.names)==2
		names=sorted(list(set(statnet.names[0])|set(statnet.names[1])))
		if pts is None:
			pts=statnet.pts
		self.func=layout_func
		self.default_names_=[names,[f'Dim {x}' for x in range(1,ndim+1)]]
		self.stat=statnet
		self.pts=pts
		self.ndim=ndim
		self.netscale=netscale
		super().__init__()
		self.compute_all(**ka)
	def default_names(self):
		return self.default_names_
	def default_label(self):
		return 'Coordinates'
	def init_pos(self,m:NDArray,sep:float=1,always:list[int]=[])->NDArray:
		"""
		Initialize node positions for the initial network.
		m:	Directed network edge matrix as numpy.array(shape=(n_node,n_node))
		sep:	Initial separation between disconnected subnetworks
		always:	List of indices to always show (even without neighbors)
		Return:	numpy.ndarray(shape=(n_node,n_dim)) as initial node positions. NAN means hidden.
		"""
		import itertools
		import networkx as nx
		g=nx.Graph()
		if self.netscale is not None:
			t1=m[(~np.isnan(m))&(m!=0)]
			if not t1.any():
				return np.ones((m.shape[0],self.ndim),dtype=float)*np.nan
			m=m*(self.netscale/np.sqrt((t1**2).mean()))			

		t1=m!=0
		t1=np.nonzero(t1|t1.T)
		g.add_edges_from(np.array(t1).T)
		nids=[np.array(list(x)) for x in nx.connected_components(g) if len(x)>1]
		#Force include always shown nodes
		t1=set(list(itertools.chain.from_iterable(nids)))
		nids+=[np.array([x]) for x in always if x not in t1]
		ans=np.ones((m.shape[0],self.ndim),dtype=float)*np.nan
		xbase=0
		for xi in nids:
			#Initialize every subnetwork
			n=len(xi)
			if n==1:
				ans[xi]=[[xbase,0]]
				xbase+=sep
				continue
			m1=m[xi][:,xi]
			posinit=np.random.rand(n,self.ndim)*np.sqrt(n)*sep
			#Compute positions
			pos0=self.func(m1,posinit)
			pos0[:,0]-=pos0[:,0].min()
			pos0[:,1]-=pos0[:,1].mean()
			ans[xi,0]=pos0[:,0]+xbase
			ans[xi,1]=pos0[:,1]
			xbase+=pos0[:,0].max()+sep
		#Reinitialize full network twice
		nids=list(itertools.chain.from_iterable(nids))
		m1=m[nids][:,nids]
		for xi in range(2):
			ans[nids]=self.func(m1,ans[nids])
		return ans
	def compute_all(self,nodrop_reg:bool=False,sabs:bool=True,scale:str='size',rand_expand:float=0.1)->None:		# noqa: C901
		"""
		Compute node locations for all time points.

		nodrop_reg:	Whether to keep regulators even if they do not have any edge
		sabs:		Whether to use the absolute value of edge strength. If not, uses raw value where negative values are weaker than positive.
		scale:		How to scale network layout coordindates. Accepts:
			none:	No scaling
			size:	Use fixed average distance from figure center
			length:	Use fixed average length of each edge
		rand_expand:	Extra distance to put newly added nodes from the existing network
		"""
		n=len(self.names[0])
		ans=[]
		tmap0=None

		stat=self.stat.compute(self.pts)
		if sabs:
			stat=np.abs(stat)
		namemap=[[self.ndict[0][x] for x in self.stat.names[y]] for y in range(2)]
		m=np.zeros((n,n,len(self.pts)),dtype=float)
		for xi in range(len(self.stat.names[0])):
			m[namemap[0][xi],namemap[1]]=stat[xi]
		if self.netscale is not None:
			t1=m[(~np.isnan(m))&(m!=0)]
			if not t1.any():
				self.pos=np.ones((m.shape[0],self.ndim,len(self.pts)),dtype=float)*np.nan
				return
			m=m*(self.netscale/np.sqrt((t1**2).mean()))

		del stat
		#Prefilter nodes to show
		t1=(m!=0).any(axis=2)
		tmap=t1.any(axis=0)|t1.any(axis=1)
		if nodrop_reg:
			t1=np.zeros_like(tmap)
			t1[namemap[0]]=True
			tmap|=t1
		self.names[0]=self.names[0][tmap]
		self.ndict[0]=dict(zip(self.names[0],range(n)))
		namemap=[[self.ndict[0][x] for x in self.stat.names[y] if x in self.ndict[0]] for y in range(2)]
		m0=m[tmap][:,tmap]
		#Initialize positions
		pos0=self.init_pos(m0[:,:,0],always=(m0!=0).any(axis=1) if nodrop_reg else [])
		#Layout for each point
		for xi in range(len(self.pts)):
			m=m0[:,:,xi]
			#Prefilter nodes to show
			t1=(m!=0)
			tmap=t1.any(axis=0)|t1.any(axis=1)
			if nodrop_reg:
				t1=np.zeros_like(tmap)
				t1[namemap[0]]=True
				tmap|=t1
			elif not tmap.any():
				ans.append(np.ones((len(tmap),self.ndim),dtype=float)*np.nan)
				pos0=None
				tmap0=None
				continue
			m=m[tmap][:,tmap]
			#Initial positions for new nodes
			posinit=np.ones((len(tmap),self.ndim),dtype=float)*np.nan
			if tmap0 is not None:
				posinit[tmap0]=pos0
				vrange=[pos0.min(axis=0),pos0.max(axis=0)]
			else:
				vrange=[np.zeros(self.ndim),np.ones(self.ndim)]
			posinit=posinit[tmap]
			vrange=[(vrange[1]+vrange[0])/2,(vrange[1]-vrange[0]).max()/2+rand_expand]
			vrange=[vrange[0]-vrange[1],vrange[0]+vrange[1]]
			t1=np.isnan(posinit).any(axis=1)
			posinit[t1]=np.random.rand(t1.sum(),self.ndim)*(vrange[1]-vrange[0])+vrange[0]
			assert not np.isnan(posinit).any()
			assert not (posinit==0).all(axis=1).any()
			#Compute positions
			pos0=self.func(m,posinit)
			tmap0=tmap
			#Format & normalize output
			ans1=pos0
			ans1=ans1-ans1.mean(axis=0)
			if scale=='none':
				pass
			elif scale=='size':
				ans1=ans1/(np.sqrt((ans1**2).mean())+1E-300)
			elif scale=='length':
				t1=np.nonzero(m)
				ans1=(ans1.T/(np.sqrt(((ans1[t1[0]]-ans1[t1[1]])**2).sum(axis=1).mean())+1E-300)).T
			else:
				raise ValueError(f'Unknown scale {scale}.')
			ans2=np.ones((len(tmap),self.ndim),dtype=float)*np.nan
			ans2[tmap]=ans1
			ans.append(ans2)
		self.pos=np.array(ans).transpose(1,2,0)
	def ptsmap(self,pts:dictys.traj.point)->NDArray:
		"""
		Maps points to known points.
		"""
		if not isinstance(pts,dictys.traj.point):
			return None
		#Map given points
		t1=[np.nonzero(x)[0] for x in pts-self.pts==0]
		t1=np.array([x[0] for x in t1 if len(x)>0])
		return t1
	def compute(self,pts):
		"""
		Computes node coordinates from the provided layout for states or known points
		pts:	Point list instance of dictys.traj.point, or state list as list of int
		Return:
		Node coordindates as numpy.array(shape=(n,n_dim,len(pts))) . Use nan to hide value or set as invalid.
		"""
		if isinstance(pts,dictys.traj.point):
			#Map given points to self.pts
			pts2=self.ptsmap(pts)
			if pts2 is None or len(pts2)!=len(pts):
				raise ValueError('Property should not be computed at any point. Use only input points or wrap with smooth instead.')
			pts=pts2
		else:
			t1=dict(zip(self.pts,range(len(self.pts))))
			pts=np.array([t1[x] for x in pts])
		return self.pos[:,:,pts]









































#
