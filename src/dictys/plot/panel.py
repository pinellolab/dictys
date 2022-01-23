#!/usr/bin/python3
# Lingfei Wang, 2020-2022. All rights reserved.

"""
Visualization panels of networks
"""

class panel_base:
	def __init__(self,ax,d,pts):
		"""Base class to visualize single panel for dynamic network
		ax:	Axis
		d:	Dataset object
		pts:Points of path to visualize network"""
		self.ax=ax
		self.d=d
		self.pts=pts
	def init(self):
		"""Draws initial canvas that doesn't change.
		Prepares artists with empty drawing for future update at each frame.
		Cannot use certain stat funtions due to unavailable point.
		"""
		raise NotImplementedError
	def draw(self,t):
		"""Draws the changing part of given frame at given trajectory location.
		Returns list of artists that were drawn/changed.
		"""
		raise NotImplementedError

class panel_statscatter_tf(panel_base):
	def __init__(self,ax,d,pts,statx,staty,namet=None,g_ann=[],lims=None,scatterka=dict(),
				staty2=None,scatterka2=dict(),**ka):
		"""Draw scatter plots from two stats of TFs.
		statx,
		staty:		Stat instances on X and Y axes 
		namet:		TFs to show as a list of gene names. Defaults to all.
		g_ann:		TFs to annotate on scatter plot as list of gene names or {gene name:annotation text}. Use 'all' for all TFs.
		lims:		Limits of X and Y axes in [[xmin,xmax],[ymin,ymax]]. Defaults to min and max values on nodes with 2% expansion on each side.
		scatterka:	Keyword arguments for ax.scatter.
		staty2,
		scatterka2:	Respective parameters to draw a second y stat. Defaults to staty2=None to disable second y axis (on right side).
		ka:			Keyword arguments passed to parent class.
		"""
		import numpy as np
		super().__init__(ax,d,pts,**ka)
		#TF names and annotations
		t1=set(statx.names[0])&set(staty.names[0])
		if staty2 is not None:
			t1&=set(staty2.names[0])
		if namet is None:
			namet=sorted(list(t1))
		else:
			assert len(namet)>0
			t2=np.nonzero([x not in t1 for x in namet])[0]
			if len(t2)>0:
				raise ValueError('Genes in namet not found in stats: '+','.join([namet[x] for x in t2]))
		if g_ann=='all':
			g_ann=list(t1)
		self.namet=namet
		self.nametdict=dict(zip(self.namet,range(len(self.namet))))
		#Genes to annotate. Defaults to none
		t1=np.nonzero([x not in self.nametdict for x in g_ann])[0]
		if len(t1)>0:
			raise ValueError('TF(s) not found: {}'.format(','.join([g_ann[x] for x in t1])))
		if type(g_ann) is dict:
			self.g_ann={self.nametdict[x]:y for x,y in g_ann.items()}
		else:
			self.g_ann={self.nametdict[x]:x for x in g_ann}
		#Search stat functions
		self.stats=[statx,staty]
		if staty2 is not None:
			#Second Y
			self.stats.append(staty2)
		assert all([isinstance(x,stat_base) for x in self.stats])
		self.scatterka=scatterka
		self.lims=lims
		if staty2 is not None:
			#Second Y
			self.ny=2
			self.scatterka2=scatterka2
		else:
			self.ny=1
	def init(self):
		#Prepare limits
		if self.lims is None:
			self.lims=[None]*(self.ny+1)
		for xi in range(self.ny+1):
			if self.lims[xi] is not None:
				continue
			self.lims[xi]=self.stats[xi].default_lims(pts=self.pts,names=[self.namet])
		assert len(self.lims)==self.ny+1 and all([len(x)==2 for x in self.lims])
		#Draw initial panel
		ans=[]
		self.objs=[]
		self.objs.append(self.ax.scatter([],[],**self.scatterka))
		self.objs+=[self.ax.text(self.lims[0][0],self.lims[1][0],'') for x in self.g_ann]
		self.ax.set_xlim(self.lims[0])
		self.ax.set_ylim(self.lims[1])
		self.ax.set_xlabel(self.stats[0].label)
		self.ax.set_ylabel(self.stats[1].label)
		if self.ny>1:
			#Second Y
			self.ax2=self.ax.twinx()
			self.objs.append(self.ax2.scatter([],[],**self.scatterka2))
			self.objs+=[self.ax2.text(self.lims[0][0],self.lims[2][0],'') for x in self.g_ann]
			self.ax2.set_ylim(self.lims[2])
			self.ax2.set_ylabel(self.stats[2].label)
		else:
			if self.ax.spines['right'].get_visible():
				self.ax.tick_params(right=True,which='both')
		if self.ax.spines['top'].get_visible():
			self.ax.tick_params(top=True,which='both')
		return ans
	def get_data(self,pts):
		#n_stat,n_t,n_pts for data and isshow
		#Compute values at point from nodes
		import numpy as np
		data0=[self.stats[x].compute_points(pts) for x in range(self.ny+1)]
		assert all([data0[x].shape==(len(self.stats[x].names[0]),pts.npt) for x in range(self.ny+1)])
		data=[data0[x][[self.stats[x].ndict[0][y] for y in self.namet]] for x in range(self.ny+1)]
		assert all([x.shape==data[0].shape for x in data[1:]])
		data=np.array(data)

		#Compute isshow from data
		isshow=[self.stats[x].compute_show(pts,data0[x]) for x in range(self.ny+1)]
		assert all([isshow[x].shape==(len(self.stats[x].names[0]),pts.npt) for x in range(self.ny+1)])
		isshow=[isshow[x][[self.stats[x].ndict[0][y] for y in self.namet]] for x in range(self.ny+1)]
		assert all([x.shape==isshow[0].shape for x in isshow[1:]])
		isshow=np.array(isshow)
		isshow[1:]&=isshow[0]
		return [data,isshow]
	def draw(self,t):
		import numpy as np
		import itertools
		objs=[]
		pts=self.pts[[t]]
		data,isshow=[x[:,:,0] for x in self.get_data(pts)]

		#Update scatter
		for xi in range(self.ny):
			self.objs[xi*(len(self.g_ann)+1)].set_offsets(data[[0,1+xi]].T)
		#Update annotation
		t1=list(self.g_ann)
		for xi,xj in itertools.product(range(self.ny),range(len(t1))):
			if isshow[xi][t1[xj]]:
				self.objs[(1+len(self.g_ann))*xi+1+xj].set_text(self.g_ann[t1[xj]])
				self.objs[(1+len(self.g_ann))*xi+1+xj].set_position(data[[0,1+xi],t1[xj]])
			else:
				self.objs[(1+len(self.g_ann))*xi+1+xj].set_text('')
		objs=list(self.objs)
		return objs

class panel_statplot_tf(panel_statscatter_tf):
	def __init__(self,*a,g_ann=[],plotka=dict(),plotka2=dict(),**ka):
		"""Draw curve plots for two stats of TFs.
		Takes same parameters as panel_statscatter_tf."""
		super().__init__(*a,**ka)
		self.plotka=plotka
		self.plotka2=plotka2
	def init(self):
		import numpy as np
		#Prepare limits
		if self.lims is None:
			self.lims=[None]*(self.ny+1)
		for xi in range(self.ny+1):
			if self.lims[xi] is not None:
				continue
			self.lims[xi]=self.stats[xi].default_lims(pts=self.pts,names=[self.namet])
		assert len(self.lims)==self.ny+1 and all([len(x)==2 for x in self.lims])		
		#Prepare curve data
		objs=[]
		self.objs=[]
		data,isshow=self.get_data(self.pts)
		#Draw curves as initial panel and prepare pointers
		for xi in np.nonzero(isshow[1].any(axis=1))[0]:
			t1=np.nonzero(isshow[1,xi])[0]
			#Hide with masked array
			t2=np.ma.array(data[1,xi])
			t2[~isshow[1,xi]]=np.ma.masked
			self.ax.plot(data[0,xi],t2,**self.plotka)
			self.ax.text(data[0,xi,t1[-1]],data[1,xi,t1[-1]],self.namet[xi])
			self.objs.append(self.ax.scatter([],[],**self.scatterka))
		self.ax.set_xlabel(self.stats[0].label)
		self.ax.set_ylabel(self.stats[1].label)
		if self.ny>1:
			#Draw second Y
			self.ax2=self.ax.twinx()
			for xi in np.nonzero(isshow[2].any(axis=1))[0]:
				t1=np.nonzero(isshow[2,xi])[0]
				t2=np.ma.array(data[2,xi])
				t2[~isshow[2,xi]]=np.ma.masked
				self.ax2.plot(data[0,xi],t2,**self.plotka2)
				self.ax2.text(data[0,xi,t1[-1]],data[2,xi,t1[-1]],self.namet[xi])
				self.objs.append(self.ax2.scatter([],[],**self.scatterka2))
			self.ax2.set_ylabel(self.stats[2].label)
		else:
			if self.ax.spines['right'].get_visible():
				self.ax.tick_params(right=True,which='both')
		if self.ax.spines['top'].get_visible():
			self.ax.tick_params(top=True,which='both')
		return objs
	def draw(self,t):
		import numpy as np
		import itertools
		objs=[]
		pts=self.pts[[t]]
		data,isshow=[x[:,:,0] for x in self.get_data(pts)]
		#Update scatter
		for xi in range(self.ny):
			for xj in np.nonzero(isshow[xi+1])[0]:
				self.objs[xi*len(self.namet)+xj].set_offsets(data[[0,1+xi],[xj]].T)
				self.objs[xi*len(self.namet)+xj].set_alpha(1)
			for xj in np.nonzero(~isshow[xi+1])[0]:
				self.objs[xi*len(self.namet)+xj].set_alpha(0)
		objs=list(self.objs)
		return objs

class panel_stathist_tf(panel_base):
	def __init__(self,ax,d,pts,stat,namet=None,g_ann=[],lims=None,scatterka=dict(),**ka):
		"""Draw scatter plots from two stats of TFs.
		statx,
		staty:		Stat instances on X and Y axes 
		namet:		TFs to show as a list of gene names. Defaults to all.
		g_ann:		TFs to annotate on scatter plot as list of gene names or {gene name:annotation text}.
		lims:		Limits of X and Y axes in [[xmin,xmax],[ymin,ymax]]. Defaults to min and max values on nodes with 2% expansion on each side.
		scatterka:	Keyword arguments for ax.scatter.
		staty2,
		scatterka2:	Respective parameters to draw a second y stat. Defaults to staty2=None to disable second y axis (on right side).
		ka:			Keyword arguments passed to parent class.
		"""
		import numpy as np
		super().__init__(ax,d,pts,**ka)
		#TF names and annotations
		if namet is None:
			namet=stat.names[0]
		self.namet=namet
		self.nametdict=dict(zip(self.namet,range(len(self.namet))))
		#Genes to annotate. Defaults to none
		assert all([x in self.nametdict for x in g_ann])
		if type(g_ann) is dict:
			self.g_ann={self.nametdict[x]:y for x,y in g_ann.items()}
		else:
			self.g_ann={self.nametdict[x]:x for x in g_ann}
		#Search stat functions
		self.stat=stat
		assert isinstance(self.stat,stat_base)
		self.scatterka=scatterka
		self.lims=lims
	def init(self):
		#Prepare limits
		if self.lims is None:
			self.lims=[None]*(self.ny+1)
		for xi in range(self.ny+1):
			if self.lims[xi] is not None:
				continue
			self.lims[xi]=self.stat.default_lims(pts=self.pts,names=[self.namet])
		assert len(self.lims)==self.ny+1 and all([len(x)==2 for x in self.lims])
		#Draw initial panel
		ans=[]
		self.objs=[]
		self.objs.append(self.ax.scatter([],[],**self.scatterka))
		self.objs+=[self.ax.text(self.lims[0][0],self.lims[1][0],'') for x in self.g_ann]
		self.ax.set_xlim(self.lims[0])
		self.ax.set_ylim(self.lims[1])
		self.ax.set_xlabel(self.stats[0].label)
		self.ax.set_ylabel(self.stats[1].label)
		if self.ny>1:
			#Second Y
			self.ax2=self.ax.twinx()
			self.objs.append(self.ax2.scatter([],[],**self.scatterka2))
			self.objs+=[self.ax2.text(self.lims[0][0],self.lims[2][0],'') for x in self.g_ann]
			self.ax2.set_ylim(self.lims[2])
			self.ax2.set_ylabel(self.stats[2].label)
		return ans
	def get_data(self,pts):
		#n_stat,n_t,n_pts for data and isshow
		#Compute values at point from nodes
		import numpy as np
		data0=[self.stats[x].compute_points(pts) for x in range(self.ny+1)]
		assert all([data0[x].shape==(len(self.stats[x].names[0]),pts.npt) for x in range(self.ny+1)])
		data=[data0[x][[self.stats[x].ndict[0][y] for y in self.namet]] for x in range(self.ny+1)]
		assert all([x.shape==data[0].shape for x in data[1:]])
		data=np.array(data)

		#Compute isshow from data
		isshow=[self.stats[x].compute_show(pts,data0[x]) for x in range(self.ny+1)]
		assert all([isshow[x].shape==(len(self.stats[x].names[0]),pts.npt) for x in range(self.ny+1)])
		isshow=[isshow[x][[self.stats[x].ndict[0][y] for y in self.namet]] for x in range(self.ny+1)]
		assert all([x.shape==isshow[0].shape for x in isshow[1:]])
		isshow=np.array(isshow)
		isshow[1:]&=isshow[0]
		return [data,isshow]
	def draw(self,t):
		import numpy as np
		import itertools
		objs=[]
		pts=self.pts[[t]]
		data,isshow=[x[:,:,0] for x in self.get_data(pts)]

		#Update scatter
		for xi in range(self.ny):
			self.objs[xi*(len(self.g_ann)+1)].set_offsets(data[[0,1+xi]].T)
		#Update annotation
		t1=list(self.g_ann)
		for xi,xj in itertools.product(range(self.ny),range(len(t1))):
			if isshow[xi][t1[xj]]:
				self.objs[(1+len(self.g_ann))*xi+1+xj].set_text(self.g_ann[t1[xj]])
				self.objs[(1+len(self.g_ann))*xi+1+xj].set_position(data[[0,1+xi],t1[xj]])
			else:
				self.objs[(1+len(self.g_ann))*xi+1+xj].set_text('')
		objs=list(self.objs)
		return objs

class panel_cellscatter(panel_base):
	#HERE: need to update
	def __init__(self,ax,d,pts,loc='mean',skeleton=False,weightfunc=('conv',[1.],dict())):
		"""Draw scatter plots of cells
		loc:		How to compute pointer location.
			mean:	Takes weighted mean of cells involved.
			given:	Uses given pseudotime for linear interpolation
		skeleton:	Whether to show trajectory skeleton
		weightfunc:	Function that weight cells for color and pointer
		"""
		from lwang.pipeline.dynet2.traj import trajc,pointc
		super().__init__(ax,d,pts)
		#Compute coordinates
		traj0=trajc.fromdist(self.d['traj-edge'],d['traj-dist'])
		point=pointc(traj0,self.d['trajsub-edge'],self.d['trajsub-loc'])
		traj_coord=self.d['traj-coord']
		trajsub_coord=point.weight_linear()@traj_coord
		self.coord_node=trajsub_coord.T
		self.coord_cell=self.d['dtumap']
		self.loc=loc
		self.weightfunc=weightfunc
		self.cellweight=d['trajsub-set']
		self.skeleton=skeleton
	def init(self):
		ans=[]
		self.objs=[]
		#Cell scatter
		self.ax.set_aspect(1.)
		self.ax.axis('off')
		self.objs.append(self.ax.scatter(*self.coord_cell,marker='.',c='k',lw=0,s=5,alpha=0))
		if self.skeleton:
			#Node scatter
			self.ax.scatter(*self.coord_node,marker='o',c='k',lw=0,s=30)
			#Edges
			for xi in self.pts.p.edges:
				self.ax.plot(*self.coord_node[:,xi],'k-');
		#Pointer
		self.objs.append(self.ax.scatter([],[],marker='o',c=[(0.4,0.,0.)],lw=0,s=50))
		return ans
	def draw(self,t):
		import numpy as np
		from lwang.pipeline.dynet2.traj import pointc
		#Compute cell weight
		objs=[]
		w=self.compute_cellweight(self.pts.edges[t],self.pts.locs[t],self.weightfunc)
		w=w.ravel()/w.max()
		#Update cell color
		c=np.repeat([[1,0,0]],len(w),axis=0).T*w
		c+=np.repeat([np.ones(3)*0.9],len(w),axis=0).T*(1-w)
		c=np.clip(c,0,1)
		self.objs[0].set_color(c.T)
		self.objs[0].set_alpha(1.)
		#Update pointer location
		if self.loc=='mean':
			loc=(self.coord_cell@w)/w.sum()
		elif self.loc=='given':
			raise NotImplementedError
		else:
			raise ValueError('Unknown value loc={} during initializatin of class {}'.format(self.loc,self.__class__.name))
		self.objs[1].set_offsets(loc.reshape(1,2))
		objs=list(self.objs)
		return objs
	def compute_nodeweight_byloc(self,edge,loc,weightfunc):
		import numpy as np
		from lwang.pipeline.dynet2.traj import pointc
		#Compute values at point
		point=pointc(self.pts.p,np.array([edge]),np.array([loc]))
		return self.compute_nodeweight_bypoint(point,weightfunc)
	def compute_nodeweight_bypoint(self,point,weightfunc):
		wfunc=getattr(point,'weight_'+weightfunc[0])
		ans=wfunc(*weightfunc[1],**weightfunc[2])
		return ans
	def compute_cellweight(self,edge,loc,weightfunc):
		w=self.compute_nodeweight_byloc(edge,loc,weightfunc)
		return w@self.cellweight.T

class panel_netheatmap(panel_base):
	def __init__(self,ax,d,pts,smoothen_func=('conv',[1.],dict()),varname='dynet',namets=None,g_ann=[[],[]],lim=None,**ka):
		"""Draw dynamic heatmap for network.
		smoothen_func:	Function to smoothen network. See lwang.pipeline.dynet2.traj.trajc.smoothened
		varname:	Variable name to represent network
		namets:		Regulator and target gene lists to visualize as [reg list, target list]. Defaults to all.
		g_ann:		TFs to annotate on scatter plot as list of gene names or {gene name:annotation text}.
		ka:			Keyword arguments passed to parent class.
		"""
		import numpy as np
		assert len(g_ann)==2
		if any([len(x)>0 for x in g_ann]):
			raise NotImplementedError
		super().__init__(ax,d,pts)
		self.varname=varname
		namets0=[]
		dynet=self.d[self.varname]
		if dynet.ndim==5:
			dynet=dynet.reshape(dynet.shape[0],dynet.shape[2],dynet.shape[4])
			for xi in range(2):
				namets0.append(self.d.dim[self.d.d[self.varname].shape[2+2*xi]])
		else:
			assert dynet.ndim==3
			for xi in range(2):
				namets0.append(self.d.dim[self.d.d[self.varname].shape[1+xi]])
		if namets is None:
			namets=[None,None]
		for xi in range(2):
			if namets[xi] is not None:
				t1=dict(zip(namets0[xi],range(len(namets0[xi]))))
				t2=np.nonzero([x not in t1 for x in namets[xi]])[0]
				if len(t2)>0:
					raise ValueError(('TF' if xi==0 else 'Gene')+'(s) not found: '+','.join([namets[xi][x] for x in t2]))
				t1=[t1[x] for x in namets[xi]]
				dynet=dynet.swapaxes(0,xi+1)[t1].swapaxes(0,xi+1)
				namets0[xi]=namets[xi]
		self.namets=namets0
		self.nametdicts=[dict(zip(x,range(len(x)))) for x in self.namets]
		if lim is None:
			lim=[dynet.min(),dynet.max()]
			lim=np.abs(lim).max()
			lim=[-lim,lim]
		self.lim=lim
		self.dynet=self.pts.p.smoothened(dynet,smoothen_func[0],*smoothen_func[1],**smoothen_func[2],axis=0)
		self.ka=ka
		#Genes to annotate. Defaults to none
		# assert all([x in self.nametdict for x in g_ann])
		# if type(g_ann) is dict:
		# 	self.g_ann={self.nametdict[x]:y for x,y in g_ann.items()}
		# else:
		# 	self.g_ann={self.nametdict[x]:x for x in g_ann}
	def init(self):
		#Draw initial panel
		ans=[]
		self.objs=[]
		net=self.dynet(self.pts[[0]])[0]
		self.objs.append(self.ax.imshow(net,vmin=self.lim[0],vmax=self.lim[1],**self.ka))
		self.ax.set_xticks(list(range(len(self.namets[1]))))
		self.ax.set_xticklabels(self.namets[1],rotation=90)
		self.ax.set_yticks(list(range(len(self.namets[0]))))
		self.ax.set_yticklabels(self.namets[0])
		self.ax.set_xlabel('Target')
		self.ax.set_ylabel('TF')
		return ans
	def draw(self,t):
		net=self.dynet(self.pts[[t]])[0]
		self.objs[0].set_array(net)
		objs=list(self.objs)
		return objs

class panel_statheatmap(panel_base):
	def __init__(self,ax,d,pts,stat,annotation=[[],[]],lims=None,swapaxes=False,**ka):
		"""Draw dynamic heatmap for stat.
		stat:		Stat instance to draw. Must output two dimensions before the time dimension. 
		annotation:	Names of X and Y dimensions to draw dynamic heatmap for. Each element as a list of names or {name:annotation text}.
					Defaults to no annotation.
		lims:		Limits of values as [vmin,vmax] for coloring. Defaults to min and max values of all times.
		swapaxes:	By default, Dimension 0/1 are shown as Y/X respectively. Set swapaxes=True to swap X and Y axes.
		ka:			Keyword arguments passed to ax.imshow.
		"""
		import numpy as np
		assert len(g_ann)==2
		assert all([len(x)>0 for x in g_ann])
		if swapaxes:
			raise NotImplementedError('swapaxes is not implemented.')

		for xi in range(2):
			if len(g_ann[xi])!=len(set(g_ann[xi])):
				raise ValueError(f'Found duplicates in g_ann[{xi}].')
			t1=set(g_ann[xi])-set(stat.names[xi])
			if len(t1)>0:
				raise ValueError('Names not found in stat: '+','.join(t1))
		tid=[dict(x) for x in stat.names]
		tid=[[x[0][y] for y in x[1]] for x in zip(tid,g_ann)]
		self.tid=tid
		self.namet=namet
		self.nametdict=dict(zip(self.namet,range(len(self.namet))))
		#Genes to annotate. Defaults to none
		assert all([x in self.nametdict for x in g_ann])
		if type(g_ann) is dict:
			self.g_ann={self.nametdict[x]:y for x,y in g_ann.items()}
		else:
			self.g_ann={self.nametdict[x]:x for x in g_ann}
		#Search stat functions
		self.stats=[statx,staty]
		if staty2 is not None:
			#Second Y
			self.stats.append(staty2)
		assert all([isinstance(x,stat_base) for x in self.stats])
		self.scatterka=scatterka
		self.lims=lims
		if staty2 is not None:
			#Second Y
			self.ny=2
			self.scatterka2=scatterka2
		else:
			self.ny=1
			
		super().__init__(ax,d,pts)
		self.varname=varname
		namets0=[]
		dynet=self.d[self.varname]
		if dynet.ndim==5:
			dynet=dynet.reshape(dynet.shape[0],dynet.shape[2],dynet.shape[4])
			for xi in range(2):
				namets0.append(self.d.dim[self.d.d[self.varname].shape[2+2*xi]])
		else:
			assert dynet.ndim==3
			for xi in range(2):
				namets0.append(self.d.dim[self.d.d[self.varname].shape[1+xi]])
		if namets is None:
			namets=[None,None]
		for xi in range(2):
			if namets[xi] is not None:
				t1=dict(zip(namets0[xi],range(len(namets0[xi]))))
				t1=[t1[x] for x in namets[xi]]
				dynet=dynet.swapaxes(0,xi+1)[t1].swapaxes(0,xi+1)
				namets0[xi]=namets[xi]
		self.namets=namets0
		self.nametdicts=[dict(zip(x,range(len(x)))) for x in self.namets]
		if lim is None:
			lim=[dynet.min(),dynet.max()]
			lim=np.abs(lim).max()
			lim=[-lim,lim]
		self.lim=lim
		self.dynet=self.pts.p.smoothened(dynet,smoothen_func[0],*smoothen_func[1],**smoothen_func[2],axis=0)
		self.ka=ka
		#Genes to annotate. Defaults to none
		# assert all([x in self.nametdict for x in g_ann])
		# if type(g_ann) is dict:
		# 	self.g_ann={self.nametdict[x]:y for x,y in g_ann.items()}
		# else:
		# 	self.g_ann={self.nametdict[x]:x for x in g_ann}
	def init(self):
		#Draw initial panel
		ans=[]
		self.objs=[]
		net=self.dynet(self.pts[[0]])[0]
		self.objs.append(self.ax.imshow(net,vmin=self.lim[0],vmax=self.lim[1],**self.ka))
		self.ax.set_xticks(list(range(len(self.namets[1]))))
		self.ax.set_xticklabels(self.namets[1],rotation=90)
		self.ax.set_yticks(list(range(len(self.namets[0]))))
		self.ax.set_yticklabels(self.namets[0])
		self.ax.set_xlabel('Target')
		self.ax.set_ylabel('TF')
		return ans
	def draw(self,t):
		net=self.dynet(self.pts[[t]])[0]
		self.objs[0].set_array(net)
		objs=list(self.objs)
		return objs

class animate_generic:
	"""Generic class to visualize animations"""
	def __init__(self,pts,fig,panels):
		"""
		"""
		assert all([isinstance(x,panel_base) for x in panels])
		self.panels=panels
		self.fig=fig
		self.hasinit=False
		self.pts=pts
	def init(self):
		from contextlib import suppress
		if self.hasinit:
			return []
		objs=[]
		for xi0 in range(len(self.panels)):
			xi=self.panels[xi0]
			objs+=xi.init()
		self.hasinit=True
		return objs
	def draw(self,t):
		objs=[]
		for xi in self.panels:
			objs+=xi.draw(t)
		return objs
	def animate(self,**ka):
		"""Draw animation.
		ka:  Keyword arguments of matplotlib.animation.FuncAnimation.
		Return:
		Animation"""
		from matplotlib.animation import FuncAnimation
		return FuncAnimation(self.fig,self.draw,init_func=self.init,frames=self.pts.npt,blit=True,**ka)









































































assert __name__ != "__main__"
