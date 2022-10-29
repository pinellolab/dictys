#!/usr/bin/python3
# Lingfei Wang, 2022. All rights reserved.

"""
Visualization panels of networks
"""

import abc
from typing import Union,Optional,Callable,Any,Tuple
import dictys
from dictys.net import stat

class base(metaclass=abc.ABCMeta):
	"""
	Abstract base class for dynamic network plotting panel
	"""
	def __init__(self,ax,pts:dictys.traj.point):
		"""
		Base class to visualize single panel for dynamic network

		Parameters
		----------
		ax:		matplotlib.pyplot.axes
			Axes to draw on
		pts:
			Points of path to visualize network
		"""
		self.ax=ax
		self.pts=pts
	@abc.abstractmethod
	def init(self)->list:
		"""
		Draws initial canvas that doesn't change.
		Prepares artists with empty drawing for future update at each frame.
		Cannot use certain stat funtions due to unavailable point.

		Returns
		----------
		list
			List of artists that may be redrawn afterwards.
		"""
	@abc.abstractmethod
	def draw(self,t:int)->list:
		"""Draws the changing part of given frame at given trajectory location.

		Parameters
		----------
		t:	int
			Frame index to draw.

		Returns
		----------
		list
			List of artists that are redrawn.
		"""

class overlay(base):
	"""
	Overlay class that draws several panel classes in the same panel.
	"""
	def __init__(self,ax,pts:dictys.traj.point,panels:list[base]):
		"""
		Overlay class that draws several panel classes in the same panel.

		Parameters
		----------
		ax:		matplotlib.pyplot.axes
			Axes to draw on
		pts:
			Points of path to visualize network
		panels:
			Panel objects to draw on the same panel
		"""
		assert len(panels)>=1
		assert all(isinstance(x,base) for x in panels)
		self.panels=panels
		super().__init__(ax,pts)
	def init(self):
		import itertools
		objs=[x.init() for x in self.panels]
		return list(itertools.chain.from_iterable(objs))
	def draw(self,t:int,**ka):
		import itertools
		objs=[x.draw(t,**ka) for x in self.panels]
		return list(itertools.chain.from_iterable(objs))

class statscatter(base):
	"""
	Draw scatter plots from two stats.
	"""
	def __init__(self,ax,pts:dictys.traj.point,statx:stat.base,staty:stat.base,names:Optional[list[str]]=None,annotate:Union[str,list[str],dict[str,str]]=[],
			annotate_err:bool=True,
			lim:Union[list[Union[tuple[float,float]]],set[str],None]=set(),scatterka:dict[str,Any]={},statka:dict[str,Any]={},
			staty2:Optional[stat.base]=None,scatterka2:dict[str,Any]={},aspect:Optional[float]=None):
		"""
		Draw scatter plots from two stats.

		Parameters
		----------
		ax:				matplotlib.pyplot.axes
			Axes to draw on
		pts:
			Points of path to visualize network
		statx:
			Stat instance for X axis. Must be one-dimensional.
		staty:
			Stat instance for Y axis. Must be one-dimensional.
		names:
			Names to show. Defaults to all.
		annotate:
			Names to annotate on scatter plot. Use 'all' for all names. Use dict to rename annotation as {name:annotation}.
		annotate_err:
			Whether to report error if any node to annotate is not found.
		lim:
			Limits of X and Y axes. Takes several formats.
			* [[xmin,xmax],[ymin,ymax],[y2min,y2max] if y2 is present]. This exactly specifies lims. Each [min,max] can be None to be determined automatically with `default_lims` method of each stat.

			* Set containing none, some, or all of the options for automatically determining limits based on `default_lims` method of each stat:

				- 'sym': Each axis use symmetric [-x,x] limits. x is determined by the max of absolute values of limits.

				- 'max': All axes use the max of max limits of all axes. If unset, all axes use their respective max limits.

				- 'min': All axes use the min of min limits of all axes. If unset, all axes use their respective min limits.

		scatterka:
			Keyword arguments for ax.scatter.
		statka:
			Keyword arguments for ax.scatter whose values are stats. Keys must be settable in matplotlib.collections.PathCollection.set.
		staty2:
			Stat instance for a second Y axis on right side. Must be one-dimensional. Defaults to staty2=None to disable.
		scatterka2:
			Keyword arguments for ax.scatter for second Y.
		"""
		import numpy as np
		from functools import reduce
		from operator import and_
		super().__init__(ax,pts)
		if not all(isinstance(x,stat.base) for x in statka.values()):
			raise TypeError('All values of statka must be a stat.')
		if any(x is not None and len(x.names)!=1 for x in [statx,staty,staty2]):
			raise ValueError('Stat locations must be 1-dimensional.')
		if any(len(x.names)<1 for x in statka.values()):
			raise ValueError('Stat parameters must be at least 1-dimensional.')
		#Dot names and annotations
		t1=set(statx.names[0])&set(staty.names[0])
		if staty2 is not None:
			t1&=set(staty2.names[0])
		if len(statka)>0:
			t1&=reduce(and_,[set(x.names[0]) for x in statka.values()])
		if names is None:
			names=sorted(list(t1))
		else:
			assert len(names)>0
			t2=np.nonzero([x not in t1 for x in names])[0]
			if len(t2)>0:
				raise ValueError('Not found in at least one of the stats: '+','.join([names[x] for x in t2]))
		assert len(names)>0
		if annotate=='all':
			annotate=list(t1)
		self.names=names
		self.namesdict=dict(zip(self.names,range(len(self.names))))
		#Dots to annotate. Defaults to none
		t1=np.nonzero([x not in self.namesdict for x in annotate])[0]
		if len(t1)>0:
			if annotate_err:
				raise ValueError('TF(s) not found: {}'.format(','.join([annotate[x] for x in t1])))
			t1=list(filter(lambda x:x in self.namesdict,annotate))
			annotate={x:annotate[x] for x in t1} if isinstance(annotate,dict) else t1
		if isinstance(annotate,dict):
			self.annotate={self.namesdict[x]:y for x,y in annotate.items()}
		else:
			self.annotate={self.namesdict[x]:x for x in annotate}
		#Search stat functions
		self.stats=[statx,staty]
		if staty2 is not None:
			#Second Y
			self.stats.append(staty2)
		assert all(isinstance(x,stat.base) for x in self.stats)
		self.scatterka=scatterka
		self.statka=statka
		self.aspect=aspect
		if staty2 is not None:
			#Second Y
			self.ny=2
			self.scatterka2=scatterka2
		else:
			self.ny=1
		if isinstance(lim,list) and len(lim)!=self.ny+1:
			raise ValueError('lim size ({}) must match the size of all axes ({}).'.format(len(lim),self.ny+1))
		self.lim=lim
	def autolimits(self):
		"""
		Determines limits of axes

		Returns
		-------
		numpy.ndarray(shape=(self.ny+1,2))
			Automatically determined limits of each axis
		"""
		import numpy as np
		#Prepare autolimits
		lim=[self.stats[x].default_lims(pts=self.pts,names=[self.names]) for x in range(self.ny+1)]
		assert len(lim)==self.ny+1 and all(len(x)==2 for x in lim)
		for xi in range(len(lim)):		# pylint: disable=C0200
			if np.isnan(lim[xi]).any() or lim[xi][0]==lim[xi][1]:
				lim[xi]=[0,1]
		lim=np.array(lim)
		#Replace/update autolimits based on user option
		if isinstance(self.lim,set):
			t1=self.lim-{'sym','min','max'}
			if len(t1)>0:
				raise ValueError('Unrecognized options for lim: '+','.join(t1))
			if 'sym' in self.lim:
				lim=np.repeat(np.abs(lim).max(axis=1).reshape(-1,1),2,axis=1)
				lim[:,0]*=-1
			lim=np.array(lim)
			if 'min' in self.lim:
				lim[:,0]=lim[:,0].min()
			if 'max' in self.lim:
				lim[:,1]=lim[:,1].max()
		else:
			for xi in range(self.ny+1):
				if self.lim[xi] is None:
					continue
				for xj in range(2):
					if self.lim[xi][xj] is not None:
						lim[xi,xj]=self.lim[xi,xj]
		return lim		
	def init(self)->list:
		lim=self.autolimits()
		self.objs=[]
		#Draw initial panel
		self.objs.append(self.ax.scatter([],[],**self.scatterka))
		self.objs+=[self.ax.text(lim[0,0],lim[1,0],'') for x in self.annotate]
		self.ax.set_xlim(lim[0])
		self.ax.set_ylim(lim[1])
		self.ax.set_xlabel(self.stats[0].label)
		self.ax.set_ylabel(self.stats[1].label)
		if self.ny>1:
			#Second Y
			self.ax2=self.ax.twinx()
			self.objs.append(self.ax2.scatter([],[],**self.scatterka2))
			self.objs+=[self.ax2.text(lim[0,0],lim[2,0],'') for x in self.annotate]
			self.ax2.set_ylim(lim[2])
			self.ax2.set_ylabel(self.stats[2].label)
		else:
			if self.ax.spines['right'].get_visible():
				self.ax.tick_params(right=True,which='both')
		if self.ax.spines['top'].get_visible():
			self.ax.tick_params(top=True,which='both')
		if self.aspect is not None:
			self.ax.set_aspect(self.aspect)
		ans=self.draw(0,force=True)
		return ans
	def get_data(self,pts:dictys.traj.point,force:bool=False):
		"""
		Computes data needed for drawing.

		Parameters
		----------
		pts:
			Points to get data for
		force:
			Whether to force getting data even if the stats are constant over time.

		Returns
		----------
		data:		numpy.ndarray(shape=(n_stat,n_scatter,len(pts))) or None
			Data for scatter locations. Only None if all stats are constant over time.
		param:		{str:numpy.ndarray(shape=(n_scatter,len(pts)))}
			Data for parameters. Only non-constant stats are included.
		"""
		#Compute values at point from nodes
		import numpy as np
		if force or not all(hasattr(x,'isconst') and x.isconst for x in self.stats):
			data=[self.stats[x].compute(pts) for x in range(self.ny+1)]
			assert all(data[x].shape==(len(self.stats[x].names[0]),len(pts)) for x in range(self.ny+1))
			data=[data[x][[self.stats[x].ndict[0][y] for y in self.names]] for x in range(self.ny+1)]
			assert all(x.shape==data[0].shape for x in data[1:])
			data=np.array(data)
		else:
			data=None

		#Compute stat based parameters
		param={x:y.compute(pts) for x,y in self.statka.items() if force or not (hasattr(y,'isconst') and y.isconst)}
		assert all(param[x].shape==tuple([len(y) for y in self.statka[x].names]+[len(pts)]) for x in param)
		param={x:param[x][[self.statka[x].ndict[0][y] for y in self.names]] for x in param}
		if data is not None:
			assert all(x.shape[0]==data[0].shape[0] for x in param.values())
		return [data,param]
	def draw(self,t:int,force:bool=False):		# pylint: disable=W0221
		"""Draws the changing part of given frame at given trajectory location.

		Parameters
		----------
		t:		int
			Frame index to draw.
		force:	bool
			Whether to force redraw even if nothing changed.

		Returns
		----------
		list
			List of artists that are redrawn.
		"""
		import numpy as np
		import itertools
		objs=[]
		pts=self.pts[[t]]
		data,param=self.get_data(pts,force=force)
		param={x:np.take(y,0,axis=-1) for x,y in param.items()}
		if data is not None:
			data=data[:,:,0]
			isshow=np.isfinite(data)
			isshow=isshow[1:]&isshow[0]
		
		#Update scatter
		for xi in range(self.ny):
			if data is not None:
				self.objs[xi*(len(self.annotate)+1)].set_offsets(data[[0,1+xi]].T)
			self.objs[xi*(len(self.annotate)+1)].set(**param)
			if data is not None or len(param)>0:
				objs.append(self.objs[xi*(len(self.annotate)+1)])
		#Update annotation
		if data is not None:
			t1=list(self.annotate)
			for xi,xj in itertools.product(range(self.ny),range(len(t1))):
				if isshow[xi][t1[xj]]:
					self.objs[(1+len(self.annotate))*xi+1+xj].set_text(self.annotate[t1[xj]])
					self.objs[(1+len(self.annotate))*xi+1+xj].set_position(data[[0,1+xi],t1[xj]])
				else:
					self.objs[(1+len(self.annotate))*xi+1+xj].set_text('')
				objs.append(self.objs[(1+len(self.annotate))*xi+1+xj])
		return objs

class statplot_static(statscatter):
	"""
	Draw constant line plots from two stats. This is a static plot and no redraw takes place.
	"""
	def __init__(self,ax,pts:dictys.traj.point,statx:stat.base,staty:stat.base,colors,plotka:dict[str,Any]={},plotka2:dict[str,Any]={},**ka):
		"""
		Draw constant line plots from two stats. This is a static plot and no redraw takes place.

		Parameters
		----------
		ax:			matplotlib.pyplot.axes
			Axes to draw on
		pts:
			Points of path to visualize network
		statx:
			Stat instance for X axis. Must be one-dimensional. Each entry is a line to draw.
		staty:
			Stat instance for Y axis. Must be one-dimensional. Each entry is a line to draw.
		colors:		list
			List of colors for each line in matplotlib format.
		plotka:
			Keyword arguments for ax.plot.
		plotka2:
			Keyword arguments for ax.plot for second Y.
		ka:
			Keyword arguments in the same format as in dictys.plot.panel.statscatter.
		"""
		super().__init__(ax,pts,statx,staty,**ka)
		self.plotka=plotka
		self.plotka2=plotka2
		self.colors=colors
	def init(self):
		import numpy as np
		# lim=self.autolimits()
		#Prepare curve data
		ans=[]
		data,_=self.get_data(self.pts,force=True)
		isshow=np.isfinite(data)
		isshow=isshow[1:]&isshow[0]
		#Draw curves as initial panel and prepare pointers
		for xi in np.nonzero(isshow[0].any(axis=1))[0]:
			t1=np.nonzero(isshow[0,xi])[0]
			#Hide with masked array
			t2=np.ma.array(data[1,xi])
			t2[~isshow[0,xi]]=np.ma.masked
			ans.append(self.ax.plot(data[0,xi],t2,c=self.colors[xi],**self.plotka))
			ans.append(self.ax.text(data[0,xi,t1[-1]],data[1,xi,t1[-1]],self.names[xi]))
		self.ax.set_xlabel(self.stats[0].label)
		self.ax.set_ylabel(self.stats[1].label)
		if self.ny>1:
			#Draw second Y
			self.ax2=self.ax.twinx()
			for xi in np.nonzero(isshow[1].any(axis=1))[0]:
				t1=np.nonzero(isshow[1,xi])[0]
				t2=np.ma.array(data[2,xi])
				t2[~isshow[1,xi]]=np.ma.masked
				ans.append(self.ax2.plot(data[0,xi],t2,c=self.colors[xi],**self.plotka2))
				ans.append(self.ax2.text(data[0,xi,t1[-1]],data[2,xi,t1[-1]],self.names[xi]))
			self.ax2.set_ylabel(self.stats[2].label)
		else:
			if self.ax.spines['right'].get_visible():
				self.ax.tick_params(right=True,which='both')
		if self.ax.spines['top'].get_visible():
			self.ax.tick_params(top=True,which='both')
		return []
	def draw(self,*a,**ka):		# pylint: disable=W0613
		return []

class statplot(overlay):
	"""
	Draw line plots from two stats by overlaying a line plot with a scatter point for pointer.
	"""
	def __init__(self,ax,pts:dictys.traj.point,statx:stat.base,staty:stat.base,names:Optional[list[str]],cmap='tab10',pointer:bool=True,plotka:dict[str,Any]={},pointerka:dict[str,Any]={}):
		"""
		Draw line plots from two stats by overlaying a line plot with a scatter point for pointer.

		Parameters
		----------
		ax:			matplotlib.pyplot.axes
			Axes to draw on
		pts:
			Points of path to visualize network
		statx:
			Stat instance for X axis. Must be one-dimensional.
		staty:
			Stat instance for Y axis. Must be one-dimensional.
		names:
			Names to show. Defaults to all.
		cmap:		str or numpy.ndarray(shape=(len(names),3 or 4))
			Colormap in matplotlib format for different names
		pointer:
			Whether to show the current scatter plot as pointers
		plotka:
			Keyword arguments for ax.plot.
		pointerka:
			Keyword arguments for ax.scatter for pointers.
		"""
		import numpy as np
		from dictys.plot import get_cmap
		
		#Keyword arguments
		pointerka_default=dict({'s':20,'zorder':99})
		pointerka_default.update(pointerka)
		#Determine colors for cell types
		if isinstance(cmap,str):
			cmap=get_cmap(cmap,len(names))
		cmap=np.array(cmap)
		panels=[]
		#Plot
		panels.append(statplot_static(ax,pts,statx,staty,cmap,names=names,plotka=plotka))
		if pointer:
			#Color stat
			statc=stat.const(cmap,[names,list('RGBA')])
			panels.append(statscatter(ax,pts,statx,staty,names=names,scatterka=pointerka_default,statka={'color':statc}))
		return super().__init__(ax,pts,panels[::-1])

class cellscatter_scatter(statscatter):
	"""
	Draw scatter plots of cells.
	"""
	def __init__(self,ax,d:dictys.net.network,pts:dictys.traj.point,statx:stat.base,staty:stat.base,statw:stat.base,
		cmap='tab10',color:str='type',alphas:Tuple[float,float]=(0.05,0.5),legend_loc:Tuple[float,float]=(1.1,1),
		legend_ka:dict[str,Any]={},aspect:float=1,**ka):
		"""
		Draw scatter plots of cells.

		Parameters
		----------
		ax:			matplotlib.pyplot.axes
			Axes to draw on
		d:
			Dynamic network
		pts:
			Points of path to visualize network
		statx:
			Stat instance for X axis with shape=(n_cell)
		staty:
			Stat instance for Y axis with shape=(n_cell)
		statw:
			Stat instance for cell weights with shape=(n_cell)
		cmap:		str or numpy.ndarray(shape=(len(names),3 or 4))
			Colormap in matplotlib format for different names
		color:
			Cell property name to use as categorical coloring
		alphas:
			Alpha/transparency for cell weight=0 and 1 respectively.
		legend_loc:
			Relative location of cell type legend
		legend_ka:
			Keyword arguments for drawing cell type legend with dictys.plot.colorlegend
		aspect:
			Aspect ratio
		ka:
			Keyword arguments scatterka for drawing cells with dictys.plot.panel.statscatter
		"""
		import numpy as np
		from dictys.plot import get_cmap
		#Determine colors for cell types
		if 'c' in d.prop and color in d.prop['c']:
			ctype=d.prop['c'][color]
		else:
			ctype=['Unknown']*d.cn
		ctypelist=sorted(list(set(ctype)))
		nctype=len(ctypelist)
		if isinstance(cmap,dict):
			assert all(x in cmap for x in ctypelist)
		elif nctype>1:
			cmap=dict(zip(ctypelist,get_cmap(cmap,nctype)))
		else:
			cmap={ctypelist[0]:[0.5]*3+[1]}
		self.cmap=cmap
		colors=np.array([cmap[x] for x in ctype])
		assert colors.shape==(len(ctype),4) and colors.min()>=0 and colors.max()<=1
		#Color stat
		statc=stat.const(colors,[d.cname,list('RGBA')])
		#Alpha stat
		stata=stat.function(lambda *x:alphas[0]+(alphas[1]-alphas[0])*x[0],[statw],names=[d.cname])
		self.d=d
		self.legend_loc=legend_loc
		self.legend_ka=legend_ka
		super().__init__(ax,pts,statx,staty,statka={'color':statc,'alpha':stata},scatterka=ka,aspect=aspect)
	def init(self):
		from dictys.plot import colorlegend
		t1=list(self.cmap)
		self.ax.axis('off')
		colorlegend(self.ax,self.legend_loc,t1,[self.cmap[x] for x in t1],**self.legend_ka)
		return super().init()

class cellscatter_pointer(statscatter):
	"""
	Draw pointer for average cell location.
	"""
	def __init__(self,ax,pts:dictys.traj.point,statx:stat.base,staty:stat.base,statw:stat.base,**ka):
		"""
		Draw pointer for average cell location.

		Parameters
		----------
		ax:			matplotlib.pyplot.axes
			Axes to draw on
		pts:
			Points of path to visualize network
		statx:
			Stat instance for X axis with shape=(n_cell)
		staty:
			Stat instance for Y axis with shape=(n_cell)
		statw:
			Stat instance for cell weights with shape=(n_cell)
		ka:
			Keyword arguments scatterka for drawing average cell pointer with dictys.plot.panel.statscatter
		"""
		statx,staty=[stat.function(lambda *x:((x[0]*x[1]).sum(axis=0)/x[1].sum(axis=0)).reshape(1,-1),[y,statw],names=[['pointer']]) for y in [statx,staty]]
		super().__init__(ax,pts,statx,staty,scatterka=ka)

class cellscatter(overlay):
	"""
	Draws overlay plot of scatter plots of cells and average cell pointer.
	"""
	def __init__(self,ax,d:dictys.net.network,pts:dictys.traj.point,fsmooth,pointer:bool=True,scatterka:dict[str,Any]={},pointerka:dict[str,Any]={}):
		"""
		Draws overlay plot of scatter plots of cells and average cell pointer.

		Parameters
		----------
		ax:			matplotlib.pyplot.axes
			Axes to draw on
		d:			
			Dynamic network object to draw cells
		pts:
			Points of path to visualize network
		fsmooth:	functools.partial
			Partial function that produces Gaussian kernel smoothened statistic on an original statistic as its parameter
		pointer:
			Whether to show the current scatter plot as pointers
		scatterka:
			Keyword arguments for ax.scatter for cells.
		pointerka:
			Keyword arguments for ax.scatter for pointers.
		"""
		#Keyword arguments
		scatterka_default=dict({'s':2,'lw':0})
		pointerka_default=dict({'color':'k','s':20,'zorder':99})
		scatterka_default.update(scatterka)
		pointerka_default.update(pointerka)
		#X&Y coordindates
		statx=stat.const(d.prop['c']['coord'][0],[d.cname],label='Dim1')
		staty=stat.const(d.prop['c']['coord'][1],[d.cname],label='Dim2')
		#Cell weight
		statw=fsmooth(stat.sprop(d,'sc','w'))
		statw=stat.function(lambda *x:(x[0]/x[0].max(axis=0)),[statw],names=[d.cname])
		panels=[]
		#Scatter
		panels.append(cellscatter_scatter(ax,d,pts,statx,staty,statw,**scatterka_default))
		if pointer:
			panels.append(cellscatter_pointer(ax,pts,statx,staty,statw,**pointerka_default))
		return super().__init__(ax,pts,panels[::-1])

class statheatmap(base):
	"""
	Draw dynamic heatmap for a single stat.
	"""
	def __init__(self,ax,pts:dictys.traj.point,stat1:stat.base,names:Optional[Tuple[Optional[list[str]],Optional[list[str]]]]=None,
		annotate:Optional[Tuple[Optional[list[str]],Optional[list[str]]]]=None,lim:Optional[Tuple[float,float]]=None,cmap_sym:bool=True,aspect='auto',**ka):
		"""
		Draw dynamic heatmap for a single stat.

		Parameters
		----------
		ax:			matplotlib.pyplot.axes
			Axes to draw on
		pts:
			Points of path to visualize network
		stat1:
			Stat instance to draw. Must be two-dimensional.
		names:
			Names to show. Defaults to all.
		annotate:
			Names to annotate. Defaults to all.
		lim:
			Limits in [min,xmax] for coloring. Defaults to min and max values.
		cmap_sym:
			Whether to use symmetric lim if lim is unspecified.
		aspect:
			Aspect ratio
		ka:
			Keyword arguments for ax.imshow.
		"""
		import numpy as np
		super().__init__(ax,pts)
		self.stat=stat1
		#Rows & columns to show
		if names is None:
			names=[None,None]
		assert len(stat1.names)==2
		for xi in range(2):
			if names[xi] is None:
				names[xi]=list(stat1.names[xi])
		t1=[np.nonzero([x not in stat1.ndict[y] for x in names[y]])[0] for y in range(2)]
		if len(t1[0])>0:
			raise ValueError('Regulator(s) not found: {}'.format(','.join([names[0][x] for x in t1[0]])))
		if len(t1[1])>0:
			raise ValueError('Target(s) not found: {}'.format(','.join([names[1][x] for x in t1[1]])))
		self.names=names
		self.namesdict=[dict(zip(x,range(len(x)))) for x in self.names]
		#Rows & columns to annotate
		if annotate is None:
			annotate=[None,None]
		assert len(annotate)==2
		annotate=list(annotate)
		for xi in range(2):
			if annotate[xi] is None:
				annotate[xi]=list(names[xi])
		t1=[np.nonzero([x not in self.namesdict[y] for x in names[y]])[0] for y in range(2)]
		if len(t1[0])>0:
			raise ValueError('Regulator(s) to annotate not found: {}'.format(','.join([annotate[0][x] for x in t1[0]])))
		if len(t1[1])>0:
			raise ValueError('Target(s) to annotate not found: {}'.format(','.join([annotate[1][x] for x in t1[1]])))
		self.annotate=annotate
		if lim is None:
			lim=stat1.default_lims(pts)
			if cmap_sym:
				lim=np.abs(lim).max()
				lim=[-lim,lim]
		self.lim=list(lim)
		ka['aspect']=aspect
		self.ka=ka
	def get_data(self,pts:dictys.traj.point):
		"""
		Obtains stat needed for heatmap

		Parameters
		----------
		pts:
			A single point to compute stat for

		Returns
		-------
		numpy.ndarray
			Stat values at the requested point
		"""
		assert len(pts)==1
		data=self.stat.compute(pts)
		data=data[[self.stat.ndict[0][x] for x in self.names[0]]][:,[self.stat.ndict[1][x] for x in self.names[1]]]
		return data
	def init(self):
		#Draw initial panel
		self.objs=[]
		data=self.get_data(self.pts[[0]])
		self.objs.append(self.ax.imshow(data,vmin=self.lim[0],vmax=self.lim[1],**self.ka))
		t1=[[x[0][y] for y in x[1]] for x in zip(self.namesdict,self.annotate)]
		self.ax.set_xticks(t1[1])
		self.ax.set_xticklabels(self.annotate[1],rotation=90)
		self.ax.set_yticks(t1[0])
		self.ax.set_yticklabels(self.annotate[0])
		self.ax.set_xlabel('Target')
		self.ax.set_ylabel('Regulator')
		return self.objs
	def draw(self,t:int):
		data=self.get_data(self.pts[[t]])
		self.objs[0].set_array(data)
		objs=list(self.objs)
		return objs

class network_node(statscatter):
	"""
	Draw network nodes with scatter plot.
	"""
	def __init__(self,ax,pts:dictys.traj.point,statloc:stat.base,*a,aspect:float=1,**ka):
		"""
		Draw network nodes with scatter plot.

		Parameters
		----------
		ax:			matplotlib.pyplot.axes
			Axes to draw on
		pts:
			Points of path to visualize network
		statloc:
			Stat instance for axes. Must be two-dimensional with shape (n_cell,2).
		a:			list
			Arguments for dictys.plot.panel.statscatter
		aspect:
			Aspect ratio
		ka:			dict
			Keyword arguments for dictys.plot.panel.statscatter
		"""
		statx=statloc[:,0]
		staty=statloc[:,1]
		super().__init__(ax,pts,statx,staty,*a,aspect=aspect,**ka)

class network_edge(base):
	"""
	Draw network edges.
	"""
	def __init__(self,ax,pts:dictys.traj.point,statloc:stat.base,statnet:stat.base,**ka):
		"""
		Draw network edges.

		Parameters
		----------
		ax:			matplotlib.pyplot.axes
			Axes to draw on
		pts:
			Points of path to visualize network
		statloc:
			Stat instance for axes. Must be two-dimensional with shape (n_cell,2).
		statnet:
			Stat instance for edge strengths. Must be two-dimensional with shape (n_reg,n_target).
		ka:
			Keyword arguments for ax.arrow
		"""
		self.statnet=statnet
		self.statloc=statloc
		self.map=[{y:self.statloc.ndict[0][self.statnet.names[x][y]] for y in range(len(self.statnet.names[x])) if self.statnet.names[x][y] in self.statloc.ndict[0]} for x in range(2)]
		self.ka=ka
		super().__init__(ax,pts)
	def init(self):
		self.objs=[]
		self.nlast=0
		return self.objs
	def draw(self,t:int):
		import numpy as np
		pts=self.pts[[t]]
		net,loc=[np.take(x.compute(pts),0,axis=-1) for x in [self.statnet,self.statloc]]
		net=np.abs(net)
		t1=np.nonzero(net)
		for xi in range(len(t1[0])):
			t3=[self.map[x][t1[x][xi]] for x in range(2)]
			if xi<len(self.objs):
				self.objs[xi].set_data(x=loc[t3[0],0],y=loc[t3[0],1],dx=loc[t3[1],0]-loc[t3[0],0],dy=loc[t3[1],1]-loc[t3[0],1])
				self.objs[xi].set_alpha(1)
			else:
				self.objs.append(self.ax.arrow(loc[t3[0],0],loc[t3[0],1],loc[t3[1],0]-loc[t3[0],0],loc[t3[1],1]-loc[t3[0],1],**self.ka))
		#Use last frame arrow count to determin objects for redraw
		if len(t1[0])<self.nlast:
			for xi in self.objs[len(t1[0]):self.nlast]:
				xi.set_alpha(0)
			t2=self.nlast
			self.nlast=len(t1[0])
			return self.objs[:t2]
		self.nlast=len(t1[0])
		return self.objs[:self.nlast]

class network(overlay):
	"""
	Draws overlay plot of network.
	"""
	def __init__(self,ax,pts:dictys.traj.point,statloc:stat.base,statnet:stat.base,nodeka:dict[str,Any]={},edgeka:dict[str,Any]={}):
		"""
		Draws overlay plot of network.

		Parameters
		----------
		ax:			matplotlib.pyplot.axes
			Axes to draw on
		pts:
			Points of path to visualize network
		statloc:
			Stat instance for node coordindates. Must be two-dimensional with shape (n_node,2).
		statnet:
			Stat instance for edge strengths. Must be two-dimensional with shape (n_reg,n_target).
		nodeka:
			Keyword arguments for drawing nodes with dictys.plot.panel.network_node
		edgeka:
			Keyword arguments for drawing edges with dictys.plot.panel.network_edge
		"""
		panels=[]
		#Scatter
		nodeka0={'scatterka':{}}
		if 'scatterka' in nodeka:
			nodeka=dict(nodeka)
			nodeka0['scatterka'].update(nodeka['scatterka'])
			del nodeka['scatterka']
		nodeka0.update(nodeka)
		panels.append(network_node(ax,pts,statloc,**nodeka0))
		panels.append(network_edge(ax,pts,statloc,statnet,**edgeka))
		return super().__init__(ax,pts,panels[::-1])
	def init(self):
		self.ax.axis('off')
		return super().init()

class animate_generic:
	"""
	Engine class to visualize animations
	"""
	def __init__(self,pts:dictys.traj.point,fig,panels:list[base]):
		"""Animation engine.

		Parameters
		----------
		pts:
			Points for drawing animation
		fig:	matplotlib.pyplot.Figure
			Figure for drawing animation
		panels:
			Panels to animate
		"""
		assert all(isinstance(x,base) for x in panels)
		self.panels=panels
		self.fig=fig
		self.hasinit=False
		self.pts=pts
	def init(self,post_init:Callable[[Any],list]=lambda _:[])->list:
		"""
		Initialize canvas.

		Returns
		----------
		list of artists
			Artists that may be redrawn in the future frames
		"""
		if self.hasinit:
			return []
		objs=[]
		for xi in self.panels:
			objs+=xi.init()
		self.hasinit=True
		objs+=post_init(self)
		return objs
	def draw(self,t:int)->list:
		"""
		Update canvas at frame t.

		Parameters
		----------
		t:
			Frame ID
		
		Returns
		----------
		list of artists
			Artists redrawn
		"""
		objs=[]
		for xi in self.panels:
			objs+=xi.draw(t)
		return objs
	def animate(self,blit:bool=True,post_init:Callable[[Any],list]=lambda _:[],**ka):
		"""
		Draw animation.

		Parameters
		----------
		blit:
			Whether to redraw only the needed parts.
		post_init:
			Function to call after initialization.
		ka: 
			Keyword arguments passed to matplotlib.animation.FuncAnimation.

		Returns
		----------
		matplotlib.animation.FuncAnimation
			Animation
		"""
		from functools import partial
		from matplotlib.animation import FuncAnimation
		return FuncAnimation(self.fig,self.draw,init_func=partial(self.init,post_init=post_init),frames=len(self.pts),blit=blit,**ka)































assert __name__ != "__main__"
