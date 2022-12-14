#!/usr/bin/python3
# Lingfei Wang, 2022. All rights reserved.

"""
Dynamic network visualization layouts
"""

import abc
from typing import Optional,Tuple
import matplotlib
import dictys

class base(metaclass=abc.ABCMeta):
	"""
	Abstract base class for dynamic network plotting layout
	"""
	def __init__(self,dist:float=1.5,nframe:int=200,dpi:float=200):
		"""
		Base class to visualize dynamic network through layout

		Parameters
		----------
		dist:
			Kernel smoothing distance
		nframe:
			Number of frames (interpolated time points), use 100 or higher for finer resolution
		dpi:
			DPI for animation
		"""
		self.dist=dist
		self.nframe=nframe
		self.dpi=dpi
	@abc.abstractmethod
	def draw(self,net:dictys.net.network,branch:Tuple[int,int])->Tuple[dictys.traj.point,matplotlib.figure.Figure,'list[dictys.plot.panel.base]',dict]:
		"""
		Draws animation for dynamic network.
		
		Parameters
		----------
		net:
			Dynamic network to draw
		branch:
			Branch to draw defined as (start node id, end node id)
			
		Returns
		----------
		pts:
			List of artists that may be redrawn afterwards.
		fig:
			matplotlib figure of the animation
		panels:
			List of panels included in the animation
		animate_ka:
			Keyword arguments passed to dictys.plot.panel.animate_generic.animate
		"""

class notch(base):
	"""
	The notch layout is like an old iphone.

	+---+       +---+
	| a |       | b |
	+---+---+---+---+
	| c | d | e | f |
	+---+---+---+---+
	| c | d | e | f |
	+---+---+---+---+
	| c | d | e | f |
	+---+---+---+---+
	| c | d | e | f |
	+---+---+---+---+
	....
	
	a and b appear only once. Each row for c,d,e,f has their own set of select TFs.
	a: Dynamic tracking of cells used for GRN inference
	b: Dynamic scatter plot for differential regulation v.s. differential expression logFCs
	c: Dynamic plot for expression level (log CPM) of select TFs as a function of pseudo-time
	d: Dynamic plot for regulatory activity (log target count) of select TFs as a function of pseudo-time
	e: Dynamic heatmap for regulation strength from select TFs to select target genes
	f: Dynamic subnetwork graph from select TF to its targets
	"""
	def draw(self,net:dictys.net.network,branch:Tuple[int,int],
		weighted:bool=False,sparsity:float=0.01,panelsize:Tuple[float,float]=(6,4),
		subplots_adjust:Optional[dict]={'wspace':0.35,'hspace':0.2},
		bcde_tfs:list[list[str]]=[],
		e_targets:list[str]=[],
		f_weight_func:Tuple[str,list,dict]=('linear',[],{}),f_tfs:list[list[str]]=[],f_stop:int=20,f_niteration:int=50,
		a_ka:dict={},b_ka:dict={},c_ka:dict={},d_ka:dict={},e_ka:dict={},f_ka:dict={},
		)->Tuple[dictys.traj.point,matplotlib.figure.Figure,'list[dictys.plot.panel.base]',dict]:
		"""
		Draws animation for dynamic network with the notch layout.
		
		Parameters
		----------
		net:
			Dynamic network to draw
		branch:
			Branch to draw defined as (start node id, end node id)
		weighted:
			Whether to use weighted instead of binarized network for panels b,d.
		sparsity:
			Overall configuration of binarized network sparsity (or equilvalent for weighted network). Affects panels b,d,f.
		panelsize:
			Size of each panel
		subplots_adjust:
			Keyword arguments for matplotlib.pyplot.subplots_adjust. If None, subplots_adjust function is skipped.
		bcde_tfs:
			Select TFs to annotate in b and in c,d,e of each row
		e_targets:
			List of target genes to annotate in all rows
		f_tfs:
			Select TFs to draw target genes in f of each row
		f_weight_func:
			Smoothing parameters for subnetwork layout coordinates. Should be left unchanged.
		f_niteration:
			Number of iterations to (re)perform subnet force-directed layout at each frame.
		f_stop:
			Early stopping iteration count for subnet force-directed layout. Smaller values than f_niteration lead to insufficient convergence.
		a_ka:
			Keyword arguments passed to a panel (dictys.plot.panel.cellscatter).
		b_ka:
			Keyword arguments passed to b panel (dictys.plot.panel.statscatter).
		c_ka:
			Keyword arguments passed to c panels (dictys.plot.panel.statplot).
		d_ka:
			Keyword arguments passed to d panels (dictys.plot.panel.statplot).
		e_ka:
			Keyword arguments passed to e panels (dictys.plot.panel.statheatmap).
		f_ka:
			Keyword arguments passed to f panels (dictys.plot.panel.network).
			
		Returns
		----------
		pts:
			List of artists that may be redrawn afterwards.
		fig:
			matplotlib figure of the animation
		panels:
			List of panels included in the animation
		animate_ka:
			Keyword arguments passed to dictys.plot.panel.animate_generic.animate
		"""
		import itertools
		from os import linesep
		from functools import partial
		import matplotlib.pyplot as plt
		from dictys.net.layout import _fruchterman_reingold
		from dictys.net import stat
		from dictys.plot import panel
		
		#Keyword arguments
		a_ka_default=dict({})
		b_ka_default=dict({'aspect':1,'lim':{'sym','min','max'}})
		c_ka_default=dict({})
		d_ka_default=dict({})
		e_ka_default=dict({'cmap':'coolwarm'})
		f_ka_default=dict({'nodeka':{'scatterka':{'s':5,'lw':0}},'edgeka':{'lw':0.05}})
		a_ka_default.update(a_ka)
		b_ka_default.update(b_ka)
		c_ka_default.update(c_ka)
		d_ka_default.update(d_ka)
		e_ka_default.update(e_ka)
		f_ka_default.update(f_ka)
		
		n=len(bcde_tfs)
		if len(f_tfs)!=n:
			raise ValueError('bcde_tfs and f_tfs should have the same length as the number of rows.')
		panelcount=(4,1+n)

		########################################
		# Prepare dynamic network properties
		########################################

		# Kernel smoothed network properties
		pts,fsmooth=net.linspace(branch[0],branch[1],self.nframe,self.dist)
		# Expression, quantified with logCPM
		stat1_lcpm=fsmooth(stat.lcpm(net,cut=0))
		# Kernel smoothed network
		stat1_net=fsmooth(stat.net(net))
		# Binarized network
		stat1_netbin=stat.fbinarize(stat1_net,sparsity=sparsity)
		# You can change network sparsity (proportion of positive edges) with:
		# stat1_netbin=stat.fbinarize(stat1_net,sparsity=0.001)
		# Regulatory activity, quantified with log target count
		if weighted:
			stat1_lntarget=stat.flnneighbor(stat1_net,weighted_sparsity=sparsity)
		else:
			stat1_lntarget=stat.flnneighbor(stat1_netbin)
		# Pseudo time
		stat1_pseudotime=stat.pseudotime(net,pts)

		# Kernel smoothed network properties for each row or panel
		b_tfs=list(set(itertools.chain.from_iterable(bcde_tfs)))
		# Selecting TF's outgoing edges as subnetwork
		stat1_subnets=[stat1_net[x] for x in f_tfs]
		stat1_subnetbins=[stat1_netbin[x] for x in f_tfs]
		stat1_subnet_truncs=[stat.function(lambda *y:y[0]*y[1],x,names=x[0].names) for x in zip(stat1_subnets,stat1_subnetbins)]
		# Performing layout with linear smoothing of node locations
		stat1_layouts=[stat.fsmooth(stat.flayout_base(x,partial(_fruchterman_reingold,stop=f_stop,iterations=f_niteration),pts=pts),pts,f_weight_func) for x in stat1_subnet_truncs]

		########################################
		# Draw animation
		########################################

		# Animation formating
		fig=plt.figure(figsize=(panelsize[0]*panelcount[0],panelsize[1]*panelcount[1]),dpi=self.dpi)
		axes=[fig.add_subplot(*panelcount[::-1],x+1) for x in range(panelcount[0]*panelcount[1])]
		for ax in axes[:3]:
			for x in ax.spines.values():
				x.set_visible(False)
			ax.set_xticks([])
			ax.set_yticks([])
		for ax in axes[3:]:
			ax.spines['top'].set_visible(False)
			ax.spines['right'].set_visible(False)
		if subplots_adjust is not None:
			plt.subplots_adjust(**subplots_adjust)
			
		def post_init(_):
			# Animation formating function after init
			for xi in itertools.product(range(1,n),range(3)):
				ax=axes[4*xi[0]+xi[1]]
				ax.set_xlabel('')
				ax.xaxis.set_major_formatter('')
			for xi in range(1,n+1):
				axes[4*xi].set_ylabel('LogCPM')
				axes[4*xi+1].set_ylabel('Log target count')
			return []

		# Draw each panel with iterator
		panels=[]
		axes_iter=iter(axes)
		# Panel a
		panels.append(panel.cellscatter(next(axes_iter),net,pts,fsmooth,**a_ka_default))
		# Create two empty panels to skip
		next(axes_iter)
		next(axes_iter)
		# Panel b
		panels.append(panel.statscatter(next(axes_iter),pts,stat.fdiff(stat1_lcpm,stat.finitial(stat1_lcpm,pts),label=f'Differential expression{linesep}logFC in CPM'),stat.fdiff(stat1_lntarget,stat.finitial(stat1_lntarget,pts),label=f'Differential regulation{linesep}logFC in target count'),annotate=b_tfs,**b_ka_default))
		# Each row
		for xi in range(n):
			# Panel c
			panels.append(panel.statplot(next(axes_iter),pts,stat1_pseudotime,stat1_lcpm,names=bcde_tfs[xi],**c_ka_default))
			# Panel d
			panels.append(panel.statplot(next(axes_iter),pts,stat1_pseudotime,stat1_lntarget,names=bcde_tfs[xi],**d_ka_default))
			# Panel e
			panels.append(panel.statheatmap(next(axes_iter),pts,stat1_net,names=[bcde_tfs[xi],e_targets],**e_ka_default))
			# Panel f
			f_ka_default['nodeka']=dict(f_ka_default['nodeka'] if 'nodeka' in f_ka_default else {})
			f_ka_default['nodeka']['annotate']=f_tfs[xi]
			panels.append(panel.network(next(axes_iter),pts,stat1_layouts[xi],stat1_subnet_truncs[xi],**f_ka_default))
		return (pts,fig,panels,{'post_init':post_init})


#
