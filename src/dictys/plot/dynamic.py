#!/usr/bin/python3
# Lingfei Wang, 2022. All rights reserved.

"""
Dynamic network visualizations
"""

from typing import Union,Tuple,Optional
import pandas as pd
import dictys
from dictys.utils.importing import matplotlib
from dictys.utils.numpy import NDArray

########################################################################
# TF discovery plot
########################################################################

def auc(dx:NDArray[float],dy:NDArray[float])->NDArray[float]:
	"""
	Computes area under the curves.

	Parameters
	----------
	dx:	numpy.ndarray(shape=(n,))
		X coordinates 
	dy:	numpy.ndarray(shape=(ny,n))
		Y coordinates, one for each y curve

	Returns
	-------
	numpy.ndarray(shape=(ny,))
		Area under the curves
	"""
	if len(dx)<2 or not (dx[1:]>dx[:-1]).all():
		raise ValueError('dx must be increasing and have at least 2 values.')
	dxdiff=dx[1:]-dx[:-1]
	dymean=(dy[:,1:]+dy[:,:-1])/2
	ans=dymean@dxdiff
	return ans

def _dynamic_network_char_terminal_logfc_(dx:NDArray[float],dy:NDArray[float])->NDArray[float]:
	"""
	Computes terminal logFC for curves.

	Parameters
	----------
	dx:	numpy.ndarray(shape=(n,))
		X coordinates 
	dy:	numpy.ndarray(shape=(ny,n))
		Y coordinates, one for each y curve

	Returns
	-------
	numpy.ndarray(shape=(ny,))
		Terminal logFCs
	"""
	if len(dx)<2 or not (dx[1:]>dx[:-1]).all():
		raise ValueError('dx must be increasing and have at least 2 values.')
	return dy[:,-1]-dy[:,0]

def _dynamic_network_char_transient_logfc_(dx:NDArray[float],dy:NDArray[float])->NDArray[float]:
	"""
	Computes transient logFC for curves.

	Parameters
	----------
	dx:	numpy.ndarray(shape=(n,))
		X coordinates 
	dy:	numpy.ndarray(shape=(ny,n))
		Y coordinates, one for each y curve

	Returns
	-------
	numpy.ndarray(shape=(ny,))
		Transient logFCs
	"""
	import numpy as np
	n=dy.shape[1]
	dx=(dx-dx[0])/(dx[-1]-dx[0])
	dy=dy-np.median([dy,np.repeat(dy[:,[0]],n,axis=1),np.repeat(dy[:,[-1]],n,axis=1)],axis=0)
	return auc(dx,dy)

def _dynamic_network_char_switching_time_(dx:NDArray[float],dy:NDArray[float])->NDArray[float]:
	"""
	Computes switching time for curves.

	Parameters
	----------
	dx:	numpy.ndarray(shape=(n,))
		X coordinates 
	dy:	numpy.ndarray(shape=(ny,n))
		Y coordinates, one for each y curve

	Returns
	-------
	numpy.ndarray(shape=(ny,))
		Switching time
	"""
	import numpy as np
	n=dy.shape[1]
	dx=(dx-dx[0])/(dx[-1]-dx[0])
	dy=np.median([dy,np.repeat(dy[:,[0]],n,axis=1),np.repeat(dy[:,[-1]],n,axis=1)],axis=0)
	return (auc(dx,(dy.T-dy[:,-1]).T))/(dy[:,0]-dy[:,-1]+1E-300)

def _compute_chars_(dcurve:pd.DataFrame,dtime:pd.DataFrame)->pd.DataFrame:
	charlist={
		'Terminal logFC':_dynamic_network_char_terminal_logfc_,
		'Transient logFC':_dynamic_network_char_transient_logfc_,
		'Switching time':_dynamic_network_char_switching_time_,
	}
	#Compute curve characteristics
	dchar={}
	for xj in charlist:
		dchar[xj]=charlist[xj](dtime.values,dcurve.values)
	dchar=pd.DataFrame.from_dict(dchar)
	dchar.set_index(dcurve.index,inplace=True,drop=True)
	return dchar

def _transform_inset_(left:float,bottom:float,width:float,height:float,figsize:float)->Tuple[float,float,float,float]:
	return (left/figsize[0],bottom/figsize[1],width/figsize[0],height/figsize[1])

def draw_discover1(dcurve:pd.DataFrame,dchar:pd.DataFrame,dtime:Optional[pd.Series]=None,
		cmap:Union[str,matplotlib.cm.ScalarMappable,dict[str,Union[str,matplotlib.cm.ScalarMappable]]]='viridis',
		curve_expand:Tuple[float,float]=(0,0.1),inset_sides:int=2,inset_lvs:int=1,inset_size:Tuple[float,float]=(0.55,0.3),inset_space0:float=0.15,inset_space:Tuple[float,float]=(0.05,0),
		heatmap_height:float=0.1,heatmap_space:float=0.05,line_space:float=0.04,fs:float=8,
		ka_curve:dict={},ka_heatmap:dict={},
	)->Tuple[matplotlib.figure.Figure,list,dict[str,matplotlib.cm.ScalarMappable]]:		# noqa: E123
	"""
	Draws a single TF discovery plot for dynamic network from given characteristic dataframes.

	Parameters
	----------
	dcurve:
		Dataframe for curves with each TF as a row and each time point as a column in chronological order
	dchar:
		Dataframe for curve characteristics with each TF as a row and each characteristic as a column
	dtime:
		Series for pseudo-time values matching columns of dcurve. Defaults to equal time interval.
	cmap:
		Color map for each characteristic. Accepts the following types:

		* str or matplotlib.cm.ScalarMappable: Same colormap for all characteristics

		* dict(str=str or matplotlib.cm.ScalarMappable): Separate colormaps for different characteristics

		For str, matplotlib names for colormap are used and normalizations are [min,max] by default.

	curve_expand:
		Relative expansion of axis limits for each curve on negative and positive sides.
	inset_sides:
		Binary: which side to put curve insets.

		* 1: below

		* 2: above

	inset_lvs:
		Number of levels of curve inset on each side. Accepts: 1,2.
	inset_size:
		Size of each curve inset in inch
	inset_space0:
		Space between curve inset and heatmap in inch
	inset_space:
		Space between curve inset themselves in inch
	heatmap_height:
		Height of each heatmap in inch
	heatmap_space:
		Space between heatmaps in inch
	line_space:
		Distance leaving blank for indication lines in inch
	fs:
		Font size
	ka_curve:
		Keyword arguments for drawing curves
	ka_heatmap:
		Keyword arguments for drawing heatmaps
	
	Returns
	-------
	fig:
		Figure object
	axes:
		List of axes objects including axes for main figure, heatmaps, and curves
	cmap:
		Dictionary of colormap by characteristic name used
	"""
	import matplotlib.pyplot as plt
	import numpy as np
	from mpl_toolkits.axes_grid1.inset_locator import inset_axes
	ka_curve_default=dict({'lw':1,'c':'k'})
	ka_heatmap_default=dict({'aspect':'auto'})
	ka_curve_default.update(ka_curve)
	ka_heatmap_default.update(ka_heatmap)

	#Match rows & columns
	t1=set(dcurve.index)&set(dchar.index)
	if len(t1)==0:
		raise ValueError('No overlapping index found for dcurve and dchar.')
	if dcurve.shape[1]<2:
		raise ValueError('At least two time points needed for dcurve.')
	dcurve,dchar=[x[x.index.isin(t1)] for x in [dcurve,dchar]]
	if (dcurve.index!=dchar.index).any():
		raise ValueError('Unmatching indices for dcurve and dchar.')
	assert np.isfinite(dcurve.values).all()
	assert np.isfinite(dchar.values).all()
	if dtime is None:
		dtime=np.linspace(0,1,dcurve.shape[1])
	else:
		if dtime.values.shape!=(dcurve.shape[1],):
			raise TypeError('dtime must be one dimensional with matching shape with dcurve.')
		if frozenset(dtime.index)!=frozenset(dcurve.columns):
			raise ValueError("dtime's indices must match dcurve's columns.")
		dtime=dtime.loc[dcurve.columns].values
		t1=dtime.argsort()
		dtime=dtime[t1]
		dcurve=dcurve[dcurve.columns[t1]]
	#Initialize color map
	if isinstance(cmap,(str,matplotlib.cm.ScalarMappable)):
		cmap={x:cmap for x in dchar.columns}
	if not isinstance(cmap,dict):
		raise TypeError('Unaccepted type for parameter cmap.')
	cmap={x:y if isinstance(y,matplotlib.cm.ScalarMappable) else matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=dchar[x].min(),vmax=dchar[x].max()),cmap=y) for x,y in cmap.items() if x in dchar.columns}

	#Draw figure
	assert line_space*2<inset_space0
	assert inset_lvs in {1,2}
	n=len(dcurve)
	inset_nside=((inset_sides&1!=0)+(inset_sides&2!=0))
	inset_nrow=inset_lvs*inset_nside
	assert n%inset_nrow==0
	inset_ncol=n//inset_nrow
	heatmap_width=(inset_size[0]+inset_space[0])/inset_nrow
	heatmap_fullwidth=heatmap_width*n
	
	ns=dchar.shape[1]
	figsize=(
		(inset_size[0]+inset_space[0])*(inset_ncol+1-1/inset_nrow),
		heatmap_height*ns+inset_size[1]*inset_nrow+heatmap_space*(ns-1)+(inset_space0+inset_space[1]*(inset_lvs-1))*inset_nside,
	)
	fig=plt.figure(figsize=figsize)
	ax0=fig.add_subplot(111)
	ax0.set_xlim([0,1])
	ax0.set_ylim([0,1])
	
	axes=[ax0]
	#Heatmaps
	if inset_sides==3:
		heatmap_bbox=(figsize[1]+heatmap_height*ns+heatmap_space*(ns-1))/2
	elif inset_sides==2:
		heatmap_bbox=0
	else:
		raise NotImplementedError('inset_sides not in {2,3}')
	heatmap_bbox=np.array([
		(figsize[0]-heatmap_fullwidth)/2,
		heatmap_bbox,
		heatmap_fullwidth,
		heatmap_height*ns+heatmap_space*(ns-1)])
	for xi in range(ns):
		column=dchar.columns[ns-xi-1]
		inset_fullsize=(heatmap_fullwidth,heatmap_height)
		inset_bbox=(heatmap_bbox[0],heatmap_bbox[1]-heatmap_bbox[3]+heatmap_height*xi+heatmap_space*xi)
		bbox=_transform_inset_(*inset_bbox,*inset_fullsize,figsize)
		ax=inset_axes(ax0,"100%","100%",loc='upper left',bbox_to_anchor=bbox,borderpad=0,bbox_transform=ax0.transAxes)
		v1=dchar[column].values
		ax.imshow(cmap[column].to_rgba(v1.reshape(1,-1)),**ka_heatmap_default)
		for xj in ax.spines.values():
			xj.set_visible(False)
		ax.set_xticks([])
		ax.set_yticks([0])
		ax.set_yticklabels([column],fontsize=fs)
		axes.append(ax)
		
	#Curves
	#Base locations for each row
	locx0=heatmap_bbox[0]+(np.arange(inset_nrow)+0.5)*heatmap_width-inset_size[0]/2
	if inset_sides==3:
		isbelow=np.array([1,0,1,0],dtype=bool)
		locy0=np.array([0,0,-1,1])
	elif inset_sides==2:
		isbelow=np.array([0,0],dtype=bool)
		locy0=np.array([0,1])
	else:
		raise NotImplementedError('inset_sides not in {2,3}')
	locy0=(1-isbelow*2)*(heatmap_bbox[3]/2+inset_space0+inset_size[1]/2)+locy0*(inset_size[1]+inset_space[1])
	locy0+=heatmap_bbox[1]-heatmap_bbox[3]/2-inset_size[1]/2
	locy0=locy0[:inset_nrow]
	#Add locations for each row
	locx1=np.ones(inset_nrow)*heatmap_width*inset_nrow
	locy1=np.zeros(inset_nrow)
	vrange=[dcurve.values.min(),dcurve.values.max()]
	for xi in range(n):
		#Column ID
		locx=xi//inset_nrow
		#Row ID
		locy=xi%inset_nrow
		inset_fullsize=inset_size
		inset_bbox=(locx0[locy]+locx1[locy]*locx,locy0[locy]+locy1[locy]*locx)
		bbox=_transform_inset_(*inset_bbox,*inset_fullsize,figsize)
		ax=inset_axes(ax0,"100%","100%",loc='upper left',bbox_to_anchor=bbox,borderpad=0,bbox_transform=ax0.transAxes)
		axes.append(ax)
		ax.plot(dtime,(dcurve.loc[dcurve.index[xi]].values-vrange[0])/(vrange[1]-vrange[0]),**ka_curve_default)
		ax.set_xlim([dtime.min(),dtime.max()])
		ax.set_ylim([-curve_expand[0],1+curve_expand[1]])
		for xj in ['top','right']:
			ax.spines[xj].set_visible(False)
		ax.set_xticks([])
		ax.set_yticks([])
		ax.set_xlabel(dcurve.index[xi],fontsize=fs)
		if not isbelow[locy]:
			ax.xaxis.set_label_position('top')
		#Draw line
		src=bbox[1]+isbelow[locy]*bbox[3]
		dst=heatmap_bbox[1]-isbelow[locy]*heatmap_bbox[3]
		dst=dst/figsize[1]
		ax0.plot([bbox[0]+bbox[2]/2]*2,[src+np.sign(dst-src)*line_space/figsize[1],dst-np.sign(dst-src)*line_space/figsize[1]],'-',lw=1,color=[0.5]*3)
	ax0.axis('off')
	return fig,axes,cmap

def fig_discover(dcurve:pd.DataFrame,dtime:pd.DataFrame,ntops:Tuple[int,int,int,int]=(8,8,4,4),
		vrange:dict[str,Union[str,Tuple[float,float]]]={
			'Terminal logFC':'symmetric',
			'Switching time':'extreme',
			'Transient logFC':'symmetric',
		},
		cmap:dict[str,str]={
			'Terminal logFC':'coolwarm',
			'Switching time':'viridis',
			'Transient logFC':'PiYG_r',
		},
		**ka)->list[Tuple[matplotlib.figure.Figure,list,dict[str,matplotlib.cm.ScalarMappable]]]:
	"""
	Draws TF discovery plots for different patterns of variation for dynamic network.

	Parameters
	----------
	dcurve:
		Dataframe for curves with each TF as a row and each time point as a column in chronological order
	dchar:
		Dataframe for curve characteristics with each TF as a row and each characteristic as a column
	ntops:
		Number of top TFs to show for activating, inactivating, transient up, and transient down patterns separately. Set to 0 to skip drawing a particular pattern.
	vrange:
		Value range dictionary for each characteristic. Supports the following options for each value:

		* (vmin,vmax): Manual specification

		* 'extreme': Automatically determined using extreme (min and max) values among those shown

		* 'symmetric':  Automatically determined using maximum absolute values among those shown

	cmap:
		Colormap dictionary for each characteristic. Values should be matplotlib colormap name.
	ka:
		Keyword arguments for dictys.plot.dynamic.draw_dynamic1

	Returns
	-------
	Tuple of returns from dictys.plot.dynamic.draw_dynamic1 for each of the patterns drawn
	"""
	import itertools
	import numpy as np

	#Compute curve characteristics
	dchar=_compute_chars_(dcurve,dtime)
	#Convert vrange to dict
	if not isinstance(vrange,dict):
		vrange={x:vrange for x in dchar}

	#Parameters for 4 plots
	ans=[]
	if ntops[0]>0:
		#Activating
		t1=dchar.sort_values('Terminal logFC',ascending=False).index[:ntops[0]]
		t1=dchar.loc[t1].sort_values('Switching time').index
		ans.append([[dcurve.loc[t1],dchar.loc[t1][['Switching time','Terminal logFC']]],dict(dtime=dtime,cmap=dict(cmap),**ka)])
	if ntops[1]>0:
		#Inactivating
		t1=dchar.sort_values('Terminal logFC',ascending=True).index[:ntops[1]]
		t1=dchar.loc[t1].sort_values('Switching time').index
		ans.append([[dcurve.loc[t1],dchar.loc[t1][['Switching time','Terminal logFC']]],dict(dtime=dtime,cmap=dict(cmap),**ka)])
	if ntops[2]>0:
		#Transient up
		t1=dchar.sort_values('Transient logFC',ascending=False).index[:ntops[2]]
		ans.append([[dcurve.loc[t1],dchar.loc[t1][['Transient logFC','Terminal logFC']]],dict(dtime=dtime,cmap=dict(cmap),**ka)])
	if ntops[3]>0:
		#Transient down
		t1=dchar.sort_values('Transient logFC',ascending=True).index[:ntops[3]]
		ans.append([[dcurve.loc[t1],dchar.loc[t1][['Transient logFC','Terminal logFC']]],dict(dtime=dtime,cmap=dict(cmap),**ka)])
	#Determining colormap normalization using rows plotted.
	for xi in vrange:
		t1=np.nonzero([xi in x[0][1].columns for x in ans])[0]
		if isinstance(vrange[xi],str):
			#Automatically determine vrange
			t2=list(set(itertools.chain.from_iterable([ans[x][0][1].index for x in t1])))
			if len(t2)==0:
				continue
			t2=dchar.loc[t2][xi].values
			if vrange[xi]=='extreme':
				vrange1=[t2.min(),t2.max()]
			elif vrange[xi]=='symmetric':
				vrange1=np.abs(t2).max()
				vrange1=[-vrange1,vrange1]
			else:
				raise ValueError(f'Unknown vrange value: {vrange[xi]}')
		else:
			vrange1=vrange[xi]
		for xj in t1:
			#Update colormap with normalization
			ans[xj][1]['cmap'].update({xi:matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=vrange1[0],vmax=vrange1[1]),cmap=ans[xj][1]['cmap'][xi])})

	ans=[draw_discover1(*x[0],**x[1]) for x in ans]
	return ans

########################################################################
# Individual regulation heatmap
########################################################################

def fig_regulation_heatmap(network:dictys.net.dynamic_network,start:int,stop:int,regulations:list[Tuple[str,str]],num:int=100,dist:float=1.5,ax:Optional[matplotlib.axes.Axes]=None,cmap:Union[str,matplotlib.cm.ScalarMappable]='coolwarm',figsize:Tuple[float,float]=(2,0.22),vmax:Optional[float]=None)->Tuple[matplotlib.pyplot.Figure,matplotlib.axes.Axes,matplotlib.cm.ScalarMappable]:
	"""Draws pseudo-time dependent heatmap of individual regulation strengths.
	
	Parameters
	----------
	network:
		Dynamic network object
	start:
		Branch starting node ID 
	stop:
		Branch ending node ID 
	regulations:
		List of regulations in (Regulator name, target name) format to draw strength
	num:
		Number of points from starting to ending nodes to draw
	dist:
		Gaussian kernel smoothing distance/length scale between cells. Larger value means stronger smoothing.
	ax:
		Axes to draw on.
	figsize:
		Figure size for each regulation as a row. Should remain unassigned when ax is assigned.
	cmap:
		Colormap in matplotlib string or matplotlib.cm.ScalarMappable
	vmax:
		Maximum value in colormap. Should remain unassigned when cmap is a matplotlib.cm.ScalarMappable instance.
	
	Returns
	-------
	fig:
		Heatmap figure
	ax:
		Heatmap axes
	cmap:
		Heatmap colormap
	"""
	import numpy as np
	import matplotlib.pyplot as plt
	from dictys.net import stat
	#Get dynamic network edge strength
	pts,fsmooth=network.linspace(start,stop,num,dist)
	stat1_net=fsmooth(stat.net(network))
	# stat1_x=stat.pseudotime(network,pts)
	#Test regulation existence
	tdict=[dict(zip(x,range(len(x)))) for x in stat1_net.names]
	t1=[[x[y] for x in regulations if x[y] not in tdict[y]] for y in range(2)]
	if len(t1[0])>0:
		raise ValueError('Regulator gene(s) {} not found in network.'.format('/'.join(t1[0])))
	if len(t1[1])>0:
		raise ValueError('Target gene(s) {} not found in network.'.format('/'.join(t1[1])))
	#Extract regulations to draw
	dnet=stat1_net.compute(pts)
	t1=np.array([[tdict[0][x[0]],tdict[1][x[1]]] for x in regulations]).T
	dnet=dnet[t1[0],t1[1]]
	#Create figure and axes
	if ax is None:
		figsize=(figsize[0],figsize[1]*dnet.shape[0])
		fig=plt.figure(figsize=figsize)
		ax=fig.add_subplot(111)
	else:
		if figsize is not None:
			raise ValueError('figsize should not be set if ax is set.')
		fig=ax.get_figure()
		figsize=fig.get_size_inches()
	aspect=(figsize[1]/dnet.shape[0])/(figsize[0]/dnet.shape[1])
	#Determine colormap
	if isinstance(cmap,str):
		if vmax is None:
			vmax=np.quantile(np.abs(dnet).ravel(),0.95)
		cmap=matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=-vmax,vmax=vmax),cmap=cmap)		# pylint: disable=E1130
	elif vmax is not None:
		raise ValueError('vmax should not be set if cmap is a matplotlib.cm.ScalarMappable.')
	#Draw figure
	ax.imshow(cmap.to_rgba(dnet),aspect=aspect,interpolation='none')
	ax.set_xticks([])
	ax.set_yticks(list(range(len(regulations))))
	ax.set_yticklabels(['-'.join(x) for x in regulations])
	return fig,ax,cmap




#
