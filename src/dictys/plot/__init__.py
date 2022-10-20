#!/usr/bin/python3
# Lingfei Wang, 2020-2022. All rights reserved.

"""
Network visualization
"""

__all__=['dynamic','layout','panel','static']

from typing import Tuple,Optional,Union
from dictys.utils.importing import matplotlib
from . import *

def get_cmap(cmap,n):
	"""
	Get discrete colors from colormap in matplotlib

	Parameters
	----------
	cmap:	str
		Matplotlib color map name
	n:		int
		Number of colors needed

	Return
	------
	numpy.ndarray(shape=(n,4))
		n colors in RGBA format
	"""
	import numpy as np
	import matplotlib.pyplot as plt
	c=plt.get_cmap(cmap)
	if n<c.N:
		return c(np.linspace(0,1,c.N)[:n])
	return c(np.linspace(0,1,n))

def colorbar(cmap:Union[str,matplotlib.cm.ScalarMappable],vmin:Optional[float]=None,vmax:Optional[float]=None,figsize:Optional[Tuple[float,float]]=(0.15,0.8),ax=None,orientation:str='vertical',title:Optional[str]=None,**ka):
	"""
	Draw simple colorbar
	
	Parameters
	----------
	cmap:
		matplotlib colormap
	vmin:
		Minimum value for colormap. Should remain unassigned when cmap is a matplotlib.cm.ScalarMappable.
	vmax:
		Maximum value for colormap. Should remain unassigned when cmap is a matplotlib.cm.ScalarMappable.
	figsize:
		Size of colorbar figure in inches. Conflicts with ax.
	ax:
		Axes to draw on. Conflicts with figsize.
	orientation:
		Colobar orientation
	title:
		Title of colorbar
	ka:
		Keyword arguments passed through, to colorbar function by default and to set_title if with prefix 'title_'
		
	Returns
	-------
	fig:
		Figure of colorbar
	ax:
		Axes of colorbar
	"""
	import matplotlib.pyplot as plt
	if isinstance(cmap,matplotlib.cm.ScalarMappable):
		if vmax is not None or vmin is not None:
			raise TypeError('vmin and vmax should be unset when cmap is a matplotlib.cm.ScalarMappable')
	else:
		cmap=matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=vmin, vmax=vmax),cmap=cmap)

	ka_title={x[len('title_'):]:y for x,y in ka.items() if x.startswith('title_')}
	ka={x:y for x,y in ka.items() if not x.startswith('title_')}

	if ax is None:
		fig,ax=plt.subplots(figsize=figsize)
	else:
		if figsize is not None:
			raise ValueError('At most one of figsize and ax should be set')
		fig=ax.get_figure()
	fig.colorbar(cmap,cax=ax,orientation=orientation,**ka)
	if title is not None:
		ax.set_title(title,loc='center',**ka_title)
	elif len(ka_title)>0:
		raise ValueError('Found keyword arguments for title when no title is specified: '+','.join(ka_title))
	return fig,ax

def colorlegend(ax,
		loc,
		labels,
		colors,
		horizontalalignment='left',
		verticalalignment='top',
		spacex=0,
		spacey=0.15,
		space=None,
		ncol=1,
		**ka):
	"""Draw legends using colored texts
	ax:		axis
	loc:	Location of first line
	labels: List of labels
	colors: List of colors
	spacex,
	spacey:	Extra proportion of space in columns/rows
	space:	Alterantive parameter name for both spacex and spacey
	ncol:	Number of columns
	horizontalalignment,
	verticalalignment,
	ka: Arguments passed to ax.text"""
	import numpy as np
	from matplotlib import transforms
	n=len(labels)
	assert len(colors) == n
	if space is not None:
		spacex=space
		spacey=space
	cid=np.array_split(np.arange(n),ncol)
	t0 = ax.transAxes
	canvas = ax.figure.canvas
	for xi in range(ncol):
		t=t0
		widthmax=0
		for xj in range(len(cid[xi])):
			text = ax.text(*loc,
				labels[cid[xi][xj]],
				color=colors[cid[xi][xj]],
				transform=t,
				horizontalalignment=horizontalalignment,
				verticalalignment=verticalalignment,
				**ka)
			# Need to draw to update the text position.
			text.draw(canvas.get_renderer())
			ex = text.get_window_extent()
			widthmax=max(widthmax,ex.width)
			if xj==0:
				text0_transform=text.get_transform()
			t = transforms.offset_copy(text.get_transform(),
					y=-ex.height * (1 + spacey),
					units='dots')
		t0 = transforms.offset_copy(text0_transform,
				x=widthmax+spacex,
				units='dots')

def heatmap(d,
		optimal_ordering=True,
		method='ward',
		metric='euclidean',
		dshow=None,
		fig=None,
		cmap='coolwarm',
		aspect=0.1,
		figscale=0.02,
		dtop=0.3,
		dright=0,
		wcolorbar=0.03,
		wedge=0.03,
		xselect=None,
		yselect=None,
		xtick=False,
		ytick=True,
		vmin=None,
		vmax=None,
		inverty=True):
	"""
	Draw 2-D hierachical clustering of pandas.DataFrame, with optional hierachical clustering on both axes.
	X/Y axis of figure corresponds to columns/rows of the dataframe.
	d:		Pandas.DataFrame 2D data with index & column names for clustering.
	optimal_ordering: passed to scipy.cluster.hierarchy.dendrogram
	method: Method of hierarchical clustering, passed to scipy.cluster.hierarchy.linkage.
			Accepts single strs or a tuple of two strings for different method options for x and y.
	metric:	Metric to compute pariwise distance for clustering, passed to scipy.spatial.distance.pdist.
			Accepts single strs or a tuple of two strings for different metric options for x and y.
	dshow:	Pandas.DataFrame 2D data with index & column names to draw heatmap. Defaults to d.
	fig:	Figure to plot on.
	cmap:	Colormap
	aspect:	Aspect ratio
	figscale:	Scale of figure compared to font.
	dtop,
	dright:	Top and right dendrogram size. Value from 0 to 1 values as proportion.
			If 0, do not cluster on given axis.
	wcolorbar: Width of colorbar. Value from 0 to 1 values as proportion.
	wedge:	Width of edges and between colorbar and main figure.
			Value from 0 to 1 values as proportion.
	xselect,
	yselect:np.array(bool) of coordinates to draw. Current only selected coordinates are used for clustering.
			TODO: All coordinates in d are used for clustering.
	xtick,
	ytick:	Whether to show ticks.
	vmin,
	vmax:	Minimum/maximum values of heatmap.
	inverty:Whether to invert direction of y.
	Return:
	figure:	plt.Figure drawn on
	x:		column IDs included
	y:		index IDs included.
	"""
	import matplotlib.pyplot as plt
	from scipy.cluster.hierarchy import dendrogram, linkage
	import numpy as np
	assert isinstance(xtick,bool) or (isinstance(xtick,list) and len(xtick) == d.shape[1])
	assert isinstance(ytick,bool) or (isinstance(ytick,list) and len(ytick) == d.shape[0])
	if isinstance(method,str):
		method=[method,method]
	if len(method)!=2:
		raise ValueError('Parameter "method" must have size 2 for x and y respectively.')
	if isinstance(metric,str):
		metric=[metric,metric]
	if metric is not None and len(metric)!=2:
		raise ValueError('Parameter "metric" must have size 2 for x and y respectively.')
	if metric is None:
		assert d.ndim==2 and d.shape[0]==d.shape[1]
		assert (d.index==d.columns).all()
		assert method[0]==method[1]
		if xselect is not None:
			assert yselect is not None
			assert (xselect==yselect).all()
		else:
			assert yselect is None
	if dshow is None:
		dshow=d
	assert dshow.shape==d.shape and (dshow.index==d.index).all() and (dshow.columns==d.columns).all()
	xt0 = d.columns if isinstance(xtick,bool) else xtick
	yt0 = d.index if isinstance(ytick,bool) else ytick
	# Genes to highlight
	d2 = d.copy()
	if xselect is not None:
		d2 = d2.loc[:, xselect]
		dshow=dshow.loc[:,xselect]
		xt0 = [xt0[x] for x in np.nonzero(xselect)[0]]
	if yselect is not None:
		d2 = d2.loc[yselect]
		dshow=dshow.loc[yselect]
		yt0 = [yt0[x] for x in np.nonzero(yselect)[0]]

	wtop = dtop / (1 + d2.shape[0] / 8)
	wright = dright / (1 + d2.shape[1] * aspect / 8)
	iscolorbar = wcolorbar > 0
	t1 = np.array(d2.T.shape)
	t1 = t1 * figscale
	t1[1] /= aspect
	t1[1] /= 1 - wedge * 2 - wtop
	t1[0] /= 1 - wedge * (2 + iscolorbar) - wright - wcolorbar
	if fig is None:
		fig = plt.figure(figsize=t1)

	d3 = dshow.copy()
	if metric is not None:
		# Right dendrogram
		if dright > 0:
			ax1 = fig.add_axes([
				1 - wedge * (1 + iscolorbar) - wright - wcolorbar, wedge, wright,
				1 - 2 * wedge - wtop])
			tl1 = linkage(d2, method=method[1], metric=metric[1], optimal_ordering=optimal_ordering)
			td1 = dendrogram(tl1, orientation='right')
			ax1.set_xticks([])
			ax1.set_yticks([])
			d3 = d3.iloc[td1['leaves'], :]
			yt0 = [yt0[x] for x in td1['leaves']]
		else:
			ax1=None

		# Top dendrogram
		if dtop > 0:
			ax2 = fig.add_axes([wedge, 1 - wedge - wtop, 1 - wedge * (2 + iscolorbar) - wright - wcolorbar, wtop])
			tl2 = linkage(d2.T, method=method[0], metric=metric[0], optimal_ordering=optimal_ordering)
			td2 = dendrogram(tl2)
			ax2.set_xticks([])
			ax2.set_yticks([])
			d3 = d3.iloc[:, td2['leaves']]
			xt0 = [xt0[x] for x in td2['leaves']]
		else:
			ax2=None
	else:
		if dright > 0 or dtop > 0:
			from scipy.spatial.distance import squareform
			tl1 = linkage(squareform(d2), method=method[0], optimal_ordering=optimal_ordering)
			# Right dendrogram
			if dright > 0:
				ax1 = fig.add_axes([
					1 - wedge * (1 + iscolorbar) - wright - wcolorbar, wedge, wright,
					1 - 2 * wedge - wtop])
				td1 = dendrogram(tl1, orientation='right')
				ax1.set_xticks([])
				ax1.set_yticks([])
			else:
				ax1=None
				td1=None
			# Top dendrogram
			if dtop > 0:
				ax2 = fig.add_axes([wedge, 1 - wedge - wtop, 1 - wedge * (2 + iscolorbar) - wright - wcolorbar, wtop])
				td2 = dendrogram(tl1)
				ax2.set_xticks([])
				ax2.set_yticks([])
			else:
				ax2=None
				td2=None
			td0=td1['leaves'] if td1 is not None else td2['leaves']
			d3 = d3.iloc[td0,:].iloc[:,td0]
			xt0,yt0 = [[y[x] for x in td0] for y in [xt0,yt0]]
	axmatrix = fig.add_axes([
		wedge, wedge, 1 - wedge * (2 + iscolorbar) - wright - wcolorbar,
		1 - 2 * wedge - wtop])
	ka = {'aspect': 1 / aspect, 'origin': 'lower', 'cmap': cmap}
	if vmin is not None:
		ka['vmin'] = vmin
	if vmax is not None:
		ka['vmax'] = vmax
	im = axmatrix.matshow(d3, **ka)
	if not isinstance(xtick,bool) or xtick:
		t1 = list(zip(range(d3.shape[1]), xt0))
		t1 = list(zip(*list(filter(lambda x: x[1] is not None, t1))))
		axmatrix.set_xticks(t1[0])
		axmatrix.set_xticklabels(t1[1], minor=False, rotation=90)
	else:
		axmatrix.set_xticks([])
	if not isinstance(ytick,bool)or ytick:
		t1 = list(zip(range(d3.shape[0]), yt0))
		t1 = list(zip(*list(filter(lambda x: x[1] is not None, t1))))
		axmatrix.set_yticks(t1[0])
		axmatrix.set_yticklabels(t1[1], minor=False)
	else:
		axmatrix.set_yticks([])

	axmatrix.tick_params(top=False,
		bottom=True,
		labeltop=False,
		labelbottom=True,
		left=True,
		labelleft=True,
		right=False,
		labelright=False)
	if inverty:
		if ax1 is not None:
			ax1.set_ylim(ax1.get_ylim()[::-1])
		axmatrix.set_ylim(axmatrix.get_ylim()[::-1])
	if wcolorbar > 0:
		cax = fig.add_axes([
			1 - wedge - wcolorbar, wedge, wcolorbar, 1 - 2 * wedge - wtop])
		fig.colorbar(im, cax=cax)

	return fig, d3.columns, d3.index

def dotplot(ds,dc,fig=None,figsize=0.2,size_transform=lambda x:x,sizes=None,vrange=None,cmap='viridis',**ka):
	"""
	Draw dotplot
	ds:				Dot size parameter as pandas.DataFrame
	dc:				Dot color parameter as pandas.DataFrame
	fig:			fig to darw on
	figsize:		Figure scale if fig is None
	size_transform:	function to transform values in ds to actual size
	sizes:			List of dot size parameters to show as legend
	vrange:			Color [min,max] range for normalizing dc
	cmap:			Colormap
	ka:				Keyword arguments passed to plt.scatter
	Return:
	fig:			plt.Figure drawn on
	ax:				matplotlib.axes._subplots.AxesSubplot drawn on
	"""
	import numpy as np
	import matplotlib.pyplot as plt
	assert ds.shape==dc.shape
	assert (ds.index==dc.index).all()
	assert (ds.columns==dc.columns).all()
	ns=[len(ds.columns),len(dc.index)]
	colorbar_width=0.5
	
	#Set color range
	if vrange is None:
		vrange=[None,None]
	assert len(vrange)==2
	if vrange[0] is None:
		vrange[0]=dc.values.min()
	if vrange[1] is None:
		vrange[1]=dc.values.max()
	#Create axis
	if fig is None:
		assert figsize is not None
		if not hasattr(figsize,'__len__'):
			figsize=[figsize*(ns[0]*2+colorbar_width),figsize*ns[1]]
		fig=plt.figure(figsize=figsize)
	ax=fig.add_gridspec(1,3,width_ratios=[ns[0],ns[0],colorbar_width])
	#Create coordinates
	dx=[np.repeat(np.arange(ns[x]).reshape(1,-1).swapaxes(x,0),ns[1-x],axis=x) for x in range(2)]
	ax1=fig.add_subplot(ax[0,0],aspect=1)
	ax1.scatter(dx[0].ravel(),dx[1].ravel(),marker='o',lw=0,s=size_transform(ds.values.ravel()),c=dc.values.ravel(),vmin=dc.values.min(),vmax=dc.values.max(),cmap=cmap,**ka)
	ax1.set_xlim([-0.5,ns[0]-0.5])
	ax1.set_ylim([-0.5,ns[1]-0.5])
	ax1.set_xticks(np.arange(ns[0]))
	ax1.set_xticklabels(ds.columns,rotation=90)
	ax1.set_yticks(np.arange(ns[1]))
	ax1.set_yticklabels(ds.index)
	ax2=fig.add_subplot(ax[0,1],sharey=ax1,sharex=ax1,aspect=1)
	ax2.set_axis_off()
	#Create coordinates
	if sizes is not None:
		n=len(sizes)
		dx2=[np.zeros(n),np.arange(n)*1.3]
		ax2.scatter(dx2[0],dx2[1],marker='o',lw=0,s=size_transform(sizes),c='k',**ka)
		for xi in range(n):
			ax2.text(dx2[0][xi]+0.7,dx2[1][xi],str(sizes[xi]),ha='left',va='center')
	ax3=fig.add_subplot(ax[0,2],aspect=15)
	matplotlib.colorbar.ColorbarBase(ax3,cmap=plt.get_cmap(cmap),norm=matplotlib.colors.Normalize(vmin=vrange[0],vmax=vrange[1]),orientation='vertical')
	return (fig,ax1)







































assert __name__ != "__main__"
