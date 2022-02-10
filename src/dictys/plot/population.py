#!/usr/bin/python3
# Lingfei Wang, 2020-2022. All rights reserved.

"""
Population level visualization of networks 
"""

def compute_reg_spec(d0,min_entropy=0.5,ncut=0.4,nmin=20,nmax_reg=10,select_state=None):
	"""
	Compute state specificity of regulators with out-degree centrality and CPM.
	Specificity is defined as value in this state / sum across all states.

	Parameters
	----------
	d0:				dictys.net.network
		Input network class
	min_entropy:	float
		State specificity entropy level required (relative to random assignment) to select regulator. Lower means more specific.
	ncut:			float
		Minimum probability required
	nmin:			int
		Minimum number of targets required
	nmax_reg:		int
		Maximum TF count for each state, selected based on probability
	select_state:	list of str
		Only compute specificity among given states

	Returns
	-------
	n:		numpy.ndarray(shape=[n_reg,n_state])
		Number of targets
	v:		numpy.ndarray(shape=[n_reg,n_state])
		Degree centrality specificity (value in this state/sum across all states)
	cpm:	numpy.ndarray(shape=[n_reg,n_state])
		CPM
	cpm_v:	numpy.ndarray(shape=[n_reg,n_state])
		CPM specificity (value in this state/sum across all states)
	reg:	numpy.ndarray
		Regulators selected from given constraints
	reg_s:	numpy.ndarray
		Each regulator's specific state corresponding to reg
	"""
	from dictys.net import stat
	import numpy as np
	import pandas as pd

	if select_state is None:
		select_state=np.arange(d0.sn)
	else:
		select_state=np.array([d0.sdict[x] for x in select_state])
	#Network mask
	mask=stat.net(d0,varname='mask')
	#Binary network
	binnet=stat.fbinarize(stat.net(d0),statmask=mask)
	#Target count
	na=stat.fcentrality_degree(binnet,roleaxis=0).compute(select_state)
	#Degree centrality rate
	dcrate=stat.fcentrality_degree(binnet,statmask=mask,roleaxis=0).compute(select_state)
	#Degree centrality specificity
	va=(dcrate.T/(dcrate.sum(axis=1)+1E-300)).T
	#CPM
	cpm=2**(stat.lcpm(d0,cut=-1,const=1).compute(select_state)[d0.nids[0]])-1

	assert (na>=0).all() and (va>=0).all() and (va<=1).all()
	assert na.shape==va.shape
	s=np.arange(va.shape[0])

	#Remove regulators not having any target
	t1=(va>0).any(axis=1)
	s=s[t1]
	v=va[t1]
	n=na[t1]
	
	#Filter regulators by entropy
	t1=(v*np.log(v+1E-300)).sum(axis=1)
	t1=t1>=np.log(1/v.shape[1])*min_entropy
	s=s[t1]
	v=v[t1]
	n=n[t1]
	
	#Filter regulators
	t1=(v>=ncut)&(n>=nmin)
	#Limit max regulator count
	t1=[np.nonzero(x)[0] for x in t1.T]
	if nmax_reg is not None:
		t1=[t1[x][v[t1[x],x]>=np.partition(v[t1[x],x],-nmax_reg)[-nmax_reg]] if len(t1[x])>nmax_reg else t1[x] for x in range(len(t1))]
	t1=np.unique(np.concatenate(t1))
	
	#Order regulators: First state, then strength within state
	t1=[t1,v[t1].argmax(axis=1)]
	t1.append(v[t1[0],t1[1]])
	t1=[list(x) for x in zip(*t1)]
	t1=[np.array(list(x)) for x in zip(*sorted(t1,key=lambda x:[x[1],-x[2]]))]
	# t1=[x[0] for x in t1]
	s=[d0.nname[d0.nids[0]][s[t1[0]]],d0.sname[select_state[t1[1]]]]
	
	#CPM specificity
	cpm_v=(cpm.T/(cpm.sum(axis=1)+1E-300)).T
	
	assert all([x.shape==na.shape for x in [va,cpm,cpm_v]])
	assert (na>=0).all()
	assert (va>=0).all() and (va<=1).all()
	assert (cpm>=0).all()
	assert (cpm_v>=0).all() and (cpm_v<=1).all()
	na,va,cpm,cpm_v=[pd.DataFrame(x,index=d0.nname[d0.nids[0]],columns=d0.sname[select_state]) for x in [na,va,cpm,cpm_v]]
	return (na,va,cpm,cpm_v,s[0],s[1])

def fig_heatmap_reg_spec(v,aspect=0.3,figscale=0.15,g_ann=None,**ka):
	"""
	Draw heatmap for regulators' state specificity.

	Parameters
	----------
	v:			pandas.DataFrame
		State specificity with states as rows and regulator genes as columns
	aspect:		float
		Aspect ratio of each entry
	figscale:	float
		Size scale of figure
	g_ann:		list of str or None
		List of regulator genes to annotate
	ka:			dict
		Keyword arguments passed to dictys.plot.heatmap

	Returns
	-------
	matplotlib.pyplot.Figure
		Figure drawn
	"""
	import numpy as np
	import matplotlib.pyplot as plt
	from dictys.plot import heatmap
	
	vrange=[v.values.min(),v.values.max()]
	#Draw heatmap
	figscale=figscale*aspect
	ka1=dict(metric=[lambda u,v:1-np.sqrt(u*v).sum()]*2,method='average',vmin=vrange[0],vmax=vrange[1],aspect=aspect,dtop=0.25,dright=0,figscale=figscale,xtick=False,ytick=True,cmap='viridis',optimal_ordering=False,colorbar=0.5/aspect/v.shape[0],wedge=0.75/aspect/v.shape[0])
	ka1.update(ka)
	g=heatmap(v.T,**ka1)
	fig1=plt.gcf()
	ax=fig1.axes[-2]
	#Annotate select regulators
	gdict=dict(zip(g[1],range(len(g[1]))))
	t1=list(filter(lambda x:x in gdict,g_ann))
	t2=[gdict[x] for x in t1]
	ax.set_xticks(t2)
	ax.set_xticklabels(t1,rotation=90)
	ax.tick_params(which='both',axis='both',top=True,labeltop=False)
	return fig1

def fig_heatmap_top(d0,selection,ntop=10,direction=0,gann=[],cmap_value='coolwarm',cmap_type='tab10',color_type=None,normalization='column',aspect=0.2,topheight=0.7,topspace=0.3,figscale=0.15):
	"""
	Draw heatmap for top targets of given regulators in given cell types/states.

	Parameters
	----------
	d0:				dictys.net.network
		Input network
	selection:		list of tuple
		List of regulator-state pair (or TF-cell type pair) to draw. Each element of list is a tuple (TF, cell type) by name.
	ntop:			int
		Number of top targets to draw
	direction:		int
		Direction of regulation to be considered for top targets. For -1,0,1, indicates repression,both,activation only.
	gann:			list of str or str
		Target genes to annotate. Accepts 'all' for all genes.
	cmap_value:		str
		Matplotlib colormap name for heatmap
	cmap_type:		str
		Matplotlib colormap name for cell state/type bar. Ignored if color_type is specified.
	color_type:		dict
		Dictionary mappping cell state name to bar color.
	normalization:	str
		Normalization of regulation strength for heatmap. Accepts:
		column:	All regulations for each heatmap column are scale normalized to -1 to 1 so the strongest is -1 or 1.
		none:	No normalization
	aspect:			float
		Aspect ratio of heatmap
	topheight:		float
		Relative height of top colorbar for cell types/states
	topspace:		float
		Relative spacing between colorbar and heatmap
	figscale:		float
		Figure size scale

	Returns
	-------
	fig:		matplotlib.pyplot.Figure
		Heatmap drawn
	colorbar:	matplotlib.pyplot.Figure
		Colorbar drawn for heatmap
	data:		pandas.DataFrame
		Data of heatmap
	"""
	from os import linesep
	import logging
	import numpy as np
	import pandas as pd
	import matplotlib
	import matplotlib.pyplot as plt

	ndict=[dict(zip(d0.nname[x],range(len(x)))) for x in d0.nids]
	names=np.array([' '+'-'.join(x) for x in selection])
	if gann=='all':
		gann=d0.nname[d0.nids[1]]
	gann=set(gann)

	#Validations
	t1=np.nonzero([x[0] not in ndict[0] for x in selection])[0]
	if len(t1)>0:
		raise ValueError('TF(s) not found in dataset: '+','.join([selection[x][0] for x in t1]))
	t1=np.nonzero([x[1] not in d0.sdict for x in selection])[0]
	if len(t1)>0:
		raise ValueError('Celltype(s) not found in dataset: '+','.join([selection[x][1] for x in t1]))
	if direction not in {-1,0,1}:
		raise ValueError('Parameter direction must be -1, 0, or 1.')

	#Cell type colors
	color_type=None
	if color_type is None:
		if cmap_type is None:
			raise ValueError('At least one of color_type and cmap_type must be assigned.')
		typelist=[]
		for xi in selection:
			if xi[1] not in typelist:
				typelist.append(xi[1])
		t1=plt.get_cmap(cmap_type)
		if hasattr(t1,'N') and t1.N<len(typelist):
			raise RuntimeError('Colormap has fewer colors than cell types.')
		t1=t1(np.linspace(0,1,t1.N)[:len(typelist)])
		color_type=dict(zip(typelist,t1))

	#Obtain network
	net=[(ndict[0][x[0]],d0.sdict[x[1]]) for x in selection]
	net=np.array([d0.prop['es']['w'][x[0],:,x[1]] for x in net])
	#Keep only requested directions for target gene selection
	net2=net.copy()
	if direction==-1:
		net2*=net2<0
	elif direction==1:
		net2*=net2>0
	net2=np.abs(net2)
	#Validations
	t1=(net2!=0).sum(axis=1)
	t2=np.nonzero(t1==0)[0]
	if len(t2)>0:
		raise ValueError('TF-cell type pair(s) have no targets in the specified direction: '+','.join(names[t1]))
	t2=np.nonzero(t1<ntop)[0]
	if len(t2)>0:
		logging.warn('TF-cell type pair(s) have <{} targets in the specified direction: '.format(ntop)+','.join(names[t1]))
	#Select target genes to show
	t1=np.argpartition(net2,-ntop,axis=1)[:,-ntop:]
	t1=np.array([t1[x,(net2 if direction!=0 else net)[x,t1[x]].argsort()[::-1]] for x in range(len(t1))])
	t2=[]
	for xi in t1.ravel():
		if xi not in t2:
			t2.append(xi)
	t1=np.array(t2)
	#Filter network and gene names
	net=net[:,t1]
	net2=net2[:,t1]
	names2=d0.nname[d0.nids[1][t1]]
	ns=np.array(net.shape)
	#Normalize network
	if normalization=='column':
		net=(net.T/net2.max(axis=1)).T
	elif normalization is None or normalization=='none':
		pass
	else:
		raise ValueError('Unknown value for parameter normalization')
	vmax=(net*direction).max() if direction!=0 else np.abs(net).max()
	
	#Figure size
	figsize1=(figscale*ns[0],figscale*aspect*ns[1])
	figsize2=(figsize1[0],figscale*(topheight+topspace))
	figsize=(figsize1[0],figsize1[1]+figsize2[1])
	fig=plt.figure(figsize=figsize)
	#Panel for cell type color
	ax=fig.add_axes([0,1-topheight*figscale/figsize[1],1, topheight*figscale/figsize[1]])
	t1=np.array([color_type[x[1]] for x in selection])
	t1=t1.reshape(1,*t1.shape)
	ax.imshow(t1,aspect=topheight);
	ax.set_xticks(np.arange(ns[0]));
	ax.set_xticklabels(names,rotation=90);
	ax.set_yticks([])
	ax.tick_params(bottom=False,top=False,left=False,labeltop=True,labelbottom=False,length=0);
	[x.set_visible(False) for x in ax.spines.values()]
	#Panel for heatmap
	ax=fig.add_axes([0,0,1, figsize1[1]/figsize[1]],sharex=ax)
	ax.imshow(net.T,cmap=cmap_value,vmin=-vmax,vmax=vmax,aspect=aspect);
	gann2=np.nonzero([x in gann for x in names2])[0]
	ax.set_yticks(gann2)
	ax.set_yticklabels(names2[gann2])
	ax.tick_params(bottom=False,top=False,left=True,labeltop=False,labelbottom=False);
	
	#Figure for colorbar
	fig2,ax=plt.subplots(figsize=(0.15, 0.8))
	norm=matplotlib.colors.Normalize(vmin=-1, vmax=1)
	fig2.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap_value),cax=ax,orientation='vertical')
	ax.set_title(f'Relative{linesep}strength',loc='center',pad=10)
	
	#Prepare data
	net=pd.DataFrame(net.T,index=names2,columns=names)
	
	return (fig,fig2,net)





































































assert __name__ != "__main__"
