#!/usr/bin/python3
# Lingfei Wang, 2020-2022. All rights reserved.

"""
Population level visualization of networks 
"""

def compute_reg_spec(d0,min_entropy=0.5,ncut=0.4,nmin=20,nmax_reg=10,select_state=None):
	"""
	Compute state specificity of regulators with out-degree centrality and CPM.
	d0:  Dataset
	min_entropy: State specificity entropy level required (relative to random assignment) to select regulator. Lower means more specific.
	ncut: Minimum probability required
	nmin: Minimum number of targets required
	nmax_reg: Maximum TF count for each cell type, selected based on probability
	Return:
	n:		Number of targets as np.array(shape=[n_reg,n_state])
	v:		Degree centrality specificity as np.array(shape=[n_reg,n_state])
	cpm:	CPM as np.array(shape=[n_reg,n_state])
	cpm_v:	CPM specificity as np.array(shape=[n_reg,n_state])
	reg:	Regulators selected from given constraints as np.array
	reg_s:	Each regulator's specific state corresponding to reg as np.array
	"""
	from dictys.net import stat
	import numpy as np
	import pandas as pd

	if select_state is None:
		select_state=np.arange(d0.sn)
	else:
		select_state=np.array([d0.sdict[x] for x in select_state])
	#Network mask
	mask=stat.stat_net(d0,varname='mask')
	#Binary network
	binnet=stat.statf_binarize(stat.stat_net(d0),mask)
	#Target count
	na=stat.statf_centrality_degree(binnet,roleaxis=0).compute(select_state)
	#Degree centrality rate
	dcrate=stat.statf_centrality_degree(binnet,statmask=mask,roleaxis=0).compute(select_state)
	#Degree centrality specificity
	va=(dcrate.T/(dcrate.sum(axis=1)+1E-300)).T
	#CPM
	cpm=2**(stat.stat_lcpm(d0,cut=-1,const=1).compute(select_state)[d0.nids[0]])-1

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
	v:		pandas.DataFrame of state specificity with states as rows and regulator genes as columns
	aspect:	Aspect ratio of each entry
	g_ann:	List of regulator genes to annotate
	ka:		Keyword arguments passed to dictys.plot.heatmap
	Return:	plt.Figure drawn on
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








































































assert __name__ != "__main__"
