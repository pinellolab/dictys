#!/usr/bin/python3
# Lingfei Wang, 2022, 2023. All rights reserved.

"""
Gene regulatory network reconstruction with TF binding network and transcriptome data
"""

import abc
from typing import Optional,Tuple,Callable,Union
from dictys.utils.file import csv_compression

_docstring2argparse_ignore_=['Trace_ELBO_site','SVI_multiloss','model_base','model_covariance','model_ou']

_model_names=frozenset(_docstring2argparse_ignore_)

def __getattr__(name):
	# Pyro classes live in dictys.network_model, imported on first use
	if name in _model_names:
		from dictys import network_model
		return getattr(network_model,name)
	raise AttributeError(f'module {__name__!r} has no attribute {name!r}')

###########################################################################
# Network reconstruction
###########################################################################

def reconstruct(fi_exp:str,fi_mask:str,fo_weight:str,fo_meanvar:str,fo_covfactor:str,fo_loss:str,fo_stats:str,lr:float=0.01,lrd:float=0.999,nstep:float=4000,npc:int=0,fi_cov:Optional[str]=None,model:str='ou',nstep_report:int=100,rseed:int=12345,device:str='cpu',dtype:str='float',loss:str='Trace_ELBO_site',nth:float=1,varmean:str='N_0val',varstd:Optional[str]=None,fo_weightz:Optional[str]=None,scale_lyapunov:float=1E5)->None:
	"""
	Reconstruct network with any pyro model in net_pyro_models that is based on covariance_model and has binary masks.

	Parameters
	----------
	fi_exp:
		Path of input tsv file of expression matrix.
	fi_mask:
		Path of input tsv file of mask matrix indicating which edges are allowed. Can be output of dictys chromatin binlinking.
	fo_weight:
		Path of output tsv file of edge weight matrix
	fo_meanvar:
		Path of output tsv file of mean and variance of each gene's relative log expression
	fo_covfactor:
		Path of output tsv file of factors for the off-diagonal component of gene covariance matrix
	fo_loss:
		Path of output tsv file of sublosses in each training step
	fo_stats:
		Path of output tsv file of basic stats of each variable during training
	lr:
		Initial iearning rate
	lrd:
		Learning rate decay
	nstep:
		Number of training steps
	npc:
		Number of unknown factors for covariance in multivariate distribution
	fi_cov:
		Path of input tsv file of covariate matrix for each cell to be included. Should have cell x covariate shape.
	model:
		Name of model to train
	nstep_report:
		Number of steps to save each parameter distribution
	rseed:
		Initial random seed
	device:
		Device for pytorch and pyro. See TBA....
	dtype:
		Data type for pytorch. See TBA...
	loss:
		Loss function name
	nth:
		Number of threads for CPU usage. When <1, use nth*(detected core count).
	varmean:
		Pyro parameter name for mean of effect size
	varstd:
		Pyro parameter name for std of effect size to compute z scores. Use None if unavailable.
	fo_weightz:
		Path of output tsv file of z score matrix of edge weights
	scale_lyapunov:
		Scale of Lyapunov equation inequality loss
	"""
	import itertools
	import logging
	import numpy as np
	import pandas as pd
	from dictys.utils.parallel import num_threads,autocount
	from dictys.utils.importing import torch,pyro
	from dictys.network_model import Trace_ELBO_site,model_ou
	#Initialize
	nth=autocount(nth)
	pyro.set_rng_seed(rseed)
	np.random.seed(rseed)
	if model=='ou':
		model=model_ou
	else:
		raise ValueError(f'Unknown network model {model}.')
	if fo_weightz is not None and varstd is None:
		raise ValueError('Output file for weight z score only available if varstd is specified.')
	dtype=getattr(torch,dtype)

	#Loading data
	logging.info(f'Reading file {fi_exp}')
	dt0=pd.read_csv(fi_exp,header=0,index_col=0,sep='\t')
	logging.info(f'Reading file {fi_mask}')
	mask=pd.read_csv(fi_mask,header=0,index_col=0,sep='\t')
	assert len(set(mask.index)-set(dt0.index))==0
	assert len(set(mask.columns)-set(dt0.index))==0
	if mask.shape[0]==0:
		raise RuntimeError('No regulator found.')
	nt,ns=dt0.shape

	#Covariates: load from file
	if fi_cov is not None:
		logging.info(f'Reading file {fi_cov}.')
		dc=pd.read_csv(fi_cov,header=0,index_col=0,sep='\t')
		if dc.shape[0]==0:
			raise ValueError(f'No covariate found in {fi_cov}. Make sure it is in cell x covariate format.')
		if len(set(dt0.columns)-set(dc.index))>0:
			raise ValueError(f'Found cells not contained in {fi_cov}. Make sure it is in cell x covariate format.')
		dc=dc.loc[dt0.columns].values
	else:
		dc=np.array([],dtype=float).reshape(dt0.shape[1],0)
	t1=min(nt,ns)-dc.shape[1]
	if npc>=t1:
		raise RuntimeError(f'Insufficient degrees of freedom {t1} compared to off-diagonal factor count {npc} in unexplained variance matrix.')

	#Match shape and order between read count matrix and mask
	namereg=mask.index
	nametarget=mask.columns
	nreg,ntarget=[len(x) for x in [namereg,nametarget]]
	t1=set(namereg)
	namet=np.r_[namereg,list(filter(lambda x:x not in t1,nametarget))]
	assert len(namet)==ntarget
	dt0=dt0.loc[namet].values
	mask=mask[namet].values.astype(bool)

	#Prepare other parameters
	if loss=='Trace_ELBO_site':
		loss=Trace_ELBO_site(model.varnames).differentiable_loss
	#Run pyro
	with num_threads(nth):
		model=model('pyromodel',np.array(np.nonzero(mask)),mask.shape,npc,dc,floattype=dtype,device=device,scale_lyapunov=scale_lyapunov)
		model.gen_params({'G_obs':model.tensor(dt0.astype(int),dtype=torch.int)})
		model.init('ClippedAdam',loss,optimizer_a=[[{'lr': lr, "betas": (0.90, 0.999),'lrd': lrd}],{}])
		model.train_svi(nstep,nstep_report=nstep_report)

	#Recover networks: mean and z score of each edge
	ans_mean=torch.zeros([nreg,mask.shape[1]],**model.tensorka)
	ans_mean[model.N_0mask[0],model.N_0mask[1]]=model.param(varmean)
	ans_mean=ans_mean.detach().cpu().numpy().astype(float)
	assert ans_mean.shape==(nreg,ntarget) and np.isfinite(ans_mean).all()
	assert (ans_mean[~mask]==0).all()
	if varstd is not None:
		ans_z=torch.zeros([nreg,mask.shape[1]],**model.tensorka)
		ans_z[model.N_0mask[0],model.N_0mask[1]]=model.param(varmean)/model.param(varstd)
		ans_z=ans_z.detach().cpu().numpy().astype(float).reshape(nreg,ntarget)
		assert ans_z.shape==(nreg,ntarget) and np.isfinite(ans_z).all()
		assert (ans_z[~mask]==0).all()
	#Recover other model parameters
	ans_gmean=model.param('mu_G').detach().cpu().numpy().astype(float)
	ans_gcov0_d=model.param('sigma0_G_d').detach().cpu().numpy().astype(float)
	ans_gcov0_nd=model.param('sigma0_G_nd').detach().cpu().numpy().astype(float).T
	ans_gcov_d=model.param('sigma_G_d').detach().cpu().numpy().astype(float)
	ans_gcov_nd=model.param('sigma_G_nd').detach().cpu().numpy().astype(float).T
	ans_reg=model.param('reg').detach().cpu().numpy().astype(float).ravel()[0]
	#Extract stats
	ans_loss=np.array(model.stat_loss,dtype=float)
	ans_lossname=np.array(model.varnames)
	ans_statsname=[]
	ans_statstype=list(model.stat_param_names)
	ans_stats=[]
	ans_steps=None
	for xi in model.p:
		if np.sum([len(x) for x in model.stat_param[xi]])==0:
			continue
		ans_stats.append([])
		ans_statsname.append(xi)
		for xj in range(len(model.stat_param[xi])):
			t1=np.array(list(model.stat_param[xi][xj].items())).T
			if ans_steps is None:
				ans_steps=t1[0]
			else:
				assert len(t1[0])==len(ans_steps) and np.all(t1[0]==ans_steps)
			ans_stats[-1].append(t1[1])
	ans_stats,ans_statsname,ans_statstype,ans_steps=[np.array(x) for x in [ans_stats,ans_statsname,ans_statstype,ans_steps]]
	ans_stats=ans_stats.astype(float)
	ans_steps=ans_steps.astype('u4')

	#Normalize to unit regularization
	ans_mean=ans_mean/ans_reg
	ans_gcov0_d=ans_gcov0_d/ans_reg
	ans_gcov0_nd=ans_gcov0_nd/np.sqrt(ans_reg)
	#Normalize to unit total expression
	ans_gmean=ans_gmean-np.log(np.exp(ans_gmean).sum())
	assert all(len(x)>0 for x in [ans_statsname,ans_statstype,ans_steps])
	assert ans_stats.shape==(len(ans_statsname),len(ans_statstype),len(ans_steps))
	assert ans_loss.ndim==2 and ans_loss.shape[1]==len(ans_lossname)
	assert (ans_steps[1:]>ans_steps[:-1]).all() and np.isfinite(ans_stats).all()
	#Reformat output
	ans_weight=pd.DataFrame(ans_mean,index=namereg,columns=namet)
	ans_meanvar=np.array([ans_gmean,ans_gcov_d]).T
	ans_meanvar=pd.DataFrame(ans_meanvar,index=namet,columns=['mean','var'])
	ans_nd=pd.DataFrame(ans_gcov_nd.T,index=namet,columns=['factor'+str(x+1) for x in range(ans_gcov_nd.shape[0])])
	ans_loss=pd.DataFrame(ans_loss,index=np.arange(ans_loss.shape[0]),columns=ans_lossname)
	ans_stats=[[ans_statsname[x[0]],ans_statstype[x[1]],ans_steps[x[2]],ans_stats[x[0],x[1],x[2]]] for x in itertools.product(*[range(y) for y in ans_stats.shape])]
	ans_stats=pd.DataFrame(ans_stats,columns=['variable','stat','step','value'])

	#Write output
	logging.info(f'Writing file {fo_weight}')
	ans_weight.to_csv(fo_weight,header=True,index=True,sep='\t',compression=csv_compression(fo_weight))
	logging.info(f'Writing file {fo_meanvar}')
	ans_meanvar.to_csv(fo_meanvar,header=True,index=True,sep='\t',compression=csv_compression(fo_meanvar))
	logging.info(f'Writing file {fo_covfactor}')
	ans_nd.to_csv(fo_covfactor,header=True,index=True,sep='\t',compression=csv_compression(fo_covfactor))
	logging.info(f'Writing file {fo_loss}')
	ans_loss.to_csv(fo_loss,header=True,index=True,sep='\t',compression=csv_compression(fo_loss))
	logging.info(f'Writing file {fo_stats}')
	ans_stats.to_csv(fo_stats,header=True,index=True,sep='\t',compression=csv_compression(fo_stats))
	if varstd is not None:
		ans_weightz=pd.DataFrame(ans_z,index=namereg,columns=namet)
		logging.info(f'Writing file {fo_weightz}')
		ans_weightz.to_csv(fo_weightz,header=True,index=True,sep='\t',compression=csv_compression(fo_weightz))

###########################################################################
# Network post processing
###########################################################################

def indirect(fi_weight:str,fi_covfactor:str,fo_iweight:str,norm:int=3,fi_meanvar:Optional[str]=None,eigmin:Optional[float]=0.2,eigmax:Optional[float]=1.8,multiplier:float=1.1,nth:float=1)->None:
	"""
	Computes steady-state indirect effect of gene perturbation from OU process.

	Performs extra regularization on network by bounding the eigenvalues of feedback loops with parameters eigmin and eigmax. Values away from 1 indicates stronger feedback loop effects. Set values closer to 1 to apply stronger regularization.

	Parameters
	----------
	fi_weight:
		Path of input tsv file of edge weight matrix
	fi_covfactor:
		Path of iput tsv file of factors for the off-diagonal component of gene covariance matrix
	fo_iweight:
		Path of output tsv file of steady-state indirect effect edge weight matrix
	norm:
		Normalizations of indirect effect matrix. Binary flag variable that accepts:
		
		1:	Divides the indirect effect size of each TF on other genes with that on itself to obtain logFC per unit logFC of TF.
		
		2:	Accounts for mRNA saturation/sampling effects.

	fi_meanvar:
		Path of iput tsv file of mean and variance of each gene's relative log expression. Needed if `norm` & 2.
	eigmin:
		Lower bound of eigenvalue for finding the minimal required extra regularization.
	eigmax:
		Upper bound of eigenvalue for finding the minimal required extra regularization.
	multiplier:
		Step size of strengthening regularization. Larger values provide faster but over-regularization. Must be greater than 1.
	nth:
		Number of threads for CPU usage. When <1, use nth*(detected core count).
	"""
	import logging
	import numpy as np
	import pandas as pd
	from dictys.utils.parallel import num_threads,autocount

	assert multiplier>1
	if norm not in {0,1,2,3}:
		raise ValueError('Invalid norm value. Choose from: 0, 1, 2, 3.')
	nth=autocount(nth)
	eigtarget=[eigmin if eigmin is not None else -np.inf,eigmax if eigmax is not None else np.inf]
	assert eigtarget[0]<1 and eigtarget[1]>1		# pylint: disable=R1716

	#Loading data
	logging.info(f'Reading file {fi_weight}')
	dw=pd.read_csv(fi_weight,header=0,index_col=0,sep='\t')
	assert dw.index.isin(dw.columns).all()
	t1=set(dw.index)
	t1=list(dw.index)+list(filter(lambda x:x not in t1,dw.columns))
	dw=dw[t1]
	d=dw.values
	ns=d.shape
	net0=np.concatenate([d,np.zeros((ns[1]-ns[0],ns[1]),dtype=float)],axis=0)
	logging.info(f'Reading file {fi_covfactor}')
	df=pd.read_csv(fi_covfactor,header=0,index_col=0,sep='\t')
	assert len(dw.columns)==len(df.index) and set(dw.columns)==set(df.index)
	df=df.loc[dw.columns]
	if norm&2 is not None:
		logging.info(f'Reading file {fi_meanvar}')
		dmv=pd.read_csv(fi_meanvar,header=0,index_col=0,sep='\t')
		assert len(dw.columns)==len(dmv.index) and set(dw.columns)==set(dmv.index)
		dmv=dmv.loc[dw.columns]
	
	with num_threads(nth):
		#Determine parameters to re-regularize OU process
		reg=1
		while True:
			net1=net0/reg
			#Compute steady state network
			beta=np.eye(ns[1])-net1.T
			betai=np.linalg.inv(beta)
			dx=betai.T
			if norm&2 is not None:
				xtot=np.exp(dmv['mean'].values).sum()
				dxtot=np.exp(dmv['mean'].values)@betai.T
				dx=dx-dxtot/xtot
			t1=np.diag(dx)[:ns[0]]
			if t1.min()>=eigtarget[0] and t1.max()<=eigtarget[1]:
				break
			reg*=multiplier
		net1=dx*reg

		if norm & 1:
			#Normalization of effect size of base rate perturbation using effect size on the perturbed TF itself.
			net1=(net1.T/np.diag(net1)).T
		net1[np.diag_indices_from(net1)]=0

	#Convert to output
	net1=pd.DataFrame(net1[:ns[0]],index=dw.index,columns=dw.columns)
	logging.info(f'Writing file {fo_iweight}.')
	net1.to_csv(fo_iweight,header=True,index=True,sep='\t',compression=csv_compression(fo_iweight))

def normalize(fi_weight:str,fi_meanvar:str,fi_covfactor:str,fo_nweight:str,norm:int=3,nth:float=1)->None:
	"""
	Normalize edge strength.

	So they are more resistant to estimation bias of true expression level variance.

	Parameters
	----------
	fi_weight:
		Path of input tsv file of edge weight matrix
	fi_meanvar:
		Path of iput tsv file of mean and variance of each gene's relative log expression
	fi_covfactor:
		Path of iput tsv file of factors for the off-diagonal component of gene covariance matrix
	fo_nweight:
		Path of output tsv file of normalized edge weight matrix
	norm:
		Type of normalization as binary flag values. Accepts:

		* 1:	Multiplying edge weight with stochastic noise std of TF
		
		* 2:	Dividing edge weight with stochastic noise std of target

	nth:
		Number of threads for CPU usage. When <1, use nth*(detected core count).
	"""
	import logging
	import numpy as np
	import pandas as pd
	from dictys.utils.parallel import num_threads,autocount

	if norm not in {1,2,3}:
		raise ValueError('Invalid norm value. Choose from: 1, 2, 3.')
	nth=autocount(nth)
	#Loading data
	logging.info(f'Reading file {fi_weight}')
	dw=pd.read_csv(fi_weight,header=0,index_col=0,sep='\t')
	logging.info(f'Reading file {fi_meanvar}')
	dmv=pd.read_csv(fi_meanvar,header=0,index_col=0,sep='\t')
	logging.info(f'Reading file {fi_covfactor}')
	df=pd.read_csv(fi_covfactor,header=0,index_col=0,sep='\t')
	assert dw.index.isin(dw.columns).all()
	assert len(dw.columns)==len(dmv.index) and set(dw.columns)==set(dmv.index)
	assert len(dw.columns)==len(df.index) and set(dw.columns)==set(df.index)
	t1=set(dw.index)
	t1=list(dw.index)+list(filter(lambda x:x not in t1,dw.columns))
	dw=dw[t1]
	dmv=dmv.loc[t1]
	df=df.loc[t1]
	d=dw.values
	#Normalize
	with num_threads(nth):
		cov=np.sqrt(dmv['var'].values+np.diag(df.values@df.values.T))
		if norm&1:
			#Edge weight*std of TF
			d=(d.T*cov[:dw.shape[0]]).T
		if norm&2:
			#Edge weight/std of target
			d=d/cov
	d=pd.DataFrame(d,index=dw.index,columns=dw.columns)
	#Saving data
	logging.info(f'Writing file {fo_nweight}')
	d.to_csv(fo_nweight,header=True,index=True,sep='\t',compression=csv_compression(fo_nweight))

###########################################################################
# Network storage
###########################################################################

def tofile(diri_data:str,diri_work:str,fi_subsets:str,fo_networks:str,dynamic:bool=False,nettype:str='n',optional:str='readcount',fi_c:Optional[str]=None)->None:
	"""
	Saving networks to a single file.

	Parameters
	-------------
	diri_data:
		Path of input data folder to load from
	diri_work:
		Path of input working folder to load from
	fi_subsets:
		Path of input txt file for cell subset names
	fo_networks:
		Path of output h5 file for all networks
	dynamic:
		Whether to load a dynamic network instead of a set of static networks
	nettype:
		Type of network. Accepts:

		* '':	Unnormalized direct network

		* 'n':	Normalized direct network

		* 'i':	Unnormalized steady-state network

		* 'in':	Normalized steady-state network

	optional:
		Optional data to include. Accepts:

		* readcount:	RNA read count for each cell

	fi_c:
		Path of input tsv file for extra property columns for each cell
	"""
	from dictys.net import network
	optional=set(filter(lambda x:len(x)>0,optional.split(',')))
	n=network.from_folders(diri_data,diri_work,fi_subsets,dynamic=dynamic,nettype=nettype,optional=optional,fi_c=fi_c)
	n.to_file(fo_networks)


































#
