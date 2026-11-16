#!/usr/bin/python3
# Lingfei Wang, 2022, 2023. All rights reserved.

"""
Pyro adaptations and pyro network models used by dictys.network.reconstruct.
Kept apart from dictys.network so that command line calls that do not train a model do not import torch and pyro.
"""

import abc
from typing import Optional,Tuple,Callable,Union
from dictys.utils.importing import torch,pyro

_docstring2argparse_ignore_=['Trace_ELBO_site','SVI_multiloss','model_base','model_covariance','model_ou']


###########################################################################
# Pyro adaptations
###########################################################################

class Trace_ELBO_site(pyro.infer.trace_elbo.Trace_ELBO):
	"""
	A modified pyro.infer.trace_elbo.Trace_ELBO that traces sublosses.
	"""
	def __init__(self,keys:list[str],*a,**ka):
		"""
		A modified pyro.infer.trace_elbo.Trace_ELBO that traces sublosses. The total loss is the sum of sublosses.

		Parameters
		----------
		keys:
			List of pyro variable/factor names to indicate each subloss.
		a:
			Arguments passed to Trace_ELBO
		ka:
			Keyword arguments passed to Trace_ELBO
		"""
		super().__init__(*a,**ka)
		self.loss_keys=keys
	def loss(self, model, guide, *args, **kwargs):
		"""
		Not needed by default.
		"""
		raise NotImplementedError
	@staticmethod
	def _differentiable_loss_particle(model_trace, guide_trace):
		"""
		Adapted to output a dictionary each for loss and surraget loss.
		"""
		from collections import defaultdict
		from pyro.infer.trace_elbo import is_identically_zero,torch_item,_compute_log_r
		loss=defaultdict(float)
		surrogate_loss=defaultdict(float)
		log_r = None

		# compute elbo and surrogate elbo
		for name, site in model_trace.nodes.items():
			if site["type"] == "sample":
				loss[name]=loss[name]-torch_item(site["log_prob_sum"])
				surrogate_loss[name]=surrogate_loss[name]-site["log_prob_sum"]

		for name, site in guide_trace.nodes.items():
			if site["type"] == "sample":
				_, score_function_term, entropy_term = site["score_parts"]

				loss[name]=loss[name]+torch_item(site["log_prob_sum"])

				if not is_identically_zero(entropy_term):
					surrogate_loss[name]=surrogate_loss[name]+entropy_term.sum()

				if not is_identically_zero(score_function_term):
					if log_r is None:
						log_r = _compute_log_r(model_trace, guide_trace)
					site = log_r.sum_to(site["cond_indep_stack"])
					surrogate_loss[name]=surrogate_loss[name]-(site * score_function_term).sum()
		loss=dict(loss)
		surrogate_loss=dict(surrogate_loss)
		return loss,surrogate_loss
	def differentiable_loss(self, model, guide, *args, **kwargs):
		"""
		Adapted to output a tuple of losses.
		"""
		from collections import defaultdict
		from pyro.infer.trace_elbo import warn_if_nan,torch_item
		loss=defaultdict(float)
		surrogate_loss=defaultdict(float)
		for model_trace, guide_trace in self._get_traces(model, guide, args, kwargs):
			loss_particle, surrogate_loss_particle = self._differentiable_loss_particle(model_trace, guide_trace)
			for xi in loss_particle:
				loss[xi]=loss[xi]+loss_particle[xi]/ self.num_particles
			for xi in surrogate_loss_particle:
				surrogate_loss[xi]=surrogate_loss[xi]+surrogate_loss_particle[xi]/ self.num_particles
		surrogate_loss=dict(surrogate_loss)
		for xi in surrogate_loss:
			warn_if_nan(surrogate_loss[xi], f"loss[{xi}]")
		for xi in surrogate_loss:
			loss[xi]=loss[xi]+(surrogate_loss[xi]-torch_item(surrogate_loss[xi]))
		loss=dict(loss)
		if len(loss)!=len(self.loss_keys):
			raise ValueError('Different size between claimed and recovered losses. Please confirm the loss names defined in model.varnames match the variables.')
		loss=tuple(loss[x] for x in self.loss_keys)
		return loss
	def loss_and_grads(self, model, guide, *args, **kwargs):
		"""
		Not needed by default.
		"""
		raise NotImplementedError

class SVI_multiloss(pyro.infer.svi.SVI):
	"""Like SVI but handles and returns multiple losses. The total loss is their sum.
	"""
	@staticmethod
	def _loss_and_grads(losses,*args,**kwargs):
		"""
		Handles multiple sublosses
		"""
		if isinstance(losses,tuple):
			loss_val = tuple(x(*args, **kwargs) for x in losses)
		elif hasattr(losses,'__call__'):
			loss_val = tuple(losses(*args, **kwargs))
		else:
			assert False
		if any(getattr(x, 'requires_grad', False) for x in loss_val):
			sum(loss_val).backward(retain_graph=True)
		return loss_val
	def __init__(self,model,guide,optim,loss,loss_and_grads=None,num_samples=0,num_steps=0,**kwargs):		# pylint: disable=W0613
		import warnings
		from functools import partial
		from pyro.infer.elbo import ELBO
		if num_steps:
			warnings.warn('The `num_steps` argument to SVI is deprecated and will be removed in a future release. Use `SVI.step` directly to control the number of iterations.', FutureWarning)
		if num_samples:
			warnings.warn('The `num_samples` argument to SVI is deprecated and will be removed in a future release. Use `pyro.infer.Predictive` class to draw samples from the posterior.', FutureWarning)

		self.model = model
		self.guide = guide
		self.optim = optim
		self.num_steps = num_steps
		self.num_samples = num_samples
		super(pyro.infer.svi.SVI,self).__init__(**kwargs)

		if not isinstance(optim, pyro.optim.PyroOptim):
			raise ValueError("Optimizer should be an instance of pyro.optim.PyroOptim class.")

		if type(loss) in [list,tuple]:
			self.loss=tuple(x.loss if isinstance(x, ELBO) else x for x in loss)
		elif hasattr(loss,'__call__'):
			self.loss=loss
		else:
			assert False
		self.loss_and_grads=partial(self._loss_and_grads,self.loss)

	def evaluate_loss(self, *args, **kwargs):
		"""
		Not needed by default.
		"""
		raise NotImplementedError

###########################################################################
# Pyro models
###########################################################################

class model_base(metaclass=abc.ABCMeta):
	"""
	Abstract base class for pyro model.
	"""
	#Stats of variables to record during model training
	stat_param_funcs=[torch.min,torch.max,torch.mean,torch.median]
	stat_param_names=['min','max','mean','median']
	def __init__(self,name:str,device:str='cpu',floattype=torch.float,suffix:str=''):
		"""
		Abstract base class for pyro model.

		Parameters
		----------
		name:
			Name of model
		device:
			Default device of pytorch variables
		floattype:
			Default float type of pytorch variables
		suffix:
			Name suffix of all pyro definitions
		"""
		self.name=name
		self.device=device
		self.dtype=floattype
		self.suffix=suffix
		self.tensorka={'device':self.device,'dtype':self.dtype}
	def tensor(self,v,device:Optional[str]=None,dtype=None,**ka):
		"""
		Creates pytorch tensor with default configuration.

		Parameters
		----------
		v:
			Data to create tensor
		device:
			Device override
		dtype:
			Pytorch data type override
		ka:
			Keyword arguments passed to torch.tensor

		Returns
		-------
		torch.tensor
			Created pytorch tensor
		"""
		if device is None:
			device=self.device
		if dtype is None:
			dtype=self.dtype
		return torch.tensor(v,device=device,dtype=dtype,**ka)
	def param(self,name:str,*a,**ka):
		"""
		Finds parameter by name for this model. Automatically handles renaming.

		Parameters
		----------
		name:
			Name of parameter

		Returns
		-------
		torch.tensor
			Parameter value
		"""
		return pyro.param(name+self.suffix,*a,**ka)
	def sample(self,name:str,dist,*a,**ka):
		"""
		Performs sample for this model. Automatically handles renaming.

		Parameters
		----------
		name:
			Name of variable to sample
		dist:
			Distribution to sample

		Returns
		-------
		torch.tensor
			Sampled value
		"""
		return pyro.sample(name+self.suffix,dist,*a,**ka)
	def factor(self,name:str,*a,**ka):
		"""
		Creates pyro loss factor by name for this model. Automatically handles renaming.

		Parameters
		----------
		name:
			Name of factor
		"""
		return pyro.factor(name+self.suffix,*a,**ka)
	def plate(self,name:str,*a,**ka):
		"""
		Creates pyro plate for this model. Automatically handles renaming.

		Parameters
		----------
		name:
			Name of plate

		Returns
		-------
		Created plate
		"""
		if name in self.subsamples and 'subsample' not in ka:
			ka=dict(ka)
			ka['subsample']=self.subsamples[name]
		return pyro.plate(name+self.suffix,*a,**ka)
	@staticmethod
	def clear_param_store():
		"""
		Clears parameter storage in pyro.
		"""
		return pyro.clear_param_store()
	@staticmethod
	def get_param_store():
		"""
		Gets parameter storage from pyro.
		"""
		return pyro.get_param_store()
	def record_stat_param(self)->None:
		"""
		Stores intermediate parameters' stats for fugure inspection.
		"""
		import itertools
		for xi,xj in itertools.product(self.p,range(len(self.stat_param_funcs))):
			if self.step not in self.stat_param[xi][xj] and torch.numel(self.param(xi).data)>0:
				self.stat_param[xi][xj][self.step]=float(self.stat_param_funcs[xj](self.param(xi).data))
	def get_pyro_loss(self,loss:str,*a,**ka):
		"""
		Obtains pyro loss instance and converts loss from name str to function.

		Parameters
		----------
		loss:
			Name of loss

		Returns
		-------
		function
			Loss function
		"""
		from functools import partial
		f0=getattr(pyro.infer,loss)(*a,**ka).differentiable_loss
		func=partial(lambda f,dof,*a2,**ka2:f(*a2,**ka2)/dof,f0,self.nrand)
		return func
	def init(self,optimizer:str,loss:Union[str,Callable,list[Callable]],subsamples:dict={},optimizer_a:Tuple[list,dict]=[[],{}],loss_a:Tuple[list,dict]=[[],{}])->None:
		"""
		Initializes model for training.

		Parameters
		----------
		optimizer:
			Name of pyro optimizer class
		loss:
			Loss function as str for single loss, function for single loss, or list of functions for multiple losses
		subsamples:	{plate_name:subsample_count}
			Data subsample for each plate. Defaults to no subsample.
		optimizer_a:
			Arguments and keyword arguments for optimizer instance constructor
		optimizer_a:
			Arguments and keyword arguments for loss instance constructor
		"""
		self.optim_name=optimizer
		self.optim_a=optimizer_a
		self.loss_a=loss_a
		self.optim=getattr(pyro.optim,optimizer)(*optimizer_a[0],**optimizer_a[1])
		if isinstance(loss,str):
			self.loss=[self.get_pyro_loss(loss,*loss_a[0],**loss_a[1])]
		elif isinstance(loss,list) or hasattr(loss,'__call__'):
			self.loss=loss
		else:
			raise ValueError('Unknown type {} for parameter loss'.format(type(loss)))
		self.subsamples={}
		self.clear_param_store()
		#Initialize model parameters in pyro
		for xi in self.p:
			ka={}
			if self.p[xi][1] is not None:
				ka['constraint']=self.p[xi][1]
			self.param(xi,self.p[xi][0].to(self.device),**ka)
		self.guide()
		self.subsamples=subsamples

		#Initialize training stats for parameters
		self.stat_param={x:[{} for _ in self.stat_param_funcs] for x in self.p}
		self.stat_loss=[]
		self.step=0
		self.record_stat_param()
	def save_param(self,path:str)->None:
		"""
		Save parameters to file.
		"""
		self.get_param_store().save(path)
	def load_param(self,path:str)->None:
		"""
		Load parameters from file.
		"""
		self.get_param_store().load(path)
	def train_svi(self,nstep:int,nstep_report:int=100):
		"""
		Train model with stochastic variational inference.

		Parameters
		----------
		nstep:
			Number of training steps
		nreportsteps:
			Number of training steps between each report of parameter stats
		"""
		if not hasattr(self,'svi') or self.svi is None:		# pylint: disable=E0203
			self.svi=SVI_multiloss(self.model,self.guide,self.optim,loss=self.loss)

		target=self.step+nstep
		while self.step<target:
			self.stat_loss.append(self.svi.step())
			self.step+=1
			if self.step%nstep_report==0:
				self.record_stat_param()
		#Store last training info
		self.last_train_params=['train_svi',[nstep],{'nstep_report':nstep_report}]
	def retrain(self):
		"""
		Train again with the same configuration as last.
		"""
		return getattr(self,self.last_train_params[0])(*self.last_train_params[1],**self.last_train_params[2])
	def draw_stats(self,lossnames:Optional[list[str]]=None)->None:
		"""
		Visualize parameter stats with matplotlib.

		Parameters
		----------
		lossnames:
			List of names for sublosses to draw. Defaults to all.

		"""
		import matplotlib.pyplot as plt
		import numpy as np
		losses2=np.array(self.stat_loss)
		if lossnames is None:
			lossnames=list(range(losses2.shape[1]))
		nstart=int(len(losses2)*0.2)
		nend=len(losses2)
		for xi in range(losses2.shape[1]):
			plt.plot(np.arange(nstart,nend),losses2[nstart:nend,xi])
			plt.ylabel('Loss {}'.format(lossnames[xi]))
			plt.show()
		plt.plot(np.arange(nstart,nend),losses2[nstart:nend].sum(axis=1))
		plt.ylabel('Total loss')
		plt.show()
		for xi in self.p:
			if np.sum([len(x) for x in self.stat_param[xi]])==0:
				continue
			for xj in range(len(self.stat_param[xi])):
				t1=np.array(list(self.stat_param[xi][xj].items())).T
				plt.plot(t1[0],t1[1],label=self.stat_param_names[xj])
			plt.legend()
			title=xi+' '+str(list(self.param(xi).shape))
			if self.p[xi][1] is not None:
				title+=' '+str(self.p[xi][1])
			plt.title(title)
			plt.show()
	def gen_params_preproc(self,observations:dict)->dict:
		"""
		Preprocessing steps to initialize model parameters to account for suffix.

		Parameters
		----------
		observations:	{name:torch.tensor}
			Observations as input of model training

		Returns
		-------
		{name:torch.tensor}
			Preprocessed observations to account for name suffix
		"""
		if len(self.suffix)==0:
			return observations
		assert all(x.endswith(self.suffix) for x in observations)
		observations={x[:-len(self.suffix)]:y for x,y in observations.items()}
		return observations
	@abc.abstractmethod
	def gen_params(self,observations:dict)->None:
		"""
		Prepare and save model parameter initialization values. This function should save them in three variables:

		* self.a:	Tuple of model hyperparameter values (including observations) which are read by model & guide functions.

		* self.p:	Dictionary of initialization values of each self.param to optimize. Must match self.param.

		* self.nrand:	Number of random variables for loss function scaling between sublosses

		Parameters
		----------
		observations:	{name:torch.tensor}
			Observations as input of model training
		"""
	@abc.abstractmethod
	def model(self):
		"""
		Pyro model function to compute likelihood or other loss. Should read in data from those placed by self.gen_params.
		"""
	@abc.abstractmethod
	def guide(self)->None:
		"""
		Pyro guide function to compute posterior likelihood or other loss. Should read in data from those placed by self.gen_params.
		"""

class model_covariance(model_base):
	"""
	Model of scRNA-seq data with a low-rank multi-variate normal distribution of true expression level, followed by binomial sampling to obtain RNA read counts.
	"""
	#Sublosses by pyro variable/factor name
	varnames=['G_0','G_obs']
	def __init__(self,name:str,npc:int,dcs,covdc_rate:float=1.,**ka):
		"""
		Model of scRNA-seq data with a low-rank multi-variate normal distribution of true expression level, followed by binomial sampling to obtain RNA read counts.
		This also accounts for linear and nonlinear technical confounders for each cell such as expressed gene count and log read count squared.
		See ... for model details.

		Parameters
		----------
		name:
			Name of model
		npc:
			Number of low rank off-diagonal degrees of freedom in multivariate normal distribution. See....
		dcs:
			Covariates to account for in shape (n_cell,n_cov)
		covdc_rate:
			Initial estimate for the ratio of contribution to offdiangonal covariance matrix from covariates v.s. low-dimensional factors. Used for model initialization.
		ka:
			Keyword arguments passed to base model
		"""
		super().__init__(name,**ka)
		assert dcs.ndim==2
		self.npc=npc
		self.dcs=self.tensor(dcs)
		#Remove single-valued covariates
		if self.dcs.shape[1]>0:
			self.dcs=self.dcs[:,self.dcs.max(axis=0)!=self.dcs.min(axis=0)]
		self.nl=self.dcs.shape[1]
		self.covdc_rate=self.tensor(covdc_rate)
	def gen_params(self,observations:dict)->None:
		from pyro.distributions import constraints
		observations=self.gen_params_preproc(observations)
		assert isinstance(observations,dict) and len(observations)==1 and 'G_obs' in observations
		# observations={x:y.to(self.device) for x,y in observations.items()}
		self.obs=observations
		dx=self.obs['G_obs']
		self.ng,self.nc=dx.shape

		##Arguments
		n_read=dx.sum(axis=0)
		obs_G_obs=dx.T
		dc=self.dcs
		if dc.shape[0]>0:
			#Orthonormalize covariates
			dc-=dc.mean(axis=0)
			dc/=(dc**2).mean(axis=0).sqrt()+1E-30
			dc=torch.pca_lowrank(dc,q=self.nl,center=False,niter=50)[0]

		##Initial values
		t1=dx.sum(axis=1)
		t1=torch.log(t1/t1.sum())
		mu_G_init=t1
		t1=torch.log(dx/dx.sum(axis=0)+1E-6).T
		t1-=t1.mean(axis=0)
		t2=torch.pca_lowrank(t1,q=self.npc,center=False,niter=50)
		sigma_G_nd_init=t2[1]*t2[2]/torch.sqrt(self.tensor(self.nc))
		if self.npc==0:
			scale_pc=torch.pca_lowrank(t1,q=1,center=False,niter=50)
			scale_pc=scale_pc[1]*scale_pc[2]/torch.sqrt(self.tensor(self.nc))
		else:
			scale_pc=sigma_G_nd_init
		scale_pc=torch.mm(scale_pc,scale_pc.T).abs().mean()
		sigma_G_d_init=(t1-torch.mm(t2[0],(t2[1]*t2[2]).T))
		sigma_G_d_init=torch.sqrt((sigma_G_d_init**2).mean(axis=0))
		G_0_mean_init=mu_G_init.expand([self.nc,self.ng])
		G_0_std_init=torch.ones([self.nc,self.ng],**self.tensorka)*self.tensor(0.5)
		self.p={
			'mu_G':[mu_G_init,constraints.interval(-20.,0.)],'sigma_G_d':[sigma_G_d_init,constraints.positive],
			'sigma_G_nd':[sigma_G_nd_init,None],'G_0_mean':[G_0_mean_init,None],
			'G_0_std':[G_0_std_init,constraints.positive],
		}
		if self.nl>0:
			gamma_init=torch.zeros([self.nl,self.ng],**self.tensorka)
			self.p['gamma']=[gamma_init,None]
		
		##Arguments
		if self.nl>0:
			#Scale covariates to same contribution level on covariance matrix
			t1=torch.sqrt(self.covdc_rate*scale_pc*self.tensor(self.nc/self.nl))
			dc*=t1
			dc=dc.T
		self.a=(n_read,obs_G_obs,dc)
		self.nrand=dx.shape[0]*dx.shape[1]
		
	_binom_cls=None
	def _binomial_cached(self,n_read,obs_G_obs,logits):
		"""
		Binomial over the full cell plate whose log_prob reuses the data-only log binomial coefficient, computed once on first use.
		Arithmetic order follows pyro.distributions.Binomial.log_prob (approx_log_prob_tol=0).

		Derived from Binomial.log_prob in Pyro 1.9.2 (pyro/distributions/torch.py), Copyright (c) 2017-2019 Uber Technologies, Inc., licensed under the Apache License, Version 2.0 (http://www.apache.org/licenses/LICENSE-2.0; the license text is in Pyro's LICENSE.md). Pyro 1.9.2 ships no NOTICE file.
		Changed from the original: the log binomial coefficient is computed once from the data and cached instead of on every call, and log_prob is defined on a subclass built on first use.
		Unless required by applicable law or agreed to in writing, Pyro is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
		"""
		import torch
		import pyro.distributions as dist
		if dist.Binomial.approx_log_prob_tol>=0.02:
			return dist.Binomial(n_read,logits=logits)
		cache=getattr(self,'_binom_cache',None)
		if cache is None or cache.shape!=logits.shape or cache.dtype!=logits.dtype:
			with torch.no_grad():
				n=dist.Binomial(n_read,logits=logits).total_count
				k=obs_G_obs.T
				n_plus_1=n+1
				cache=n_plus_1.lgamma()-(k+1).lgamma()-(n_plus_1-k).lgamma()
			self._binom_cache=cache
		cls=self._binom_cls
		if cls is None:
			class cls(dist.Binomial):
				def log_prob(self,value):
					if self._validate_args:
						self._validate_sample(value)
					n=self.total_count
					k=value
					x=self.logits
					# Same as pyro's _clamp_by_zero
					clamp=(x.clamp(min=0)+x-x.clamp(max=0))/2
					normalize_term=n*(clamp+x.abs().neg().exp().log1p())
					return k*x-normalize_term+type(self)._cache
			self._binom_cls=cls
		cls._cache=cache
		return cls(n_read,logits=logits)
	def model(self)->None:
		import pyro.distributions as dist
		n_read,obs_G_obs,dc=self.a
		
		mu_G=self.param("mu_G")
		sigma_G_d=self.param("sigma_G_d")
		sigma_G_nd=self.param("sigma_G_nd")
		if self.nl>0:
			gamma=self.param("gamma")
		with self.plate('cell1',self.nc,dim=-1) as ind_c:
			t1=dist.LowRankMultivariateNormal(mu_G,sigma_G_nd,sigma_G_d)
			G_0=self.sample("G_0",t1)
			if self.nl>0:
				G_C=torch.mm(dc.T[ind_c],gamma)
				G_1=G_0+G_C
			else:
				G_1=G_0
			with self.plate('lv2',self.ng,dim=-2) as ind_g:
				if len(ind_c)==self.nc:
					t1=self._binomial_cached(n_read,obs_G_obs,G_1.T)
				else:
					t1=dist.Binomial(n_read[ind_c],logits=G_1.T)
				G_obs=self.sample("G_obs",t1,obs=obs_G_obs[ind_c].T)
				assert G_0.shape in {(len(ind_c),len(ind_g)),(self.nc,len(ind_g))}
				assert G_obs.shape==(len(ind_g),len(ind_c))

	def guide(self)->None:
		import pyro.distributions as dist
		G_0_mean=self.param('G_0_mean')
		G_0_std=self.param('G_0_std')
		with self.plate('cell1',self.nc,dim=-1) as ind_c:
			t1=dist.Normal(G_0_mean[ind_c],G_0_std[ind_c]).to_event(1)
			self.sample("G_0",t1)

class model_ou(model_covariance):
	"""
	On top of model_covariance, also models the source of covariance matrix with the steady-state distribution of Ornstein-Uhlenbeck process.
	"""
	varnames=['G_0','sigma_G_obs','G_obs']
	def __init__(self,name:str,N_0mask,fullsize:Tuple[int,int],npc0:int,dcs,scale_lyapunov:float=1000,npc:Optional[int]=None,**ka):
		"""
		On top of model_covariance, also models the source of covariance matrix with the steady-state distribution of Ornstein-Uhlenbeck process.
		See ... for model details.

		Parameters
		----------
		name:
			Name of model
		N_0mask:
			Mask of network (N_0) that allows nonzeros values of edge strength in shape [2,n_edge]
		fullsize:
			Full size of network as [n_reg,n_target]
		npc0:
			Number of low rank off-diagonal degrees of freedom in stochastic differential equation
		scale_lyapunov:
			Scale of Lyapunov equation inequality loss.
		npc:
			Number of low-dimensional factors for true expression covariance matrix sigma_G. Defaults to min(npc+n_reg,n_target-1). This is in theory sufficient for the full covariance structure generated by network SDE steady-state distribution.
		ka:
			Keyword arguments passed to covariance model
		"""
		self.ntf,self.ng=fullsize
		assert N_0mask.ndim==2 and N_0mask.shape[0]==2 and N_0mask.shape[1]>0
		assert (N_0mask[0]!=N_0mask[1]).all()
		self.npc0=npc0
		npcmax=min(self.ntf+self.npc0,self.ng-1)
		if npc is None:
			npc=npcmax
		elif npc>npcmax:
			raise ValueError('Need npc>=npc0+number of TFs.')
		if npc<npc0:
			raise ValueError('Not enough target genes for specified {} principal components'.format(npc0))
		super().__init__(name,npc,dcs,**ka)
		self.N_0mask=self.tensor(N_0mask,dtype=torch.long)
		self.nedge=self.N_0mask.shape[1]
		self.scale_lyapunov=scale_lyapunov
	def gen_params(self,observations:dict)->None:
		from pyro.distributions import constraints
		super().gen_params(observations)
		observations=self.gen_params_preproc(observations)
		assert observations['G_obs'].ndim==2 and observations['G_obs'].shape[0]==self.ng

		N_0val_init=torch.normal(0,0.0001,size=[self.nedge],**self.tensorka)
		reg_init=self.tensor(1)
		self.p['N_0val']=[N_0val_init,constraints.interval(-5.,5.)]
		self.p['reg']=[reg_init,constraints.interval(-10,10)]
		#Reuse top PCs for initialization
		self.p['sigma0_G_d']=list(self.p['sigma_G_d'])
		self.p['sigma0_G_nd']=[self.p['sigma_G_nd'][0][:,:self.npc0],None]
	def lyapunov_lhs0(self,N_0val):
		"""
		Left hand side of Lyapunov equation.
		"""
		sigma_G_d,sigma_G_nd,reg=[self.param(x) for x in 'sigma_G_d,sigma_G_nd,reg'.split(',')]
		N_0=torch.zeros((self.ng,self.ng),**self.tensorka)
		N_0[self.N_0mask[1],self.N_0mask[0]]=N_0val
		matA=torch.eye(self.ng,**self.tensorka)*reg-N_0
		matX=torch.diag_embed(sigma_G_d)+sigma_G_nd@sigma_G_nd.T
		ans=matA@matX
		return ans+ans.T
	def lyapunov_lhs(self):
		"""
		Left hand side of Lyapunov equation.
		"""
		return self.lyapunov_lhs0(self.param('N_0val'))
	def lyapunov_rhs(self):
		"""
		Right hand side of Lyapunov equation.
		"""
		sigma0_G_d,sigma0_G_nd=[self.param(x) for x in 'sigma0_G_d,sigma0_G_nd'.split(',')]
		return torch.diag_embed(sigma0_G_d)+sigma0_G_nd@sigma0_G_nd.T
	def model(self)->None:
		super().model()
		#MSE loss for lyapunov
		lossl=self.lyapunov_lhs()-self.lyapunov_rhs()
		lossl=(lossl**2).sum()*((self.scale_lyapunov*self.nc)/(4*self.ng))
		self.factor("sigma_G_obs",-lossl)
