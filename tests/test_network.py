"""
Tests for dictys.network: the pyro loss and SVI adaptations, the pyro network models, reconstruct, indirect,
normalize and tofile, and the `dictys network` subcommands.

Expected values come from closed forms (SGD on quadratic losses, the first Adam step), scipy densities, numpy
re-derivations of the formulas in the source and docstrings, and the outputs of a reconstruct run with lr=0,
which leaves every parameter at its initial value.
"""

import os

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from dictys import network

###########################################################################
# Shared helpers and fixtures
###########################################################################

OUTPUTS = ['weight', 'meanvar', 'covfactor', 'loss', 'stats']
GENES = ['tf', 'up', 'down', 'ind', 'tf2']
REGS = ['tf', 'tf2']
ORDER = ['tf', 'tf2', 'up', 'down', 'ind']
STATS = ['min', 'max', 'mean', 'median']


def read_tsv(path):
	"""Read a tab-separated table with a header row and an index column."""
	return pd.read_csv(path, header=0, index_col=0, sep='\t')


def write_lines(path, lines):
	"""Write one item per line."""
	with open(path, 'w') as f:
		f.write('\n'.join(lines) + '\n')


def lower_median(x):
	"""Median as torch.median computes it: the lower of the two middle values for an even count."""
	x = np.sort(np.ravel(x))
	return x[(len(x) - 1) // 2]


def summary(x):
	"""The four statistics model_base records for a parameter."""
	x = np.ravel(x)
	return np.array([x.min(), x.max(), x.mean(), lower_median(x)])


def pca_init(counts, k):
	"""
	Initial mu_G, sigma_G_d and sigma_G_nd of model_covariance recomputed with numpy's exact SVD from a
	genes x cells count matrix.
	"""
	counts = np.asarray(counts, dtype=float)
	tot = counts.sum(axis=1)
	mu = np.log(tot / tot.sum())
	x = np.log(counts / counts.sum(axis=0) + 1E-6).T
	x = x - x.mean(axis=0)
	u, s, vt = np.linalg.svd(x, full_matrices=False)
	factor = vt[:k].T * s[:k] / np.sqrt(x.shape[0])
	res = x - (u[:, :k] * s[:k]) @ vt[:k]
	return mu, np.sqrt((res ** 2).mean(axis=0)), factor


@pytest.fixture
def tp():
	"""torch and pyro modules, skipping without them; pyro's parameter store is cleared around the test."""
	torch = pytest.importorskip('torch')
	pyro = pytest.importorskip('pyro')
	pyro.clear_param_store()
	pyro.set_rng_seed(0)
	yield torch, pyro
	pyro.clear_param_store()


@pytest.fixture
def threads():
	"""Skip unless joblib and threadpoolctl, which network functions use to limit threads, are installed."""
	pytest.importorskip('joblib')
	pytest.importorskip('threadpoolctl')


def simulate_counts():
	"""
	Read counts of 5 genes in 40 cells: 'up' follows the log expression of 'tf', 'down' mirrors it, 'ind' and
	'tf2' are independent.
	"""
	rng = np.random.default_rng(7)
	nc = 40
	tf = rng.normal(0, 0.8, nc)
	logx = np.array([tf, tf + rng.normal(0, 0.2, nc), -tf + rng.normal(0, 0.2, nc), rng.normal(0, 0.8, nc), rng.normal(0, 0.8, nc)])
	p = np.exp(logx)
	p /= p.sum(axis=0)
	return pd.DataFrame(rng.poisson(2000 * p), index=GENES, columns=['c' + str(x) for x in range(nc)])


@pytest.fixture(scope='module')
def sim(tmp_path_factory):
	"""Simulated expression and mask files for reconstruct (mask columns in GENES order, regulators tf and tf2)."""
	pytest.importorskip('torch')
	pytest.importorskip('pyro')
	pytest.importorskip('joblib')
	pytest.importorskip('threadpoolctl')
	d = str(tmp_path_factory.mktemp('sim'))
	expr = simulate_counts()
	mask = pd.DataFrame(0, index=REGS, columns=GENES)
	mask.loc['tf', ['up', 'down', 'ind']] = 1
	mask.loc['tf2', ['ind', 'up']] = 1
	ans = {'exp': os.path.join(d, 'expression.tsv'), 'mask': os.path.join(d, 'mask.tsv'), 'expr': expr, 'maskdf': mask}
	expr.to_csv(ans['exp'], sep='\t')
	mask.to_csv(ans['mask'], sep='\t')
	ans['allowed'] = mask.loc[REGS, ORDER].values.astype(bool)
	return ans


def output_paths(folder, prefix=''):
	"""Paths of reconstruct's five outputs in a folder."""
	return [os.path.join(str(folder), prefix + x + '.tsv') for x in OUTPUTS]


def run_reconstruct(sim, folder, prefix='', exp=None, **ka):
	"""Run reconstruct on the simulated data (nth=1 unless given) and return its five output tables."""
	paths = output_paths(folder, prefix)
	ka.setdefault('nth', 1)
	network.reconstruct(sim['exp'] if exp is None else exp, sim['mask'], *paths, **ka)
	return dict(zip(OUTPUTS, [read_tsv(x) for x in paths]))


@pytest.fixture(scope='module')
def sim_lr0(sim, tmp_path_factory):
	"""reconstruct outputs after 3 steps with lr=0, i.e. the initial parameter values (npc=0, default seed)."""
	return run_reconstruct(sim, tmp_path_factory.mktemp('lr0'), nstep=3, lr=0)


def assert_outputs_equal(a, b):
	"""All five reconstruct output tables are identical."""
	for xi in OUTPUTS:
		pd.testing.assert_frame_equal(a[xi], b[xi], check_exact=True)


def lyapunov_loss(weight, var, factor, npc, scale, ncell):
	"""
	model_ou's Lyapunov loss sum((A X + X A^T - Q)^2) * scale * ncell / (4 * ngene), with A = reg*I - N_0
	(reg=1, N_0[target, regulator] = weight), X = diag(var) + F F^T and Q = diag(var) + F[:, :npc] F[:, :npc]^T.
	"""
	ng = weight.shape[1]
	n0 = np.zeros((ng, ng))
	n0[:weight.shape[0]] = weight
	a = np.eye(ng) - n0.T
	x = np.diag(var) + factor @ factor.T
	q = np.diag(var) + factor[:, :npc] @ factor[:, :npc].T
	lhs = a @ x + x @ a.T
	return ((lhs - q) ** 2).sum() * scale * ncell / (4 * ng)


###########################################################################
# Trace_ELBO_site
###########################################################################

def delta_model(torch, pyro, y0=2.):
	"""
	z ~ Normal(0,1), y ~ Normal(z,1) observed at y0, with a Delta(m) guide at m=0.5: the per-site losses are
	-log N(m; 0, 1) for z and -log N(y0; m, 1) for y, and d(total)/dm = 2m - y0.
	"""
	import pyro.distributions as dist

	def model():
		z = pyro.sample('z', dist.Normal(torch.tensor(0.), torch.tensor(1.)))
		pyro.sample('y', dist.Normal(z, torch.tensor(1.)), obs=torch.tensor(y0))

	def guide():
		m = pyro.param('m', torch.tensor(0.5))
		pyro.sample('z', dist.Delta(m))
	return model, guide


def test_trace_elbo_site_returns_per_site_losses_in_key_order(tp):
	"""differentiable_loss returns one loss per site in the order of keys, matching scipy's normal densities."""
	from scipy.stats import norm
	torch, pyro = tp
	model, guide = delta_model(torch, pyro)
	ans = network.Trace_ELBO_site(['y', 'z']).differentiable_loss(model, guide)
	assert isinstance(ans, tuple) and len(ans) == 2
	assert_allclose([float(x.detach()) for x in ans], [-norm.logpdf(2., loc=0.5), -norm.logpdf(0.5)], rtol=1E-6)
	assert all(x.requires_grad for x in ans)
	ans2 = network.Trace_ELBO_site(['z', 'y']).differentiable_loss(model, guide)
	assert_allclose([float(x.detach()) for x in ans2], [float(x.detach()) for x in ans[::-1]], rtol=1E-7)


def test_trace_elbo_site_forwards_elbo_arguments(tp):
	"""Extra arguments reach Trace_ELBO: num_particles is set and a deterministic loss is unchanged."""
	from scipy.stats import norm
	torch, pyro = tp
	model, guide = delta_model(torch, pyro)
	elbo = network.Trace_ELBO_site(['y', 'z'], num_particles=3)
	assert elbo.num_particles == 3
	ans = elbo.differentiable_loss(model, guide)
	assert_allclose([float(x.detach()) for x in ans], [-norm.logpdf(2., loc=0.5), -norm.logpdf(0.5)], rtol=1E-6)


def test_trace_elbo_site_rejects_wrong_key_count(tp):
	"""A key list whose length differs from the number of sites raises ValueError."""
	torch, pyro = tp
	model, guide = delta_model(torch, pyro)
	with pytest.raises(ValueError):
		network.Trace_ELBO_site(['z']).differentiable_loss(model, guide)
	with pytest.raises(ValueError):
		network.Trace_ELBO_site(['y', 'z', 'w']).differentiable_loss(model, guide)


def test_trace_elbo_site_loss_methods_not_implemented(tp):
	"""loss and loss_and_grads raise NotImplementedError."""
	torch, pyro = tp
	model, guide = delta_model(torch, pyro)
	elbo = network.Trace_ELBO_site(['y', 'z'])
	with pytest.raises(NotImplementedError):
		elbo.loss(model, guide)
	with pytest.raises(NotImplementedError):
		elbo.loss_and_grads(model, guide)


###########################################################################
# SVI_multiloss
###########################################################################

def test_svi_multiloss_step_with_site_losses(tp):
	"""A step returns the per-site losses as floats and moves m by -lr * d(sum)/dm = -0.1 * (2*0.5 - 2)."""
	from scipy.stats import norm
	torch, pyro = tp
	model, guide = delta_model(torch, pyro)
	svi = network.SVI_multiloss(model, guide, pyro.optim.SGD({'lr': 0.1}), loss=network.Trace_ELBO_site(['y', 'z']).differentiable_loss)
	ans = svi.step()
	assert isinstance(ans, tuple) and all(isinstance(x, float) for x in ans)
	assert_allclose(ans, [-norm.logpdf(2., loc=0.5), -norm.logpdf(0.5)], rtol=1E-6)
	assert_allclose(float(pyro.param('m').detach()), 0.6, rtol=1E-6)


def test_svi_multiloss_step_with_tuple_of_losses(tp):
	"""A tuple of loss functions gives one value each, and the gradient step uses their sum."""
	torch, pyro = tp
	model, guide = delta_model(torch, pyro)

	def f1(model, guide):
		return (pyro.param('m', torch.tensor(0.5)) - 1) ** 2

	def f2(model, guide):
		return (pyro.param('m', torch.tensor(0.5)) - 3) ** 2
	svi = network.SVI_multiloss(model, guide, pyro.optim.SGD({'lr': 0.1}), loss=(f1, f2))
	assert svi.loss == (f1, f2)
	assert_allclose(svi.step(), [0.25, 6.25], rtol=1E-6)
	# d/dm[(m-1)^2+(m-3)^2] at 0.5 is -6
	assert_allclose(float(pyro.param('m').detach()), 1.1, rtol=1E-6)


def test_svi_multiloss_uses_loss_method_of_elbo_instances(tp):
	"""An ELBO in the loss list is replaced by its .loss, whose float value is returned without a parameter update."""
	from scipy.stats import norm
	torch, pyro = tp
	model, guide = delta_model(torch, pyro)
	elbo = pyro.infer.Trace_ELBO()
	svi = network.SVI_multiloss(model, guide, pyro.optim.SGD({'lr': 0.1}), loss=[elbo])
	assert svi.loss == (elbo.loss,)
	ans = svi.step()
	assert_allclose(ans, [-norm.logpdf(2., loc=0.5) - norm.logpdf(0.5)], rtol=1E-6)
	assert_allclose(float(pyro.param('m').detach()), 0.5)


def test_svi_multiloss_rejects_invalid_arguments(tp):
	"""A non-pyro optimizer raises ValueError, a loss that is not callable or a sequence raises AssertionError."""
	torch, pyro = tp
	model, guide = delta_model(torch, pyro)
	loss = network.Trace_ELBO_site(['y', 'z']).differentiable_loss
	with pytest.raises(ValueError):
		network.SVI_multiloss(model, guide, torch.optim.SGD([torch.zeros(1, requires_grad=True)], lr=0.1), loss=loss)
	with pytest.raises(AssertionError):
		network.SVI_multiloss(model, guide, pyro.optim.SGD({'lr': 0.1}), loss=5)


def test_svi_multiloss_deprecated_arguments_warn(tp):
	"""num_steps and num_samples give FutureWarning and are stored."""
	torch, pyro = tp
	model, guide = delta_model(torch, pyro)
	loss = network.Trace_ELBO_site(['y', 'z']).differentiable_loss
	with pytest.warns(FutureWarning):
		svi = network.SVI_multiloss(model, guide, pyro.optim.SGD({'lr': 0.1}), loss=loss, num_steps=3)
	assert svi.num_steps == 3
	with pytest.warns(FutureWarning):
		svi = network.SVI_multiloss(model, guide, pyro.optim.SGD({'lr': 0.1}), loss=loss, num_samples=2)
	assert svi.num_samples == 2


def test_svi_multiloss_evaluate_loss_not_implemented(tp):
	"""evaluate_loss raises NotImplementedError."""
	torch, pyro = tp
	model, guide = delta_model(torch, pyro)
	svi = network.SVI_multiloss(model, guide, pyro.optim.SGD({'lr': 0.1}), loss=network.Trace_ELBO_site(['y', 'z']).differentiable_loss)
	with pytest.raises(NotImplementedError):
		svi.evaluate_loss()


###########################################################################
# model_base, through a quadratic toy model
###########################################################################

TOY_Y = [1., 2., 3., 4.]


def toy_class(torch):
	"""A model_base subclass whose loss is the factor f = sum((mu - y)^2) over a plate 'data' of 4 entries."""
	from pyro.distributions import constraints

	class Toy(network.model_base):
		"""
		Toy model: 'mu' (4 values from 0) is trained, 's' (positive) and 'empty' (no values) are never used.
		"""
		def gen_params(self, observations):
			observations = self.gen_params_preproc(observations)
			self.y = observations['y']
			self.p = {
				'mu': [torch.zeros(4, **self.tensorka), None],
				's': [self.tensor([1., 2.]), constraints.positive],
				'empty': [torch.zeros(0, **self.tensorka), None]}
			self.a = ()
			self.nrand = 2

		def model(self):
			mu = self.param('mu')
			with self.plate('data', 4) as ind:
				self.factor('f', -(mu[ind] - self.y[ind]) ** 2)

		def guide(self):
			pass
	return Toy


def make_toy(torch, suffix='', loss=None, lr=0.1, **ka):
	"""A fresh toy model with parameters initialised for SGD at the given learning rate."""
	m = toy_class(torch)('toy', suffix=suffix)
	m.gen_params({'y' + suffix: m.tensor(TOY_Y)})
	if loss is None:
		loss = network.Trace_ELBO_site(['f' + suffix]).differentiable_loss
	m.init('SGD', loss, optimizer_a=[[{'lr': lr}], {}], **ka)
	return m


def sgd_mu(step, lr=0.1):
	"""mu after `step` SGD steps on sum((mu-y)^2) from 0: y * (1 - (1-2 lr)^step)."""
	return np.array(TOY_Y) * (1 - (1 - 2 * lr) ** step)


def sgd_loss(step, lr=0.1):
	"""sum((mu-y)^2) before SGD step `step`: 30 * (1-2 lr)^(2 step)."""
	return 30 * (1 - 2 * lr) ** (2 * step)


def recorded_losses(m):
	"""The per-step losses train_svi recorded (what reconstruct writes to fo_loss), as a steps x losses array."""
	return np.array(m.stat_loss, dtype=float)


def recorded_stats(m, name):
	"""{step: [min, max, mean, median]} recorded for a parameter (what reconstruct writes to fo_stats)."""
	steps = sorted(m.stat_param[name][0])
	return {x: [m.stat_param[name][y][x] for y in range(len(STATS))] for x in steps}


def test_model_base_is_abstract(tp):
	"""model_base and subclasses missing model/guide/gen_params cannot be instantiated."""
	with pytest.raises(TypeError):
		network.model_base('x')

	class Partial(network.model_base):
		"""Subclass without guide."""
		def gen_params(self, observations):
			pass

		def model(self):
			pass
	with pytest.raises(TypeError):
		Partial('x')


def test_model_base_tensor_defaults_and_overrides(tp):
	"""tensor() uses the model's dtype and device unless overridden, and passes other keywords to torch.tensor."""
	torch, pyro = tp
	m = toy_class(torch)('toy', floattype=torch.float64)
	t = m.tensor([1, 2])
	assert t.dtype == torch.float64 and t.device.type == 'cpu'
	assert_array_equal(t.numpy(), [1., 2.])
	assert m.tensor([1, 2], dtype=torch.int32).dtype == torch.int32
	assert m.tensor([1.], device='cpu', requires_grad=True).requires_grad


def test_model_base_gen_params_preproc(tp):
	"""Observation names lose the suffix; without a suffix the dict is returned as is; a missing suffix fails."""
	torch, pyro = tp
	Toy = toy_class(torch)
	obs = {'y': 1}
	assert Toy('toy').gen_params_preproc(obs) is obs
	assert Toy('toy', suffix='_s').gen_params_preproc({'y_s': 1, 'G_s': 2}) == {'y': 1, 'G': 2}
	with pytest.raises(AssertionError):
		Toy('toy', suffix='_s').gen_params_preproc({'y': 1})


def test_model_base_suffix_names_pyro_sites(tp):
	"""param, factor and plate append the suffix to pyro names; param reads the store by the suffixed name."""
	torch, pyro = tp
	from pyro import poutine
	m = make_toy(torch, suffix='_s')
	assert set(m.get_param_store().keys()) == {'mu_s', 's_s', 'empty_s'}
	assert tuple(m.param('mu').shape) == (4,)
	tr = poutine.trace(m.model).get_trace()
	assert 'f_s' in tr.nodes and 'data_s' in tr.nodes and 'mu_s' in tr.nodes
	assert 'f' not in tr.nodes


def test_model_base_sample_uses_suffix(tp):
	"""sample() draws under the suffixed name and returns the drawn value."""
	torch, pyro = tp
	from pyro import poutine
	import pyro.distributions as dist
	m = toy_class(torch)('toy', suffix='_s')
	tr = poutine.trace(lambda: m.sample('z', dist.Normal(torch.tensor(0.), torch.tensor(1.)))).get_trace()
	assert 'z_s' in tr.nodes and 'z' not in tr.nodes
	assert tr.nodes['_RETURN']['value'] is tr.nodes['z_s']['value']


def test_model_base_init_registers_parameters_and_stats(tp):
	"""init registers every parameter with its constraint and records the step-0 statistics of non-empty ones."""
	torch, pyro = tp
	from pyro.distributions import constraints
	m = make_toy(torch)
	store = m.get_param_store()
	assert store is pyro.get_param_store()
	assert set(store.keys()) == {'mu', 's', 'empty'}
	assert store.get_state()['constraints']['s'] is constraints.positive
	assert recorded_losses(m).size == 0
	assert recorded_stats(m, 'mu') == {0: [0, 0, 0, 0]}
	assert list(recorded_stats(m, 's')) == [0]
	assert_allclose(recorded_stats(m, 's')[0], [1, 2, 1.5, 1], rtol=1E-6)
	assert recorded_stats(m, 'empty') == {}


def test_model_base_train_svi_follows_sgd(tp):
	"""Losses, parameters and recorded statistics follow SGD on sum((mu-y)^2): mu_t = y (1 - 0.8^t)."""
	torch, pyro = tp
	m = make_toy(torch)
	m.train_svi(5, nstep_report=2)
	losses = recorded_losses(m)
	assert losses.shape == (5, 1)
	assert_allclose(losses[:, 0], [sgd_loss(x) for x in range(5)], rtol=1E-5)
	assert_allclose(m.param('mu').detach().numpy(), sgd_mu(5), rtol=1E-5)
	assert_allclose(m.param('s').detach().numpy(), [1, 2], rtol=1E-6)
	stats = recorded_stats(m, 'mu')
	assert list(stats) == [0, 2, 4]
	for step in [2, 4]:
		assert_allclose(stats[step], summary(sgd_mu(step)), rtol=1E-5)


def test_model_base_retrain_continues_training(tp):
	"""retrain runs the same number of steps again, continuing the SGD trajectory and the step count."""
	torch, pyro = tp
	m = make_toy(torch)
	m.train_svi(3, nstep_report=2)
	m.retrain()
	assert_allclose(recorded_losses(m)[:, 0], [sgd_loss(x) for x in range(6)], rtol=1E-5)
	assert_allclose(m.param('mu').detach().numpy(), sgd_mu(6), rtol=1E-5)
	assert list(recorded_stats(m, 'mu')) == [0, 2, 4, 6]


def test_model_base_record_stat_param(tp):
	"""record_stat_param keeps the first record of a step and records the current values at a new step."""
	torch, pyro = tp
	m = make_toy(torch)
	m.get_param_store()['mu'] = torch.tensor([5., 6., 7., 9.])
	m.record_stat_param()
	assert recorded_stats(m, 'mu') == {0: [0, 0, 0, 0]}
	m.train_svi(3, nstep_report=100)
	assert list(recorded_stats(m, 'mu')) == [0]
	m.get_param_store()['mu'] = torch.tensor([5., 6., 7., 9.])
	m.record_stat_param()
	assert_allclose(recorded_stats(m, 'mu')[3], [5, 9, 6.75, 6])
	assert recorded_stats(m, 'empty') == {}


def test_model_base_get_pyro_loss_divides_by_nrand(tp):
	"""A pyro loss name gives Trace_ELBO's loss divided by nrand=2, and SGD then follows (mu-y)^2/2."""
	torch, pyro = tp
	m = make_toy(torch, loss='Trace_ELBO')
	assert isinstance(m.loss, list) and len(m.loss) == 1
	assert_allclose(float(m.loss[0](m.model, m.guide).detach()), 15., rtol=1E-6)
	f = m.get_pyro_loss('Trace_ELBO')
	assert_allclose(float(f(m.model, m.guide).detach()), 15., rtol=1E-6)
	m.train_svi(2)
	assert_allclose(recorded_losses(m)[:, 0], [15., 15. * 0.81], rtol=1E-5)
	assert_allclose(m.param('mu').detach().numpy(), np.array(TOY_Y) * (1 - 0.9 ** 2), rtol=1E-5)


def test_model_base_init_forwards_loss_arguments(tp):
	"""loss_a reaches the pyro loss constructor: an unknown keyword raises TypeError."""
	torch, pyro = tp
	m = toy_class(torch)('toy')
	m.gen_params({'y': m.tensor(TOY_Y)})
	m.init('SGD', 'Trace_ELBO', optimizer_a=[[{'lr': 0.1}], {}], loss_a=[[], {'num_particles': 2}])
	assert_allclose(float(m.loss[0](m.model, m.guide).detach()), 15., rtol=1E-6)
	with pytest.raises(TypeError):
		m.init('SGD', 'Trace_ELBO', optimizer_a=[[{'lr': 0.1}], {}], loss_a=[[], {'not_an_argument': 1}])


def test_model_base_init_accepts_list_of_losses(tp):
	"""A list of losses gives one value per loss per step, and the update follows their sum."""
	torch, pyro = tp
	m = toy_class(torch)('toy')
	m.gen_params({'y': m.tensor(TOY_Y)})
	m.init('SGD', [m.get_pyro_loss('Trace_ELBO'), m.get_pyro_loss('Trace_ELBO')], optimizer_a=[[{'lr': 0.1}], {}])
	m.train_svi(1)
	assert_allclose(recorded_losses(m)[0], [15., 15.], rtol=1E-6)
	# gradient of 2 * (mu-y)^2/2 is the same as for the single quadratic loss
	assert_allclose(m.param('mu').detach().numpy(), sgd_mu(1), rtol=1E-5)


def test_model_base_init_rejects_invalid_loss(tp):
	"""A loss that is neither a name, a function nor a list raises ValueError."""
	torch, pyro = tp
	m = toy_class(torch)('toy')
	m.gen_params({'y': m.tensor(TOY_Y)})
	with pytest.raises(ValueError):
		m.init('SGD', 5, optimizer_a=[[{'lr': 0.1}], {}])


def test_model_base_plate_subsamples(tp):
	"""A subsample index tensor restricts the plate; pyro scales the factor by 4/2, so the loss is 2*(1+9)."""
	torch, pyro = tp
	m = make_toy(torch, subsamples={'data': torch.tensor([0, 2])})
	m.train_svi(1)
	assert_allclose(recorded_losses(m)[0], [20.], rtol=1E-6)
	# gradient of 2*(mu-y)^2 on entries 0 and 2 only
	assert_allclose(m.param('mu').detach().numpy(), [0.4, 0, 1.2, 0], rtol=1E-5, atol=1E-7)


def test_model_base_save_and_load_param(tp, tmp_path):
	"""save_param and load_param round-trip the parameter store through a file."""
	torch, pyro = tp
	m = make_toy(torch)
	m.train_svi(2)
	path = str(tmp_path / 'params.pt')
	m.save_param(path)
	assert os.path.getsize(path) > 0
	m.clear_param_store()
	assert len(list(m.get_param_store().keys())) == 0
	m.load_param(path)
	assert set(m.get_param_store().keys()) == {'mu', 's', 'empty'}
	assert_allclose(m.param('mu').detach().numpy(), sgd_mu(2), rtol=1E-5)
	assert_allclose(m.param('s').detach().numpy(), [1, 2], rtol=1E-6)


def test_model_base_draw_stats(tp, monkeypatch):
	"""draw_stats plots each subloss and the total from 20% of the steps on, then each parameter's statistics."""
	torch, pyro = tp
	pytest.importorskip('matplotlib')
	import matplotlib.pyplot as plt
	shown = []

	def show(*a, **ka):
		ax = plt.gca()
		shown.append((ax.get_ylabel(), ax.get_title(), [(x.get_label(), np.asarray(x.get_xdata()), np.asarray(x.get_ydata())) for x in ax.get_lines()]))
		plt.close('all')
	monkeypatch.setattr(plt, 'show', show)
	m = make_toy(torch)
	m.train_svi(10, nstep_report=5)
	try:
		m.draw_stats(lossnames=['quad'])
	finally:
		plt.close('all')
	assert len(shown) == 4
	losses = np.array([sgd_loss(x) for x in range(10)])
	assert shown[0][0] == 'Loss quad' and len(shown[0][2]) == 1
	assert_array_equal(shown[0][2][0][1], np.arange(2, 10))
	assert_allclose(shown[0][2][0][2], losses[2:], rtol=1E-5)
	assert shown[1][0] == 'Total loss'
	assert_allclose(shown[1][2][0][2], losses[2:], rtol=1E-5)
	assert shown[2][1] == 'mu [4]'
	assert [x[0] for x in shown[2][2]] == STATS
	assert_array_equal(shown[2][2][0][1], [0, 5, 10])
	assert_allclose(np.array([x[2] for x in shown[2][2]]).T, [summary(sgd_mu(x)) for x in [0, 5, 10]], rtol=1E-5)
	assert shown[3][1].startswith('s [2] ')
	m.draw_stats()
	assert shown[4][0] == 'Loss 0'


def test_clear_param_store(tp):
	"""clear_param_store empties pyro's global parameter store."""
	torch, pyro = tp
	pyro.param('x', torch.tensor(1.))
	assert 'x' in network.model_base.get_param_store()
	network.model_base.clear_param_store()
	assert len(list(pyro.get_param_store().keys())) == 0


###########################################################################
# model_covariance
###########################################################################

def small_counts():
	"""Read counts of 4 genes in 12 cells."""
	return np.random.default_rng(5).poisson(30, size=(4, 12)) + 1


def init_covariance(torch, m, counts, suffix=''):
	"""Generate initial values for m from counts and register them for training."""
	m.gen_params({'G_obs' + suffix: m.tensor(counts, dtype=torch.int)})
	m.init('ClippedAdam', network.Trace_ELBO_site([x + suffix for x in m.varnames]).differentiable_loss, optimizer_a=[[{'lr': 0.01}], {}])
	return m


def traces(m):
	"""A guide trace and the model trace replayed on it, with log probabilities computed."""
	from pyro import poutine
	tg = poutine.trace(m.guide).get_trace()
	tm = poutine.trace(poutine.replay(m.model, trace=tg)).get_trace()
	tg.compute_log_prob()
	tm.compute_log_prob()
	return tg, tm


def pval(pyro, name):
	"""Current value of a pyro parameter as numpy float64."""
	return pyro.param(name).detach().numpy().astype(float)


def test_model_covariance_init(tp):
	"""npc sets the factor count, an (n_cell, 0) covariate matrix adds no covariate effect, tensors default to float32, and a 1-D matrix fails."""
	torch, pyro = tp
	m = network.model_covariance('m', 2, np.zeros((12, 0)), covdc_rate=3.)
	assert m.tensor([1.]).dtype == torch.float32
	assert m.varnames == ['G_0', 'G_obs']
	m.gen_params({'G_obs': m.tensor(small_counts(), dtype=torch.int)})
	assert tuple(m.p['sigma_G_nd'][0].shape) == (4, 2)
	assert 'gamma' not in m.p
	with pytest.raises(AssertionError):
		network.model_covariance('m', 2, np.zeros(12))


@pytest.mark.parametrize('npc', [0, 2])
def test_model_covariance_gen_params_initial_values(tp, npc):
	"""Initial values: log expression shares, a rank-npc PCA of centred log shares and its residual RMS."""
	torch, pyro = tp
	from pyro.distributions import constraints
	counts = small_counts()
	m = network.model_covariance('m', npc, np.zeros((12, 0)))
	m.gen_params({'G_obs': m.tensor(counts, dtype=torch.int)})
	assert m.nrand == 48
	assert set(m.p) == {'mu_G', 'sigma_G_d', 'sigma_G_nd', 'G_0_mean', 'G_0_std'}
	assert (m.p['mu_G'][1].lower_bound, m.p['mu_G'][1].upper_bound) == (-20., 0.)
	assert m.p['sigma_G_d'][1] is constraints.positive and m.p['G_0_std'][1] is constraints.positive
	assert m.p['sigma_G_nd'][1] is None and m.p['G_0_mean'][1] is None
	mu, sd, factor = pca_init(counts, npc)
	assert_allclose(m.p['mu_G'][0].numpy(), mu, atol=1E-5)
	assert_allclose(m.p['sigma_G_d'][0].numpy(), sd, rtol=1E-4, atol=1E-6)
	f = m.p['sigma_G_nd'][0].numpy()
	assert f.shape == (4, npc)
	assert_allclose(f @ f.T, factor @ factor.T, atol=1E-5)
	assert_allclose(np.abs(f), np.abs(factor), atol=1E-4)
	assert_allclose(m.p['G_0_mean'][0].numpy(), np.repeat(mu[None], 12, axis=0), atol=1E-5)
	assert_allclose(m.p['G_0_std'][0].numpy(), np.full((12, 4), 0.5))


def test_model_covariance_gen_params_checks_observations(tp):
	"""Observations must be exactly {'G_obs'+suffix: counts}; the suffix is stripped and the counts set the initial mean."""
	torch, pyro = tp
	counts = small_counts()
	m = network.model_covariance('m', 1, np.zeros((12, 0)), suffix='_x')
	m.gen_params({'G_obs_x': m.tensor(counts, dtype=torch.int)})
	assert_allclose(m.p['mu_G'][0].numpy(), np.log(counts.sum(axis=1) / counts.sum()), atol=1E-5)
	m = network.model_covariance('m', 1, np.zeros((12, 0)))
	with pytest.raises(AssertionError):
		m.gen_params({'G_obs': m.tensor(counts, dtype=torch.int), 'other': m.tensor(counts)})


def test_model_covariance_log_probabilities(tp):
	"""Model and guide log probabilities match scipy: low-rank normal G_0, binomial counts with logits G_0, normal guide."""
	torch, pyro = tp
	from scipy.special import expit
	from scipy.stats import binom, multivariate_normal, norm
	counts = small_counts()
	m = init_covariance(torch, network.model_covariance('m', 2, np.zeros((12, 0))), counts)
	m.get_param_store()['G_0_std'] = torch.full((12, 4), 0.3)
	tg, tm = traces(m)
	g0 = tg.nodes['G_0']['value'].detach().numpy().astype(float)
	assert g0.shape == (12, 4)
	expected = norm.logpdf(g0, pval(pyro, 'G_0_mean'), 0.3).sum()
	assert_allclose(float(tg.nodes['G_0']['log_prob_sum'].detach()), expected, rtol=1E-5)
	f = pval(pyro, 'sigma_G_nd')
	cov = f @ f.T + np.diag(pval(pyro, 'sigma_G_d'))
	expected = multivariate_normal(pval(pyro, 'mu_G'), cov).logpdf(g0).sum()
	assert_allclose(float(tm.nodes['G_0']['log_prob_sum'].detach()), expected, rtol=1E-4)
	assert_array_equal(tm.nodes['G_obs']['value'].numpy(), counts)
	expected = binom.logpmf(counts, counts.sum(axis=0)[None, :], expit(g0.T)).sum()
	assert_allclose(float(tm.nodes['G_obs']['log_prob_sum'].detach()), expected, rtol=1E-4)


@pytest.mark.parametrize('covdc_rate', [1., 4.])
def test_model_covariance_covariates(tp, covdc_rate):
	"""Covariates are orthonormalised, scaled so dc dc^T = covdc_rate*mean|S S^T|*n_cell/n_cov * I, and shift the count logits by dc^T gamma."""
	torch, pyro = tp
	from scipy.special import expit
	from scipy.stats import binom
	counts = small_counts()
	dcs = np.random.default_rng(2).normal(size=(12, 2))
	m = network.model_covariance('m', 2, np.zeros((12, 0)), covdc_rate=covdc_rate)
	# the constructor cannot keep covariate columns (torch max/min return named tuples), so set them here;
	# the orthonormal basis pca_lowrank picks is arbitrary, so the scaled covariates are read from m.a
	m.dcs = m.tensor(dcs)
	m.nl = 2
	init_covariance(torch, m, counts)
	assert m.p['gamma'][1] is None
	assert_array_equal(m.p['gamma'][0].detach().numpy(), np.zeros((2, 4)))
	dc = m.a[2].numpy().astype(float)
	assert dc.shape == (2, 12)
	s = m.p['sigma_G_nd'][0].detach().numpy().astype(float)
	assert_allclose(dc @ dc.T, covdc_rate * np.abs(s @ s.T).mean() * 12 / 2 * np.eye(2), rtol=1E-4, atol=1E-5)
	z = (dcs - dcs.mean(axis=0)) / np.sqrt(((dcs - dcs.mean(axis=0)) ** 2).mean(axis=0))
	res = z - dc.T @ np.linalg.lstsq(dc.T, z, rcond=None)[0]
	assert_allclose(res, 0, atol=1E-4)
	gamma = np.array([[0.1, -0.2, 0.3, 0.], [0.05, 0.1, -0.1, 0.2]])
	m.get_param_store()['gamma'] = torch.tensor(gamma, dtype=torch.float)
	tg, tm = traces(m)
	g0 = tg.nodes['G_0']['value'].detach().numpy().astype(float)
	logits = (g0 + dc.T @ gamma).T
	expected = binom.logpmf(counts, counts.sum(axis=0)[None, :], expit(logits)).sum()
	assert_allclose(float(tm.nodes['G_obs']['log_prob_sum'].detach()), expected, rtol=1E-4)


###########################################################################
# model_ou
###########################################################################

OU_MASK = np.array([[0, 0, 1], [1, 2, 3]])


def make_ou(torch, npc0=1, npc=None, scale_lyapunov=1000.):
	"""model_ou with regulators 0 and 1 over 4 genes (edges 0->1, 0->2, 1->3), initialised on small_counts."""
	m = network.model_ou('m', OU_MASK, (2, 4), npc0, np.zeros((12, 0)), scale_lyapunov=scale_lyapunov, npc=npc)
	return init_covariance(torch, m, small_counts())


def test_model_ou_init(tp):
	"""One edge value per mask column, npc0 noise factors, and npc expression factors defaulting to min(npc0 + n_reg, n_gene - 1)."""
	torch, pyro = tp

	def shapes(npc0, ngene, **ka):
		m = network.model_ou('m', OU_MASK, (2, ngene), npc0, np.zeros((12, 0)), **ka)
		m.gen_params({'G_obs': m.tensor(np.random.default_rng(5).poisson(30, size=(ngene, 12)) + 1, dtype=torch.int)})
		return tuple(tuple(m.p[x][0].shape) for x in ['N_0val', 'sigma0_G_nd', 'sigma_G_nd'])
	assert shapes(1, 4) == ((3,), (4, 1), (4, 3))
	assert shapes(1, 6) == ((3,), (6, 1), (6, 3))
	assert shapes(0, 4, npc=1) == ((3,), (4, 0), (4, 1))
	assert network.model_ou.varnames == ['G_0', 'sigma_G_obs', 'G_obs']


def test_model_ou_init_rejects_invalid_input(tp):
	"""Masks with self-edges, no edges or wrong shape fail; npc above its maximum or below npc0 raises ValueError."""
	torch, pyro = tp
	dc = np.zeros((12, 0))
	for mask in [np.array([[0], [0]]), np.zeros((2, 0), dtype=int), np.array([[0, 1, 1]])]:
		with pytest.raises(AssertionError):
			network.model_ou('m', mask, (2, 4), 1, dc)
	with pytest.raises(ValueError):
		network.model_ou('m', OU_MASK, (2, 4), 1, dc, npc=4)
	with pytest.raises(ValueError):
		network.model_ou('m', OU_MASK, (2, 4), 2, dc, npc=1)


def test_model_ou_gen_params(tp):
	"""Edges start near zero in (-5, 5), reg at 1 in (-10, 10), and sigma0 copies sigma_G_d and npc0 factors."""
	torch, pyro = tp
	m = make_ou(torch, npc0=1)
	assert set(m.p) == {'mu_G', 'sigma_G_d', 'sigma_G_nd', 'G_0_mean', 'G_0_std', 'N_0val', 'reg', 'sigma0_G_d', 'sigma0_G_nd'}
	n0 = m.p['N_0val'][0].numpy()
	assert n0.shape == (3,) and np.abs(n0).max() < 1E-3
	assert (m.p['N_0val'][1].lower_bound, m.p['N_0val'][1].upper_bound) == (-5., 5.)
	assert float(m.p['reg'][0]) == 1. and (m.p['reg'][1].lower_bound, m.p['reg'][1].upper_bound) == (-10, 10)
	assert_array_equal(m.p['sigma0_G_d'][0].detach().numpy(), m.p['sigma_G_d'][0].detach().numpy())
	assert m.p['sigma0_G_d'][1] is m.p['sigma_G_d'][1]
	assert m.p['sigma0_G_nd'][0].shape == (4, 1) and m.p['sigma_G_nd'][0].shape == (4, 3)
	assert_array_equal(m.p['sigma0_G_nd'][0].detach().numpy(), m.p['sigma_G_nd'][0].detach().numpy()[:, :1])


def set_ou_params(torch, m, rng):
	"""Set model_ou's parameters to random values; return them as float64 arrays."""
	vals = {
		'sigma_G_d': rng.uniform(0.5, 1.5, 4), 'sigma_G_nd': rng.normal(0, 0.5, (4, 3)), 'reg': np.array(1.7),
		'N_0val': np.array([0.4, -0.6, 0.3]), 'sigma0_G_d': rng.uniform(0.5, 1.5, 4), 'sigma0_G_nd': rng.normal(0, 0.5, (4, 1))}
	for xi, xj in vals.items():
		m.get_param_store()[xi] = torch.tensor(xj, dtype=torch.float)
	return {x: m.param(x).detach().numpy().astype(float) for x in vals}


def ou_matrices(v, n0val=None):
	"""A = reg*I - N_0 with N_0[target, regulator] = edge value, and X = diag(sigma_G_d) + sigma_G_nd sigma_G_nd^T."""
	n0 = np.zeros((4, 4))
	n0[OU_MASK[1], OU_MASK[0]] = v['N_0val'] if n0val is None else n0val
	a = np.eye(4) * v['reg'] - n0
	x = np.diag(v['sigma_G_d']) + v['sigma_G_nd'] @ v['sigma_G_nd'].T
	return a, x


def test_model_ou_lyapunov_sides(tp):
	"""lyapunov_lhs is A X + X A^T (scipy's solver recovers X from it), lhs0 takes other edges, rhs is diag + low rank."""
	torch, pyro = tp
	from scipy.linalg import solve_continuous_lyapunov
	m = make_ou(torch)
	v = set_ou_params(torch, m, np.random.default_rng(3))
	a, x = ou_matrices(v)
	lhs = m.lyapunov_lhs().detach().numpy().astype(float)
	assert_allclose(lhs, a @ x + x @ a.T, rtol=1E-5, atol=1E-6)
	assert_allclose(solve_continuous_lyapunov(a, lhs), x, rtol=1E-4, atol=1E-5)
	other = np.array([-1., 0.2, 0.5])
	a2, _ = ou_matrices(v, other)
	lhs0 = m.lyapunov_lhs0(torch.tensor(other, dtype=torch.float)).detach().numpy().astype(float)
	assert_allclose(lhs0, a2 @ x + x @ a2.T, rtol=1E-5, atol=1E-6)
	rhs = m.lyapunov_rhs().detach().numpy().astype(float)
	assert_allclose(rhs, np.diag(v['sigma0_G_d']) + v['sigma0_G_nd'] @ v['sigma0_G_nd'].T, rtol=1E-5, atol=1E-6)


def test_model_ou_model_adds_lyapunov_factor(tp):
	"""The model's 'sigma_G_obs' factor is -sum((lhs - rhs)^2) * scale * n_cell / (4 n_gene), beside the covariance sites."""
	torch, pyro = tp
	m = make_ou(torch, scale_lyapunov=30.)
	v = set_ou_params(torch, m, np.random.default_rng(4))
	a, x = ou_matrices(v)
	q = np.diag(v['sigma0_G_d']) + v['sigma0_G_nd'] @ v['sigma0_G_nd'].T
	expected = -(((a @ x + x @ a.T) - q) ** 2).sum() * 30. * 12 / (4 * 4)
	tg, tm = traces(m)
	assert {'G_0', 'G_obs', 'sigma_G_obs'} <= set(tm.nodes)
	assert_allclose(float(tm.nodes['sigma_G_obs']['log_prob_sum'].detach()), expected, rtol=1E-5)
	losses = network.Trace_ELBO_site(m.varnames).differentiable_loss(m.model, m.guide)
	assert len(losses) == 3
	assert_allclose(float(losses[1].detach()), -expected, rtol=1E-5)


###########################################################################
# reconstruct
###########################################################################

def test_reconstruct_output_tables(sim, tmp_path):
	"""Outputs have the documented rows and columns, regulators first, and finite values."""
	out = run_reconstruct(sim, tmp_path, nstep=6, nstep_report=2)
	assert list(out['weight'].index) == REGS and list(out['weight'].columns) == ORDER
	assert list(out['meanvar'].index) == ORDER and list(out['meanvar'].columns) == ['mean', 'var']
	assert list(out['covfactor'].index) == ORDER
	# min(n_reg + npc, n_gene - 1) factors
	assert list(out['covfactor'].columns) == ['factor1', 'factor2']
	assert list(out['loss'].columns) == ['G_0', 'sigma_G_obs', 'G_obs']
	assert list(out['loss'].index) == list(range(6))
	st = out['stats']
	assert list(st.columns) == ['variable', 'stat', 'step', 'value']
	assert set(st['variable']) == {'mu_G', 'sigma_G_d', 'sigma_G_nd', 'G_0_mean', 'G_0_std', 'N_0val', 'reg', 'sigma0_G_d'}
	assert set(st['stat']) == set(STATS)
	assert sorted(set(st['step'])) == [0, 2, 4, 6] and len(st) == 8 * 4 * 4
	for xi in ['weight', 'meanvar', 'covfactor', 'loss']:
		assert np.isfinite(out[xi].values).all()
	assert np.isfinite(st['value'].values).all()
	assert (out['meanvar']['var'] > 0).all()
	assert_allclose(np.exp(out['meanvar']['mean']).sum(), 1, rtol=1E-6)


def test_reconstruct_respects_mask(sim, tmp_path):
	"""Edges outside the mask, including self-edges, are exactly zero; allowed edges are nonzero."""
	w = run_reconstruct(sim, tmp_path, nstep=20)['weight'].values
	assert (w[~sim['allowed']] == 0).all()
	assert (w[sim['allowed']] != 0).all()


def test_reconstruct_is_deterministic_for_a_seed(sim, tmp_path):
	"""The same rseed gives identical outputs; another rseed gives other weights."""
	a = run_reconstruct(sim, tmp_path, 'a', nstep=10, rseed=3)
	b = run_reconstruct(sim, tmp_path, 'b', nstep=10, rseed=3)
	c = run_reconstruct(sim, tmp_path, 'c', nstep=10, rseed=4)
	assert_outputs_equal(a, b)
	assert not np.allclose(a['weight'].values[sim['allowed']], c['weight'].values[sim['allowed']])


def test_reconstruct_lr0_outputs_initial_values(sim, sim_lr0):
	"""With lr=0: mean is the log expression share, factors and var are a rank-2 PCA of log shares and its residual RMS, reg stays 1."""
	out = sim_lr0
	mu, sd, factor = pca_init(sim['expr'].loc[ORDER].values, 2)
	assert_allclose(out['meanvar']['mean'].values, mu, atol=2E-5)
	assert_allclose(out['meanvar']['var'].values, sd, rtol=1E-4, atol=1E-6)
	f = out['covfactor'].values
	assert_allclose(f @ f.T, factor @ factor.T, atol=1E-5)
	assert_allclose(np.abs(f), np.abs(factor), atol=1E-4)
	st = out['stats']
	assert_allclose(st[st.variable == 'reg']['value'].values, 1, atol=1E-5)
	w = out['weight'].values
	assert np.abs(w).max() < 1E-3
	assert (out['loss']['sigma_G_obs'] == out['loss']['sigma_G_obs'][0]).all()


@pytest.mark.parametrize('npc', [0, 1])
def test_reconstruct_lyapunov_loss_at_initial_values(sim, tmp_path, npc):
	"""The sigma_G_obs loss equals the Lyapunov loss recomputed from lr=0 outputs, also in the first step of training."""
	base = run_reconstruct(sim, tmp_path, 'a', nstep=2, lr=0, npc=npc, scale_lyapunov=1000)
	expected = lyapunov_loss(base['weight'].values, base['meanvar']['var'].values, base['covfactor'].values, npc, 1000, 40)
	assert_allclose(base['loss']['sigma_G_obs'].values, expected, rtol=1E-4)
	trained = run_reconstruct(sim, tmp_path, 'b', nstep=2, npc=npc, scale_lyapunov=1000)
	assert_allclose(trained['loss']['sigma_G_obs'].values[0], expected, rtol=1E-4)


def test_reconstruct_scale_lyapunov_scales_lyapunov_loss(sim, tmp_path):
	"""At the shared initial values a 100x scale_lyapunov gives a 100x sigma_G_obs loss and the same other losses."""
	a = run_reconstruct(sim, tmp_path, 'a', nstep=1, scale_lyapunov=10)
	b = run_reconstruct(sim, tmp_path, 'b', nstep=1, scale_lyapunov=1000)
	assert_allclose(b['loss']['sigma_G_obs'].values, 100 * a['loss']['sigma_G_obs'].values, rtol=1E-5)
	pd.testing.assert_frame_equal(a['loss'][['G_0', 'G_obs']], b['loss'][['G_0', 'G_obs']], check_exact=True)


@pytest.mark.parametrize('lr,lrd', [(0.01, 0.5), (0.002, 1.)])
def test_reconstruct_first_step_moves_factors_by_lr_times_lrd(sim, sim_lr0, tmp_path, lr, lrd):
	"""One ClippedAdam step moves each covariance factor (unconstrained) by lr*lrd: the first Adam step is lr*lrd*sign(gradient)."""
	out = run_reconstruct(sim, tmp_path, nstep=1, lr=lr, lrd=lrd)
	assert_allclose(np.abs(out['covfactor'].values - sim_lr0['covfactor'].values), lr * lrd, rtol=1E-3)


def test_reconstruct_lrd_zero_freezes_parameters(sim, sim_lr0, tmp_path):
	"""ClippedAdam multiplies lr by lrd before each step, so lrd=0 gives the same outputs as lr=0."""
	assert_outputs_equal(run_reconstruct(sim, tmp_path, nstep=3, lrd=0), sim_lr0)


def test_reconstruct_stats_match_final_parameters(sim, tmp_path):
	"""Stats are recorded at multiples of nstep_report; at the last step they summarise the written outputs."""
	st = run_reconstruct(sim, tmp_path, 'a', nstep=7, nstep_report=3)['stats']
	assert sorted(set(st['step'])) == [0, 3, 6]
	out = run_reconstruct(sim, tmp_path, 'b', nstep=6, nstep_report=3)
	st = out['stats']

	def stat(name, step):
		return st[(st.variable == name) & (st.step == step)].set_index('stat').loc[STATS, 'value'].values
	reg = stat('reg', 6)
	assert np.ptp(reg) == 0 and reg[0] != 1
	n0 = out['weight'].values[sim['allowed']] * reg[0]
	assert_allclose(stat('N_0val', 6), summary(n0), rtol=1E-5, atol=1E-7)
	assert_allclose(stat('sigma_G_d', 6), summary(out['meanvar']['var'].values), rtol=1E-6)
	assert_allclose(stat('sigma_G_nd', 6), summary(out['covfactor'].values), rtol=1E-6, atol=1E-8)
	mu = stat('mu_G', 6)
	mean = out['meanvar']['mean'].values
	assert_allclose(mu[1] - mu[0], mean.max() - mean.min(), rtol=1E-5)
	assert_allclose(stat('reg', 0), 1, atol=1E-5)
	assert_allclose(stat('G_0_std', 0), 0.5, rtol=1E-6)
	assert_allclose(stat('sigma0_G_d', 0), stat('sigma_G_d', 0), rtol=1E-6)


def test_reconstruct_recovers_edge_signs(sim, tmp_path):
	"""On data where 'up' follows 'tf' and 'down' mirrors it, 50 steps give a positive and a negative edge and the loss falls."""
	out = run_reconstruct(sim, tmp_path, nstep=50)
	w = out['weight']
	assert w.loc['tf', 'up'] > 0.5
	assert w.loc['tf', 'down'] < -0.5
	assert abs(w.loc['tf', 'ind']) < 0.3
	total = out['loss'].sum(axis=1).values
	assert total[-5:].max() < 0.05 * total[0]


@pytest.mark.parametrize('npc,nfactor', [(0, 2), (1, 3), (2, 4), (4, 4)])
def test_reconstruct_npc_sets_factor_count(sim, tmp_path, npc, nfactor):
	"""Factors number min(n_reg + npc, n_gene - 1); sigma0_G_nd holds the first npc initial factors."""
	out = run_reconstruct(sim, tmp_path, nstep=2, lr=0, npc=npc)
	f = out['covfactor'].values
	assert f.shape == (5, nfactor)
	st = out['stats']
	assert ('sigma0_G_nd' in set(st['variable'])) == (npc > 0)
	if npc > 0:
		got = st[(st.variable == 'sigma0_G_nd') & (st.step == 0)].set_index('stat').loc[STATS, 'value'].values
		assert_allclose(got, summary(f[:, :npc]), rtol=1E-5, atol=1E-8)


def test_reconstruct_varstd_writes_z_scores(sim, tmp_path):
	"""fo_weightz holds varmean/varstd on allowed edges: N_0val/N_0val is 1, and N_0val/reg is the weight."""
	z1 = str(tmp_path / 'z1.tsv')
	run_reconstruct(sim, tmp_path, 'a', nstep=5, varstd='N_0val', fo_weightz=z1)
	z = read_tsv(z1)
	assert list(z.index) == REGS and list(z.columns) == ORDER
	assert_allclose(z.values[sim['allowed']], 1, rtol=1E-6)
	assert (z.values[~sim['allowed']] == 0).all()
	z2 = str(tmp_path / 'z2.tsv')
	out = run_reconstruct(sim, tmp_path, 'b', nstep=5, varstd='reg', fo_weightz=z2)
	assert_allclose(read_tsv(z2).values, out['weight'].values, rtol=1E-6)


def test_reconstruct_varmean_selects_parameter(sim, tmp_path):
	"""varmean='reg' puts the scalar reg on every allowed edge, which is 1 after division by reg."""
	w = run_reconstruct(sim, tmp_path, nstep=5, varmean='reg')['weight'].values
	assert_allclose(w[sim['allowed']], 1, rtol=1E-6)
	assert (w[~sim['allowed']] == 0).all()


def test_reconstruct_weightz_requires_varstd(sim, tmp_path):
	"""fo_weightz without varstd raises ValueError."""
	with pytest.raises(ValueError):
		run_reconstruct(sim, tmp_path, nstep=2, fo_weightz=str(tmp_path / 'z.tsv'))


def test_reconstruct_fi_cov_without_columns_changes_nothing(sim, tmp_path):
	"""A covariate table listing the cells (in another order) but no covariates gives the same outputs as none."""
	path = str(tmp_path / 'cov.tsv')
	pd.DataFrame(index=pd.Index(sim['expr'].columns[::-1], name='cell')).to_csv(path, sep='\t')
	a = run_reconstruct(sim, tmp_path, 'a', nstep=5, fi_cov=path)
	b = run_reconstruct(sim, tmp_path, 'b', nstep=5)
	assert_outputs_equal(a, b)


def test_reconstruct_fi_cov_errors(sim, tmp_path):
	"""Covariate tables missing cells or rows raise ValueError; covariates count against the npc budget."""
	cells = sim['expr'].columns
	path = str(tmp_path / 'cov.tsv')
	pd.DataFrame({'a': np.arange(len(cells) - 1.)}, index=cells[1:]).to_csv(path, sep='\t')
	with pytest.raises(ValueError):
		run_reconstruct(sim, tmp_path, nstep=2, fi_cov=path)
	with open(path, 'w') as f:
		f.write('cell\ta\n')
	with pytest.raises(ValueError):
		run_reconstruct(sim, tmp_path, nstep=2, fi_cov=path)
	# 5 genes - 3 covariates leaves 2 degrees of freedom, so npc=2 is too many
	pd.DataFrame(np.random.default_rng(0).normal(size=(len(cells), 3)), index=cells, columns=['a', 'b', 'c']).to_csv(path, sep='\t')
	with pytest.raises(RuntimeError):
		run_reconstruct(sim, tmp_path, nstep=2, fi_cov=path, npc=2)


def test_reconstruct_matches_genes_by_name(sim, tmp_path):
	"""Expression rows are matched to the mask by name: shuffled rows and an extra gene give the same outputs."""
	expr = sim['expr'].iloc[[3, 0, 4, 2, 1]].copy()
	expr.loc['extra'] = np.arange(expr.shape[1]) + 5
	path = str(tmp_path / 'expr.tsv')
	expr.to_csv(path, sep='\t')
	a = run_reconstruct(sim, tmp_path, 'a', exp=path, nstep=5)
	b = run_reconstruct(sim, tmp_path, 'b', nstep=5)
	assert_outputs_equal(a, b)


def test_reconstruct_rejects_invalid_input(sim, tmp_path):
	"""Unknown models raise ValueError; a mask without regulators or npc beyond the genes raises RuntimeError."""
	with pytest.raises(ValueError):
		run_reconstruct(sim, tmp_path, nstep=2, model='cov')
	with pytest.raises(RuntimeError):
		run_reconstruct(sim, tmp_path, nstep=2, npc=5)
	path = str(tmp_path / 'mask.tsv')
	with open(path, 'w') as f:
		f.write('\t' + '\t'.join(GENES) + '\n')
	with pytest.raises(RuntimeError):
		network.reconstruct(sim['exp'], path, *output_paths(tmp_path), nstep=2, nth=1)


def test_reconstruct_dtype_float32_is_default(sim, tmp_path):
	"""dtype='float32' gives the same outputs as the default 'float'."""
	assert_outputs_equal(run_reconstruct(sim, tmp_path, 'a', nstep=5, dtype='float32'), run_reconstruct(sim, tmp_path, 'b', nstep=5))


def test_reconstruct_thread_count_does_not_change_results(sim, tmp_path):
	"""nth=2 gives the same outputs as nth=1 up to float rounding accumulated over 10 steps."""
	a = run_reconstruct(sim, tmp_path, 'a', nstep=10, nth=2)
	b = run_reconstruct(sim, tmp_path, 'b', nstep=10, nth=1)
	for xi in ['weight', 'meanvar', 'covfactor', 'loss']:
		assert_allclose(a[xi].values, b[xi].values, rtol=1E-4, atol=1E-6)


def test_cli_reconstruct_all_options(sim, tmp_path, run_dictys):
	"""`dictys network reconstruct` with every option writes the same files as the function call."""
	cov = str(tmp_path / 'cov.tsv')
	pd.DataFrame(index=pd.Index(sim['expr'].columns, name='cell')).to_csv(cov, sep='\t')
	ka = dict(lr=0.02, lrd=0.99, nstep=6, npc=1, fi_cov=cov, model='ou', nstep_report=3, rseed=5, device='cpu', dtype='float',
		loss='Trace_ELBO_site', nth=1, varmean='N_0val', varstd='reg', scale_lyapunov=2000.)
	cli = output_paths(tmp_path, 'cli_')
	zcli = str(tmp_path / 'cli_z.tsv')
	opts = []
	for xi, xj in ka.items():
		opts += ['--' + xi, xj]
	run_dictys('network', 'reconstruct', sim['exp'], sim['mask'], *cli, *opts, '--fo_weightz', zcli)
	fun = output_paths(tmp_path, 'fun_')
	zfun = str(tmp_path / 'fun_z.tsv')
	network.reconstruct(sim['exp'], sim['mask'], *fun, fo_weightz=zfun, **ka)
	for xi, xj in zip(cli + [zcli], fun + [zfun]):
		with open(xi) as f1, open(xj) as f2:
			assert f1.read() == f2.read()
	out = dict(zip(OUTPUTS, [read_tsv(x) for x in cli]))
	assert out['loss'].shape == (6, 3) and out['covfactor'].shape == (5, 3)
	assert sorted(set(out['stats']['step'])) == [0, 3, 6]
	assert_allclose(read_tsv(zcli).values, out['weight'].values, rtol=1E-6)


###########################################################################
# indirect
###########################################################################

IREGS = ['r1', 'r2', 'r3']
IGENES = ['r1', 'r2', 'r3', 't1', 't2']
ICOLS = ['t1', 'r2', 'r1', 't2', 'r3']


def zero_sum_network(scale=1.):
	"""
	Weights of 3 regulators on 5 genes with feedback among the regulators, whose rows and columns all sum to
	zero. With equal gene means the saturation term is then the same constant 1/5 for every entry, whichever
	axis it is taken along.
	"""
	b = np.array([[0, .9, -.5], [.7, 0, .5], [-.7, -.9, 0]])
	t = np.array([.4, -.6, .2])
	return np.c_[b, t, -b.sum(axis=1) - t] * scale


def indirect_reference(w, mean, norm, eigmin=0.2, eigmax=1.8, multiplier=1.1):
	"""
	Steady-state effects dx = inv(I - W^T/reg)^T, minus the expression-weighted change sum_k e_k dx[i,k]/sum(e)
	when norm & 2, with reg multiplied until the regulators' diagonal lies in [eigmin, eigmax]; then times reg,
	rows divided by their diagonal when norm & 1, and the diagonal set to zero. Returns (effects, reg).
	"""
	nreg, ng = w.shape
	net0 = np.vstack([w, np.zeros((ng - nreg, ng))])
	e = np.exp(mean)
	lo = -np.inf if eigmin is None else eigmin
	hi = np.inf if eigmax is None else eigmax
	reg = 1.
	for _ in range(1000):
		dx = np.linalg.inv(np.eye(ng) - net0.T / reg).T
		if norm & 2:
			dx = dx - (dx @ e)[:, None] / e.sum()
		d = np.diag(dx)[:nreg]
		if d.min() >= lo and d.max() <= hi:
			break
		reg *= multiplier
	else:
		raise AssertionError('bounds not reached')
	ans = dx * reg
	if norm & 1:
		ans = ans / np.diag(ans)[:, None]
	np.fill_diagonal(ans, 0)
	return ans[:nreg], reg


def write_indirect_inputs(folder, w, mean, factors=None):
	"""Weight (columns shuffled), covariance factor and mean/var files for IGENES; returns their paths."""
	folder = str(folder)
	if factors is None:
		factors = np.arange(10.).reshape(5, 2) / 10
	paths = [os.path.join(folder, x) for x in ['weight.tsv', 'covfactor.tsv', 'meanvar.tsv']]
	pd.DataFrame(w, index=IREGS, columns=IGENES)[ICOLS].to_csv(paths[0], sep='\t')
	pd.DataFrame(factors, index=IGENES, columns=['factor1', 'factor2']).iloc[::-1].to_csv(paths[1], sep='\t')
	pd.DataFrame({'mean': mean, 'var': np.ones(5)}, index=IGENES).iloc[[2, 4, 0, 1, 3]].to_csv(paths[2], sep='\t')
	return paths


EQUAL_MEANS = np.log(np.full(5, 0.2))


def test_indirect_reference_setup():
	"""The zero-sum network makes the regularization loop iterate: reg is 1.1^3 for the default bounds."""
	_, reg = indirect_reference(zero_sum_network(), EQUAL_MEANS, 3)
	assert_allclose(reg, 1.1 ** 3)
	_, reg = indirect_reference(zero_sum_network(), EQUAL_MEANS, 3, 0.6, 1.1, 1.5)
	assert_allclose(reg, 1.5 ** 2)


@pytest.mark.parametrize('norm', [2, 3])
@pytest.mark.parametrize('bounds', [(0.2, 1.8, 1.1), (0.6, 1.1, 1.5), (None, None, 1.1), (0.5, None, 1.2), (None, 1.2, 1.3)])
def test_indirect_matches_reference(threads, tmp_path, norm, bounds):
	"""Steady-state effects match a numpy re-derivation for each bound and multiplier, with outputs in regulator-first order."""
	fw, fc, fm = write_indirect_inputs(tmp_path, zero_sum_network(), EQUAL_MEANS)
	fo = str(tmp_path / 'out.tsv')
	network.indirect(fw, fc, fo, norm=norm, fi_meanvar=fm, eigmin=bounds[0], eigmax=bounds[1], multiplier=bounds[2], nth=1)
	out = read_tsv(fo)
	assert list(out.index) == IREGS and list(out.columns) == IGENES
	expected, _ = indirect_reference(zero_sum_network(), EQUAL_MEANS, norm, *bounds)
	assert_allclose(out.values, expected, rtol=1E-9, atol=1E-12)
	assert (np.diag(out.values[:, :3]) == 0).all()


def test_indirect_unequal_means_matched_by_name(threads, tmp_path):
	"""Gene means are matched by name; reference for a network whose rows and columns sum to the same value per gene."""
	# Each gene's row sum equals its column sum s_j and e_j = 1 - s_j, so the saturation term is constant
	w = np.zeros((3, 5))
	w[0, 1] = w[1, 2] = w[2, 0] = 0.3
	w[0, 1] += 0.2
	w[1, 0] = 0.2
	mean = np.log(1 - np.r_[w.sum(axis=1), 0, 0])
	assert_allclose(np.r_[w.sum(axis=1), 0, 0], w.sum(axis=0))
	fw, fc, fm = write_indirect_inputs(tmp_path, w, mean)
	fo = str(tmp_path / 'out.tsv')
	network.indirect(fw, fc, fo, norm=2, fi_meanvar=fm, eigmin=None, eigmax=None)
	expected, _ = indirect_reference(w, mean, 2, None, None)
	assert_allclose(read_tsv(fo).values, expected, rtol=1E-9, atol=1E-12)


def test_indirect_bounds_change_regularization(threads, tmp_path):
	"""Tighter bounds or a larger multiplier change the output, as the reference's reg does."""
	fw, fc, fm = write_indirect_inputs(tmp_path, zero_sum_network(), EQUAL_MEANS)
	outs = []
	for xi, bounds in enumerate([(None, None, 1.1), (0.2, 1.8, 1.1), (0.2, 1.8, 2.)]):
		fo = str(tmp_path / f'out{xi}.tsv')
		network.indirect(fw, fc, fo, norm=2, fi_meanvar=fm, eigmin=bounds[0], eigmax=bounds[1], multiplier=bounds[2])
		outs.append(read_tsv(fo).values)
	assert not np.allclose(outs[0], outs[1]) and not np.allclose(outs[1], outs[2])
	# with no bounds reg stays 1: dx - 1/5 off the diagonal
	dx = np.linalg.inv(np.eye(5) - np.vstack([zero_sum_network(), np.zeros((2, 5))]).T).T - 0.2
	np.fill_diagonal(dx, 0)
	assert_allclose(outs[0], dx[:3], rtol=1E-9, atol=1E-12)


def regulator_diagonal(reg):
	"""Diagonal of inv(I - W^T/reg)^T - 1/5 for the zero-sum network's regulators, as indirect's bound check sees it."""
	dx = np.linalg.inv(np.eye(5) - np.vstack([zero_sum_network(), np.zeros((2, 5))]).T / reg).T - 0.2
	return np.diag(dx)[:3]


@pytest.mark.parametrize('side', ['max', 'min'])
def test_indirect_stops_at_first_reg_within_bounds(threads, tmp_path, side):
	"""With a bound 1e-9 beyond the diagonal at reg=1.1^k, indirect stops at exactly that reg (k=3 for eigmax, k=4 for eigmin)."""
	if side == 'max':
		k = 3
		bounds = (None, regulator_diagonal(1.1 ** k).max() + 1E-9, 1.1)
		assert regulator_diagonal(1.1 ** (k - 1)).max() > bounds[1] + 0.1
	else:
		k = 4
		bounds = (regulator_diagonal(1.1 ** k).min() - 1E-9, None, 1.1)
		assert regulator_diagonal(1.1 ** (k - 1)).min() < bounds[0] - 0.01
	expected, reg = indirect_reference(zero_sum_network(), EQUAL_MEANS, 2, *bounds)
	assert_allclose(reg, 1.1 ** k)
	fw, fc, fm = write_indirect_inputs(tmp_path, zero_sum_network(), EQUAL_MEANS)
	fo = str(tmp_path / 'out.tsv')
	network.indirect(fw, fc, fo, norm=2, fi_meanvar=fm, eigmin=bounds[0], eigmax=bounds[1], multiplier=bounds[2])
	assert_allclose(read_tsv(fo).values, expected, rtol=1E-9, atol=1E-12)


def test_indirect_ignores_covfactor_values(threads, tmp_path):
	"""Covariance factor values do not enter the result; their gene names must match the weight columns."""
	fw, fc, fm = write_indirect_inputs(tmp_path, zero_sum_network(), EQUAL_MEANS)
	network.indirect(fw, fc, str(tmp_path / 'a.tsv'), fi_meanvar=fm)
	pd.DataFrame(np.random.default_rng(0).normal(size=(5, 3)), index=IGENES).to_csv(fc, sep='\t')
	network.indirect(fw, fc, str(tmp_path / 'b.tsv'), fi_meanvar=fm)
	pd.testing.assert_frame_equal(read_tsv(str(tmp_path / 'a.tsv')), read_tsv(str(tmp_path / 'b.tsv')), check_exact=True)
	pd.DataFrame(np.zeros((4, 2)), index=IGENES[:4]).to_csv(fc, sep='\t')
	with pytest.raises(AssertionError):
		network.indirect(fw, fc, str(tmp_path / 'c.tsv'), fi_meanvar=fm)


def test_indirect_rejects_invalid_arguments(threads, tmp_path):
	"""Invalid norm raises ValueError; multiplier <= 1, bounds not around 1, or a regulator absent from the columns fail."""
	fw, fc, fm = write_indirect_inputs(tmp_path, zero_sum_network(), EQUAL_MEANS)
	fo = str(tmp_path / 'out.tsv')
	with pytest.raises(ValueError):
		network.indirect(fw, fc, fo, norm=4, fi_meanvar=fm)
	for ka in [{'multiplier': 1.}, {'eigmin': 1.}, {'eigmax': 0.9}]:
		with pytest.raises(AssertionError):
			network.indirect(fw, fc, fo, fi_meanvar=fm, **ka)
	pd.DataFrame(zero_sum_network()[:, :4], index=IREGS, columns=['t1', 'r1', 'r2', 'x']).to_csv(fw, sep='\t')
	with pytest.raises(AssertionError):
		network.indirect(fw, fc, fo, fi_meanvar=fm)
	assert not os.path.exists(fo)


def test_indirect_thread_count_does_not_change_results(threads, tmp_path):
	"""nth=2 gives the same output as nth=1."""
	fw, fc, fm = write_indirect_inputs(tmp_path, zero_sum_network(2.), EQUAL_MEANS)
	network.indirect(fw, fc, str(tmp_path / 'a.tsv'), fi_meanvar=fm, nth=1)
	network.indirect(fw, fc, str(tmp_path / 'b.tsv'), fi_meanvar=fm, nth=2)
	assert_allclose(read_tsv(str(tmp_path / 'a.tsv')).values, read_tsv(str(tmp_path / 'b.tsv')).values, rtol=1E-12)


def test_cli_indirect_all_options(threads, tmp_path, run_dictys):
	"""`dictys network indirect` with every option matches the reference."""
	fw, fc, fm = write_indirect_inputs(tmp_path, zero_sum_network(), EQUAL_MEANS)
	fo = str(tmp_path / 'out.tsv')
	run_dictys('network', 'indirect', fw, fc, fo, '--norm', 3, '--fi_meanvar', fm, '--eigmin', 0.6, '--eigmax', 1.1, '--multiplier', 1.5, '--nth', 1)
	expected, _ = indirect_reference(zero_sum_network(), EQUAL_MEANS, 3, 0.6, 1.1, 1.5)
	assert_allclose(read_tsv(fo).values, expected, rtol=1E-9, atol=1E-12)


###########################################################################
# normalize
###########################################################################

def write_normalize_inputs(folder):
	"""Weights of r1, r2 on 4 genes (columns not regulator-first), shuffled mean/var and factor tables."""
	rng = np.random.default_rng(9)
	genes = ['r1', 'r2', 't1', 't2']
	w = pd.DataFrame(rng.normal(size=(2, 4)), index=['r1', 'r2'], columns=genes)[['t1', 'r2', 't2', 'r1']]
	mv = pd.DataFrame({'mean': rng.normal(size=4), 'var': rng.uniform(0.5, 2, 4)}, index=genes)
	f = pd.DataFrame(rng.normal(size=(4, 2)), index=genes, columns=['factor1', 'factor2'])
	paths = [os.path.join(str(folder), x) for x in ['weight.tsv', 'meanvar.tsv', 'covfactor.tsv']]
	w.to_csv(paths[0], sep='\t')
	mv.iloc[[3, 1, 0, 2]].to_csv(paths[1], sep='\t')
	f.iloc[[2, 0, 3, 1]].to_csv(paths[2], sep='\t')
	return paths, w, mv, f


@pytest.mark.parametrize('norm', [1, 2, 3])
def test_normalize_matches_formula(threads, tmp_path, norm):
	"""Weights times the regulator's std (norm&1) and divided by the target's std (norm&2), std = sqrt(var + sum of squared factors)."""
	paths, w, mv, f = write_normalize_inputs(tmp_path)
	fo = str(tmp_path / 'out.tsv')
	network.normalize(*paths, fo, norm=norm, nth=1)
	out = read_tsv(fo)
	assert list(out.index) == ['r1', 'r2'] and list(out.columns) == ['r1', 'r2', 't1', 't2']
	std = np.sqrt(mv['var'] + (f ** 2).sum(axis=1))
	expected = w.copy()
	if norm & 1:
		expected = expected.mul(std[expected.index], axis=0)
	if norm & 2:
		expected = expected / std[expected.columns]
	assert_allclose(out.values, expected[out.columns].values, rtol=1E-12)


def test_normalize_rejects_invalid_input(threads, tmp_path):
	"""norm outside 1-3 raises ValueError; gene sets that disagree between files fail."""
	paths, w, mv, f = write_normalize_inputs(tmp_path)
	fo = str(tmp_path / 'out.tsv')
	for norm in [0, 4]:
		with pytest.raises(ValueError):
			network.normalize(*paths, fo, norm=norm)
	mv.iloc[:3].to_csv(paths[1], sep='\t')
	with pytest.raises(AssertionError):
		network.normalize(*paths, fo)
	assert not os.path.exists(fo)


def test_normalize_thread_count_does_not_change_results(threads, tmp_path):
	"""nth=2 gives the same output as nth=1."""
	paths = write_normalize_inputs(tmp_path)[0]
	network.normalize(*paths, str(tmp_path / 'a.tsv'), nth=1)
	network.normalize(*paths, str(tmp_path / 'b.tsv'), nth=2)
	assert_allclose(read_tsv(str(tmp_path / 'a.tsv')).values, read_tsv(str(tmp_path / 'b.tsv')).values, rtol=1E-12)


def test_cli_normalize_all_options(threads, tmp_path, run_dictys):
	"""`dictys network normalize` with every option applies the regulator-std scaling only for --norm 1."""
	paths, w, mv, f = write_normalize_inputs(tmp_path)
	fo = str(tmp_path / 'out.tsv')
	run_dictys('network', 'normalize', *paths, fo, '--norm', 1, '--nth', 1)
	std = np.sqrt(mv['var'] + (f ** 2).sum(axis=1))
	expected = w.mul(std[w.index], axis=0)
	out = read_tsv(fo)
	assert_allclose(out.values, expected[out.columns].values, rtol=1E-12)


###########################################################################
# tofile
###########################################################################

def h5_tree(path):
	"""All datasets of an h5 file as {path: numpy array}, with fixed and variable-length byte strings decoded."""
	import h5py
	ans = {}

	def visit(name, obj):
		if isinstance(obj, h5py.Dataset):
			v = np.array(obj)
			if v.dtype.char == 'S':
				v = v.astype(str)
			elif v.dtype == object:
				v = np.array([x.decode() if isinstance(x, bytes) else x for x in v.ravel()]).reshape(v.shape)
			ans[name] = v
	with h5py.File(path, 'r') as f:
		f.visititems(visit)
	return ans


TF_NETS = {
	's1': pd.DataFrame([[0, 1.5, -2., 0], [3., 0, 0, 4.]], index=['a', 'b'], columns=['a', 'b', 'c', 'd']),
	's2': pd.DataFrame([[0, 5., 6., 7.]], index=['a'], columns=['a', 'b', 'c', 'e'])}
TF_INETS = {'s1': pd.DataFrame([[0, -1., -2.5, .5], [-3., 0, 1., -4.]], index=['a', 'b'], columns=['a', 'b', 'c', 'd'])}
TF_GENES = ['a', 'b', 'c', 'd', 'e']


def build_tofile_layout(root, subsets, nettype='n'):
	"""
	Data folder with expression.tsv.gz (genes a-e and z, cells c0-c7) and working folder with subsets.txt,
	<subset>/names_rna.txt and the networks in TF_NETS (direct) and TF_INETS (steady-state).
	"""
	data = os.path.join(str(root), 'data')
	work = os.path.join(str(root), 'work')
	os.makedirs(data)
	os.makedirs(work)
	counts = np.random.default_rng(11).integers(0, 1000, size=(6, 8))
	counts[0, 0] = 999
	expr = pd.DataFrame(counts, index=TF_GENES + ['z'], columns=['c' + str(x) for x in range(8)])
	expr.to_csv(os.path.join(data, 'expression.tsv.gz'), sep='\t')
	write_lines(os.path.join(work, 'subsets.txt'), list(subsets))
	for xi, xj in subsets.items():
		os.makedirs(os.path.join(work, xi))
		write_lines(os.path.join(work, xi, 'names_rna.txt'), xj)
		if xi in TF_NETS:
			TF_NETS[xi].to_csv(os.path.join(work, xi, f'net_{nettype}weight.tsv.gz'), sep='\t')
		if xi in TF_INETS:
			TF_INETS[xi].to_csv(os.path.join(work, xi, f'net_i{nettype}weight.tsv.gz'), sep='\t')
	return data, work, os.path.join(work, 'subsets.txt'), expr


def expected_edges(nets, snames):
	"""Edge weights and masks (regulator a/b x gene a-e x subset) placed by name."""
	w = np.zeros((2, 5, len(snames)))
	m = np.zeros((2, 5, len(snames)), dtype=bool)
	for xk, s in enumerate(snames):
		if s not in nets:
			continue
		d = nets[s]
		for xi in d.index:
			for xj in d.columns:
				w[['a', 'b'].index(xi), TF_GENES.index(xj), xk] = d.loc[xi, xj]
				m[['a', 'b'].index(xi), TF_GENES.index(xj), xk] = True
	return w, m


STATIC_SUBSETS = {'s1': ['c0', 'c1', 'c2'], 's2': ['c3', 'c4'], 's3': ['c5', 'c6']}


def test_tofile_static_networks(tmp_path):
	"""The h5 file holds names, memberships, read counts, CPM and edge weights placed by name; a subset without networks is dropped."""
	data, work, fsub, expr = build_tofile_layout(tmp_path, STATIC_SUBSETS)
	fo = str(tmp_path / 'net.h5')
	network.tofile(data, work, fsub, fo)
	h = h5_tree(fo)
	cells = ['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6']
	assert_array_equal(h['_version_'], [0, 1, 0])
	assert_array_equal(h['cname'], cells)
	assert_array_equal(h['sname'], ['s1', 's2'])
	assert_array_equal(h['nname'], TF_GENES)
	assert_array_equal(h['nids1'], [0, 1])
	assert_array_equal(h['nids2'], [0, 1, 2, 3, 4])
	member = np.array([[x in STATIC_SUBSETS[y] for x in cells] for y in ['s1', 's2']], dtype=float)
	assert_array_equal(h['prop/sc/w'], member)
	rc = expr.loc[TF_GENES, cells].values
	assert_array_equal(h['prop/nc/readcount'], rc)
	assert h['prop/nc/readcount'].dtype == np.int16
	pseudo = rc @ member.T
	assert_allclose(h['prop/ns/cpm'], 1E6 * pseudo / pseudo.sum(axis=0), rtol=1E-12)
	w, m = expected_edges(TF_NETS, ['s1', 's2'])
	assert_array_equal(h['prop/es/w_n'], w)
	assert_array_equal(h['prop/es/mask_n'], m)
	assert_array_equal(h['prop/es/w'], w)
	assert_array_equal(h['prop/es/mask'], m)
	w, m = expected_edges(TF_INETS, ['s1', 's2'])
	assert_array_equal(h['prop/es/w_in'], w)
	assert_array_equal(h['prop/es/mask_in'], m)
	assert not any(x.startswith('prop/c/') or x.startswith('traj') or x.startswith('point') for x in h)


def test_tofile_reads_network_files_of_nettype(tmp_path):
	"""nettype '' reads net_weight/net_iweight; nettype 'in' uses the steady-state networks and drops subsets without them."""
	data, work, fsub, expr = build_tofile_layout(tmp_path / 'u', STATIC_SUBSETS, nettype='')
	fo = str(tmp_path / 'u.h5')
	network.tofile(data, work, fsub, fo, nettype='')
	h = h5_tree(fo)
	w, m = expected_edges(TF_NETS, ['s1', 's2'])
	assert_array_equal(h['prop/es/w_'], w)
	assert_array_equal(h['prop/es/w'], w)
	assert_array_equal(h['prop/es/mask_i'], expected_edges(TF_INETS, ['s1', 's2'])[1])
	data, work, fsub, expr = build_tofile_layout(tmp_path / 'n', STATIC_SUBSETS)
	fo = str(tmp_path / 'in.h5')
	network.tofile(data, work, fsub, fo, nettype='in')
	h = h5_tree(fo)
	assert_array_equal(h['sname'], ['s1'])
	w, m = expected_edges(TF_INETS, ['s1'])
	assert_array_equal(h['prop/es/w'], w[:, :4])
	assert_array_equal(h['prop/es/mask'], m[:, :4])
	assert_array_equal(h['nname'], ['a', 'b', 'c', 'd'])


def test_tofile_optional_and_cell_properties(tmp_path):
	"""optional='' leaves out read counts; fi_c adds one cell property per column, ordered by cell name."""
	data, work, fsub, expr = build_tofile_layout(tmp_path, STATIC_SUBSETS)
	fc = str(tmp_path / 'cells.tsv')
	cells = ['c' + str(x) for x in range(8)]
	props = pd.DataFrame({'depth': np.arange(8.) * 10, 'batch': ['x', 'y'] * 4}, index=cells).iloc[::-1]
	props.to_csv(fc, sep='\t')
	fo = str(tmp_path / 'net.h5')
	network.tofile(data, work, fsub, fo, optional='', fi_c=fc)
	h = h5_tree(fo)
	assert 'prop/nc/readcount' not in h
	assert_allclose(h['prop/c/depth'], np.arange(7.) * 10)
	assert_array_equal(h['prop/c/batch'], (['x', 'y'] * 4)[:7])


def test_tofile_rejects_invalid_input(tmp_path):
	"""Unknown optional data or cells missing from fi_c raise ValueError; no subsets or no networks raise RuntimeError."""
	data, work, fsub, expr = build_tofile_layout(tmp_path / 'a', STATIC_SUBSETS)
	fo = str(tmp_path / 'net.h5')
	with pytest.raises(ValueError):
		network.tofile(data, work, fsub, fo, optional='readcount,unknown')
	fc = str(tmp_path / 'cells.tsv')
	pd.DataFrame({'depth': [1., 2.]}, index=['c0', 'c1']).to_csv(fc, sep='\t')
	with pytest.raises(ValueError):
		network.tofile(data, work, fsub, fo, fi_c=fc)
	with pytest.raises(RuntimeError):
		network.tofile(data, work, os.path.join(work, 'absent.txt'), fo)
	data, work, fsub, expr = build_tofile_layout(tmp_path / 'b', {'s3': ['c5', 'c6']})
	with pytest.raises(RuntimeError):
		network.tofile(data, work, fsub, fo)


DYNAMIC_SUBSETS = {'s1': ['c0', 'c1', 'c2', 'c3'], 's2': ['c4', 'c2', 'c3', 'c5']}


def build_dynamic_layout(root):
	"""tofile layout plus the trajectory, cell and state locations, UMAP coordinates, cell types and state neighbours."""
	import h5py
	data, work, fsub, expr = build_tofile_layout(root, DYNAMIC_SUBSETS)
	TF_INETS['s1'].to_csv(os.path.join(work, 's2', 'net_inweight.tsv.gz'), sep='\t')
	cells = ['c' + str(x) for x in range(8)]
	coord = pd.DataFrame({'x': np.arange(8.), 'y': -np.arange(8.) / 2}, index=cells).iloc[[7, 3, 1, 0, 2, 6, 5, 4]]
	coord.to_csv(os.path.join(data, 'coord_rna.tsv.gz'), sep='\t')
	write_lines(os.path.join(data, 'subsets.txt'), ['typeA', 'typeB'])
	for xi, xj in [('typeA', cells[:3]), ('typeB', cells[3:])]:
		os.makedirs(os.path.join(data, 'subsets', xi))
		write_lines(os.path.join(data, 'subsets', xi, 'names_rna.txt'), xj)
	pd.DataFrame([[0, 1], [1, 0]], index=['s1', 's2'], columns=['s1', 's2']).to_csv(os.path.join(work, 'subset_edges.tsv.gz'), sep='\t')
	pts = {
		os.path.join(data, 'traj_node.h5'): {'edges': np.array([[0, 1], [1, 2]]), 'lens': np.array([1., 2.])},
		os.path.join(data, 'traj_cell_rna.h5'): {'edges': np.array([0, 0, 0, 1, 1, 1]), 'locs': np.array([.1, .5, .9, .2, 1., 1.8])},
		os.path.join(work, 'subset_locs.h5'): {'edges': np.array([0, 1]), 'locs': np.array([.5, 1.])}}
	for xi, xj in pts.items():
		with h5py.File(xi, 'w') as f:
			for xk, xl in xj.items():
				f.create_dataset(xk, data=xl)
	return data, work, fsub, expr, coord


def test_tofile_dynamic_network(tmp_path):
	"""dynamic=True stores the trajectory, cell and state points, coordinates, cell types and state neighbours."""
	data, work, fsub, expr, coord = build_dynamic_layout(tmp_path)
	fo = str(tmp_path / 'net.h5')
	network.tofile(data, work, fsub, fo, dynamic=True)
	h = h5_tree(fo)
	cells = ['c0', 'c1', 'c2', 'c3', 'c4', 'c5']
	assert_array_equal(h['cname'], cells)
	assert_array_equal(h['traj/edges'], [[0, 1], [1, 2]])
	assert_array_equal(h['traj/lens'], [1., 2.])
	assert_array_equal(h['point/c/edges'], [0, 0, 0, 1, 1, 1])
	assert_allclose(h['point/c/locs'], [.1, .5, .9, .2, 1., 1.8])
	assert_array_equal(h['point/s/edges'], [0, 1])
	assert_allclose(h['point/s/locs'], [.5, 1.])
	assert_allclose(h['prop/c/coord'], coord.loc[cells].values.T)
	assert_array_equal(h['prop/c/type'], ['typeA'] * 3 + ['typeB'] * 3)
	assert_array_equal(h['prop/ss/traj-neighbor'], [[0, 1], [1, 0]])
	member = np.array([[x in DYNAMIC_SUBSETS[y] for x in cells] for y in ['s1', 's2']], dtype=float)
	assert_array_equal(h['prop/sc/w'], member)
	w, m = expected_edges({'s1': TF_INETS['s1'], 's2': TF_INETS['s1']}, ['s1', 's2'])
	assert_array_equal(h['prop/es/w_in'], w)


def test_cli_tofile_all_options(tmp_path, run_dictys):
	"""`dictys network tofile` with every option writes the dynamic network with cell properties and without read counts."""
	data, work, fsub, expr, coord = build_dynamic_layout(tmp_path)
	fc = str(tmp_path / 'cells.tsv')
	pd.DataFrame({'depth': np.arange(8.)}, index=['c' + str(x) for x in range(8)]).to_csv(fc, sep='\t')
	fo = str(tmp_path / 'net.h5')
	run_dictys('network', 'tofile', data, work, fsub, fo, '--dynamic', '--nettype', 'n', '--optional', '', '--fi_c', fc)
	h = h5_tree(fo)
	assert 'traj/edges' in h and 'prop/nc/readcount' not in h
	assert_allclose(h['prop/c/depth'], np.arange(6.))
	w, m = expected_edges(TF_NETS, ['s1', 's2'])
	assert_array_equal(h['prop/es/w'], w)
	data, work, fsub, expr = build_tofile_layout(tmp_path / 'static', STATIC_SUBSETS, nettype='')
	fo = str(tmp_path / 'static.h5')
	run_dictys('network', 'tofile', data, work, fsub, fo, '--nettype', '')
	h = h5_tree(fo)
	assert 'traj/edges' not in h
	assert 'prop/es/w_n' not in h
	assert_array_equal(h['prop/es/w_'], w)
	assert_array_equal(h['prop/nc/readcount'], expr.loc[TF_GENES, ['c' + str(x) for x in range(7)]].values)
