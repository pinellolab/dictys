"""
Tests for dictys.net.stat: statistics computed from networks at states or
points on a trajectory.

Expected values are worked out by hand on a three-state network, or with
numpy and networkx reference computations inside each test.
"""

import operator

import networkx as nx
import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from dictys.net import network, stat
from dictys.traj import trajectory, point

CNAME = np.array(['c0', 'c1', 'c2', 'c3', 'c4', 'c5'])
SNAME = np.array(['s0', 's1', 's2'])
NNAME = np.array(['TF1', 'TF2', 'G3', 'G4', 'G5', 'G6'])
STATES = np.arange(3)


def edge_w():
	"""Edge weights (2 regulators x 6 targets x 3 states); G6 never has an edge."""
	w = np.zeros((2, 6, 3))
	w[:, :, 0] = [[0, 0.1, 0.9, -0.4, 0, 0], [0.7, 0, -0.2, 0, 0.3, 0]]
	w[:, :, 1] = [[0, 0.5, 0, 0, 0, 0], [0, 0, 0, 0.6, 0, 0]]
	w[:, :, 2] = [[0, -0.8, 0.2, 0, 0.05, 0], [0, 0, 0, 0, 0, 0]]
	return w


def edge_mask():
	"""Tested edges: all but TF2->G5 at s0, all at s1, TF1 to TF1..G4 at s2."""
	m = np.ones((2, 6, 3), dtype=bool)
	m[1, 4, 0] = False
	m[1, :, 2] = False
	m[0, 4:, 2] = False
	return m


CPM = np.array([
	[1e5, 2e5, 4e5],
	[3e5, 0, 1e5],
	[2e5, 3e5, 1e5],
	[4e5, 4e5, 3e5],
	[0, 1e5, 1e5],
	[5e-3, 1, 2],
])


def make_net(w=None, **ka):
	w = edge_w() if w is None else w
	params = dict(
		cname=CNAME.copy(), sname=SNAME.copy(), nname=NNAME.copy(),
		nids=[np.array([0, 1]), np.arange(6)],
		esprop={'w': w, 'mask': edge_mask()},
		nsprop={'cpm': CPM.copy(), 'x': np.arange(36, dtype=float).reshape(2, 6, 3)},
		scprop={'w': np.repeat(np.eye(3), 2, axis=1)},
		cprop={'coord': np.arange(12, dtype=float).reshape(2, 6)},
	)
	params.update(ka)
	return network(**params)


@pytest.fixture
def d():
	"""Static network with regulators TF1, TF2, targets TF1..G6 and states s0..s2."""
	return make_net()


@pytest.fixture
def chain():
	"""Trajectory 0-1-2 with unit edges; its nodes stand for states s0..s2."""
	return trajectory(np.array([[0, 1], [1, 2]]), np.array([1.0, 1.0]))


def const2():
	"""A 3 x 2 constant stat with rows a, b, c and columns x, y."""
	return stat.const(np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]), [['a', 'b', 'c'], ['x', 'y']], label='K')


def identity_layout(m, pos):
	return pos.copy()


# base


def test_base_is_abstract():
	"""base cannot be instantiated without its abstract methods."""
	with pytest.raises(TypeError):
		stat.base()


def test_default_lims():
	"""default_lims spans finite values expanded by the relative amount on each side."""
	s = stat.const(np.array([[1.0, np.nan], [3.0, 4.0]]), [['a', 'b'], ['x', 'y']])
	assert_allclose(s.default_lims(STATES), [1 - 0.06, 4 + 0.06])
	assert_allclose(s.default_lims(STATES, expansion=0), [1, 4])
	assert_allclose(s.default_lims(STATES, names=[['b'], ['x', 'y']], expansion=0.5), [3 - 0.5, 4 + 0.5])
	allnan = stat.const(np.array([np.nan, np.nan]), [['a', 'b']])
	assert np.isnan(allnan.default_lims(STATES)).all()
	with pytest.raises(NotImplementedError):
		s.default_lims()


@pytest.mark.parametrize('op,sym', [
	(operator.add, '+'), (operator.sub, '-'), (operator.mul, '*'), (operator.truediv, '/'),
	(operator.lt, '<'), (operator.le, '<='), (operator.eq, '=='), (operator.ne, '!='),
	(operator.gt, '>'), (operator.ge, '>='),
])
def test_arithmetic_between_stats(op, sym):
	"""Operators between two stats compute elementwise and label the result '(a)op(b)'."""
	a = stat.const(np.array([1.0, 2.0, 3.0]), [['x', 'y', 'z']], label='A')
	b = stat.const(np.array([3.0, 2.0, 1.0]), [['x', 'y', 'z']], label='B')
	s = op(a, b)
	assert isinstance(s, stat.function)
	assert s.label == f'(A){sym}(B)'
	assert s.isconst
	assert_array_equal(s.compute(STATES[:2]), np.repeat(op(np.array([1.0, 2.0, 3.0]), np.array([3.0, 2.0, 1.0]))[:, None], 2, axis=1))


@pytest.mark.parametrize('op', [operator.and_, operator.or_])
def test_logical_between_stats(op):
	"""& and | combine boolean stats elementwise."""
	a = stat.const(np.array([True, True, False, False]), [list('wxyz')])
	b = stat.const(np.array([True, False, True, False]), [list('wxyz')])
	assert_array_equal(op(a, b).compute(STATES)[:, 0], op(np.array([True, True, False, False]), np.array([True, False, True, False])))


@pytest.mark.parametrize('op', [operator.add, operator.sub, operator.mul, operator.truediv, operator.gt])
def test_arithmetic_with_number(op):
	"""Operators between a stat and a number apply the number to every value."""
	a = stat.const(np.array([1.0, 2.0, 4.0]), [['x', 'y', 'z']])
	assert_array_equal(op(a, 2).compute(STATES)[:, 1], op(np.array([1.0, 2.0, 4.0]), 2))


def test_getitem_by_names_slices_and_indices():
	"""Indexing subsets each axis by name lists, slices or index lists and drops axes indexed by an int."""
	s = const2()
	val = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
	t = s[['c', 'a']]
	assert [list(x) for x in t.names] == [['c', 'a'], ['x', 'y']]
	assert t.label == 'K'
	assert_array_equal(t.compute(STATES[:2]), np.repeat(val[[2, 0]][:, :, None], 2, axis=2))
	t = s[0:2]
	assert [list(x) for x in t.names] == [['a', 'b'], ['x', 'y']]
	assert_array_equal(t.compute(STATES[:1])[:, :, 0], val[:2])
	t = s[[1], ['y']]
	assert [list(x) for x in t.names] == [['b'], ['y']]
	assert_array_equal(t.compute(STATES)[:, :, 2], [[4.0]])
	t = s[1]
	assert [list(x) for x in t.names] == [['x', 'y']]
	assert_array_equal(t.compute(STATES)[:, 0], [3.0, 4.0])
	t = s[:, 0]
	assert [list(x) for x in t.names] == [['a', 'b', 'c']]
	assert_array_equal(t.compute(STATES)[:, 0], [1.0, 3.0, 5.0])


def test_getitem_errors():
	"""Indexing with too many keys raises KeyError and with float keys raises TypeError."""
	s = const2()
	with pytest.raises(KeyError):
		s[0, 0, 0]
	with pytest.raises(TypeError):
		s[[0.5]]


# const


def test_const():
	"""const repeats its value along a last axis of len(pts) and keeps its names and label."""
	val = np.array([[1.0, 2.0], [3.0, 4.0]])
	s = stat.const(val, [['a', 'b'], ['x', 'y']])
	assert s.label == s.default_label() == 'Constant'
	assert [list(x) for x in s.default_names()] == [['a', 'b'], ['x', 'y']]
	assert s.ndict == [{'a': 0, 'b': 1}, {'x': 0, 'y': 1}]
	assert s.isconst
	ans = s.compute(np.array([4, 0, 1]))
	assert ans.shape == (2, 2, 3)
	assert all((ans[:, :, x] == val).all() for x in range(3))
	assert stat.const(val, [['a', 'b'], ['x', 'y']], label='L').label == 'L'


# function


def test_function_compute_and_extra_arguments():
	"""function applies func to the computed stats followed by its extra arguments."""
	a = stat.const(np.array([1.0, 2.0]), [['x', 'y']])
	b = stat.const(np.array([10.0, 20.0]), [['x', 'y']])
	s = stat.function(lambda u, v, k: (u + v) * k, [a, b], 3)
	assert s.label == s.default_label() == 'Function'
	assert [list(x) for x in s.names] == [['x', 'y']]
	assert_array_equal(s.compute(STATES[:2]), [[33.0, 33.0], [66.0, 66.0]])


def test_function_default_names_are_sorted_intersection():
	"""Without names, function uses the sorted intersection of its stats' names on every axis."""
	a = stat.const(np.zeros(3), [['z', 'y', 'x']])
	b = stat.const(np.zeros(3), [['y', 'w', 'z']])
	s = stat.function(np.add, [a, b], label='S')
	assert [list(x) for x in s.names] == [['y', 'z']]
	assert s.label == 'S'
	s = stat.function(np.add, [a, b], names=[['q', 'r', 's']])
	assert [list(x) for x in s.names] == [['q', 'r', 's']]
	with pytest.raises(NotImplementedError):
		stat.function(np.add, [a, const2()])


def test_function_isconst(d):
	"""function is constant only when all its stats are, unless isconst overrides it."""
	a = stat.const(CPM[:, 0], [NNAME])
	b = stat.lcpm(d, cut=-1)
	assert stat.function(np.add, [a, a]).isconst
	assert not stat.function(np.add, [a, b]).isconst
	assert stat.function(np.add, [a, b], isconst=True).isconst
	assert not stat.function(np.add, [a, a], isconst=False).isconst


# fsinglestat, finitial, fdiff, fmasked


def test_fsinglestat(d):
	"""fsinglestat turns a function of the values at given points into a constant stat."""
	s0 = stat.sprop(d, 'ns', 'cpm')
	s = stat.fsinglestat(lambda x: x.max(axis=-1), s0, np.array([0, 2]))
	assert isinstance(s, stat.const)
	assert s.label == 'Single value of cpm'
	assert_array_equal(s.compute(np.arange(4))[:, 3], np.maximum(CPM[:, 0], CPM[:, 2]))
	assert [list(x) for x in s.names] == [list(NNAME)]
	assert stat.fsinglestat(lambda x: x.min(axis=-1), s0, STATES, label='low').label == 'low'


def test_finitial(d):
	"""finitial holds the value at the first of the given points for every point."""
	s0 = stat.sprop(d, 'ns', 'cpm')
	s = stat.finitial(s0, np.array([2, 0, 1]))
	assert_array_equal(s.compute(STATES), np.repeat(CPM[:, [2]], 3, axis=1))
	assert s.label == 'Single value of cpm'
	assert stat.finitial(s0, STATES, label='init').label == 'init'


def sorted_lcpm(d):
	"""lcpm(cut=-1) with its genes put in sorted order by indexing, and the row order of CPM that matches."""
	# function stats name their rows by sorted names without reordering values (a bug in the base), so the
	# tests that build them start from stats whose names are already sorted.
	order = np.argsort(NNAME)
	return stat.lcpm(d, cut=-1)[list(NNAME[order])], order


def test_fdiff(d):
	"""fdiff subtracts a base stat and labels the result as a difference unless given a label."""
	s0, order = sorted_lcpm(d)
	s = stat.fdiff(s0, stat.finitial(s0, STATES))
	assert list(s.names[0]) == sorted(NNAME)
	expected = np.log2(CPM[order] + 1)
	assert_allclose(s.compute(STATES), expected - expected[:, [0]], atol=1e-12)
	assert s.label == 'Difference in Log2 CPM'
	assert stat.fdiff(s0, s0, label='zero').label == 'zero'
	assert_array_equal(stat.fdiff(s0, s0).compute(STATES), 0)


def test_fmasked(d):
	"""fmasked replaces values where the mask stat is False by nan."""
	s0, order = sorted_lcpm(d)
	s = stat.fmasked(s0, s0 > 18, label='high')
	assert list(s.names[0]) == sorted(NNAME)
	expected = np.log2(CPM[order] + 1)
	expected[expected <= 18] = np.nan
	assert_allclose(s.compute(STATES), expected)
	assert s.label == 'high'


# fsmooth


def test_fsmooth_linear_on_trajectory(d, chain):
	"""Linear smoothing over trajectory nodes interpolates node values along each edge."""
	s = stat.fsmooth(stat.sprop(d, 'ns', 'cpm'), chain, ['linear', [], {}])
	assert s.label == 'cpm'
	assert [list(x) for x in s.names] == [list(NNAME)]
	pts = point(chain, np.array([0, 0, 1]), np.array([0.25, 1.0, 0.5]))
	expected = np.array([0.75 * CPM[:, 0] + 0.25 * CPM[:, 1], CPM[:, 1], 0.5 * CPM[:, 1] + 0.5 * CPM[:, 2]]).T
	assert_allclose(s.compute(pts), expected, rtol=1e-12)


def test_fsmooth_gaussian_and_node_subset(d, chain):
	"""Gaussian smoothing weights nodes by exp(-(distance/radius)^2/2); nodes_path limits the nodes used."""
	pts = point(chain, np.array([0, 1]), np.array([0.3, 0.9]))
	pos = np.array([0.3, 1.9])
	dist = np.abs(pos[None, :] - np.array([0, 1, 2])[:, None])
	k = np.exp(-(dist / 0.7) ** 2 / 2)
	s = stat.fsmooth(stat.sprop(d, 'ns', 'cpm'), chain, ['conv', [0.7], {}], label='smooth')
	assert s.label == 'smooth'
	assert_allclose(s.compute(pts), CPM @ (k / k.sum(axis=0)), rtol=1e-12)
	s = stat.fsmooth(stat.sprop(d, 'ns', 'cpm'), chain, ['conv', [0.7], {'nodes_path': [0, 1], 'criterion_path': 'strict'}])
	assert_allclose(s.compute(pts), CPM[:, :2] @ (k[:2] / k[:2].sum(axis=0)), rtol=1e-12)


def test_fsmooth_on_points(d, chain):
	"""Smoothing from a point list interpolates values computed at those points."""
	base = chain.linspace(0, 2, 3)
	s = stat.fsmooth(stat.pseudotime(d, base), base, ['linear', [], {}])
	pts = point(chain, np.array([0, 1]), np.array([0.25, 0.6]))
	assert_allclose(s.compute(pts), np.repeat([[0.25, 1.6]], 6, axis=0), rtol=1e-12)


def test_fsmooth_errors(d, chain):
	"""fsmooth needs a trajectory or point list, and computes only at points."""
	with pytest.raises(TypeError):
		stat.fsmooth(stat.sprop(d, 'ns', 'cpm'), STATES, ['linear', [], {}])
	s = stat.fsmooth(stat.sprop(d, 'ns', 'cpm'), chain, ['linear', [], {}])
	with pytest.raises(TypeError):
		s.compute(STATES)


# fbinarize


def test_fbinarize_top_edges_by_absolute_value(d):
	"""fbinarize marks the top int(sparsity*12) edges by absolute weight at each state."""
	s = stat.fbinarize(stat.net(d), sparsity=0.25)
	assert s.label == 'Network edge strength'
	assert [list(x) for x in s.names] == [['TF1', 'TF2'], list(NNAME)]
	ans = s.compute(np.array([0, 2]))
	assert ans.dtype == bool
	expected = np.zeros((2, 6, 2), dtype=bool)
	# s0: 0.9 TF1->G3, 0.7 TF2->TF1, |-0.4| TF1->G4; s2: |-0.8| TF1->TF2, 0.2 TF1->G3, 0.05 TF1->G5
	expected[[0, 1, 0], [2, 0, 3], 0] = True
	expected[[0, 0, 0], [1, 2, 4], 1] = True
	assert_array_equal(ans, expected)


def test_fbinarize_unsigned(d):
	"""With signed=False the top edges are the largest raw weights, so negative edges rank last."""
	ans = stat.fbinarize(stat.net(d), sparsity=0.25, signed=False).compute(np.array([0]))
	expected = np.zeros((2, 6, 1), dtype=bool)
	expected[[0, 1, 1], [2, 0, 4], 0] = True
	assert_array_equal(ans, expected)


def test_fbinarize_without_sparsity(d):
	"""sparsity=None marks every non-zero edge (signed) or every positive edge (unsigned)."""
	w = edge_w()
	assert_array_equal(stat.fbinarize(stat.net(d), sparsity=None).compute(STATES), w != 0)
	assert_array_equal(stat.fbinarize(stat.net(d), sparsity=None, signed=False).compute(STATES), w > 0)


def test_fbinarize_with_mask():
	"""With statmask the edge count is int(sparsity * tested edges) at each state; untested edges have zero weight."""
	d = make_net(w=edge_w() * edge_mask())
	s = stat.fbinarize(stat.net(d), statmask=stat.netmask(d), sparsity=0.45, label='bin')
	assert s.label == 'bin'
	ans = s.compute(np.array([0, 2]))
	# s0: 11 tested edges -> top 4 (12 edges would give 5); s2: 4 tested edges -> top 1 (12 would give 5)
	expected = np.zeros((2, 6, 2), dtype=bool)
	expected[[0, 1, 0, 1], [2, 0, 3, 2], 0] = True
	expected[0, 1, 1] = True
	assert_array_equal(ans, expected)


# pseudotime


def test_pseudotime(d, chain):
	"""pseudotime is the trajectory distance from the first reference point, the same for every gene."""
	pts = chain.linspace(0, 2, 5)
	s = stat.pseudotime(d, pts)
	assert s.label == 'Pseudo-time'
	assert [list(x) for x in s.names] == [list(NNAME)]
	assert_allclose(s.compute(pts), np.repeat([[0, 0.5, 1, 1.5, 2]], 6, axis=0), atol=1e-12)
	assert_allclose(s.compute(np.array([2, 0, 1])), np.repeat([[2, 0, 1]], 6, axis=0), atol=1e-12)


def test_pseudotime_reverse_and_explicit_trajectory(d, chain):
	"""Distances are measured from pts[0]; the traj argument converts state indices to nodes."""
	pts = chain.linspace(2, 0, 3)
	s = stat.pseudotime(d, pts, traj=chain, names=[['g']], label='T')
	assert s.label == 'T'
	assert_allclose(s.compute(np.array([0, 1, 2])), [[2, 1, 0]], atol=1e-12)


# lcpm


def test_lcpm(d):
	"""lcpm is log2(CPM + constant) with CPM below the cut hidden as nan."""
	s = stat.lcpm(d)
	assert s.label == 'Log2 CPM'
	assert [list(x) for x in s.names] == [list(NNAME)]
	expected = np.log2(CPM + 1)
	expected[CPM < 0.01] = np.nan
	assert_allclose(s.compute(STATES), expected)
	assert_allclose(stat.lcpm(d, cut=-1, constant=0.5).compute(np.array([2, 1])), np.log2(CPM[:, [2, 1]] + 0.5))
	expected = np.log2(CPM + 1)
	expected[CPM < 2.5] = np.nan
	assert_allclose(stat.lcpm(d, cut=2.5).compute(STATES), expected)


def test_lcpm_errors(d):
	"""lcpm refuses points, unknown genes and networks without CPM."""
	chain = trajectory(np.array([[0, 1]]), np.array([1.0]))
	with pytest.raises(TypeError):
		stat.lcpm(d).compute(chain.topoint())
	with pytest.raises(ValueError):
		stat.lcpm(d, names=[['TF1', 'missing']]).compute(STATES)
	d2 = make_net(nsprop={'x': np.zeros((6, 3))})
	with pytest.raises(ValueError, match='CPM results not found'):
		stat.lcpm(d2).compute(STATES)


# sprop, net, netmask


def test_sprop_state_cell_property(d):
	"""sprop moves the state axis last: a state-cell property gives (cells, states)."""
	s = stat.sprop(d, 'sc', 'w')
	assert s.label == 'w'
	assert [list(x) for x in s.names] == [list(CNAME)]
	expected = np.repeat(np.eye(3), 2, axis=1)
	assert_array_equal(s.compute(np.array([2, 0])), expected[[2, 0]].T)


def test_sprop_name_prefix(d):
	"""names_pref names the leading axes of a property with extra dimensions."""
	s = stat.sprop(d, 'ns', 'x', names_pref=[['a', 'b']], label='X')
	assert s.label == 'X'
	assert [list(x) for x in s.names] == [['a', 'b'], list(NNAME)]


def test_sprop_errors(d, chain):
	"""sprop rejects state-free properties, wrong prefixes, and points."""
	with pytest.raises(ValueError, match='only for state-dependent'):
		stat.sprop(d, 'c', 'coord')
	with pytest.raises(ValueError, match='Incorrect prefix dimension'):
		stat.sprop(d, 'ns', 'x')
	with pytest.raises(ValueError, match='Incorrect prefix shape'):
		stat.sprop(d, 'ns', 'x', names_pref=[['a', 'b', 'c']])
	with pytest.raises(ValueError):
		stat.sprop(d, 'ns', 'cpm').compute(chain.topoint())


def test_net_and_netmask(d):
	"""net and netmask read edge weights and masks with regulator and target names."""
	s = stat.net(d)
	assert s.label == 'Network edge strength'
	assert [list(x) for x in s.names] == [['TF1', 'TF2'], list(NNAME)]
	assert_array_equal(s.compute(np.array([1, 0])), edge_w()[:, :, [1, 0]])
	m = stat.netmask(d)
	assert m.label == 'Network edge mask'
	assert_array_equal(m.compute(STATES), edge_mask())
	assert_array_equal(stat.net(d, varname='mask').compute(STATES), edge_mask())
	assert_array_equal(stat.netmask(d, varname='w').compute(STATES), edge_w())


# fcentrality_base and its wrappers


# fcentrality_base builds its graph on regulator and target indices as if regulator i were target i
# (a bug in the base). The fixture's regulators TF1, TF2 are its first two targets, so the mix-up has no
# effect there (expected_centrality checks this); the reference graphs are built on gene names.


def named_graph(d, w, directed=False):
	"""networkx graph of the non-zero edges of one regulator x target matrix, nodes named by gene."""
	g = (nx.DiGraph if directed else nx.Graph)()
	for i, j in zip(*np.nonzero(w)):
		g.add_edge(d.nname[d.nids[0][i]], d.nname[d.nids[1][j]])
	return g


def expected_centrality(d, func, names, directed=False):
	assert_array_equal(d.nids[0], d.nids[1][:len(d.nids[0])])
	ans = np.full((len(names), 3), np.nan)
	for k in range(3):
		c = func(named_graph(d, edge_w()[:, :, k], directed=directed))
		for i, x in enumerate(names):
			if x in c:
				ans[i, k] = c[x]
	return ans


@pytest.mark.parametrize('wrapper,func,label', [
	(stat.fcentrality_eigenvector, nx.eigenvector_centrality, 'Eigenvalue centrality'),
	(stat.fcentrality_betweenness, nx.betweenness_centrality, 'Betweenness centrality'),
	(stat.fcentrality_closeness, nx.closeness_centrality, 'Closeness centrality'),
])
def test_centrality_wrappers(d, wrapper, func, label):
	"""Regulator centralities match networkx on the graph of non-zero edges at each state."""
	s = wrapper(stat.net(d))
	assert s.label == label
	assert [list(x) for x in s.names] == [['TF1', 'TF2']]
	assert_allclose(s.compute(STATES), expected_centrality(d, func, ['TF1', 'TF2']), rtol=1e-4)
	assert wrapper(stat.net(d), label='C').label == 'C'


@pytest.mark.parametrize('func', [nx.betweenness_centrality, nx.closeness_centrality])
def test_centrality_targets_directed(d, func):
	"""roleaxis=1 and directed=True give target centralities on the directed graph; absent nodes are nan."""
	s = stat.fcentrality_base(stat.net(d), func, directed=True, roleaxis=1)
	assert s.label == 'Centrality'
	assert [list(x) for x in s.names] == [list(NNAME)]
	ans = s.compute(STATES)
	assert_allclose(ans, expected_centrality(d, func, NNAME, directed=True), rtol=1e-10)
	assert np.isnan(ans[5]).all()


# fcentrality_degree and flnneighbor


def test_degree_centrality(d):
	"""Degree centrality sums absolute weights over targets (outdegree) or regulators (indegree) plus a constant."""
	a = np.abs(edge_w())
	s = stat.fcentrality_degree(stat.net(d))
	assert s.label == 'Outdegree centrality'
	assert [list(x) for x in s.names] == [['TF1', 'TF2']]
	assert_allclose(s.compute(STATES), a.sum(axis=1))
	s = stat.fcentrality_degree(stat.net(d), roleaxis=1, constant=2)
	assert s.label == 'Indegree centrality'
	assert [list(x) for x in s.names] == [list(NNAME)]
	assert_allclose(s.compute(STATES), a.sum(axis=0) + 2)


def test_degree_centrality_weighted_sparsity(d):
	"""weighted_sparsity rescales weights so their total equals sparsity times the number of edges."""
	a = np.abs(edge_w())
	ans = stat.fcentrality_degree(stat.net(d), weighted_sparsity=0.1).compute(STATES)
	assert_allclose(ans.sum(axis=0), 0.1 * 12, rtol=1e-12)
	assert_allclose(ans, a.sum(axis=1) * 0.1 / a.mean(axis=(0, 1)), rtol=1e-12)


def test_degree_centrality_rate_with_mask(d):
	"""With statmask the stat is (masked degree + constant) / (number of tested partner nodes + constant)."""
	a = np.abs(edge_w())
	m = edge_mask()
	ans = stat.fcentrality_degree(stat.net(d), statmask=stat.netmask(d), constant=1, weighted_sparsity=5).compute(STATES)
	# Targets tested by any regulator: 6, 6, 4
	assert_allclose(ans, ((a * m).sum(axis=1) + 1) / (np.array([6, 6, 4]) + 1), rtol=1e-12)
	ans = stat.fcentrality_degree(stat.net(d), roleaxis=1, statmask=stat.netmask(d)).compute(STATES)
	# Regulators testing any target: 2, 2, 1
	assert_allclose(ans, (a * m).sum(axis=0) / np.array([2, 2, 1]), rtol=1e-12)


def test_flnneighbor(d):
	"""flnneighbor is log2(outdegree + constant), with optional relative rate and weighting."""
	binary = stat.fbinarize(stat.net(d), sparsity=None)
	s = stat.flnneighbor(binary)
	assert s.label == 'Log2 (Outdegree + 1)'
	assert [list(x) for x in s.names] == [['TF1', 'TF2']]
	nonzero = (edge_w() != 0).sum(axis=1)
	assert_allclose(s.compute(STATES), np.log2(nonzero + 1))
	assert_allclose(stat.flnneighbor(binary, constant=3, label='L').compute(STATES), np.log2(nonzero + 3))
	assert stat.flnneighbor(binary, label='L').label == 'L'
	masked = ((edge_w() != 0) & edge_mask()).sum(axis=1)
	ans = stat.flnneighbor(binary, statmask=stat.netmask(d)).compute(STATES)
	assert_allclose(ans, np.log2(masked + 1) - np.log2(np.array([6, 6, 4]) + 1), rtol=1e-12)
	a = np.abs(edge_w())
	ans = stat.flnneighbor(stat.net(d), weighted_sparsity=0.2).compute(STATES)
	assert_allclose(ans, np.log2(a.sum(axis=1) * 0.2 / a.mean(axis=(0, 1)) + 1), rtol=1e-12)
	# Targets sorted first: flnneighbor is a function stat, which names its rows in sorted order (a bug in the base)
	order = np.argsort(NNAME)
	s = stat.flnneighbor(binary[:, list(NNAME[order])], roleaxis=1)
	assert list(s.names[0]) == sorted(NNAME)
	assert_allclose(s.compute(STATES), np.log2((edge_w() != 0).sum(axis=0) + 1)[order])


# flayout_base


def shown_nodes(w):
	"""Sorted names of genes with an edge in the regulator x target matrix w."""
	t = w != 0
	nodes = set(NNAME[:2][t.any(axis=1)]) | set(NNAME[t.any(axis=0)])
	return sorted(nodes)


def test_layout_names_and_size_scaling(d):
	"""Layout keeps genes with any edge, gives nan to genes without an edge at a state, and scales to unit RMS."""
	np.random.seed(0)
	s = stat.flayout_base(stat.net(d), identity_layout, pts=STATES)
	assert s.label == 'Coordinates'
	names = list(s.names[0])
	assert names == ['G3', 'G4', 'G5', 'TF1', 'TF2']
	assert list(s.names[1]) == ['Dim 1', 'Dim 2']
	pos = s.compute(STATES)
	assert pos.shape == (5, 2, 3)
	for k in range(3):
		shown = [names.index(x) for x in shown_nodes(edge_w()[:, :, k])]
		hidden = [x for x in range(5) if x not in shown]
		assert np.isnan(pos[hidden, :, k]).all()
		p = pos[shown, :, k]
		assert not np.isnan(p).any()
		assert_allclose(p.mean(axis=0), 0, atol=1e-12)
		assert_allclose(np.sqrt((p ** 2).mean()), 1, rtol=1e-12)


def test_layout_keeps_previous_positions(d):
	"""Nodes shown at consecutive states keep their positions; new nodes start near the existing ones."""
	np.random.seed(1)
	s = stat.flayout_base(stat.net(d), identity_layout, pts=STATES, scale='none', rand_expand=0.1)
	names = list(s.names[0])
	pos = s.compute(STATES)
	tf1, tf2, g4 = [names.index(x) for x in ['TF1', 'TF2', 'G4']]
	assert_allclose(pos[tf1, :, 1] - pos[tf2, :, 1], pos[tf1, :, 0] - pos[tf2, :, 0], atol=1e-12)
	assert_allclose(pos[tf1, :, 1] - pos[g4, :, 1], pos[tf1, :, 0] - pos[g4, :, 0], atol=1e-12)
	assert_allclose(pos[tf1, :, 2] - pos[tf2, :, 2], pos[tf1, :, 1] - pos[tf2, :, 1], atol=1e-12)
	# G3 and G5 are new at s2: they start in the square around the s1 nodes expanded by rand_expand
	old = pos[[tf1, tf2, g4], :, 1]
	shift = pos[tf1, :, 2] - pos[tf1, :, 1]
	center = (old.min(axis=0) + old.max(axis=0)) / 2 + shift
	half = (old.max(axis=0) - old.min(axis=0)).max() / 2 + 0.1
	for x in ['G3', 'G5']:
		assert (np.abs(pos[names.index(x), :, 2] - center) <= half + 1e-12).all()


def test_layout_edge_scaling_and_three_dimensions():
	"""scale='length' gives unit edge length, and ndim sets the number of coordinates."""
	w = np.zeros((2, 6, 3))
	w[0, 2, 0] = 2.0
	w[1, 3, 1] = -0.5
	w[0, 1, 2] = 1.0
	d1 = make_net(w=w)
	np.random.seed(2)
	s = stat.flayout_base(stat.net(d1), identity_layout, ndim=3, pts=STATES, scale='length')
	assert list(s.names[1]) == ['Dim 1', 'Dim 2', 'Dim 3']
	names = list(s.names[0])
	assert names == ['G3', 'G4', 'TF1', 'TF2']
	pos = s.compute(STATES)
	assert pos.shape == (4, 3, 3)
	for k, (a, b) in enumerate([('TF1', 'G3'), ('TF2', 'G4'), ('TF1', 'TF2')]):
		assert_allclose(np.linalg.norm(pos[names.index(a), :, k] - pos[names.index(b), :, k]), 1, rtol=1e-12)


def test_layout_passes_scaled_network_to_layout_function(d):
	"""The layout function gets |weights| scaled to RMS netscale, or raw signed weights with sabs=False and netscale=None."""
	calls = []

	def record(m, pos):
		calls.append(m.copy())
		return pos.copy()

	np.random.seed(3)
	s = stat.flayout_base(stat.net(d), record, pts=STATES, netscale=2.0)
	names = list(s.names[0])
	last = calls[-3:]
	nonzero = np.concatenate([x[x != 0] for x in last])
	assert_allclose(np.sqrt((nonzero ** 2).mean()), 2.0, rtol=1e-12)
	assert (nonzero > 0).all()
	calls.clear()
	s = stat.flayout_base(stat.net(d), record, pts=STATES, netscale=None, sabs=False)
	m = calls[-1]
	shown = [names.index(x) for x in shown_nodes(edge_w()[:, :, 2])]
	assert m.shape == (len(shown), len(shown))
	expected = np.zeros((5, 5))
	for i in range(2):
		for j in range(6):
			if edge_w()[i, j, 2] != 0:
				expected[names.index(NNAME[i]), names.index(NNAME[j])] = edge_w()[i, j, 2]
	assert_allclose(m, expected[np.ix_(shown, shown)])


def test_layout_empty_network_and_unknown_scale():
	"""A network without edges gives nan coordinates, and an unknown scale raises ValueError."""
	d0 = make_net(w=np.zeros((2, 6, 3)))
	s = stat.flayout_base(stat.net(d0), identity_layout, pts=STATES)
	assert np.isnan(s.pos).all()
	np.random.seed(4)
	with pytest.raises(ValueError, match='Unknown scale'):
		stat.flayout_base(stat.net(make_net()), identity_layout, pts=STATES, scale='area')


def test_layout_compute_at_states_and_points(d, chain):
	"""compute maps states or known points to their layouts and rejects unknown points."""
	np.random.seed(5)
	s = stat.flayout_base(stat.net(d), identity_layout, pts=np.array([2, 0, 1]))
	assert_array_equal(s.compute(np.array([0, 2])), s.pos[:, :, [1, 0]])
	assert s.ptsmap(np.array([0])) is None
	smooth = stat.fsmooth(stat.net(d), chain, ['linear', [], {}])
	pts = chain.linspace(0, 2, 3)
	s = stat.flayout_base(smooth, identity_layout, pts=pts)
	assert_array_equal(s.ptsmap(pts[[2, 0]]), [2, 0])
	assert_array_equal(s.compute(pts[[2, 1]]), s.pos[:, :, [2, 1]])
	for k in range(3):
		names = list(s.names[0])
		shown = [names.index(x) for x in shown_nodes(edge_w()[:, :, k])]
		assert not np.isnan(s.pos[shown, :, k]).any()
	with pytest.raises(ValueError):
		s.compute(point(chain, np.array([0]), np.array([0.5])))


def test_layout_init_pos(d):
	"""init_pos lays out connected components side by side, sep apart, with always-shown nodes on the x axis."""
	np.random.seed(6)
	s = stat.flayout_base(stat.net(d), identity_layout, pts=STATES)
	m = np.zeros((6, 6))
	m[0, 1] = 1.0
	m[2, 3] = 2.0
	m[3, 4] = 1.0
	pos = s.init_pos(m, sep=0.5, always=[5])
	assert pos.shape == (6, 2)
	assert not np.isnan(pos).any()
	comps = sorted([[0, 1], [2, 3, 4]], key=lambda x: pos[x, 0].min())
	assert_allclose(pos[comps[0], 0].min(), 0, atol=1e-12)
	assert_allclose(pos[comps[1], 0].min(), pos[comps[0], 0].max() + 0.5, rtol=1e-12)
	for c in comps:
		assert_allclose(pos[c, 1].mean(), 0, atol=1e-12)
	assert_allclose(pos[5], [pos[comps[1], 0].max() + 0.5, 0], rtol=1e-12)
	pos = s.init_pos(m)
	assert np.isnan(pos[5]).all() and not np.isnan(pos[:5]).any()
	assert np.isnan(s.init_pos(np.zeros((3, 3)))).all()
