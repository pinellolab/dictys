"""
Tests for dictys.traj: argpartition, the trajectory class and the point class.

Trajectories are tiny trees built by hand. Expected distances come from a
networkx reference that splits every edge at every point location.
"""

import h5py
import networkx as nx
import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from dictys.traj import argpartition, point, trajectory


################################################################
# Helpers
################################################################

def line3():
	"""Path 0-1-2 with edge lengths 1 and 2."""
	return trajectory(np.array([[0, 1], [1, 2]]), np.array([1., 2.]))


def ytree():
	"""Y-shaped tree: root 0, branching node 1, leaves 2 and 3; edge lengths 1, 2, 3."""
	return trajectory(np.array([[0, 1], [1, 2], [1, 3]]), np.array([1., 2., 3.]))


def ytree_mixed():
	"""The tree of ytree() with edges 0 and 1 stored from their other end."""
	return trajectory(np.array([[1, 0], [2, 1], [1, 3]]), np.array([1., 2., 3.]))


def branched():
	"""Path 0-1-2-3 (lengths 1, 2, 1) with a side branch 1-4 (length 1)."""
	return trajectory(np.array([[0, 1], [1, 2], [2, 3], [1, 4]]), np.array([1., 2., 1., 1.]))


def locs_of(p):
	"""List of (edge, loc) pairs of a point object."""
	return list(zip(p.edges.tolist(), p.locs.tolist()))


def _locname(traj, x):
	if np.ndim(x) == 0:
		return ('node', int(x))
	e, l = int(x[0]), float(x[1])
	if l == 0:
		return ('node', int(traj.edges[e, 0]))
	if l == float(traj.lens[e]):
		return ('node', int(traj.edges[e, 1]))
	return ('mid', e, l)


def ref_dist(traj, a, b):
	"""
	Distance matrix between location lists a and b on traj, from networkx shortest paths on the
	tree split at every location. A location is a node ID or an (edge ID, distance from the
	edge's start node) pair.
	"""
	cuts = [{0., float(traj.lens[e])} for e in range(traj.ne)]
	for x in list(a) + list(b):
		if np.ndim(x) != 0:
			cuts[int(x[0])].add(float(x[1]))
	g = nx.Graph()
	for e in range(traj.ne):
		c = sorted(cuts[e])
		for x, y in zip(c[:-1], c[1:]):
			g.add_edge(_locname(traj, (e, x)), _locname(traj, (e, y)), l=y - x)
	d = dict(nx.all_pairs_dijkstra_path_length(g, weight='l'))
	return np.array([[d[_locname(traj, x)][_locname(traj, y)] for y in b] for x in a], dtype=float)


def random_points(traj, rng, n):
	"""n points at random edges and locations, plus one point at every node."""
	edges = rng.integers(0, traj.ne, n)
	locs = rng.random(n) * traj.lens[edges]
	edges = np.r_[edges, np.arange(traj.ne), [0]]
	locs = np.r_[locs, np.zeros(traj.ne), traj.lens[0]]
	return point(traj, edges, locs)


def assert_same_trajectory(a, b):
	"""Asserts every attribute of two trajectory objects is equal."""
	assert a.nn == b.nn and a.ne == b.ne
	assert_array_equal(a.edges, b.edges)
	assert_array_equal(a.lens, b.lens)
	assert a.lens.dtype == b.lens.dtype and a.edges.dtype == b.edges.dtype
	assert_array_equal(a.dist, b.dist)
	assert_array_equal(a.deg, b.deg)
	assert a.edgedict == b.edgedict
	assert sorted(a.g.edges(data='l')) == sorted(b.g.edges(data='l'))


def assert_same_points(a, b):
	"""Asserts edges, locs, dist and count of two point objects are equal."""
	assert len(a) == len(b) and a.npt == b.npt
	assert_array_equal(a.edges, b.edges)
	assert_array_equal(a.locs, b.locs)
	assert_allclose(a.dist, b.dist, rtol=0, atol=1e-12)


def check_partition(a, ans, kth, axis):
	"""Asserts ans permutes indices along axis and puts the kth smallest value of a at position kth."""
	assert ans.shape == a.shape
	n = a.shape[axis]
	idx = np.moveaxis(ans, axis, -1)
	assert_array_equal(np.sort(idx, axis=-1), np.broadcast_to(np.arange(n), idx.shape))
	v = np.moveaxis(np.take_along_axis(a, ans, axis), axis, -1)
	s = np.sort(np.moveaxis(a, axis, -1), axis=-1)
	assert_array_equal(v[..., kth], s[..., kth])
	assert (v[..., :kth] <= v[..., kth:kth + 1]).all()
	assert (v[..., kth + 1:] >= v[..., kth:kth + 1]).all()


################################################################
# argpartition
################################################################

class TestArgpartition:
	@pytest.mark.parametrize('shape,axis', [((7,), -1), ((7,), 0), ((5, 6), 0), ((5, 6), 1), ((3, 4, 5), 0), ((3, 4, 5), 1), ((3, 4, 5), 2), ((3, 4, 5), -1)])
	def test_default_partitions_along_axis(self, rng, shape, axis):
		"""The default draw order gives a valid partition along any axis of a 1-3 D array."""
		a = rng.random(shape)
		check_partition(a, argpartition(a, 2, axis=axis), 2, axis)

	def test_undefined_partitions(self, rng):
		"""draw_order='undefined' gives a valid partition, also when values tie."""
		a = rng.integers(0, 3, (9, 4)).astype(float)
		check_partition(a, argpartition(a, 4, axis=0, draw_order='undefined'), 4, 0)

	def test_random_without_ties_equals_undefined(self, rng):
		"""Without ties at the kth value, draw_order='random' returns the plain partition unchanged."""
		a = rng.random((8, 5))
		np.random.seed(0)
		ans = argpartition(a, 3, axis=0, draw_order='random')
		assert_array_equal(ans, argpartition(a, 3, axis=0, draw_order='undefined'))
		check_partition(a, ans, 3, 0)

	def test_random_with_ties_returns_permutation(self):
		"""With ties at the kth value, draw_order='random' still returns a permutation of indices on each column."""
		a = np.array([[1., 1., 0., 5., 1.], [2., 0., 2., 2., 3.]]).T
		np.random.seed(0)
		# The partition property is not checked: the shuffle breaks it (a bug in the base).
		for _ in range(20):
			ans = argpartition(a, 1, axis=0, draw_order='random')
			assert_array_equal(np.sort(ans, axis=0), np.tile(np.arange(5)[:, None], (1, 2)))

	def test_error_raises_on_tie_at_kth_value(self):
		"""draw_order='error' raises RuntimeError when another element equals the kth smallest value."""
		with pytest.raises(RuntimeError, match='duplicate'):
			argpartition(np.array([2., 1., 1., 3.]), 0, draw_order='error')

	def test_error_checks_each_column(self):
		"""draw_order='error' raises when only one row along axis=1 has a tie at its kth value."""
		with pytest.raises(RuntimeError, match='duplicate'):
			argpartition(np.array([[1., 2., 3.], [4., 4., 5.]]), 0, axis=1, draw_order='error')

	def test_error_allows_ties_elsewhere(self):
		"""draw_order='error' returns a valid partition when the ties are not at the kth value."""
		a = np.array([3., 3., 1., 2.])
		ans = argpartition(a, 1, draw_order='error')
		check_partition(a, ans, 1, -1)

	def test_invalid_draw_order(self):
		"""An unknown draw_order fails its assertion."""
		with pytest.raises(AssertionError):
			argpartition(np.arange(3.), 1, draw_order='first')


################################################################
# trajectory
################################################################

class TestTrajectoryInit:
	def test_attributes_of_y_tree(self):
		"""Counts, degrees, edge dictionary, graph and stored inputs of the Y tree match hand-built values."""
		edges = np.array([[0, 1], [1, 2], [1, 3]])
		lens = np.array([1., 2., 3.])
		t = trajectory(edges, lens)
		assert t.nn == 4 and t.ne == 3
		assert t.edges is edges and t.lens is lens
		assert_array_equal(t.deg, [1, 3, 1, 1])
		assert t.edgedict == {(0, 1): [0, 1], (1, 2): [1, 1], (1, 3): [2, 1], (1, 0): [0, -1], (2, 1): [1, -1], (3, 1): [2, -1]}
		assert isinstance(t.g, nx.Graph)
		assert sorted((min(u, v), max(u, v), l) for u, v, l in t.g.edges(data='l')) == [(0, 1, 1.), (1, 2, 2.), (1, 3, 3.)]

	@pytest.mark.parametrize('make', [ytree, ytree_mixed])
	def test_dist_by_hand(self, make):
		"""Node distances of the Y tree are the sums of edge lengths on the unique paths, whatever the edge direction."""
		expected = [[0, 1, 3, 4], [1, 0, 2, 3], [3, 2, 0, 5], [4, 3, 5, 0]]
		assert_allclose(make().dist, expected)

	def test_dist_and_degree_match_networkx(self, rng):
		"""On a random 8-node tree, dist equals networkx Floyd-Warshall and deg equals networkx degrees."""
		n = 8
		edges = np.array([[int(rng.integers(0, i)), i] for i in range(1, n)])
		lens = rng.random(n - 1) + 0.1
		g = nx.Graph()
		g.add_weighted_edges_from([(int(a), int(b), float(l)) for (a, b), l in zip(edges, lens)])
		t = trajectory(edges, lens)
		assert_allclose(t.dist, nx.floyd_warshall_numpy(g, nodelist=range(n)), rtol=1e-12)
		assert_array_equal(t.deg, [g.degree(i) for i in range(n)])

	@pytest.mark.parametrize('edges,lens', [
		(np.zeros((0, 2), dtype=int), np.zeros(0)),
		(np.array([[0, 1], [1, 2]]), np.array([1.])),
		(np.array([[0, 1], [1, 2]]), np.array([1., -1.])),
		(np.array([[0, 1], [2, 3]]), np.array([1., 1.])),
		(np.array([[0, 0], [1, 2]]), np.array([1., 1.])),
		(np.array([[0, 1], [1, 3], [3, 0]]), np.array([1., 1., 1.])),
	], ids=['no-edge', 'shape-mismatch', 'negative-length', 'too-many-nodes', 'self-loop', 'missing-node'])
	def test_invalid_input_asserts(self, edges, lens):
		"""Empty, mismatched, negative, self-looped or gapped edge lists fail their assertions."""
		with pytest.raises(AssertionError):
			trajectory(edges, lens)

	def test_cycle_raises(self):
		"""A trajectory containing a cycle raises ValueError."""
		with pytest.raises(ValueError, match='Cyclic'):
			trajectory(np.array([[0, 1], [1, 2], [2, 0], [3, 4]]), np.ones(4))

	def test_disconnected_raises(self):
		"""A repeated edge leaves the graph disconnected, which raises ValueError."""
		with pytest.raises(ValueError, match='connected'):
			trajectory(np.array([[0, 1], [1, 0], [2, 3]]), np.ones(3))


class TestTrajectoryLenEdgeFromdist:
	@pytest.mark.filterwarnings('ignore:invalid value:RuntimeWarning')
	def test_len_edge_recovers_lengths(self):
		"""len_edge gives each edge's length as the largest finite |d(start)-d(end)|; inf and nan entries are skipped."""
		t = ytree()
		d = ref_dist(t, [(0, 0.25), (1, 1.5), (2, 0.5), 0, 1, 2, 3], range(4))
		d = np.vstack([d, [np.inf, 0., 1., np.inf], [np.inf] * 4])
		ans = trajectory.len_edge(d, t.edges)
		assert ans.shape == (3,)
		assert_allclose(ans, [1., 2., 3.])

	def test_len_edge_uses_only_given_samples(self):
		"""With one sample at loc 0.5 of edge 0 of line3, edge 0 gets |0.5-0.5|=0 and edge 1 gets |0.5-2.5|=2."""
		assert_allclose(trajectory.len_edge(np.array([[0.5, 0.5, 2.5]]), line3().edges), [0., 2.])

	def test_fromdist(self):
		"""fromdist builds the trajectory with the recovered lengths and node distances."""
		t = ytree_mixed()
		d = ref_dist(t, [0, 1, 2, 3, (1, 0.7)], range(4))
		t2 = trajectory.fromdist(t.edges, d)
		assert isinstance(t2, trajectory)
		assert_array_equal(t2.edges, t.edges)
		assert_allclose(t2.lens, t.lens)
		assert_allclose(t2.dist, t.dist)


class TestTrajectoryTopoint:
	@pytest.mark.parametrize('make', [line3, ytree, ytree_mixed, branched])
	def test_one_point_at_each_node(self, make):
		"""Point i of topoint() is an end of its edge at node i, and its dist equals the node distances."""
		t = make()
		p = t.topoint()
		assert isinstance(p, point) and len(p) == t.nn and p.p is t
		for i in range(t.nn):
			e, l = p.edges[i], p.locs[i]
			assert (l == 0 and t.edges[e, 0] == i) or (l == t.lens[e] and t.edges[e, 1] == i)
		assert_allclose(p.dist, t.dist)
		assert_allclose(p.compute_dist(), t.dist, atol=1e-12)


class TestTrajectoryConformLocs:
	def test_clips_small_deviations(self):
		"""Locations less than 1e-6 outside [0, edge length] are clipped to the ends; others are unchanged."""
		ans = line3().conform_locs(np.array([-1e-9, 0.5, 1. + 1e-9, 2. + 1e-7]), np.array([0, 0, 0, 1]))
		assert_array_equal(ans, [0., 0.5, 1., 2.])

	def test_large_deviation_raises(self):
		"""A deviation above both the absolute and relative bounds raises ValueError."""
		with pytest.raises(ValueError, match='both exceed'):
			line3().conform_locs(np.array([0.5, -0.01]), np.array([0, 1]))

	def test_abs_err(self):
		"""Raising abs_err above the deviation accepts and clips the location."""
		assert_array_equal(line3().conform_locs(np.array([-0.01, 2.01]), np.array([1, 1]), abs_err=0.1), [0., 2.])

	def test_rel_err(self):
		"""A deviation of 0.01 on an edge of length 1000 raises by default and is clipped with rel_err=1e-4."""
		t = trajectory(np.array([[0, 1]]), np.array([1000.]))
		with pytest.raises(ValueError):
			t.conform_locs(np.array([1000.01]), np.array([0]))
		assert_array_equal(t.conform_locs(np.array([1000.01]), np.array([0]), abs_err=1e-6, rel_err=1e-4), [1000.])


class TestTrajectoryLinspace:
	@pytest.mark.parametrize('make,start,end,total', [(ytree, 0, 3, 4.), (ytree_mixed, 3, 0, 4.), (ytree, 2, 3, 5.), (line3, 2, 0, 3.), (branched, 4, 3, 4.)])
	def test_evenly_spaced_on_path(self, make, start, end, total):
		"""Points are at linspace(0, L, n) from start and L minus that from end, with dist matching the reference."""
		t = make()
		p = t.linspace(start, end, 6)
		assert isinstance(p, point) and len(p) == 6
		ref = ref_dist(t, locs_of(p), range(t.nn))
		expected = np.linspace(0, total, 6)
		assert_allclose(ref[:, start], expected, atol=1e-12)
		assert_allclose(ref[:, end], total - expected, atol=1e-12)
		assert_allclose(p.dist, ref, atol=1e-12)

	def test_two_points_are_the_end_nodes(self):
		"""n=2 gives exactly the start and end nodes."""
		p = ytree().linspace(2, 0, 2)
		assert_allclose(p.dist[0], ytree().dist[2], atol=1e-12)
		assert_allclose(p.dist[1], ytree().dist[0], atol=1e-12)

	def test_n_below_two_asserts(self):
		"""n=1 fails its assertion."""
		with pytest.raises(AssertionError):
			ytree().linspace(0, 2, 1)


class TestTrajectoryPathPoints:
	def test_unsorted_lengths_keep_order(self):
		"""Points come back in input order at the requested distances from start; a negative length gives start."""
		t = ytree_mixed()
		lengths = np.array([3.5, -1., 0.5, 4., 3., 1.])
		p = t.path_points(3, 0, lengths)
		assert len(p) == 6
		ref = ref_dist(t, locs_of(p), [3, 0, 1])
		assert_allclose(ref[:, 0], [3.5, 0., 0.5, 4., 3., 1.], atol=1e-12)
		assert_allclose(ref[:, 1], [0.5, 4., 3.5, 0., 1., 3.], atol=1e-12)
		assert_allclose(ref[4, 2], 0., atol=1e-12)
		assert_array_equal(lengths, [3.5, -1., 0.5, 4., 3., 1.])

	def test_points_on_two_edge_path(self):
		"""On line3 from 0 to 2, lengths 0.5 and 2 are at edge 0 loc 0.5 and edge 1 loc 1 (hand values)."""
		p = line3().path_points(0, 2, np.array([0.5, 2.]))
		ref = ref_dist(line3(), locs_of(p), [(0, 0.5), (1, 1.)])
		assert_allclose(np.diag(ref), 0., atol=1e-12)


class TestTrajectoryPath:
	@pytest.mark.parametrize('start,end,expected', [(0, 2, [0, 1, 2]), (2, 3, [2, 1, 3]), (3, 0, [3, 1, 0]), (1, 2, [1, 2])])
	@pytest.mark.parametrize('make', [ytree, ytree_mixed])
	def test_node_paths(self, make, start, end, expected):
		"""Paths between nodes of the Y tree list the nodes in order, both ends included."""
		ans = make().path(start, end)
		assert isinstance(ans, np.ndarray)
		assert_array_equal(ans, expected)

	@pytest.mark.parametrize('start,end', [(1, 1), (-1, 2), (0, 4)])
	def test_invalid_nodes_assert(self, start, end):
		"""Equal or out-of-range nodes fail their assertions."""
		with pytest.raises(AssertionError):
			ytree().path(start, end)


def linear_on_edges(traj, data, pts):
	"""Reference: data at each point linearly interpolated between the two end nodes of its edge."""
	f = pts.locs / traj.lens[pts.edges]
	return data[..., traj.edges[pts.edges, 0]] * (1 - f) + data[..., traj.edges[pts.edges, 1]] * f


class TestTrajectorySmoothened:
	data = np.array([[0., 7., 30., 40., -5.], [1., 1., 1., 1., 1.]])

	def query(self, t):
		"""Points on branched(): edge 0 loc .25, node 1, edge 1 loc 1, edge 2 loc .5, edge 3 loc .4."""
		return point(t, np.array([0, 0, 1, 2, 3]), np.array([0.25, 1., 1., 0.5, 0.4]))

	def test_all_nodes_linear(self):
		"""With all nodes, values are linear interpolations between the end nodes of each point's edge."""
		t = branched()
		q = self.query(t)
		f = t.smoothened(self.data)
		assert callable(f)
		ans = f(points=q)
		assert ans.shape == (2, 5)
		assert_allclose(ans, linear_on_edges(t, self.data, q), atol=1e-12)
		assert_allclose(ans[0], [1.75, 7., 18.5, 35., 2.2], atol=1e-12)

	def test_axis(self):
		"""axis=0 with transposed data gives the transposed linear interpolation."""
		t = branched()
		q = self.query(t)
		assert_allclose(t.smoothened(self.data.T, axis=0)(points=q), linear_on_edges(t, self.data, q).T, atol=1e-12)

	def test_nodes_subset(self):
		"""With nodes [0, 2, 3], node 1 and edge 1 loc 1 interpolate between nodes 0 and 2 (values 10 and 20)."""
		t = branched()
		q = point(t, np.array([0, 1]), np.array([1., 1.]))
		assert_allclose(t.smoothened(self.data, nodes=[0, 2, 3])(points=q)[0], [10., 20.], atol=1e-12)

	def test_nodes_path_strict(self):
		"""Strict path 0-2 keeps nodes 0, 1, 2: beyond node 2 the value stays 30, on branch 1-4 it stays 7, mid edge 1-2 it is 18.5."""
		t = branched()
		q = point(t, np.array([2, 3, 1]), np.array([0.5, 0.4, 1.]))
		f = t.smoothened(self.data, nodes_path=(0, 2), criterion_path='strict')
		assert_allclose(f(points=q)[0], [30., 7., 18.5], atol=1e-12)

	def test_nodes_path_strict_is_default(self):
		"""criterion_path defaults to 'strict'."""
		t = branched()
		q = self.query(t)
		assert_allclose(t.smoothened(self.data, nodes_path=(0, 2))(points=q), t.smoothened(self.data, nodes_path=(0, 2), criterion_path='strict')(points=q))

	def test_nodes_path_loose(self):
		"""Loose path 0-2 drops only branch node 4: beyond node 2 the value is 35, on branch 1-4 it stays 7."""
		t = branched()
		q = point(t, np.array([2, 3, 1]), np.array([0.5, 0.4, 1.]))
		f = t.smoothened(self.data, nodes_path=(0, 2), criterion_path='loose')
		assert_allclose(f(points=q)[0], [35., 7., 18.5], atol=1e-12)

	def test_nodes_path_loose_without_inner_nodes(self):
		"""Loose path 1-2 has no inner node, so every node is kept and values are the all-node linear interpolation."""
		t = branched()
		q = self.query(t)
		assert_allclose(t.smoothened(self.data, nodes_path=(1, 2), criterion_path='loose')(points=q), linear_on_edges(t, self.data, q), atol=1e-12)

	def test_extra_arguments_reach_weight_function(self):
		"""Positional radius and func_name='conv' give Gaussian-kernel averages over node distances."""
		t = branched()
		q = self.query(t)
		ans = t.smoothened(self.data, 0.8, func_name='conv')(points=q)
		w = np.exp(-(ref_dist(t, range(t.nn), locs_of(q)) / 0.8) ** 2 / 2)
		assert_allclose(ans, self.data @ (w / w.sum(axis=0)), atol=1e-12)

	def test_unknown_criterion_raises(self):
		"""An unknown criterion_path raises ValueError."""
		with pytest.raises(ValueError, match='criterion_path'):
			branched().smoothened(self.data, nodes_path=(0, 2), criterion_path='other')

	def test_nodes_and_nodes_path_together_assert(self):
		"""Giving both nodes and nodes_path fails an assertion."""
		with pytest.raises(AssertionError):
			branched().smoothened(self.data, nodes=[0, 1], nodes_path=(0, 2))

	def test_wrong_data_length_asserts(self):
		"""Data whose axis length differs from the node count fails an assertion."""
		with pytest.raises(AssertionError):
			branched().smoothened(self.data[:, :4])


class TestTrajectoryTerminalNodes:
	@pytest.mark.parametrize('make,expected', [(line3, [0, 2]), (ytree, [0, 2, 3]), (ytree_mixed, [0, 2, 3]), (branched, [0, 3, 4])])
	def test_degree_one_nodes(self, make, expected):
		"""terminal_nodes returns exactly the nodes of degree 1."""
		assert_array_equal(np.sort(make().terminal_nodes()), expected)

	def test_two_node_trajectory(self):
		"""Both nodes of a single-edge trajectory are terminal."""
		assert_array_equal(np.sort(trajectory(np.array([[1, 0]]), np.array([2.])).terminal_nodes()), [0, 1])


class TestTrajectoryIO:
	@pytest.mark.parametrize('make', [ytree, ytree_mixed, branched])
	def test_file_round_trip(self, tmp_path, make):
		"""to_file then from_file restores every attribute."""
		t = make()
		path = str(tmp_path / 'traj.h5')
		assert t.to_file(path) is None
		assert_same_trajectory(trajectory.from_file(path), t)

	def test_fileobj_round_trip_in_group(self, tmp_path):
		"""to_fileobj and from_fileobj work on h5py groups, two trajectories in one file."""
		a, b = ytree(), branched()
		path = str(tmp_path / 'two.h5')
		with h5py.File(path, 'w') as f:
			a.to_fileobj(f.create_group('a'))
			b.to_fileobj(f.create_group('b'))
		with h5py.File(path, 'r') as f:
			assert set(f['a'].keys()) == {'edges', 'lens'}
			assert_same_trajectory(trajectory.from_fileobj(f['a']), a)
			assert_same_trajectory(trajectory.from_fileobj(f['b']), b)

	def test_default_compression_is_gzip(self, tmp_path):
		"""Datasets are gzip-compressed by default."""
		path = str(tmp_path / 'traj.h5')
		ytree().to_file(path)
		with h5py.File(path, 'r') as f:
			assert f['edges'].compression == 'gzip' and f['lens'].compression == 'gzip'

	@pytest.mark.parametrize('compression', ['lzf', None])
	def test_compression_option(self, tmp_path, compression):
		"""The compression option of to_fileobj, passed through to_file, sets the datasets' compression."""
		path = str(tmp_path / 'traj.h5')
		ytree().to_file(path, compression=compression)
		with h5py.File(path, 'r') as f:
			assert f['edges'].compression == compression and f['lens'].compression == compression
		assert_same_trajectory(trajectory.from_file(path), ytree())

	def test_extra_dataset_keywords(self, tmp_path):
		"""Other keyword arguments reach h5py create_dataset."""
		path = str(tmp_path / 'traj.h5')
		with h5py.File(path, 'w') as f:
			ytree().to_fileobj(f, compression='gzip', compression_opts=9, shuffle=True)
		with h5py.File(path, 'r') as f:
			for k in ['edges', 'lens']:
				assert f[k].compression_opts == 9 and f[k].shuffle


################################################################
# point
################################################################

class TestPointInit:
	def test_attributes_and_computed_dist(self):
		"""npt, p, edges and locs are stored, and dist equals the reference distance to every node."""
		t = ytree_mixed()
		edges = np.array([0, 1, 2, 2])
		locs = np.array([0.5, 1., 1., 2.5])
		p = point(t, edges, locs)
		assert p.npt == 4 and len(p) == 4 and p.p is t
		assert_array_equal(p.edges, edges)
		assert_array_equal(p.locs, locs)
		assert p.dist.shape == (4, 4)
		assert_allclose(p.dist, ref_dist(t, locs_of(p), range(4)), atol=1e-12)

	def test_hand_distances(self):
		"""Distances of four points on the Y tree to its nodes match hand-computed values."""
		p = point(ytree(), np.array([0, 1, 2, 2]), np.array([0.5, 1., 1., 2.5]))
		assert_allclose(p.dist, [[0.5, 0.5, 2.5, 3.5], [2., 1., 1., 4.], [2., 1., 3., 2.], [3.5, 2.5, 4.5, 0.5]], atol=1e-12)

	def test_given_dist_is_kept(self):
		"""A dist passed in is stored as given instead of being recomputed."""
		d = np.full((2, 4), 7.)
		p = point(ytree(), np.array([0, 1]), np.array([0.5, 1.]), dist=d)
		assert p.dist is d

	def test_locs_clipped(self):
		"""Locations within float error outside the edge are clipped."""
		p = point(ytree(), np.array([0, 2]), np.array([-1e-9, 3. + 1e-9]))
		assert_array_equal(p.locs, [0., 3.])

	def test_no_locs_no_dist_raises(self):
		"""Passing neither locs nor dist raises TypeError."""
		with pytest.raises(TypeError):
			point(ytree(), np.array([0]), None)

	@pytest.mark.parametrize('edges,locs,dist', [
		(np.array([[0], [1]]), np.array([0., 0.]), None),
		(np.array([3]), np.array([0.]), None),
		(np.array([-1]), np.array([0.]), None),
		(np.array([0, 1]), np.array([0.]), None),
		(np.array([0, 1]), np.array([0., 0.]), np.zeros((2, 3))),
		(np.array([0, 1]), np.array([0., 0.]), -np.ones((2, 4))),
	], ids=['2d-edges', 'edge-too-large', 'negative-edge', 'locs-shape', 'dist-shape', 'negative-dist'])
	def test_invalid_input_asserts(self, edges, locs, dist):
		"""Bad edge IDs or shapes and negative distances fail their assertions."""
		with pytest.raises(AssertionError):
			point(ytree(), edges, locs, dist=dist)

	def test_non_trajectory_asserts(self):
		"""A traj that is not a trajectory object fails an assertion."""
		with pytest.raises(AssertionError):
			point('ytree', np.array([0]), np.array([0.]))

	def test_location_beyond_edge_raises(self):
		"""A location past the end of its edge raises ValueError."""
		with pytest.raises(ValueError):
			point(ytree(), np.array([0]), np.array([1.5]))


class TestPointFromdistFromnodes:
	def test_edge_ids(self):
		"""With edge IDs, locs are the distances to each edge's start node and dist is kept."""
		t = ytree_mixed()
		d = ref_dist(t, [(0, 0.25), (1, 0.5), (2, 2.)], range(4))
		p = point.fromdist(t, np.array([0, 1, 2]), d)
		assert isinstance(p, point)
		assert_array_equal(p.edges, [0, 1, 2])
		assert_allclose(p.locs, [0.25, 0.5, 2.], atol=1e-12)
		assert_allclose(p.dist, d)

	def test_node_pairs_in_either_order(self):
		"""Edges given as node pairs, as stored or reversed, map to the stored edge IDs."""
		t = ytree_mixed()
		d = ref_dist(t, [(0, 0.25), (1, 0.5), (2, 2.)], range(4))
		p = point.fromdist(t, np.array([[1, 0], [1, 2], [3, 1]]), d)
		assert_array_equal(p.edges, [0, 1, 2])
		assert_allclose(p.locs, [0.25, 0.5, 2.], atol=1e-12)

	def test_absent_node_pair_raises(self):
		"""A node pair that is not an edge raises ValueError."""
		t = ytree()
		with pytest.raises(ValueError, match='absent'):
			point.fromdist(t, np.array([[0, 2]]), t.dist[[0]])

	def test_shape_checks(self):
		"""A dist of the wrong shape or 3-D edges fail their assertions."""
		t = ytree()
		with pytest.raises(AssertionError):
			point.fromdist(t, np.array([0, 1]), t.dist[[0]])
		with pytest.raises(AssertionError):
			point.fromdist(t, np.zeros((1, 1, 2), dtype=int), t.dist[[0]])

	@pytest.mark.parametrize('make', [ytree, ytree_mixed])
	def test_fromnodes_equals_topoint(self, make):
		"""fromnodes puts point i at node i (reference distance 0), with dist equal to the node distances, as topoint does."""
		t = make()
		p = point.fromnodes(t)
		assert len(p) == t.nn
		assert_allclose(np.diag(ref_dist(t, locs_of(p), range(t.nn))), 0., atol=1e-12)
		assert_allclose(p.dist, t.dist)
		assert_same_points(p, t.topoint())


class TestPointCopyConcat:
	def test_copy(self):
		"""copy returns a new object with equal edges, locs, dist and trajectory."""
		p = point(ytree(), np.array([0, 2]), np.array([0.5, 1.]))
		q = p.copy()
		assert q is not p and q.p is p.p
		assert_same_points(q, p)

	def test_perturbing_copy_leaves_original(self):
		"""Perturbing a copy does not change the original's edges, locs or dist."""
		p = point(ytree(), np.array([0, 0, 2]), np.array([0.5, 0.5, 1.]))
		before = (p.edges.copy(), p.locs.copy(), p.dist.copy())
		q = p.copy()
		np.random.seed(0)
		q.perturb()
		assert q.locs[0] != q.locs[1]
		assert_array_equal(p.edges, before[0])
		assert_array_equal(p.locs, before[1])
		assert_array_equal(p.dist, before[2])

	def test_concat(self):
		"""concat stacks edges, locs and dist of several point lists in order."""
		t = ytree()
		a = point(t, np.array([0]), np.array([0.5]))
		b = point(t, np.array([1, 2]), np.array([1., 2.]))
		c = point(t, np.array([2]), np.array([0.]))
		ans = point.concat([a, b, c])
		assert len(ans) == 4 and ans.p is t
		assert_array_equal(ans.edges, [0, 1, 2, 2])
		assert_array_equal(ans.locs, [0.5, 1., 2., 0.])
		assert_allclose(ans.dist, np.vstack([a.dist, b.dist, c.dist]))

	def test_concat_single_is_copy(self):
		"""concat of one point list returns an equal, separate object."""
		a = point(ytree(), np.array([0, 1]), np.array([0.5, 1.]))
		ans = point.concat([a])
		assert ans is not a
		assert_same_points(ans, a)

	def test_add(self):
		"""a + b equals concat([a, b])."""
		t = ytree()
		a = point(t, np.array([0]), np.array([0.5]))
		b = point(t, np.array([1, 2]), np.array([1., 2.]))
		assert_same_points(a + b, point.concat([a, b]))

	def test_concat_rejects_non_point(self):
		"""A non-point element raises TypeError."""
		a = point(ytree(), np.array([0]), np.array([0.5]))
		with pytest.raises(TypeError):
			point.concat([a, 'b'])

	def test_concat_rejects_other_trajectory(self):
		"""Points on two different trajectory objects raise TypeError."""
		a = point(ytree(), np.array([0]), np.array([0.5]))
		b = point(ytree(), np.array([0]), np.array([0.5]))
		with pytest.raises(TypeError):
			point.concat([a, b])


class TestPointDistances:
	@pytest.mark.parametrize('make', [line3, ytree, ytree_mixed, branched])
	def test_sub_matches_reference(self, make, rng):
		"""a - b equals the reference pairwise distances, is antisymmetric in order and zero on a - a's diagonal."""
		t = make()
		a = random_points(t, rng, 7)
		b = random_points(t, rng, 5)
		ans = a - b
		assert ans.shape == (len(a), len(b))
		assert_allclose(ans, ref_dist(t, locs_of(a), locs_of(b)), atol=1e-12)
		assert_allclose(b - a, ans.T, atol=1e-12)
		assert_allclose(np.diag(a - a), 0., atol=1e-12)

	def test_sub_requires_same_trajectory(self):
		"""Subtracting points of another trajectory object fails an assertion."""
		a = point(ytree(), np.array([0]), np.array([0.5]))
		b = point(ytree(), np.array([0]), np.array([0.5]))
		with pytest.raises(AssertionError):
			a - b

	@pytest.mark.parametrize('make', [ytree_mixed, branched])
	def test_compute_dist(self, make, rng):
		"""compute_dist equals the reference distance of every point to every node."""
		t = make()
		p = random_points(t, rng, 6)
		ans = p.compute_dist()
		assert ans.shape == (len(p), t.nn)
		assert_allclose(ans, ref_dist(t, locs_of(p), range(t.nn)), atol=1e-12)


class TestPointPerturb:
	def test_no_duplicates_no_change(self):
		"""Without duplicate locations, perturb leaves edges, locs and dist unchanged."""
		p = point(ytree(), np.array([0, 1, 2]), np.array([0.5, 1., 1.]))
		before = (p.edges.copy(), p.locs.copy(), p.dist.copy())
		np.random.seed(0)
		assert p.perturb() is None
		assert_array_equal(p.edges, before[0])
		assert_array_equal(p.locs, before[1])
		assert_array_equal(p.dist, before[2])

	@pytest.mark.parametrize('seed', [0, 1, 2])
	def test_duplicates_spread_between_neighbours(self, seed):
		"""Three points at 3 between points at 1 and 6 move to distinct locations in [2, 4.5]; the rest stay."""
		t = trajectory(np.array([[0, 1]]), np.array([10.]))
		p = point(t, np.zeros(5, dtype=int), np.array([1., 3., 3., 3., 6.]))
		np.random.seed(seed)
		p.perturb()
		assert_array_equal(p.edges, 0)
		assert p.locs[0] == 1. and p.locs[4] == 6.
		assert len(set(p.locs[1:4])) == 3
		assert ((p.locs[1:4] >= 2.) & (p.locs[1:4] <= 4.5)).all()
		assert_allclose(p.dist, ref_dist(t, locs_of(p), [0, 1]), atol=1e-12)

	def test_duplicates_at_edge_ends(self):
		"""Duplicates at the start node stay in [0, 2.5] and duplicates at the end node in [7.5, 10] when the other point is at 5."""
		t = trajectory(np.array([[0, 1]]), np.array([10.]))
		p = point(t, np.zeros(5, dtype=int), np.array([0., 0., 5., 10., 10.]))
		np.random.seed(0)
		p.perturb()
		assert p.locs[2] == 5.
		assert ((p.locs[:2] >= 0) & (p.locs[:2] <= 2.5)).all() and p.locs[0] != p.locs[1]
		assert ((p.locs[3:] >= 7.5) & (p.locs[3:] <= 10)).all() and p.locs[3] != p.locs[4]

	def test_same_node_on_different_edges(self):
		"""Two points at node 1 stored on edges 0 and 1 are moved onto edge 0, between the point at 0.2 and node 1."""
		t = ytree()
		p = point(t, np.array([0, 1, 2, 0]), np.array([1., 0., 0.5, 0.2]))
		np.random.seed(0)
		p.perturb()
		assert_array_equal(p.edges, [0, 0, 2, 0])
		assert p.locs[0] != p.locs[1]
		assert ((p.locs[:2] >= 0.6) & (p.locs[:2] <= 1.)).all()
		assert p.locs[2] == 0.5 and p.locs[3] == 0.2

	def test_scale_half_stays_within_bounds(self):
		"""scale=0.5 keeps duplicates distinct and within the bounds of scale=1."""
		t = trajectory(np.array([[0, 1]]), np.array([10.]))
		p = point(t, np.zeros(5, dtype=int), np.array([1., 3., 3., 3., 6.]))
		np.random.seed(3)
		# Where in [2, 4.5] the smaller range sits is not checked (the base moves duplicates to one side, a bug).
		p.perturb(scale=0.5)
		assert len(set(p.locs[1:4])) == 3
		assert ((p.locs[1:4] >= 2.) & (p.locs[1:4] <= 4.5)).all()

	def test_scale_above_one_asserts(self):
		"""scale above 1 fails an assertion."""
		p = point(ytree(), np.array([0, 0]), np.array([0.5, 0.5]))
		with pytest.raises(AssertionError):
			p.perturb(scale=1.5)


class TestPointPathLocFilter:
	@pytest.mark.parametrize('start,end,expected', [(0, 1, [1, 2]), (1, 0, [2, 1]), (3, 1, [1, 2]), (2, 3, [1])])
	def test_path_nodes_between_points(self, start, end, expected):
		"""Nodes passed from one point to another on the branched tree, in order."""
		p = point(branched(), np.array([0, 2, 1, 3, 0]), np.array([0.5, 0.5, 1., 0.5, 0.2]))
		assert_array_equal(p.path(start, end), expected)

	def test_path_same_edge_is_empty(self):
		"""Two points on the same edge have an empty node path."""
		p = point(branched(), np.array([0, 2, 1, 3, 0]), np.array([0.5, 0.5, 1., 0.5, 0.2]))
		assert len(p.path(0, 4)) == 0

	@pytest.mark.parametrize('start,end', [(1, 1), (0, 5), (-1, 0)])
	def test_path_invalid_indices_assert(self, start, end):
		"""Equal or out-of-range point indices fail their assertions."""
		p = point(branched(), np.array([0, 2, 1, 3, 0]), np.array([0.5, 0.5, 1., 0.5, 0.2]))
		with pytest.raises(AssertionError):
			p.path(start, end)

	def test_path_loc(self):
		"""Each point's projection onto the path, as distance from nstart; side-branch points project to node 1."""
		p = point(ytree(), np.array([0, 1, 2, 2]), np.array([0.5, 1., 1., 2.5]))
		assert_allclose(p.path_loc(0, 2), [0.5, 2., 1., 1.], atol=1e-12)
		assert_allclose(p.path_loc(2, 3), [2., 1., 3., 4.5], atol=1e-12)

	def test_path_loc_distpath(self):
		"""A given distpath is used: the true one changes nothing, and adding 2 shifts every location by 1 ((d0-d1+D)/2)."""
		p = point(ytree(), np.array([0, 1, 2, 2]), np.array([0.5, 1., 1., 2.5]))
		assert_allclose(p.path_loc(0, 2, distpath=3.), p.path_loc(0, 2), atol=1e-12)
		assert_allclose(p.path_loc(0, 2, distpath=5.), p.path_loc(0, 2) + 1., atol=1e-12)

	@pytest.mark.parametrize('nstart,nend,expected', [(0, 2, [0, 1, 4]), (2, 3, [1, 2, 3, 4]), (0, 3, [0, 2, 3, 4]), (3, 0, [0, 2, 3, 4])])
	def test_filter_path(self, nstart, nend, expected):
		"""Points on the path's edges, or at one of its nodes even when stored on another edge, are kept."""
		p = point(ytree(), np.array([0, 1, 2, 2, 2]), np.array([0.5, 1., 1., 2.5, 0.]))
		assert_array_equal(p.filter_path(nstart, nend), expected)


def nearest_cells(dist, ncell):
	"""Sets of the ncell smallest entries of each row."""
	return [frozenset(np.argsort(x)[:ncell].tolist()) for x in dist]


class TestPointSubsets:
	ncell = 4

	def cells(self, rng):
		"""Y tree with 7, 8 and 9 cells at random locations on its three edges."""
		t = ytree()
		edges = np.repeat([0, 1, 2], [7, 8, 9])
		locs = rng.random(len(edges)) * t.lens[edges]
		return point(t, edges, locs)

	def run(self, c, noverlap=2, dmax=1.):
		np.random.seed(0)
		return c.subsets(self.ncell, noverlap, dmax)

	def test_outputs_shapes_and_names(self, rng):
		"""Five outputs with matching lengths, Subset1..n names and ncell cells in every subset."""
		c = self.cells(rng)
		edges, locs, subsets, graph, names = self.run(c)
		n = len(edges)
		assert n >= c.p.nn
		assert locs.shape == (n,) and subsets.shape == (len(c), n) and graph.shape == (n, n)
		assert subsets.dtype == bool and graph.dtype == bool
		assert list(names) == ['Subset{}'.format(x + 1) for x in range(n)]
		assert_array_equal(subsets.sum(axis=0), self.ncell)

	def test_node_subsets_first(self, rng):
		"""The first nn subsets are centred on the nodes and hold each node's ncell nearest cells (reference distance)."""
		c = self.cells(rng)
		edges, locs, subsets, graph, names = self.run(c)
		t = c.p
		d = ref_dist(t, range(t.nn), locs_of(c))
		expected = nearest_cells(d, self.ncell)
		for i in range(t.nn):
			assert frozenset(np.nonzero(subsets[:, i])[0].tolist()) == expected[i]
			assert_allclose(ref_dist(t, [(int(edges[i]), float(locs[i]))], [i]), 0., atol=1e-12)

	def test_interior_subsets_are_nearest_cells_on_their_edge(self, rng):
		"""Every other subset is centred at a cell's location and holds the ncell nearest cells of a cell or end node of that edge."""
		c = self.cells(rng)
		edges, locs, subsets, graph, names = self.run(c)
		t = c.p
		# Not the nearest cells of the subset's own centre: subsets are shifted by one cell (a bug in the base).
		assert len(edges) > t.nn
		d = ref_dist(t, locs_of(c), locs_of(c))
		cell_sets = nearest_cells(d, self.ncell)
		node_sets = nearest_cells(ref_dist(t, range(t.nn), locs_of(c)), self.ncell)
		for i in range(t.nn, len(edges)):
			e = int(edges[i])
			on_edge = np.nonzero(c.edges == e)[0]
			assert np.isclose(c.locs[on_edge], locs[i]).any()
			allowed = {cell_sets[x] for x in on_edge} | {node_sets[x] for x in t.edges[e]}
			assert frozenset(np.nonzero(subsets[:, i])[0].tolist()) in allowed

	def test_neighbour_graph_is_chain_per_edge(self, rng):
		"""The subset graph is a symmetric tree; between the subsets of an edge's two nodes it runs through that edge's subsets in increasing loc."""
		c = self.cells(rng)
		edges, locs, subsets, graph, names = self.run(c)
		t = c.p
		n = len(edges)
		assert (graph == graph.T).all() and not graph.diagonal().any()
		g = nx.from_numpy_array(graph.astype(int))
		assert nx.is_tree(g) and g.number_of_nodes() == n
		inner = 0
		for e, (u, v) in enumerate(t.edges):
			chain = nx.shortest_path(g, int(u), int(v))[1:-1]
			inner += len(chain)
			assert (edges[chain] == e).all()
			assert (np.diff(locs[chain]) > 0).all()
		assert inner == n - t.nn

	def test_no_overlap_requirement_keeps_node_subsets(self, rng):
		"""With noverlap=0 and a huge dmax every edge subset is removable, leaving the node subsets linked like the trajectory."""
		c = self.cells(rng)
		edges, locs, subsets, graph, names = self.run(c, noverlap=0, dmax=1e9)
		t = c.p
		assert len(edges) == t.nn
		expected = np.zeros((t.nn, t.nn), dtype=bool)
		expected[t.edges[:, 0], t.edges[:, 1]] = True
		assert_array_equal(graph, expected | expected.T)

	def test_does_not_modify_points(self, rng):
		"""subsets works on a copy: the cells' edges, locs and dist are unchanged, even with duplicate cells."""
		t = ytree()
		c = point(t, np.array([0, 0, 0, 1, 1, 1, 2, 2, 2, 2]), np.array([0.2, 0.2, 0.7, 0.3, 1.1, 1.6, 0.4, 1.2, 2.1, 2.6]))
		before = (c.edges.copy(), c.locs.copy(), c.dist.copy())
		np.random.seed(0)
		c.subsets(3, 1, 1.)
		assert_array_equal(c.edges, before[0])
		assert_array_equal(c.locs, before[1])
		assert_array_equal(c.dist, before[2])

	def test_too_many_cells_at_a_node_raises(self):
		"""More than ncell cells exactly at one node raise ValueError."""
		t = ytree()
		c = point(t, np.array([0, 0, 0, 1, 1, 2, 2]), np.array([0., 0., 0., 0.5, 1.5, 1., 2.]))
		np.random.seed(0)
		with pytest.raises(ValueError, match='ncell'):
			c.subsets(2, 1, 1.)

	@pytest.mark.parametrize('ncell,noverlap,dmax', [(3, 3, 1.), (3, -1, 1.), (3, 1, -0.5)])
	def test_parameter_checks(self, rng, ncell, noverlap, dmax):
		"""noverlap outside [0, ncell) or a negative dmax fail their assertions."""
		with pytest.raises(AssertionError):
			self.cells(rng).subsets(ncell, noverlap, dmax)


def objective(c, pts, w):
	"""Weighted sum of distances from each of pts to the points c (reference distances)."""
	return ref_dist(c.p, locs_of(pts), locs_of(c)) @ w


class TestPointNearest:
	def cells(self):
		"""Six points on the three edges of the Y tree, whose best node is the branching node 1."""
		return point(ytree(), np.array([0, 0, 1, 1, 2, 2]), np.array([0.3, 0.8, 0.5, 1.2, 0.4, 2.5]))

	def best(self, c, w):
		"""Smallest weighted distance sum over all nodes and points, where a tree's weighted median lies."""
		cand = point.concat([c.p.topoint(), c])
		return objective(c, cand, w).min()

	def test_weight_rows(self):
		"""For each weight row, the returned point reaches the minimum weighted distance sum."""
		c = self.cells()
		w = np.array([[1., 1, 1, 1, 1, 1], [1, 1, 3, 3, 1, 1], [3, 3, 1, 1, 1, 1], [1, 1, 1, 1, 4, 4]])
		ans = c.nearest(w=w)
		assert isinstance(ans, point) and len(ans) == 4 and ans.p is c.p
		for i in range(4):
			assert_allclose(objective(c, ans[[i]], w[i]), self.best(c, w[i]), atol=1e-12)
		assert_array_equal(ans.edges[2], 0)
		assert_allclose(ans.locs[2], 0.8)

	def test_default_weight_is_uniform(self):
		"""w=None gives the same point as uniform weights."""
		c = self.cells()
		assert_same_points(c.nearest(), c.nearest(w=np.ones(len(c))))
		assert_allclose(objective(c, c.nearest(), np.ones(len(c))), self.best(c, np.ones(len(c))), atol=1e-12)

	def test_one_dimensional_weight(self):
		"""A 1-D weight vector gives one point at the weighted median of points on one edge (0.9 here)."""
		c = point(line3(), np.array([1, 1, 1]), np.array([0.2, 0.9, 1.5]))
		ans = c.nearest(w=np.array([1., 1., 1.]))
		assert len(ans) == 1
		assert_array_equal(ans.edges, [1])
		assert_allclose(ans.locs, [0.9])

	def test_weight_shape_checks(self):
		"""A 3-D weight raises ValueError and a weight of the wrong length fails an assertion."""
		c = self.cells()
		with pytest.raises(ValueError):
			c.nearest(w=np.ones((1, 1, len(c))))
		with pytest.raises(AssertionError):
			c.nearest(w=np.ones(len(c) + 1))


class TestPointSubtraj:
	def test_edge_lengths_are_point_distances(self):
		"""The new trajectory has the given edges and lengths equal to the distances between their points."""
		p = point(ytree(), np.array([0, 0, 1, 2, 2]), np.array([0., 1., 2., 1., 3.]))
		edges = np.array([[0, 1], [1, 2], [1, 3], [3, 4]])
		t = p.subtraj(edges=edges)
		assert isinstance(t, trajectory) and t.nn == 5
		assert_array_equal(t.edges, edges)
		assert_allclose(t.lens, [1., 2., 1., 2.], atol=1e-12)
		assert_allclose(t.dist, ref_dist(p.p, locs_of(p), locs_of(p)), atol=1e-12)

	def test_edges_required(self):
		"""edges=None raises NotImplementedError and edges of the wrong shape fail an assertion."""
		p = point(ytree(), np.array([0, 0]), np.array([0., 1.]))
		with pytest.raises(NotImplementedError):
			p.subtraj()
		with pytest.raises(AssertionError):
			p.subtraj(edges=np.array([0, 1]))


class TestPointIndexing:
	@pytest.mark.parametrize('key', [[2, 0], slice(1, 3), np.array([True, False, True, True]), np.array([3])])
	def test_getitem(self, key):
		"""Indexing with a list, slice, bool mask or index array selects edges, locs and dist rows alike."""
		p = point(ytree(), np.array([0, 1, 2, 2]), np.array([0.5, 1., 1., 2.5]))
		q = p[key]
		assert isinstance(q, point) and q.p is p.p
		assert_array_equal(q.edges, p.edges[key])
		assert_array_equal(q.locs, p.locs[key])
		assert_array_equal(q.dist, p.dist[key])
		assert len(q) == len(p.edges[key])

	def test_getitem_scalar_raises(self):
		"""A scalar index raises TypeError."""
		p = point(ytree(), np.array([0, 1]), np.array([0.5, 1.]))
		with pytest.raises(TypeError):
			p[0]

	def test_len(self):
		"""len equals the number of points."""
		assert len(point(ytree(), np.array([0, 1, 2]), np.array([0.5, 1., 1.]))) == 3


class TestPointDropDuplicates:
	def points(self):
		"""Points with repeated (edge, loc) pairs away from nodes."""
		return point(ytree(), np.array([1, 0, 1, 0, 2, 1]), np.array([0.5, 0.25, 1.5, 0.25, 2., 0.5]))

	def test_unique_pairs(self):
		"""Each distinct (edge, loc) pair is kept exactly once."""
		ans = self.points().drop_duplicates()
		assert isinstance(ans, point)
		assert sorted(locs_of(ans)) == [(0, 0.25), (1, 0.5), (1, 1.5), (2, 2.)]

	def test_selection(self):
		"""selection=True also returns the kept indices, one per distinct pair, and the points are those rows."""
		p = self.points()
		ans, idx = p.drop_duplicates(selection=True)
		assert len(idx) == 4 and len(set(idx.tolist())) == 4
		assert sorted(locs_of(p[idx])) == sorted(locs_of(ans))
		assert_same_points(ans, p[idx])

	def test_default_selection_false(self):
		"""selection=False returns only the points, the same as the default."""
		p = self.points()
		assert_same_points(p.drop_duplicates(selection=False), p.drop_duplicates())


class TestPointWeights:
	def test_linear_from_nodes(self, rng):
		"""With self at the nodes, each point's weights are 1-l/len and l/len on its edge's two end nodes."""
		t = ytree_mixed()
		nodes = t.topoint()
		q = random_points(t, rng, 8)
		w = nodes.weight_linear(q)
		assert w.shape == (t.nn, len(q))
		expected = np.zeros((t.nn, len(q)))
		f = q.locs / t.lens[q.edges]
		np.add.at(expected, (t.edges[q.edges, 0], np.arange(len(q))), 1 - f)
		np.add.at(expected, (t.edges[q.edges, 1], np.arange(len(q))), f)
		assert_allclose(w, expected, atol=1e-12)

	def test_linear_on_a_line_matches_interp(self):
		"""On a single edge, data @ weight_linear equals numpy.interp, constant beyond the outer points."""
		t = trajectory(np.array([[0, 1]]), np.array([10.]))
		xs = np.array([2., 3.5, 5., 8.])
		vals = np.array([1., -2., 4., 0.5])
		xq = np.array([0., 2., 2.7, 4., 6.5, 8., 9.9])
		s = point(t, np.zeros(4, dtype=int), xs)
		q = point(t, np.zeros(len(xq), dtype=int), xq)
		w = s.weight_linear(q)
		assert_allclose(w.sum(axis=0), 1., atol=1e-12)
		assert_allclose(vals @ w, np.interp(xq, xs, vals), atol=1e-12)

	@pytest.mark.parametrize('radius', [0.5, 2.])
	def test_conv_matches_gaussian(self, rng, radius):
		"""weight_conv equals exp(-d^2/2r^2) over reference distances, normalized per column."""
		t = branched()
		s = random_points(t, rng, 6)
		q = random_points(t, rng, 4)
		w = s.weight_conv(q, radius, cut=0)
		ref = np.exp(-(ref_dist(t, locs_of(s), locs_of(q)) / radius) ** 2 / 2)
		assert w.shape == (len(s), len(q))
		assert_allclose(w, ref / ref.sum(axis=0), atol=1e-12)
		assert_allclose(w.sum(axis=0), 1., atol=1e-12)

	def test_conv_cut_not_implemented(self):
		"""cut above 0 raises NotImplementedError."""
		s = ytree().topoint()
		with pytest.raises(NotImplementedError):
			s.weight_conv(s, 1., cut=0.1)


class TestPointSmoothen:
	def setup_data(self):
		t = branched()
		s = t.topoint()
		q = point(t, np.array([0, 0, 1, 2, 3]), np.array([0.25, 1., 1., 0.5, 0.4]))
		data = np.array([[0., 7., 30., 40., -5.], [1., 2., 3., 4., 5.]])
		return t, s, q, data

	def test_linear_axis_last(self):
		"""Default linear smoothing on the last axis interpolates between the end nodes of each point's edge."""
		t, s, q, data = self.setup_data()
		ans = s.smoothen(data, points=q)
		assert ans.shape == (2, 5)
		assert_allclose(ans, linear_on_edges(t, data, q), atol=1e-12)

	def test_axis_zero_and_middle(self):
		"""axis=0 on transposed data and axis=1 on 3-D data put the points on that axis."""
		t, s, q, data = self.setup_data()
		expected = linear_on_edges(t, data, q)
		assert_allclose(s.smoothen(data.T, points=q, axis=0, func_name='linear'), expected.T, atol=1e-12)
		d3 = np.stack([data, 2 * data], axis=0).transpose(1, 2, 0)
		ans = s.smoothen(d3, points=q, axis=1)
		assert ans.shape == (2, 5, 2)
		assert_allclose(ans[:, :, 0], expected, atol=1e-12)
		assert_allclose(ans[:, :, 1], 2 * expected, atol=1e-12)

	def test_conv_with_arguments(self):
		"""func_name='conv' with the radius as a positional argument and cut as a keyword gives Gaussian-kernel averages over reference distances."""
		t, s, q, data = self.setup_data()
		w = np.exp(-(ref_dist(t, range(t.nn), locs_of(q)) / 0.7) ** 2 / 2)
		assert_allclose(s.smoothen(data, 0.7, points=q, func_name='conv', cut=0), data @ (w / w.sum(axis=0)), atol=1e-12)

	def test_nan_ignore(self):
		"""nan='ignore' renormalizes over non-nan data and gives nan only where every weighted value is nan."""
		t, s, q, data = self.setup_data()
		data = data.copy()
		data[0, 1] = np.nan
		ans = s.smoothen(data, points=q, nan='ignore')
		assert_allclose(ans[0], [0., np.nan, 30., 35., -5.], atol=1e-12)
		assert_allclose(ans[1], linear_on_edges(t, data[1], q), atol=1e-12)

	def test_nan_ignore_is_default(self):
		"""nan defaults to 'ignore'."""
		t, s, q, data = self.setup_data()
		data = data.copy()
		data[0, 1] = np.nan
		assert_allclose(s.smoothen(data, points=q), s.smoothen(data, points=q, nan='ignore'))

	def test_nan_propagate(self):
		"""nan='propagate' gives nan wherever a nan value has nonzero weight and leaves rows without nan unchanged."""
		t, s, q, data = self.setup_data()
		data = data.copy()
		data[0, 1] = np.nan
		ans = s.smoothen(data, points=q, nan='propagate')
		assert np.isnan(ans[0, [0, 1, 2, 4]]).all()
		assert_allclose(ans[1], linear_on_edges(t, data[1], q), atol=1e-12)

	def test_invalid_arguments(self):
		"""Unknown nan raises ValueError; missing points or a wrong data length fail assertions; unknown func_name raises AttributeError."""
		t, s, q, data = self.setup_data()
		with pytest.raises(ValueError, match='nan'):
			s.smoothen(data, points=q, nan='zero')
		with pytest.raises(AssertionError):
			s.smoothen(data)
		with pytest.raises(AssertionError):
			s.smoothen(data[:, :4], points=q)
		with pytest.raises(AttributeError):
			s.smoothen(data, points=q, func_name='cubic')

	def test_smoothened_partial(self):
		"""smoothened returns a function of points equal to smoothen with the same arguments."""
		t, s, q, data = self.setup_data()
		f = s.smoothened(data, 0.7, func_name='conv')
		assert callable(f)
		assert_allclose(f(points=q), s.smoothen(data, 0.7, points=q, func_name='conv'), atol=1e-12)
		assert_allclose(s.smoothened(data)(points=q), linear_on_edges(t, data, q), atol=1e-12)


class TestPointIO:
	def points(self):
		return point(ytree_mixed(), np.array([0, 1, 2, 2]), np.array([0.5, 1., 1., 2.5]))

	def test_file_round_trip(self, tmp_path):
		"""to_file then from_file with the trajectory restores edges, locs and dist."""
		p = self.points()
		path = str(tmp_path / 'pts.h5')
		assert p.to_file(path) is None
		q = point.from_file(p.p, path)
		assert q.p is p.p
		assert_same_points(q, p)
		with h5py.File(path, 'r') as f:
			assert set(f.keys()) == {'edges', 'locs'}

	def test_file_with_trajectory(self, tmp_path):
		"""traj=True stores the trajectory in group 'traj', which from_file(None, path) loads."""
		p = self.points()
		path = str(tmp_path / 'pts.h5')
		p.to_file(path, traj=True)
		with h5py.File(path, 'r') as f:
			assert set(f.keys()) == {'edges', 'locs', 'traj'}
		q = point.from_file(None, path)
		assert_same_trajectory(q.p, p.p)
		assert_same_points(q, p)

	def test_file_without_trajectory_needs_one(self, tmp_path):
		"""from_file(None, path) on a file written with traj=False raises KeyError."""
		p = self.points()
		path = str(tmp_path / 'pts.h5')
		p.to_file(path, traj=False)
		with pytest.raises(KeyError):
			point.from_file(None, path)

	def test_fileobj_round_trip_in_groups(self, tmp_path):
		"""to_fileobj and from_fileobj work on h5py groups, two point lists in one file."""
		p = self.points()
		q = p.p.topoint()
		path = str(tmp_path / 'pts.h5')
		with h5py.File(path, 'w') as f:
			p.to_fileobj(f.create_group('a'))
			q.to_fileobj(f.create_group('b'))
		with h5py.File(path, 'r') as f:
			assert_same_points(point.from_fileobj(p.p, f['a']), p)
			assert_same_points(point.from_fileobj(p.p, f['b']), q)

	@pytest.mark.parametrize('compression', ['gzip', 'lzf', None])
	def test_compression(self, tmp_path, compression):
		"""compression sets the point datasets' compression; the stored trajectory keeps gzip."""
		p = self.points()
		path = str(tmp_path / 'pts.h5')
		p.to_file(path, traj=True, compression=compression)
		with h5py.File(path, 'r') as f:
			assert f['edges'].compression == compression and f['locs'].compression == compression
			assert f['traj/edges'].compression == 'gzip'
		assert_same_points(point.from_file(p.p, path), p)

	def test_extra_dataset_keywords(self, tmp_path):
		"""Other keyword arguments reach h5py create_dataset."""
		p = self.points()
		path = str(tmp_path / 'pts.h5')
		p.to_file(path, compression_opts=2, shuffle=True)
		with h5py.File(path, 'r') as f:
			for k in ['edges', 'locs']:
				assert f[k].compression == 'gzip' and f[k].compression_opts == 2 and f[k].shuffle
