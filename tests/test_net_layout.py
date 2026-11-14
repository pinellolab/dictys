"""
Tests for dictys.net.layout: the Fruchterman-Reingold layout adapted from
networkx, checked against networkx's spring_layout and a short reference
implementation of the force and cooling rules.
"""

import networkx as nx
import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from dictys.net import layout


def small_graph():
	"""Weighted graph on five nodes: a path a-b-c-d plus a chord a-c and a pendant e on d."""
	g = nx.Graph()
	g.add_weighted_edges_from([('a', 'b', 1.0), ('b', 'c', 2.0), ('c', 'd', 0.5), ('a', 'c', 1.5), ('d', 'e', 1.0)])
	return g


def initial_positions(g, seed=0):
	rng = np.random.default_rng(seed)
	return {x: rng.uniform(0, 1, size=2) for x in g}


def reference_steps(A, pos, k, t0, dt, nsteps, fixed=()):
	"""Fruchterman-Reingold updates: repulsion k^2/d, attraction A*d/k, step capped by a linearly cooling temperature."""
	pos = pos.astype(float).copy()
	t = t0
	for _ in range(nsteps):
		delta = pos[:, None, :] - pos[None, :, :]
		dist = np.maximum(np.linalg.norm(delta, axis=-1), 0.01)
		disp = (delta * (k * k / dist ** 2 - A * dist / k)[:, :, None]).sum(axis=1)
		length = np.linalg.norm(disp, axis=-1)
		length = np.where(length < 0.01, 0.1, length)
		step = disp * (t / length)[:, None]
		step[list(fixed)] = 0
		pos += step
		t -= dt
	return pos


# fruchterman_reingold_layout


def test_layout_of_empty_graph():
	"""An empty graph has an empty layout."""
	assert layout.fruchterman_reingold_layout(nx.Graph()) == {}


def test_layout_of_single_node():
	"""A single node is placed at the center, which defaults to the origin in dim dimensions."""
	g = nx.Graph()
	g.add_node('a')
	ans = layout.fruchterman_reingold_layout(g)
	assert list(ans) == ['a']
	assert_array_equal(ans['a'], [0, 0])
	assert_array_equal(layout.fruchterman_reingold_layout(g, center=(1.5, -2))['a'], [1.5, -2])
	assert_array_equal(layout.fruchterman_reingold_layout(g, dim=3, seed=1)['a'], [0, 0, 0])
	assert_array_equal(layout.fruchterman_reingold_layout(['x'])['x'], [0, 0])


def test_layout_fixed_needs_positions():
	"""Fixing nodes without giving their positions raises ValueError."""
	g = small_graph()
	with pytest.raises(ValueError):
		layout.fruchterman_reingold_layout(g, fixed=['a'])
	with pytest.raises(ValueError):
		layout.fruchterman_reingold_layout(g, fixed=['a'], pos={'b': (0, 0)})


# _fruchterman_reingold: the layout function that dictys.net.stat.flayout_base takes


def test_kernel_matches_networkx_spring_layout():
	"""With default step sizes the kernel reproduces networkx spring_layout from the same start without rescaling."""
	g = small_graph()
	pos = initial_positions(g)
	A = nx.to_numpy_array(g)
	pos_arr = np.array([pos[x] for x in g])
	ans = layout._fruchterman_reingold(A, pos_arr, iterations=30, threshold=1e-6)
	ref = nx.spring_layout(g, pos=pos, iterations=30, threshold=1e-6, scale=None)
	assert_allclose(ans, np.array([ref[x] for x in g]), rtol=1e-8, atol=1e-10)
	assert_array_equal(pos_arr, np.array([pos[x] for x in g]))


def test_kernel_fixed_nodes_match_networkx():
	"""Fixed node indices keep their initial positions, as in networkx spring_layout with fixed nodes."""
	g = small_graph()
	pos = initial_positions(g, seed=1)
	A = nx.to_numpy_array(g)
	pos_arr = np.array([pos[x] for x in g])
	ans = layout._fruchterman_reingold(A, pos_arr, k=0.4, fixed=np.array([0, 3]), iterations=25)
	assert_array_equal(ans[[0, 3]], pos_arr[[0, 3]])
	assert not np.allclose(ans[[1, 2, 4]], pos_arr[[1, 2, 4]])
	ref = nx.spring_layout(g, k=0.4, pos=pos, fixed=['a', 'd'], iterations=25)
	assert_allclose(ans, np.array([ref[x] for x in g]), rtol=1e-8, atol=1e-10)


def test_kernel_three_dimensions():
	"""In three dimensions the kernel still agrees with networkx spring_layout."""
	g = small_graph()
	rng = np.random.default_rng(2)
	pos = {x: rng.uniform(0, 1, size=3) for x in g}
	pos_arr = np.array([pos[x] for x in g])
	ans = layout._fruchterman_reingold(nx.to_numpy_array(g), pos_arr, iterations=20, dim=3)
	assert ans.shape == (5, 3)
	ref = nx.spring_layout(g, pos=pos, iterations=20, scale=None, dim=3)
	assert_allclose(ans, np.array([ref[x] for x in g]), rtol=1e-8, atol=1e-10)


@pytest.mark.parametrize('stop', [1, 2, 5])
def test_kernel_stop_and_pretend_iterations(stop):
	"""stop ends after that many steps; step s moves every node by t0-(s-1)*dt, t0=0.1*span*pretendIterations/iterations, dt=t0/(iterations+1)."""
	g = small_graph()
	A = nx.to_numpy_array(g)
	pos_arr = np.array(list(initial_positions(g, seed=3).values()))
	k = np.sqrt(1 / 5)
	span = (pos_arr.max(axis=0) - pos_arr.min(axis=0)).max()
	t0 = span * 0.1 * 10 / 40
	ans = layout._fruchterman_reingold(A, pos_arr, iterations=40, threshold=0, pretendIterations=10, stop=stop)
	prev = layout._fruchterman_reingold(A, pos_arr, iterations=40, threshold=0, pretendIterations=10, stop=stop - 1)
	assert_allclose(np.linalg.norm(ans - prev, axis=1), t0 - (stop - 1) * t0 / 41, rtol=1e-10)
	assert_allclose(ans, reference_steps(A, pos_arr, k, t0, t0 / 41, stop), rtol=1e-10, atol=1e-12)


def test_kernel_defaults_and_determinism():
	"""stop defaults to iterations, the result is deterministic, and a large threshold stops after one step."""
	g = small_graph()
	A = nx.to_numpy_array(g)
	pos_arr = np.array(list(initial_positions(g, seed=4).values()))
	full = layout._fruchterman_reingold(A, pos_arr, iterations=15, threshold=0)
	assert full.shape == (5, 2)
	assert_array_equal(full, layout._fruchterman_reingold(A, pos_arr, iterations=15, threshold=0, stop=15))
	assert_array_equal(full, layout._fruchterman_reingold(A, pos_arr, iterations=15, threshold=0))
	early = layout._fruchterman_reingold(A, pos_arr, iterations=15, threshold=1e6)
	assert_allclose(early, layout._fruchterman_reingold(A, pos_arr, iterations=15, threshold=0, stop=1), rtol=1e-12)
	span = (pos_arr.max(axis=0) - pos_arr.min(axis=0)).max()
	assert_allclose(full, reference_steps(A, pos_arr, np.sqrt(1 / 5), 0.1 * span, 0.1 * span / 16, 15), rtol=1e-10, atol=1e-12)


def test_kernel_rejects_non_matrix():
	"""The kernel raises NetworkXError when the adjacency matrix is not an array."""
	with pytest.raises(nx.NetworkXError):
		layout._fruchterman_reingold([[0, 1], [1, 0]], np.zeros((2, 2)))


def test_sparse_kernel_agrees_with_dense():
	"""The sparse kernel gives the dense kernel's layout, with and without fixed nodes."""
	g = small_graph()
	A = nx.to_numpy_array(g)
	pos_arr = np.array(list(initial_positions(g, seed=5).values()))
	dense = layout._fruchterman_reingold(A, pos_arr, iterations=20, pretendIterations=10, stop=12)
	sparse = layout._sparse_fruchterman_reingold(A, pos_arr, iterations=20, pretendIterations=10, stop=12)
	assert_allclose(sparse, dense, rtol=1e-10, atol=1e-12)
	dense = layout._fruchterman_reingold(A, pos_arr, k=0.3, fixed=np.array([1]), iterations=10)
	sparse = layout._sparse_fruchterman_reingold(A, pos_arr, k=0.3, fixed=[1], iterations=10)
	assert_allclose(sparse, dense, rtol=1e-10, atol=1e-12)
	assert_array_equal(sparse[1], pos_arr[1])
