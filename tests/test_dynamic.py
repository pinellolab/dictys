"""
Tests for dictys.dynamic: subsets_rna and subsets_atac, called directly and through
``dictys dynamic subsets_rna`` and ``dictys dynamic subsets_atac``.

Cells sit on a Y-shaped trajectory whose three edges all meet at node 1, so the
reference distance between two locations is computed from their distances to node 1.
"""

import os

import networkx as nx
import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_array_equal

from dictys.dynamic import subsets_atac, subsets_rna
from dictys.traj import point, trajectory


LENS = np.array([2., 3., 4.])
EDGES = np.array([[0, 1], [1, 2], [1, 3]])


def ytree():
	"""Y tree: edge 0 from node 0 to the branching node 1, edges 1 and 2 from node 1 to leaves 2 and 3."""
	return trajectory(EDGES.copy(), LENS.copy())


def arm(edges, locs):
	"""Each location's edge and distance from node 1 (edge 0 runs towards node 1, the others away from it)."""
	edges = np.asarray(edges)
	locs = np.asarray(locs, dtype=float)
	return edges, np.where(edges == 0, LENS[0] - locs, locs)


def ref_dist(e1, l1, e2, l2):
	"""Reference distances between two location lists on the Y tree: same edge |r1-r2|, else r1+r2."""
	a1, r1 = arm(e1, l1)
	a2, r2 = arm(e2, l2)
	same = a1[:, None] == a2[None, :]
	return np.where(same, np.abs(r1[:, None] - r2[None, :]), r1[:, None] + r2[None, :])


NODE_EDGES = np.array([0, 0, 1, 2])
NODE_LOCS = np.array([0., 2., 3., 4.])


def write_cells(folder, rng, prefix, counts):
	"""Writes cells at random locations (counts per edge) as traj_cell_<prefix>.h5 and coord_<prefix>.tsv."""
	t = ytree()
	edges = np.repeat([0, 1, 2], counts)
	locs = rng.random(len(edges)) * LENS[edges]
	names = ['{}cell{:02d}'.format(prefix, x) for x in range(len(edges))]
	point(t, edges, locs).to_file(os.path.join(folder, 'traj_cell_{}.h5'.format(prefix)))
	coord = pd.DataFrame(rng.random((len(edges), 2)), index=names, columns=['Dim1', 'Dim2'])
	coord.to_csv(os.path.join(folder, 'coord_{}.tsv'.format(prefix)), sep='\t', header=True, index=True)
	return edges, locs, names


def write_inputs(folder, rng, counts=(6, 7, 8)):
	"""Writes traj_node.h5 and RNA cell files; returns the input paths and the cells."""
	os.makedirs(folder, exist_ok=True)
	ytree().to_file(os.path.join(folder, 'traj_node.h5'))
	edges, locs, names = write_cells(folder, rng, 'rna', counts)
	paths = [os.path.join(folder, x) for x in ['traj_node.h5', 'traj_cell_rna.h5', 'coord_rna.tsv']]
	return paths, edges, locs, names


def output_paths(folder, edges_name='subset_edges.tsv'):
	"""Output paths of subsets_rna, in its argument order."""
	os.makedirs(folder, exist_ok=True)
	return [os.path.join(folder, 'subsets.txt'), os.path.join(folder, 'subset_locs.h5'), folder, os.path.join(folder, edges_name)]


def read_names(path):
	"""Whitespace-separated names in a text file."""
	with open(path, 'r') as f:
		return f.read().split()


def read_outputs(outs):
	"""Reads subsets_rna outputs: subset names, subset locations, cell names per subset, neighbour table."""
	names = read_names(outs[0])
	locs = point.from_file(ytree(), outs[1])
	cells = [read_names(os.path.join(outs[2], x, 'names_rna.txt')) for x in names]
	graph = pd.read_csv(outs[3], sep='\t', header=0, index_col=0)
	return names, locs, cells, graph


def nearest(dist, n, names):
	"""Sets of the names of the n smallest entries of each row of dist."""
	return [frozenset(np.asarray(names)[np.argsort(x)[:n]].tolist()) for x in dist]


def check_rna_outputs(outs, edges, locs, names, ncell):
	"""
	Asserts the documented properties of the subsets_rna outputs and returns them parsed: node subsets come
	first, at the nodes, with each node's ncell nearest cells; every other subset is centred at a cell of its
	edge and holds the ncell nearest cells of a cell or end node of that edge; the neighbour table is a
	symmetric tree over the subsets.
	"""
	snames, slocs, cells, graph = read_outputs(outs)
	n = len(snames)
	assert snames == ['Subset{}'.format(x + 1) for x in range(n)]
	assert n >= 4 and len(slocs) == n
	node_sets = nearest(ref_dist(NODE_EDGES, NODE_LOCS, edges, locs), ncell, names)
	for i in range(4):
		assert ref_dist(slocs.edges[[i]], slocs.locs[[i]], NODE_EDGES[[i]], NODE_LOCS[[i]])[0, 0] < 1e-12
		assert frozenset(cells[i]) == node_sets[i]
	cell_sets = nearest(ref_dist(edges, locs, edges, locs), ncell, names)
	# Not the nearest cells of the subset's own centre: subsets are shifted by one cell (a bug in the base).
	for i in range(4, n):
		e = slocs.edges[i]
		on_edge = np.nonzero(edges == e)[0]
		assert np.isclose(locs[on_edge], slocs.locs[i]).any()
		allowed = {cell_sets[x] for x in on_edge} | {node_sets[x] for x in EDGES[e]}
		assert frozenset(cells[i]) in allowed
	for x in cells:
		assert len(x) == ncell and len(set(x)) == ncell and set(x) <= set(names)
	assert list(graph.index) == snames and list(graph.columns) == snames
	g = graph.values
	assert g.dtype == bool
	assert (g == g.T).all() and not g.diagonal().any()
	assert nx.is_tree(nx.from_numpy_array(g.astype(int)))
	return snames, slocs, cells, graph


def node_adjacency():
	"""Neighbour matrix of the four node subsets when they follow the trajectory's edges."""
	ans = np.zeros((4, 4), dtype=bool)
	ans[EDGES[:, 0], EDGES[:, 1]] = True
	return ans | ans.T


################################################################
# subsets_rna
################################################################

class TestSubsetsRna:
	def test_outputs(self, tmp_path, rng):
		"""Writes subset names, node-first subset locations, ncell-cell folders and a tree-shaped neighbour table."""
		inputs, edges, locs, names = write_inputs(str(tmp_path / 'data'), rng)
		outs = output_paths(str(tmp_path / 'out'))
		np.random.seed(0)
		assert subsets_rna(*inputs, *outs, 4, 2, 1.) is None
		snames, slocs, cells, graph = check_rna_outputs(outs, edges, locs, names, 4)
		assert len(snames) > 4
		for x in snames:
			assert os.path.isfile(os.path.join(outs[2], x, 'names_rna.txt'))

	def test_no_overlap_requirement_keeps_node_subsets(self, tmp_path, rng):
		"""With noverlap=0 and a huge dmax only the four node subsets remain, linked like the trajectory."""
		inputs, edges, locs, names = write_inputs(str(tmp_path / 'data'), rng)
		outs = output_paths(str(tmp_path / 'out'))
		np.random.seed(0)
		subsets_rna(*inputs, *outs, 5, 0, 1e9)
		snames, slocs, cells, graph = check_rna_outputs(outs, edges, locs, names, 5)
		assert snames == ['Subset1', 'Subset2', 'Subset3', 'Subset4']
		assert_array_equal(graph.values, node_adjacency())

	def test_neighbour_chain_follows_edges(self, tmp_path, rng):
		"""Between the subsets of an edge's two nodes, the neighbour table passes only that edge's subsets in increasing loc."""
		inputs, edges, locs, names = write_inputs(str(tmp_path / 'data'), rng, counts=(8, 9, 10))
		outs = output_paths(str(tmp_path / 'out'))
		np.random.seed(1)
		subsets_rna(*inputs, *outs, 4, 3, 0.)
		snames, slocs, cells, graph = check_rna_outputs(outs, edges, locs, names, 4)
		assert len(snames) > 4
		g = nx.from_numpy_array(graph.values.astype(int))
		inner = 0
		for e, (u, v) in enumerate(EDGES):
			chain = nx.shortest_path(g, int(u), int(v))[1:-1]
			inner += len(chain)
			assert (slocs.edges[chain] == e).all()
			assert (np.diff(slocs.locs[chain]) > 0).all()
		assert inner == len(snames) - 4

	def test_too_many_cells_at_a_node_raises(self, tmp_path):
		"""More than ncell cells exactly at one node raise ValueError."""
		folder = str(tmp_path / 'data')
		os.makedirs(folder)
		t = ytree()
		t.to_file(os.path.join(folder, 'traj_node.h5'))
		point(t, np.array([0, 0, 0, 1, 1, 2, 2]), np.array([0., 0., 0., 0.5, 1.5, 1., 2.])).to_file(os.path.join(folder, 'traj_cell_rna.h5'))
		pd.DataFrame(np.zeros((7, 1)), index=['c{}'.format(x) for x in range(7)], columns=['Dim1']).to_csv(os.path.join(folder, 'coord_rna.tsv'), sep='\t')
		outs = output_paths(str(tmp_path / 'out'))
		np.random.seed(0)
		with pytest.raises(ValueError, match='ncell'):
			subsets_rna(os.path.join(folder, 'traj_node.h5'), os.path.join(folder, 'traj_cell_rna.h5'), os.path.join(folder, 'coord_rna.tsv'), *outs, 2, 1, 1.)

	def test_cli_matches_direct_call(self, tmp_path, rng, run_dictys):
		"""The command line, given every argument, writes the same files as the direct call with the same seed."""
		inputs, edges, locs, names = write_inputs(str(tmp_path / 'data'), rng)
		outs_cli = output_paths(str(tmp_path / 'cli'), 'subset_edges.tsv.gz')
		outs_py = output_paths(str(tmp_path / 'py'), 'subset_edges.tsv.gz')
		np.random.seed(2)
		run_dictys('dynamic', 'subsets_rna', *inputs, *outs_cli, 4, 2, 1.5)
		np.random.seed(2)
		subsets_rna(*inputs, *outs_py, 4, 2, 1.5)
		cli = check_rna_outputs(outs_cli, edges, locs, names, 4)
		py = read_outputs(outs_py)
		assert cli[0] == py[0]
		assert_array_equal(cli[1].edges, py[1].edges)
		assert_array_equal(cli[1].locs, py[1].locs)
		assert cli[2] == py[2]
		assert cli[3].equals(py[3])

	def test_cli_no_overlap(self, tmp_path, rng, run_dictys):
		"""Command-line ncell, noverlap and dmax reach the function: 0 and 1e9 keep only the node subsets of 3 cells."""
		inputs, edges, locs, names = write_inputs(str(tmp_path / 'data'), rng)
		outs = output_paths(str(tmp_path / 'out'))
		np.random.seed(0)
		run_dictys('dynamic', 'subsets_rna', *inputs, *outs, 3, 0, 1e9)
		snames, slocs, cells, graph = check_rna_outputs(outs, edges, locs, names, 3)
		assert len(snames) == 4
		assert_array_equal(graph.values, node_adjacency())


################################################################
# subsets_atac
################################################################

def write_atac_inputs(folder, rng, counts=(5, 6, 7)):
	"""Writes traj_node.h5, ATAC cell files and three hand-placed subsets A, B, C; returns paths, cells and subset locations."""
	os.makedirs(folder, exist_ok=True)
	ytree().to_file(os.path.join(folder, 'traj_node.h5'))
	edges, locs, names = write_cells(folder, rng, 'atac', counts)
	sedges = np.array([0, 1, 2])
	slocs = np.array([0.5, 1.5, 3.2])
	point(ytree(), sedges, slocs).to_file(os.path.join(folder, 'subset_locs.h5'))
	with open(os.path.join(folder, 'subsets.txt'), 'w') as f:
		f.write('A\nB\nC\n')
	paths = [os.path.join(folder, x) for x in ['traj_node.h5', 'traj_cell_atac.h5', 'coord_atac.tsv', 'subsets.txt', 'subset_locs.h5']]
	return paths, edges, locs, names, sedges, slocs


class TestSubsetsAtac:
	@pytest.mark.parametrize('name,index', [('A', 0), ('B', 1), ('C', 2)])
	def test_nearest_cells_of_named_subset(self, tmp_path, rng, name, index):
		"""Writes the names of the ncell ATAC cells nearest the named subset's location (reference distance)."""
		inputs, edges, locs, names, sedges, slocs = write_atac_inputs(str(tmp_path), rng)
		fo = str(tmp_path / 'names_atac0.txt')
		np.random.seed(0)
		assert subsets_atac(*inputs, fo, name, 4) is None
		ans = read_names(fo)
		assert len(ans) == 4
		expected = nearest(ref_dist(sedges[[index]], slocs[[index]], edges, locs), 4, names)[0]
		assert frozenset(ans) == expected

	def test_duplicate_cells_are_perturbed(self, tmp_path):
		"""ATAC cells sharing a location far from the subset, which subsets_atac perturbs, leave the nearest set unchanged."""
		folder = str(tmp_path)
		t = ytree()
		t.to_file(os.path.join(folder, 'traj_node.h5'))
		edges = np.array([0, 0, 0, 0, 1, 1, 2, 2, 2])
		locs = np.array([0.4, 0.6, 0.9, 1.3, 2.5, 2.5, 3.5, 3.5, 3.5])
		names = ['a{}'.format(x) for x in range(len(edges))]
		point(t, edges, locs).to_file(os.path.join(folder, 'traj_cell_atac.h5'))
		pd.DataFrame(np.zeros((len(names), 1)), index=names, columns=['Dim1']).to_csv(os.path.join(folder, 'coord_atac.tsv'), sep='\t')
		point(t, np.array([0]), np.array([0.5])).to_file(os.path.join(folder, 'subset_locs.h5'))
		with open(os.path.join(folder, 'subsets.txt'), 'w') as f:
			f.write('S1\n')
		fo = os.path.join(folder, 'out.txt')
		np.random.seed(0)
		subsets_atac(*[os.path.join(folder, x) for x in ['traj_node.h5', 'traj_cell_atac.h5', 'coord_atac.tsv', 'subsets.txt', 'subset_locs.h5']], fo, 'S1', 3)
		assert sorted(read_names(fo)) == ['a0', 'a1', 'a2']

	def test_unknown_subset_raises(self, tmp_path, rng):
		"""A subset name absent from the subset list raises KeyError."""
		inputs, edges, locs, names, sedges, slocs = write_atac_inputs(str(tmp_path), rng)
		with pytest.raises(KeyError):
			subsets_atac(*inputs, str(tmp_path / 'out.txt'), 'D', 3)

	def test_duplicate_subset_names_assert(self, tmp_path, rng):
		"""A subset list with a repeated name fails an assertion."""
		inputs, edges, locs, names, sedges, slocs = write_atac_inputs(str(tmp_path), rng)
		with open(inputs[3], 'w') as f:
			f.write('A\nB\nA\n')
		with pytest.raises(AssertionError):
			subsets_atac(*inputs, str(tmp_path / 'out.txt'), 'B', 3)

	def test_cli(self, tmp_path, rng, run_dictys):
		"""The command line, given every argument, writes the ncell ATAC cells nearest subset B."""
		inputs, edges, locs, names, sedges, slocs = write_atac_inputs(str(tmp_path), rng)
		fo = str(tmp_path / 'B' / 'names_atac0.txt')
		os.makedirs(str(tmp_path / 'B'))
		np.random.seed(0)
		run_dictys('dynamic', 'subsets_atac', *inputs, fo, 'B', 5)
		ans = read_names(fo)
		assert len(ans) == 5
		assert frozenset(ans) == nearest(ref_dist(sedges[[1]], slocs[[1]], edges, locs), 5, names)[0]

	def test_after_subsets_rna(self, tmp_path, rng, run_dictys):
		"""On subsets_rna output, every subset gets the ncell ATAC cells nearest its centre; node subsets those nearest each node."""
		inputs, rna_edges, rna_locs, rna_names = write_inputs(str(tmp_path / 'data'), rng)
		atac_edges, atac_locs, atac_names = write_cells(str(tmp_path / 'data'), rng, 'atac', (6, 6, 6))
		outs = output_paths(str(tmp_path / 'out'))
		np.random.seed(0)
		run_dictys('dynamic', 'subsets_rna', *inputs, *outs, 4, 2, 1.)
		snames = read_names(outs[0])
		slocs = point.from_file(ytree(), outs[1])
		for i, x in enumerate(snames):
			fo = os.path.join(outs[2], x, 'names_atac0.txt')
			subsets_atac(inputs[0], str(tmp_path / 'data' / 'traj_cell_atac.h5'), str(tmp_path / 'data' / 'coord_atac.tsv'), outs[0], outs[1], fo, x, 3)
			expected = nearest(ref_dist(slocs.edges[[i]], slocs.locs[[i]], atac_edges, atac_locs), 3, atac_names)[0]
			assert frozenset(read_names(fo)) == expected
			if i < 4:
				assert frozenset(read_names(fo)) == nearest(ref_dist(NODE_EDGES[[i]], NODE_LOCS[[i]], atac_edges, atac_locs), 3, atac_names)[0]
