"""
Tests for dictys.net: the network and dynamic_network classes.

Networks are built by hand through network.__init__, and the pipeline
folder layout that network.from_folders reads is written to tmp_path.
"""

import os
import sys

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose, assert_array_equal

import h5py
from dictys.net import network, dynamic_network
from dictys.traj import trajectory, point

CNAME = np.array(['c0', 'c1', 'c2', 'c3', 'c4', 'c5'])
SNAME = np.array(['s0', 's1', 's2'])
NNAME = np.array(['TF1', 'TF2', 'G3', 'G4', 'G5'])
CPM = np.array([
	[1e5, 2e5, 4e5],
	[3e5, 0, 1e5],
	[2e5, 3e5, 1e5],
	[4e5, 4e5, 3e5],
	[0, 1e5, 1e5],
])


@pytest.fixture(autouse=True)
def _close_figures():
	yield
	if 'matplotlib.pyplot' in sys.modules:
		sys.modules['matplotlib.pyplot'].close('all')


def static_w():
	"""Edge weights (2 regulators, 5 targets, 3 states); state s2 has no edge."""
	w = np.zeros((2, 5, 3))
	w[:, :, 0] = [[0, 0.1, 0.9, -0.4, 0], [0.7, 0, -0.2, 0, 0.3]]
	w[:, :, 1] = [[0, 0.5, 0, 0, 0], [0, 0, 0, 0.6, 0]]
	return w


def static_params(**ka):
	w = static_w()
	params = dict(
		cname=CNAME.copy(), sname=SNAME.copy(), nname=NNAME.copy(),
		nids=[np.array([0, 1]), np.arange(5)],
		esprop={'w': w, 'mask': w != 0},
		nsprop={'cpm': CPM.copy()},
		scprop={'w': np.repeat(np.eye(3), 2, axis=1)},
		ncprop={'readcount': np.arange(30).reshape(5, 6)},
		cprop={'type': np.array(['a', 'a', 'b', 'b', 'c', 'c'])},
	)
	params.update(ka)
	return params


def make_static(**ka):
	return network(**static_params(**ka))


def build_dynamic(traj_edges, traj_lens, state_edges, state_locs, neighbors, seed=0, cls=dynamic_network):
	"""Dynamic network with two cells sitting at each state's location."""
	rng = np.random.default_rng(seed)
	traj = trajectory(np.array(traj_edges), np.array(traj_lens, dtype=float))
	se = np.array(state_edges)
	sl = np.array(state_locs, dtype=float)
	ns = len(se)
	nb = np.zeros((ns, ns), dtype=bool)
	for a, b in neighbors:
		nb[a, b] = nb[b, a] = True
	w = rng.normal(size=(2, 5, ns))
	cpm = rng.uniform(1e4, 2e5, size=(5, ns))
	return cls(
		cname=np.array(['c%02d' % x for x in range(2 * ns)]),
		sname=np.array(['s%d' % x for x in range(ns)]),
		nname=NNAME.copy(), nids=[np.array([0, 1]), np.arange(5)],
		esprop={'w': w, 'mask': np.ones(w.shape, dtype=bool)},
		nsprop={'cpm': cpm},
		scprop={'w': np.repeat(np.eye(ns), 2, axis=1)},
		ssprop={'traj-neighbor': nb},
		traj=traj,
		point={'c': point(traj, np.repeat(se, 2), np.repeat(sl, 2)), 's': point(traj, se, sl)},
	)


STAR = dict(
	traj_edges=[[0, 1], [1, 2], [1, 3]], traj_lens=[1, 1, 1],
	state_edges=[0, 0, 0, 1, 1, 2, 2], state_locs=[0, 0.5, 1, 0.5, 1, 0.5, 1],
	neighbors=[(0, 1), (1, 2), (2, 3), (3, 4), (2, 5), (5, 6)],
)
CHAIN = dict(
	traj_edges=[[0, 1], [1, 2]], traj_lens=[1, 1],
	state_edges=[0, 0, 1], state_locs=[0, 1, 1],
	neighbors=[(0, 1), (1, 2)],
)


@pytest.fixture
def star():
	"""
	Dynamic network on a star trajectory 0-1, 1-2, 1-3 (unit lengths). States
	s0..s4 lie on the path 0->2 at distances 0, 0.5, 1, 1.5, 2; s5 and s6 lie
	on the side branch 1->3.
	"""
	return build_dynamic(**STAR)


def path_smoother(num, dist):
	"""Gaussian kernel weights of states s0..s4 at num even points on path 0->2, and their pseudotime."""
	x = np.array([0, 0.5, 1, 1.5, 2])
	t = np.linspace(0, 2, num)
	w = np.exp(-((t[None, :] - x[:, None]) / dist) ** 2 / 2)
	return w / w.sum(axis=0), t


def binarize_top(v, k):
	"""Top-k entries by absolute value in each slice v[:,:,i]."""
	a = np.abs(v)
	ans = np.zeros(v.shape, dtype=bool)
	for i in range(v.shape[2]):
		order = np.argsort(a[:, :, i].ravel())[::-1][:k]
		ans[:, :, i].flat[order] = True
	return ans


def curve_chars(y, t):
	"""
	Curve characteristics from their definitions: end minus start; area of the
	excursion beyond the range of the start and end values; and the area of the
	start-end-clipped curve above its end value relative to the total change.
	"""
	x = (t - t[0]) / (t[-1] - t[0])

	def area(v):
		return ((v[:, 1:] + v[:, :-1]) / 2) @ np.diff(x)

	lo = np.minimum(y[:, [0]], y[:, [-1]])
	hi = np.maximum(y[:, [0]], y[:, [-1]])
	clipped = np.clip(y, lo, hi)
	change = y[:, 0] - y[:, -1]
	switching = area(clipped - y[:, [-1]]) / np.where(change == 0, 1, change)
	return pd.DataFrame({
		'Terminal logFC': y[:, -1] - y[:, 0],
		'Transient logFC': area(y - clipped),
		'Switching time': switching,
	})


def read_tsv(path, **ka):
	return pd.read_csv(path, sep='\t', header=0, index_col=0, **ka)


# network.__init__


def test_init_counts_names_and_props():
	"""__init__ fills counts, name dictionaries, source/target counts and the per-owner property dicts."""
	d = make_static()
	assert (d.cn, d.sn, d.nn) == (6, 3, 5)
	assert_array_equal(d.nns, [2, 5])
	assert d.cdict == {'c%d' % x: x for x in range(6)}
	assert d.sdict == {'s0': 0, 's1': 1, 's2': 2}
	assert d.ndict == {x: i for i, x in enumerate(NNAME)}
	assert set(d.prop) == {'es', 'ns', 'sc', 'nc', 'c'}
	assert_array_equal(d.prop['es']['w'], static_w())
	assert_array_equal(d.prop['c']['type'], ['a', 'a', 'b', 'b', 'c', 'c'])
	assert d.traj is None
	assert d.point == {}


def test_init_prop_argument_replaces_property_dict():
	"""A 'prop' argument is stored as the whole property dict."""
	params = static_params()
	prop = {'ns': {'cpm': CPM.copy()}}
	for x in ['esprop', 'nsprop', 'scprop', 'ncprop', 'cprop']:
		del params[x]
	d = network(prop=prop, **params)
	assert d.prop is prop


def test_init_rejects_unknown_parameter():
	"""__init__ raises ValueError for a keyword it does not know."""
	with pytest.raises(ValueError, match='Unknown parameter'):
		make_static(gname=np.array(['x']))


@pytest.mark.parametrize('name,dim', [('cname', 'cell'), ('sname', 'state'), ('nname', 'gene')])
def test_init_rejects_empty_dimension(name, dim):
	"""__init__ raises ValueError when there is no cell, state or gene."""
	with pytest.raises(ValueError, match=f'No {dim} included'):
		network(**static_params(**{name: np.array([], dtype=str)}, esprop={}, nsprop={}, scprop={}, ncprop={}, cprop={}))


def test_init_rejects_duplicate_names():
	"""__init__ raises ValueError naming the duplicated cell."""
	with pytest.raises(ValueError, match='duplicate cell names. First three: c1'):
		make_static(cname=np.array(['c0', 'c1', 'c1', 'c3', 'c4', 'c5']))


def test_check_rejects_name_count_mismatch():
	"""check raises ValueError when a name array no longer matches its count."""
	d = make_static()
	d.sname = np.array(['s0', 's1'])
	with pytest.raises(ValueError, match='Expecting 3 state names, actual count: 2'):
		d.check()


@pytest.mark.parametrize('key,value', [
	('esprop', {'w': np.zeros((5, 2, 3))}),
	('nsprop', {'cpm': np.zeros((3, 5))}),
	('scprop', {'w': np.zeros((6, 3))}),
	('cprop', {'coord': np.zeros((2, 5))}),
])
def test_init_rejects_property_with_wrong_trailing_shape(key, value):
	"""__init__ raises ValueError when a property's trailing axes do not match its owner type."""
	with pytest.raises(ValueError, match='incorrect size'):
		make_static(**{key: value})


def test_init_accepts_leading_axes_and_multi_owner_properties():
	"""Extra leading axes are allowed, and an 'ss' property needs trailing shape (sn, sn)."""
	d = make_static(cprop={'coord': np.zeros((2, 6))}, ssprop={'adj': np.eye(3)})
	assert d.prop['c']['coord'].shape == (2, 6)
	assert d.prop['ss']['adj'].shape == (3, 3)


# network.check


def test_check_detects_property_added_later():
	"""check raises ValueError for a property added after construction with a wrong shape."""
	d = make_static()
	d.check()
	d.prop['es']['bad'] = np.zeros((2, 4, 3))
	with pytest.raises(ValueError, match='Properties of type es'):
		d.check()


def test_check_dynamic_network_points(star):
	"""check accepts cell/state points of matching sizes and rejects a state point list of another length."""
	star.check()
	assert len(star.point['c']) == star.cn
	star.point['s'] = star.point['s'][np.arange(star.sn - 1)]
	with pytest.raises(AssertionError):
		star.check()


# network.rename


def test_rename_with_list():
	"""rename with a list replaces every node name and the name-to-index dict."""
	d = make_static()
	d.rename('n', ['A', 'B', 'C', 'D', 'E'])
	assert_array_equal(d.nname, ['A', 'B', 'C', 'D', 'E'])
	assert d.ndict == {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4}
	assert d.nn == 5


def test_rename_with_partial_dict():
	"""rename with a dict renames only the listed cells or states."""
	d = make_static()
	d.rename('c', {'c0': 'x', 'c4': 'y'})
	assert_array_equal(d.cname, ['x', 'c1', 'c2', 'c3', 'y', 'c5'])
	assert d.cdict['y'] == 4 and 'c0' not in d.cdict
	d.rename('s', {'s2': 'late'})
	assert_array_equal(d.sname, ['s0', 's1', 'late'])
	assert d.sdict == {'s0': 0, 's1': 1, 'late': 2}


def test_rename_errors():
	"""rename rejects unknown dimensions, wrong list length, other value types and duplicate results."""
	d = make_static()
	with pytest.raises(ValueError):
		d.rename('x', ['a'])
	with pytest.raises(TypeError):
		d.rename('s', ['a', 'b'])
	with pytest.raises(TypeError):
		d.rename('s', ('a', 'b', 'c'))
	with pytest.raises(ValueError, match='Duplicate'):
		d.rename('s', {'s0': 's1'})
	assert_array_equal(d.sname, SNAME)


# network.to_file and network.from_file


def assert_same_network(a, b):
	for x in ['cname', 'sname', 'nname']:
		assert_array_equal(getattr(a, x), getattr(b, x))
		assert getattr(b, x).dtype.kind == 'U'
	for x, y in zip(a.nids, b.nids):
		assert_array_equal(x, y)
	assert set(a.prop) == set(b.prop)
	for x in a.prop:
		assert set(a.prop[x]) == set(b.prop[x])
		for y in a.prop[x]:
			assert_array_equal(a.prop[x][y], b.prop[x][y])


def test_to_file_from_file_round_trip(tmp_path):
	"""A static network written by to_file reads back identically, string properties included."""
	d = make_static()
	path = os.path.join(tmp_path, 'net.h5')
	d.to_file(path)
	d2 = network.from_file(path)
	assert type(d2) is network
	assert_same_network(d, d2)
	assert d2.traj is None and d2.point == {}
	with h5py.File(path, 'r') as f:
		assert_array_equal(f['_version_'], [0, 1, 0])
		assert f['cname'].compression == 'gzip'
		assert f['prop']['es']['w'].compression == 'gzip'


def test_to_file_from_file_dynamic_round_trip(tmp_path, star):
	"""A dynamic network keeps its trajectory and cell/state points through to_file and from_file."""
	path = os.path.join(tmp_path, 'dyn.h5')
	star.to_file(path)
	d2 = dynamic_network.from_file(path)
	assert type(d2) is dynamic_network
	assert_same_network(star, d2)
	assert_array_equal(d2.traj.edges, star.traj.edges)
	assert_array_equal(d2.traj.lens, star.traj.lens)
	for x in ['c', 's']:
		assert_array_equal(d2.point[x].edges, star.point[x].edges)
		assert_allclose(d2.point[x].locs, star.point[x].locs)
	pts1 = star.linspace(0, 2, 4, 0.5)[0]
	pts2 = d2.linspace(0, 2, 4, 0.5)[0]
	assert_allclose(pts2.locs, pts1.locs)


def test_to_file_options(tmp_path):
	"""to_file applies compression and extra dataset keywords, and links repeated arrays only when links=True."""
	w = static_w()
	d = make_static(esprop={'w': w, 'w_n': w})
	p1 = os.path.join(tmp_path, 'a.h5')
	p2 = os.path.join(tmp_path, 'b.h5')
	d.to_file(p1, compression='lzf', links=True, fletcher32=True)
	d.to_file(p2, compression=None, links=False)
	with h5py.File(p1, 'r') as f:
		assert f['nname'].compression == 'lzf'
		assert f['prop']['ns']['cpm'].compression == 'lzf'
		assert f['nname'].fletcher32 and f['prop']['ns']['cpm'].fletcher32
		assert f['prop']['es']['w'].id == f['prop']['es']['w_n'].id
	with h5py.File(p2, 'r') as f:
		assert f['nname'].compression is None
		assert not f['nname'].fletcher32
		assert f['prop']['es']['w'].id != f['prop']['es']['w_n'].id
	for p in [p1, p2]:
		d2 = network.from_file(p)
		assert_array_equal(d2.prop['es']['w_n'], w)


def test_from_file_without_version_reads_as_0_1_0(tmp_path):
	"""A file without a _version_ entry is read as version 0.1.0."""
	path = os.path.join(tmp_path, 'net.h5')
	make_static().to_file(path)
	with h5py.File(path, 'a') as f:
		del f['_version_']
	assert_same_network(make_static(), network.from_file(path))


def test_from_file_rejects_unknown_version(tmp_path):
	"""from_file raises TypeError for a file version it does not know."""
	path = os.path.join(tmp_path, 'net.h5')
	make_static().to_file(path)
	with h5py.File(path, 'a') as f:
		f['_version_'][:] = [0, 2, 0]
	with pytest.raises(TypeError, match='Unknown network version 0.2.0'):
		network.from_file(path)


# network.from_folders


def write_tsv(path, df):
	df.to_csv(path, sep='\t', header=True, index=True)


def write_lines(path, lines):
	with open(path, 'w') as f:
		f.write('\n'.join(lines) + '\n')


STATIC_CELLS = {'A': ['a1', 'a2', 'a3'], 'B': ['b1', 'b2']}
STATIC_NETS = {
	'A': pd.DataFrame([[0.1, 0.2, 0.3, 0.4], [0.5, 0.6, 0.7, 0.8]], index=['TF1', 'TF2'], columns=['TF1', 'TF2', 'G3', 'G4']),
	'B': pd.DataFrame([[1.1, 1.2, 1.3]], index=['TF1'], columns=['TF1', 'G3', 'G5']),
}
STATIC_GENES = ['G3', 'G4', 'G5', 'TF1', 'TF2']
STATIC_COUNTS = pd.DataFrame(
	[[1, 0, 3, 2, 2], [4, 1, 0, 0, 1], [0, 2, 2, 5, 0], [3, 3, 1, 1, 4], [2, 0, 6, 0, 300]],
	index=STATIC_GENES, columns=['a1', 'a2', 'a3', 'b1', 'b2'])


def make_static_folders(root, nettype='n', skip=(), empty=()):
	"""
	Pipeline layout for two static subsets A and B: data/expression.tsv.gz,
	work/subsets.txt, work/<subset>/names_rna.txt and the direct (scaled by 1)
	and steady-state (scaled by -1) network files. (subset, prefix) pairs in
	skip are not written and those in empty are written with no rows.
	"""
	data = os.path.join(root, 'data')
	work = os.path.join(root, 'work')
	os.makedirs(data)
	os.makedirs(work)
	write_tsv(os.path.join(data, 'expression.tsv.gz'), STATIC_COUNTS)
	fi_subsets = os.path.join(work, 'subsets.txt')
	write_lines(fi_subsets, list(STATIC_CELLS))
	for s, cells in STATIC_CELLS.items():
		os.makedirs(os.path.join(work, s))
		write_lines(os.path.join(work, s, 'names_rna.txt'), cells)
		for prefix, scale in [(nettype, 1), ('i' + nettype, -1)]:
			if (s, prefix) in skip:
				continue
			net = STATIC_NETS[s] * scale
			if (s, prefix) in empty:
				net = net.iloc[:0]
			write_tsv(os.path.join(work, s, f'net_{prefix}weight.tsv.gz'), net)
	return data, work, fi_subsets


def expected_edges(states, scale=1):
	"""Edge weights and masks over sorted regulators x sorted targets x states, filled from STATIC_NETS."""
	regs = ['TF1', 'TF2']
	w = np.zeros((2, 5, len(states)))
	mask = np.zeros(w.shape, dtype=bool)
	for k, s in enumerate(states):
		net = STATIC_NETS[s]
		for r in net.index:
			for t in net.columns:
				w[regs.index(r), STATIC_GENES.index(t), k] = net.loc[r, t] * scale
				mask[regs.index(r), STATIC_GENES.index(t), k] = True
	return w, mask


def test_from_folders_static(tmp_path):
	"""from_folders reads names, membership, read counts, pseudobulk CPM and both network types of static subsets."""
	data, work, fi_subsets = make_static_folders(str(tmp_path))
	d = network.from_folders(data, work, fi_subsets)
	assert type(d) is network
	assert_array_equal(d.sname, ['A', 'B'])
	assert_array_equal(d.cname, ['a1', 'a2', 'a3', 'b1', 'b2'])
	assert_array_equal(d.nname, STATIC_GENES)
	assert_array_equal(d.nids[0], [3, 4])
	assert_array_equal(d.nids[1], [0, 1, 2, 3, 4])
	assert_array_equal(d.prop['sc']['w'], [[1, 1, 1, 0, 0], [0, 0, 0, 1, 1]])
	counts = STATIC_COUNTS.values
	assert_array_equal(d.prop['nc']['readcount'], counts)
	assert np.issubdtype(d.prop['nc']['readcount'].dtype, np.integer)
	pseudobulk = np.array([counts[:, :3].sum(axis=1), counts[:, 3:].sum(axis=1)]).T
	assert_allclose(d.prop['ns']['cpm'], 1e6 * pseudobulk / pseudobulk.sum(axis=0))
	w, mask = expected_edges(['A', 'B'])
	assert_allclose(d.prop['es']['w_n'], w)
	assert_array_equal(d.prop['es']['mask_n'], mask)
	assert_allclose(d.prop['es']['w_in'], -w)
	assert_array_equal(d.prop['es']['mask_in'], mask)
	assert_allclose(d.prop['es']['w'], w)
	assert_array_equal(d.prop['es']['mask'], mask)
	assert d.traj is None


def test_from_folders_unnormalized_type_and_no_readcount(tmp_path):
	"""nettype='' reads net_weight/net_iweight files, and optional=set() leaves out read counts."""
	data, work, fi_subsets = make_static_folders(str(tmp_path), nettype='')
	d = network.from_folders(data, work, fi_subsets, dynamic=False, nettype='', optional=set())
	assert set(d.prop['es']) == {'w_', 'mask_', 'w_i', 'mask_i', 'w', 'mask'}
	w, mask = expected_edges(['A', 'B'])
	assert_allclose(d.prop['es']['w_'], w)
	assert_allclose(d.prop['es']['w_i'], -w)
	assert 'readcount' not in d.prop['nc']


def test_from_folders_extra_cell_properties(tmp_path):
	"""fi_c adds one cell property per column, reordered to the network's cells."""
	data, work, fi_subsets = make_static_folders(str(tmp_path))
	fi_c = os.path.join(str(tmp_path), 'cells.tsv')
	write_tsv(fi_c, pd.DataFrame({'batch': [2, 1, 1, 2, 2, 3], 'label': list('qrstuv')}, index=['b2', 'a1', 'a2', 'a3', 'b1', 'extra']))
	d = network.from_folders(data, work, fi_subsets, fi_c=fi_c)
	assert_array_equal(d.prop['c']['batch'], [1, 1, 2, 2, 2])
	assert_array_equal(d.prop['c']['label'], ['r', 's', 't', 'u', 'q'])


def test_from_folders_extra_cell_properties_missing_cell(tmp_path):
	"""fi_c that misses a known cell raises ValueError."""
	data, work, fi_subsets = make_static_folders(str(tmp_path))
	fi_c = os.path.join(str(tmp_path), 'cells.tsv')
	write_tsv(fi_c, pd.DataFrame({'batch': [1, 2]}, index=['a1', 'a2']))
	with pytest.raises(ValueError, match='Not all known cells'):
		network.from_folders(data, work, fi_subsets, fi_c=fi_c)


def test_from_folders_drops_subset_without_networks(tmp_path):
	"""A subset whose network files are all missing is removed from every state axis; its cells stay."""
	data, work, fi_subsets = make_static_folders(str(tmp_path), skip={('B', 'n'), ('B', 'in')})
	d = network.from_folders(data, work, fi_subsets)
	assert_array_equal(d.sname, ['A'])
	assert_array_equal(d.cname, ['a1', 'a2', 'a3', 'b1', 'b2'])
	assert_array_equal(d.prop['sc']['w'], [[1, 1, 1, 0, 0]])
	assert_array_equal(d.nname, ['G3', 'G4', 'TF1', 'TF2'])
	assert d.prop['ns']['cpm'].shape == (4, 1)
	assert d.prop['es']['w'].shape == (2, 4, 1)
	assert_allclose(d.prop['es']['w'][:, :, 0], STATIC_NETS['A'][['G3', 'G4', 'TF1', 'TF2']].values)


def test_from_folders_keeps_subset_with_one_network_type(tmp_path):
	"""A subset missing one network type, or with an empty network file, keeps an all-False mask for it."""
	data, work, fi_subsets = make_static_folders(str(tmp_path), skip={('B', 'in')}, empty={('A', 'in')})
	d = network.from_folders(data, work, fi_subsets)
	assert_array_equal(d.sname, ['A', 'B'])
	w, mask = expected_edges(['A', 'B'])
	assert_allclose(d.prop['es']['w_n'], w)
	assert not d.prop['es']['mask_in'].any()
	assert not d.prop['es']['w_in'].any()


def test_from_folders_errors(tmp_path):
	"""from_folders rejects unknown optional data, a missing subset list and folders without any network."""
	data, work, fi_subsets = make_static_folders(str(tmp_path), skip={('A', 'n'), ('A', 'in'), ('B', 'n'), ('B', 'in')})
	with pytest.raises(ValueError, match='Unknown optional data types: bad'):
		network.from_folders(data, work, fi_subsets, optional={'bad'})
	with pytest.raises(RuntimeError, match='Could not find cell subsets'):
		network.from_folders(data, work, os.path.join(work, 'missing.txt'))
	with pytest.raises(RuntimeError, match='No network found'):
		network.from_folders(data, work, fi_subsets)


def make_dynamic_folders(root, cells_extra=None):
	"""
	Pipeline layout for a dynamic network on the star trajectory: the
	trajectory, cell and subset locations, subset neighbours, low-dimensional
	coordinates, cell types and one network pair per subset.
	"""
	rng = np.random.default_rng(3)
	data = os.path.join(root, 'data')
	work = os.path.join(root, 'work')
	os.makedirs(data)
	os.makedirs(work)
	traj = trajectory(np.array(STAR['traj_edges']), np.array(STAR['traj_lens'], dtype=float))
	se = np.array(STAR['state_edges'])
	sl = np.array(STAR['state_locs'], dtype=float)
	ns = len(se)
	cells = ['c%02d' % x for x in range(2 * ns)]
	snames = ['Subset%d' % (x + 1) for x in range(ns)]
	traj.to_file(os.path.join(data, 'traj_node.h5'))
	point(traj, np.repeat(se, 2), np.repeat(sl, 2)).to_file(os.path.join(data, 'traj_cell_rna.h5'))
	point(traj, se, sl).to_file(os.path.join(work, 'subset_locs.h5'))
	nb = np.zeros((ns, ns), dtype=bool)
	for a, b in STAR['neighbors']:
		nb[a, b] = nb[b, a] = True
	write_tsv(os.path.join(work, 'subset_edges.tsv.gz'), pd.DataFrame(nb, index=snames, columns=snames))
	coord = pd.DataFrame(rng.normal(size=(2 * ns, 2)), index=cells, columns=['x', 'y'])
	write_tsv(os.path.join(data, 'coord_rna.tsv.gz'), coord)
	types = {'early': cells[:6], 'late': cells[6:]}
	write_lines(os.path.join(data, 'subsets.txt'), list(types))
	for x, y in types.items():
		os.makedirs(os.path.join(data, 'subsets', x))
		write_lines(os.path.join(data, 'subsets', x, 'names_rna.txt'), y)
	counts = pd.DataFrame(rng.integers(0, 300, size=(5, 2 * ns)), index=NNAME, columns=cells)
	write_tsv(os.path.join(data, 'expression.tsv.gz'), counts)
	write_lines(os.path.join(work, 'subsets.txt'), snames)
	members = [cells[2 * x:2 * x + 2] for x in range(ns)]
	members[1] = members[1] + [cells[0]]
	nets = []
	for x in range(ns):
		os.makedirs(os.path.join(work, snames[x]))
		write_lines(os.path.join(work, snames[x], 'names_rna.txt'), members[x])
		net = pd.DataFrame(rng.normal(size=(2, 5)), index=['TF1', 'TF2'], columns=['TF1', 'TF2', 'G3', 'G4', 'G5'])
		nets.append(net)
		write_tsv(os.path.join(work, snames[x], 'net_nweight.tsv.gz'), net)
		write_tsv(os.path.join(work, snames[x], 'net_inweight.tsv.gz'), net * 2)
	return data, work, os.path.join(work, 'subsets.txt'), dict(traj=traj, nb=nb, coord=coord, members=members, nets=nets, cells=cells, snames=snames)


def test_from_folders_dynamic(tmp_path):
	"""dynamic=True also reads the trajectory, points, subset neighbours, coordinates and cell types."""
	data, work, fi_subsets, ref = make_dynamic_folders(str(tmp_path))
	d = dynamic_network.from_folders(data, work, fi_subsets, dynamic=True)
	assert type(d) is dynamic_network
	assert_array_equal(d.cname, sorted(ref['cells']))
	assert_array_equal(d.sname, ref['snames'])
	assert_array_equal(d.traj.edges, ref['traj'].edges)
	assert_allclose(d.traj.lens, ref['traj'].lens)
	assert_array_equal(d.point['s'].edges, STAR['state_edges'])
	assert_allclose(d.point['s'].locs, STAR['state_locs'])
	assert_array_equal(d.point['c'].edges, np.repeat(STAR['state_edges'], 2))
	assert_array_equal(d.prop['ss']['traj-neighbor'], ref['nb'])
	assert_allclose(d.prop['c']['coord'], ref['coord'].loc[d.cname].values.T)
	assert_array_equal(d.prop['c']['type'], ['early'] * 6 + ['late'] * 8)
	membership = np.array([[x in m for x in d.cname] for m in ref['members']], dtype=float)
	assert_array_equal(d.prop['sc']['w'], membership)
	assert membership[1].sum() == 3
	for x in range(len(ref['nets'])):
		assert_allclose(d.prop['es']['w'][:, :, x], ref['nets'][x][list(d.nname)].values)
		assert_allclose(d.prop['es']['w_in'][:, :, x], 2 * ref['nets'][x][list(d.nname)].values)
	pts, _ = d.linspace(0, 2, 5, 0.5)
	assert len(pts) == 5


def test_from_folders_dynamic_without_cell_types(tmp_path):
	"""Without data/subsets.txt a dynamic network loads with no cell type property."""
	data, work, fi_subsets, _ = make_dynamic_folders(str(tmp_path))
	os.remove(os.path.join(data, 'subsets.txt'))
	d = dynamic_network.from_folders(data, work, fi_subsets, dynamic=True)
	assert 'type' not in d.prop['c']
	assert d.prop['c']['coord'].shape == (2, d.cn)


def test_from_folders_dynamic_drops_failed_subset(tmp_path):
	"""A failed subset of a dynamic network is removed from the state points and state-state properties too."""
	data, work, fi_subsets, ref = make_dynamic_folders(str(tmp_path))
	for x in ['net_nweight.tsv.gz', 'net_inweight.tsv.gz']:
		os.remove(os.path.join(work, ref['snames'][6], x))
	d = dynamic_network.from_folders(data, work, fi_subsets, dynamic=True)
	assert_array_equal(d.sname, ref['snames'][:6])
	assert len(d.point['s']) == 6
	assert_array_equal(d.prop['ss']['traj-neighbor'], ref['nb'][:6, :6])
	assert d.prop['es']['w'].shape == (2, 5, 6)


# network.export


def test_export_full_network_and_cpm(tmp_path):
	"""export with no sparsity writes log2(CPM+1) and each state's non-empty rows and columns of edge weights."""
	d = make_static()
	out = os.path.join(str(tmp_path), 'out')
	d.export(out)
	assert sorted(os.listdir(out)) == ['Full', 'cpm.tsv.gz']
	cpm = read_tsv(os.path.join(out, 'cpm.tsv.gz'))
	assert list(cpm.index) == list(NNAME)
	assert list(cpm.columns) == list(SNAME)
	assert_allclose(cpm.values, np.log2(CPM + 1))
	assert sorted(os.listdir(os.path.join(out, 'Full'))) == ['s0.tsv.gz', 's1.tsv.gz', 's2.tsv.gz']
	w = static_w()
	s0 = read_tsv(os.path.join(out, 'Full', 's0.tsv.gz'))
	assert list(s0.index) == ['TF1', 'TF2'] and list(s0.columns) == list(NNAME)
	assert_allclose(s0.values, w[:, :, 0])
	s1 = read_tsv(os.path.join(out, 'Full', 's1.tsv.gz'))
	assert list(s1.index) == ['TF1', 'TF2'] and list(s1.columns) == ['TF2', 'G4']
	assert_allclose(s1.values, [[0.5, 0], [0, 0.6]])
	assert read_tsv(os.path.join(out, 'Full', 's2.tsv.gz')).shape == (0, 0)


def test_export_binary_networks(tmp_path):
	"""
	export with sparsities keeps the top int(sparsity*10) edges by absolute
	weight in each state, and writes no binary network for the state without edges.
	"""
	d = make_static()
	out = os.path.join(str(tmp_path), 'out')
	d.export(out, sparsities=[None, 0.2])
	assert sorted(os.listdir(out)) == ['Full', 'cpm.tsv.gz', 'sparsity=0.2']
	assert sorted(os.listdir(os.path.join(out, 'sparsity=0.2'))) == ['s0.tsv.gz', 's1.tsv.gz']
	s0 = read_tsv(os.path.join(out, 'sparsity=0.2', 's0.tsv.gz'))
	# Top 2 of 10 edges at s0: TF1->G3 (0.9) and TF2->TF1 (0.7)
	assert list(s0.index) == ['TF1', 'TF2'] and list(s0.columns) == ['TF1', 'G3']
	assert_array_equal(s0.values, [[False, True], [True, False]])
	s1 = read_tsv(os.path.join(out, 'sparsity=0.2', 's1.tsv.gz'))
	assert list(s1.columns) == ['TF2', 'G4']
	assert_array_equal(s1.values, [[True, False], [False, True]])


def test_export_single_float_sparsity_and_existing_folder(tmp_path):
	"""A float sparsity writes only its own folder, and an existing output folder raises FileExistsError."""
	d = make_static()
	out = os.path.join(str(tmp_path), 'out')
	d.export(out, sparsities=0.2)
	assert sorted(os.listdir(out)) == ['cpm.tsv.gz', 'sparsity=0.2']
	with pytest.raises(FileExistsError):
		d.export(out)


# dynamic_network.linspace


def test_linspace_points_and_smoother(star):
	"""
	linspace returns num even points on the state path between two nodes, and
	a smoother that Gaussian-averages states on that path, ignoring side branches.
	"""
	from dictys.net import stat
	pts, fsmooth = star.linspace(0, 2, 6, 0.5)
	assert len(pts) == 6
	assert pts.p.nn == star.sn
	assert_allclose(pts.p.lens, 0.5)
	start = pts.p.topoint()[[0]]
	assert_allclose((pts - start).ravel(), np.linspace(0, 2, 6), atol=1e-12)
	smoothed = fsmooth(stat.lcpm(star, cut=0))
	assert isinstance(smoothed, stat.fsmooth)
	w, _ = path_smoother(6, 0.5)
	assert_allclose(smoothed.compute(pts), np.log2(star.prop['ns']['cpm'][:, :5] + 1) @ w, rtol=1e-10)


def test_linspace_reverse_direction(star):
	"""Swapping start and stop reverses the points along the path."""
	pts, fsmooth = star.linspace(2, 0, 5, 0.5)
	start = pts.p.topoint()[[0]]
	assert_allclose((pts - start).ravel(), np.linspace(2, 0, 5), atol=1e-12)


def test_linspace_check_rejects_non_branching_node():
	"""With check=True a start node of degree 2 raises ValueError; check=False allows it."""
	d = build_dynamic(**CHAIN)
	with pytest.raises(ValueError, match='not a branching or terminal node'):
		d.linspace(1, 2, 3, 0.5)
	with pytest.raises(ValueError, match='not a branching or terminal node'):
		d.linspace(0, 1, 3, 0.5, check=True)
	pts, _ = d.linspace(1, 2, 3, 0.5, check=False)
	start = pts.p.topoint()[[1]]
	assert_allclose((pts - start).ravel(), [0, 0.5, 1], atol=1e-12)


# dynamic_network.compute_chars


def test_compute_chars_expression(star):
	"""Expression characteristics match those of independently smoothed log2(CPM+1) curves."""
	ans = star.compute_chars(0, 2, num=7, dist=0.4, mode='expression')
	w, t = path_smoother(7, 0.4)
	y = np.log2(star.prop['ns']['cpm'][:, :5] + 1) @ w
	expected = curve_chars(y, t)
	assert list(ans.index) == list(NNAME)
	assert list(ans.columns) == ['Terminal logFC', 'Transient logFC', 'Switching time']
	assert_allclose(ans.values, expected[ans.columns].values, rtol=1e-8, atol=1e-10)


def test_compute_chars_regulation(star):
	"""Regulation characteristics use log2(1 + target count) of the top int(sparsity*10) smoothed edges."""
	ans = star.compute_chars(0, 2, num=6, dist=0.6, mode='regulation', sparsity=0.3)
	w, t = path_smoother(6, 0.6)
	net = star.prop['es']['w'][:, :, :5] @ w
	y = np.log2(binarize_top(net, 3).sum(axis=1) + 1)
	assert list(ans.index) == ['TF1', 'TF2']
	assert_allclose(ans.values, curve_chars(y, t)[ans.columns].values, atol=1e-10)


def test_compute_chars_defaults_and_unknown_mode(star):
	"""compute_chars defaults to 100 points, distance 1.5 and regulation mode, and rejects an unknown mode."""
	ans = star.compute_chars(0, 2, sparsity=0.4)
	w, t = path_smoother(100, 1.5)
	net = star.prop['es']['w'][:, :, :5] @ w
	y = np.log2(binarize_top(net, 4).sum(axis=1) + 1)
	assert_allclose(ans.values, curve_chars(y, t)[ans.columns].values, atol=1e-10)
	with pytest.raises(ValueError, match='Unknown mode'):
		star.compute_chars(0, 2, num=5, mode='other')


# dynamic_network.draw_discover


def test_draw_discover_expression(star):
	"""draw_discover draws, per requested pattern, the top gene's normalised curve against pseudotime."""
	pytest.importorskip('matplotlib')
	ans = star.draw_discover(0, 2, num=6, dist=0.5, mode='expression', sparsity=0.5, ntops=(1, 1, 0, 0))
	assert len(ans) == 2
	w, t = path_smoother(6, 0.5)
	y = np.log2(star.prop['ns']['cpm'][:, :5] + 1) @ w
	terminal = y[:, -1] - y[:, 0]
	shown = [terminal.argmax(), terminal.argmin()]
	vmax = np.abs(terminal[shown]).max()
	for (fig, axes, cmap), k in zip(ans, shown):
		assert len(fig.axes) == len(axes) == 4
		assert {axes[x].get_yticklabels()[0].get_text() for x in [1, 2]} == {'Switching time', 'Terminal logFC'}
		assert axes[3].get_xlabel() == NNAME[k]
		line = axes[3].lines[0]
		assert_allclose(line.get_xdata(), t, atol=1e-12)
		assert_allclose(line.get_ydata(), (y[k] - y[k].min()) / (y[k].max() - y[k].min()), atol=1e-10)
		assert set(cmap) == {'Switching time', 'Terminal logFC'}
		assert_allclose(cmap['Terminal logFC'].norm.vmax, vmax, rtol=1e-10)
		assert_allclose(cmap['Terminal logFC'].norm.vmin, -vmax, rtol=1e-10)


@pytest.mark.parametrize('mode', ['regulation', 'weighted_regulation'])
def test_draw_discover_regulation_modes(star, mode):
	"""In regulation modes draw_discover draws four patterns of both regulators' curves."""
	pytest.importorskip('matplotlib')
	ans = star.draw_discover(0, 2, num=5, dist=0.5, mode=mode, sparsity=0.3)
	assert len(ans) == 4
	w, t = path_smoother(5, 0.5)
	net = star.prop['es']['w'][:, :, :5] @ w
	if mode == 'regulation':
		deg = binarize_top(net, 3).sum(axis=1)
	else:
		a = np.abs(net)
		deg = (a * (0.3 / a.mean(axis=(0, 1)))).sum(axis=1)
	y = np.log2(deg + 1)
	for fig, axes, cmap in ans:
		assert len(axes) == 5
		names = [x.get_xlabel() for x in axes[3:]]
		assert sorted(names) == ['TF1', 'TF2']
		rows = y[[list(NNAME).index(x) for x in names]]
		for ax, row in zip(axes[3:], rows):
			assert_allclose(ax.lines[0].get_ydata(), (row - rows.min()) / (rows.max() - rows.min()), atol=1e-10)


def test_draw_discover_unknown_mode(star):
	"""draw_discover rejects an unknown mode."""
	pytest.importorskip('matplotlib')
	with pytest.raises(ValueError, match='Unknown mode'):
		star.draw_discover(0, 2, num=5, mode='other')


# dynamic_network.draw_regulation_heatmap


def test_draw_regulation_heatmap_default_colormap(star):
	"""The heatmap shows smoothed edge strengths on a symmetric colormap capped at the 95% quantile."""
	pytest.importorskip('matplotlib')
	regs = [('TF1', 'G3'), ('TF2', 'G4'), ('TF2', 'TF1')]
	fig, ax, cmap = star.draw_regulation_heatmap(0, 2, regs, num=8, dist=0.5)
	w, _ = path_smoother(8, 0.5)
	net = star.prop['es']['w'][:, :, :5] @ w
	dnet = net[[0, 1, 1], [2, 3, 0]]
	vmax = np.quantile(np.abs(dnet).ravel(), 0.95)
	assert_allclose(cmap.norm.vmax, vmax, rtol=1e-10)
	assert_allclose(cmap.norm.vmin, -vmax, rtol=1e-10)
	assert cmap.get_cmap().name == 'coolwarm'
	assert_allclose(ax.images[0].get_array(), cmap.to_rgba(dnet), atol=1e-12)
	assert [x.get_text() for x in ax.get_yticklabels()] == ['TF1-G3', 'TF2-G4', 'TF2-TF1']
	assert_allclose(fig.get_size_inches(), [2, 0.22 * 3])


def test_draw_regulation_heatmap_options(star):
	"""vmax, a named colormap, figsize, an existing axes and a ScalarMappable colormap are honoured."""
	mpl = pytest.importorskip('matplotlib')
	import matplotlib.pyplot as plt
	regs = [('TF1', 'G5')]
	fig, ax, cmap = star.draw_regulation_heatmap(0, 2, regs, num=4, dist=1.0, cmap='viridis', figsize=(3, 0.5), vmax=2.0)
	assert (cmap.norm.vmin, cmap.norm.vmax) == (-2.0, 2.0)
	assert cmap.get_cmap().name == 'viridis'
	assert_allclose(fig.get_size_inches(), [3, 0.5])
	fig2, ax2 = plt.subplots(figsize=(4, 1))
	sm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=-1, vmax=1), cmap='bwr')
	fig3, ax3, cmap3 = star.draw_regulation_heatmap(0, 2, regs, num=4, dist=1.0, ax=ax2, cmap=sm, figsize=None)
	assert fig3 is fig2 and ax3 is ax2 and cmap3 is sm
	w, _ = path_smoother(4, 1.0)
	dnet = (star.prop['es']['w'][:, :, :5] @ w)[[0], [4]]
	assert_allclose(ax2.images[0].get_array(), sm.to_rgba(dnet), atol=1e-12)


def test_draw_regulation_heatmap_errors(star):
	"""Unknown genes, figsize with ax, and vmax with a ScalarMappable raise ValueError."""
	mpl = pytest.importorskip('matplotlib')
	import matplotlib.pyplot as plt
	with pytest.raises(ValueError, match='Regulator gene'):
		star.draw_regulation_heatmap(0, 2, [('G3', 'TF1')], num=3)
	with pytest.raises(ValueError, match='Target gene'):
		star.draw_regulation_heatmap(0, 2, [('TF1', 'X')], num=3)
	_, ax = plt.subplots()
	with pytest.raises(ValueError, match='figsize'):
		star.draw_regulation_heatmap(0, 2, [('TF1', 'G3')], num=3, ax=ax)
	sm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=-1, vmax=1), cmap='bwr')
	with pytest.raises(ValueError, match='vmax'):
		star.draw_regulation_heatmap(0, 2, [('TF1', 'G3')], num=3, cmap=sm, vmax=1.0)


# dynamic_network.export


def test_dynamic_export(tmp_path, star):
	"""export writes pseudotime, smoothed log2(CPM+1) and per-point full and binarized networks."""
	out = os.path.join(str(tmp_path), 'out')
	star.export(out, 0, 2, 4, 0.5, sparsities=[None, 0.5])
	assert sorted(os.listdir(out)) == ['Full', 'cpm.tsv.gz', 'pseudotime.tsv.gz', 'sparsity=0.5']
	w, t = path_smoother(4, 0.5)
	pt = pd.read_csv(os.path.join(out, 'pseudotime.tsv.gz'), sep='\t')
	assert list(pt.columns) == ['State_index', 'Pseudo-time']
	assert_array_equal(pt['State_index'], [1, 2, 3, 4])
	assert_allclose(pt['Pseudo-time'], t, atol=1e-12)
	cpm = read_tsv(os.path.join(out, 'cpm.tsv.gz'))
	assert list(cpm.index) == list(NNAME)
	assert list(cpm.columns) == ['1', '2', '3', '4']
	assert_allclose(cpm.values, np.log2(star.prop['ns']['cpm'][:, :5] + 1) @ w, rtol=1e-10)
	net = star.prop['es']['w'][:, :, :5] @ w
	binary = binarize_top(net, 5)
	for k in range(4):
		full = read_tsv(os.path.join(out, 'Full', f'{k + 1}.tsv.gz'))
		assert list(full.index) == ['TF1', 'TF2'] and list(full.columns) == list(NNAME)
		assert_allclose(full.values, net[:, :, k], rtol=1e-10)
		cols = binary[:, :, k].any(axis=0)
		rows = binary[:, :, k].any(axis=1)
		assert not binary[:, :, k][rows][:, cols].all()
		b = read_tsv(os.path.join(out, 'sparsity=0.5', f'{k + 1}.tsv.gz'))
		assert list(b.columns) == list(NNAME[cols])
		assert_array_equal(b.values, binary[:, :, k][rows][:, cols])


def test_dynamic_export_float_sparsity_and_existing_folder(tmp_path, star):
	"""A float sparsity writes only its own network folder; an existing output folder raises FileExistsError."""
	out = os.path.join(str(tmp_path), 'out')
	star.export(out, 0, 2, 3, 0.5, sparsities=0.5)
	assert sorted(os.listdir(out)) == ['cpm.tsv.gz', 'pseudotime.tsv.gz', 'sparsity=0.5']
	assert sorted(os.listdir(os.path.join(out, 'sparsity=0.5'))) == ['1.tsv.gz', '2.tsv.gz', '3.tsv.gz']
	with pytest.raises(FileExistsError):
		star.export(out, 0, 2, 3, 0.5)
