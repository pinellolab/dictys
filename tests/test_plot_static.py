"""
Tests for dictys.plot.static: compute_spec, compute_reg_spec,
fig_heatmap_reg_spec, fig_heatmap_top, fig_diff_scatter, fig_diff_rank and
fig_subnet, on a small hand-built dictys.net.network.
"""

import logging
import os

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose

matplotlib = pytest.importorskip('matplotlib')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import ticker  # noqa: E402
from matplotlib.collections import QuadMesh  # noqa: E402
from matplotlib.colors import ListedColormap, to_rgba, to_rgba_array  # noqa: E402

from dictys.net import network  # noqa: E402
from dictys.plot import static  # noqa: E402

REGS = ['TF1', 'TF2', 'TF3', 'TF4']
TARGETS = ['TG1', 'TG2', 'TG3', 'TG4', 'TG5', 'TG6']
STATES = ['S1', 'S2', 'S3']


@pytest.fixture(autouse=True)
def close_figures():
	"""Close every matplotlib figure a test opened."""
	yield
	plt.close('all')


@pytest.fixture
def adjusttext():
	"""Skip unless adjustText, imported unconditionally by fig_diff_scatter and fig_diff_rank, is installed."""
	return pytest.importorskip('adjustText')


def _weights():
	"""
	Edge weights (regulator, target, state). Every edge has a small distinct
	positive background weight; at sparsity 0.25 (6 of 24 edges per state) only
	the designed strong edges are positive:
	S1: TF1->TG1..TG6 (10..15); S2: TF2->TG1..TG4 (10..13), TF4->TG5,TG6 (-20,-21);
	S3: TF3->TG1..TG3 (10..12), TF4->TG4..TG6 (30..32).
	"""
	w = 0.01 * (1 + np.arange(72).reshape(4, 6, 3) / 100)
	w[0, :, 0] = 10 + np.arange(6)
	w[1, :4, 1] = 10 + np.arange(4)
	w[3, 4:, 1] = -(20 + np.arange(2))
	w[2, :3, 2] = 10 + np.arange(3)
	w[3, 3:, 2] = 30 + np.arange(3)
	return w


# CPM of each gene (rows: TF1..TF4, TG1..TG6) in S1..S3
CPM = np.array([
	[5, 40, 0.5],
	[3, 3, 3],
	[0.2, 0.3, 0.1],
	[7, 1, 0],
	[100, 0, 0],
	[1, 2, 3],
	[1, 2, 3],
	[1, 2, 3],
	[1, 2, 3],
	[1, 2, 3],
], dtype=float)


def _network(w, cpm=CPM, mask=None):
	"""Static network over sorted gene names TF1..TF4 (regulators) and TG1..TG6 (targets)."""
	if mask is None:
		mask = np.ones(w.shape, dtype=bool)
	return network(
		cname=np.array(['c%d' % x for x in range(6)]),
		sname=np.array(STATES),
		nname=np.array(REGS + TARGETS),
		nids=[np.arange(4), np.arange(4, 10)],
		esprop={'w': w, 'mask': mask},
		nsprop={'cpm': cpm})


@pytest.fixture
def d0():
	"""The designed static network."""
	return _network(_weights())


def _row_normalize(v):
	s = v.sum(axis=1, keepdims=True)
	return np.where(s > 0, v / np.where(s > 0, s, 1), 0)


########################################################################
# compute_spec and compute_reg_spec
########################################################################

def test_compute_spec_tables(d0):
	"""Target counts, their specificity, CPM and CPM specificity match hand-derived tables."""
	deg, deg_spec, cpm, cpm_spec, _, _ = static.compute_spec(d0, sparsity=0.25, vmin=1)
	expected_deg = np.array([[6, 0, 0], [0, 4, 0], [0, 0, 3], [0, 2, 3]])
	for x in [deg, deg_spec]:
		assert list(x.index) == REGS and list(x.columns) == STATES
	assert_allclose(deg.values, expected_deg)
	# Every target is tested, so the rate is count/6 and specificity is count/row sum
	assert_allclose(deg_spec.values, [[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0.4, 0.6]])
	for x in [cpm, cpm_spec]:
		assert list(x.index) == REGS + TARGETS and list(x.columns) == STATES
	assert_allclose(cpm.values, CPM, rtol=1e-12, atol=1e-12)
	assert_allclose(cpm_spec.values, CPM / CPM.sum(axis=1, keepdims=True), rtol=1e-12)


@pytest.mark.parametrize('ka,markers,states', [
	# TF4's specificity (0, 0.4, 0.6) has entropy 0.673 > 0.5*ln(3)
	({'vmin': 1}, ['TF1', 'TF2', 'TF3'], ['S1', 'S2', 'S3']),
	# ... but passes 0.7*ln(3)=0.769, and ranks after TF3 (specificity 1) in S3
	({'vmin': 1, 'min_entropy': 0.7}, ['TF1', 'TF2', 'TF3', 'TF4'], ['S1', 'S2', 'S3', 'S3']),
	({'vmin': 1, 'min_entropy': 0.7, 'nmax': 1}, ['TF1', 'TF2', 'TF3'], ['S1', 'S2', 'S3']),
	({'vmin': 1, 'min_entropy': 0.7, 'ncut': 0.7}, ['TF1', 'TF2', 'TF3'], ['S1', 'S2', 'S3']),
	({'vmin': 1, 'min_entropy': 0.7, 'nmax': None}, ['TF1', 'TF2', 'TF3', 'TF4'], ['S1', 'S2', 'S3', 'S3']),
	# TF3 and TF4 have 3 targets only
	({'vmin': 4, 'min_entropy': 0.7}, ['TF1', 'TF2'], ['S1', 'S2']),
])
def test_compute_spec_regulation_markers(d0, ka, markers, states):
	"""Regulation markers follow the entropy, ncut, vmin and nmax filters, ordered by state then specificity."""
	ans = static.compute_spec(d0, sparsity=0.25, **ka)
	assert list(ans[4]) == markers
	assert list(ans[5]) == states


@pytest.mark.parametrize('ka,markers,states', [
	# TG1 (100,0,0) and TF1 (5,40,0.5) are specific; TF4 (7,1,0) has CPM 7 < 20
	({}, ['TG1', 'TF1'], ['S1', 'S2']),
	({'vmin': 5}, ['TG1', 'TF4', 'TF1'], ['S1', 'S1', 'S2']),
	({'vmin': 5, 'nmax': 1}, ['TG1', 'TF1'], ['S1', 'S2']),
])
def test_compute_spec_expression_markers(d0, ka, markers, states):
	"""base='cpm' selects expression markers among all genes from CPM specificity."""
	ans = static.compute_spec(d0, base='cpm', sparsity=0.25, **ka)
	assert list(ans[4]) == markers
	assert list(ans[5]) == states


def test_compute_spec_select_state(d0):
	"""select_state restricts and orders the states, and specificity is computed among them only."""
	deg, deg_spec, cpm, cpm_spec, markers, states = static.compute_spec(d0, sparsity=0.25, vmin=1, select_state=['S3', 'S1'])
	for x in [deg, deg_spec, cpm, cpm_spec]:
		assert list(x.columns) == ['S3', 'S1']
	assert_allclose(deg.values, [[0, 6], [0, 0], [3, 0], [3, 0]])
	assert_allclose(deg_spec.values, [[0, 1], [0, 0], [1, 0], [1, 0]])
	assert_allclose(cpm.values, CPM[:, [2, 0]], atol=1e-12)
	assert_allclose(cpm_spec.values, _row_normalize(CPM[:, [2, 0]]), atol=1e-12)
	assert list(markers) == ['TF3', 'TF4', 'TF1']
	assert list(states) == ['S3', 'S3', 'S1']


def test_compute_spec_weighted(d0):
	"""weighted=True uses weighted outdegree scaled so each state's total equals sparsity*edge count."""
	w = np.abs(_weights())
	deg, deg_spec, _, _, markers, states = static.compute_spec(d0, sparsity=0.25, vmin=1, weighted=True)
	expected = 0.25 * 24 * w.sum(axis=1) / w.sum(axis=(0, 1))
	assert_allclose(deg.values, expected, rtol=1e-10)
	assert_allclose(deg.values.sum(axis=0), [0.25 * 24] * 3, rtol=1e-10)
	assert_allclose(deg_spec.values, _row_normalize(w.sum(axis=1) / 6), rtol=1e-10)
	# TF4's weighted specificity (~0, 0.31, 0.69) still fails min_entropy=0.5
	assert list(markers) == ['TF1', 'TF2', 'TF3']
	assert list(states) == ['S1', 'S2', 'S3']


def test_compute_spec_partial_mask_reference(rng):
	"""With untested edges, binarization keeps sparsity*tested edges and the rate divides by tested targets (numpy reference)."""
	w = rng.normal(size=(4, 6, 3))
	mask = rng.random(size=w.shape) < 0.7
	mask[:, 5, 1] = False
	w[~mask] = 0
	d = _network(w, mask=mask)
	sparsity = 0.3
	deg, deg_spec, _, _, _, _ = static.compute_spec(d, sparsity=sparsity, vmin=0, min_entropy=1)
	exp_deg = np.zeros((4, 3))
	exp_rate = np.zeros((4, 3))
	for xi in range(3):
		a = np.abs(w[:, :, xi])
		k = int(sparsity * mask[:, :, xi].sum())
		thr = np.sort(a.ravel())[::-1][k - 1]
		b = a >= thr
		exp_deg[:, xi] = b.sum(axis=1)
		exp_rate[:, xi] = (b & mask[:, :, xi]).sum(axis=1) / mask[:, :, xi].any(axis=0).sum()
	assert_allclose(deg.values, exp_deg)
	assert_allclose(deg_spec.values, _row_normalize(exp_rate), rtol=1e-12)


def test_compute_spec_unknown_base_raises(d0):
	"""An unknown base raises ValueError."""
	with pytest.raises(ValueError):
		static.compute_spec(d0, base='tpm', sparsity=0.25, vmin=1)


def test_compute_reg_spec_equals_compute_spec(d0):
	"""compute_reg_spec returns compute_spec's regulation-marker results and fixes base."""
	ka = dict(sparsity=0.25, vmin=1, min_entropy=0.7)
	a = static.compute_reg_spec(d0, **ka)
	b = static.compute_spec(d0, base='deg', **ka)
	for x, y in zip(a[:4], b[:4]):
		pd.testing.assert_frame_equal(x, y)
	assert list(a[4]) == ['TF1', 'TF2', 'TF3', 'TF4']
	assert list(a[5]) == list(b[5])
	with pytest.raises(TypeError):
		static.compute_reg_spec(d0, base='cpm', **ka)


########################################################################
# fig_heatmap_reg_spec
########################################################################

def _bhattacharyya(u, v):
	return 1 - np.sqrt(u * v).sum()


@pytest.fixture
def spec():
	"""Specificity of 12 regulators (rows) across 3 states (columns): 4 groups specific to S1, S2, S3 or shared."""
	base = {
		'S1': [0.9, 0.05, 0.05],
		'S2': [0.05, 0.9, 0.05],
		'S3': [0.05, 0.05, 0.9],
		'all': [0.34, 0.33, 0.33],
	}
	rows = []
	names = []
	for xi, x in enumerate(['S1', 'S2', 'S3', 'all']):
		for xj in range(3):
			v = np.array(base[x]) + 0.01 * (xj + 1) * np.array([1, -0.5, -0.5])[np.roll(np.arange(3), xi)]
			rows.append(v / v.sum())
			names.append(f'R{xi * 3 + xj + 1}')
	return pd.DataFrame(np.array(rows), index=names, columns=STATES)


def test_fig_heatmap_reg_spec_layout(spec):
	"""Regulators are clustered by Bhattacharyya distance into columns, states are rows, and only g_ann regulators are labelled."""
	from scipy.cluster.hierarchy import dendrogram, linkage
	fig = static.fig_heatmap_reg_spec(spec, g_ann=['R2', 'R7', 'missing'])
	leaves = dendrogram(linkage(spec.values, method='average', metric=_bhattacharyya, optimal_ordering=False), no_plot=True)['leaves']
	order = list(spec.index[leaves])
	assert len(fig.axes) == 3
	ax = fig.axes[1]
	im = ax.images[0]
	assert_allclose(np.asarray(im.get_array()), spec.T.loc[:, order].values)
	assert im.get_cmap().name == 'viridis'
	assert_allclose(im.get_clim(), [spec.values.min(), spec.values.max()])
	assert [x.get_text() for x in ax.get_yticklabels()] == STATES
	assert list(ax.get_xticks()) == [order.index('R2'), order.index('R7')]
	assert [x.get_text() for x in ax.get_xticklabels()] == ['R2', 'R7']
	# The 4 groups of 3 regulators are leaves next to each other
	groups = [{order.index(f'R{x * 3 + y + 1}') for y in range(3)} for x in range(4)]
	assert all(max(x) - min(x) == 2 for x in groups)


def test_fig_heatmap_reg_spec_size_and_overrides(spec):
	"""Cells are figscale*aspect wide and figscale tall; extra keywords override the heatmap defaults."""
	fig = static.fig_heatmap_reg_spec(spec, aspect=0.5, figscale=0.2, g_ann=[], cmap='magma', dtop=0)
	assert len(fig.axes) == 2
	ax = fig.axes[0]
	pos = ax.get_position()
	size = fig.get_size_inches()
	assert_allclose([pos.width * size[0], pos.height * size[1]], [12 * 0.2 * 0.5, 3 * 0.2])
	assert ax.images[0].get_cmap().name == 'magma'
	assert_allclose(np.asarray(ax.images[0].get_array()), spec.T.values)
	assert list(ax.get_xticks()) == []


########################################################################
# fig_heatmap_top
########################################################################

def _expected_top(w, selection, ntop, normalize=True):
	"""Reference: top |w| targets per pair (descending), union in order, each pair scaled by its max |w|."""
	order = []
	cols = {}
	for reg, state in selection:
		v = w[REGS.index(reg), :, STATES.index(state)]
		for x in np.argsort(-np.abs(v))[:ntop]:
			if TARGETS[x] not in order:
				order.append(TARGETS[x])
		cols[' ' + reg + '-' + state] = v / (np.abs(v).max() if normalize else 1)
	return pd.DataFrame({x: [y[TARGETS.index(z)] for z in order] for x, y in cols.items()}, index=order)


def test_fig_heatmap_top_default(d0):
	"""Top targets of each pair are stacked in a column-normalized heatmap under a bar of state colors."""
	sel = [('TF1', 'S1'), ('TF2', 'S2'), ('TF4', 'S3')]
	fig, fig2, net = static.fig_heatmap_top(d0, sel, ntop=3, gann=['TG1', 'TG6'])
	expected = _expected_top(_weights(), sel, 3)
	assert list(expected.index) == ['TG6', 'TG5', 'TG4', 'TG3', 'TG2']
	pd.testing.assert_frame_equal(net, expected, check_exact=False, rtol=1e-12)
	assert len(fig.axes) == 2
	ax_bar, ax_heat = fig.axes
	assert_allclose(np.asarray(ax_bar.images[0].get_array())[0], to_rgba_array(plt.get_cmap('tab10').colors[:3]))
	assert [x.get_text() for x in ax_bar.get_xticklabels()] == [' TF1-S1', ' TF2-S2', ' TF4-S3']
	im = ax_heat.images[0]
	assert_allclose(np.asarray(im.get_array()), expected.values)
	assert im.get_cmap().name == 'coolwarm'
	assert_allclose(im.get_clim(), [-1, 1])
	assert list(ax_heat.get_yticks()) == [0]
	assert [x.get_text() for x in ax_heat.get_yticklabels()] == ['TG6']
	# Width figscale per pair; height figscale*(aspect per target + topheight + topspace)
	assert_allclose(fig.get_size_inches(), [0.15 * 3, 0.15 * 0.2 * 5 + 0.15 * (0.7 + 0.3)])
	cax = fig2.axes[0]
	assert cax.get_title() == 'Relative' + os.linesep + 'strength'
	assert_allclose(cax.get_ylim(), [-1, 1])
	assert [x for x in cax.collections if isinstance(x, QuadMesh)][0].get_cmap().name == 'coolwarm'


def test_fig_heatmap_top_shared_state_color(d0):
	"""Pairs from the same state share one bar color, assigned in order of first appearance."""
	sel = [('TF4', 'S3'), ('TF3', 'S3'), ('TF2', 'S2')]
	fig, _, net = static.fig_heatmap_top(d0, sel, ntop=2, cmap_type='Set1')
	colors = to_rgba_array(plt.get_cmap('Set1').colors[:2])
	assert_allclose(np.asarray(fig.axes[0].images[0].get_array())[0], colors[[0, 0, 1]])
	pd.testing.assert_frame_equal(net, _expected_top(_weights(), sel, 2), check_exact=False, rtol=1e-12)


def test_fig_heatmap_top_repression_only(d0, caplog):
	"""direction=-1 keeps only negative targets, scales the strongest to -1 and warns when fewer than ntop exist."""
	sel = [('TF4', 'S2')]
	fig, _, net = static.fig_heatmap_top(d0, sel, ntop=2, direction=-1)
	assert list(net.index) == ['TG6', 'TG5']
	assert_allclose(net.values[:, 0], [-1, -20 / 21])
	assert_allclose(fig.axes[1].images[0].get_clim(), [-1, 1])
	plt.close('all')
	with caplog.at_level(logging.WARNING):
		_, _, net = static.fig_heatmap_top(d0, sel, ntop=3, direction=-1)
	assert any('<3 targets' in x.getMessage() for x in caplog.records)
	assert list(net.index[:2]) == ['TG6', 'TG5']


def test_fig_heatmap_top_activation_only(d0):
	"""direction=1 ranks targets by positive weight only and sets the color range from the largest positive value."""
	sel = [('TF4', 'S2')]
	_, _, net = static.fig_heatmap_top(d0, sel, ntop=2, direction=1, normalization='none')
	w = _weights()[3, :, 1]
	assert list(net.index) == ['TG4', 'TG3']
	assert_allclose(net.values[:, 0], w[[3, 2]])


@pytest.mark.parametrize('normalization', ['none', None])
def test_fig_heatmap_top_no_normalization(d0, normalization):
	"""normalization 'none' or None keeps raw weights and the symmetric color range is the max |weight|."""
	sel = [('TF1', 'S1'), ('TF4', 'S2')]
	fig, _, net = static.fig_heatmap_top(d0, sel, ntop=2, normalization=normalization)
	expected = _expected_top(_weights(), sel, 2, normalize=False)
	pd.testing.assert_frame_equal(net, expected, check_exact=False, rtol=1e-12)
	assert_allclose(fig.axes[1].images[0].get_clim(), [-21, 21])


def test_fig_heatmap_top_options(d0):
	"""gann='all', cmap_value, aspect, topheight, topspace and figscale set labels, colormaps and figure size."""
	sel = [('TF1', 'S1'), ('TF2', 'S2')]
	fig, fig2, net = static.fig_heatmap_top(d0, sel, ntop=2, gann='all', cmap_value='PiYG', aspect=0.5, topheight=1, topspace=0.5, figscale=0.4)
	ax_bar, ax_heat = fig.axes
	assert [x.get_text() for x in ax_heat.get_yticklabels()] == list(net.index)
	assert ax_heat.images[0].get_cmap().name == 'PiYG'
	assert [x for x in fig2.axes[0].collections if isinstance(x, QuadMesh)][0].get_cmap().name == 'PiYG'
	height = 0.4 * 0.5 * len(net) + 0.4 * 1.5
	assert_allclose(fig.get_size_inches(), [0.4 * 2, height])
	assert_allclose(ax_bar.get_position().height, 0.4 * 1 / height)


@pytest.mark.parametrize('sel,ka,error', [
	([('XX', 'S1')], {}, ValueError),
	([('TF1', 'S9')], {}, ValueError),
	([('TF1', 'S1')], {'direction': 2}, ValueError),
	([('TF1', 'S1')], {'direction': -1}, ValueError),
	([('TF1', 'S1')], {'normalization': 'row'}, ValueError),
	([('TF1', 'S1')], {'cmap_type': None}, ValueError),
	([('TF1', 'S1'), ('TF2', 'S2'), ('TF3', 'S3')], {'cmap_type': ListedColormap(['r', 'g'])}, RuntimeError),
])
def test_fig_heatmap_top_invalid_input_raises(d0, sel, ka, error):
	"""Unknown TFs or states, bad direction/normalization, no targets in the direction or too few colors raise."""
	with pytest.raises(error):
		static.fig_heatmap_top(d0, sel, ntop=2, **ka)


########################################################################
# fig_diff_scatter
########################################################################

def _log2fc(a, b):
	return np.log2(b + 1) - np.log2(a + 1)


def _sym_lim(*vals):
	"""Symmetric shared limit: max |.| over each axis' finite [min,max] expanded by 2% of its range."""
	ans = []
	for v in vals:
		v = np.asarray(v, dtype=float)
		v = v[np.isfinite(v)]
		r = v.max() - v.min()
		ans += [v.min() - 0.02 * r, v.max() + 0.02 * r]
	return np.abs(ans).max()


def test_fig_diff_scatter_values(d0, adjusttext):
	"""Returned and plotted logFCs of CPM and target count from S1 to S2 match hand-derived values, NaN rows dropped."""
	fig = plt.figure()
	ax = fig.add_subplot(111)
	data = static.fig_diff_scatter(d0, ax, ('S1', 'S2'), annotate=['TF1', 'TF2', 'TF3'], sparsity=0.25, cut_ntarget=1, ka_adjust_text=None)
	# CPM: TF3 is below 1 CPM in both states; targets: TF3 has none in both
	de = {'TF1': _log2fc(5, 40), 'TF2': 0.0, 'TF4': _log2fc(7, 1)}
	dr = {'TF1': _log2fc(6, 0), 'TF2': _log2fc(0, 4), 'TF4': _log2fc(0, 2)}
	assert list(data.columns) == ['DE_logFC', 'DR_logFC']
	assert list(data.index) == ['TF1', 'TF2', 'TF4']
	assert_allclose(data['DE_logFC'].values, [de[x] for x in data.index], atol=1e-12)
	assert_allclose(data['DR_logFC'].values, [dr[x] for x in data.index], atol=1e-12)
	coll = ax.collections[0]
	offsets = np.ma.filled(np.ma.asarray(coll.get_offsets()).astype(float), np.nan)
	assert offsets.shape == (4, 2)
	assert_allclose(offsets[[0, 1, 3]], data.values, atol=1e-12)
	assert np.isnan(offsets[2]).any()
	assert_allclose(coll.get_sizes(), [20])
	assert coll.get_alpha() == 0.7
	texts = {x.get_text(): x.get_position() for x in ax.texts if x.get_text() != ''}
	assert set(texts) == {'TF1', 'TF2'}
	for x in texts:
		assert_allclose(texts[x], data.loc[x].values, atol=1e-12)
	lim = _sym_lim(list(de.values()), list(dr.values()))
	assert_allclose(ax.get_xlim(), [-lim, lim])
	assert_allclose(ax.get_ylim(), [-lim, lim])
	assert ax.get_xlabel() == 'Differential expression logFC'
	assert ax.get_ylabel() == 'Differential regulation logFC'
	assert ax.get_aspect() == 1
	assert isinstance(ax.xaxis.get_major_locator(), ticker.MultipleLocator)
	# axes_alpha=0.4 draws the x=0 and y=0 lines
	assert len(ax.lines) == 2
	assert_allclose(ax.lines[0].get_ydata(), [0, 0])
	assert_allclose(ax.lines[1].get_xdata(), [0, 0])
	assert all(x.get_alpha() == 0.4 for x in ax.lines)


def test_fig_diff_scatter_default_cuts(d0, adjusttext):
	"""With the default 5-target cut only TF1 (6 targets in S1) has a differential regulation value."""
	fig = plt.figure()
	ax = fig.add_subplot(111)
	data = static.fig_diff_scatter(d0, ax, ('S1', 'S2'), sparsity=0.25)
	assert list(data.index) == ['TF1']
	assert_allclose(data.values, [[_log2fc(5, 40), _log2fc(6, 0)]], atol=1e-12)
	assert [x.get_text() for x in ax.texts] == []


def test_fig_diff_scatter_zero_cuts_and_reversed_states(d0, adjusttext):
	"""Lower cut_cpm and cut_ntarget keep TF3; swapping the states negates every logFC."""
	fig = plt.figure()
	ax = fig.add_subplot(111)
	data = static.fig_diff_scatter(d0, ax, ('S2', 'S1'), sparsity=0.25, cut_cpm=0.25, cut_ntarget=0)
	assert list(data.index) == REGS
	assert_allclose(data.loc['TF3'].values, [_log2fc(0.3, 0.2), 0], atol=1e-12)
	assert_allclose(data.loc['TF1'].values, [_log2fc(40, 5), _log2fc(0, 6)], atol=1e-12)


def test_fig_diff_scatter_weighted_and_options(d0, adjusttext):
	"""weighted=True uses scaled weighted outdegree; axes_alpha=1, aspect, lim=set() and scatter keywords take effect."""
	w = np.abs(_weights())
	deg = 0.25 * 24 * w.sum(axis=1) / w.sum(axis=(0, 1))
	fig = plt.figure()
	ax = fig.add_subplot(111)
	data = static.fig_diff_scatter(d0, ax, ('S1', 'S2'), annotate=['TF1', 'TF2'], weighted=True, sparsity=0.25, cut_ntarget=1, axes_alpha=1, aspect=2, lim=set(), ka_adjust_text=None, s=7)
	dr = np.log2(deg[:, 1] + 1) - np.log2(deg[:, 0] + 1)
	keep = (deg[:, :2] >= 1).any(axis=1) & np.array([True, True, False, True])
	assert list(data.index) == list(np.array(REGS)[keep])
	assert_allclose(data['DR_logFC'].values, dr[keep], rtol=1e-10)
	assert len(ax.lines) == 0
	assert ax.get_aspect() == 2
	assert_allclose(ax.collections[0].get_sizes(), [7])
	# lim=set(): each axis uses its own data range expanded by 2%
	for lim, v in [(ax.get_xlim(), data['DE_logFC'].values), (ax.get_ylim(), dr[[0, 1, 3]])]:
		r = v.max() - v.min()
		assert_allclose(lim, [v.min() - 0.02 * r, v.max() + 0.02 * r], rtol=1e-10)


def test_fig_diff_scatter_annotate_all_with_adjust_text(d0, adjusttext):
	"""annotate='all' labels every gene shown, and adjustText runs on the labels."""
	fig = plt.figure()
	ax = fig.add_subplot(111)
	data = static.fig_diff_scatter(d0, ax, ('S1', 'S2'), annotate='all', sparsity=0.25, cut_ntarget=1, ka_adjust_text={'arrowprops': {'arrowstyle': '-', 'color': 'r', 'lw': 0.5}})
	labels = {x.get_text() for x in ax.texts} - {''}
	assert labels == set(data.index)


########################################################################
# fig_diff_rank
########################################################################

@pytest.fixture
def diffdata():
	"""Differential expression and regulation logFCs of three genes."""
	return pd.DataFrame({'DE_logFC': [2.0, -1.0, 0.5], 'DR_logFC': [-3.0, 1.0, 0.0]}, index=['A', 'B', 'C'])


def test_fig_diff_rank_bars(diffdata, adjusttext):
	"""Three panels (mean, DE, DR) show values sorted descending, colored by value/max|value| on coolwarm."""
	data0 = diffdata.copy()
	fig = static.fig_diff_rank(diffdata)
	pd.testing.assert_frame_equal(diffdata, data0)
	assert len(fig.axes) == 3
	assert [x.get_ylabel() for x in fig.axes] == ['Mean logFC', 'DE logFC', 'DR logFC']
	expected = [[0.25, 0, -0.5], [2, 0.5, -1], [1, 0, -3]]
	cmap = plt.get_cmap('coolwarm')
	for ax, v in zip(fig.axes, expected):
		v = np.array(v, dtype=float)
		assert len(ax.patches) == 3
		assert_allclose([x.get_height() for x in ax.patches], v)
		assert_allclose([x.get_x() + x.get_width() / 2 for x in ax.patches], [0, 1, 2])
		assert_allclose([x.get_width() for x in ax.patches], [1, 1, 1])
		assert_allclose([x.get_facecolor() for x in ax.patches], cmap((v / np.abs(v).max() + 1) / 2))
		assert_allclose(ax.get_xlim(), [-0.53, 2.53])
		assert list(ax.get_xticks()) == []
		assert len(ax.texts) == 0
	assert_allclose(fig.get_size_inches(), [0.015 * 3, 2 * 3])


def test_fig_diff_rank_annotations(diffdata, adjusttext):
	"""Annotated genes are labelled at their bar tops with ka_text applied; ka_adjust_text=None keeps positions."""
	fig = static.fig_diff_rank(diffdata, annotate=['A', 'C'], ka_text={'fontsize': 5}, ka_adjust_text=None)
	positions = [{'A': (2, -0.5), 'C': (0, 0.25)}, {'A': (0, 2), 'C': (1, 0.5)}, {'A': (2, -3), 'C': (1, 0)}]
	for ax, pos in zip(fig.axes, positions):
		assert [x.get_text() for x in ax.texts] == ['A', 'C']
		for x in ax.texts:
			assert_allclose(x.get_position(), pos[x.get_text()])
			assert x.get_fontsize() == 5
			assert x.get_rotation() == 90


def test_fig_diff_rank_options(diffdata, adjusttext):
	"""figsize, cmap and bar keywords take effect; adjustText runs on annotations."""
	fig = static.fig_diff_rank(diffdata, figsize=(0.5, 1), annotate=['B'], cmap='viridis', ka_adjust_text={'arrowprops': {'arrowstyle': '-', 'color': 'r', 'lw': 0.5}}, width=0.8)
	assert_allclose(fig.get_size_inches(), [1.5, 3])
	cmap = plt.get_cmap('viridis')
	v = np.array([2, 0.5, -1])
	ax = fig.axes[1]
	assert_allclose([x.get_width() for x in ax.patches], [0.8] * 3)
	assert_allclose([x.get_facecolor() for x in ax.patches], cmap((v / 2 + 1) / 2))
	assert all([x.get_text() for x in a.texts] == ['B'] for a in fig.axes)


########################################################################
# fig_subnet
########################################################################

def _check_layout(nodes, ax):
	"""Nodes are centered with unit RMS distance and drawn at the returned coordinates."""
	assert list(nodes.columns) == ['x', 'y']
	assert np.isfinite(nodes.values).all()
	assert_allclose(nodes.values.mean(axis=0), [0, 0], atol=1e-10)
	assert_allclose(np.sqrt((nodes.values ** 2).mean()), 1, rtol=1e-10)
	assert_allclose(np.asarray(ax.collections[0].get_offsets()), nodes.values)
	assert not ax.axison


def test_fig_subnet_regulators(d0):
	"""Outgoing binarized edges of the chosen regulators in S2 are drawn as arrows with annotated nodes."""
	np.random.seed(0)
	fig = plt.figure()
	ax = fig.add_subplot(111)
	nodes, edges = static.fig_subnet(d0, ax, 'S2', regulators=['TF2', 'TF4'], sparsity=0.25, annotate=['TF2', 'TG1', 'missing'])
	assert list(nodes.index) == ['TF2', 'TF4'] + TARGETS
	_check_layout(nodes, ax)
	assert list(edges.columns) == ['Regulator', 'Target', 'Strength']
	assert list(edges['Regulator']) == ['TF2'] * 4 + ['TF4'] * 2
	assert list(edges['Target']) == TARGETS
	assert_allclose(edges['Strength'].astype(float).values, [10, 11, 12, 13, -20, -21])
	assert len(ax.patches) == 6
	assert all(x.get_linewidth() == pytest.approx(0.05) for x in ax.patches)
	assert_allclose(ax.collections[0].get_sizes(), [5])
	texts = {x.get_text(): x.get_position() for x in ax.texts if x.get_text() != ''}
	assert set(texts) == {'TF2', 'TG1'}
	for x in texts:
		assert_allclose(texts[x], nodes.loc[x].values)


def test_fig_subnet_targets_annotate_all(d0):
	"""targets restricts to incoming edges of those genes; annotate='all' labels every node shown."""
	np.random.seed(1)
	fig = plt.figure()
	ax = fig.add_subplot(111)
	nodes, edges = static.fig_subnet(d0, ax, 'S2', targets=['TG1', 'TG2'], sparsity=0.25, annotate='all')
	assert list(nodes.index) == ['TF2', 'TG1', 'TG2']
	_check_layout(nodes, ax)
	assert list(edges['Regulator']) == ['TF2', 'TF2']
	assert list(edges['Target']) == ['TG1', 'TG2']
	assert_allclose(edges['Strength'].astype(float).values, [10, 11])
	assert len(ax.patches) == 2
	assert {x.get_text() for x in ax.texts} == {'TF2', 'TG1', 'TG2'}


def test_fig_subnet_whole_network_options(d0):
	"""Without restrictions all binarized edges of S3 are shown; transformation, ka_node and ka_edge take effect."""
	np.random.seed(2)
	fig = plt.figure()
	ax = fig.add_subplot(111)
	nodes, edges = static.fig_subnet(d0, ax, 'S3', sparsity=0.25, transformation=lambda x: 2 * x, ka_node={'s': 30}, ka_edge={'lw': 0.5})
	assert list(nodes.index) == ['TF3', 'TF4'] + TARGETS
	_check_layout(nodes, ax)
	assert list(edges['Regulator']) == ['TF3'] * 3 + ['TF4'] * 3
	assert list(edges['Target']) == TARGETS
	assert_allclose(edges['Strength'].astype(float).values, 2 * np.array([10, 11, 12, 30, 31, 32]))
	assert_allclose(ax.collections[0].get_sizes(), [30])
	assert len(ax.patches) == 6
	assert all(x.get_linewidth() == pytest.approx(0.5) for x in ax.patches)
	assert [x.get_text() for x in ax.texts] == []
