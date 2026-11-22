"""
Tests for dictys.plot.dynamic: auc, draw_discover1, fig_discover and
fig_regulation_heatmap (on a hand-built dictys.net.dynamic_network).
"""

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose

matplotlib = pytest.importorskip('matplotlib')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402

from dictys.plot import dynamic  # noqa: E402


@pytest.fixture(autouse=True)
def close_figures():
	"""Close every matplotlib figure a test opened."""
	yield
	plt.close('all')


########################################################################
# auc
########################################################################

def test_auc_hand_values():
	"""Trapezoid areas on uneven x match hand-computed values."""
	dx = np.array([0., 1., 3.])
	dy = np.array([[1., 1., 1.], [0., 2., 2.], [1., 3., 0.]])
	assert_allclose(dynamic.auc(dx, dy), [3, 5, 5])


def test_auc_linear_function_is_exact():
	"""The area under 2x+1 on [0,2] is 6 for any increasing sampling including the ends."""
	dx = np.array([0., 0.3, 1.1, 2.])
	assert_allclose(dynamic.auc(dx, (2 * dx + 1).reshape(1, -1)), [6])


def test_auc_matches_scipy(rng):
	"""Random curves agree with scipy.integrate.trapezoid."""
	from scipy.integrate import trapezoid
	dx = np.cumsum(rng.random(7) + 0.1)
	dy = rng.normal(size=(4, 7))
	assert_allclose(dynamic.auc(dx, dy), trapezoid(dy, dx, axis=1), rtol=1e-12)


@pytest.mark.parametrize('dx', [np.array([0.]), np.array([0., 1., 1.]), np.array([0., 2., 1.])])
def test_auc_invalid_x_raises(dx):
	"""Fewer than two x values or non-increasing x raise ValueError."""
	with pytest.raises(ValueError):
		dynamic.auc(dx, np.ones((2, len(dx))))


########################################################################
# draw_discover1
########################################################################

@pytest.fixture
def curves():
	"""Curves of 4 TFs over 5 time points and two characteristics per TF."""
	dcurve = pd.DataFrame([
		[0, 1, 2, 3, 4],
		[4, 3, 2, 1, 0],
		[0, 2, 4, 2, 0],
		[1, 1, 1, 1, 2],
	], index=['A', 'B', 'C', 'D'], columns=['t0', 't1', 't2', 't3', 't4'], dtype=float)
	dchar = pd.DataFrame({'c1': [1.0, 2.0, 3.0, 5.0], 'c2': [-1.0, 0.0, 0.5, 1.0]}, index=dcurve.index)
	return dcurve, dchar


def _figsize(n, ns, inset_sides=2, inset_lvs=1, inset_size=(0.55, 0.3), inset_space0=0.15, inset_space=(0.05, 0), heatmap_height=0.1, heatmap_space=0.05):
	"""
	Figure size from the documented component sizes in inches. Height: the
	stacked heatmaps, nrow inset rows and the gaps between them. Width: the n
	insets are staggered across nrow rows with a pitch of
	(inset_size[0]+inset_space[0])/nrow, so they span pitch*(n-1)+inset_size[0],
	plus inset_space[0]/2 of margin on each side.
	"""
	nside = (inset_sides & 1 != 0) + (inset_sides & 2 != 0)
	nrow = inset_lvs * nside
	pitch = (inset_size[0] + inset_space[0]) / nrow
	return (
		pitch * (n - 1) + inset_size[0] + inset_space[0],
		heatmap_height * ns + heatmap_space * (ns - 1) + inset_size[1] * nrow + (inset_space0 + inset_space[1] * (inset_lvs - 1)) * nside,
	)


def _layout_boxes(fig, axes):
	"""
	(left, bottom, width, height) of each heatmap and curve inset, in inches of
	the layout: fractions of the main axes times the figure size, the units
	draw_discover1 places its insets in.
	"""
	fig.canvas.draw()
	size = fig.get_size_inches()
	p0 = axes[0].get_position()
	ans = []
	for ax in axes[1:]:
		p = ax.get_position()
		ans.append([(p.x0 - p0.x0) / p0.width * size[0], (p.y0 - p0.y0) / p0.height * size[1], p.width / p0.width * size[0], p.height / p0.height * size[1]])
	return np.array(ans)


def test_draw_discover1_layout(curves):
	"""One heatmap row per characteristic (last first) and one curve inset per TF, normalized to the global range."""
	dcurve, dchar = curves
	fig, axes, cmap = dynamic.draw_discover1(dcurve, dchar)
	assert len(axes) == 1 + 2 + 4
	assert not axes[0].axison
	assert set(cmap) == {'c1', 'c2'}
	vir = plt.get_cmap('viridis')
	for ax, col in zip(axes[1:3], ['c2', 'c1']):
		assert [x.get_text() for x in ax.get_yticklabels()] == [col]
		assert list(ax.get_xticks()) == []
		v = dchar[col].values
		assert (cmap[col].norm.vmin, cmap[col].norm.vmax) == (v.min(), v.max())
		img = np.asarray(ax.images[0].get_array())
		assert img.shape == (1, 4, 4)
		assert_allclose(img[0], vir((v - v.min()) / (v.max() - v.min())))
		assert ax.get_aspect() == 'auto'
	for ax, name in zip(axes[3:], dcurve.index):
		assert ax.get_xlabel() == name
		assert ax.xaxis.get_label_position() == 'top'
		x, y = ax.lines[0].get_data()
		assert_allclose(x, np.linspace(0, 1, 5))
		assert_allclose(y, dcurve.loc[name].values / 4)
		assert_allclose(ax.get_xlim(), [0, 1])
		assert_allclose(ax.get_ylim(), [0, 1.1])
		assert ax.lines[0].get_color() == 'k'
	# One indicator line per curve on the main axes
	assert len(axes[0].lines) == 4
	assert_allclose(fig.get_size_inches(), _figsize(4, 2))


def test_draw_discover1_component_geometry(curves):
	"""Heatmaps and insets have the documented sizes and gaps, and together span exactly the figure size."""
	dcurve, dchar = curves
	fig, axes, _ = dynamic.draw_discover1(dcurve, dchar)
	boxes = _layout_boxes(fig, axes)
	heat, insets = boxes[:2], boxes[2:]
	# Heatmap rows (c2 below c1): 0.1 high, 0.05 apart, as wide as the 4 inset pitches of 0.55+0.05
	assert_allclose(heat[:, 2:], [[2.4, 0.1]] * 2, atol=1e-9)
	assert_allclose(heat[1, 1] - (heat[0, 1] + 0.1), 0.05, atol=1e-9)
	# Insets: 0.55 by 0.3, one pitch apart, all 0.15 above the top heatmap
	assert_allclose(insets[:, 2:], [[0.55, 0.3]] * 4, atol=1e-9)
	assert_allclose(np.diff(insets[:, 0]), [0.6] * 3, atol=1e-9)
	assert_allclose(insets[:, 1], heat[1, 1] + 0.1 + 0.15, atol=1e-9)
	# The heatmaps span the full width; the insets leave inset_space[0]/2 on each side
	size = fig.get_size_inches()
	assert_allclose([insets[0, 0] - heat[0, 0], heat[0, 0] + heat[0, 2] - insets[-1, 0] - 0.55], [0.025, 0.025], atol=1e-9)
	assert_allclose(heat[0, 2], size[0], atol=1e-9)
	assert_allclose(insets[0, 1] + 0.3 - heat[0, 1], size[1], atol=1e-9)


def test_draw_discover1_dtime_sorts_columns(curves):
	"""Columns of dcurve are reordered by dtime, which is matched by index rather than position."""
	dcurve, dchar = curves
	dtime = pd.Series([3.0, 0.0, 1.0, 10.0, 2.0], index=['t4', 't0', 't1', 't2', 't3'])
	_, axes, _ = dynamic.draw_discover1(dcurve, dchar, dtime=dtime)
	x, y = axes[3].lines[0].get_data()
	# Times of t0..t4 are 0,1,10,2,3 so the order is t0,t1,t3,t4,t2
	assert_allclose(x, [0, 1, 2, 3, 10])
	assert_allclose(y, np.array([0, 1, 3, 4, 2]) / 4)
	assert_allclose(axes[3].get_xlim(), [0, 10])


def test_draw_discover1_matches_rows_by_index(curves):
	"""Only TFs present in both tables are drawn."""
	dcurve, dchar = curves
	_, axes, _ = dynamic.draw_discover1(dcurve, dchar.iloc[:2])
	assert [x.get_xlabel() for x in axes[3:]] == ['A', 'B']
	_, axes, _ = dynamic.draw_discover1(dcurve.iloc[1:3], dchar)
	assert [x.get_xlabel() for x in axes[3:]] == ['B', 'C']


def test_draw_discover1_cmap_dict(curves):
	"""A dict of colormaps keeps given ScalarMappables, normalizes names to the data range and drops unknown keys."""
	dcurve, dchar = curves
	sm = matplotlib.cm.ScalarMappable(norm=Normalize(-5, 5), cmap='coolwarm')
	_, axes, cmap = dynamic.draw_discover1(dcurve, dchar, cmap={'c1': sm, 'c2': 'magma', 'other': 'viridis'})
	assert set(cmap) == {'c1', 'c2'}
	assert cmap['c1'] is sm
	assert cmap['c2'].get_cmap().name == 'magma'
	assert (cmap['c2'].norm.vmin, cmap['c2'].norm.vmax) == (-1, 1)
	assert_allclose(np.asarray(axes[2].images[0].get_array())[0], plt.get_cmap('coolwarm')((dchar['c1'].values + 5) / 10))
	_, _, cmap = dynamic.draw_discover1(dcurve, dchar, cmap=sm)
	assert cmap['c1'] is sm and cmap['c2'] is sm


@pytest.mark.parametrize('inset_sides,inset_lvs,positions', [
	(3, 1, ['bottom', 'top', 'bottom', 'top']),
	(2, 2, ['top', 'top', 'top', 'top']),
	(3, 2, ['bottom', 'top', 'bottom', 'top']),
])
def test_draw_discover1_inset_arrangements(curves, inset_sides, inset_lvs, positions):
	"""inset_sides=3 alternates insets below and above; inset_lvs adds rows; figure size follows the documented parts."""
	dcurve, dchar = curves
	fig, axes, _ = dynamic.draw_discover1(dcurve, dchar, inset_sides=inset_sides, inset_lvs=inset_lvs, inset_space=(0.1, 0.2))
	assert len(axes) == 7
	assert [x.xaxis.get_label_position() for x in axes[3:]] == positions
	assert_allclose(fig.get_size_inches(), _figsize(4, 2, inset_sides=inset_sides, inset_lvs=inset_lvs, inset_space=(0.1, 0.2)))


def test_draw_discover1_options(curves):
	"""Sizes, spacing, font size, curve expansion and drawing keywords take effect."""
	dcurve, dchar = curves
	ka = dict(inset_size=(0.8, 0.5), inset_space0=0.3, inset_space=(0.1, 0), heatmap_height=0.2, heatmap_space=0.1)
	fig, axes, _ = dynamic.draw_discover1(dcurve, dchar, curve_expand=(0.2, 0.3), line_space=0.1, fs=5, ka_curve={'c': 'r', 'lw': 2}, ka_heatmap={'interpolation': 'nearest'}, **ka)
	assert_allclose(fig.get_size_inches(), _figsize(4, 2, **ka))
	# Indicator lines span the inset-heatmap gap minus line_space at each end, in the layout's
	# inches (main-axes fraction times figure size; the base keeps subplot margins, a bug)
	height = fig.get_size_inches()[1]
	assert len(axes[0].lines) == 4
	for line in axes[0].lines:
		assert_allclose(np.abs(np.diff(line.get_ydata())) * height, [0.3 - 2 * 0.1])
	for ax in axes[1:3]:
		assert ax.images[0].get_interpolation() == 'nearest'
		assert ax.get_yticklabels()[0].get_fontsize() == 5
	for ax in axes[3:]:
		assert_allclose(ax.get_ylim(), [-0.2, 1.3])
		assert ax.lines[0].get_color() == 'r'
		assert ax.lines[0].get_linewidth() == 2
		assert ax.xaxis.label.get_fontsize() == 5


def test_draw_discover1_invalid_input_raises(curves):
	"""Mismatched tables, too few time points, bad dtime, cmap type or layout values raise."""
	dcurve, dchar = curves
	with pytest.raises(ValueError):
		dynamic.draw_discover1(dcurve, dchar.rename(index=lambda x: x + 'x'))
	with pytest.raises(ValueError):
		dynamic.draw_discover1(dcurve.iloc[:, :1], dchar)
	with pytest.raises(ValueError):
		dynamic.draw_discover1(dcurve, dchar.iloc[::-1])
	with pytest.raises(TypeError):
		dynamic.draw_discover1(dcurve, dchar, dtime=pd.Series([0., 1.], index=['t0', 't1']))
	with pytest.raises(ValueError):
		dynamic.draw_discover1(dcurve, dchar, dtime=pd.Series(np.arange(5.), index=['t0', 't1', 't2', 't3', 'x']))
	with pytest.raises(TypeError):
		dynamic.draw_discover1(dcurve, dchar, cmap=['viridis'])
	with pytest.raises(AssertionError):
		dynamic.draw_discover1(dcurve.iloc[:3], dchar, inset_sides=3)
	with pytest.raises(AssertionError):
		dynamic.draw_discover1(dcurve, dchar, inset_lvs=3)
	with pytest.raises(AssertionError):
		dynamic.draw_discover1(dcurve, dchar, line_space=0.1)


########################################################################
# fig_discover
########################################################################

@pytest.fixture
def discover():
	"""
	Curves at 5 evenly spaced times with hand-derived characteristics
	(terminal logFC, transient logFC, switching time):
	act_early 4,0,0.125; act_lin 2,0,0.5; act_late 3,0,0.875; inact_lin -4,0,0.5;
	inact_early -3,0,0.125; bump 0,2,0; dip 0,-1,0; flat 0,0,0.
	"""
	dcurve = pd.DataFrame({
		'act_early': [0, 4, 4, 4, 4],
		'act_lin': [0, 0.5, 1, 1.5, 2],
		'act_late': [0, 0, 0, 0, 3],
		'inact_lin': [4, 3, 2, 1, 0],
		'inact_early': [3, 0, 0, 0, 0],
		'bump': [0, 2, 4, 2, 0],
		'dip': [1, 0, -1, 0, 1],
		'flat': [1, 1, 1, 1, 1],
	}, index=['t0', 't1', 't2', 't3', 't4']).T.astype(float)
	dtime = pd.Series(np.linspace(0, 1, 5), index=dcurve.columns)
	return dcurve, dtime


def _heat(axes):
	"""Characteristic name and color array of each heatmap row of one discovery figure."""
	return {ax.get_yticklabels()[0].get_text(): np.asarray(ax.images[0].get_array())[0] for ax in axes[1:3]}


def test_fig_discover_patterns(discover):
	"""Activating, inactivating, transient up and down TFs are selected and ordered from the hand-derived characteristics."""
	dcurve, dtime = discover
	ans = dynamic.fig_discover(dcurve, dtime, ntops=(3, 2, 1, 1))
	assert len(ans) == 4
	names = [[x.get_xlabel() for x in y[1][3:]] for y in ans]
	assert names == [['act_early', 'act_lin', 'act_late'], ['inact_early', 'inact_lin'], ['bump'], ['dip']]
	for fig, axes, cmap in ans[:2]:
		assert set(cmap) == {'Switching time', 'Terminal logFC'}
	for fig, axes, cmap in ans[2:]:
		assert set(cmap) == {'Transient logFC', 'Terminal logFC'}
	# Default ranges: terminal symmetric over all shown, switching extreme, transient symmetric
	norms = {}
	for _, _, cmap in ans:
		for x, y in cmap.items():
			norms.setdefault(x, set()).add((y.norm.vmin, y.norm.vmax, y.get_cmap().name))
	assert norms == {
		'Terminal logFC': {(-4, 4, 'coolwarm')},
		'Switching time': {(0.125, 0.875, 'viridis')},
		'Transient logFC': {(-2, 2, 'PiYG_r')},
	}
	heat = _heat(ans[0][1])
	assert_allclose(heat['Switching time'], plt.get_cmap('viridis')([0, 0.5, 1]), atol=1e-12)
	assert_allclose(heat['Terminal logFC'], plt.get_cmap('coolwarm')((np.array([4, 2, 3]) + 4) / 8), atol=1e-12)
	heat = _heat(ans[3][1])
	assert_allclose(heat['Transient logFC'], plt.get_cmap('PiYG_r')([(-1 + 2) / 4]), atol=1e-12)
	# Each figure normalizes its own curves to [0,1]
	x, y = ans[2][1][3].lines[0].get_data()
	assert_allclose(x, dtime.values)
	assert_allclose(y, [0, 0.5, 1, 0.5, 0])


def test_fig_discover_uneven_time(discover):
	"""Characteristics use time rescaled to [0,1]: switching times on times 0,1,2,4,8 are 0.0625 and 0.75."""
	dcurve, _ = discover
	dtime = pd.Series([0., 1., 2., 4., 8.], index=dcurve.columns)
	ans = dynamic.fig_discover(dcurve, dtime, ntops=(2, 0, 0, 0))
	assert len(ans) == 1
	_, axes, cmap = ans[0]
	assert [x.get_xlabel() for x in axes[3:]] == ['act_early', 'act_late']
	assert_allclose([cmap['Switching time'].norm.vmin, cmap['Switching time'].norm.vmax], [0.0625, 0.75])
	assert_allclose(axes[3].lines[0].get_xdata(), dtime.values)


def test_fig_discover_vrange_and_cmap(discover):
	"""Manual, symmetric and extreme ranges and per-characteristic colormaps are applied to the TFs shown."""
	dcurve, dtime = discover
	vrange = {'Terminal logFC': (-10, 10), 'Switching time': 'symmetric', 'Transient logFC': 'extreme'}
	cmap = {'Terminal logFC': 'RdBu', 'Switching time': 'magma', 'Transient logFC': 'BrBG'}
	ans = dynamic.fig_discover(dcurve, dtime, ntops=(2, 0, 1, 1), vrange=vrange, cmap=cmap)
	assert len(ans) == 3
	assert [[x.get_xlabel() for x in y[1][3:]] for y in ans] == [['act_early', 'act_late'], ['bump'], ['dip']]
	got = {}
	for _, _, c in ans:
		for x, y in c.items():
			got.setdefault(x, set()).add((y.norm.vmin, y.norm.vmax, y.get_cmap().name))
	assert got == {
		'Terminal logFC': {(-10, 10, 'RdBu')},
		'Switching time': {(-0.875, 0.875, 'magma')},
		'Transient logFC': {(-1, 2, 'BrBG')},
	}


def test_fig_discover_single_vrange_and_keywords(discover):
	"""A non-dict vrange applies to every characteristic; extra keywords reach draw_discover1."""
	dcurve, dtime = discover
	ans = dynamic.fig_discover(dcurve, dtime, ntops=(2, 2, 0, 0), vrange='extreme', fs=6, inset_sides=3)
	assert len(ans) == 2
	terminal = {(c['Terminal logFC'].norm.vmin, c['Terminal logFC'].norm.vmax) for _, _, c in ans}
	assert terminal == {(-4, 4)}
	switching = {(c['Switching time'].norm.vmin, c['Switching time'].norm.vmax) for _, _, c in ans}
	assert switching == {(0.125, 0.875)}
	for _, axes, _ in ans:
		assert all(x.xaxis.label.get_fontsize() == 6 for x in axes[3:])
		assert [x.xaxis.get_label_position() for x in axes[3:]] == ['bottom', 'top']


def test_fig_discover_no_patterns(discover):
	"""ntops of all zeros draws nothing."""
	dcurve, dtime = discover
	assert dynamic.fig_discover(dcurve, dtime, ntops=(0, 0, 0, 0)) == []


def test_fig_discover_unknown_vrange_raises(discover):
	"""An unknown automatic range name raises ValueError."""
	dcurve, dtime = discover
	with pytest.raises(ValueError):
		dynamic.fig_discover(dcurve, dtime, ntops=(1, 0, 0, 0), vrange={'Terminal logFC': 'widest'})


########################################################################
# fig_regulation_heatmap
########################################################################

# Network edge weights (regulator TF1/TF2, target TG1..TG3) at states s0, s1, s2
EDGE_W = np.array([
	[[0.5, -1.0, 2.0], [1.5, 0.2, -0.3], [-2.0, 0.0, 1.0]],
	[[-0.4, 0.8, 1.2], [0.0, -1.5, 0.6], [0.9, 0.9, -0.9]],
])


@pytest.fixture
def dnet():
	"""
	Dynamic network on a one-edge trajectory of length 2. States s0, s1, s2 and
	their two cells each sit at positions 0, 1 and 2, and the states form the
	chain s0-s1-s2, so the smoothed path runs over state positions 0..2.
	"""
	from dictys.net import dynamic_network
	from dictys.traj import trajectory, point
	traj = trajectory(np.array([[0, 1]]), np.array([2.0]))
	pc = point(traj, np.zeros(6, dtype=int), np.array([0., 0., 1., 1., 2., 2.]))
	ps = point(traj, np.zeros(3, dtype=int), np.array([0., 1., 2.]))
	membership = np.zeros((3, 6))
	for xi in range(3):
		membership[xi, 2 * xi:2 * xi + 2] = 1
	adj = np.zeros((3, 3))
	adj[0, 1] = adj[1, 0] = adj[1, 2] = adj[2, 1] = 1
	return dynamic_network(
		cname=np.array(['c%d' % x for x in range(6)]),
		sname=np.array(['s0', 's1', 's2']),
		nname=np.array(['TF1', 'TF2', 'TG1', 'TG2', 'TG3']),
		nids=[np.array([0, 1]), np.array([2, 3, 4])],
		traj=traj,
		point={'c': pc, 's': ps},
		scprop={'w': membership},
		ssprop={'traj-neighbor': adj},
		esprop={'w': EDGE_W},
		nsprop={'cpm': np.ones((5, 3))})


def _smoothed(reg_targets, num, dist):
	"""Reference: Gaussian kernel (sigma=dist) average over states at 0,1,2 of each weight, at num evenly spaced positions on [0,2]."""
	x = np.linspace(0, 2, num)
	k = np.exp(-((np.arange(3)[:, None] - x[None, :]) / dist) ** 2 / 2)
	k /= k.sum(axis=0)
	return np.array([EDGE_W[r, t] @ k for r, t in reg_targets])


@pytest.mark.parametrize('num,dist', [(5, 0.5), (7, 2.0)])
def test_fig_regulation_heatmap_values(dnet, num, dist):
	"""Each row is the kernel-smoothed regulation strength along the branch, colored on a range of +-95% quantile of |strength|."""
	regulations = [('TF1', 'TG2'), ('TF2', 'TG1')]
	fig, ax, cmap = dynamic.fig_regulation_heatmap(dnet, 0, 1, regulations, num=num, dist=dist)
	expected = _smoothed([(0, 1), (1, 0)], num, dist)
	q = np.quantile(np.abs(expected).ravel(), 0.95)
	assert_allclose([cmap.norm.vmin, cmap.norm.vmax], [-q, q], rtol=1e-10)
	assert cmap.get_cmap().name == 'coolwarm'
	img = np.asarray(ax.images[0].get_array())
	assert img.shape == (2, num, 4)
	assert_allclose(img, plt.get_cmap('coolwarm')((expected + q) / (2 * q)), atol=1e-6)
	assert [x.get_text() for x in ax.get_yticklabels()] == ['TF1-TG2', 'TF2-TG1']
	assert list(ax.get_xticks()) == []
	assert_allclose(fig.get_size_inches(), [2, 0.22 * 2])
	assert ax.get_aspect() == pytest.approx((0.22 * 2 / 2) / (2 / num))
	assert fig.axes == [ax]


def test_fig_regulation_heatmap_reverse_branch(dnet):
	"""Going from node 1 to node 0 reverses the columns."""
	regulations = [('TF1', 'TG3'), ('TF2', 'TG2'), ('TF1', 'TG1')]
	sm = matplotlib.cm.ScalarMappable(norm=Normalize(-3, 3), cmap='viridis')
	_, ax1, _ = dynamic.fig_regulation_heatmap(dnet, 0, 1, regulations, num=6, cmap=sm)
	_, ax2, _ = dynamic.fig_regulation_heatmap(dnet, 1, 0, regulations, num=6, cmap=sm)
	a1 = np.asarray(ax1.images[0].get_array())
	a2 = np.asarray(ax2.images[0].get_array())
	assert_allclose(a2, a1[:, ::-1], atol=1e-6)


def test_fig_regulation_heatmap_colormaps(dnet):
	"""A ScalarMappable is used as is; a name with vmax gives the range [-vmax,vmax]."""
	regulations = [('TF2', 'TG3')]
	expected = _smoothed([(1, 2)], 4, 1.5)
	sm = matplotlib.cm.ScalarMappable(norm=Normalize(-2, 1), cmap='viridis')
	_, ax, cmap = dynamic.fig_regulation_heatmap(dnet, 0, 1, regulations, num=4, cmap=sm)
	assert cmap is sm
	assert_allclose(np.asarray(ax.images[0].get_array()), plt.get_cmap('viridis')((expected + 2) / 3), atol=1e-6)
	fig, ax, cmap = dynamic.fig_regulation_heatmap(dnet, 0, 1, regulations, num=4, cmap='PiYG', vmax=4, figsize=(3, 0.5))
	assert_allclose(fig.get_size_inches(), [3, 0.5])
	assert ax.get_aspect() == pytest.approx(0.5 / (3 / 4))
	assert (cmap.norm.vmin, cmap.norm.vmax) == (-4, 4)
	assert cmap.get_cmap().name == 'PiYG'
	assert_allclose(np.asarray(ax.images[0].get_array()), plt.get_cmap('PiYG')((expected + 4) / 8), atol=1e-6)


def test_fig_regulation_heatmap_given_axes(dnet):
	"""On given axes (figsize=None) the cell aspect follows that figure's size."""
	fig0 = plt.figure(figsize=(4, 1))
	ax0 = fig0.add_subplot(111)
	fig, ax, _ = dynamic.fig_regulation_heatmap(dnet, 0, 1, [('TF1', 'TG1'), ('TF2', 'TG2')], num=8, ax=ax0, figsize=None)
	assert fig is fig0 and ax is ax0
	assert ax.get_aspect() == pytest.approx((1 / 2) / (4 / 8))
	assert np.asarray(ax.images[0].get_array()).shape == (2, 8, 4)


def test_fig_regulation_heatmap_invalid_input_raises(dnet):
	"""Unknown regulators or targets, ax with figsize, or vmax with a ScalarMappable raise ValueError."""
	with pytest.raises(ValueError):
		dynamic.fig_regulation_heatmap(dnet, 0, 1, [('TG1', 'TG2')], num=4)
	with pytest.raises(ValueError):
		dynamic.fig_regulation_heatmap(dnet, 0, 1, [('TF1', 'XX')], num=4)
	fig = plt.figure()
	ax = fig.add_subplot(111)
	with pytest.raises(ValueError):
		dynamic.fig_regulation_heatmap(dnet, 0, 1, [('TF1', 'TG1')], num=4, ax=ax, figsize=(2, 0.2))
	sm = matplotlib.cm.ScalarMappable(norm=Normalize(-2, 1), cmap='viridis')
	with pytest.raises(ValueError):
		dynamic.fig_regulation_heatmap(dnet, 0, 1, [('TF1', 'TG1')], num=4, cmap=sm, vmax=1)
