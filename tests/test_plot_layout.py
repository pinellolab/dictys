"""
Tests for dictys.plot.layout: the abstract layout base and the notch layout,
which arranges dictys.plot.panel panels into an animation of a dynamic network.

The dynamic network has 6 cells, 3 states at the nodes of the path trajectory
0-1-2 (edges of length 1) and 4 genes, of which G1 and G2 are regulators.
Smoothed statistics are checked against a numpy Gaussian-kernel reference
over the three states.
"""

import os

import numpy as np
import pytest
from numpy.testing import assert_allclose

matplotlib = pytest.importorskip('matplotlib')
import matplotlib.animation  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import StrMethodFormatter  # noqa: E402
from matplotlib.colors import to_rgba_array  # noqa: E402

from dictys.traj import trajectory, point  # noqa: E402
from dictys.net import dynamic_network  # noqa: E402
from dictys.plot import layout, panel  # noqa: E402

CELL_W = np.array([[1, .5, 0, 0, 0, 0], [0, 0, 1, .25, 0, 0], [0, 0, 0, 0, .5, 1]])
COORD = np.array([[0., 1, 2, 3, 4, 5], [0., 1, 0, -1, 0, 1]])
CPM = np.array([[1., 3, 7], [15, 7, 3], [0, 1, 3], [2, 2, 2]])
# Edge weights (regulator, target, state) for regulators G1, G2 and targets G1-G4
W = np.array([
	[[0, .9, .5, .1], [.2, 0, -.7, .3]],
	[[0, .6, -.4, .8], [.5, 0, .1, -.2]],
	[[0, .1, .3, .2], [.9, 0, .6, -.8]],
]).transpose(1, 2, 0)
GENES = ['G1', 'G2', 'G3', 'G4']
DIST = 0.8
NFRAME = 5
S = np.linspace(0, 2, NFRAME)


def make_network():
	"""Dynamic network whose states sit at the nodes of the path 0-1-2, each the weighted centre of two cells."""
	traj = trajectory(np.array([[0, 1], [1, 2]]), np.array([1., 1.]))
	cells = point(traj, np.array([0, 0, 0, 0, 1, 1]), np.array([0., 0, 1, 1, 1, 1]))
	return dynamic_network(
		cname=np.array(['C{}'.format(x) for x in range(6)]), sname=np.array(['S0', 'S1', 'S2']),
		nname=np.array(GENES), nids=[np.array([0, 1]), np.array([0, 1, 2, 3])],
		traj=traj, point={'c': cells, 's': point.fromnodes(traj)},
		cprop={'coord': COORD.copy(), 'type': np.array(['T2', 'T1', 'T1', 'T3', 'T2', 'T2'])},
		scprop={'w': CELL_W.copy()}, ssprop={'traj-neighbor': np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])},
		nsprop={'cpm': CPM.copy()}, esprop={'w': W.copy()})


def draw_notch(nframe=NFRAME, dist=DIST, branch=(0, 2), **ka):
	"""notch.draw on the test network with two rows of panels unless overridden."""
	args = dict(sparsity=0.5, panelsize=(1, 1), bcde_tfs=[['G1'], ['G1', 'G2']], e_targets=['G3', 'G2'],
		f_tfs=[['G1'], ['G2']], f_stop=5, f_niteration=10)
	args.update(ka)
	return layout.notch(dist=dist, nframe=nframe, dpi=20).draw(make_network(), branch, **args)


def gauss(values, s, dist=DIST):
	"""Gaussian-kernel smoothing of values at the 3 states (last axis) onto path positions s."""
	d = np.abs(np.asarray(s, dtype=float)[None, :] - np.arange(3)[:, None])
	w = np.exp(-(d / dist) ** 2 / 2)
	return values @ (w / w.sum(axis=0))


def ref_lcpm(s, dist=DIST):
	return gauss(np.log2(CPM + 1), s, dist)


def ref_bin(s, sparsity):
	"""Binarized smoothed network: entries at least as strong as the int(sparsity*8)-th strongest per point."""
	a = np.abs(gauss(W, s))
	k = int(sparsity * a.shape[0] * a.shape[1])
	cut = np.sort(a.reshape(-1, a.shape[2]), axis=0)[-k]
	return a >= cut


def ref_lntarget(s, sparsity, weighted):
	"""Log2 (outdegree + 1) of the binarized network, or of the network rescaled to the sparsity when weighted."""
	if weighted:
		a = np.abs(gauss(W, s))
		a = a * (sparsity / a.mean(axis=(0, 1)))
		return np.log2(a.sum(axis=1) + 1)
	return np.log2(ref_bin(s, sparsity).sum(axis=1) + 1)


def expanded_lims(v):
	"""Finite min and max of v expanded by 2% of the range on each side."""
	lo, hi = np.nanmin(v), np.nanmax(v)
	return np.array([lo - 0.02 * (hi - lo), hi + 0.02 * (hi - lo)])


@pytest.fixture(autouse=True)
def close_figures():
	"""Close every figure after each test."""
	yield
	plt.close('all')


@pytest.fixture(autouse=True)
def seeded_global_rng():
	"""flayout_base draws initial node positions from numpy's global random state."""
	state = np.random.get_state()
	np.random.seed(0)
	yield
	np.random.set_state(state)


@pytest.fixture
def drawn():
	"""notch.draw output with its animation engine initialised (post_init applied)."""
	pts, fig, panels, animate_ka = draw_notch()
	ca = panel.animate_generic(pts, fig, panels)
	ca.init(**animate_ka)
	return pts, fig, panels, ca


# base

class minimal(layout.base):
	def draw(self, net, branch):
		return (None, None, [], {})


def test_base_is_abstract():
	"""layout.base cannot be instantiated without draw."""
	with pytest.raises(TypeError):
		layout.base()


def test_base_stores_settings():
	"""base keeps the smoothing distance, frame count and dpi, defaulting to 1.5, 200 and 200."""
	x = minimal()
	assert (x.dist, x.nframe, x.dpi) == (1.5, 200, 200)
	x = minimal(dist=0.5, nframe=7, dpi=30)
	assert (x.dist, x.nframe, x.dpi) == (0.5, 7, 30)
	assert minimal().draw(None, (0, 1)) == (None, None, [], {})


def test_notch_is_a_layout():
	"""notch subclasses layout.base and takes the same settings."""
	x = layout.notch(dist=0.3, nframe=9, dpi=11)
	assert isinstance(x, layout.base)
	assert (x.dist, x.nframe, x.dpi) == (0.3, 9, 11)
	assert (layout.notch().dist, layout.notch().nframe, layout.notch().dpi) == (1.5, 200, 200)


# notch.draw: structure

def test_notch_rejects_mismatched_rows():
	"""bcde_tfs and f_tfs must have one entry per row."""
	with pytest.raises(ValueError):
		draw_notch(bcde_tfs=[['G1'], ['G2']], f_tfs=[['G1']])


@pytest.mark.parametrize('branch', [(0, 2), (2, 0)])
def test_notch_points_follow_branch(branch):
	"""The returned points are nframe evenly spaced points from the branch start to its end."""
	pts = draw_notch(branch=branch)[0]
	assert len(pts) == NFRAME
	nodes = point.fromnodes(pts.p)
	assert_allclose((pts - nodes[[branch[0]]]).ravel(), S)
	assert_allclose((pts - nodes[[branch[1]]]).ravel(), S[::-1])


@pytest.mark.parametrize('panelsize,nrow', [((1, 1), 2), ((2, 1.5), 1)])
def test_notch_figure_grid(panelsize, nrow):
	"""The figure has 4 columns and 1+n rows of panels of the given size, at the layout's dpi."""
	pts, fig, panels, animate_ka = draw_notch(panelsize=panelsize, bcde_tfs=[['G1']] * nrow, f_tfs=[['G1']] * nrow)
	assert_allclose(fig.get_size_inches(), [panelsize[0] * 4, panelsize[1] * (1 + nrow)])
	assert fig.dpi == 20
	axes = fig.axes
	assert len(axes) == 4 * (1 + nrow)
	for xi, ax in enumerate(axes):
		spec = ax.get_subplotspec()
		assert (spec.rowspan.start, spec.colspan.start) == (xi // 4, xi % 4)
	for ax in axes[:3]:
		assert not any(x.get_visible() for x in ax.spines.values())
		assert len(ax.get_xticks()) == 0 and len(ax.get_yticks()) == 0
	for ax in axes[3:]:
		assert not ax.spines['top'].get_visible() and not ax.spines['right'].get_visible()
		assert ax.spines['left'].get_visible() and ax.spines['bottom'].get_visible()
	assert list(animate_ka) == ['post_init']


@pytest.mark.parametrize('adjust,expected', [
	('default', {'wspace': 0.35, 'hspace': 0.2}),
	({'wspace': 0.5, 'hspace': 0.1, 'left': 0.2}, {'wspace': 0.5, 'hspace': 0.1, 'left': 0.2}),
	(None, None),
])
def test_notch_subplots_adjust(adjust, expected):
	"""subplots_adjust is applied to the figure, and skipped when None."""
	ka = {} if adjust == 'default' else {'subplots_adjust': adjust}
	fig = draw_notch(**ka)[1]
	if expected is None:
		expected = {x: matplotlib.rcParams['figure.subplot.' + x] for x in ['wspace', 'hspace', 'left']}
	for key, value in expected.items():
		assert getattr(fig.subplotpars, key) == pytest.approx(value)


def test_notch_panel_types_and_axes():
	"""Panels are a (cells) and b (scatter) on the top row, then c, d, e, f on each row, each on its own axes."""
	pts, fig, panels, _ = draw_notch()
	types = [panel.cellscatter, panel.statscatter] + [panel.statplot, panel.statplot, panel.statheatmap, panel.network] * 2
	assert [type(x) for x in panels] == types
	axes = fig.axes
	assert [x.ax for x in panels] == [axes[0], axes[3]] + axes[4:]
	assert all(x.pts is pts for x in panels)


def test_notch_row_selection():
	"""Each row's c, d, e, f panels show that row's TFs; e shows e_targets; b annotates every TF of every row."""
	panels = draw_notch()[2]
	assert panels[1].names == ['G1', 'G2']
	assert sorted(panels[1].annotate.values()) == ['G1', 'G2']
	assert panels[2].panels[1].names == ['G1'] and panels[3].panels[1].names == ['G1']
	assert panels[6].panels[1].names == ['G1', 'G2'] and panels[7].panels[1].names == ['G1', 'G2']
	assert panels[4].names == [['G1'], ['G3', 'G2']]
	assert panels[8].names == [['G1', 'G2'], ['G3', 'G2']]
	assert list(panels[5].panels[1].annotate.values()) == ['G1']
	assert list(panels[9].panels[1].annotate.values()) == ['G2']


def test_notch_without_rows():
	"""With no rows only the cell and scatter panels are drawn on a single row."""
	pts, fig, panels, animate_ka = draw_notch(bcde_tfs=[], f_tfs=[])
	assert [type(x) for x in panels] == [panel.cellscatter, panel.statscatter]
	assert len(fig.axes) == 4
	assert panels[1].annotate == {}
	ca = panel.animate_generic(pts, fig, panels)
	assert len(ca.init(**animate_ka)) > 0


# notch.draw: drawn values

@pytest.mark.parametrize('dist,nframe', [(DIST, NFRAME), (0.3, 3)])
def test_notch_expression_panel(dist, nframe):
	"""Panel c plots Gaussian-smoothed log2 (CPM+1) of each row TF against pseudotime, with a pointer per frame."""
	pts, fig, panels, animate_ka = draw_notch(dist=dist, nframe=nframe)
	ca = panel.animate_generic(pts, fig, panels)
	ca.init(**animate_ka)
	s = np.linspace(0, 2, nframe)
	expected = ref_lcpm(s, dist)
	for ax, rows in [(fig.axes[4], [0]), (fig.axes[8], [0, 1])]:
		lines = ax.get_lines()
		assert len(lines) == len(rows)
		for line, row in zip(lines, rows):
			assert_allclose(line.get_xdata(), s, atol=1e-12)
			assert_allclose(line.get_ydata(), expected[row])
	ca.draw(nframe - 1)
	assert_allclose(fig.axes[8].collections[0].get_offsets(), [[2, expected[0, -1]], [2, expected[1, -1]]])
	assert_allclose(fig.axes[8].collections[0].get_facecolors(), to_rgba_array(matplotlib.colormaps['tab10'].colors[:2]))


@pytest.mark.parametrize('weighted,sparsity', [(False, 0.5), (False, 0.25), (True, 0.5)])
def test_notch_target_count_panel(weighted, sparsity):
	"""Panel d plots log2 (1 + outdegree) of the binarized network, or of the sparsity-rescaled weighted network."""
	pts, fig, panels, animate_ka = draw_notch(weighted=weighted, sparsity=sparsity)
	panel.animate_generic(pts, fig, panels).init(**animate_ka)
	expected = ref_lntarget(S, sparsity, weighted)
	assert_allclose(fig.axes[5].get_lines()[0].get_ydata(), expected[0])
	lines = fig.axes[9].get_lines()
	assert len(lines) == 2
	for line, row in zip(lines, expected):
		assert_allclose(line.get_xdata(), S, atol=1e-12)
		assert_allclose(line.get_ydata(), row)


@pytest.mark.parametrize('weighted', [False, True])
def test_notch_differential_panel(weighted):
	"""Panel b scatters each TF's change in log CPM against its change in log target count, on shared symmetric limits."""
	pts, fig, panels, animate_ka = draw_notch(weighted=weighted)
	ca = panel.animate_generic(pts, fig, panels)
	ca.init(**animate_ka)
	dx = ref_lcpm(S)[:2] - ref_lcpm(S[:1])[:2]
	lnt = ref_lntarget(S, 0.5, weighted)
	dy = lnt - lnt[:, :1]
	ax = fig.axes[3]
	assert ax.get_xlabel() == 'Differential expression{}logFC in CPM'.format(os.linesep)
	assert ax.get_ylabel() == 'Differential regulation{}logFC in target count'.format(os.linesep)
	assert ax.get_aspect() == 1
	lim = max(np.abs(expanded_lims(dx)).max(), np.abs(expanded_lims(dy)).max())
	assert_allclose(ax.get_xlim(), [-lim, lim])
	assert_allclose(ax.get_ylim(), [-lim, lim])
	for t in [2, 4]:
		ca.draw(t)
		assert_allclose(ax.collections[0].get_offsets(), np.array([dx[:, t], dy[:, t]]).T, atol=1e-12)


def test_notch_heatmap_panel(drawn):
	"""Panel e shows the smoothed edge weights from the row TFs to e_targets, with symmetric color limits."""
	pts, fig, panels, ca = drawn
	net = gauss(W, S)
	vmax = np.abs(expanded_lims(net)).max()
	for ax, rows in [(fig.axes[6], [0]), (fig.axes[10], [0, 1])]:
		im = ax.images[0]
		assert_allclose(im.get_clim(), [-vmax, vmax])
		assert im.get_cmap().name == 'coolwarm'
		assert_allclose(np.reshape(im.get_array(), (len(rows), 2)), net[rows][:, [2, 1], 0])
		assert_allclose(ax.get_xticks(), [0, 1])
	# post_init hides the target names above the last row
	assert [x.get_text() for x in fig.axes[10].get_xticklabels()] == ['G3', 'G2']
	assert [x.get_text() for x in fig.axes[6].get_xticklabels()] == ['', '']
	ca.draw(3)
	assert_allclose(np.squeeze(fig.axes[10].images[0].get_array()), net[:, [2, 1], 3])


def test_notch_cell_panel(drawn):
	"""Panel a draws every cell at its coordinates, with alpha from its max-normalised smoothed state weight."""
	pts, fig, panels, ca = drawn
	ax = fig.axes[0]
	cells = panels[0].panels[1].objs[0]
	pointer = panels[0].panels[0].objs[0]
	assert not ax.axison
	assert [x.get_text() for x in ax.texts] == ['T1', 'T2', 'T3']
	assert_allclose(cells.get_offsets(), COORD.T)
	w = gauss(CELL_W.T, S)
	w = w / w.max(axis=0)
	for t in [0, 3]:
		ca.draw(t)
		assert_allclose(cells.get_facecolors()[:, 3], 0.05 + 0.45 * w[:, t])
		assert_allclose(pointer.get_offsets(), [(COORD * w[:, t]).sum(axis=1) / w[:, t].sum()])


def test_notch_subnetwork_panel(drawn):
	"""Panel f draws the strongest edges of the row TF as arrows and shows only nodes with an edge at each frame."""
	pts, fig, panels, ca = drawn
	strong = ref_bin(S, 0.5) & (gauss(W, S) != 0)
	ax = fig.axes[7]
	assert not ax.axison
	node = panels[5].panels[1]
	keep = strong[0].any(axis=1)
	keep[0] |= strong[0].any()
	assert node.names == [x for x, y in zip(GENES, keep) if y]
	assert node.scatterka == {'s': 5, 'lw': 0}
	assert panels[5].panels[0].ka == {'lw': 0.05}
	for t in [0, 4, 2]:
		ca.draw(t)
		assert sum(x.get_alpha() != 0 for x in ax.patches) == strong[0, :, t].sum()
		shown = strong[0, :, t].copy()
		shown[0] |= shown.any()
		shown = shown[keep]
		assert (np.isfinite(ax.collections[0].get_offsets()).all(axis=1) == shown).all()
		assert all(x.get_linewidth() == pytest.approx(0.05) for x in ax.patches)


def test_notch_post_init_labels(drawn):
	"""post_init clears x labels above the last row and sets the c and d y labels on every row."""
	axes = drawn[1].axes
	for ax in axes[4:7]:
		assert ax.get_xlabel() == ''
		fmt = ax.xaxis.get_major_formatter()
		assert isinstance(fmt, StrMethodFormatter) and fmt.fmt == ''
	assert axes[8].get_xlabel() == 'Pseudo-time'
	assert axes[10].get_xlabel() == 'Target'
	assert not isinstance(axes[8].xaxis.get_major_formatter(), StrMethodFormatter)
	assert [axes[x].get_ylabel() for x in [4, 8]] == ['LogCPM', 'LogCPM']
	assert [axes[x].get_ylabel() for x in [5, 9]] == ['Log target count'] * 2


# notch.draw: keyword arguments

def test_notch_panel_keywords():
	"""a_ka to f_ka update each panel's defaults without changing the caller's dicts."""
	f_ka = {'edgeka': {'lw': 0.5}, 'nodeka': {'scatterka': {'s': 3}}}
	panels = draw_notch(
		a_ka={'pointer': False, 'scatterka': {'legend_loc': (0.6, 1)}},
		b_ka={'lim': {'sym'}, 'aspect': 2},
		c_ka={'pointer': False},
		d_ka={'cmap': 'Set1'},
		e_ka={'lim': [-0.02, 0.02], 'cmap': 'bwr'},
		f_ka=f_ka)[2]
	assert f_ka == {'edgeka': {'lw': 0.5}, 'nodeka': {'scatterka': {'s': 3}}}
	assert [type(x) for x in panels[0].panels] == [panel.cellscatter_scatter]
	assert panels[0].panels[0].legend_loc == (0.6, 1)
	assert panels[1].lim == {'sym'} and panels[1].aspect == 2
	assert [type(x) for x in panels[2].panels] == [panel.statplot_static]
	assert_allclose(panels[3].panels[1].colors, to_rgba_array(matplotlib.colormaps['Set1'].colors[:1]))
	assert panels[4].lim == [-0.02, 0.02] and panels[4].ka == {'cmap': 'bwr', 'aspect': 'auto'}
	assert panels[5].panels[0].ka == {'lw': 0.5}
	assert panels[5].panels[1].scatterka == {'s': 3}
	assert list(panels[9].panels[1].annotate.values()) == ['G2']


def test_notch_default_panel_keywords():
	"""Without keywords, b uses aspect 1 with 'sym','min','max' limits and e uses the coolwarm colormap."""
	panels = draw_notch()[2]
	assert panels[1].aspect == 1 and panels[1].lim == {'sym', 'min', 'max'}
	assert panels[4].ka == {'cmap': 'coolwarm', 'aspect': 'auto'}
	assert [type(x) for x in panels[0].panels] == [panel.cellscatter_pointer, panel.cellscatter_scatter]


def test_notch_subnetwork_layout_options():
	"""f_stop and f_niteration reach the force-directed layout; f_weight_func sets how layouts are smoothed over points."""
	pts, fig, panels, _ = draw_notch(f_stop=3, f_niteration=7)
	smooth = panels[5].panels[1].stats[0].stats[0]
	assert smooth.stat.func.keywords == {'stop': 3, 'iterations': 7}
	# The default linear interpolation at the layout's own points returns the layouts unchanged.
	assert_allclose(smooth.compute(pts), smooth.stat.compute(pts))
	pts, fig, panels, _ = draw_notch(f_weight_func=('conv', [100.], {}))
	smooth = panels[5].panels[1].stats[0].stats[0]
	raw = smooth.stat.compute(pts)
	expected = np.repeat(np.nanmean(raw, axis=2, keepdims=True), len(pts), axis=2)
	assert_allclose(smooth.compute(pts), expected, atol=1e-3)


def test_notch_animation_end_to_end(tmp_path):
	"""The panels animate through every frame: saving a GIF leaves the last frame and post_init labels drawn."""
	pytest.importorskip('PIL')
	pts, fig, panels, animate_ka = draw_notch()
	anim = panel.animate_generic(pts, fig, panels).animate(**animate_ka)
	anim.save(str(tmp_path / 'notch.gif'), writer=matplotlib.animation.PillowWriter(fps=5))
	assert os.path.getsize(str(tmp_path / 'notch.gif')) > 0
	assert fig.axes[4].get_ylabel() == 'LogCPM'
	assert_allclose(np.squeeze(fig.axes[10].images[0].get_array()), gauss(W, S)[:, [2, 1], -1])
