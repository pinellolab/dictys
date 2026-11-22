"""
Tests for dictys.plot.panel: the panels that draw statistics of a dynamic
network over pseudotime on matplotlib axes, and the animate_generic engine.

Most tests use a path trajectory 0-1-2 with two edges of length 1 and five
points at path positions 0, 0.5, 1, 1.5 and 2. Hand-made stats compute their
values from each point's path position, so expected values follow directly
from the formulas below.
"""

import functools

import numpy as np
import pytest
from numpy.testing import assert_allclose

matplotlib = pytest.importorskip('matplotlib')
import matplotlib.animation  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import to_rgba, to_rgba_array  # noqa: E402

from dictys.traj import trajectory, point  # noqa: E402
from dictys.net import stat, dynamic_network  # noqa: E402
from dictys.plot import panel  # noqa: E402

S = np.linspace(0, 2, 5)
NAMES = ['a', 'b', 'c']
TAB10 = matplotlib.colormaps['tab10'].colors


def fx(s):
	return np.array([s, 2 * s, -s])


def fy(s):
	return np.array([s ** 2, 1 + s, np.where(s > 1, np.nan, s)])


def fy2(s):
	return np.array([3 * s, -2 * s, s + 1])


def pathpos(pts):
	"""Distance of each point from node 0 along the path 0-1-2."""
	return pts.edges + pts.locs


class posstat(stat.base):
	"""Stat whose values are a function of each point's path position."""
	def __init__(self, func, names, label='Position stat'):
		self.func = func
		self.label0 = label
		super().__init__(names=names)

	def default_names(self):
		raise NotImplementedError

	def default_label(self):
		return self.label0

	def compute(self, pts):
		return np.asarray(self.func(pathpos(pts)), dtype=float)


class stub(panel.base):
	"""Panel that records its calls and returns placeholder artists."""
	def __init__(self, ax, pts, name, log):
		super().__init__(ax, pts)
		self.name = name
		self.log = log

	def init(self):
		self.log.append((self.name, 'init'))
		return [self.name + '-init']

	def draw(self, t, **ka):
		self.log.append((self.name, t, ka))
		return [self.name + '-draw{}'.format(t)]


class recorder(panel.base):
	"""Panel with one real line artist that records its calls."""
	def __init__(self, ax, pts, log):
		super().__init__(ax, pts)
		self.log = log
		self.line = None

	def init(self):
		self.log.append('init')
		self.line, = self.ax.plot([], [])
		return [self.line]

	def draw(self, t):
		self.log.append(t)
		self.line.set_data([0, t], [0, t])
		return [self.line]


CELL_W = np.array([[1, .5, 0, 0, 0, 0], [0, 0, 1, .25, 0, 0], [0, 0, 0, 0, .5, 1]])
COORD = np.array([[0., 1, 2, 3, 4, 5], [0., 1, 0, -1, 0, 1]])
CTYPE = np.array(['T2', 'T1', 'T1', 'T3', 'T2', 'T2'])


def make_network(cprop_extra=None, types=True):
	"""Dynamic network with 6 cells, 3 states at the nodes of the path 0-1-2, and 4 genes (2 regulators)."""
	traj = trajectory(np.array([[0, 1], [1, 2]]), np.array([1., 1.]))
	cells = point(traj, np.array([0, 0, 0, 0, 1, 1]), np.array([0., 0, 1, 1, 1, 1]))
	cprop = {'coord': COORD.copy()}
	if types:
		cprop['type'] = CTYPE.copy()
	cprop.update(cprop_extra or {})
	rng = np.random.default_rng(0)
	return dynamic_network(
		cname=np.array(['C{}'.format(x) for x in range(6)]), sname=np.array(['S0', 'S1', 'S2']),
		nname=np.array(['G1', 'G2', 'G3', 'G4']), nids=[np.array([0, 1]), np.array([0, 1, 2, 3])],
		traj=traj, point={'c': cells, 's': point.fromnodes(traj)}, cprop=cprop,
		scprop={'w': CELL_W.copy()}, ssprop={'traj-neighbor': np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])},
		nsprop={'cpm': rng.random((4, 3))}, esprop={'w': rng.standard_normal((2, 4, 3))})


def line_xy(line):
	"""x and y data of a Line2D with masked entries as nan."""
	return [np.ma.filled(np.ma.asarray(x, dtype=float), np.nan) for x in [line.get_xdata(), line.get_ydata()]]


def tick2_visible(axis):
	"""Whether the top (x axis) or right (y axis) major ticks are shown."""
	return axis.majorTicks[0].tick2line.get_visible()


def arrow_segment(arrow):
	"""(tail, tip) of a FancyArrow drawn with zero width and a small head; the tip is the end nearer the head corners."""
	xy = np.unique(np.round(np.asarray(arrow.get_xy(), dtype=float), 9), axis=0)
	dist = np.linalg.norm(xy[:, None] - xy[None], axis=2)
	i, j = np.unravel_index(dist.argmax(), dist.shape)
	a, b = xy[i], xy[j]
	u = (b - a) / np.linalg.norm(b - a)
	corner = xy[np.abs((xy - a) @ np.array([-u[1], u[0]])).argmax()]
	if np.linalg.norm(corner - b) < np.linalg.norm(corner - a):
		return a, b
	return b, a


@pytest.fixture(autouse=True)
def close_figures():
	"""Close every figure after each test."""
	yield
	plt.close('all')


@pytest.fixture
def traj():
	return trajectory(np.array([[0, 1], [1, 2]]), np.array([1., 1.]))


@pytest.fixture
def pts(traj):
	return traj.linspace(0, 2, 5)


@pytest.fixture
def ax():
	return plt.subplots()[1]


@pytest.fixture
def statx():
	return posstat(fx, [NAMES], 'X label')


@pytest.fixture
def staty():
	return posstat(fy, [NAMES], 'Y label')


@pytest.fixture
def staty2():
	return posstat(fy2, [NAMES], 'Y2 label')


def test_fixture_points_positions(pts):
	"""The test points lie at path positions 0, 0.5, 1, 1.5 and 2."""
	assert_allclose(pathpos(pts), S)


# base

def test_base_is_abstract(ax, pts):
	"""panel.base cannot be instantiated without init and draw."""
	with pytest.raises(TypeError):
		panel.base(ax, pts)


def test_base_stores_axes_and_points(ax, pts):
	"""A concrete subclass keeps the axes and points it was given."""
	p = stub(ax, pts, 'x', [])
	assert p.ax is ax
	assert p.pts is pts


# overlay

def test_overlay_rejects_empty_or_non_panel(ax, pts):
	"""overlay needs at least one panel and only panel.base instances."""
	with pytest.raises(AssertionError):
		panel.overlay(ax, pts, [])
	with pytest.raises(AssertionError):
		panel.overlay(ax, pts, [stub(ax, pts, 'a', []), object()])


def test_overlay_init_and_draw_concatenate_in_order(ax, pts):
	"""overlay.init and draw call each panel in order, pass t and keywords on, and chain the artist lists."""
	log = []
	p = panel.overlay(ax, pts, [stub(ax, pts, 'a', log), stub(ax, pts, 'b', log)])
	assert p.ax is ax and p.pts is pts and len(p.panels) == 2
	assert p.init() == ['a-init', 'b-init']
	assert p.draw(3, force=True) == ['a-draw3', 'b-draw3']
	assert log == [('a', 'init'), ('b', 'init'), ('a', 3, {'force': True}), ('b', 3, {'force': True})]


# statscatter: construction

def test_statscatter_default_names_are_sorted_intersection(ax, pts):
	"""Without names, statscatter shows the sorted names shared by statx, staty, staty2 and every statka stat."""
	sx = stat.const(np.arange(4.), [['c', 'a', 'b', 'd']])
	sy = stat.const(np.arange(3.), [['b', 'c', 'a']])
	sy2 = stat.const(np.arange(4.), [['a', 'b', 'c', 'e']])
	sc = stat.const(np.zeros((3, 4)), [['b', 'a', 'x'], list('RGBA')])
	p = panel.statscatter(ax, pts, sx, sy, staty2=sy2, statka={'color': sc})
	assert p.names == ['a', 'b']
	assert p.namesdict == {'a': 0, 'b': 1}
	assert p.ny == 2
	assert p.stats == [sx, sy, sy2]
	p = panel.statscatter(ax, pts, sx, sy)
	assert p.names == ['a', 'b', 'c']
	assert p.ny == 1
	assert p.stats == [sx, sy]


def test_statscatter_explicit_names(ax, pts, statx, staty):
	"""Explicit names keep their order; names absent from a stat raise ValueError and an empty list is rejected."""
	p = panel.statscatter(ax, pts, statx, staty, names=['c', 'a'])
	assert p.names == ['c', 'a']
	assert p.namesdict == {'c': 0, 'a': 1}
	with pytest.raises(ValueError):
		panel.statscatter(ax, pts, statx, staty, names=['a', 'zz'])
	with pytest.raises(AssertionError):
		panel.statscatter(ax, pts, statx, staty, names=[])


def test_statscatter_annotate_forms(ax, pts, statx, staty):
	"""annotate as a list maps name index to the name, as a dict to the given label, and 'all' covers every name."""
	assert panel.statscatter(ax, pts, statx, staty, annotate=['c']).annotate == {2: 'c'}
	assert panel.statscatter(ax, pts, statx, staty, annotate={'b': 'B', 'a': 'A'}).annotate == {1: 'B', 0: 'A'}
	assert panel.statscatter(ax, pts, statx, staty, names=['c', 'a'], annotate=['a']).annotate == {1: 'a'}
	assert panel.statscatter(ax, pts, statx, staty, annotate='all').annotate == {0: 'a', 1: 'b', 2: 'c'}
	assert panel.statscatter(ax, pts, statx, staty).annotate == {}


def test_statscatter_annotate_missing_names(ax, pts, statx, staty):
	"""Unknown names to annotate raise ValueError, or are dropped with annotate_err=False."""
	with pytest.raises(ValueError, match='zz'):
		panel.statscatter(ax, pts, statx, staty, annotate=['a', 'zz'])
	p = panel.statscatter(ax, pts, statx, staty, annotate=['a', 'zz'], annotate_err=False)
	assert p.annotate == {0: 'a'}
	p = panel.statscatter(ax, pts, statx, staty, annotate={'zz': 'Z', 'c': 'C'}, annotate_err=False)
	assert p.annotate == {2: 'C'}


def test_statscatter_rejects_invalid_stats(ax, pts, statx, staty):
	"""Non-stat statka values raise TypeError; multi-dimensional location stats and statka without names raise ValueError."""
	with pytest.raises(TypeError):
		panel.statscatter(ax, pts, statx, staty, statka={'color': 'red'})
	s2d = stat.const(np.zeros((3, 2)), [NAMES, ['u', 'v']])
	with pytest.raises(ValueError):
		panel.statscatter(ax, pts, s2d, staty)
	with pytest.raises(ValueError):
		panel.statscatter(ax, pts, statx, s2d)
	with pytest.raises(ValueError):
		panel.statscatter(ax, pts, statx, staty, staty2=s2d)
	empty = stat.const(np.zeros(3), [NAMES])
	empty.names = []
	with pytest.raises(ValueError):
		panel.statscatter(ax, pts, statx, staty, statka={'alpha': empty})


def test_statscatter_lim_list_must_match_axis_count(ax, pts, statx, staty, staty2):
	"""A list lim needs one entry per axis: two without staty2 and three with it."""
	with pytest.raises(ValueError):
		panel.statscatter(ax, pts, statx, staty, lim=[None])
	with pytest.raises(ValueError):
		panel.statscatter(ax, pts, statx, staty, staty2=staty2, lim=[None, None])
	p = panel.statscatter(ax, pts, statx, staty, staty2=staty2, lim=[None, None, None])
	assert p.lim == [None, None, None]


def test_statscatter_stores_options(ax, pts, statx, staty, staty2):
	"""scatterka, scatterka2, statka and aspect are kept as given."""
	sa = posstat(lambda s: np.array([s / 2] * 3), [NAMES])
	p = panel.statscatter(ax, pts, statx, staty, staty2=staty2, scatterka={'s': 3}, scatterka2={'s': 4},
		statka={'alpha': sa}, aspect=2)
	assert p.scatterka == {'s': 3}
	assert p.scatterka2 == {'s': 4}
	assert p.statka == {'alpha': sa}
	assert p.aspect == 2


# statscatter: limits

def test_statscatter_autolimits_default(ax, pts, statx, staty):
	"""Default limits span each stat's finite min and max over all points, expanded by 2% of the range per side."""
	# x: values s, 2s, -s for s in [0,2] give [-2,4]; y: s^2, 1+s and s (s<=1) give [0,4]
	p = panel.statscatter(ax, pts, statx, staty)
	assert_allclose(p.autolimits(), [[-2.12, 4.12], [-0.08, 4.08]])


def test_statscatter_autolimits_use_shown_names(ax, pts, statx, staty):
	"""Limits only consider the names shown."""
	p = panel.statscatter(ax, pts, statx, staty, names=['a'])
	assert_allclose(p.autolimits(), [[-0.04, 2.04], [-0.08, 4.08]])


def test_statscatter_autolimits_fall_back_to_unit_range(ax, pts, statx):
	"""A stat with a single value or with no finite value gets limits [0,1]."""
	flat = stat.const(np.ones(3), [NAMES])
	allnan = posstat(lambda s: np.full((3, len(s)), np.nan), [NAMES])
	p = panel.statscatter(ax, pts, flat, allnan)
	assert_allclose(p.autolimits(), [[0, 1], [0, 1]])
	p = panel.statscatter(ax, pts, statx, flat)
	assert_allclose(p.autolimits(), [[-2.12, 4.12], [0, 1]])


@pytest.mark.parametrize('lim,expected', [
	(set(), [[-2.12, 4.12], [-0.08, 4.08]]),
	({'sym'}, [[-4.12, 4.12], [-4.08, 4.08]]),
	({'min'}, [[-2.12, 4.12], [-2.12, 4.08]]),
	({'max'}, [[-2.12, 4.12], [-0.08, 4.12]]),
	({'min', 'max'}, [[-2.12, 4.12], [-2.12, 4.12]]),
	({'sym', 'min', 'max'}, [[-4.12, 4.12], [-4.12, 4.12]]),
])
def test_statscatter_autolimits_set_options(ax, pts, statx, staty, lim, expected):
	"""'sym' makes each axis symmetric about 0, 'min' and 'max' share the lowest minimum and highest maximum across axes."""
	p = panel.statscatter(ax, pts, statx, staty, lim=lim)
	assert_allclose(p.autolimits(), expected)


def test_statscatter_autolimits_set_options_with_second_y(ax, pts, statx, staty, staty2):
	"""The set options cover the second Y axis too; its default range [-4,6] expands to [-4.2,6.2]."""
	p = panel.statscatter(ax, pts, statx, staty, staty2=staty2)
	assert_allclose(p.autolimits(), [[-2.12, 4.12], [-0.08, 4.08], [-4.2, 6.2]])
	p = panel.statscatter(ax, pts, statx, staty, staty2=staty2, lim={'sym', 'min', 'max'})
	assert_allclose(p.autolimits(), [[-6.2, 6.2]] * 3)


def test_statscatter_autolimits_unknown_option(ax, pts, statx, staty):
	"""An unknown option in a set lim raises ValueError."""
	p = panel.statscatter(ax, pts, statx, staty, lim={'sym', 'wide'})
	with pytest.raises(ValueError, match='wide'):
		p.autolimits()


def test_statscatter_autolimits_explicit_values(ax, pts, statx, staty):
	"""Explicit limits override the automatic ones entry by entry; None entries stay automatic."""
	p = panel.statscatter(ax, pts, statx, staty, lim=np.array([[0, 1], [None, 2]], dtype=object))
	assert_allclose(p.autolimits(), [[0, 1], [-0.08, 2]])
	for lim in [[None, None], [[None, None], None]]:
		p = panel.statscatter(ax, pts, statx, staty, lim=lim)
		assert_allclose(p.autolimits(), [[-2.12, 4.12], [-0.08, 4.08]])


# statscatter: drawing

def test_statscatter_init_sets_up_axes(ax, pts, statx, staty):
	"""init sets limits and labels, adds an empty text per annotation, turns on top/right ticks and draws frame 0."""
	p = panel.statscatter(ax, pts, statx, staty, annotate=['a', 'c'], scatterka={'s': 7, 'zorder': 5}, aspect=2)
	objs = p.init()
	assert_allclose(ax.get_xlim(), [-2.12, 4.12])
	assert_allclose(ax.get_ylim(), [-0.08, 4.08])
	assert ax.get_xlabel() == 'X label'
	assert ax.get_ylabel() == 'Y label'
	assert len(ax.collections) == 1 and len(ax.texts) == 2
	sc = ax.collections[0]
	assert objs == [sc] + list(ax.texts)
	assert_allclose(sc.get_sizes(), [7])
	assert sc.get_zorder() == 5
	assert_allclose(sc.get_offsets(), [[0, 0], [0, 1], [0, 0]])
	assert [x.get_text() for x in ax.texts] == ['a', 'c']
	assert ax.get_aspect() == 2
	assert tick2_visible(ax.xaxis) and tick2_visible(ax.yaxis)


def test_statscatter_init_keeps_ticks_off_hidden_spines(ax, pts, statx, staty):
	"""Top and right ticks stay off when those spines are hidden, and aspect stays automatic by default."""
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	panel.statscatter(ax, pts, statx, staty).init()
	assert not tick2_visible(ax.xaxis)
	assert not tick2_visible(ax.yaxis)
	assert ax.get_aspect() == 'auto'


def test_statscatter_draw_positions_and_annotations(ax, pts, statx, staty):
	"""draw(t) moves dots to (x,y) at frame t and shows each annotation only where both coordinates are finite."""
	p = panel.statscatter(ax, pts, statx, staty, annotate={'a': 'A', 'c': 'C'})
	p.init()
	sc = ax.collections[0]
	objs = p.draw(1)
	assert objs == [sc] + list(ax.texts)
	assert_allclose(sc.get_offsets(), [[.5, .25], [1, 1.5], [-.5, .5]])
	assert [x.get_text() for x in ax.texts] == ['A', 'C']
	assert_allclose([x.get_position() for x in ax.texts], [[.5, .25], [-.5, .5]])
	p.draw(4)
	assert_allclose(sc.get_offsets(), [[2, 4], [4, 3], [-2, np.nan]])
	assert [x.get_text() for x in ax.texts] == ['A', '']
	assert_allclose(ax.texts[0].get_position(), [2, 4])


def test_statscatter_draw_follows_names_order(ax, pts, statx, staty):
	"""Dots follow the order of names."""
	p = panel.statscatter(ax, pts, statx, staty, names=['c', 'a'])
	p.init()
	p.draw(2)
	assert_allclose(ax.collections[0].get_offsets(), [[-1, 1], [1, 1]])


def test_statscatter_draw_sets_stat_parameters(ax, pts, statx, staty):
	"""statka stats set scatter properties per dot; constant ones only when forced, at init."""
	colors = to_rgba_array(['red', 'green', 'blue'])
	sc_color = stat.const(colors, [NAMES, list('RGBA')])
	sc_alpha = posstat(lambda s: np.array([s / 2, s / 4, 1 - s / 2]), [NAMES])
	p = panel.statscatter(ax, pts, statx, staty, statka={'color': sc_color, 'alpha': sc_alpha})
	p.init()
	sc = ax.collections[0]
	assert_allclose(sc.get_facecolors()[:, :3], colors[:, :3])
	assert_allclose(sc.get_facecolors()[:, 3], [0, 0, 1])
	assert p.draw(2) == [sc]
	assert_allclose(sc.get_facecolors()[:, :3], colors[:, :3])
	assert_allclose(sc.get_facecolors()[:, 3], [.5, .25, .5])


def test_statscatter_draw_skips_constant_stats(ax, pts):
	"""With only constant stats, draw redraws nothing unless forced."""
	sx = stat.const(np.array([1., 2]), [['a', 'b']])
	sy = stat.const(np.array([3., 5]), [['a', 'b']])
	p = panel.statscatter(ax, pts, sx, sy, annotate=['b'])
	objs = p.init()
	assert len(objs) == 2
	assert_allclose(ax.collections[0].get_offsets(), [[1, 3], [2, 5]])
	assert p.draw(3) == []
	assert len(p.draw(3, force=True)) == 2
	assert_allclose(ax.texts[0].get_position(), [2, 5])


def test_statscatter_second_y_axis(ax, pts, statx, staty, staty2):
	"""staty2 draws a second scatter with its own annotations on a twin axis with its own limits and label."""
	p = panel.statscatter(ax, pts, statx, staty, staty2=staty2, annotate=['b'], scatterka2={'s': 11})
	objs = p.init()
	assert p.ax2 is not ax and p.ax2 in ax.figure.axes
	assert_allclose(p.ax2.get_ylim(), [-4.2, 6.2])
	assert p.ax2.get_ylabel() == 'Y2 label'
	sc2 = p.ax2.collections[0]
	assert_allclose(sc2.get_sizes(), [11])
	assert objs == [ax.collections[0], sc2, ax.texts[0], p.ax2.texts[0]]
	p.draw(2)
	assert_allclose(sc2.get_offsets(), [[1, 3], [2, -2], [-1, 2]])
	assert p.ax2.texts[0].get_text() == 'b'
	assert_allclose(p.ax2.texts[0].get_position(), [2, -2])
	assert_allclose(ax.texts[0].get_position(), [2, 2])


def test_statscatter_get_data(ax, pts, statx, staty):
	"""get_data stacks location stats as (stat, name, point) in names order, and statka values as (name, ..., point)."""
	sc_alpha = posstat(lambda s: np.array([s / 2, s / 4, 1 - s / 2]), [NAMES])
	p = panel.statscatter(ax, pts, statx, staty, names=['c', 'a'], statka={'alpha': sc_alpha})
	data, param = p.get_data(pts)
	assert data.shape == (2, 2, 5)
	assert_allclose(data[0], fx(S)[[2, 0]])
	assert_allclose(data[1], fy(S)[[2, 0]])
	assert list(param) == ['alpha']
	assert_allclose(param['alpha'], np.array([S / 2, S / 4, 1 - S / 2])[[2, 0]])
	data, param = p.get_data(pts[[3]])
	assert_allclose(data[:, :, 0], [[-1.5, 1.5], [np.nan, 2.25]])


def test_statscatter_get_data_constant_stats(ax, pts):
	"""get_data returns no data and no constant parameters unless forced."""
	sx = stat.const(np.array([1., 2]), [['a', 'b']])
	sy = stat.const(np.array([3., 5]), [['a', 'b']])
	sc = stat.const(np.array([[1., 0, 0, 1], [0, 0, 1, 1]]), [['b', 'a'], list('RGBA')])
	p = panel.statscatter(ax, pts, sx, sy, statka={'color': sc})
	data, param = p.get_data(pts)
	assert data is None
	assert param == {}
	data, param = p.get_data(pts[[0, 1]], force=True)
	assert_allclose(data, [[[1, 1], [2, 2]], [[3, 3], [5, 5]]])
	assert param['color'].shape == (2, 4, 2)
	assert_allclose(param['color'][:, :, 0], [[0, 0, 1, 1], [1, 0, 0, 1]])


# statplot_static

def test_statplot_static_init_draws_lines(ax, pts, statx, staty):
	"""init draws one line per name over all points, masking non-finite values, labelled at its last shown point."""
	colors = ['red', 'green', 'blue']
	p = panel.statplot_static(ax, pts, statx, staty, colors)
	assert p.colors is colors
	assert p.init() == []
	lines = ax.get_lines()
	assert len(lines) == 3
	for xi in range(3):
		x, y = line_xy(lines[xi])
		assert_allclose(x, fx(S)[xi])
		assert_allclose(y, fy(S)[xi])
		assert to_rgba(lines[xi].get_color()) == to_rgba(colors[xi])
	assert [x.get_text() for x in ax.texts] == NAMES
	assert_allclose([x.get_position() for x in ax.texts], [[2, 4], [4, 3], [-1, 1]])
	assert ax.get_xlabel() == 'X label' and ax.get_ylabel() == 'Y label'
	assert tick2_visible(ax.xaxis) and tick2_visible(ax.yaxis)


def test_statplot_static_skips_hidden_names_and_keeps_color_index(ax, pts, statx):
	"""Names with no finite point get no line; colors stay indexed by position in names."""
	sy = posstat(lambda s: np.array([s, np.full_like(s, np.nan), np.where(s < 1, np.nan, s)]), [NAMES])
	panel.statplot_static(ax, pts, statx, sy, ['red', 'green', 'blue']).init()
	lines = ax.get_lines()
	assert len(lines) == 2
	assert [to_rgba(x.get_color()) for x in lines] == [to_rgba('red'), to_rgba('blue')]
	assert [x.get_text() for x in ax.texts] == ['a', 'c']
	assert_allclose(line_xy(lines[1])[1], [np.nan, np.nan, 1, 1.5, 2])


def test_statplot_static_names_plot_kwargs_and_second_y(ax, pts, statx, staty, staty2):
	"""names, plotka and the second Y axis with plotka2 are honoured."""
	p = panel.statplot_static(ax, pts, statx, staty, ['red', 'green'], names=['c', 'a'], staty2=staty2,
		plotka={'lw': 3}, plotka2={'lw': 0.5})
	assert p.plotka == {'lw': 3} and p.plotka2 == {'lw': 0.5}
	p.init()
	lines = ax.get_lines()
	assert len(lines) == 2
	assert_allclose(line_xy(lines[0])[0], fx(S)[2])
	assert [x.get_linewidth() for x in lines] == [3, 3]
	assert [to_rgba(x.get_color()) for x in lines] == [to_rgba('red'), to_rgba('green')]
	lines2 = p.ax2.get_lines()
	assert len(lines2) == 2
	assert_allclose(line_xy(lines2[1])[1], fy2(S)[0])
	assert [x.get_linewidth() for x in lines2] == [0.5, 0.5]
	assert [x.get_text() for x in p.ax2.texts] == ['c', 'a']
	assert_allclose([x.get_position() for x in p.ax2.texts], [[-2, 3], [2, 6]])
	assert p.ax2.get_ylabel() == 'Y2 label'


def test_statplot_static_draw_changes_nothing(ax, pts, statx, staty):
	"""draw returns no artists and leaves the lines as drawn by init."""
	p = panel.statplot_static(ax, pts, statx, staty, ['red', 'green', 'blue'])
	p.init()
	before = [line_xy(x) for x in ax.get_lines()]
	assert p.draw(3) == []
	assert_allclose([line_xy(x) for x in ax.get_lines()], before)


# statplot

def test_statplot_builds_pointer_and_lines(ax, pts, statx, staty):
	"""statplot overlays a pointer scatter on static lines, coloured by the first colors of tab10."""
	p = panel.statplot(ax, pts, statx, staty, names=NAMES)
	assert [type(x) for x in p.panels] == [panel.statscatter, panel.statplot_static]
	assert_allclose(p.panels[1].colors, to_rgba_array(TAB10[:3]))
	assert p.panels[0].scatterka == {'s': 20, 'zorder': 99}
	assert p.panels[1].plotka == {}
	p = panel.statplot(ax, pts, statx, staty, names=NAMES, pointer=False, plotka={'lw': 2})
	assert [type(x) for x in p.panels] == [panel.statplot_static]
	assert p.panels[0].plotka == {'lw': 2}


def test_statplot_pointer_keywords_and_cmap(ax, pts, statx, staty):
	"""pointerka updates the pointer defaults and an array cmap gives the colors directly."""
	cmap = to_rgba_array(['red', 'black'])
	p = panel.statplot(ax, pts, statx, staty, names=['b', 'a'], cmap=cmap, pointerka={'s': 5, 'marker': 'x'})
	assert p.panels[0].scatterka == {'s': 5, 'zorder': 99, 'marker': 'x'}
	assert_allclose(p.panels[1].colors, cmap)
	p = panel.statplot(ax, pts, statx, staty, names=['b', 'a'], cmap='Set1')
	assert_allclose(p.panels[1].colors, to_rgba_array(matplotlib.colormaps['Set1'].colors[:2]))


def test_statplot_draws_lines_and_moving_pointer(ax, pts, statx, staty):
	"""After init, each name has a line, and draw(t) moves its pointer of the same color to frame t."""
	p = panel.statplot(ax, pts, statx, staty, names=NAMES)
	objs = p.init()
	sc = ax.collections[0]
	assert objs == [sc]
	assert len(ax.get_lines()) == 3
	assert_allclose(sc.get_sizes(), [20])
	assert sc.get_zorder() == 99
	assert_allclose(sc.get_facecolors(), to_rgba_array(TAB10[:3]))
	assert_allclose(sc.get_offsets(), [[0, 0], [0, 1], [0, 0]])
	assert p.draw(3) == [sc]
	assert_allclose(sc.get_offsets(), [[1.5, 2.25], [3, 2.5], [-1.5, np.nan]])


# cellscatter_scatter, cellscatter_pointer, cellscatter

@pytest.fixture
def cellstats():
	"""Constant cell coordinates and a cell weight stat that changes along the path."""
	d = make_network()
	sx = stat.const(COORD[0], [d.cname], label='Dim1')
	sy = stat.const(COORD[1], [d.cname], label='Dim2')
	sw = posstat(lambda s: np.array([s / 2, 1 - s / 2, 0 * s + 1, s / 4, 0 * s, 0 * s + .5]), [d.cname])
	return d, sx, sy, sw


def cell_weights(s):
	return np.array([s / 2, 1 - s / 2, 1, s / 4, 0, .5])


def test_cellscatter_scatter_colors_by_type(ax, pts, cellstats):
	"""Cells are coloured by sorted cell type with tab10, with alpha from weights, and a colour legend is drawn."""
	d, sx, sy, sw = cellstats
	p = panel.cellscatter_scatter(ax, d, pts, sx, sy, sw, alphas=(0.2, 0.8), legend_loc=(0.3, 0.9),
		legend_ka={'fontsize': 6}, s=3)
	assert list(p.cmap) == ['T1', 'T2', 'T3']
	assert p.d is d and p.legend_loc == (0.3, 0.9) and p.legend_ka == {'fontsize': 6}
	objs = p.init()
	assert not ax.axison
	assert [x.get_text() for x in ax.texts] == ['T1', 'T2', 'T3']
	assert [to_rgba(x.get_color()) for x in ax.texts] == [to_rgba(x) for x in TAB10[:3]]
	assert all(x.get_fontsize() == 6 for x in ax.texts)
	assert_allclose(ax.texts[0].get_position(), (0.3, 0.9))
	sc = ax.collections[0]
	assert objs == [sc]
	assert_allclose(sc.get_sizes(), [3])
	assert_allclose(sc.get_offsets(), COORD.T)
	expected = to_rgba_array([TAB10[{'T1': 0, 'T2': 1, 'T3': 2}[x]] for x in CTYPE])
	assert_allclose(sc.get_facecolors()[:, :3], expected[:, :3])
	assert_allclose(sc.get_facecolors()[:, 3], 0.2 + 0.6 * cell_weights(0))
	assert p.draw(4) == [sc]
	assert_allclose(sc.get_facecolors()[:, 3], 0.2 + 0.6 * cell_weights(2))
	assert ax.get_aspect() == 1


def test_cellscatter_scatter_default_alphas_and_aspect(ax, pts, cellstats):
	"""Default alphas map weight 0 to 0.05 and weight 1 to 0.5; aspect is configurable."""
	d, sx, sy, sw = cellstats
	panel.cellscatter_scatter(ax, d, pts, sx, sy, sw, aspect=2).init()
	assert_allclose(ax.collections[0].get_facecolors()[:, 3], 0.05 + 0.45 * cell_weights(0))
	assert ax.get_aspect() == 2


def test_cellscatter_scatter_without_types_uses_gray(ax, pts, cellstats):
	"""Without the color property, or with a single type, every cell is gray."""
	d, sx, sy, sw = cellstats
	p = panel.cellscatter_scatter(ax, d, pts, sx, sy, sw, color='missing')
	assert list(p.cmap) == ['Unknown']
	assert_allclose(p.cmap['Unknown'], [.5, .5, .5, 1])
	p.init()
	assert [x.get_text() for x in ax.texts] == ['Unknown']
	assert_allclose(ax.collections[0].get_facecolors()[:, :3], np.full((6, 3), .5))
	d1 = make_network(cprop_extra={'group': np.array(['X'] * 6)}, types=False)
	p = panel.cellscatter_scatter(plt.subplots()[1], d1, pts, sx, sy, sw, color='group')
	assert list(p.cmap) == ['X']
	assert_allclose(p.cmap['X'], [.5, .5, .5, 1])


def test_cellscatter_scatter_cmap_options(ax, pts, cellstats):
	"""A dict cmap gives each type its color and must cover every type; a named cmap is sampled per type."""
	d, sx, sy, sw = cellstats
	cmap = {'T1': [1, 0, 0, 1], 'T2': [0, 1, 0, 1], 'T3': [0, 0, 1, 1]}
	p = panel.cellscatter_scatter(ax, d, pts, sx, sy, sw, cmap=cmap)
	p.init()
	assert_allclose(ax.collections[0].get_facecolors()[:, :3], [cmap[x][:3] for x in CTYPE])
	with pytest.raises(AssertionError):
		panel.cellscatter_scatter(ax, d, pts, sx, sy, sw, cmap={'T1': [1, 0, 0, 1]})
	p = panel.cellscatter_scatter(ax, d, pts, sx, sy, sw, cmap='Set1')
	assert_allclose([p.cmap[x] for x in ['T1', 'T2', 'T3']], to_rgba_array(matplotlib.colormaps['Set1'].colors[:3]))


def test_cellscatter_pointer_is_weighted_mean(ax, pts, cellstats):
	"""The pointer sits at the weight-averaged cell coordinates of each frame."""
	d, sx, sy, sw = cellstats
	p = panel.cellscatter_pointer(ax, pts, sx, sy, sw, color='k', s=15)
	assert p.names == ['pointer']
	assert p.scatterka == {'color': 'k', 's': 15}
	p.init()
	sc = ax.collections[0]
	assert_allclose(sc.get_sizes(), [15])
	for t in [0, 4]:
		p.draw(t)
		w = cell_weights(S[t])
		assert_allclose(sc.get_offsets(), [(COORD * w).sum(axis=1) / w.sum()])


def test_cellscatter_overlay_with_smoothed_weights(ax):
	"""cellscatter draws cells with alpha from max-normalised smoothed state weights and a pointer at their weighted mean."""
	d = make_network()
	pts = d.traj.linspace(0, 2, 5)
	fsmooth = functools.partial(stat.fsmooth, pts=d.traj, smoothen_func=('linear', [], {}))
	p = panel.cellscatter(ax, d, pts, fsmooth)
	assert [type(x) for x in p.panels] == [panel.cellscatter_pointer, panel.cellscatter_scatter]
	assert p.panels[1].scatterka == {'s': 2, 'lw': 0}
	assert p.panels[0].scatterka == {'color': 'k', 's': 20, 'zorder': 99}
	objs = p.init()
	pointer, cells = p.panels[0].objs[0], p.panels[1].objs[0]
	assert objs == [pointer, cells]
	assert_allclose(ax.get_xlim(), [-0.1, 5.1])
	assert_allclose(ax.get_ylim(), [-1.04, 1.04])
	assert_allclose(cells.get_offsets(), COORD.T)
	assert_allclose(cells.get_linewidths(), [0])
	assert_allclose(pointer.get_facecolors(), [to_rgba('k')])
	assert pointer.get_zorder() == 99
	# Point 1 is midway between states 0 and 1, so its weights are their average; point 4 is state 2.
	for t, w in [(1, (CELL_W[0] + CELL_W[1]) / 2), (4, CELL_W[2])]:
		p.draw(t)
		wn = w / w.max()
		assert_allclose(cells.get_facecolors()[:, 3], 0.05 + 0.45 * wn)
		assert_allclose(pointer.get_offsets(), [(COORD * wn).sum(axis=1) / wn.sum()])


def test_cellscatter_keywords_and_no_pointer(ax):
	"""scatterka and pointerka update the defaults (including cellscatter_scatter options); pointer=False drops the pointer."""
	d = make_network()
	pts = d.traj.linspace(0, 2, 5)
	fsmooth = functools.partial(stat.fsmooth, pts=d.traj, smoothen_func=('linear', [], {}))
	p = panel.cellscatter(ax, d, pts, fsmooth, scatterka={'s': 6, 'legend_loc': (0.6, 1)}, pointerka={'color': 'r'})
	assert p.panels[1].scatterka == {'s': 6, 'lw': 0}
	assert p.panels[1].legend_loc == (0.6, 1)
	assert p.panels[0].scatterka == {'color': 'r', 's': 20, 'zorder': 99}
	p = panel.cellscatter(ax, d, pts, fsmooth, pointer=False)
	assert [type(x) for x in p.panels] == [panel.cellscatter_scatter]


# statheatmap

def fheat(s):
	return np.array([[s, -s, 2 * s], [1 + 0 * s, s - 1, 3 * s]])


@pytest.fixture
def statnet():
	return posstat(fheat, [['r1', 'r2'], ['t1', 't2', 't3']], 'Net')


def test_statheatmap_default_names_annotate_and_lims(ax, pts, statnet):
	"""Defaults show and annotate every name, with symmetric color limits from the stat's expanded range."""
	# All values over the points span [-2,6]; 2% expansion gives [-2.16,6.16].
	h = panel.statheatmap(ax, pts, statnet)
	assert h.names == [['r1', 'r2'], ['t1', 't2', 't3']]
	assert h.namesdict == [{'r1': 0, 'r2': 1}, {'t1': 0, 't2': 1, 't3': 2}]
	assert h.annotate == h.names
	assert_allclose(h.lim, [-6.16, 6.16])
	assert h.ka == {'aspect': 'auto'}
	h = panel.statheatmap(ax, pts, statnet, cmap_sym=False)
	assert_allclose(h.lim, [-2.16, 6.16])
	h = panel.statheatmap(ax, pts, statnet, lim=(-1, 3), aspect=2, cmap='coolwarm')
	assert h.lim == [-1, 3]
	assert h.ka == {'aspect': 2, 'cmap': 'coolwarm'}


def test_statheatmap_rejects_unknown_names(ax, pts, statnet):
	"""Unknown regulators or targets raise ValueError; a stat that is not two-dimensional is rejected."""
	with pytest.raises(ValueError, match='Regulator'):
		panel.statheatmap(ax, pts, statnet, names=[['zz'], None])
	with pytest.raises(ValueError, match='Target'):
		panel.statheatmap(ax, pts, statnet, names=[None, ['t1', 'zz']])
	with pytest.raises(AssertionError):
		panel.statheatmap(ax, pts, posstat(fx, [NAMES]))


def test_statheatmap_get_data(ax, pts, statnet):
	"""get_data returns the stat at one point for the chosen rows and columns, in names order."""
	h = panel.statheatmap(ax, pts, statnet, names=[['r2', 'r1'], ['t3', 't1']])
	data = h.get_data(pts[[1]])
	assert data.size == 4
	assert_allclose(np.reshape(data, (2, 2)), [[1.5, 1], [1, .5]])
	with pytest.raises(AssertionError):
		h.get_data(pts[[0, 1]])


def test_statheatmap_init_and_draw(ax, pts, statnet):
	"""init shows frame 0 with the color limits, ticks and labels; draw(t) replaces the image data."""
	h = panel.statheatmap(ax, pts, statnet, names=[['r2', 'r1'], ['t3', 't1', 't2']], annotate=[['r1'], ['t2', 't3']],
		lim=[-2, 2], cmap='coolwarm', interpolation='nearest')
	objs = h.init()
	assert len(objs) == 1 and objs[0] is ax.images[0]
	im = objs[0]
	assert_allclose(np.squeeze(im.get_array()), fheat(0)[[1, 0]][:, [2, 0, 1]])
	assert_allclose(im.get_clim(), [-2, 2])
	assert im.get_cmap().name == 'coolwarm'
	assert im.get_interpolation() == 'nearest'
	assert ax.get_aspect() == 'auto'
	assert_allclose(ax.get_xticks(), [2, 0])
	assert [x.get_text() for x in ax.get_xticklabels()] == ['t2', 't3']
	assert all(x.get_rotation() == 90 for x in ax.get_xticklabels())
	assert_allclose(ax.get_yticks(), [1])
	assert [x.get_text() for x in ax.get_yticklabels()] == ['r1']
	assert ax.get_xlabel() == 'Target' and ax.get_ylabel() == 'Regulator'
	assert h.draw(3) == [im]
	assert_allclose(np.squeeze(im.get_array()), fheat(1.5)[[1, 0]][:, [2, 0, 1]])


# network_node, network_edge, network

def floc(s):
	return np.array([[s, 0 * s], [0 * s, s + 1], [1 + 0 * s, -s]])


def fnet(s):
	z = 0 * s
	e13 = (s > 0.75) * 1.
	e21 = -1. * (s < 1.25)
	e23 = (s < 1.75) * 1.
	return np.array([[z, z + 1, e13, z], [e21, z, e23, z]])


@pytest.fixture
def statloc():
	return posstat(floc, [['n1', 'n2', 'n3'], ['Dim 1', 'Dim 2']], 'Coordinates')


@pytest.fixture
def statedge():
	return posstat(fnet, [['n1', 'n2'], ['n1', 'n2', 'n3', 'n4']], 'Edges')


ARROWKA = {'width': 0, 'head_width': 0.1, 'head_length': 0.1, 'length_includes_head': True}


def test_network_node_draws_coordinates(ax, pts, statloc):
	"""network_node places a dot per node at its two coordinates, with equal aspect by default."""
	p = panel.network_node(ax, pts, statloc, annotate=['n2'], scatterka={'s': 4})
	assert p.names == ['n1', 'n2', 'n3']
	assert p.aspect == 1
	p.init()
	assert ax.get_aspect() == 1
	p.draw(2)
	assert_allclose(ax.collections[0].get_offsets(), floc(1.))
	assert_allclose(ax.collections[0].get_sizes(), [4])
	assert ax.texts[0].get_text() == 'n2'
	assert_allclose(ax.texts[0].get_position(), [0, 2])


def test_network_node_positional_names_and_aspect(ax, pts, statloc):
	"""Positional arguments go to statscatter (names first) and aspect can be changed."""
	p = panel.network_node(ax, pts, statloc, ['n3', 'n1'], aspect=3)
	assert p.names == ['n3', 'n1']
	p.init()
	assert ax.get_aspect() == 3
	assert_allclose(ax.collections[0].get_offsets(), [[1, 0], [0, 0]])


def test_network_edge_maps_names_to_locations(ax, pts, statloc, statedge):
	"""Edge rows and columns map to node locations; names without a location are left out."""
	p = panel.network_edge(ax, pts, statloc, statedge, color='g')
	assert p.map == [{0: 0, 1: 1}, {0: 0, 1: 1, 2: 2}]
	assert p.ka == {'color': 'g'}
	assert p.init() == []
	assert p.nlast == 0


def test_network_edge_draws_reuses_and_hides_arrows(ax, pts, statloc, statedge):
	"""draw(t) draws an arrow per nonzero edge, reuses arrows across frames and hides the unused ones."""
	p = panel.network_edge(ax, pts, statloc, statedge, color='g', **ARROWKA)
	p.init()

	def check(objs, nvisible, segs):
		assert all(x.get_alpha() != 0 for x in p.objs[:nvisible])
		assert all(x.get_alpha() == 0 for x in p.objs[nvisible:])
		for arrow, (tail, tip) in zip(p.objs, segs):
			t1 = arrow_segment(arrow)
			assert_allclose(t1[0], tail, atol=1e-9)
			assert_allclose(t1[1], tip, atol=1e-9)

	# s=0: n1=(0,0), n2=(0,1), n3=(1,0); edges n1->n2, n2->n1, n2->n3
	objs = p.draw(0)
	assert len(p.objs) == 3 and objs == p.objs and p.nlast == 3
	assert all(to_rgba(x.get_facecolor()) == to_rgba('g') for x in p.objs)
	check(objs, 3, [((0, 0), (0, 1)), ((0, 1), (0, 0)), ((0, 1), (1, 0))])
	# s=1: n1=(1,0), n2=(0,2), n3=(1,-1); edges n1->n2, n1->n3, n2->n1, n2->n3
	objs = p.draw(2)
	assert len(p.objs) == 4 and objs == p.objs and p.nlast == 4
	check(objs, 4, [((1, 0), (0, 2)), ((1, 0), (1, -1)), ((0, 2), (1, 0)), ((0, 2), (1, -1))])
	# s=2: n1=(2,0), n2=(0,3), n3=(1,-2); edges n1->n2, n1->n3
	objs = p.draw(4)
	assert len(p.objs) == 4 and objs == p.objs and p.nlast == 2
	check(objs, 2, [((2, 0), (0, 3)), ((2, 0), (1, -2))])
	# s=1.5: n1=(1.5,0), n2=(0,2.5), n3=(1,-1.5); edges n1->n2, n1->n3, n2->n3
	objs = p.draw(3)
	assert objs == p.objs[:3] and p.nlast == 3
	assert all(x.get_alpha() == 1 for x in p.objs[:3])
	check(objs, 3, [((1.5, 0), (0, 2.5)), ((1.5, 0), (1, -1.5)), ((0, 2.5), (1, -1.5))])
	assert len(ax.patches) == 4


def test_network_overlay(ax, pts, statloc, statedge):
	"""network overlays nodes over edges, merges nodeka scatterka without changing the caller's dict, and hides the axis."""
	nodeka = {'scatterka': {'s': 9}, 'annotate': ['n1']}
	p = panel.network(ax, pts, statloc, statedge, nodeka=nodeka, edgeka=dict(ARROWKA))
	assert nodeka == {'scatterka': {'s': 9}, 'annotate': ['n1']}
	assert [type(x) for x in p.panels] == [panel.network_edge, panel.network_node]
	assert p.panels[0].ka == ARROWKA
	assert p.panels[1].scatterka == {'s': 9}
	assert p.panels[1].annotate == {0: 'n1'}
	objs = p.init()
	assert not ax.axison
	assert objs == [ax.collections[0], ax.texts[0]]
	objs = p.draw(0)
	assert len(objs) == 5
	assert objs[:3] == ax.patches[:3]
	assert_allclose(ax.collections[0].get_offsets(), floc(0.))
	p = panel.network(ax, pts, statloc, statedge)
	assert p.panels[1].scatterka == {}
	assert p.panels[0].ka == {}


# animate_generic

def test_animate_generic_rejects_non_panels(pts):
	"""animate_generic only accepts panel.base instances."""
	fig = plt.figure()
	with pytest.raises(AssertionError):
		panel.animate_generic(pts, fig, [object()])


def test_animate_generic_init_runs_once(pts):
	"""init calls each panel's init once, appends post_init's artists, and returns nothing on later calls."""
	fig, ax = plt.subplots()
	log = []
	ca = panel.animate_generic(pts, fig, [stub(ax, pts, 'a', log), stub(ax, pts, 'b', log)])
	assert ca.pts is pts and ca.fig is fig and not ca.hasinit
	seen = []

	def post_init(x):
		seen.append(x)
		return ['post']

	assert ca.init(post_init=post_init) == ['a-init', 'b-init', 'post']
	assert seen == [ca] and ca.hasinit
	assert ca.init(post_init=post_init) == []
	assert log == [('a', 'init'), ('b', 'init')]
	assert seen == [ca]
	ca2 = panel.animate_generic(pts, fig, [stub(ax, pts, 'c', [])])
	assert ca2.init() == ['c-init']


def test_animate_generic_draw_chains_panels(pts):
	"""draw(t) chains the artists each panel redraws at frame t."""
	fig, ax = plt.subplots()
	log = []
	ca = panel.animate_generic(pts, fig, [stub(ax, pts, 'a', log), stub(ax, pts, 'b', log)])
	assert ca.draw(2) == ['a-draw2', 'b-draw2']
	assert log == [('a', 2, {}), ('b', 2, {})]


def test_animate_generic_animate_arguments(pts, monkeypatch):
	"""animate builds a FuncAnimation over one frame per point with draw, an init wrapper with post_init, blit and keywords."""
	calls = []

	class fake:
		def __init__(self, *a, **ka):
			calls.append((a, ka))

	monkeypatch.setattr(matplotlib.animation, 'FuncAnimation', fake)
	fig, ax = plt.subplots()
	log = []
	ca = panel.animate_generic(pts, fig, [stub(ax, pts, 'a', log)])
	anim = ca.animate(post_init=lambda x: ['post'], interval=40)
	assert isinstance(anim, fake)
	(a, ka), = calls
	assert a == (fig, ca.draw)
	assert ka['frames'] == 5 and ka['blit'] is True and ka['interval'] == 40
	assert ka['init_func']() == ['a-init', 'post']
	ca.animate(blit=False)
	assert calls[1][1]['blit'] is False
	assert calls[1][1]['init_func']() == []


def test_animate_generic_animation_saves_every_frame(pts, tmp_path):
	"""The animation runs init once and draws every frame in order when saved as a GIF with Pillow."""
	pytest.importorskip('PIL')
	from PIL import Image
	fig, ax = plt.subplots(figsize=(1, 1), dpi=20)
	log = []
	ca = panel.animate_generic(pts, fig, [recorder(ax, pts, log)])
	anim = ca.animate(interval=10)
	assert list(anim.new_frame_seq()) == [0, 1, 2, 3, 4]
	fo = str(tmp_path / 'anim.gif')
	anim.save(fo, writer=matplotlib.animation.PillowWriter(fps=5))
	assert log == ['init', 0, 1, 2, 3, 4]
	with Image.open(fo) as im:
		assert im.format == 'GIF'
