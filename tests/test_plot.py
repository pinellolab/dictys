"""
Tests for the helpers in dictys.plot: get_cmap, colorbar, colorlegend,
heatmap and dotplot.
"""

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose

matplotlib = pytest.importorskip('matplotlib')
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import QuadMesh  # noqa: E402
from matplotlib.colors import ListedColormap, Normalize, to_rgba, to_rgba_array  # noqa: E402

import dictys.plot as dplot  # noqa: E402


@pytest.fixture(autouse=True)
def close_figures():
	"""Close every matplotlib figure a test opened."""
	yield
	plt.close('all')


def _ticklabels(labels):
	return [x.get_text() for x in labels]


def _quadmesh(ax):
	return [x for x in ax.collections if isinstance(x, QuadMesh)][0]


def _size_inches(fig, ax):
	"""Width and height of an axes in inches."""
	pos = ax.get_position()
	size = fig.get_size_inches()
	return np.array([pos.width * size[0], pos.height * size[1]])


def _leaves(data, method, metric, optimal_ordering):
	"""Reference dendrogram leaf order of the rows of data, from scipy."""
	from scipy.cluster.hierarchy import dendrogram, linkage
	return dendrogram(linkage(data, method=method, metric=metric, optimal_ordering=optimal_ordering), no_plot=True)['leaves']


########################################################################
# get_cmap
########################################################################

def test_get_cmap_fewer_than_listed_colors():
	"""Asking for fewer colors than a qualitative colormap holds returns its first colors."""
	ans = dplot.get_cmap('tab10', 3)
	assert ans.shape == (3, 4)
	assert_allclose(ans, to_rgba_array(plt.get_cmap('tab10').colors[:3]))


def test_get_cmap_all_listed_colors():
	"""Asking for exactly as many colors as a listed colormap holds returns all of them in order."""
	ans = dplot.get_cmap('tab10', 10)
	assert_allclose(ans, to_rgba_array(plt.get_cmap('tab10').colors))


def test_get_cmap_more_than_listed_colors():
	"""Asking for more colors than listed samples the map evenly, repeating colors (index floor(k/(n-1)*N), last clipped)."""
	n = 12
	ans = dplot.get_cmap('tab10', n)
	idx = np.minimum(np.floor(np.arange(n) / (n - 1) * 10 + 1e-9).astype(int), 9)
	assert list(idx) == [0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9]
	assert_allclose(ans, to_rgba_array(np.array(plt.get_cmap('tab10').colors)[idx]))


def test_get_cmap_colormap_object():
	"""A Colormap object is accepted in place of a name."""
	cmap = ListedColormap(['red', 'green', 'blue'])
	ans = dplot.get_cmap(cmap, 2)
	assert_allclose(ans, to_rgba_array(['red', 'green']))


def test_get_cmap_many_colors_spans_colormap():
	"""Asking for more colors than a 256-entry map holds spans it from its first to its last color."""
	cmap = plt.get_cmap('coolwarm')
	ans = dplot.get_cmap('coolwarm', 300)
	assert ans.shape == (300, 4)
	assert_allclose(ans[0], cmap(0.0))
	assert_allclose(ans[-1], cmap(1.0))


########################################################################
# colorbar
########################################################################

def test_colorbar_from_name_vertical():
	"""A colormap name with vmin/vmax draws a vertical colorbar spanning [vmin,vmax] in a new figure."""
	fig, ax = dplot.colorbar('viridis', -2, 3)
	assert fig.axes == [ax]
	assert_allclose(fig.get_size_inches(), [0.15, 0.8])
	assert_allclose(ax.get_ylim(), [-2, 3])
	assert _quadmesh(ax).get_cmap().name == 'viridis'
	assert ax.get_title() == ''


def test_colorbar_horizontal_title_and_passthrough():
	"""Orientation, figsize, title, title_* and colorbar keyword arguments all take effect."""
	fig, ax = dplot.colorbar('magma', vmin=-2, vmax=3, figsize=(2, 0.3), orientation='horizontal', title='Strength', title_fontsize=7, label='lab', ticks=[-2, 0, 3])
	assert_allclose(fig.get_size_inches(), [2, 0.3])
	assert_allclose(ax.get_xlim(), [-2, 3])
	assert ax.get_title(loc='center') == 'Strength'
	assert ax.title.get_fontsize() == 7
	assert ax.get_xlabel() == 'lab'
	assert_allclose(ax.get_xticks(), [-2, 0, 3])
	assert _quadmesh(ax).get_cmap().name == 'magma'


def test_colorbar_scalarmappable_on_given_axes():
	"""A ScalarMappable is drawn as is on given axes, returning that axes' figure."""
	sm = matplotlib.cm.ScalarMappable(norm=Normalize(1, 5), cmap='coolwarm')
	fig0 = plt.figure()
	ax0 = fig0.add_axes([0, 0, 0.1, 1])
	fig, ax = dplot.colorbar(sm, ax=ax0, figsize=None)
	assert fig is fig0 and ax is ax0
	assert_allclose(ax.get_ylim(), [1, 5])
	assert _quadmesh(ax).get_cmap().name == 'coolwarm'


@pytest.mark.parametrize('ka', [{'vmin': 0}, {'vmax': 1}])
def test_colorbar_scalarmappable_with_limits_raises(ka):
	"""vmin or vmax together with a ScalarMappable raises TypeError."""
	sm = matplotlib.cm.ScalarMappable(norm=Normalize(1, 5), cmap='coolwarm')
	with pytest.raises(TypeError):
		dplot.colorbar(sm, **ka)


def test_colorbar_axes_and_figsize_raises():
	"""Passing both ax and a figsize raises ValueError."""
	fig = plt.figure()
	ax = fig.add_subplot(111)
	with pytest.raises(ValueError):
		dplot.colorbar('viridis', 0, 1, ax=ax, figsize=(1, 1))


def test_colorbar_title_keywords_without_title_raises():
	"""title_* keyword arguments without a title raise ValueError."""
	with pytest.raises(ValueError):
		dplot.colorbar('viridis', 0, 1, title_fontsize=5)


########################################################################
# colorlegend
########################################################################

def _legend(ka, labels=('a', 'bb', 'ccc', 'd', 'e'), colors=('r', 'g', 'b', 'k', 'm'), loc=(0.1, 0.9)):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ans = dplot.colorlegend(ax, loc, list(labels), list(colors), **ka)
	renderer = fig.canvas.get_renderer()
	return ans, ax, [x.get_window_extent(renderer) for x in ax.texts]


def test_colorlegend_texts_and_colors():
	"""Each label becomes one text in its own color, with extra keyword arguments passed to ax.text."""
	ans, ax, _ = _legend({'fontsize': 13})
	assert ans is None
	assert [x.get_text() for x in ax.texts] == ['a', 'bb', 'ccc', 'd', 'e']
	assert [to_rgba(x.get_color()) for x in ax.texts] == [to_rgba(x) for x in 'rgbkm']
	assert all(x.get_fontsize() == 13 for x in ax.texts)


def test_colorlegend_single_column_spacing():
	"""In one column the first text's top-left corner sits at loc and each next text is lower by height*(1+spacey)."""
	_, ax, ext = _legend({'spacey': 0.5})
	x0, y0 = ax.transAxes.transform((0.1, 0.9))
	assert_allclose([ext[0].x0, ext[0].y1], [x0, y0], atol=1e-6)
	for xi in range(1, 5):
		assert_allclose(ext[xi].x0, x0, atol=1e-6)
		assert_allclose(ext[xi].y1, ext[xi - 1].y1 - ext[xi - 1].height * 1.5, atol=1e-6)


def test_colorlegend_columns():
	"""ncol=2 splits 5 labels 3+2; the second column starts right of the widest first-column text by spacex dots."""
	_, ax, ext = _legend({'ncol': 2, 'spacex': 10, 'spacey': 0.2})
	widthmax = max(x.width for x in ext[:3])
	assert_allclose(ext[3].x0, ext[0].x0 + widthmax + 10, atol=1e-6)
	assert_allclose(ext[3].y1, ext[0].y1, atol=1e-6)
	assert_allclose(ext[4].x0, ext[3].x0, atol=1e-6)
	assert_allclose(ext[4].y1, ext[3].y1 - ext[3].height * 1.2, atol=1e-6)
	assert_allclose(ext[2].y1, ext[1].y1 - ext[1].height * 1.2, atol=1e-6)


def test_colorlegend_space_sets_both_spacings():
	"""space overrides spacex and spacey with the same value."""
	_, _, ext1 = _legend({'ncol': 2, 'space': 3})
	_, _, ext2 = _legend({'ncol': 2, 'spacex': 3, 'spacey': 3, 'space': None})
	_, _, ext3 = _legend({'ncol': 2, 'spacex': 0, 'spacey': 0.15})
	assert_allclose([[x.x0, x.y1] for x in ext1], [[x.x0, x.y1] for x in ext2], atol=1e-6)
	assert not np.allclose([[x.x0, x.y1] for x in ext1], [[x.x0, x.y1] for x in ext3])


def test_colorlegend_alignment():
	"""horizontalalignment and verticalalignment anchor the first text's right/bottom edges at loc."""
	_, ax, ext = _legend({'horizontalalignment': 'right', 'verticalalignment': 'bottom', 'spacey': 0})
	x0, y0 = ax.transAxes.transform((0.1, 0.9))
	assert all(x.get_horizontalalignment() == 'right' for x in ax.texts)
	assert all(x.get_verticalalignment() == 'bottom' for x in ax.texts)
	assert_allclose([ext[0].x1, ext[0].y0], [x0, y0], atol=1e-6)
	assert_allclose(ext[1].y0, ext[0].y0 - ext[0].height, atol=1e-6)


def test_colorlegend_length_mismatch_raises():
	"""Different numbers of labels and colors fail the length assertion."""
	fig = plt.figure()
	ax = fig.add_subplot(111)
	with pytest.raises(AssertionError):
		dplot.colorlegend(ax, (0, 1), ['a', 'b'], ['r'])


########################################################################
# heatmap
########################################################################

@pytest.fixture
def hdata():
	"""5x4 table whose columns c1/c3 and c2/c4 are near-duplicates, and rows r1/r4 and r2/r5 likewise."""
	v = np.array([
		[0.0, 5.0, 0.1, 5.1],
		[3.0, 1.0, 3.1, 1.2],
		[6.0, 6.0, 6.2, 6.1],
		[0.2, 5.1, 0.3, 5.0],
		[3.1, 1.1, 3.0, 1.0],
	])
	return pd.DataFrame(v, index=['r1', 'r2', 'r3', 'r4', 'r5'], columns=['c1', 'c2', 'c3', 'c4'])


def _adjacent(order, a, b):
	order = list(order)
	return abs(order.index(a) - order.index(b)) == 1


def test_heatmap_default_clusters_columns(hdata):
	"""By default only columns are clustered (ward/euclidean, optimal ordering) and drawn reordered with a top dendrogram and colorbar."""
	fig, cols, rows = dplot.heatmap(hdata)
	leaves = _leaves(hdata.T.values, 'ward', 'euclidean', True)
	assert list(cols) == list(hdata.columns[leaves])
	assert _adjacent(cols, 'c1', 'c3') and _adjacent(cols, 'c2', 'c4')
	assert list(rows) == list(hdata.index)
	assert len(fig.axes) == 3
	axmatrix = fig.axes[1]
	assert len(axmatrix.images) == 1
	im = axmatrix.images[0]
	assert_allclose(np.asarray(im.get_array()), hdata.loc[rows, cols].values)
	assert im.get_cmap().name == 'coolwarm'
	assert_allclose(im.get_clim(), [hdata.values.min(), hdata.values.max()])
	assert list(axmatrix.get_xticks()) == []
	assert _ticklabels(axmatrix.get_yticklabels()) == list(rows)
	assert axmatrix.get_ylim()[0] > axmatrix.get_ylim()[1]
	assert fig.axes[0].get_xticks().size == 0 and fig.axes[0].get_yticks().size == 0
	assert_allclose(fig.axes[2].get_ylim(), im.get_clim())


def test_heatmap_both_axes_separate_options(hdata):
	"""dright clusters rows; tuples of method/metric apply to x then y; optimal_ordering is passed to scipy."""
	fig, cols, rows = dplot.heatmap(hdata, dright=0.3, method=('average', 'single'), metric=('cityblock', 'euclidean'), optimal_ordering=False, xtick=True)
	assert list(cols) == list(hdata.columns[_leaves(hdata.T.values, 'average', 'cityblock', False)])
	assert list(rows) == list(hdata.index[_leaves(hdata.values, 'single', 'euclidean', False)])
	assert _adjacent(rows, 'r1', 'r4') and _adjacent(rows, 'r2', 'r5')
	assert len(fig.axes) == 4
	ax1, ax2, axmatrix, cax = fig.axes
	assert_allclose(np.asarray(axmatrix.images[0].get_array()), hdata.loc[rows, cols].values)
	assert _ticklabels(axmatrix.get_xticklabels()) == list(cols)
	assert all(x.get_rotation() == 90 for x in axmatrix.get_xticklabels())
	assert _ticklabels(axmatrix.get_yticklabels()) == list(rows)
	# inverty flips the right dendrogram together with the matrix
	assert ax1.get_ylim()[0] > ax1.get_ylim()[1]


def test_heatmap_inverty_false(hdata):
	"""inverty=False keeps the matrix and right dendrogram y axes increasing."""
	fig, _, _ = dplot.heatmap(hdata, dright=0.3, inverty=False)
	ax1, _, axmatrix, _ = fig.axes
	assert axmatrix.get_ylim()[0] < axmatrix.get_ylim()[1]
	assert ax1.get_ylim()[0] < ax1.get_ylim()[1]


def test_heatmap_no_clustering_no_colorbar(hdata):
	"""dtop=0, dright=0 and wcolorbar=0 draw only the matrix, in the input order."""
	fig, cols, rows = dplot.heatmap(hdata, dtop=0, dright=0, wcolorbar=0, ytick=False)
	assert list(cols) == list(hdata.columns)
	assert list(rows) == list(hdata.index)
	assert len(fig.axes) == 1
	axmatrix = fig.axes[0]
	assert_allclose(np.asarray(axmatrix.images[0].get_array()), hdata.values)
	assert list(axmatrix.get_yticks()) == []


def test_heatmap_square_distance_matrix():
	"""metric=None treats d as a distance matrix and orders rows and columns by the same leaves."""
	from scipy.spatial.distance import squareform
	names = ['a', 'b', 'c', 'd']
	dist = np.array([
		[0, 5, 1, 6],
		[5, 0, 6, 1.5],
		[1, 6, 0, 5.5],
		[6, 1.5, 5.5, 0],
	], dtype=float)
	d = pd.DataFrame(dist, index=names, columns=names)
	from scipy.cluster.hierarchy import dendrogram, linkage
	leaves = dendrogram(linkage(squareform(dist), method='average', optimal_ordering=True), no_plot=True)['leaves']
	order = [names[x] for x in leaves]
	assert _adjacent(order, 'a', 'c') and _adjacent(order, 'b', 'd')
	for dtop, dright in [(0.3, 0.3), (0.3, 0), (0, 0.3)]:
		fig, cols, rows = dplot.heatmap(d, metric=None, method='average', dtop=dtop, dright=dright, xtick=True)
		assert list(cols) == order and list(rows) == order
		axmatrix = [x for x in fig.axes if len(x.images) > 0][0]
		assert_allclose(np.asarray(axmatrix.images[0].get_array()), d.loc[order, order].values)
		assert _ticklabels(axmatrix.get_xticklabels()) == order
		plt.close('all')


def test_heatmap_dshow_limits_and_cmap(hdata):
	"""dshow is drawn in the order clustered from d, with the given vmin, vmax and cmap."""
	dshow = -hdata * 2
	fig, cols, rows = dplot.heatmap(hdata, dshow=dshow, vmin=-3, vmax=4, cmap='viridis')
	assert list(cols) == list(hdata.columns[_leaves(hdata.T.values, 'ward', 'euclidean', True)])
	im = fig.axes[1].images[0]
	assert_allclose(np.asarray(im.get_array()), dshow.loc[rows, cols].values)
	assert_allclose(im.get_clim(), [-3, 4])
	assert im.get_cmap().name == 'viridis'


def test_heatmap_select_and_tick_lists(hdata):
	"""xselect/yselect keep the chosen columns/rows; tick lists relabel them and None hides a label."""
	xselect = np.array([True, True, True, False])
	yselect = np.array([True, False, True, True, True])
	xtick = ['X1', None, 'X3', 'X4']
	ytick = ['Y1', 'Y2', 'Y3', None, 'Y5']
	fig, cols, rows = dplot.heatmap(hdata, xselect=xselect, yselect=yselect, xtick=xtick, ytick=ytick, dtop=0)
	assert list(cols) == ['c1', 'c2', 'c3']
	assert list(rows) == ['r1', 'r3', 'r4', 'r5']
	axmatrix = fig.axes[0]
	assert_allclose(np.asarray(axmatrix.images[0].get_array()), hdata.loc[rows, cols].values)
	assert list(axmatrix.get_xticks()) == [0, 2]
	assert _ticklabels(axmatrix.get_xticklabels()) == ['X1', 'X3']
	assert list(axmatrix.get_yticks()) == [0, 1, 3]
	assert _ticklabels(axmatrix.get_yticklabels()) == ['Y1', 'Y3', 'Y5']


def test_heatmap_cell_size_from_figscale_and_aspect(hdata):
	"""Each cell is figscale inches wide and figscale/aspect inches tall, whatever the dendrogram and edge widths."""
	fig, _, _ = dplot.heatmap(hdata, figscale=0.3, aspect=0.5, dtop=0.4, dright=0.2, wcolorbar=0.05, wedge=0.02)
	axmatrix = fig.axes[2]
	assert_allclose(_size_inches(fig, axmatrix), [4 * 0.3, 5 * 0.3 / 0.5])
	assert axmatrix.get_aspect() == pytest.approx(2)
	# wedge is the margin left and below the matrix and right of the colorbar
	assert_allclose([axmatrix.get_position().x0, axmatrix.get_position().y0], [0.02, 0.02])
	cax = fig.axes[3].get_position()
	assert_allclose([cax.x0, cax.width], [1 - 0.02 - 0.05, 0.05])


def test_heatmap_draws_on_given_figure(hdata):
	"""A given figure is drawn on and returned without resizing."""
	fig0 = plt.figure(figsize=(3, 2))
	fig, _, _ = dplot.heatmap(hdata, fig=fig0)
	assert fig is fig0
	assert_allclose(fig.get_size_inches(), [3, 2])
	assert len(fig.axes) == 3


@pytest.mark.parametrize('ka', [{'method': ('ward',)}, {'method': ('ward', 'single', 'average')}, {'metric': ('euclidean',)}])
def test_heatmap_wrong_option_length_raises(hdata, ka):
	"""method or metric tuples not of length 2 raise ValueError."""
	with pytest.raises(ValueError):
		dplot.heatmap(hdata, **ka)


########################################################################
# dotplot
########################################################################

@pytest.fixture
def ddata():
	"""Dot size and color tables of 3 rows by 4 columns."""
	ds = pd.DataFrame(np.arange(1, 13, dtype=float).reshape(3, 4), index=['r1', 'r2', 'r3'], columns=['a', 'b', 'c', 'd'])
	dc = pd.DataFrame(np.array([[0.5, 0.1, 0.9, 0.3], [0.2, 0.8, 0.4, 0.6], [0.7, 0.0, 1.0, 0.2]]), index=ds.index, columns=ds.columns)
	return ds, dc


def test_dotplot_dots(ddata):
	"""Each entry is a dot at (column, row) with size from ds and color from dc, ticks named by columns and index."""
	ds, dc = ddata
	fig, ax = dplot.dotplot(ds, dc)
	assert len(fig.axes) == 3 and fig.axes[0] is ax
	coll = ax.collections[0]
	offsets = np.asarray(coll.get_offsets())
	rr, cc = np.meshgrid(np.arange(3), np.arange(4), indexing='ij')
	assert_allclose(offsets, np.c_[cc.ravel(), rr.ravel()])
	assert_allclose(coll.get_sizes(), ds.values.ravel())
	assert_allclose(np.asarray(coll.get_array()), dc.values.ravel())
	assert_allclose(coll.get_clim(), [0, 1])
	assert coll.get_cmap().name == 'viridis'
	assert _ticklabels(ax.get_xticklabels()) == ['a', 'b', 'c', 'd']
	assert _ticklabels(ax.get_yticklabels()) == ['r1', 'r2', 'r3']
	assert_allclose(ax.get_xlim(), [-0.5, 3.5])
	assert_allclose(ax.get_ylim(), [-0.5, 2.5])
	assert ax.get_aspect() == 1
	assert_allclose(fig.axes[2].get_ylim(), [0, 1])
	assert not fig.axes[1].axison
	assert len(fig.axes[1].collections) == 0


def test_dotplot_size_legend(ddata):
	"""size_transform maps ds and the sizes legend, drawn in black at y=0,1.3,... with value labels."""
	ds, dc = ddata
	fig, ax = dplot.dotplot(ds, dc, size_transform=lambda x: 2 * np.asarray(x), sizes=[1, 5, 10])
	assert_allclose(ax.collections[0].get_sizes(), 2 * ds.values.ravel())
	ax2 = fig.axes[1]
	coll = ax2.collections[0]
	assert_allclose(np.asarray(coll.get_offsets()), [[0, 0], [0, 1.3], [0, 2.6]])
	assert_allclose(coll.get_sizes(), [2, 10, 20])
	assert_allclose(coll.get_facecolors()[0], to_rgba('k'))
	assert [x.get_text() for x in ax2.texts] == ['1', '5', '10']
	assert_allclose([x.get_position() for x in ax2.texts], [[0.7, 0], [0.7, 1.3], [0.7, 2.6]])


def test_dotplot_vrange_sets_colorbar(ddata):
	"""vrange sets the colorbar range, and a None end falls back to the min/max of dc."""
	ds, dc = ddata
	fig, _ = dplot.dotplot(ds, dc, vrange=[-1, 2])
	assert_allclose(fig.axes[2].get_ylim(), [-1, 2])
	plt.close('all')
	fig, _ = dplot.dotplot(ds, dc, vrange=[None, 3])
	assert_allclose(fig.axes[2].get_ylim(), [0, 3])


def test_dotplot_figsize(ddata):
	"""A scalar figsize scales to (figsize*(2*ncol+0.5), figsize*nrow); a pair is used as is."""
	ds, dc = ddata
	fig, _ = dplot.dotplot(ds, dc, figsize=0.5)
	assert_allclose(fig.get_size_inches(), [0.5 * 8.5, 0.5 * 3])
	fig, _ = dplot.dotplot(ds, dc, figsize=(4, 2))
	assert_allclose(fig.get_size_inches(), [4, 2])


def test_dotplot_given_figure_cmap_and_scatter_keywords(ddata):
	"""A given figure is drawn on; cmap colors dots and colorbar; extra keywords reach scatter."""
	ds, dc = ddata
	fig0 = plt.figure(figsize=(5, 3))
	fig, ax = dplot.dotplot(ds, dc, fig=fig0, cmap='magma', alpha=0.5)
	assert fig is fig0
	assert_allclose(fig.get_size_inches(), [5, 3])
	coll = ax.collections[0]
	assert coll.get_cmap().name == 'magma'
	assert coll.get_alpha() == 0.5
	assert _quadmesh(fig.axes[2]).get_cmap().name == 'magma'


def test_dotplot_mismatched_tables_raise(ddata):
	"""ds and dc with different shapes or labels fail the input assertions."""
	ds, dc = ddata
	with pytest.raises(AssertionError):
		dplot.dotplot(ds, dc.iloc[:2])
	with pytest.raises(AssertionError):
		dplot.dotplot(ds, dc.rename(index={'r1': 'x'}))
	with pytest.raises(AssertionError):
		dplot.dotplot(ds, dc.rename(columns={'a': 'x'}))
