"""Tests for dictys.utils.numpy."""

import numpy as np
import pytest

from dictys.utils import numpy as dnp
from dictys.utils.numpy import dtype_min, groupby, median


# groupby

def test_groupby_indices_per_value():
	"""Each distinct value maps to the positions where it occurs, in increasing order."""
	ans = groupby(np.array([3, 1, 3, 2, 1, 3]))
	assert set(ans) == {1, 2, 3}
	np.testing.assert_array_equal(ans[3], [0, 2, 5])
	np.testing.assert_array_equal(ans[1], [1, 4])
	np.testing.assert_array_equal(ans[2], [3])
	assert all(isinstance(x, np.ndarray) for x in ans.values())


def test_groupby_accepts_lists_of_str():
	"""Works on a plain list of str, as chromatin.binding calls it."""
	ans = groupby(['a', 'b', 'a'])
	assert set(ans) == {'a', 'b'}
	np.testing.assert_array_equal(ans['a'], [0, 2])
	np.testing.assert_array_equal(ans['b'], [1])


def test_groupby_partition_covers_every_index(rng):
	"""The groups partition range(n) and every index in a group holds that group's value."""
	v = rng.integers(0, 5, size=50)
	ans = groupby(v)
	allidx = np.sort(np.concatenate(list(ans.values())))
	np.testing.assert_array_equal(allidx, np.arange(50))
	for key, idx in ans.items():
		assert (v[idx] == key).all()


def test_groupby_empty():
	"""An empty array gives an empty dict."""
	assert groupby(np.array([])) == {}


# median

def test_median_unweighted_uses_numpy_median():
	"""Without weights the result is numpy.median, keyword arguments included."""
	v = np.array([[1., 5., 3.], [2., 8., 4.]])
	assert median(np.array([4., 1., 3., 2.])) == 2.5
	np.testing.assert_array_equal(median(v, axis=1), [3., 4.])
	np.testing.assert_array_equal(median(v, axis=0), [1.5, 6.5, 3.5])


@pytest.mark.parametrize('v,w,expected', [
	# cumulative weights 1,2,3; half of total 1.5 is first reached at the 2nd sorted value
	([3., 1., 2.], [1., 1., 1.], 2.),
	# cumulative weights 1,2,3,4; half 2 is reached exactly at the 2nd sorted value (lower median)
	([4., 3., 2., 1.], [1., 1., 1., 1.], 2.),
	# one heavy value dominates
	([1., 2., 3., 4.], [1., 0., 0., 10.], 4.),
	([10., 20., 30.], [5., 1., 1.], 10.),
	# zero weight on the middle values is skipped: cumulative 2,2,2,4, half 2 reached at sorted value 1
	([1., 2., 3., 4.], [2., 0., 0., 2.], 1.),
])
def test_median_weighted_by_hand(v, w, expected):
	"""Weighted median is the first sorted value whose cumulative weight reaches half the total."""
	assert median(np.array(v), np.array(w)) == expected


def test_median_weighted_equal_weights_matches_lower_median(rng):
	"""With equal weights the weighted median is the lower median of the values."""
	v = rng.normal(size=11)
	assert median(v, np.ones(11)) == np.sort(v)[5]
	v = rng.normal(size=10)
	assert median(v, np.full(10, 2.5)) == np.sort(v)[4]


def test_median_weighted_integer_weights_match_repetition(rng):
	"""Integer weights act like repeating each value that many times (lower median)."""
	v = rng.permutation(np.arange(7.))
	w = rng.integers(1, 4, size=7)
	rep = np.sort(np.repeat(v, w))
	lower = rep[(len(rep) - 1) // 2]
	assert median(v, w.astype(float)) == lower


def test_median_weighted_rejects_keywords():
	"""Keyword arguments are only accepted without weights."""
	with pytest.raises(TypeError):
		median(np.array([1., 2.]), np.array([1., 1.]), axis=0)


@pytest.mark.parametrize('v,w', [
	([1., np.nan], [1., 1.]),
	([1., 2.], [1., np.inf]),
	([1., 2.], [1., -1.]),
	([1., 2., 3.], [1., 1.]),
])
def test_median_weighted_rejects_invalid_input(v, w):
	"""Non-finite values or weights, negative weights and mismatched shapes are rejected."""
	with pytest.raises(AssertionError):
		median(np.array(v), np.array(w))


def test_median_weighted_rejects_2d():
	"""Weighted median needs 1-D input."""
	with pytest.raises(AssertionError):
		median(np.ones((2, 2)), np.ones((2, 2)))


# dtype_min

@pytest.mark.parametrize('values,dtype,expected', [
	([0, 200], 'i8', 'i2'),
	([-129, 5], 'i8', 'i2'),
	([0, 32767], 'i8', 'i2'),
	([0, 32768], 'i8', 'i4'),
	([-2**31, 0], 'i8', 'i4'),
	([0, 2**31], 'i8', 'i8'),
	([0, 256], 'u8', 'u2'),
	([0, 65535], 'u4', 'u2'),
	([0, 65536], 'u8', 'u4'),
	([2**32], 'u8', 'u8'),
	([1000], 'i4', 'i2'),
])
def test_dtype_min_smallest_fitting_type(values, dtype, expected):
	"""The result is the smallest integer type of the same signedness, at least 2 bytes here, that holds the range."""
	v = np.array(values, dtype=dtype)
	ans = dtype_min(v)
	assert ans == np.dtype(expected)
	info = np.iinfo(ans)
	assert info.min <= v.min() and info.max >= v.max()


def test_dtype_min_keeps_type_when_nothing_smaller_fits():
	"""A range that needs the input's own size returns the input's size."""
	assert dtype_min(np.array([-2**40, 2**40], dtype='i8')) == np.dtype('i8')
	assert dtype_min(np.array([0, 2**17], dtype='i4')) == np.dtype('i4')


@pytest.mark.parametrize('values', [[1.5, 2.], [True, False]])
def test_dtype_min_rejects_non_integer(values):
	"""Float and bool arrays raise TypeError."""
	with pytest.raises(TypeError):
		dtype_min(np.array(values))


def test_typing_aliases():
	"""NDArray and ArrayLike are numpy.typing's aliases when numpy.typing is available."""
	import numpy.typing as npt
	assert dnp.NDArray is npt.NDArray
	assert dnp.ArrayLike is npt.ArrayLike
