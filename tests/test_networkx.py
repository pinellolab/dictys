"""Tests for dictys.utils.networkx."""

import numpy as np

from dictys.utils.networkx import random_state


def test_random_state_turns_int_seed_into_generator():
	"""The decorator converts an int seed argument into numpy's RandomState seeded with it."""
	@random_state('seed')
	def draw(seed=None):
		return seed.random(3)

	np.testing.assert_array_equal(draw(seed=7), np.random.RandomState(7).random(3))
	np.testing.assert_array_equal(draw(seed=7), draw(seed=7))
	assert not np.array_equal(draw(seed=7), draw(seed=8))


def test_random_state_by_position():
	"""The seed can be named by argument position."""
	@random_state(1)
	def draw(n, seed):
		return seed.randint(0, 1000, size=n)

	np.testing.assert_array_equal(draw(4, 3), np.random.RandomState(3).randint(0, 1000, size=4))


def test_random_state_passes_generator_through():
	"""A RandomState passed in is used as is."""
	@random_state('seed')
	def same(seed=None):
		return seed

	r = np.random.RandomState(1)
	assert same(seed=r) is r
