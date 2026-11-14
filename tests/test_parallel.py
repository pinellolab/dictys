"""Tests for dictys.utils.parallel."""

import logging
import os
from multiprocessing import cpu_count

import pytest

from dictys.utils.parallel import autocount, num_threads, recover_num_threads, set_num_threads

KEYS = ['OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'OMP_MAX_THREADS',
	'MKL_MAX_THREADS', 'NUMEXPR_MAX_THREADS', 'OPENBLAS_MAX_THREADS', 'VECLIB_MAXIMUM_THREADS']


@pytest.fixture
def clean_env(monkeypatch):
	"""Remove every thread variable so each test starts from a known environment; monkeypatch restores them."""
	for x in KEYS:
		monkeypatch.delenv(x, raising=False)
	return monkeypatch


# autocount

def test_autocount_zero_means_all_cpus():
	"""n=0 returns the CPU count."""
	assert autocount(0) == cpu_count()


@pytest.mark.parametrize('n,expected', [(1, 1), (3, 3), (2.7, 2), (64, 64)])
def test_autocount_integer_part_of_n(n, expected):
	"""n>=1 returns int(n), independent of the CPU count."""
	assert autocount(n) == expected


def test_autocount_fraction_of_cpus(monkeypatch):
	"""0<n<1 returns max(1, int(n*cpu_count))."""
	import multiprocessing
	monkeypatch.setattr(multiprocessing, 'cpu_count', lambda: 10)
	assert autocount(0.5) == 5
	assert autocount(0.25) == 2
	assert autocount(0.01) == 1
	assert autocount(0) == 10


def test_autocount_rejects_negative():
	"""Negative n is rejected."""
	with pytest.raises(AssertionError):
		autocount(-1)


def test_autocount_returns_int():
	"""The result is an int for every accepted input."""
	assert all(isinstance(autocount(x), int) for x in [0, 0.5, 1, 2.0])


# set_num_threads / recover_num_threads

def test_set_num_threads_sets_every_key_and_returns_previous(clean_env):
	"""Every thread variable becomes str(n); the return records previous values, None when unset."""
	clean_env.setenv('OMP_NUM_THREADS', '7')
	old = set_num_threads(3)
	assert set(old) == set(KEYS)
	assert old['OMP_NUM_THREADS'] == '7'
	assert all(old[x] is None for x in KEYS if x != 'OMP_NUM_THREADS')
	assert all(os.environ[x] == '3' for x in KEYS)


def test_recover_num_threads_restores_previous_state(clean_env):
	"""Recovering restores set values and removes variables that were unset."""
	clean_env.setenv('MKL_NUM_THREADS', '5')
	old = set_num_threads(2)
	recover_num_threads(old)
	assert os.environ['MKL_NUM_THREADS'] == '5'
	assert all(x not in os.environ for x in KEYS if x != 'MKL_NUM_THREADS')


def test_recover_num_threads_warns_on_missing_variable(clean_env, caplog):
	"""A variable removed in between is reported with a warning and left unset."""
	old = set_num_threads(2)
	del os.environ['OPENBLAS_NUM_THREADS']
	with caplog.at_level(logging.WARNING):
		recover_num_threads(old)
	assert 'OPENBLAS_NUM_THREADS' in caplog.text
	assert all(x not in os.environ for x in KEYS)


def test_recover_num_threads_missing_variable_with_old_value(clean_env, caplog):
	"""A variable that had a value but was removed in between is not recreated."""
	clean_env.setenv('OMP_NUM_THREADS', '9')
	old = set_num_threads(2)
	del os.environ['OMP_NUM_THREADS']
	with caplog.at_level(logging.WARNING):
		recover_num_threads(old)
	assert 'OMP_NUM_THREADS' not in os.environ
	assert 'OMP_NUM_THREADS' in caplog.text


# num_threads

def test_num_threads_sets_environment_and_backends(clean_env):
	"""Inside the context the variables, joblib's backend and the thread pools use n threads; all are restored after."""
	pytest.importorskip('joblib')
	threadpoolctl = pytest.importorskip('threadpoolctl')
	from joblib.parallel import get_active_backend
	import numpy as np
	np.ones(2) @ np.ones(2)
	before = [x['num_threads'] for x in threadpoolctl.threadpool_info()]
	with num_threads(2):
		assert all(os.environ[x] == '2' for x in KEYS)
		backend, n_jobs = get_active_backend()
		assert n_jobs == 2
		assert type(backend).__name__ == 'ThreadingBackend'
		assert all(x['num_threads'] <= 2 for x in threadpoolctl.threadpool_info())
	assert all(x not in os.environ for x in KEYS)
	assert [x['num_threads'] for x in threadpoolctl.threadpool_info()] == before


def test_num_threads_restores_on_exception(clean_env):
	"""The environment is restored when the body raises."""
	pytest.importorskip('joblib')
	pytest.importorskip('threadpoolctl')
	clean_env.setenv('OMP_NUM_THREADS', '6')
	with pytest.raises(KeyError):
		with num_threads(1):
			assert os.environ['OMP_NUM_THREADS'] == '1'
			raise KeyError('x')
	assert os.environ['OMP_NUM_THREADS'] == '6'
	assert 'MKL_NUM_THREADS' not in os.environ


def test_num_threads_fraction_uses_autocount(clean_env, monkeypatch):
	"""A fractional n is converted with autocount before use."""
	pytest.importorskip('joblib')
	pytest.importorskip('threadpoolctl')
	import multiprocessing
	monkeypatch.setattr(multiprocessing, 'cpu_count', lambda: 8)
	with num_threads(0.5):
		assert os.environ['OMP_NUM_THREADS'] == '4'


def test_num_threads_rejects_negative(clean_env):
	"""Negative n is rejected before anything changes."""
	pytest.importorskip('joblib')
	pytest.importorskip('threadpoolctl')
	with pytest.raises(AssertionError):
		with num_threads(-2):
			pass
	assert all(x not in os.environ for x in KEYS)
