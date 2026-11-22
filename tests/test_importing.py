"""Tests for dictys.utils.importing."""

import runpy
import sys

import pytest

import dictys.utils.importing as importing
from dictys.utils.importing import dummy_class


# dummy_class

def test_dummy_attribute_access_returns_itself():
	"""Any attribute, at any depth, is the dummy itself."""
	d = dummy_class()
	assert d.anything is d
	assert d.nn.Module.forward is d


def test_dummy_call_returns_itself():
	"""Calling the dummy, with any arguments, returns the dummy."""
	d = dummy_class()
	assert d() is d
	assert d(1, 2, key='v') is d
	assert d.optim.Adam(lr=0.1) is d


def test_dummy_can_be_a_base_class():
	"""A class can name the dummy as a base: __mro_entries__ removes it, leaving object."""
	d = dummy_class()

	class child(d.nn.Module):
		value = 3

	assert child.__bases__ == (object,)
	assert child().value == 3


def test_dummy_mro_entries_returns_empty_tuple():
	"""__mro_entries__ is a callable giving an empty tuple."""
	d = dummy_class()
	assert d.__mro_entries__((d,)) == tuple()


# module-level imports

def test_installed_modules_are_the_real_ones():
	"""With torch, pyro and matplotlib installed the module exposes the real modules."""
	torch = pytest.importorskip('torch')
	pyro = pytest.importorskip('pyro')
	matplotlib = pytest.importorskip('matplotlib')
	import matplotlib.pyplot
	import matplotlib.figure
	assert importing.torch is torch
	assert importing.pyro is pyro
	assert importing.matplotlib is matplotlib
	assert importing.matplotlib.pyplot is matplotlib.pyplot
	assert importing.matplotlib.figure is matplotlib.figure


@pytest.mark.parametrize('blocked,name', [
	(['torch'], 'torch'),
	(['pyro'], 'pyro'),
	(['matplotlib', 'matplotlib.pyplot', 'matplotlib.figure'], 'matplotlib'),
])
def test_missing_module_becomes_dummy(monkeypatch, blocked, name):
	"""A module that cannot be imported is replaced by a dummy_class instance."""
	for x in blocked:
		monkeypatch.setitem(sys.modules, x, None)
	g = runpy.run_path(importing.__file__)
	assert isinstance(g[name], g['dummy_class'])
	# The other names are unaffected
	for other in {'torch', 'pyro', 'matplotlib'} - {name}:
		if importing.__dict__[other].__class__.__name__ != 'dummy_class':
			assert not isinstance(g[other], g['dummy_class'])


def test_missing_matplotlib_submodule_becomes_dummy(monkeypatch):
	"""matplotlib is a dummy when only matplotlib.pyplot cannot be imported."""
	pytest.importorskip('matplotlib')
	monkeypatch.setitem(sys.modules, 'matplotlib.pyplot', None)
	g = runpy.run_path(importing.__file__)
	assert isinstance(g['matplotlib'], g['dummy_class'])
