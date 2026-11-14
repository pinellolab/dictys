"""
Shared fixtures for the dictys unit tests.

The tests import dictys and docstring2argparse from an installed copy or a
PYTHONPATH when one exists, and from this checkout's ``src`` folder otherwise.
Command-line calls run in this process
through the same argparse parser ``python -m dictys`` builds, and are appended
to ``$TESTS_PR_CALL_LOG`` (one JSON line per call) when that variable is set.
"""

import importlib.util
import io
import json
import os
import runpy
import shutil
import subprocess
import sys

import pytest

os.environ.setdefault('MPLBACKEND', 'Agg')

TESTS = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(TESTS)
SRC = os.path.join(REPO, 'src')
# Test modules import helpers from each other; this keeps that working when pytest
# runs with --import-mode=importlib, which does not put this folder on sys.path.
if TESTS not in sys.path:
	sys.path.insert(0, TESTS)
# Use this checkout's src only when dictys is not importable otherwise (an install,
# or a PYTHONPATH naming a worktree), so an installed copy is the one tested.
if importlib.util.find_spec('dictys') is None and os.path.isdir(os.path.join(SRC, 'dictys')):
	sys.path.insert(0, SRC)


def _log_call(program, argv):
	path = os.environ.get('TESTS_PR_CALL_LOG')
	if not path:
		return
	with open(path, 'a') as f:
		f.write(json.dumps({'program': program, 'argv': [str(x) for x in argv]}) + '\n')


def scripts_dir():
	"""Folder of the shell and Python scripts shipped inside the dictys package."""
	import dictys
	return os.path.join(os.path.dirname(os.path.abspath(dictys.__file__)), 'scripts')


def subprocess_env(**extra):
	"""Environment for a child process that imports the same dictys as the tests do."""
	import dictys
	root = os.path.dirname(os.path.dirname(os.path.abspath(dictys.__file__)))
	env = dict(os.environ)
	env['PYTHONPATH'] = root + (os.pathsep + env['PYTHONPATH'] if env.get('PYTHONPATH') else '')
	env.update({x: str(y) for x, y in extra.items()})
	return env


def pytest_configure(config):
	# docutils, which docstring2argparse uses, warns once per option it defines
	config.addinivalue_line('filterwarnings', 'ignore:The frontend.Option:DeprecationWarning')
	config.addinivalue_line('filterwarnings', 'ignore:The frontend.OptionParser:DeprecationWarning')


@pytest.fixture
def rng():
	"""A seeded numpy random generator."""
	import numpy as np
	return np.random.default_rng(12345)


@pytest.fixture
def run_dictys(monkeypatch):
	"""
	Run the ``dictys`` command line in this process, as ``python -m dictys``
	does: ``run_dictys('preproc', 'selects_rna', a, b, c)``. argparse errors
	raise SystemExit.
	"""
	def run(*args):
		from docstring2argparse import docstringrunner
		argv = [str(x) for x in args]
		_log_call('dictys', argv)
		monkeypatch.setattr(sys, 'argv', ['dictys'] + argv)
		return docstringrunner('dictys')
	return run


@pytest.fixture
def run_helper(monkeypatch):
	"""
	Run a Python helper script under ``dictys/scripts`` in this process, as
	``dictys_helper <name> args...`` would: ``run_helper('makefile_update.py',
	path, json)``. ``stdin`` is a str fed to ``sys.stdin``; ``cwd`` changes
	the working directory for the call. Returns the script's globals.
	"""
	def run(name, *args, stdin=None, cwd=None, folder='helper'):
		path = os.path.join(scripts_dir(), folder, name) if folder else os.path.join(scripts_dir(), name)
		argv = [str(x) for x in args]
		_log_call('dictys_helper' if folder == 'helper' else name, ([name] if folder == 'helper' else []) + argv)
		monkeypatch.setattr(sys, 'argv', [path] + argv)
		if stdin is not None:
			monkeypatch.setattr(sys, 'stdin', io.StringIO(stdin))
		if cwd is not None:
			monkeypatch.chdir(cwd)
		return runpy.run_path(path, run_name='__main__')
	return run


def require_programs(*names):
	"""Skip the calling test unless every named program is on PATH."""
	missing = [x for x in names if shutil.which(x) is None]
	if missing:
		pytest.skip('needs program(s) not on PATH: ' + ', '.join(missing))


@pytest.fixture
def run_shell():
	"""
	Run a bash script under ``dictys/scripts`` (``folder`` is a subfolder such
	as 'helper', or '' for the top level) and return the CompletedProcess.
	Skips when bash is not available.
	"""
	def run(name, *args, cwd=None, folder='helper', env=None, timeout=120, input=None):
		require_programs('bash')
		path = os.path.join(scripts_dir(), folder, name) if folder else os.path.join(scripts_dir(), name)
		argv = [str(x) for x in args]
		_log_call('dictys_helper' if folder == 'helper' else name, ([name] if folder == 'helper' else []) + argv)
		return subprocess.run(['bash', path] + argv, cwd=cwd, env=subprocess_env(**(env or {})),
			capture_output=True, text=True, timeout=timeout, input=input)
	return run


@pytest.fixture
def dictys_scripts():
	"""Path of the ``scripts`` folder inside the dictys package."""
	return scripts_dir()


@pytest.fixture
def require():
	"""``require('samtools', 'bedtools')`` skips the test unless all are on PATH."""
	return require_programs


@pytest.fixture
def child_env():
	"""``child_env(KEY=value)``: environment for a child process importing this checkout."""
	return subprocess_env
