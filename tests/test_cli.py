"""Tests for the dictys command line as a whole: python -m dictys, bin/dictys and the generated parser.

Each subcommand's own behaviour is tested beside its function (test_chromatin.py, test_network.py, ...).
"""

import argparse
import inspect
import os
import runpy
import shutil
import subprocess
import sys
import typing

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Written by hand from the functions each module exposes; dictys.net, dictys.plot, dictys.traj and
# dictys.utils are excluded through _docstring2argparse_ignore_.
SUBCOMMANDS = {
	'network': ['reconstruct', 'indirect', 'normalize', 'tofile'],
	'chromatin': ['macs2', 'wellington', 'homer', 'binding', 'tssdist', 'linking', 'binlinking'],
	'preproc': ['selects_rna', 'selects_atac', 'qc_reads'],
	'dynamic': ['subsets_rna', 'subsets_atac'],
}
ALL = [(m, f) for m, fs in SUBCOMMANDS.items() for f in fs]


@pytest.fixture(scope='module')
def parser():
	from docstring2argparse import docstringparser, function_parser_union, function_parser_signature, function_parser_docstring_numpy
	return docstringparser('dictys', function_parser_union([function_parser_signature(), function_parser_docstring_numpy()]))


def subparser_actions(p):
	return [x for x in p._actions if isinstance(x, argparse._SubParsersAction)][0]


@pytest.fixture
def run_module(child_env):
	"""Run `python -m dictys args...` in a child process."""
	def run(*args):
		return subprocess.run([sys.executable, '-m', 'dictys'] + list(args), env=child_env(),
			capture_output=True, text=True, timeout=120)
	return run


# python -m dictys

def test_module_without_arguments_prints_help_and_fails(run_module):
	"""`python -m dictys` alone prints the usage to stderr and exits 1."""
	ans = run_module()
	assert ans.returncode == 1
	assert 'usage: dictys' in ans.stderr
	for module in SUBCOMMANDS:
		assert module in ans.stderr


def test_module_help_lists_modules(run_module):
	"""`python -m dictys -h` exits 0 and names the four modules with their docstrings."""
	ans = run_module('-h')
	assert ans.returncode == 0
	for module in SUBCOMMANDS:
		assert module in ans.stdout
	assert 'Chromatin accessibility analyses' in ans.stdout


def test_module_subcommand_help_as_in_ci(run_module):
	"""`python -m dictys preproc selects_rna -h`, the repository's install check, exits 0."""
	ans = run_module('preproc', 'selects_rna', '-h')
	assert ans.returncode == 0
	assert 'fi_reads' in ans.stdout and 'fo_reads' in ans.stdout


def test_module_unknown_subcommand_fails(run_module):
	"""An unknown module or function exits 2 with an argparse error."""
	assert run_module('nothing').returncode == 2
	assert run_module('network', 'nothing').returncode == 2


def test_main_module_in_process(monkeypatch, capsys):
	"""dictys.__main__ run as __main__ builds the parser and prints help for -h."""
	monkeypatch.setattr(sys, 'argv', ['dictys', 'network', '-h'])
	with pytest.raises(SystemExit) as e:
		runpy.run_module('dictys', run_name='__main__')
	assert e.value.code == 0
	out = capsys.readouterr().out
	for name in SUBCOMMANDS['network']:
		assert name in out


def test_main_module_import_does_not_run(monkeypatch, capsys):
	"""Importing dictys.__main__ under its own name runs nothing."""
	monkeypatch.setattr(sys, 'argv', ['dictys'])
	runpy.run_module('dictys.__main__', run_name='dictys.__main__')
	assert capsys.readouterr().err == ''


# bin/dictys

def test_bin_dictys_wrapper(tmp_path, child_env):
	"""bin/dictys runs `python3 -m dictys` with single-threaded BLAS variables."""
	if shutil.which('bash') is None:
		pytest.skip('needs bash')
	path = os.path.join(REPO, 'bin', 'dictys')
	# a python3 on PATH that reports the thread variables, then hands over to this interpreter
	fake = tmp_path / 'python3'
	fake.write_text('#!/bin/bash\necho "OMP=$OMP_NUM_THREADS MKL=$MKL_NUM_THREADS OPENBLAS=$OPENBLAS_NUM_THREADS NUMEXPR=$NUMEXPR_NUM_THREADS"\n'
		'exec "' + sys.executable + '" "$@"\n')
	fake.chmod(0o755)
	env = child_env(PATH=str(tmp_path) + os.pathsep + os.environ.get('PATH', ''))
	ans = subprocess.run(['bash', path, '-h'], env=env, capture_output=True, text=True, timeout=120)
	assert ans.returncode == 0
	assert 'MKL=1 OPENBLAS=1 NUMEXPR=1' in ans.stdout
	assert 'usage: dictys' in ans.stdout


# The generated parser

def test_parser_modules_and_functions(parser):
	"""The parser exposes exactly the expected modules and functions."""
	p, funcs = parser
	assert set(funcs) == {'dictys.{}.{}'.format(m, f) for m, f in ALL}
	modules = subparser_actions(p).choices
	assert set(modules) == set(SUBCOMMANDS)
	for module, fs in SUBCOMMANDS.items():
		assert set(subparser_actions(modules[module]).choices) == set(fs)


@pytest.mark.parametrize('module,func', ALL)
def test_parser_arguments_follow_signature(parser, module, func):
	"""Parameters without defaults are positionals in order; the others are --options with the function's default (bool ones are store_true flags, so their default must be False)."""
	import importlib
	p, _ = parser
	f = getattr(importlib.import_module('dictys.' + module), func)
	sig = inspect.signature(f)
	sub = subparser_actions(subparser_actions(p).choices[module]).choices[func]
	actions = [x for x in sub._actions if not isinstance(x, argparse._HelpAction)]
	positional = [x.dest for x in actions if not x.option_strings]
	options = {x.dest: x for x in actions if x.option_strings}
	assert positional == [x for x, y in sig.parameters.items() if y.default is inspect.Parameter.empty]
	assert set(options) == {x for x, y in sig.parameters.items() if y.default is not inspect.Parameter.empty}
	for name, action in options.items():
		param = sig.parameters[name]
		assert action.option_strings == ['--' + name]
		annot = param.annotation
		if typing.get_origin(annot) is typing.Union:
			annot = [x for x in typing.get_args(annot) if x is not type(None)][0]
		if annot is bool:
			# A store_true option always defaults to False, so the function's own default must be False too
			assert param.default is False
			assert isinstance(action, argparse._StoreTrueAction)
			assert action.default is False
		else:
			assert action.default == param.default
			assert action.type is annot
		assert action.help, 'option --{} of {} {} has no description'.format(name, module, func)
	for action in actions:
		if not action.option_strings:
			assert action.help, 'argument {} of {} {} has no description'.format(action.dest, module, func)


@pytest.mark.parametrize('module,func', ALL)
def test_subcommand_help(run_dictys, capsys, module, func):
	"""`dictys <module> <function> -h` exits 0 and lists every parameter."""
	import importlib
	f = getattr(importlib.import_module('dictys.' + module), func)
	with pytest.raises(SystemExit) as e:
		run_dictys(module, func, '-h')
	assert e.value.code == 0
	out = capsys.readouterr().out
	for name, param in inspect.signature(f).parameters.items():
		assert (name if param.default is inspect.Parameter.empty else '--' + name) in out


@pytest.mark.parametrize('module,func', ALL)
def test_subcommand_missing_arguments(run_dictys, module, func):
	"""A subcommand without its positional arguments exits 2."""
	with pytest.raises(SystemExit) as e:
		run_dictys(module, func)
	assert e.value.code == 2


def test_type_conversion_errors(run_dictys):
	"""An option value of the wrong type exits 2 before the function runs."""
	with pytest.raises(SystemExit) as e:
		run_dictys('chromatin', 'binlinking', 'a', 'b', 'not_an_int')
	assert e.value.code == 2
	with pytest.raises(SystemExit) as e:
		run_dictys('network', 'reconstruct', 'a', 'b', 'c', 'd', 'e', 'f', 'g', '--npc', '1.5')
	assert e.value.code == 2
