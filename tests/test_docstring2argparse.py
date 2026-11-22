"""Tests for docstring2argparse, which builds the dictys command line from signatures and numpy docstrings."""

import argparse
import sys
import textwrap
import typing
from inspect import _empty

import pytest

import docstring2argparse as d2a

PKG = 'd2a_fixture_pkg'

INIT = '''
"""Fixture package docstring"""
__all__ = ['tools', 'hidden', 'extra']
_docstring2argparse_ignore_ = ['hidden']
from . import *
'''

TOOLS = '''
"""Tools module doc"""
from typing import Optional

def add(a: int, b: int, scale: float = 1.0, label: Optional[str] = None, flag: bool = False) -> None:
	"""
	Add two numbers.

	Writes the scaled sum to the record.

	Parameters
	----------
	a:
		First number
	b:
		Second number
	scale:
		Multiplier
	label:
		Optional label
	flag:
		Switch
	"""
	RECORD.append((a, b, scale, label, flag))
	return (a + b) * scale

def _private(x: int) -> None:
	"""Private function"""

RECORD = []
'''

HIDDEN = '''
"""Hidden module doc"""
def secret(x: int) -> None:
	"""Secret function"""
'''

EXTRA = '''
"""Extra module doc"""
from .tools import add as imported_add

def skip_me(x: int) -> None:
	"""
	Skipped by func_ignore.

	Parameters
	----------
	x:
		Value
	"""

def double(x: float) -> None:
	"""
	Double a number.

	Parameters
	----------
	x:
		Value
	"""
	return 2 * x
'''


@pytest.fixture
def fixture_pkg(tmp_path, monkeypatch):
	"""A small package on sys.path: tools.add, extra.double, extra.skip_me and an ignored module hidden."""
	d = tmp_path / PKG
	d.mkdir()
	for name, text in [('__init__.py', INIT), ('tools.py', TOOLS), ('hidden.py', HIDDEN), ('extra.py', EXTRA)]:
		(d / name).write_text(textwrap.dedent(text))
	monkeypatch.syspath_prepend(str(tmp_path))
	for x in [x for x in sys.modules if x == PKG or x.startswith(PKG + '.')]:
		monkeypatch.delitem(sys.modules, x)
	yield PKG
	for x in [x for x in sys.modules if x == PKG or x.startswith(PKG + '.')]:
		del sys.modules[x]


def full_parser():
	return d2a.function_parser_union([d2a.function_parser_signature(), d2a.function_parser_docstring_numpy()])


def sample(a: int, b: str = 'x', c: typing.Optional[float] = None) -> None:
	"""
	Short line.

	Longer description
	over two lines.

	Parameters
	----------
	a:
		The a value
	b:
		The b value
	c:
		The c value
	"""


def typed_doc(a: int, b: str = '') -> None:
	"""
	Short.

	Parameters
	----------
	a:		int
		The a value
	b:
		The b value
	"""


def no_return(a: int):
	"""
	Short.

	Parameters
	----------
	a:
		A
	"""


# function_parser_signature

def test_signature_parser_reads_names_types_defaults():
	"""Each parameter gives (name, annotation, None, (has default, default)); the return annotation is kept."""
	ans = d2a.function_parser_signature().parse(sample)
	assert ans[0] is None and ans[1] is None
	assert ans[2] == [('a', int, None, (False, None)), ('b', str, None, (True, 'x')),
		('c', typing.Optional[float], None, (True, None))]
	assert ans[3] is None


def test_signature_parser_call_equals_parse():
	"""Calling the parser is the same as parse."""
	p = d2a.function_parser_signature()
	assert p(sample) == p.parse(sample)


def test_signature_parser_requires_return_annotation():
	"""A function without a return annotation is rejected."""
	with pytest.raises(TypeError):
		d2a.function_parser_signature().parse(no_return)


# function_parser_docstring_numpy

def test_docstring_parser_descriptions_and_parameters():
	"""The first paragraph is the short description, the remaining text the long one; parameters carry their descriptions."""
	ans = d2a.function_parser_docstring_numpy().parse(sample)
	assert ans[0] == 'Short line.'
	assert ans[1] == 'Longer description\nover two lines.'
	assert ans[2] == [('a', None, 'The a value', None), ('b', None, 'The b value', None), ('c', None, 'The c value', None)]
	assert ans[3] is None


def test_docstring_parser_type_after_colon():
	"""Text after 'name:' in a parameter term is kept as the type text; None when empty."""
	ans = d2a.function_parser_docstring_numpy().parse(typed_doc)
	assert ans[2] == [('a', 'int', 'The a value', None), ('b', None, 'The b value', None)]


def test_docstring_parser_short_description_only():
	"""A docstring of one paragraph has a short description only (checked before the format check)."""
	def f() -> None:
		"""Only a line."""
	ans = d2a.function_parser_docstring_numpy()._parse(f)
	assert ans == ['Only a line.', None, None, None]


def test_docstring_parser_drops_yields_section():
	"""Yields and Receives sections are removed from the long description."""
	def f(x: int) -> None:
		"""
		Short.

		Kept text.

		Parameters
		----------
		x:
			Value

		Yields
		------
		int
			Numbers
		"""
	ans = d2a.function_parser_docstring_numpy().parse(f)
	assert ans[1] == 'Kept text.'


def test_docstring_parser_results_section_not_implemented():
	"""A Results section is not supported."""
	def f() -> None:
		"""
		Short.

		Results
		-------
		x:
			Value
		"""
	with pytest.raises(NotImplementedError):
		d2a.function_parser_docstring_numpy().parse(f)


def test_docstring_parser_parameter_needs_colon():
	"""A parameter term without 'name:' is rejected."""
	def f(a: int) -> None:
		"""
		Short.

		Parameters
		----------
		a
			Value
		"""
	with pytest.raises(AssertionError):
		d2a.function_parser_docstring_numpy().parse(f)


# function_parser_union

def test_union_combines_signature_and_docstring():
	"""Signature types and defaults are combined with docstring descriptions."""
	ans = full_parser().parse(sample)
	assert ans[0] == 'Short line.'
	assert ans[1] == 'Longer description\nover two lines.'
	assert ans[2] == [('a', int, 'The a value', (False, None)), ('b', str, 'The b value', (True, 'x')),
		('c', typing.Optional[float], 'The c value', (True, None))]


def test_union_rejects_conflicting_types():
	"""A docstring type text conflicts with a signature annotation of the same parameter."""
	with pytest.raises(AssertionError):
		full_parser().parse(typed_doc)


def test_union_rejects_non_parsers():
	"""Only function_parser_base instances can be combined."""
	with pytest.raises(AssertionError):
		d2a.function_parser_union([d2a.function_parser_signature(), object()])


@pytest.mark.parametrize('values,expected', [
	([None, None], None),
	([None, 3], 3),
	([3, 3], 3),
	([(1, None), (None, 2)], (1, 2)),
	([[None, 'a'], ['b', None]], ['b', 'a']),
	([(1, (None, 5)), (None, (4, None))], (1, (4, 5))),
])
def test_union_merges_values(values, expected):
	"""None is filled from the other value, equal values pass through, sequences merge element-wise."""
	assert d2a.function_parser_union.union(values) == expected


def test_union_different_lengths_raise():
	"""Sequences of different lengths cannot be merged."""
	with pytest.raises(ValueError):
		d2a.function_parser_union.union([[1, 2], [1, 2, 3]])


def test_union_conflict_raises():
	"""Two different non-None scalars of the same type conflict."""
	with pytest.raises(AssertionError):
		d2a.function_parser_union.union([1, 2])


def test_union_different_types_raise():
	"""Values of different types conflict."""
	with pytest.raises(AssertionError):
		d2a.function_parser_union.union([1, 'a'])


def test_parser_base_parse_not_implemented():
	"""The base class has no parser."""
	with pytest.raises(NotImplementedError):
		d2a.function_parser_base().parse(sample)


# function_parser_base.check

def test_check_accepts_valid_and_returns_input():
	"""A well-formed result is returned unchanged."""
	ans = ['s', 'l', [('a', int, 'd', (True, 1)), None], ('r', int, 'doc')]
	assert d2a.function_parser_base.check(ans, 'f') is ans


@pytest.mark.parametrize('ans', [
	['s', 'l', []],
	[1, None, None, None],
	[None, 1, None, None],
	[None, None, 'x', None],
	[None, None, [('a', int, 'd', (False, 3))], None],
	[None, None, [('a', int, 'd')], None],
	[None, None, [(1, int, 'd', (True, 1))], None],
])
def test_check_rejects_malformed(ans):
	"""Wrong lengths or types, or a default for a required parameter, are rejected."""
	with pytest.raises(AssertionError):
		d2a.function_parser_base.check(ans, 'f')


def test_check_rejects_missing_return_annotation():
	"""An empty return annotation raises TypeError naming the function."""
	with pytest.raises(TypeError, match='myfunc'):
		d2a.function_parser_base.check([None, None, [], _empty], 'myfunc')


# totype

def test_totype():
	"""A type passes through; Optional[X] becomes X; other annotations are rejected."""
	assert d2a.totype(int) is int
	assert d2a.totype(str) is str
	assert d2a.totype(typing.Optional[float]) is float
	assert d2a.totype(typing.Union[None, int]) is int
	with pytest.raises(TypeError):
		d2a.totype(typing.Union[int, str])
	with pytest.raises(TypeError):
		d2a.totype(typing.Union[int, str, None])


# get_functions

def test_get_functions_finds_public_functions(fixture_pkg):
	"""Public functions defined in submodules are found; ignored modules, private names and imports are skipped."""
	f, m = d2a.get_functions(fixture_pkg, full_parser())
	assert set(f) == {PKG + '.tools.add', PKG + '.extra.double', PKG + '.extra.skip_me'}
	assert f[PKG + '.tools.add'][0] == 'Add two numbers.'
	assert f[PKG + '.tools.add'][1] == 'Writes the scaled sum to the record.'
	assert [x[0] for x in f[PKG + '.tools.add'][2]] == ['a', 'b', 'scale', 'label', 'flag']
	assert m == {PKG: 'Fixture package docstring', PKG + '.tools': 'Tools module doc', PKG + '.extra': 'Extra module doc'}


def test_get_functions_custom_ignore(fixture_pkg):
	"""func_ignore and varname_ignore control what is skipped."""
	f, _ = d2a.get_functions(fixture_pkg, full_parser(), func_ignore=lambda name, obj: name.startswith('_') or name == 'skip_me')
	assert set(f) == {PKG + '.tools.add', PKG + '.extra.double'}
	f, m = d2a.get_functions(fixture_pkg, d2a.function_parser_signature(), varname_ignore=None)
	assert PKG + '.hidden.secret' in f
	assert m[PKG + '.hidden'] == 'Hidden module doc'


# docstringparser and run_args

def test_docstringparser_builds_nested_subcommands(fixture_pkg):
	"""Subcommands follow module and function names; options follow the signature."""
	p, f = d2a.docstringparser(fixture_pkg, full_parser())
	assert isinstance(p, argparse.ArgumentParser)
	args = p.parse_args(['tools', 'add', '2', '3'])
	assert args.subcommand_1 == 'tools' and args.subcommand_2 == 'add'
	assert (args.a, args.b, args.scale, args.label, args.flag) == (2, 3, 1.0, None, False)
	args = p.parse_args(['tools', 'add', '2', '3', '--scale', '0.5', '--label', 'L', '--flag'])
	assert (args.scale, args.label, args.flag) == (0.5, 'L', True)
	assert isinstance(args.scale, float)
	args = p.parse_args(['extra', 'double', '1.5'])
	assert args.x == 1.5


def test_docstringparser_rejects_bad_input(fixture_pkg, capsys):
	"""Missing positionals, wrong types and unknown subcommands exit with code 2."""
	p, _ = d2a.docstringparser(fixture_pkg, full_parser())
	for argv in [['tools', 'add', '2'], ['tools', 'add', 'x', '3'], ['tools', 'nothing'], [], ['tools']]:
		with pytest.raises(SystemExit) as e:
			p.parse_args(argv)
		assert e.value.code == 2


def test_docstringparser_help_shows_docstrings(fixture_pkg, capsys):
	"""Help shows the package doc, the function's short description and parameter descriptions with defaults."""
	p, _ = d2a.docstringparser(fixture_pkg, full_parser())
	with pytest.raises(SystemExit):
		p.parse_args(['-h'])
	out = capsys.readouterr().out
	assert 'Fixture package docstring' in out
	with pytest.raises(SystemExit):
		p.parse_args(['tools', 'add', '-h'])
	out = capsys.readouterr().out
	assert 'Add two numbers.' in out and 'Writes the scaled sum' in out
	assert 'Multiplier' in out and 'Second number' in out
	assert '1.0' in out


def test_run_args_calls_function(fixture_pkg):
	"""run_args passes positionals and options to the function and returns its value."""
	p, f = d2a.docstringparser(fixture_pkg, full_parser())
	import importlib
	tools = importlib.import_module(PKG + '.tools')
	args = p.parse_args(['tools', 'add', '2', '3', '--scale', '2', '--flag'])
	assert d2a.run_args(fixture_pkg, f, args) == 10.0
	assert tools.RECORD[-1] == (2, 3, 2.0, None, True)
	args = p.parse_args(['extra', 'double', '4'])
	assert d2a.run_args(fixture_pkg, f, args) == 8.0


def test_docstringrunner_runs_from_argv(fixture_pkg, monkeypatch):
	"""docstringrunner parses sys.argv and calls the function."""
	import importlib
	tools = importlib.import_module(PKG + '.tools')
	monkeypatch.setattr(sys, 'argv', ['prog', 'tools', 'add', '1', '4', '--label', 'z'])
	d2a.docstringrunner(fixture_pkg)
	assert tools.RECORD[-1] == (1, 4, 1.0, 'z', False)


def test_docstringrunner_without_arguments_prints_help(fixture_pkg, monkeypatch, capsys):
	"""Without arguments the help goes to stderr and the exit code is 1."""
	monkeypatch.setattr(sys, 'argv', ['prog'])
	with pytest.raises(SystemExit) as e:
		d2a.docstringrunner(fixture_pkg)
	assert e.value.code == 1
	assert 'tools' in capsys.readouterr().err
