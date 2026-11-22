"""Tests for dictys.utils.file."""

import os

import pytest

from dictys.utils.file import read_columns, read_index, read_txt


def _write(path, text):
	with open(path, 'w') as f:
		f.write(text)
	return str(path)


# read_txt

def test_read_txt_strips_and_drops_empty_lines(tmp_path):
	"""Lines are stripped of surrounding white space and empty lines are removed, order kept."""
	fi = _write(tmp_path / 'a.txt', '  gene1 \n\n\tgene2\n   \ngene3\n')
	assert read_txt(fi) == ['gene1', 'gene2', 'gene3']


def test_read_txt_keeps_duplicates_by_default(tmp_path):
	"""Without unique=True duplicate lines are all returned."""
	fi = _write(tmp_path / 'a.txt', 'x\ny\nx\n')
	assert read_txt(fi) == ['x', 'y', 'x']


def test_read_txt_unique_accepts_distinct_lines(tmp_path):
	"""unique=True returns the lines when none repeats."""
	fi = _write(tmp_path / 'a.txt', 'x\ny\nz')
	assert read_txt(fi, unique=True) == ['x', 'y', 'z']


def test_read_txt_unique_rejects_duplicates(tmp_path):
	"""unique=True raises when a line repeats after stripping."""
	fi = _write(tmp_path / 'a.txt', 'x\ny\n x \n')
	with pytest.raises(AssertionError):
		read_txt(fi, unique=True)


def test_read_txt_empty_file(tmp_path):
	"""A file of blank lines gives an empty list."""
	fi = _write(tmp_path / 'a.txt', '\n \n')
	assert read_txt(fi) == []


def test_read_txt_missing_file(tmp_path):
	"""A missing file raises FileNotFoundError."""
	with pytest.raises(FileNotFoundError):
		read_txt(os.path.join(str(tmp_path), 'missing.txt'))


# read_columns

def test_read_columns_excludes_index_column(tmp_path):
	"""Column names after the index column, as str, in file order."""
	fi = _write(tmp_path / 't.tsv', 'gene\tc1\tc2\tc3\ng1\t1\t2\t3\ng2\t4\t5\t6\n')
	assert read_columns(fi) == ['c1', 'c2', 'c3']


def test_read_columns_numeric_names_become_str(tmp_path):
	"""Numeric column names are returned as str."""
	fi = _write(tmp_path / 't.tsv', 'x\t1\t2\na\t0\t0\n')
	assert read_columns(fi) == ['1', '2']


def test_read_columns_custom_separator_and_no_index(tmp_path):
	"""sep and index_col are passed to pandas.read_csv."""
	fi = _write(tmp_path / 't.csv', 'a,b,c\n1,2,3\n')
	assert read_columns(fi, sep=',', index_col=None) == ['a', 'b', 'c']


def test_read_columns_passes_extra_keywords(tmp_path):
	"""Extra keyword arguments reach pandas.read_csv (here skiprows)."""
	fi = _write(tmp_path / 't.tsv', '# comment line\nid\tu\tv\nr\t1\t2\n')
	assert read_columns(fi, skiprows=1) == ['u', 'v']


def test_read_columns_unique(tmp_path):
	"""unique=True returns the names when they are distinct."""
	fi = _write(tmp_path / 't.tsv', 'id\tu\tv\nr\t1\t2\n')
	assert read_columns(fi, unique=True) == ['u', 'v']


def test_read_columns_without_header(tmp_path):
	"""header=None numbers the columns; the index column 0 is excluded."""
	fi = _write(tmp_path / 't.tsv', 'r\tx\ty\n')
	assert read_columns(fi, header=None) == ['1', '2']


def test_read_columns_reads_only_first_row(tmp_path):
	"""Rows after the first data row do not need to parse: only the header and one row are read."""
	import pandas as pd
	# The third row has more fields than the header, which pandas rejects when it reads the whole file
	fi = _write(tmp_path / 't.tsv', 'id\ta\tb\nr1\t1\t2\nr2\t1\t2\t3\t4\t5\n')
	with pytest.raises(pd.errors.ParserError):
		pd.read_csv(fi, sep='\t', header=0, index_col=0)
	assert read_columns(fi) == ['a', 'b']


# read_index

def test_read_index_returns_row_names(tmp_path):
	"""Row names from the first column, as str, in file order."""
	fi = _write(tmp_path / 't.tsv', 'id\tc1\tc2\ngA\t1\t2\ngB\t3\t4\n7\t5\t6\n')
	assert read_index(fi) == ['gA', 'gB', '7']


def test_read_index_unique(tmp_path):
	"""unique=True rejects a repeated row name."""
	fi = _write(tmp_path / 't.tsv', 'id\tc\ngA\t1\ngA\t2\n')
	assert read_index(fi) == ['gA', 'gA']
	with pytest.raises(AssertionError):
		read_index(fi, unique=True)


def test_read_index_custom_separator(tmp_path):
	"""sep is passed to pandas.read_csv."""
	fi = _write(tmp_path / 't.csv', 'id,c\nr1,1\nr2,2\n')
	assert read_index(fi, sep=',') == ['r1', 'r2']


def test_read_index_without_header(tmp_path):
	"""header=None treats the first line as a row."""
	fi = _write(tmp_path / 't.tsv', 'r1\t1\nr2\t2\n')
	assert read_index(fi, header=None) == ['r1', 'r2']


def test_read_index_other_index_column_not_implemented(tmp_path):
	"""Only index_col=0 is supported."""
	fi = _write(tmp_path / 't.tsv', 'id\tc\nr1\t1\n')
	with pytest.raises(NotImplementedError):
		read_index(fi, index_col=1)
