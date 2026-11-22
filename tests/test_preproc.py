"""
Tests for dictys.preproc: selects_rna, selects_atac and qc_reads, called
directly and through ``dictys preproc ...``.
"""

import os

import numpy as np
import pandas as pd
import pytest

from dictys import preproc


def _write_matrix(path, values, genes, cells):
	pd.DataFrame(np.asarray(values), index=genes, columns=cells).to_csv(path, sep='\t')
	return str(path)


def _write_lines(path, lines):
	with open(path, 'w') as f:
		f.write('\n'.join(lines) + '\n')
	return str(path)


def _read_matrix(path):
	return pd.read_csv(path, sep='\t', header=0, index_col=0)


def _read_names(path):
	with open(path) as f:
		return [x.strip() for x in f.read().splitlines() if x.strip()]


################################################################
# selects_rna
################################################################

@pytest.fixture
def rna_matrix(tmp_path):
	"""Four genes by four cells; g2 is zero in c1 and c3, g4 is zero everywhere but c2."""
	values = [
		[1, 0, 2, 0],
		[0, 5, 0, 3],
		[4, 4, 4, 4],
		[0, 7, 0, 0],
	]
	return _write_matrix(tmp_path / 'reads.tsv', values, ['g1', 'g2', 'g3', 'g4'], ['c1', 'c2', 'c3', 'c4'])


def test_selects_rna_keeps_named_cells_and_expressed_genes(tmp_path, rna_matrix):
	"""Selecting c3 and c1 keeps those columns in matrix order and drops genes all-zero in them (g2, g4)."""
	names = _write_lines(tmp_path / 'names.txt', ['c3', 'cX', 'c1'])
	fo = tmp_path / 'out.tsv'
	preproc.selects_rna(rna_matrix, names, str(fo))
	d = _read_matrix(fo)
	assert list(d.columns) == ['c1', 'c3']
	assert list(d.index) == ['g1', 'g3']
	np.testing.assert_array_equal(d.values, [[1, 2], [4, 4]])


def test_selects_rna_keeps_all_cells_and_genes(tmp_path, rna_matrix):
	"""Selecting every cell keeps the whole matrix, because every gene is expressed in some cell."""
	names = _write_lines(tmp_path / 'names.txt', ['c1', '', '  c2  ', 'c3', 'c4'])
	fo = tmp_path / 'out.tsv'
	preproc.selects_rna(rna_matrix, names, str(fo))
	pd.testing.assert_frame_equal(_read_matrix(fo), _read_matrix(rna_matrix))


def test_selects_rna_cli(tmp_path, rna_matrix, run_dictys):
	"""``dictys preproc selects_rna`` writes the same selection as the function: c2 and c4 keep g2, g3, g4."""
	names = _write_lines(tmp_path / 'names.txt', ['c2', 'c4'])
	fo = tmp_path / 'out.tsv'
	run_dictys('preproc', 'selects_rna', rna_matrix, names, fo)
	d = _read_matrix(fo)
	assert list(d.columns) == ['c2', 'c4']
	assert list(d.index) == ['g2', 'g3', 'g4']
	np.testing.assert_array_equal(d.values, [[5, 3], [4, 4], [7, 0]])


################################################################
# selects_atac
################################################################

def test_selects_atac_filters_by_expression_columns(tmp_path, rna_matrix):
	"""Names absent from the expression columns and blank lines are dropped; list order is kept."""
	fi_list = _write_lines(tmp_path / 'list.txt', ['c3', 'cX', '', 'c1', 'c9'])
	fo = tmp_path / 'out.txt'
	preproc.selects_atac(rna_matrix, fi_list, str(fo))
	assert _read_names(fo) == ['c3', 'c1']


def test_selects_atac_none_selected_raises(tmp_path, rna_matrix):
	"""A list sharing no name with the expression columns raises RuntimeError and writes nothing."""
	fi_list = _write_lines(tmp_path / 'list.txt', ['cX', 'cY'])
	fo = tmp_path / 'out.txt'
	with pytest.raises(RuntimeError):
		preproc.selects_atac(rna_matrix, fi_list, str(fo))
	assert not os.path.exists(fo)


def test_selects_atac_cli(tmp_path, rna_matrix, run_dictys):
	"""``dictys preproc selects_atac`` keeps the listed cells present in the expression matrix."""
	fi_list = _write_lines(tmp_path / 'list.txt', ['c4', 'c2', 'cZ'])
	fo = tmp_path / 'out.txt'
	run_dictys('preproc', 'selects_atac', rna_matrix, fi_list, fo)
	assert _read_names(fo) == ['c4', 'c2']


################################################################
# qc_reads
################################################################

# Genes by cells. Gene totals: g1 10, g2 1, g3 16, g4 2; expressed cells: 3, 1, 4, 1.
# Cell totals: c1 9, c2 6, c3 8, c4 6; expressed genes: 2, 2, 3, 2.
QC_VALUES = [
	[5, 0, 3, 2],
	[0, 0, 1, 0],
	[4, 4, 4, 4],
	[0, 2, 0, 0],
]
QC_GENES = ['g1', 'g2', 'g3', 'g4']
QC_CELLS = ['c1', 'c2', 'c3', 'c4']


@pytest.fixture
def qc_matrix(tmp_path):
	"""The QC_VALUES matrix written as a tsv file."""
	return _write_matrix(tmp_path / 'reads.tsv', QC_VALUES, QC_GENES, QC_CELLS)


def _expected_qc(genes, cells):
	d = pd.DataFrame(QC_VALUES, index=QC_GENES, columns=QC_CELLS)
	return d.loc[genes, cells]


@pytest.mark.parametrize('params,genes,cells', [
	# all criteria disabled
	(dict(), QC_GENES, QC_CELLS),
	# gene total >= 2 drops g2 (total 1)
	(dict(n_gene=2), ['g1', 'g3', 'g4'], QC_CELLS),
	# gene expressed in >= 2 cells drops g2 and g4
	(dict(nc_gene=2), ['g1', 'g3'], QC_CELLS),
	# gene expressed in >= 0.75 * 4 = 3 cells keeps g1 and g3
	(dict(ncp_gene=0.75), ['g1', 'g3'], QC_CELLS),
	# cell total >= 7 keeps c1 (9) and c3 (8)
	(dict(n_cell=7), QC_GENES, ['c1', 'c3']),
	# cell expressing >= 3 genes keeps c3 only
	(dict(nt_cell=3), QC_GENES, ['c3']),
	# cell expressing >= 0.6 * 4 = 2.4 genes keeps c3 only
	(dict(ntp_cell=0.6), QC_GENES, ['c3']),
	# cell expressing >= 0.5 * 4 = 2 genes keeps every cell
	(dict(ntp_cell=0.5), QC_GENES, QC_CELLS),
])
def test_qc_reads_single_criterion(tmp_path, qc_matrix, params, genes, cells):
	"""Each QC bound alone removes the genes or cells worked out by hand from the matrix totals."""
	ka = dict(n_gene=0, nc_gene=0, ncp_gene=0, n_cell=0, nt_cell=0, ntp_cell=0)
	ka.update(params)
	fo = tmp_path / 'out.tsv'
	preproc.qc_reads(qc_matrix, str(fo), **ka)
	pd.testing.assert_frame_equal(_read_matrix(fo), _expected_qc(genes, cells))


def test_qc_reads_iterates_until_stable(tmp_path):
	"""g2 passes n_gene=3 at first but fails once c2 (total 4 < 5) is removed, so it goes in the second round."""
	values = [[5, 0, 3], [0, 3, 0], [4, 1, 4]]
	fi = _write_matrix(tmp_path / 'reads.tsv', values, ['g1', 'g2', 'g3'], ['c1', 'c2', 'c3'])
	fo = tmp_path / 'out.tsv'
	preproc.qc_reads(fi, str(fo), 3, 0, 0, 5, 0, 0)
	d = _read_matrix(fo)
	assert list(d.index) == ['g1', 'g3']
	assert list(d.columns) == ['c1', 'c3']
	np.testing.assert_array_equal(d.values, [[5, 3], [4, 4]])


def test_qc_reads_ncp_gene_uses_remaining_cells(tmp_path):
	"""ncp_gene is a proportion of the cells left: after c3 and c4 go (nt_cell=2), g2 has 1 of 2 cells and passes 0.5."""
	values = [[1, 1, 0, 0], [1, 0, 1, 0], [0, 1, 0, 1]]
	fi = _write_matrix(tmp_path / 'reads.tsv', values, ['g1', 'g2', 'g3'], ['c1', 'c2', 'c3', 'c4'])
	fo = tmp_path / 'out.tsv'
	# Round 1 on 4 cells: genes need >= 2 expressed cells (all pass); cells need >= 2 genes: c1, c2.
	# Round 2 on 2 cells: genes need >= 1: g1 (2), g2 (1), g3 (1) all pass; cells c1, c2 have 2 genes each.
	preproc.qc_reads(fi, str(fo), 0, 0, 0.5, 0, 2, 0)
	d = _read_matrix(fo)
	assert list(d.index) == ['g1', 'g2', 'g3']
	assert list(d.columns) == ['c1', 'c2']


def test_qc_reads_all_genes_removed(tmp_path, qc_matrix):
	"""A gene bound above every gene total raises RuntimeError."""
	with pytest.raises(RuntimeError):
		preproc.qc_reads(qc_matrix, str(tmp_path / 'out.tsv'), 100, 0, 0, 0, 0, 0)


def test_qc_reads_all_cells_removed(tmp_path, qc_matrix):
	"""A cell bound above every cell total raises RuntimeError."""
	with pytest.raises(RuntimeError):
		preproc.qc_reads(qc_matrix, str(tmp_path / 'out.tsv'), 0, 0, 0, 100, 0, 0)


@pytest.mark.parametrize('index', range(6))
def test_qc_reads_negative_parameter(tmp_path, qc_matrix, index):
	"""Any negative QC parameter raises ValueError."""
	args = [0, 0, 0, 0, 0, 0]
	args[index] = -1
	with pytest.raises(ValueError):
		preproc.qc_reads(qc_matrix, str(tmp_path / 'out.tsv'), *args)


@pytest.mark.parametrize('index', [2, 5])
def test_qc_reads_proportion_above_one(tmp_path, qc_matrix, index):
	"""ncp_gene or ntp_cell above 1 raises ValueError."""
	args = [0, 0, 0, 0, 0, 0]
	args[index] = 1.5
	with pytest.raises(ValueError):
		preproc.qc_reads(qc_matrix, str(tmp_path / 'out.tsv'), *args)


def test_qc_reads_cli(tmp_path, qc_matrix, run_dictys):
	"""``dictys preproc qc_reads`` with every bound set: n_gene=2, nc_gene=1, ncp_gene=0.25, n_cell=7, nt_cell=2, ntp_cell=0.5."""
	fo = tmp_path / 'out.tsv'
	# Round 1: genes need total >= 2 and >= 1 cell: g2 goes. Cells need total >= 7 and >= 2 genes: c1, c3.
	# Round 2 on g1, g3, g4 x c1, c3 = [[5, 3], [4, 4], [0, 0]]: g4 total 0 goes; c1 and c3 keep 2 genes.
	# Round 3 on g1, g3 x c1, c3: stable.
	run_dictys('preproc', 'qc_reads', qc_matrix, fo, 2, 1, 0.25, 7, 2, 0.5)
	pd.testing.assert_frame_equal(_read_matrix(fo), _expected_qc(['g1', 'g3'], ['c1', 'c3']))
