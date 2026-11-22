"""
Tests for the Python helper scripts under dictys/scripts/helper:
expression_mtx.py, expression.py, makefile_update.py, makefile_check.py and
split_bam_text.py. Each script runs in this process through ``run_helper``.
"""

import gzip
import logging
import os
import re
import shutil

import numpy as np
import pandas as pd
import pytest


def read_lines(path):
	"""Non-empty lines of a text file, independent of the platform's line separator."""
	with open(path, 'r') as f:
		return [x for x in f.read().splitlines() if x]


def write_lines(path, lines):
	"""Write each string in ``lines`` as one line."""
	with open(path, 'w') as f:
		f.write(''.join(x + '\n' for x in lines))


######################################################################
# expression_mtx.py
######################################################################

def write_mtx_folder(path, mat, genes, barcodes, nbarcode_cols=1):
	"""Write a 10x-style folder: matrix.mtx.gz (genes x cells), features.tsv.gz and barcodes.tsv.gz."""
	import scipy.io
	import scipy.sparse
	os.makedirs(path, exist_ok=True)
	fraw = os.path.join(path, 'matrix.mtx')
	scipy.io.mmwrite(fraw, scipy.sparse.coo_matrix(np.asarray(mat)))
	with open(fraw, 'rb') as f, gzip.open(os.path.join(path, 'matrix.mtx.gz'), 'wb') as fo:
		fo.write(f.read())
	os.remove(fraw)
	with gzip.open(os.path.join(path, 'features.tsv.gz'), 'wt') as f:
		for i, g in enumerate(genes):
			f.write('ENSG{:04d}\t{}\tGene Expression\n'.format(i, g))
	with gzip.open(os.path.join(path, 'barcodes.tsv.gz'), 'wt') as f:
		for b in barcodes:
			f.write('\t'.join([b] * nbarcode_cols) + '\n')


class TestExpressionMtx:
	def test_converts_dedups_filters_and_sorts(self, run_helper, tmp_path, caplog):
		"""Output keeps first occurrences, drops names with ':' or '.', sorts both axes (numpy indexing reference)."""
		mat = np.arange(24).reshape(6, 4) % 5
		genes = ['Zeb', 'Abc', 'Zeb', 'mt:X', 'Gene.1', 'Mid']
		cells = ['c3', 'c1', 'c2', 'c1']
		din = str(tmp_path / 'in')
		write_mtx_folder(din, mat, genes, cells)
		fo = str(tmp_path / 'out.tsv.gz')
		caplog.set_level(logging.WARNING)
		run_helper('expression_mtx.py', din, fo)
		with open(fo, 'rb') as f:
			assert f.read(2) == b'\x1f\x8b'
		d = pd.read_csv(fo, header=0, index_col=0, sep='\t')
		assert list(d.index) == ['Abc', 'Mid', 'Zeb']
		assert list(d.columns) == ['c1', 'c2', 'c3']
		# rows Abc=1, Mid=5, Zeb=0 (first occurrence); columns c1=1, c2=2, c3=0
		np.testing.assert_array_equal(d.values, mat[[1, 5, 0]][:, [1, 2, 0]])
		assert 'Skipped duplicate occurrence of gene names: Zeb' in caplog.text
		assert 'Skipped duplicate occurrence of cell names: c1' in caplog.text

	def test_no_duplicates_no_warning(self, run_helper, tmp_path, caplog):
		"""Unique sorted names pass through unchanged with no warning."""
		mat = np.array([[1, 0], [0, 2], [3, 4]])
		din = str(tmp_path / 'in')
		write_mtx_folder(din, mat, ['A', 'B', 'C'], ['c1', 'c2'])
		fo = str(tmp_path / 'out.tsv.gz')
		caplog.set_level(logging.WARNING)
		run_helper('expression_mtx.py', din, fo)
		d = pd.read_csv(fo, header=0, index_col=0, sep='\t')
		assert list(d.index) == ['A', 'B', 'C']
		assert list(d.columns) == ['c1', 'c2']
		np.testing.assert_array_equal(d.values, mat)
		assert 'Skipped duplicate' not in caplog.text

	def test_barcodes_with_two_columns_rejected(self, run_helper, tmp_path):
		"""barcodes.tsv.gz must have a single column."""
		din = str(tmp_path / 'in')
		write_mtx_folder(din, np.eye(2), ['A', 'B'], ['c1', 'c2'], nbarcode_cols=2)
		with pytest.raises(AssertionError):
			run_helper('expression_mtx.py', din, str(tmp_path / 'out.tsv.gz'))

	def test_shape_mismatch_rejected(self, run_helper, tmp_path):
		"""A matrix whose shape differs from the feature and barcode counts is rejected."""
		din = str(tmp_path / 'in')
		write_mtx_folder(din, np.eye(2), ['A', 'B', 'C'], ['c1', 'c2'])
		with pytest.raises(AssertionError):
			run_helper('expression_mtx.py', din, str(tmp_path / 'out.tsv.gz'))

	def test_missing_output_argument(self, run_helper, tmp_path):
		"""The output file is a required positional argument."""
		with pytest.raises(SystemExit):
			run_helper('expression_mtx.py', str(tmp_path))


######################################################################
# expression.py
######################################################################

def test_expression_is_deprecated(run_helper):
	"""The deprecated expression.py helper raises and points to expression_mtx.py."""
	with pytest.raises(RuntimeError, match='expression_mtx.py'):
		run_helper('expression.py')


######################################################################
# makefile_update.py
######################################################################

class TestMakefileUpdate:
	def test_replace_keeps_operator_indent_and_newline(self, run_helper, tmp_path):
		"""Replacing a value keeps the line's '=' or ':=' operator, its indentation and its line end."""
		f = tmp_path / 'a.mk'
		f.write_text('#comment\nA=1\n  B:=x y\nC=keep\n')
		run_helper('makefile_update.py', str(f), '{"A":"2","B":"z"}')
		assert f.read_text() == '#comment\nA=2\n  B:=z\nC=keep\n'

	def test_append_with_plus_suffix(self, run_helper, tmp_path):
		"""A '+' suffix on the variable name appends the value to the current one."""
		f = tmp_path / 'a.mk'
		f.write_text('KP:=--nth 4\t \nX=1\n')
		run_helper('makefile_update.py', str(f), '{"KP+":" --device cpu"}')
		assert f.read_text() == 'KP:=--nth 4 --device cpu\t \nX=1\n'

	def test_last_assignment_is_updated(self, run_helper, tmp_path):
		"""Only the last assignment of a variable is changed."""
		f = tmp_path / 'a.mk'
		f.write_text('A=1\nB=2\nA:=3\n')
		run_helper('makefile_update.py', str(f), '{"A":"9"}')
		assert f.read_text() == 'A=1\nB=2\nA:=9\n'

	def test_prefix_names_untouched(self, run_helper, tmp_path):
		"""Variables whose names start with the target name are not matched."""
		f = tmp_path / 'a.mk'
		f.write_text('NTH=1\nNTH2=1\nNTH-X:=1\n')
		run_helper('makefile_update.py', str(f), '{"NTH":"8"}')
		assert f.read_text() == 'NTH=8\nNTH2=1\nNTH-X:=1\n'

	def test_variable_not_found(self, run_helper, tmp_path):
		"""Updating a variable that is never assigned raises ValueError and leaves the file unchanged."""
		f = tmp_path / 'a.mk'
		f.write_text('A=1\n')
		with pytest.raises(ValueError, match='Variable not found: Z'):
			run_helper('makefile_update.py', str(f), '{"Z":"1"}')
		assert f.read_text() == 'A=1\n'

	def test_final_rule_append_rejected(self, run_helper, tmp_path):
		"""A variable whose last occurrence is a '+=' line raises ValueError."""
		f = tmp_path / 'a.mk'
		f.write_text('A=1\nA+=2\n')
		with pytest.raises(ValueError, match='Final rule is not assignment: A'):
			run_helper('makefile_update.py', str(f), '{"A":"3"}')

	def test_invalid_json(self, run_helper, tmp_path):
		"""A malformed JSON update string raises ValueError."""
		f = tmp_path / 'a.mk'
		f.write_text('A=1\n')
		with pytest.raises(ValueError):
			run_helper('makefile_update.py', str(f), '{A:1')

	def test_missing_argument(self, run_helper, tmp_path):
		"""Both positional arguments are required."""
		with pytest.raises(SystemExit):
			run_helper('makefile_update.py', str(tmp_path / 'a.mk'))

	def test_package_config(self, run_helper, dictys_scripts, tmp_path):
		"""On the shipped config.mk only the JOINT and DEVICE lines change."""
		src = os.path.join(dictys_scripts, 'makefiles', 'config.mk')
		f = tmp_path / 'config.mk'
		shutil.copy(src, str(f))
		before = f.read_text().splitlines(True)
		run_helper('makefile_update.py', str(f), '{"JOINT":"0","DEVICE":"cuda:0"}')
		after = f.read_text().splitlines(True)
		assert len(after) == len(before)
		changed = [(x, y) for x, y in zip(before, after) if x != y]
		assert sorted(y for _, y in changed) == ['DEVICE=cuda:0\n', 'JOINT=0\n']


######################################################################
# makefile_check.py
######################################################################

def set_joint(dir_makefiles, value):
	"""Rewrite the JOINT= line of config.mk."""
	f = os.path.join(dir_makefiles, 'config.mk')
	with open(f) as fh:
		s = fh.read()
	s = re.sub(r'(?m)^JOINT=.*$', 'JOINT=' + value, s)
	with open(f, 'w') as fh:
		fh.write(s)


def write_expression(path, genes, cells):
	"""Write a genes x cells expression.tsv.gz with small integer counts."""
	d = pd.DataFrame(np.arange(len(genes) * len(cells)).reshape(len(genes), len(cells)) % 7, index=genes, columns=cells)
	d.to_csv(path, sep='\t')


def write_motifs(path, names):
	"""Write a HOMER motif file with one 1-position motif per name."""
	lines = []
	for i, x in enumerate(names):
		lines += ['>ACGT\t{}\t5.0'.format(x), '0.25\t0.25\t0.25\t0.25']
	write_lines(path, lines)


def build_pipeline(root, scripts, joint=True, ncell=100, ngene=100, ntf=10, static=True, dynamic=True,
		covariate=True, blacklist=True, dir_data='data', dir_makefiles='makefiles'):
	"""
	Write a minimal pipeline folder that makefile_check.py accepts: the package makefiles,
	expression of ngene genes in ncell cells, one empty bam per ATAC cell, motifs for the first
	ntf genes plus one motif for TF 'TFX' absent from the data and one shared by G000 and G001,
	gene.bed, a genome, and optionally covariates, a blacklist, two subsets and a one-edge trajectory.
	"""
	root = str(root)
	dmk = os.path.join(root, dir_makefiles)
	os.makedirs(dmk)
	src = os.path.join(scripts, 'makefiles')
	for x in os.listdir(src):
		if x.endswith('.mk'):
			shutil.copy(os.path.join(src, x), dmk)
	if not static:
		os.remove(os.path.join(dmk, 'static.mk'))
	if not dynamic:
		os.remove(os.path.join(dmk, 'dynamic.mk'))
	set_joint(dmk, '1' if joint else '0')

	dd = os.path.join(root, dir_data)
	os.makedirs(dd)
	genes = ['G{:03d}'.format(x) for x in range(ngene)]
	cells_rna = ['R{:03d}'.format(x) for x in range(ncell)]
	cells_atac = cells_rna if joint else ['A{:03d}'.format(x) for x in range(ncell)]
	write_expression(os.path.join(dd, 'expression.tsv.gz'), genes, cells_rna)
	os.makedirs(os.path.join(dd, 'bams'))
	for x in cells_atac:
		open(os.path.join(dd, 'bams', x + '.bam'), 'w').close()
	open(os.path.join(dd, 'bams', 'README.txt'), 'w').close()
	motifs = ['{}_m{}'.format(x, i) for i, x in enumerate(genes[:ntf])] + ['TFX_extra', 'G000,G001_shared']
	write_motifs(os.path.join(dd, 'motifs.motif'), motifs)
	bed = ['chr1\t{}\t{}\t{}\t.\t+'.format(1000 * i, 1000 * i + 500, x) for i, x in enumerate(genes)]
	write_lines(os.path.join(dd, 'gene.bed'), bed + ['chr2\t0\t10\tNOTEXPR\t.\t-'])
	os.makedirs(os.path.join(dd, 'genome'))
	write_lines(os.path.join(dd, 'genome', 'chr1.fa'), ['>chr1', 'ACGT'])
	if blacklist:
		write_lines(os.path.join(dd, 'blacklist.bed'), ['chr1\t0\t100'])
	if covariate:
		pd.DataFrame(np.ones((ncell, 2)), index=cells_rna, columns=['cov1', 'cov2']).to_csv(
			os.path.join(dd, 'covariate.tsv.gz'), sep='\t')
	if static:
		write_lines(os.path.join(dd, 'subsets.txt'), ['S1', 'S2'])
		half = ncell // 2
		for name, sl in [('S1', slice(0, half)), ('S2', slice(half, None))]:
			os.makedirs(os.path.join(dd, 'subsets', name))
			write_lines(os.path.join(dd, 'subsets', name, 'names_rna.txt'), cells_rna[sl])
			write_lines(os.path.join(dd, 'subsets', name, 'names_atac.txt'), cells_atac[sl])
	if dynamic:
		from dictys.traj import trajectory, point
		traj = trajectory(np.array([[0, 1]]), np.array([1.0]))
		traj.to_file(os.path.join(dd, 'traj_node.h5'))
		modes = [('rna', cells_rna)] + ([] if joint else [('atac', cells_atac)])
		for mode, cells in modes:
			point(traj, np.zeros(len(cells), dtype=int), np.linspace(0, 1, len(cells))).to_file(
				os.path.join(dd, 'traj_cell_{}.h5'.format(mode)))
			pd.DataFrame(np.ones((len(cells), 2)), index=cells, columns=['x', 'y']).to_csv(
				os.path.join(dd, 'coord_{}.tsv.gz'.format(mode)), sep='\t')
	return {'root': root, 'data': dd, 'makefiles': dmk, 'genes': genes, 'cells_rna': cells_rna,
		'cells_atac': cells_atac}


def read_tsv(path):
	return pd.read_csv(path, header=0, index_col=0, sep='\t')


def mut_dup_gene_expression(p):
	f = os.path.join(p['data'], 'expression.tsv.gz')
	d = read_tsv(f)
	pd.concat([d, d.iloc[[50]]]).to_csv(f, sep='\t')


def mut_no_bams(p):
	shutil.rmtree(os.path.join(p['data'], 'bams'))


def mut_missing_bam(p):
	os.remove(os.path.join(p['data'], 'bams', 'R000.bam'))


def mut_dup_motif(p):
	with open(os.path.join(p['data'], 'motifs.motif'), 'a') as f:
		f.write('>CCCC\tG005_m5\t5.0\n0.25\t0.25\t0.25\t0.25\n')


def mut_empty_tf(p):
	with open(os.path.join(p['data'], 'motifs.motif'), 'a') as f:
		f.write('>CCCC\t_bad\t5.0\n0.25\t0.25\t0.25\t0.25\n')


def mut_bed_columns(p):
	with open(os.path.join(p['data'], 'gene.bed'), 'a') as f:
		f.write('chr2\t20\t30\tNOTEXPR2\t.\n')


def mut_bed_dup(p):
	with open(os.path.join(p['data'], 'gene.bed'), 'a') as f:
		f.write('chr1\t500000\t500100\tG000\t.\t+\n')


def mut_bed_few(p):
	f = os.path.join(p['data'], 'gene.bed')
	write_lines(f, read_lines(f)[1:])


def mut_cov_dup_cell(p):
	f = os.path.join(p['data'], 'covariate.tsv.gz')
	d = read_tsv(f)
	pd.concat([d, d.iloc[[3]]]).to_csv(f, sep='\t')


def mut_cov_missing_cell(p):
	f = os.path.join(p['data'], 'covariate.tsv.gz')
	read_tsv(f).iloc[1:].to_csv(f, sep='\t')


def mut_subset_rna_unknown(p):
	for x in ['names_rna.txt', 'names_atac.txt']:
		with open(os.path.join(p['data'], 'subsets', 'S1', x), 'a') as f:
			f.write('NOPE\n')


def mut_subset_joint_mismatch(p):
	f = os.path.join(p['data'], 'subsets', 'S1', 'names_atac.txt')
	write_lines(f, read_lines(f)[1:])


def mut_subset_atac_unknown(p):
	with open(os.path.join(p['data'], 'subsets', 'S1', 'names_atac.txt'), 'a') as f:
		f.write('NOPE\n')


def mut_coord_unknown(mode):
	def f(p):
		fn = os.path.join(p['data'], 'coord_{}.tsv.gz'.format(mode))
		d = read_tsv(fn)
		d.index = ['NOPE'] + list(d.index[1:])
		d.to_csv(fn, sep='\t')
	return f


def mut_coord_count(mode):
	def f(p):
		fn = os.path.join(p['data'], 'coord_{}.tsv.gz'.format(mode))
		read_tsv(fn).iloc[1:].to_csv(fn, sep='\t')
	return f


# (id, build keyword arguments, mutation, exception without -c, error count with -c, logged message)
CHECK_FAILURES = [
	('dup_gene_expression', {}, mut_dup_gene_expression, ValueError, 1, 'Duplicate gene names found in expression.tsv.gz'),
	('few_cells', {'ncell': 99}, None, ValueError, 1, '<100 cells found with RNA profile'),
	('few_genes', {'ngene': 99}, None, ValueError, 2, '<100 genes found with RNA profile'),
	('no_bams', {}, mut_no_bams, FileNotFoundError, 1, 'Using RNA cell names for ATAC cell names'),
	('missing_bam_joint', {}, mut_missing_bam, FileNotFoundError, 1, 'Not all cells with RNA profile'),
	('dup_motif', {}, mut_dup_motif, ValueError, 1, 'Duplicate motif names found'),
	('empty_tf', {}, mut_empty_tf, ValueError, 1, 'TF name "" found in motif file'),
	('few_tfs', {'ntf': 9}, None, ValueError, 1, '<10 TFs found in motifs.motif file'),
	('bed_columns', {}, mut_bed_columns, ValueError, 1, 'Unequal number of columns in gene.bed'),
	('bed_dup', {}, mut_bed_dup, ValueError, 1, 'Duplicate gene names found in gene.bed'),
	('bed_few', {}, mut_bed_few, ValueError, 1, '<100 overlapping genes with TSS information'),
	('cov_dup_cell', {}, mut_cov_dup_cell, ValueError, 1, 'Duplicate cell names found in covariate.tsv.gz'),
	('cov_missing_cell', {}, mut_cov_missing_cell, ValueError, 1, 'Covariate information cannot be found'),
	('subset_rna_unknown', {}, mut_subset_rna_unknown, ValueError, 1, 'Subset RNA cells in subsets/*/names_rna.txt could not be found'),
	('subset_joint_mismatch', {}, mut_subset_joint_mismatch, ValueError, 1, 'must be identical in each subset for joint profiles'),
	('subset_atac_unknown', {'joint': False}, mut_subset_atac_unknown, ValueError, 1, 'Subset ATAC cells in subsets/*/names_atac.txt could not be found'),
	('coord_rna_unknown', {}, mut_coord_unknown('rna'), ValueError, 1, 'Low dimensional RNA cells in coord_rna.tsv.gz cannot be found'),
	('coord_rna_count', {}, mut_coord_count('rna'), ValueError, 1, 'Inconsistent RNA cell count'),
	('coord_atac_unknown', {'joint': False}, mut_coord_unknown('atac'), ValueError, 1, 'Low dimensional ATAC cells in coord_atac.tsv.gz cannot be found'),
	('coord_atac_count', {'joint': False}, mut_coord_count('atac'), ValueError, 1, 'Inconsistent ATAC cell count'),
]


def valid_output(joint):
	"""Lines makefile_check.py prints for the default build_pipeline layout, counted from what it writes."""
	return [
		'Joint profile: {}'.format(joint),
		'Found 100 cells with RNA profile',
		'Found 100 genes with RNA profile',
		'Found 100 cells with ATAC profile',
		'Found 12 motifs',
		'Found 11 TFs',
		'Found 10 TFs in current dataset',
		'Missing 1 TFs in current dataset: TFX',
		'Found 100 genes with TSS information',
		'Found 2 covariates',
	]


class TestMakefileCheck:
	@pytest.mark.parametrize('joint', [True, False])
	def test_valid_layout_passes(self, run_helper, dictys_scripts, tmp_path, capsys, caplog, joint):
		"""A complete valid layout passes and the printed counts match the files written."""
		p = build_pipeline(tmp_path, dictys_scripts, joint=joint)
		caplog.set_level(logging.WARNING)
		run_helper('makefile_check.py', cwd=tmp_path)
		assert capsys.readouterr().out.splitlines() == valid_output(joint)
		assert not [x for x in caplog.records if x.levelno >= logging.ERROR]
		assert 'Skipping' not in caplog.text

	def test_valid_layout_passes_with_continue(self, run_helper, dictys_scripts, tmp_path, caplog):
		"""With -c a valid layout raises nothing and logs no error."""
		build_pipeline(tmp_path, dictys_scripts)
		caplog.set_level(logging.WARNING)
		g = run_helper('makefile_check.py', '-c', cwd=tmp_path)
		assert g['nerr'] == 0
		assert not [x for x in caplog.records if x.levelno >= logging.ERROR]

	def test_custom_folders(self, run_helper, dictys_scripts, tmp_path, capsys):
		"""--dir_data and --dir_makefiles locate folders with other names, and every check reads from them."""
		build_pipeline(tmp_path, dictys_scripts, joint=False, dir_data='mydata', dir_makefiles='mymk')
		run_helper('makefile_check.py', '--dir_data', str(tmp_path / 'mydata'), '--dir_makefiles',
			str(tmp_path / 'mymk'), cwd=tmp_path)
		assert capsys.readouterr().out.splitlines() == valid_output(False)

	def test_missing_makefiles_folder(self, run_helper, dictys_scripts, tmp_path):
		"""A missing makefiles folder raises FileNotFoundError even with -c."""
		build_pipeline(tmp_path, dictys_scripts, dir_makefiles='mymk')
		with pytest.raises(FileNotFoundError, match='--dir_makefiles'):
			run_helper('makefile_check.py', '-c', cwd=tmp_path)

	def test_missing_data_folder(self, run_helper, dictys_scripts, tmp_path):
		"""A missing data folder raises FileNotFoundError even with -c."""
		build_pipeline(tmp_path, dictys_scripts, dir_data='mydata')
		with pytest.raises(FileNotFoundError, match='--dir_data'):
			run_helper('makefile_check.py', '-c', cwd=tmp_path)

	@pytest.mark.parametrize('cont', [False, True])
	def test_invalid_joint(self, run_helper, dictys_scripts, tmp_path, cont):
		"""A JOINT value other than 0 or 1 raises ValueError, with or without -c."""
		p = build_pipeline(tmp_path, dictys_scripts)
		set_joint(p['makefiles'], '2')
		with pytest.raises(ValueError, match='Invalid JOINT variable'):
			run_helper('makefile_check.py', *(['-c'] if cont else []), cwd=tmp_path)

	@pytest.mark.parametrize('case', CHECK_FAILURES, ids=[x[0] for x in CHECK_FAILURES])
	def test_failure_stops(self, run_helper, dictys_scripts, tmp_path, case):
		"""Without -c each documented problem raises its exception."""
		_, ka, mut, exc, _, _ = case
		p = build_pipeline(tmp_path, dictys_scripts, **ka)
		if mut is not None:
			mut(p)
		with pytest.raises(exc):
			run_helper('makefile_check.py', cwd=tmp_path)

	@pytest.mark.parametrize('case', CHECK_FAILURES, ids=[x[0] for x in CHECK_FAILURES])
	def test_failure_counted_with_continue(self, run_helper, dictys_scripts, tmp_path, caplog, case):
		"""With -c each documented problem is logged and counted in the final RuntimeError."""
		_, ka, mut, _, nerr, msg = case
		p = build_pipeline(tmp_path, dictys_scripts, **ka)
		if mut is not None:
			mut(p)
		caplog.set_level(logging.WARNING)
		with pytest.raises(RuntimeError, match=r'^Found {} error\(s\) in total\.$'.format(nerr)):
			run_helper('makefile_check.py', '-c', cwd=tmp_path)
		assert msg in caplog.text

	def test_errors_accumulate_with_continue(self, run_helper, dictys_scripts, tmp_path, caplog):
		"""Independent problems add up under -c: duplicate motif, too few TFs and a duplicate bed gene."""
		p = build_pipeline(tmp_path, dictys_scripts, ntf=9)
		mut_dup_motif(p)
		mut_bed_dup(p)
		caplog.set_level(logging.WARNING)
		with pytest.raises(RuntimeError, match=r'Found 3 error\(s\) in total'):
			run_helper('makefile_check.py', '-c', cwd=tmp_path)
		assert len([x for x in caplog.records if x.levelno == logging.ERROR]) == 3

	def test_warning_no_strand(self, run_helper, dictys_scripts, tmp_path, caplog):
		"""A five-column gene.bed passes with a warning about missing strands."""
		p = build_pipeline(tmp_path, dictys_scripts)
		f = os.path.join(p['data'], 'gene.bed')
		write_lines(f, ['\t'.join(x.split('\t')[:5]) for x in read_lines(f)])
		caplog.set_level(logging.WARNING)
		run_helper('makefile_check.py', cwd=tmp_path)
		assert 'No strand information in gene.bed' in caplog.text

	def test_warning_chromosome_not_in_genome(self, run_helper, dictys_scripts, tmp_path, caplog):
		"""gene.bed chromosomes absent from genome/*.fa are reported as a warning only."""
		p = build_pipeline(tmp_path, dictys_scripts)
		write_lines(os.path.join(p['data'], 'genome', 'chr1.fa'), ['>chr5', 'ACGT'])
		caplog.set_level(logging.WARNING)
		run_helper('makefile_check.py', cwd=tmp_path)
		assert 'Chromosomes not found in reference genome: chr1' in caplog.text

	def test_warning_many_motifs(self, run_helper, dictys_scripts, tmp_path, caplog, capsys):
		"""More than 10000 motifs trigger the slow-inference warning."""
		p = build_pipeline(tmp_path, dictys_scripts)
		with open(os.path.join(p['data'], 'motifs.motif'), 'a') as f:
			for i in range(10000):
				f.write('>A\tEXTRA{}_x\t1\n0.25\t0.25\t0.25\t0.25\n'.format(i))
		caplog.set_level(logging.WARNING)
		run_helper('makefile_check.py', cwd=tmp_path)
		assert 'Found >10000 motifs' in caplog.text
		assert 'Found 10012 motifs' in capsys.readouterr().out

	@pytest.mark.parametrize('missing', ['static.mk', 'subsets.txt'])
	def test_static_checks_skipped(self, run_helper, dictys_scripts, tmp_path, caplog, missing):
		"""Without static.mk or subsets.txt the static checks are skipped, so bad subsets pass."""
		p = build_pipeline(tmp_path, dictys_scripts)
		mut_subset_rna_unknown(p)
		os.remove(os.path.join(p['makefiles'] if missing == 'static.mk' else p['data'], missing))
		caplog.set_level(logging.WARNING)
		run_helper('makefile_check.py', cwd=tmp_path)
		assert 'Skipping static network inference checks' in caplog.text

	@pytest.mark.parametrize('missing', ['dynamic.mk', 'traj_node.h5'])
	def test_dynamic_checks_skipped(self, run_helper, dictys_scripts, tmp_path, caplog, missing):
		"""Without dynamic.mk or traj_node.h5 the dynamic checks are skipped, so bad coordinates pass."""
		p = build_pipeline(tmp_path, dictys_scripts)
		mut_coord_count('rna')(p)
		os.remove(os.path.join(p['makefiles'] if missing == 'dynamic.mk' else p['data'], missing))
		caplog.set_level(logging.WARNING)
		run_helper('makefile_check.py', cwd=tmp_path)
		assert 'Skipping dynamic network inference checks' in caplog.text

	@pytest.mark.parametrize('joint', [True, False])
	def test_warning_cells_in_two_subsets(self, run_helper, dictys_scripts, tmp_path, caplog, joint):
		"""A cell listed in two subsets is a warning naming both subsets."""
		p = build_pipeline(tmp_path, dictys_scripts, joint=joint)
		for x in ['names_rna.txt', 'names_atac.txt']:
			first = read_lines(os.path.join(p['data'], 'subsets', 'S1', x))[0]
			with open(os.path.join(p['data'], 'subsets', 'S2', x), 'a') as f:
				f.write(first + '\n')
		caplog.set_level(logging.WARNING)
		run_helper('makefile_check.py', cwd=tmp_path)
		assert 'Found RNA cells in appearing multiple times' in caplog.text
		assert 'R000: S1,S2' in caplog.text
		if not joint:
			assert 'Found ATAC cells appearing multiple times' in caplog.text
			assert 'A000: S1,S2' in caplog.text

	@pytest.mark.parametrize('joint', [True, False])
	def test_warning_unassigned_cells(self, run_helper, dictys_scripts, tmp_path, caplog, joint):
		"""Cells in no subset are a warning, for RNA and (separate profiles) for ATAC."""
		p = build_pipeline(tmp_path, dictys_scripts, joint=joint)
		for x in ['names_rna.txt', 'names_atac.txt']:
			f = os.path.join(p['data'], 'subsets', 'S1', x)
			write_lines(f, read_lines(f)[1:])
		caplog.set_level(logging.WARNING)
		run_helper('makefile_check.py', cwd=tmp_path)
		assert 'Found RNA cells in expression.tsv.gz not assigned to any subset' in caplog.text
		assert ('Found ATAC cells in bams folder not assigned to any subset' in caplog.text) == (not joint)


######################################################################
# split_bam_text.py
######################################################################

def sam_line(name, *tags):
	"""One SAM text line (as from `samtools view`) with the given optional fields."""
	return '\t'.join([name, '0', 'chr1', '100', '60', '4M', '*', '0', '0', 'ACGT', 'IIII'] + list(tags))


def sam_input(cells):
	"""SAM text for reads tagged CB:Z:<cell>; None gives an untagged read. Returns (text, expected dict)."""
	lines = []
	expected = {}
	for i, c in enumerate(cells):
		line = sam_line('r{}'.format(i), *([] if c is None else ['CB:Z:' + c]), 'NM:i:0')
		lines.append(line)
		expected.setdefault(c, []).append(line)
	return ''.join(x + '\n' for x in lines), expected


CELLS = ['A', 'B', None, 'A', 'C', 'B', None, 'A', 'C', 'A']


class TestSplitBamText:
	def test_split_by_section(self, run_helper, tmp_path):
		"""Reads go to one file per barcode in input order; reads without a barcode are dropped."""
		text, expected = sam_input(CELLS)
		dout = tmp_path / 'out'
		run_helper('split_bam_text.py', '--section', 'CB:Z:', str(dout), stdin=text)
		assert sorted(os.listdir(str(dout))) == ['A', 'B', 'C']
		for c in 'ABC':
			assert read_lines(str(dout / c)) == expected[c]

	@pytest.mark.parametrize('buffer_size', [1, 3, 10000])
	def test_output_unknown_and_buffer_size(self, run_helper, tmp_path, buffer_size):
		"""--output_unknown collects untagged reads; output is the same for any --buffer_size."""
		text, expected = sam_input(CELLS)
		dout = tmp_path / 'out'
		funk = tmp_path / 'unknown.txt'
		run_helper('split_bam_text.py', '--section', 'CB:Z:', '--output_unknown', str(funk), '--buffer_size',
			buffer_size, str(dout), stdin=text)
		for c in 'ABC':
			assert read_lines(str(dout / c)) == expected[c]
		assert read_lines(str(funk)) == expected[None]

	def test_default_section(self, run_helper, tmp_path):
		"""The default section 'RG:' takes the barcode as the rest of a field 'RG:<barcode>' (not valid SAM; valid SAM gives a barcode starting with Z:, a bug in the base)."""
		lines = [sam_line('r1', 'RG:c1', 'NM:i:0'), sam_line('r2', 'RG:c2', 'NM:i:0'), sam_line('r3', 'RG:c1', 'NM:i:1')]
		dout = tmp_path / 'out'
		run_helper('split_bam_text.py', str(dout), stdin=''.join(x + '\n' for x in lines))
		assert read_lines(str(dout / 'c1')) == [lines[0], lines[2]]
		assert read_lines(str(dout / 'c2')) == [lines[1]]

	@pytest.mark.parametrize('buffer_size', [1, 10000])
	def test_ref_expression(self, run_helper, tmp_path, buffer_size):
		"""Barcodes absent from the --ref_expression columns count as unknown."""
		text, expected = sam_input(CELLS)
		fref = str(tmp_path / 'expression.tsv.gz')
		write_expression(fref, ['g1', 'g2'], ['A', 'C'])
		dout = tmp_path / 'out'
		funk = tmp_path / 'unknown.txt'
		run_helper('split_bam_text.py', '--section', 'CB:Z:', '--ref_expression', fref, '--output_unknown',
			str(funk), '--buffer_size', buffer_size, str(dout), stdin=text)
		assert sorted(os.listdir(str(dout))) == ['A', 'C']
		assert read_lines(str(dout / 'A')) == expected['A']
		assert sorted(read_lines(str(funk))) == sorted(expected['B'] + expected[None])

	@pytest.mark.parametrize('buffer_size', [1, 10000])
	def test_namemap(self, run_helper, tmp_path, buffer_size):
		"""--namemap renames ATAC barcodes to RNA names; unmapped barcodes count as unknown."""
		text, expected = sam_input(CELLS)
		fmap = tmp_path / 'map.tsv'
		fmap.write_text('x\trnaA\tA\nx\trnaB\tB\n')
		dout = tmp_path / 'out'
		funk = tmp_path / 'unknown.txt'
		run_helper('split_bam_text.py', '--section', 'CB:Z:', '--namemap', '{},1,2'.format(fmap), '--output_unknown',
			str(funk), '--buffer_size', buffer_size, str(dout), stdin=text)
		assert sorted(os.listdir(str(dout))) == ['rnaA', 'rnaB']
		assert read_lines(str(dout / 'rnaA')) == expected['A']
		assert read_lines(str(dout / 'rnaB')) == expected['B']
		assert sorted(read_lines(str(funk))) == sorted(expected['C'] + expected[None])

	def test_namemap_with_ref_expression(self, run_helper, tmp_path):
		"""With both options, a barcode is kept only when its mapped RNA name is in the reference."""
		text, expected = sam_input(CELLS)
		fmap = tmp_path / 'map.tsv'
		fmap.write_text('rnaA\tA\nrnaB\tB\nrnaC\tC\n')
		fref = str(tmp_path / 'expression.tsv.gz')
		write_expression(fref, ['g1'], ['rnaA', 'rnaC'])
		dout = tmp_path / 'out'
		funk = tmp_path / 'unknown.txt'
		run_helper('split_bam_text.py', '--section', 'CB:Z:', '--namemap', '{},0,1'.format(fmap), '--ref_expression',
			fref, '--output_unknown', str(funk), str(dout), stdin=text)
		assert sorted(os.listdir(str(dout))) == ['rnaA', 'rnaC']
		assert read_lines(str(dout / 'rnaC')) == expected['C']
		assert sorted(read_lines(str(funk))) == sorted(expected['B'] + expected[None])

	def test_namemap_needs_three_parts(self, run_helper, tmp_path):
		"""--namemap must be file,RNA_column,ATAC_column."""
		with pytest.raises(AssertionError):
			run_helper('split_bam_text.py', '--namemap', 'a.tsv,1', str(tmp_path / 'out'), stdin='')

	def test_existing_outputs_replaced(self, run_helper, tmp_path):
		"""Old files in the output folder and an old unknown file are removed before splitting."""
		dout = tmp_path / 'out'
		dout.mkdir()
		(dout / 'stale').write_text('x\n')
		funk = tmp_path / 'unknown.txt'
		funk.write_text('old\n')
		text, expected = sam_input(['A', None])
		run_helper('split_bam_text.py', '--section', 'CB:Z:', '--output_unknown', str(funk), str(dout), stdin=text)
		assert os.listdir(str(dout)) == ['A']
		assert read_lines(str(funk)) == expected[None]

	def test_multiple_matching_fields(self, run_helper, tmp_path):
		"""A read with two fields matching --section raises RuntimeError and removes all outputs."""
		dout = tmp_path / 'out'
		funk = tmp_path / 'unknown.txt'
		# NM:i:0 last, so the barcode does not end in the line's newline (otherwise the barcode keeps the newline, a bug in the base)
		text = sam_line('r1', 'CB:Z:A', 'NM:i:0') + '\n' + sam_line('r2', 'CB:Z:A', 'CB:Z:B', 'NM:i:0') + '\n'
		with pytest.raises(RuntimeError, match='multiple fields'):
			run_helper('split_bam_text.py', '--section', 'CB:Z:', '--output_unknown', str(funk), '--buffer_size', 1,
				str(dout), stdin=text)
		assert not dout.exists()
		assert not funk.exists()

	def test_missing_parent_folder(self, run_helper, tmp_path):
		"""An output folder whose parent does not exist raises FileNotFoundError."""
		with pytest.raises(FileNotFoundError):
			run_helper('split_bam_text.py', str(tmp_path / 'no' / 'out'), stdin='')

	def test_empty_input(self, run_helper, tmp_path):
		"""Empty input creates an empty output folder."""
		dout = tmp_path / 'out'
		run_helper('split_bam_text.py', str(dout), stdin='')
		assert os.listdir(str(dout)) == []


class TestSplitBamTextCleanup:
	@pytest.fixture
	def cleanup(self, run_helper, tmp_path):
		"""The script's module-level cleanup function, taken from one run on empty input."""
		return run_helper('split_bam_text.py', str(tmp_path / 'first'), stdin='')['cleanup']

	def test_removes_folder_and_file(self, cleanup, tmp_path):
		"""cleanup removes the output folder with its contents and the unknown-reads file."""
		d = tmp_path / 'd'
		(d / 'sub').mkdir(parents=True)
		(d / 'sub' / 'x').write_text('1')
		f = tmp_path / 'u.txt'
		f.write_text('1')
		cleanup(str(d), str(f))
		assert not d.exists()
		assert not f.exists()

	def test_missing_paths_ignored(self, cleanup, tmp_path):
		"""cleanup ignores a missing folder and a missing file."""
		cleanup(str(tmp_path / 'none'), str(tmp_path / 'none.txt'))
		assert sorted(os.listdir(str(tmp_path))) == ['first']

	def test_no_unknown_file(self, cleanup, tmp_path):
		"""With funknown=None only the folder is removed."""
		d = tmp_path / 'd'
		d.mkdir()
		keep = tmp_path / 'keep.txt'
		keep.write_text('1')
		cleanup(str(d), None)
		assert not d.exists()
		assert keep.exists()
