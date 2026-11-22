"""
Tests for dictys.chromatin (macs2, wellington, homer, binding, tssdist,
linking, binlinking), the scripts macs2, wellington and homer run, and the
``dictys chromatin ...`` subcommands.

macs2, wellington and homer call bash scripts that call samtools, macs2,
wellington_footprints.py and findMotifsGenome.pl. Their tests put small
stand-in programs first on PATH and check what the script and dictys do with
the stand-ins' output. These tests need bash and a POSIX shell.
"""

import logging
import os
import sys
import tempfile

import numpy as np
import pandas as pd
import pytest

from dictys import chromatin


def _write_text(path, text):
	with open(path, 'w', newline='\n') as f:
		f.write(text)
	return str(path)


def _write_rows(path, rows):
	return _write_text(path, ''.join('\t'.join(str(y) for y in x) + '\n' for x in rows))


def _read_bed(path):
	with open(path) as f:
		return [x.split('\t') for x in f.read().splitlines() if x]


def _write_exp(path, genes):
	"""Expression matrix with the given genes as rows and two cells."""
	pd.DataFrame(np.ones((len(genes), 2)), index=genes, columns=['c1', 'c2']).to_csv(path, sep='\t')
	return str(path)


################################################################
# Stand-in programs
################################################################

class _Stubs:
	"""Writes stand-in programs into a folder that is first on PATH; each logs its argv to $STUB_LOG."""
	def __init__(self, bindir, log):
		self.bindir = bindir
		self.log = log

	def add(self, name, body):
		path = os.path.join(self.bindir, name)
		_write_text(path, '#!/bin/sh\necho "' + name + ' $*" >> "$STUB_LOG"\n' + body)
		os.chmod(path, 0o755)

	def calls(self, name):
		"""Argument lists of every logged call of program ``name``."""
		if not os.path.exists(self.log):
			return []
		with open(self.log) as f:
			lines = f.read().splitlines()
		return [x.split()[1:] for x in lines if x.split()[:1] == [name]]


@pytest.fixture
def stubs(tmp_path, monkeypatch, require):
	"""A _Stubs folder first on PATH. Skips without bash or on Windows, where dictys cannot start its scripts."""
	if sys.platform == 'win32':
		pytest.skip('dictys starts its bash scripts through a POSIX shell')
	require('bash')
	bindir = tmp_path / 'stubbin'
	bindir.mkdir()
	monkeypatch.setenv('PATH', str(bindir) + os.pathsep + os.environ.get('PATH', ''))
	monkeypatch.setenv('STUB_LOG', str(tmp_path / 'stub.log'))
	return _Stubs(str(bindir), str(tmp_path / 'stub.log'))


# samtools stand-in: "BAM" files are SAM text. view/sort copy their input to -o or stdout; index writes the index file.
SAMTOOLS_STUB = '''sub=$1; shift
out=; a=; b=
while [ $# -gt 0 ]; do
	case "$1" in
		-o) out=$2; shift;;
		-@|-l) shift;;
		-?*) ;;
		*) if [ -z "$a" ]; then a=$1; else b=$1; fi;;
	esac
	shift
done
case "$sub" in
	index) echo index > "$b";;
	*) if [ -n "$out" ]; then cat "$a" > "$out"; else cat "$a"; fi;;
esac
'''


################################################################
# macs2
################################################################

SAM_HEADER = ['@HD\tVN:1.6', '@SQ\tSN:chr1\tLN:1000', '@PG\tID:bwa']
SAM_READ = '\t0\t{}\t{}\t60\t4M\t*\t0\t0\tACGT\tIIII'
SAM_CELL1 = SAM_HEADER + ['r1' + SAM_READ.format('chr1', 100), 'r2' + SAM_READ.format('chrM', 100)]
SAM_CELL2 = ['r3' + SAM_READ.format('chr1', 200), 'r4' + SAM_READ.format('chrUn_KI270302v1', 5),
	'r5' + SAM_READ.format('chr1_KI270706v1_random', 5), 'r6' + SAM_READ.format('GL000220.1', 5)]
# narrowPeak rows; column 8 (0-based) holds 5, 9, 1, 7, 3, 8
PEAKS = [['chr1', 10 + 20 * i, 20 + 20 * i, 'p' + str(i + 1), 50, '.', 2.0, 6.0, q, 5] for i, q in enumerate([5, 9, 1, 7, 3, 8])]


@pytest.fixture
def macs2_files(tmp_path):
	"""Names file, folder of two cells' SAM-text "BAM" files, and output paths for macs2."""
	bamdir = tmp_path / 'bams'
	bamdir.mkdir()
	_write_text(bamdir / 'cell1.bam', '\n'.join(SAM_CELL1) + '\n')
	_write_text(bamdir / 'cell2.bam', '\n'.join(SAM_CELL2) + '\n')
	names = _write_text(tmp_path / 'names.txt', 'cell1\n\ncell2\n')
	out = tmp_path / 'out'
	out.mkdir()
	return dict(fi_names=names, fi_bam=str(bamdir), fo_bam=str(out / 'all.bam'), fo_bai=str(out / 'all.bai'),
		fo_bed=str(out / 'peaks.bed'))


@pytest.fixture
def macs2_stubs(stubs, tmp_path, monkeypatch):
	"""samtools stand-in plus a macs2 stand-in that writes the PEAKS rows as its narrowPeak output."""
	stubs.add('samtools', SAMTOOLS_STUB)
	stubs.add('macs2', 'cp "$STUB_PEAKS" 04_peaks.narrowPeak\n')
	monkeypatch.setenv('STUB_PEAKS', _write_rows(tmp_path / 'peaks.narrowPeak', PEAKS))
	return stubs


def _macs2_args(files):
	return [files[x] for x in ['fi_names', 'fi_bam', 'fo_bam', 'fo_bai', 'fo_bed']]


def test_macs2_merges_and_filters_reads(macs2_files, macs2_stubs):
	"""The merged BAM keeps cell1's @SQ header and the chr1 reads of both cells, dropping @HD, @PG, chrM, chrUn_, *random and GL00 reads."""
	chromatin.macs2(*_macs2_args(macs2_files), 'mm', nmax=0)
	with open(macs2_files['fo_bam']) as f:
		assert f.read().splitlines() == [SAM_HEADER[1], SAM_CELL1[3], SAM_CELL2[0]]
	assert os.path.isfile(macs2_files['fo_bai'])


def test_macs2_passes_arguments(macs2_files, macs2_stubs):
	"""genome_size and qcut reach macs2 callpeak (-g, -q) on the merged BAM; nth=3 reaches samtools as -@ 3."""
	chromatin.macs2(*_macs2_args(macs2_files), 'hs', qcut=0.01, nth=3, nmax=0)
	(callpeak,) = macs2_stubs.calls('macs2')
	assert callpeak[0] == 'callpeak'
	for flag, value in [('-t', macs2_files['fo_bam']), ('-g', 'hs'), ('-q', '0.01')]:
		assert callpeak[callpeak.index(flag) + 1] == value
	calls = macs2_stubs.calls('samtools')
	assert {x[0] for x in calls} == {'view', 'sort', 'index'}
	assert all(x[x.index('-@') + 1] == '3' for x in calls)


def test_macs2_single_thread_omits_samtools_threads(macs2_files, macs2_stubs):
	"""With nth=1 no samtools call gets a -@ option."""
	chromatin.macs2(*_macs2_args(macs2_files), 'hs', nth=1, nmax=0)
	calls = macs2_stubs.calls('samtools')
	assert len(calls) >= 4
	assert not any('-@' in x for x in calls)


@pytest.mark.parametrize('nmax', [0, 6, 100])
def test_macs2_nmax_keeps_all(macs2_files, macs2_stubs, tmp_path, nmax):
	"""nmax=0, or nmax at least the peak count (6), leaves macs2's peak file unchanged."""
	chromatin.macs2(*_macs2_args(macs2_files), 'hs', nmax=nmax)
	with open(macs2_files['fo_bed']) as f, open(tmp_path / 'peaks.narrowPeak') as g:
		assert f.read() == g.read()


@pytest.mark.parametrize('nmax,names', [(3, ['p2', 'p4', 'p6']), (1, ['p2']), (5, ['p1', 'p2', 'p4', 'p5', 'p6'])])
def test_macs2_nmax_truncates_by_column_8(macs2_files, macs2_stubs, nmax, names):
	"""nmax keeps the peaks with the nmax largest column-8 values (5, 9, 1, 7, 3, 8), in their original order."""
	chromatin.macs2(*_macs2_args(macs2_files), 'hs', nmax=nmax)
	d = _read_bed(macs2_files['fo_bed'])
	assert [x[3] for x in d] == names
	expected = {x[3]: x[8] for x in PEAKS}
	assert [float(x[8]) for x in d] == [expected[x] for x in names]


def test_macs2_program_failure(macs2_files, stubs):
	"""A failing macs2 program makes macs2 raise RuntimeError."""
	stubs.add('samtools', SAMTOOLS_STUB)
	stubs.add('macs2', 'exit 1\n')
	with pytest.raises(RuntimeError):
		chromatin.macs2(*_macs2_args(macs2_files), 'hs')


def test_macs2_cli(macs2_files, macs2_stubs, run_dictys):
	"""``dictys chromatin macs2`` with --qcut, --nth and --nmax passes them on and truncates to the top 3 peaks."""
	run_dictys('chromatin', 'macs2', *_macs2_args(macs2_files), 'hs', '--qcut', 0.2, '--nth', 2, '--nmax', 3)
	assert [x[3] for x in _read_bed(macs2_files['fo_bed'])] == ['p2', 'p4', 'p6']
	(callpeak,) = macs2_stubs.calls('macs2')
	assert callpeak[callpeak.index('-q') + 1] == '0.2'
	assert all(x[x.index('-@') + 1] == '2' for x in macs2_stubs.calls('samtools'))


@pytest.mark.parametrize('qcut', [0, 1, -0.1, 1.5])
def test_macs2_qcut_out_of_range(macs2_files, qcut):
	"""qcut outside (0, 1) raises ValueError."""
	with pytest.raises(ValueError):
		chromatin.macs2(*_macs2_args(macs2_files), 'hs', qcut=qcut)


def test_macs2_missing_names_file(macs2_files, tmp_path):
	"""A missing fi_names raises FileNotFoundError."""
	macs2_files['fi_names'] = str(tmp_path / 'missing.txt')
	with pytest.raises(FileNotFoundError):
		chromatin.macs2(*_macs2_args(macs2_files), 'hs')


@pytest.mark.parametrize('kind', ['missing', 'file'])
def test_macs2_bam_folder_not_a_folder(macs2_files, tmp_path, kind):
	"""fi_bam that is missing or is a file rather than a folder raises FileNotFoundError."""
	macs2_files['fi_bam'] = str(tmp_path / 'missing') if kind == 'missing' else macs2_files['fi_names']
	with pytest.raises(FileNotFoundError):
		chromatin.macs2(*_macs2_args(macs2_files), 'hs')


def test_macs2_no_names(macs2_files):
	"""A names file with only blank lines raises ValueError."""
	_write_text(macs2_files['fi_names'], '\n  \n')
	with pytest.raises(ValueError):
		chromatin.macs2(*_macs2_args(macs2_files), 'hs')


def test_macs2_cli_validates(macs2_files, run_dictys):
	"""``dictys chromatin macs2 --qcut 1.5`` raises ValueError before any program runs."""
	with pytest.raises(ValueError):
		run_dictys('chromatin', 'macs2', *_macs2_args(macs2_files), 'hs', '--qcut', 1.5, '--nth', 2, '--nmax', 10)


################################################################
# wellington
################################################################

# Footprints the wellington stand-in reports: chr, start, end, name, log p-value
FOOTPRINTS = [
	['chr1', 110, 120, 'fp1', -12],
	['chr1', 130, 140, 'fp2', -8],
	['chr1', 150, 160, 'fp3', -30],
	['chr2', 10, 20, 'fp4', -15],
	['chr2', 30, 40, 'fp5', -11],
]
FOOTPRINTS_OUT = {x[3]: [x[0], str(x[1]), str(x[2]), '{}:{}:{}:{}:{}'.format(x[3], *x[:3], x[4]), str(x[4])] for x in FOOTPRINTS}


@pytest.fixture
def wellington_files(tmp_path):
	"""Dummy BAM and BAI files, a peak bed with two overlapping chr1 peaks, a chr2 blacklist and an output path."""
	return dict(
		fi_bam=_write_text(tmp_path / 'all.bam', 'bam'),
		fi_bai=_write_text(tmp_path / 'all.bai', 'bai'),
		fi_bed=_write_rows(tmp_path / 'peaks.bed', [['chr1', 100, 200], ['chr1', 150, 250], ['chr2', 0, 50]]),
		fo_bed=str(tmp_path / 'footprints.bed'),
		blacklist=_write_rows(tmp_path / 'blacklist.bed', [['chr2', 0, 100]]))


@pytest.fixture
def wellington_stubs(stubs, require, tmp_path, monkeypatch):
	"""Real bedtools plus a wellington_footprints.py stand-in that copies its peak input to $STUB_LOG.peaks and reports FOOTPRINTS."""
	require('bedtools')
	stubs.add('wellington_footprints.py', 'cp "$1" "$STUB_LOG.peaks"\nmkdir -p "$3/p value cutoffs"\ncp "$STUB_FP" "$3/p value cutoffs/fp$9.bed"\n')
	monkeypatch.setenv('STUB_FP', _write_rows(tmp_path / 'fp.bed', FOOTPRINTS))
	return stubs


def _wellington_args(files):
	return [files[x] for x in ['fi_bam', 'fi_bai', 'fi_bed', 'fo_bed']]


def test_wellington_reformats_footprints(wellington_files, wellington_stubs):
	"""With cut=10 and nmax=0, footprints below -10 are kept in order and named name:chr:start:end:score."""
	chromatin.wellington(*_wellington_args(wellington_files), cut=10, nth=1, nmax=0)
	assert _read_bed(wellington_files['fo_bed']) == [FOOTPRINTS_OUT[x] for x in ['fp1', 'fp3', 'fp4', 'fp5']]
	assert _read_bed(wellington_stubs.log + '.peaks') == [['chr1', '100', '250'], ['chr2', '0', '50']]
	(call,) = wellington_stubs.calls('wellington_footprints.py')
	assert call[1] == wellington_files['fi_bam']
	assert call[call.index('-pv') + 1] == '-10'
	assert call[call.index('-p') + 1] == '1'


def test_wellington_cut(wellington_files, wellington_stubs):
	"""cut=11.5 drops fp5 (-11) as well as fp2 (-8)."""
	chromatin.wellington(*_wellington_args(wellington_files), cut=11.5, nmax=0)
	assert _read_bed(wellington_files['fo_bed']) == [FOOTPRINTS_OUT[x] for x in ['fp1', 'fp3', 'fp4']]
	(call,) = wellington_stubs.calls('wellington_footprints.py')
	assert call[call.index('-pv') + 1] == '-11.5'


def test_wellington_nmax(wellington_files, wellington_stubs):
	"""nmax=2 keeps the two strongest footprints, fp3 (-30) then fp4 (-15)."""
	chromatin.wellington(*_wellington_args(wellington_files), nmax=2)
	assert _read_bed(wellington_files['fo_bed']) == [FOOTPRINTS_OUT['fp3'], FOOTPRINTS_OUT['fp4']]


def test_wellington_program_failure(wellington_files, stubs, require):
	"""A failing wellington_footprints.py makes wellington raise RuntimeError."""
	require('bedtools')
	stubs.add('wellington_footprints.py', 'exit 1\n')
	with pytest.raises(RuntimeError):
		chromatin.wellington(*_wellington_args(wellington_files))


def test_wellington_cli(wellington_files, wellington_stubs, run_dictys):
	"""``dictys chromatin wellington`` with every option: the blacklist removes the chr2 peak, nth=3 gives -p 2, nmax=2 keeps fp3 and fp4."""
	run_dictys('chromatin', 'wellington', *_wellington_args(wellington_files), '--fi_blacklist', wellington_files['blacklist'],
		'--cut', 11.5, '--nth', 3, '--nmax', 2)
	assert _read_bed(wellington_files['fo_bed']) == [FOOTPRINTS_OUT['fp3'], FOOTPRINTS_OUT['fp4']]
	assert _read_bed(wellington_stubs.log + '.peaks') == [['chr1', '100', '250']]
	(call,) = wellington_stubs.calls('wellington_footprints.py')
	assert call[call.index('-p') + 1] == '2'


def test_wellington_missing_blacklist(wellington_files, tmp_path):
	"""A fi_blacklist that does not exist raises FileNotFoundError."""
	with pytest.raises(FileNotFoundError):
		chromatin.wellington(*_wellington_args(wellington_files), fi_blacklist=str(tmp_path / 'missing.bed'))


@pytest.mark.parametrize('cut', [0, -1])
def test_wellington_cut_not_positive(wellington_files, cut):
	"""cut <= 0 raises ValueError."""
	with pytest.raises(ValueError):
		chromatin.wellington(*_wellington_args(wellington_files), cut=cut)


@pytest.mark.parametrize('key', ['fi_bam', 'fi_bai', 'fi_bed'])
def test_wellington_missing_input(wellington_files, tmp_path, key):
	"""A missing BAM, BAI or peak bed raises FileNotFoundError."""
	wellington_files[key] = str(tmp_path / 'missing')
	with pytest.raises(FileNotFoundError):
		chromatin.wellington(*_wellington_args(wellington_files))


def test_wellington_cli_validates(wellington_files, run_dictys):
	"""``dictys chromatin wellington --cut 0`` raises ValueError before any program runs."""
	with pytest.raises(ValueError):
		run_dictys('chromatin', 'wellington', *_wellington_args(wellington_files), '--fi_blacklist', wellington_files['blacklist'],
			'--cut', 0, '--nth', 2, '--nmax', 5)


################################################################
# homer, chromatin_homer.py and _motif_postproc
################################################################

HOMER_HEADER = ['PositionID', 'Offset', 'Sequence', 'Motif Name', 'Strand', 'MotifScore']
# Regions in the wellington output format: name column is footprint:chr:start:end:score
HOMER_REGIONS = [
	['chr1', 110, 120, 'fp1:chr1:110:120:-12', -12],
	['chr1', 150, 160, 'fp3:chr1:150:160:-30', -30],
	['chr2', 10, 20, 'fp4:chr2:10:20:-15', -15],
]
HOMER_HITS = [
	['fp1:chr1:110:120:-12', 0, 'ACGT', 'GATA1_m1', '+', 8.5],
	['fp1:chr1:110:120:-12', 2, 'ACGT', 'GATA1_m1', '-', 9.5],
	['fp1:chr1:110:120:-12', 1, 'ACGT', 'TAL1,TAL2_m2', '+', 7.0],
	['fp3:chr1:150:160:-30', 0, 'ACGT', 'FOXA1_m3', '+', 6.0],
	['fp4:chr2:10:20:-15', 0, 'ACGT', 'TAL1,TAL2_m2', '+', 5.0],
	['fp4:chr2:10:20:-15', 3, 'ACGT', 'GATA1_m1', '+', 0],
	['fpX:chr9:1:2:-50', 0, 'ACGT', 'GATA1_m1', '+', 4.0],
]
# 16-long.bed expected from the hits above: MotifScore 0 dropped, the GATA1_m1 duplicate at fp1 keeps 9.5,
# w is |score| from the region name; rows sorted by chr, start, end, motif. fpX is not in HOMER_REGIONS.
HOMER_LONG = [
	['chr1', 110, 120, 'chr1:110:120', 12.0, 'GATA1_m1', 9.5],
	['chr1', 110, 120, 'chr1:110:120', 12.0, 'TAL1,TAL2_m2', 7.0],
	['chr1', 150, 160, 'chr1:150:160', 30.0, 'FOXA1_m3', 6.0],
	['chr2', 10, 20, 'chr2:10:20', 15.0, 'TAL1,TAL2_m2', 5.0],
]
HOMER_GENES = ['GATA1', 'TAL1', 'MYC']


def _read_long(path):
	return [[x[0], int(x[1]), int(x[2]), x[3], float(x[4]), x[5], float(x[6])] for x in _read_bed(path)]


def _read_table(path):
	return pd.read_csv(path, sep='\t', header=0, index_col=0)


def _assert_table(path, index, columns, values):
	d = _read_table(path)
	assert list(d.index) == index
	assert list(d.columns) == columns
	np.testing.assert_allclose(d.values, np.array(values, dtype=float))


def _assert_homer_outputs(files):
	"""The outputs of homer for HOMER_HITS: FOXA1 is not in HOMER_GENES, so its motif and chr1:150:160 go; TAL2 is dropped from the TAL1,TAL2 name."""
	assert _read_long(files['fo_bed']) == HOMER_LONG
	_assert_table(files['fo_wellington'], ['chr1:110:120', 'chr2:10:20'], ['GATA1_m1', 'TAL1_m2'], [[12, 12], [0, 15]])
	_assert_table(files['fo_homer'], ['chr1:110:120', 'chr2:10:20'], ['GATA1_m1', 'TAL1_m2'], [[9.5, 7.0], [0, 5.0]])


def test_chromatin_homer_script(tmp_path, run_helper):
	"""chromatin_homer.py turns 15-tf.bed into 16-long.bed and region-by-motif wellington and homer tables, NaN where no hit."""
	_write_rows(tmp_path / '15-tf.bed', [HOMER_HEADER] + HOMER_HITS[:6])
	run_helper('chromatin_homer.py', folder='', cwd=tmp_path)
	assert _read_long(tmp_path / '16-long.bed') == HOMER_LONG
	index = ['chr1:110:120', 'chr1:150:160', 'chr2:10:20']
	columns = ['FOXA1_m3', 'GATA1_m1', 'TAL1,TAL2_m2']
	nan = np.nan
	_assert_table(tmp_path / '19-w.tsv', index, columns, [[nan, 12, 12], [30, nan, nan], [nan, nan, 15]])
	_assert_table(tmp_path / '19-h.tsv', index, columns, [[nan, 9.5, 7.0], [6.0, nan, nan], [nan, nan, 5.0]])


@pytest.fixture
def homer_files(tmp_path):
	"""Inputs and output paths for homer: HOMER_REGIONS, a motif file, a genome folder and an expression matrix of HOMER_GENES."""
	genome = tmp_path / 'genome'
	genome.mkdir()
	return dict(
		fi_bed=_write_rows(tmp_path / 'regions.bed', HOMER_REGIONS),
		fi_motif=_write_text(tmp_path / 'motifs.motif', '>ACGT\tGATA1_m1\t5.0\n0.97\t0.01\t0.01\t0.01\n'),
		dirio_genome=str(genome),
		fi_exp=_write_exp(tmp_path / 'exp.tsv', HOMER_GENES),
		fo_bed=str(tmp_path / 'motifs.bed'),
		fo_wellington=str(tmp_path / 'wellington.tsv'),
		fo_homer=str(tmp_path / 'homer.tsv'))


@pytest.fixture
def homer_stubs(stubs, tmp_path, monkeypatch):
	"""findMotifsGenome.pl stand-in printing the HOMER_HITS rows of its input regions; python3 runs this interpreter."""
	stubs.add('findMotifsGenome.pl', "printf 'PositionID\\tOffset\\tSequence\\tMotif Name\\tStrand\\tMotifScore\\n'\n"
		"awk -F'\\t' 'NR==FNR{k[$4]=1;next} ($1 in k)' \"$1\" \"$STUB_HITS\"\n")
	stubs.add('python3', 'exec "' + sys.executable + '" "$@"\n')
	monkeypatch.setenv('STUB_HITS', _write_rows(tmp_path / 'hits.tsv', HOMER_HITS))
	return stubs


def _homer_args(files):
	return [files[x] for x in ['fi_bed', 'fi_motif', 'dirio_genome', 'fi_exp', 'fo_bed', 'fo_wellington', 'fo_homer']]


def test_homer_outputs(homer_files, homer_stubs):
	"""homer with one thread scans all regions in one findMotifsGenome.pl call and writes the filtered motif tables."""
	chromatin.homer(*_homer_args(homer_files))
	(call,) = homer_stubs.calls('findMotifsGenome.pl')
	assert call[1] == homer_files['dirio_genome']
	assert call[3:] == ['-size', 'given', '-mask', '-find', homer_files['fi_motif']]
	_assert_homer_outputs(homer_files)


def test_homer_cli_threads(homer_files, homer_stubs, run_dictys):
	"""``dictys chromatin homer --nth 2`` splits the 3 regions into 2 chunks and merges them to the same outputs."""
	run_dictys('chromatin', 'homer', *_homer_args(homer_files), '--nth', 2)
	assert len(homer_stubs.calls('findMotifsGenome.pl')) == 2
	_assert_homer_outputs(homer_files)


def test_homer_genome_file(homer_files, homer_stubs, tmp_path):
	"""dirio_genome may be a file, which is passed to findMotifsGenome.pl as given."""
	homer_files['dirio_genome'] = _write_text(tmp_path / 'genome.fa', '>chr1\nACGT\n')
	chromatin.homer(*_homer_args(homer_files))
	(call,) = homer_stubs.calls('findMotifsGenome.pl')
	assert call[1] == homer_files['dirio_genome']
	_assert_homer_outputs(homer_files)


def test_homer_no_hits(homer_files, homer_stubs, tmp_path, monkeypatch):
	"""When findMotifsGenome.pl finds nothing, the scripts still finish and homer raises dictys's 'No motif found.'."""
	monkeypatch.setenv('STUB_HITS', _write_text(tmp_path / 'nohits.tsv', ''))
	with pytest.raises(RuntimeError, match='No motif found'):
		chromatin.homer(*_homer_args(homer_files))


@pytest.mark.parametrize('key', ['fi_bed', 'fi_motif', 'dirio_genome'])
def test_homer_missing_input(homer_files, tmp_path, key):
	"""A missing region bed, motif file or genome raises FileNotFoundError."""
	homer_files[key] = str(tmp_path / 'missing')
	with pytest.raises(FileNotFoundError):
		chromatin.homer(*_homer_args(homer_files))


def test_homer_cli_validates(homer_files, tmp_path, run_dictys):
	"""``dictys chromatin homer --nth 2`` with a missing genome raises FileNotFoundError."""
	homer_files['dirio_genome'] = str(tmp_path / 'missing')
	with pytest.raises(FileNotFoundError):
		run_dictys('chromatin', 'homer', *_homer_args(homer_files), '--nth', 2)


def _tsv_bytes(d):
	return d.to_csv(None, sep='\t', header=True, index=True).encode()


@pytest.fixture
def postproc_out(tmp_path):
	"""Output paths for _motif_postproc."""
	return [str(tmp_path / x) for x in ['motifs.bed', 'w.tsv', 'h.tsv']]


def test_motif_postproc_nan_and_names(tmp_path, postproc_out):
	"""A NaN in either table zeroes both; TFs outside fi_exp leave motif names; regions left all-zero are dropped; the bed is copied."""
	nan = np.nan
	columns = ['GATA1,GATA2,TAL1_a', 'MYC_b', 'TAL1_c_x']
	dw = pd.DataFrame([[1.0, nan, 2.0], [nan, 4.0, 3.0], [5.0, 6.0, 2.0]], index=['r1', 'r2', 'r3'], columns=columns)
	dh = pd.DataFrame([[7.0, 1.0, nan], [8.0, 9.0, nan], [nan, 2.0, 3.0]], index=['r1', 'r2', 'r3'], columns=columns)
	bed = b'chr1\t1\t2\tr1\t1.0\tGATA1,GATA2,TAL1_a\t7.0\n'
	fi_exp = _write_exp(tmp_path / 'exp.tsv', ['TAL1', 'GATA2'])
	chromatin._motif_postproc((_tsv_bytes(dw), _tsv_bytes(dh), bed), fi_exp, *postproc_out)
	with open(postproc_out[0], 'rb') as f:
		assert f.read() == bed
	# After zeroing NaN pairs: r1 [1, 0, 0], r2 [0, 4, 0], r3 [0, 6, 2] (wellington). MYC_b goes (MYC not in
	# fi_exp), which leaves r2 all zero, so it goes too.
	_assert_table(postproc_out[1], ['r1', 'r3'], ['GATA2,TAL1_a', 'TAL1_c_x'], [[1.0, 0], [0, 2.0]])
	_assert_table(postproc_out[2], ['r1', 'r3'], ['GATA2,TAL1_a', 'TAL1_c_x'], [[7.0, 0], [0, 3.0]])


def test_motif_postproc_keeps_region_with_homer_score_only(tmp_path, postproc_out):
	"""A region is kept when only its homer score is positive."""
	dw = pd.DataFrame([[0.0], [2.0]], index=['r1', 'r2'], columns=['TAL1_a'])
	dh = pd.DataFrame([[3.0], [0.0]], index=['r1', 'r2'], columns=['TAL1_a'])
	fi_exp = _write_exp(tmp_path / 'exp.tsv', ['TAL1'])
	chromatin._motif_postproc((_tsv_bytes(dw), _tsv_bytes(dh), b'x\t1\t2\ty\t1\tz\t1\n'), fi_exp, *postproc_out)
	_assert_table(postproc_out[1], ['r1', 'r2'], ['TAL1_a'], [[0.0], [2.0]])
	_assert_table(postproc_out[2], ['r1', 'r2'], ['TAL1_a'], [[3.0], [0.0]])


def test_motif_postproc_duplicate_names(tmp_path, postproc_out):
	"""Motifs that become identical after dropping TFs outside fi_exp raise ValueError."""
	d = pd.DataFrame([[1.0, 2.0]], index=['r1'], columns=['GATA1_m', 'GATA1,GATA9_m'])
	fi_exp = _write_exp(tmp_path / 'exp.tsv', ['GATA1'])
	with pytest.raises(ValueError):
		chromatin._motif_postproc((_tsv_bytes(d), _tsv_bytes(d), b'x\t1\t2\ty\t1\tz\t1\n'), fi_exp, *postproc_out)


@pytest.mark.parametrize('dret', [None, (b'a', b'b')])
def test_motif_postproc_failed_call(tmp_path, postproc_out, dret):
	"""A missing or incomplete homer result raises RuntimeError."""
	fi_exp = _write_exp(tmp_path / 'exp.tsv', ['GATA1'])
	with pytest.raises(RuntimeError):
		chromatin._motif_postproc(dret, fi_exp, *postproc_out)


def test_motif_postproc_no_motif(tmp_path, postproc_out):
	"""An empty motif bed raises RuntimeError."""
	d = pd.DataFrame([[1.0]], index=['r1'], columns=['GATA1_m'])
	fi_exp = _write_exp(tmp_path / 'exp.tsv', ['GATA1'])
	with pytest.raises(RuntimeError):
		chromatin._motif_postproc((_tsv_bytes(d), _tsv_bytes(d), b''), fi_exp, *postproc_out)


################################################################
# binding
################################################################

BIND_REGIONS = ['chr1:100:200', 'chr1:300:400', 'chr2:50:80']
BIND_MOTIFS = ['GATA1_m1', 'GATA1_m2', 'TAL1_m1', 'GATA1,GATA2_m3']
BIND_W = [[2, 4, 0, 1], [0, 3, 5, 0], [6, 0, 2, 3]]
BIND_H = [[3, 1, 2, 0], [1, 2, 4, 0], [2, 0, 0, 5]]
R1, R2, R3 = BIND_REGIONS
L = np.log


@pytest.fixture
def binding_files(tmp_path):
	"""Wellington and homer score tables for three regions and four motifs (BIND_W, BIND_H)."""
	fw = tmp_path / 'w.tsv'
	fh = tmp_path / 'h.tsv'
	pd.DataFrame(BIND_W, index=BIND_REGIONS, columns=BIND_MOTIFS).to_csv(fw, sep='\t')
	pd.DataFrame(BIND_H, index=BIND_REGIONS, columns=BIND_MOTIFS).to_csv(fh, sep='\t')
	return str(fw), str(fh), str(tmp_path / 'bind.tsv')


def _read_binding(path):
	d = pd.read_csv(path, sep='\t', header=0, index_col=None)
	assert list(d.columns) == ['TF', 'loc', 'score']
	return {(x, y): z for x, y, z in zip(d['TF'], d['loc'], d['score'])}


def _assert_binding(path, expected):
	ans = _read_binding(path)
	assert set(ans) == set(expected)
	keys = sorted(expected)
	np.testing.assert_allclose([ans[x] for x in keys], [expected[x] for x in keys], rtol=1e-12)


# Pairs with both scores > 0: GATA1 r1 (m1 2*3, m2 4*1), r2 (m2 3*2), r3 (m1 6*2); TAL1 r2 (5*4); GATA1,GATA2 r3 (3*5).
@pytest.mark.parametrize('ka,expected', [
	(dict(), {('GATA1', R1): L(6), ('GATA1', R2): L(6), ('GATA1', R3): L(12), ('TAL1', R2): L(20), ('GATA1,GATA2', R3): L(15)}),
	(dict(combine='sum'), {('GATA1', R1): L(24), ('GATA1', R2): L(6), ('GATA1', R3): L(12), ('TAL1', R2): L(20), ('GATA1,GATA2', R3): L(15)}),
	(dict(combine='mean'), {('GATA1', R1): L(24) / 2, ('GATA1', R2): L(6), ('GATA1', R3): L(12), ('TAL1', R2): L(20), ('GATA1,GATA2', R3): L(15)}),
	(dict(mode=1), {('GATA1', R1): L(4), ('GATA1', R2): L(3), ('GATA1', R3): L(6), ('TAL1', R2): L(5), ('GATA1,GATA2', R3): L(3)}),
	(dict(mode=2), {('GATA1', R1): L(3), ('GATA1', R2): L(2), ('GATA1', R3): L(2), ('TAL1', R2): L(4), ('GATA1,GATA2', R3): L(5)}),
	# wellington > 2.5 removes GATA1_m1 at r1
	(dict(cutw=2.5), {('GATA1', R1): L(4), ('GATA1', R2): L(6), ('GATA1', R3): L(12), ('TAL1', R2): L(20), ('GATA1,GATA2', R3): L(15)}),
	# homer > 2.5 keeps GATA1_m1 at r1, TAL1 at r2 and GATA1,GATA2_m3 at r3
	(dict(cuth=2.5), {('GATA1', R1): L(6), ('TAL1', R2): L(20), ('GATA1,GATA2', R3): L(15)}),
	# score > log(10) keeps log(12), log(20) and log(15)
	(dict(cut=float(L(10))), {('GATA1', R3): L(12), ('TAL1', R2): L(20), ('GATA1,GATA2', R3): L(15)}),
	# cut applies to each motif before combining: GATA1_m2 at r1 (log 4 < 1.5) is not added to the sum
	(dict(cut=1.5, combine='sum'), {('GATA1', R1): L(6), ('GATA1', R2): L(6), ('GATA1', R3): L(12), ('TAL1', R2): L(20), ('GATA1,GATA2', R3): L(15)}),
])
def test_binding_scores(binding_files, ka, expected):
	"""Binding scores per TF and region, worked out by hand from BIND_W and BIND_H for each option."""
	chromatin.binding(*binding_files, **ka)
	_assert_binding(binding_files[2], expected)


def test_binding_nothing_left(binding_files):
	"""A wellington cutoff above every score raises RuntimeError."""
	with pytest.raises(RuntimeError):
		chromatin.binding(*binding_files, cutw=100)


def test_binding_unknown_combine(binding_files):
	"""An unknown combine method raises ValueError."""
	with pytest.raises(ValueError):
		chromatin.binding(*binding_files, combine='median')


def test_binding_cli(binding_files, run_dictys):
	"""``dictys chromatin binding`` with --cuth 0.5 --cutw 1.5 --cut 1 --combine sum --mode 1, worked out by hand."""
	run_dictys('chromatin', 'binding', *binding_files, '--cuth', 0.5, '--cutw', 1.5, '--cut', 1.0, '--combine', 'sum', '--mode', 1)
	# GATA1_m1 at r1 (log 2 < 1) is cut; TAL1 at r3 fails homer > 0.5.
	_assert_binding(binding_files[2], {('GATA1', R1): L(4), ('GATA1', R2): L(3), ('GATA1', R3): L(6), ('TAL1', R2): L(5), ('GATA1,GATA2', R3): L(3)})


################################################################
# tssdist
################################################################

TSS_REGIONS = ['chr1:1500:1600', 'chr1:300:400', 'chr1:900:1100', 'chr1:2500:2600', 'chr2:5300:5400', 'chr2:4200:4300', 'chr3:100:200']
TSS_BED = [['chr1', 1000, 3000, 'G1', 0, '+'], ['chr2', 4000, 5000, 'G2', 0, '-'], ['chr3', 100, 200, 'GX', 0, '+']]
# Within 1000 of G1 (TSS 1000, +) and G2 (TSS 5000, -). Positive is upstream, 0 overlaps the TSS.
TSS_EXPECTED = {('chr1:1500:1600', 'G1'): -500, ('chr1:300:400', 'G1'): 600, ('chr1:900:1100', 'G1'): 0,
	('chr2:5300:5400', 'G2'): 300, ('chr2:4200:4300', 'G2'): -700}


@pytest.fixture
def tss_env(require, tmp_path, monkeypatch):
	"""Skips without bedtools or pybedtools; keeps pybedtools' temporary files in tmp_path."""
	require('bedtools')
	pytest.importorskip('pybedtools')
	monkeypatch.setattr(tempfile, 'tempdir', str(tmp_path))


def _tss_files(tmp_path, regions=TSS_REGIONS, bed=TSS_BED, genes=('G1', 'G2', 'G3')):
	fw = tmp_path / 'w.tsv'
	pd.DataFrame(np.ones((len(regions), 1)), index=regions, columns=['m']).to_csv(fw, sep='\t')
	return _write_exp(tmp_path / 'exp.tsv', list(genes)), str(fw), _write_rows(tmp_path / 'tss.bed', bed), str(tmp_path / 'dist.tsv')


def _read_dist(path):
	d = pd.read_csv(path, sep='\t', header=0, index_col=None)
	assert list(d.columns) == ['region', 'target', 'dist']
	ans = {(x, y): z for x, y, z in zip(d['region'], d['target'], d['dist'])}
	assert len(ans) == len(d)
	return ans


def test_tssdist_signed(tss_env, tmp_path):
	"""Regions within 1000 bp of a TSS get signed distances; GX (not in fi_exp) and far regions are left out."""
	files = _tss_files(tmp_path)
	chromatin.tssdist(*files, cut=1000)
	assert _read_dist(files[3]) == TSS_EXPECTED


def test_tssdist_cut(tss_env, tmp_path):
	"""cut=550 drops the regions 600 and 700 bp away."""
	files = _tss_files(tmp_path)
	chromatin.tssdist(*files, cut=550)
	assert _read_dist(files[3]) == {x: y for x, y in TSS_EXPECTED.items() if abs(y) < 550}


@pytest.mark.parametrize('ncol', [4, 5])
def test_tssdist_unsigned(tss_env, tmp_path, ncol):
	"""A bed without strand puts every TSS at the start: G2's TSS is 4000, so chr2:4200:4300 is 200 bp downstream and chr2:5300:5400 is too far."""
	files = _tss_files(tmp_path, bed=[x[:ncol] for x in TSS_BED])
	chromatin.tssdist(*files, cut=1000)
	expected = {x: y for x, y in TSS_EXPECTED.items() if x[1] == 'G1'}
	expected['chr2:4200:4300', 'G2'] = -200
	assert _read_dist(files[3]) == expected


def test_tssdist_nmin_nmax(tss_env, tmp_path):
	"""nmin and nmax bound the number of relations (5 here): equal bounds pass, tighter bounds raise RuntimeError."""
	files = _tss_files(tmp_path)
	chromatin.tssdist(*files, cut=1000, nmin=5, nmax=5)
	assert len(_read_dist(files[3])) == 5
	with pytest.raises(RuntimeError):
		chromatin.tssdist(*files, cut=1000, nmin=6)
	with pytest.raises(RuntimeError):
		chromatin.tssdist(*files, cut=1000, nmax=4)


def test_tssdist_nmin_above_nmax(tss_env, tmp_path):
	"""nmin greater than nmax raises ValueError."""
	with pytest.raises(ValueError):
		chromatin.tssdist(*_tss_files(tmp_path), nmin=3, nmax=2)


def test_tssdist_no_relation(tss_env, tmp_path):
	"""Regions on a chromosome with no gene raise RuntimeError."""
	files = _tss_files(tmp_path, regions=['chr5:100:200', 'chr5:300:400'])
	with pytest.raises(RuntimeError):
		chromatin.tssdist(*files)


@pytest.mark.parametrize('nplus,nminus,ncol,warns', [
	(100, 0, 6, True),
	(99, 0, 6, False),
	(100, 0, 4, False),
	(100, 100, 6, False),
])
def test_tssdist_strand_warning(tss_env, tmp_path, caplog, nplus, nminus, ncol, warns):
	"""A strand-signed bed with at least 100 relations, over 90% on one strand, logs a warning."""
	regions = ['chr1:{}:{}'.format(1000 + 10 * i, 1005 + 10 * i) for i in range(nplus)]
	regions += ['chr2:{}:{}'.format(1000 + 10 * i, 1005 + 10 * i) for i in range(nminus)]
	bed = [['chr1', 1000, 3000, 'G1', 0, '+'][:ncol], ['chr2', 1000, 3000, 'G2', 0, '-'][:ncol]]
	files = _tss_files(tmp_path, regions=regions, bed=bed)
	with caplog.at_level(logging.WARNING):
		chromatin.tssdist(*files)
	assert len(_read_dist(files[3])) == nplus + nminus
	assert any(x.getMessage() == 'Heavily unbalanced strand found.' for x in caplog.records) == warns


def test_tssdist_cli(tss_env, tmp_path, run_dictys):
	"""``dictys chromatin tssdist`` with --cut 550 --nmin 3 --nmax 3 writes the three relations within 550 bp."""
	files = _tss_files(tmp_path)
	run_dictys('chromatin', 'tssdist', *files, '--cut', 550, '--nmin', 3, '--nmax', 3)
	assert _read_dist(files[3]) == {x: y for x, y in TSS_EXPECTED.items() if abs(y) < 550}


################################################################
# linking
################################################################

LINK_BINDING = [['A', 'chr1:100:200', 2.0], ['A', 'chr1:300:400', 1.0], ['B', 'chr1:300:400', 3.0], ['A,C', 'chr2:100:200', 0.5]]
LINK_DIST = [['chr1:100:200', 'G1', 1000], ['chr1:300:400', 'G1', -2000], ['chr1:300:400', 'B', 500],
	['chr2:100:200', 'G2', 0], ['chr3:1:2', 'G3', 10]]
# Distance penalty per bp for mode 4
P = np.log(10) / 1E6
NINF = -np.inf


@pytest.fixture
def linking_files(tmp_path):
	"""Binding table (LINK_BINDING), distance table (LINK_DIST) and output path."""
	fb = tmp_path / 'binding.tsv'
	fd = tmp_path / 'dist.tsv'
	pd.DataFrame(LINK_BINDING, columns=['TF', 'loc', 'score']).to_csv(fb, sep='\t', index=False)
	pd.DataFrame(LINK_DIST, columns=['region', 'target', 'dist']).to_csv(fd, sep='\t', index=False)
	return str(fb), str(fd), str(tmp_path / 'linking.tsv')


def _assert_linking(path, expected):
	"""``expected`` maps (TF, target) to the finite scores; every other cell must be -inf."""
	d = _read_table(path)
	tfs = sorted({x[0] for x in expected})
	assert sorted(d.index) == tfs
	# Columns are the targets with a score plus every TF
	assert sorted(d.columns) == sorted({x[1] for x in expected} | set(tfs))
	d = d.loc[tfs, sorted(d.columns)]
	want = np.full(d.shape, NINF)
	for (tf, target), value in expected.items():
		want[tfs.index(tf), list(d.columns).index(target)] = value
	np.testing.assert_allclose(d.values, want, rtol=1e-12)


# TF -> target through a shared region: binding score minus the distance penalty. B -> B is a self link and goes.
@pytest.mark.parametrize('ka,expected', [
	(dict(), {('A', 'B'): 1 - 500 * P, ('A', 'G1'): 2 - 1000 * P, ('A', 'G2'): 0.5, ('B', 'G1'): 3 - 2000 * P, ('C', 'G2'): 0.5}),
	(dict(mode=7), {('A', 'B'): 1 - 500 * P, ('A', 'G1'): 2 - 1000 * P, ('A', 'G2'): 0.5, ('B', 'G1'): 3 - 2000 * P, ('C', 'G2'): 0.5}),
	(dict(combine='sum'), {('A', 'B'): 1 - 500 * P, ('A', 'G1'): 3 - 3000 * P, ('A', 'G2'): 0.5, ('B', 'G1'): 3 - 2000 * P, ('C', 'G2'): 0.5}),
	(dict(combine='mean'), {('A', 'B'): 1 - 500 * P, ('A', 'G1'): (3 - 3000 * P) / 2, ('A', 'G2'): 0.5, ('B', 'G1'): 3 - 2000 * P, ('C', 'G2'): 0.5}),
	(dict(mode=0), {('A', 'B'): 1.0, ('A', 'G1'): 2.0, ('A', 'G2'): 0.5, ('B', 'G1'): 3.0, ('C', 'G2'): 0.5}),
	(dict(mode=3, combine='sum'), {('A', 'B'): 1.0, ('A', 'G1'): 3.0, ('A', 'G2'): 0.5, ('B', 'G1'): 3.0, ('C', 'G2'): 0.5}),
])
def test_linking_scores(linking_files, ka, expected):
	"""TF-target scores chain binding and distance scores, worked out by hand for each combine and mode."""
	chromatin.linking(*linking_files, **ka)
	_assert_linking(linking_files[2], expected)


@pytest.fixture
def whitelist(tmp_path):
	"""Whitelist: chr1 150-250 for G1 (named G1_x) and chr1 0-1000 for B."""
	return _write_rows(tmp_path / 'whitelist.bed', [['chr1', 150, 250, 'G1_x'], ['chr1', 0, 1000, 'B']])


def test_linking_whitelist_intersect(linking_files, whitelist):
	"""intersect keeps chr1:100:200 for G1 (overlaps 150-250) and chr1:300:400 for B, dropping everything else."""
	chromatin.linking(*linking_files, fi_whitelist=whitelist)
	_assert_linking(linking_files[2], {('A', 'B'): 1 - 500 * P, ('A', 'G1'): 2 - 1000 * P})


def test_linking_whitelist_within(linking_files, whitelist):
	"""within drops chr1:100:200 for G1, which only partly lies in 150-250."""
	chromatin.linking(*linking_files, fi_whitelist=whitelist, whitelist_mode='within')
	_assert_linking(linking_files[2], {('A', 'B'): 1 - 500 * P})


def test_linking_whitelist_within_inside(linking_files, tmp_path):
	"""within keeps chr1:100:200 for G1 when the whitelist interval 50-250 contains it."""
	wl = _write_rows(tmp_path / 'wl.bed', [['chr1', 50, 250, 'G1_y']])
	chromatin.linking(*linking_files, fi_whitelist=wl, whitelist_mode='within')
	_assert_linking(linking_files[2], {('A', 'G1'): 2 - 1000 * P})


def test_linking_whitelist_removes_all(linking_files, tmp_path):
	"""A whitelist matching no target's regions raises RuntimeError."""
	wl = _write_rows(tmp_path / 'wl.bed', [['chr1', 5000, 6000, 'G1'], ['chr9', 0, 1000, 'G2']])
	with pytest.raises(RuntimeError):
		chromatin.linking(*linking_files, fi_whitelist=wl)


def test_linking_unknown_whitelist_mode(linking_files, whitelist):
	"""An unknown whitelist_mode raises ValueError."""
	with pytest.raises(ValueError):
		chromatin.linking(*linking_files, fi_whitelist=whitelist, whitelist_mode='near')


def test_linking_unknown_combine(linking_files):
	"""An unknown combine method raises ValueError."""
	with pytest.raises(ValueError):
		chromatin.linking(*linking_files, combine='median')


def test_linking_no_shared_region(tmp_path):
	"""Binding and distance tables without a shared region raise RuntimeError."""
	fb = tmp_path / 'binding.tsv'
	fd = tmp_path / 'dist.tsv'
	pd.DataFrame([['A', 'chr1:1:2', 1.0]], columns=['TF', 'loc', 'score']).to_csv(fb, sep='\t', index=False)
	pd.DataFrame([['chr1:5:6', 'G1', 0]], columns=['region', 'target', 'dist']).to_csv(fd, sep='\t', index=False)
	with pytest.raises(RuntimeError):
		chromatin.linking(str(fb), str(fd), str(tmp_path / 'out.tsv'))


def test_linking_cli(linking_files, whitelist, run_dictys):
	"""``dictys chromatin linking`` with every option: intersect whitelist, sum, mode 0 (no distance penalty)."""
	run_dictys('chromatin', 'linking', *linking_files, '--fi_whitelist', whitelist, '--whitelist_mode', 'intersect',
		'--combine', 'sum', '--mode', 0)
	_assert_linking(linking_files[2], {('A', 'B'): 1.0, ('A', 'G1'): 2.0})


################################################################
# binlinking
################################################################

BIN_INDEX = ['A', 'B', 'C']
BIN_COLUMNS = ['B', 'G1', 'G2', 'G3', 'A', 'C']
BIN_VALUES = [
	[1.0, 5.0, 0.5, NINF, NINF, NINF],
	[NINF, 4.0, 2.0, 1.5, NINF, NINF],
	[0.2, 3.0, NINF, NINF, NINF, NINF],
]


@pytest.fixture
def binlinking_files(tmp_path):
	"""Link score matrix BIN_VALUES (TFs A, B, C by targets) and output path."""
	fi = tmp_path / 'linking.tsv'
	pd.DataFrame(BIN_VALUES, index=BIN_INDEX, columns=BIN_COLUMNS).to_csv(fi, sep='\t')
	return str(fi), str(tmp_path / 'binlinking.tsv')


def _assert_binlinking(path, links):
	"""``links`` is the set of selected (TF, target) pairs. Rows are TFs with a link; columns are targets with a link plus those TFs."""
	d = _read_table(path)
	tfs = {x[0] for x in links}
	assert set(d.index) == tfs
	assert set(d.columns) == tfs | {x[1] for x in links}
	got = {(x, y) for x in d.index for y in d.columns if d.loc[x, y]}
	assert d.values.dtype == bool
	assert got == links


ALL_FINITE = {('A', 'B'), ('A', 'G1'), ('A', 'G2'), ('B', 'G1'), ('B', 'G2'), ('B', 'G3'), ('C', 'B'), ('C', 'G1')}


@pytest.mark.parametrize('n,ka,links', [
	# strongest regulator of each target; A and C have no finite regulator
	(1, dict(), {('A', 'B'), ('A', 'G1'), ('B', 'G2'), ('B', 'G3')}),
	# two strongest; G3 has only one finite regulator
	(2, dict(axis=1), {('A', 'B'), ('C', 'B'), ('A', 'G1'), ('B', 'G1'), ('A', 'G2'), ('B', 'G2'), ('B', 'G3')}),
	# strongest target of each TF
	(1, dict(axis=0), {('A', 'G1'), ('B', 'G1'), ('C', 'G1')}),
	# two strongest links overall: 5.0 and 4.0
	(2, dict(axis=None), {('A', 'G1'), ('B', 'G1')}),
	(4, dict(axis=None), {('A', 'G1'), ('B', 'G1'), ('C', 'G1'), ('B', 'G2')}),
	(-1, dict(), ALL_FINITE),
	(-1, dict(axis=None), ALL_FINITE),
	(100, dict(axis=0), ALL_FINITE),
])
def test_binlinking_selection(binlinking_files, n, ka, links):
	"""The selected links match the strongest finite scores of BIN_VALUES along the chosen axis, picked by hand."""
	chromatin.binlinking(*binlinking_files, n, **ka)
	_assert_binlinking(binlinking_files[1], links)


def test_binlinking_self_regulation(tmp_path):
	"""A finite self-regulation score raises ValueError."""
	fi = tmp_path / 'linking.tsv'
	pd.DataFrame([[1.0, 2.0]], index=['A'], columns=['A', 'G1']).to_csv(fi, sep='\t')
	with pytest.raises(ValueError):
		chromatin.binlinking(str(fi), str(tmp_path / 'out.tsv'), 1)


def test_binlinking_unknown_selfreg(binlinking_files):
	"""An unknown selfreg option raises ValueError."""
	with pytest.raises(ValueError):
		chromatin.binlinking(*binlinking_files, 1, selfreg='ignore')


def test_binlinking_no_link(tmp_path):
	"""A matrix of only -inf scores raises RuntimeError."""
	fi = tmp_path / 'linking.tsv'
	pd.DataFrame([[NINF, NINF]], index=['A'], columns=['A', 'G1']).to_csv(fi, sep='\t')
	with pytest.raises(RuntimeError):
		chromatin.binlinking(str(fi), str(tmp_path / 'out.tsv'), 1)


def test_binlinking_cli(binlinking_files, run_dictys):
	"""``dictys chromatin binlinking 2 --axis 0 --selfreg error --inf never`` keeps each TF's two strongest targets."""
	run_dictys('chromatin', 'binlinking', *binlinking_files, 2, '--axis', 0, '--selfreg', 'error', '--inf', 'never')
	_assert_binlinking(binlinking_files[1], {('A', 'G1'), ('A', 'B'), ('B', 'G1'), ('B', 'G2'), ('C', 'G1'), ('C', 'B')})
