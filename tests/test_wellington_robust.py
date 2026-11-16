"""
Tests that the wellington step fails or starts cleanly: an env shebang on the
pyDNase script, and a BAM that pysam cannot read. They need bash and a POSIX
shell; pyDNase is replaced by a stand-in package, so it need not be installed.
"""

import os
import subprocess
import sys

import pytest

from dictys import chromatin
from test_chromatin import stubs, _write_rows, _write_text, _wellington_args, _read_bed, wellington_files, FOOTPRINTS, FOOTPRINTS_OUT  # noqa: F401

SCRIPTS = os.path.join(os.path.dirname(chromatin.__file__), 'scripts')


@pytest.fixture
def env_shebang_stubs(stubs, require, tmp_path, monkeypatch):
	"""pyDNase's script behind an ``#!/usr/bin/env python3`` shebang, and a python3 stand-in that reports FOOTPRINTS."""
	require('bedtools')
	_write_text(os.path.join(stubs.bindir, 'wellington_footprints.py'), '#!/usr/bin/env python3\nimport pyDNase\n')
	os.chmod(os.path.join(stubs.bindir, 'wellington_footprints.py'), 0o755)
	stubs.add('python3', 'mkdir -p "$4/p value cutoffs"\ncp "$STUB_FP" "$4/p value cutoffs/fp${10}.bed"\n')
	monkeypatch.setenv('STUB_FP', _write_rows(tmp_path / 'fp.bed', FOOTPRINTS))
	return stubs


def test_wellington_env_shebang(wellington_files, env_shebang_stubs):
	"""A pyDNase script with an env shebang still runs dictys's driver through python3, not through env itself (exit 126)."""
	chromatin.wellington(*_wellington_args(wellington_files), cut=10, nth=3, nmax=0)
	assert _read_bed(wellington_files['fo_bed']) == [FOOTPRINTS_OUT[x] for x in ['fp1', 'fp3', 'fp4', 'fp5']]
	(call,) = env_shebang_stubs.calls('python3')
	assert call[0] == os.path.join(SCRIPTS, 'wellington_footprints.py')
	assert call[call.index('-p') + 1] == '2'


FAKE_PYDNASE = '''
class BAMHandler:
	def __init__(self, *args, **kwargs):
		raise OSError('could not open alignment file')

class GenomicInterval:
	def __init__(self, chromosome, startbp, endbp, strand='+'):
		self.chromosome, self.startbp, self.endbp, self.strand = chromosome, startbp, endbp, strand

class GenomicIntervalSet:
	def __init__(self, path):
		self.intervals = {'chr1': [GenomicInterval('chr1', 0, 100)]}
'''


def test_wellington_unreadable_bam(tmp_path):
	"""A BAM that cannot be opened ends the driver at once with exit 1 and a message, instead of respawning workers forever."""
	pkg = tmp_path / 'fake' / 'pyDNase'
	pkg.mkdir(parents=True)
	_write_text(pkg / '__init__.py', FAKE_PYDNASE)
	_write_text(pkg / 'footprinting.py', 'wellington = wellington1D = None\n')
	bam = _write_text(tmp_path / 'bad.bam', 'not a bam file\n')
	bed = _write_rows(tmp_path / 'peaks.bed', [['chr1', 0, 100]])
	out = tmp_path / 'out'
	out.mkdir()
	env = dict(os.environ, PYTHONPATH=str(tmp_path / 'fake') + os.pathsep + os.path.join(os.path.dirname(chromatin.__file__), '..'))
	try:
		proc = subprocess.run([sys.executable, os.path.join(SCRIPTS, 'wellington_footprints.py'), '-p', '2', bed, bam, str(out)],
			env=env, capture_output=True, text=True, timeout=60)
	except subprocess.TimeoutExpired as e:
		pytest.fail('driver hangs on an unreadable BAM')
	assert proc.returncode == 1
	assert 'could not open alignment file' in proc.stderr
