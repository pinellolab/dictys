"""
Tests that the macs2 step still succeeds where the serial script did: a chunk
whose reads are all filtered out, and a thread count of 0. They need bash and
a POSIX shell; samtools and macs2 are replaced by stand-ins.
"""

import pytest

from dictys import chromatin
from test_chromatin import stubs, macs2_stubs, macs2_files, _macs2_args, _write_text, SAM_HEADER, SAM_CELL1, SAM_READ  # noqa: F401


def test_macs2_all_filtered_chunk(macs2_files, macs2_stubs):
	"""A chunk whose every read is removed by the filters (one cell per chunk at nth=2) does not fail the step."""
	_write_text(macs2_files['fi_bam'] + '/cell2.bam', 'r3' + SAM_READ.format('chrM', 5) + '\nr4' + SAM_READ.format('chrUn_KI270302v1', 5) + '\n')
	chromatin.macs2(*_macs2_args(macs2_files), 'hs', nth=2, nmax=0)
	with open(macs2_files['fo_bam']) as f:
		assert f.read().splitlines() == [SAM_HEADER[1], SAM_CELL1[3]]


def test_macs2_zero_threads(macs2_files, macs2_stubs, capfd):
	"""nth=0 is accepted as before, without a shell arithmetic error (division by 0 in the chunk size); samtools gets -@ 0."""
	chromatin.macs2(*_macs2_args(macs2_files), 'hs', nth=0, nmax=0)
	out, err = capfd.readouterr()
	assert 'division' not in out + err
	with open(macs2_files['fo_bam']) as f:
		assert f.read().splitlines()[0] == SAM_HEADER[1]
	calls = macs2_stubs.calls('samtools')
	assert all(x[x.index('-@') + 1] == '0' for x in calls)
