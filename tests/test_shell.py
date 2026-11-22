"""Tests for dictys.utils.shell.

Commands run the current Python interpreter through the shell, so they hold
under sh and cmd.exe alike. The Python code in a command avoids braces, which
cmdfile formats with the temporary folder's path.
"""

import os
import subprocess
import sys
import threading

import pytest

from dictys.utils.shell import cmdfile, run_timed, tmpfolder_context

PY = '"' + sys.executable + '"'


def pycmd(code):
	"""Shell command running the given Python code (no double quotes or braces inside)."""
	return PY + ' -c "' + code + '"'


# tmpfolder_context

def test_tmpfolder_context_creates_and_removes(tmp_path):
	"""A new folder exists inside the context and is removed after it."""
	with tmpfolder_context(dir=str(tmp_path)) as d:
		assert os.path.isdir(d)
		assert os.path.dirname(d) == os.path.realpath(str(tmp_path))
		with open(os.path.join(d, 'f.txt'), 'w') as f:
			f.write('x')
	assert not os.path.exists(d)


def test_tmpfolder_context_persist_keeps_folder(tmp_path):
	"""persist=True keeps the folder and its files."""
	with tmpfolder_context(persist=True, dir=str(tmp_path), prefix='keep_') as d:
		with open(os.path.join(d, 'f.txt'), 'w') as f:
			f.write('x')
	assert os.path.basename(d).startswith('keep_')
	assert os.path.isfile(os.path.join(d, 'f.txt'))


def test_tmpfolder_context_existing_folder(tmp_path):
	"""d names an existing folder to use; it is returned as a real path and removed unless persist."""
	d0 = tmp_path / 'given'
	d0.mkdir()
	with tmpfolder_context(d=str(d0), persist=True) as d:
		assert d == os.path.realpath(str(d0))
	assert d0.is_dir()
	with tmpfolder_context(d=str(d0)) as d:
		assert d == os.path.realpath(str(d0))
	assert not d0.exists()


def test_tmpfolder_context_removes_on_exception(tmp_path):
	"""The folder is removed when the body raises."""
	with pytest.raises(ValueError):
		with tmpfolder_context(dir=str(tmp_path)) as d:
			raise ValueError
	assert not os.path.exists(d)


# run_timed

def test_run_timed_returns_process_and_output():
	"""run returns [Popen, stdout, stderr] with the process finished."""
	ans = run_timed.run(pycmd('import sys; sys.stdout.write(chr(65)); sys.stderr.write(chr(66))'), 30,
		shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	assert len(ans) == 3
	assert isinstance(ans[0], subprocess.Popen)
	assert ans[0].returncode == 0
	assert ans[1] == b'A'
	assert ans[2] == b'B'


def test_run_timed_nonzero_exit_is_returned():
	"""A failing command is not an error for run_timed; its return code is kept."""
	ans = run_timed.run(pycmd('import sys; sys.exit(4)'), 30, shell=True)
	assert ans[0].returncode == 4
	assert ans[1] is None and ans[2] is None


def test_run_timed_object_records_command():
	"""The instance keeps the command and fills process and data after run_func."""
	t = run_timed(pycmd('print(1)'))
	assert t.process is None and t.data is None
	ans = t.run_func(30, shell=True, stdout=subprocess.PIPE)
	assert t.process is ans[0]
	assert ans[1].strip() == b'1'


@pytest.mark.filterwarnings('ignore::pytest.PytestUnhandledThreadExceptionWarning')
def test_run_timed_timeout_raises():
	"""A command running past the timeout raises TimeoutError."""
	with pytest.raises(TimeoutError):
		run_timed.run(pycmd('import time; time.sleep(5)'), 0.5, shell=True)


# cmdfile

def test_cmdfile_reads_output_file():
	"""outfile as str returns the bytes of that file in the temporary folder."""
	ans = cmdfile(pycmd("open('out.txt','w').write('hello')"), 'out.txt', cd=True)
	assert ans == b'hello'


def test_cmdfile_folder_placeholder_without_cd():
	"""Without cd the command reaches the temporary folder through the format placeholder."""
	ans = cmdfile(PY + ' -c "import sys; open(sys.argv[1],\'w\').write(\'abc\')" "{}' + os.sep + 'o.txt"', 'o.txt')
	assert ans == b'abc'


def test_cmdfile_input_files_are_written():
	"""infiles are written (str encoded, bytes as is) before the command runs."""
	code = "a=open('a.txt').read(); b=open('b.bin','rb').read(); open('out.txt','w').write(a+str(len(b)))"
	ans = cmdfile(pycmd(code), 'out.txt', infiles={'a.txt': 'text', 'b.bin': b'\x00\x01\x02'}, cd=True)
	assert ans == b'text3'


def test_cmdfile_list_of_outputs_and_missing_file():
	"""outfile as list returns one entry per name, None for a file not written."""
	code = "open('x.txt','w').write('1'); open('y.txt','w').write('22')"
	ans = cmdfile(pycmd(code), ['x.txt', 'missing.txt', 'y.txt'], cd=True)
	assert ans == [b'1', None, b'22']


def test_cmdfile_missing_single_output_is_none():
	"""A single outfile that the command did not write gives None."""
	assert cmdfile(pycmd('pass'), 'none.txt', cd=True) is None


def test_cmdfile_empty_output_list():
	"""An empty outfile list returns an empty list, as chromatin.macs2 expects on success."""
	assert cmdfile(pycmd('pass'), [], cd=True) == []


def test_cmdfile_stdout_and_stderr_when_no_outfile():
	"""outfile=None returns stdout with stderr merged, as str."""
	code = "import sys; sys.stdout.write('o'); sys.stdout.flush(); sys.stderr.write('e')"
	assert cmdfile(pycmd(code), None) == 'oe'


def test_cmdfile_nostderr_drops_stderr():
	"""nostderr=True returns stdout only."""
	code = "import sys; sys.stdout.write('o'); sys.stderr.write('e')"
	assert cmdfile(pycmd(code), None, nostderr=True) == 'o'


def test_cmdfile_quiet_false_prints_output(capsys):
	"""quiet=False prints the command's merged output after it finishes."""
	code = "open('f','w').write('z'); print('visible')"
	ans = cmdfile(pycmd(code), 'f', quiet=False, cd=True)
	assert ans == b'z'
	assert 'visible' in capsys.readouterr().out


def test_cmdfile_quiet_true_prints_nothing(capfd):
	"""quiet=True with an outfile sends the command's stdout and stderr nowhere (checked at the file descriptor level)."""
	code = "import sys; print('hidden'); sys.stderr.write('alsohidden'); open('f','w').write('z')"
	assert cmdfile(pycmd(code), 'f', quiet=True, cd=True) == b'z'
	out, err = capfd.readouterr()
	assert 'hidden' not in out and 'hidden' not in err


def test_cmdfile_quiet_integer_one_equals_true(capfd):
	"""quiet accepts 0 and 1 as integers: 1 suppresses the output, 0 prints it."""
	code = "print('qq'); open('f','w').write('q')"
	assert cmdfile(pycmd(code), 'f', quiet=1, cd=True) == b'q'
	assert 'qq' not in capfd.readouterr().out
	assert cmdfile(pycmd(code), 'f', quiet=0, cd=True) == b'q'
	assert 'qq' in capfd.readouterr().out


def test_cmdfile_quiet_two_not_supported():
	"""quiet=2 (real-time output) is not implemented."""
	with pytest.raises(AssertionError):
		cmdfile(pycmd('pass'), None, quiet=2)


def test_cmdfile_failure_raises():
	"""A non-zero exit raises RuntimeError by default."""
	with pytest.raises(RuntimeError):
		cmdfile(pycmd('import sys; sys.exit(2)'), 'x')


def test_cmdfile_sizelimit():
	"""An output file larger than sizelimit raises BufferError; one within it is returned; None means no limit."""
	code = "open('big','w').write('a'*100)"
	with pytest.raises(BufferError):
		cmdfile(pycmd(code), 'big', cd=True, sizelimit=99)
	assert cmdfile(pycmd(code), 'big', cd=True, sizelimit=100) == b'a' * 100
	assert cmdfile(pycmd(code), 'big', cd=True, sizelimit=None) == b'a' * 100


def test_cmdfile_tmpfolder_prefix_and_persist(tmp_path):
	"""tmpfolder with isprefix=True names a prefix for a new folder, kept with persist=True."""
	prefix = os.path.join(str(tmp_path), 'run_')
	ans = cmdfile(pycmd("open('o','w').write('k')"), 'o', tmpfolder=prefix, persist=True, cd=True)
	assert ans == b'k'
	made = [x for x in os.listdir(str(tmp_path)) if x.startswith('run_')]
	assert len(made) == 1
	with open(os.path.join(str(tmp_path), made[0], 'o')) as f:
		assert f.read() == 'k'


def test_cmdfile_tmpfolder_prefix_removed_without_persist(tmp_path):
	"""Without persist the prefixed temporary folder is removed."""
	prefix = os.path.join(str(tmp_path), 'run_')
	cmdfile(pycmd("open('o','w').write('k')"), 'o', tmpfolder=prefix, cd=True)
	assert [x for x in os.listdir(str(tmp_path)) if x.startswith('run_')] == []


def test_cmdfile_tmpfolder_full_path(tmp_path):
	"""isprefix=False runs in the given folder itself."""
	d = tmp_path / 'work'
	d.mkdir()
	ans = cmdfile(pycmd("open('o','w').write('w')"), 'o', tmpfolder=str(d), isprefix=False, persist=True, cd=True)
	assert ans == b'w'
	assert (d / 'o').read_text() == 'w'


@pytest.mark.filterwarnings('ignore::pytest.PytestUnhandledThreadExceptionWarning')
def test_cmdfile_timeout_raises():
	"""A command running past the timeout raises TimeoutError."""
	with pytest.raises(TimeoutError):
		cmdfile(pycmd('import time; time.sleep(5)'), None, timeout=0.5)


def test_cmdfile_negative_timeout_runs_in_thread(tmp_path):
	"""A negative timeout returns a started thread that runs the command."""
	d = tmp_path / 'bg'
	d.mkdir()
	t = cmdfile(pycmd("open('o','w').write('bg')"), 'o', tmpfolder=str(d), isprefix=False, persist=True,
		cd=True, timeout=-30)
	assert isinstance(t, threading.Thread)
	t.join(60)
	assert not t.is_alive()
	assert (d / 'o').read_text() == 'bg'


@pytest.mark.parametrize('ka', [
	dict(cmd=1, outfile=None),
	dict(cmd='x', outfile=1),
	dict(cmd='x', outfile=None, infiles=[]),
	dict(cmd='x', outfile=None, persist=1),
	dict(cmd='x', outfile=None, noerror=1),
	dict(cmd='x', outfile=None, sizelimit=-1),
	dict(cmd='x', outfile=None, sizelimit=1.5),
	dict(cmd='x', outfile=None, cd=1),
])
def test_cmdfile_argument_checks(ka):
	"""Arguments of the wrong type are rejected before anything runs."""
	cmd = ka.pop('cmd')
	outfile = ka.pop('outfile')
	with pytest.raises(AssertionError):
		cmdfile(cmd, outfile, **ka)
