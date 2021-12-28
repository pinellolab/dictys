# Lingfei Wang, 2018. All rights reserved.
"""
Shell script wrapper
"""

assert __name__ != "__main__"

from contextlib import contextmanager


@contextmanager
def tmpfolder_context(persist=False, d=None,**ka):
	"""Create temporary folder
	persist:	Whether to keep folder after exiting context
	d:			Path of an existing folder to use instead.
	ka:			Keyword args passed to tempfile.mkdtemp if d is None"""
	from tempfile import mkdtemp
	from os.path import realpath
	if d is None:
		d = mkdtemp(**ka)
	d = realpath(d)
	try:
		yield d
	finally:
		if not persist:
			from shutil import rmtree
			rmtree(d)


class run_timed(object):
	"""Class to run process with time out"""
	def __init__(self, cmd):
		"""Class to run process with time out
		cmd: command to execute"""
		self.cmd = cmd
		self.process = None
		self.data = None

	@classmethod
	def run(cls, cmd, timeout, **ka):
		"""Function wrapper to run process with time out.
		cmd: 		command to execute
		timeout:	Timeout
		ka:			Keyword arguments passed to subprocess.Popen
		Return:		[subprocess.Popen(),stdout,stderr (both as b'')]
		Raises TimeoutError on time-out."""
		t = cls(cmd)
		return t.run_func(timeout, **ka)

	def run_func(self, timeout, **ka):
		"""Run process with time-out
		timeout:	Timeout
		ka:			Keyword arguments passed to subprocess.Popen
		Return:		[subprocess.Popen(),stdout,stderr (both as b'')]
		Raises TimeoutError on time-out."""
		import threading, logging

		def target(cmd, timeout, **ka):
			import subprocess
			self.process = subprocess.Popen(cmd, **ka)
			try:
				self.data = self.process.communicate(timeout=timeout)
			except subprocess.TimeoutExpired:
				raise TimeoutError('Command failed due to time out: ' + self.cmd)

		thread = threading.Thread(target=target, args=(self.cmd, timeout), kwargs=ka)
		thread.start()
		thread.join(timeout)
		if thread.is_alive():
			if self.process is not None:
				self.process.kill()
			thread.join()
			raise TimeoutError('Command failed due to time out: ' + self.cmd)
		elif self.data is None:
			raise TimeoutError('Command failed due to time out: ' + self.cmd)
		return [self.process] + list(self.data)


def cmdfile(cmd,
			outfile,
			infiles=dict(),
			tmpfolder=None,
			isprefix=True,
			persist=False,
			noerror=False,
			timeout=None,
			sizelimit=1000000000,
			quiet=True,
			cd=False,
			nostderr=False):
	"""Run command and return the content of file.
	cmd:		string of command to run. The actual command is cmd.format(tmpfolder).
	outfile:	Relative file path to return contents as str or file paths as list of str.
				The absolute file path is os.path.join(tmpfolder,file).
				If None, return stdout and stderr instead.
	infiles:	Dict of input file relative path:contents to write before command execution.
	tmpfolder:	Specify temporary folder path or prefix.
	isprefix:	Whether tmpfolder indicates prefix of tmp folder, or full path.
	persist:	Whether to keep tmp folder after execution.
	noerror:	Whether to ignore nonzero return value from command. If not, an RuntimeError will be thrown.
	timeout:	Time out for command in seconds.
			Default:	Wait till return
			<0:		Immediate return. No wait. Command executed in separate thread with abs(timeout) time-out.
	sizelimit:	The size limit of output. If file is bigger, BufferError will be raised.
				Set to None for unlimited.
	quiet:		How to handle stdout and stderr
			0, False:	Both to stdout after cmd finishes
			1, True:	Suppressed
			2:			Both straight to stdout in realtime. Not Implemented.
	cd:			Change working directory to temp dir before execution.
	nostderr:	Whether to ignore stderr in output. Only used when outfile is None.
	Return:		For timeout>=0 or None: Contents of the designated file as bytes or list of bytes for outfile list.
				Or None if file does not exist.
				For timeout<0: thread object"""
	assert type(cmd) is str
	assert outfile is None or type(outfile) is str or type(outfile) is list
	assert type(infiles) is dict
	assert type(persist) is bool and type(noerror) is bool
	assert sizelimit is None or (type(sizelimit) is int and sizelimit >= 0)
	if type(quiet) is bool:
		quiet = int(quiet)
	assert type(quiet) is int and quiet in {0, 1}
	assert type(cd) is bool
	if timeout is not None and timeout < 0:
		import numpy as np
		from threading import Thread
		t1 = Thread(target=cmdfile,
					args=(cmd, outfile),
					kwargs=dict(infiles=infiles,
								tmpfolder=tmpfolder,
								isprefix=isprefix,
								persist=persist,
								noerror=noerror,
								timeout=-timeout,
								sizelimit=sizelimit,
								quiet=quiet,
								cd=cd))
		t1.start()
		return t1
	import logging
	from os.path import join as pjoin
	import subprocess
	import os

	# Create tmp folder
	ka = dict()
	if tmpfolder is not None:
		assert type(tmpfolder) is str and type(isprefix) is bool
		ka['prefix' if isprefix else 'd'] = tmpfolder

	with tmpfolder_context(persist=persist, **ka) as td:
		logging.debug('Using temporary folder: ' + td)
		# Prepare input files
		for fn in infiles:
			assert type(fn) is str and (type(infiles[fn]) is str or
										type(infiles[fn]) is bytes)
			fi = pjoin(td, fn)
			logging.debug('Writing temporary input file: ' + fi)
			f = open(fi, 'wb')
			f.write(infiles[fn].encode()if type(infiles[fn]) is str else infiles[fn])
			f.close()

		# Run command
		cmda = cmd.format(td)
		ka = dict()
		if quiet == 1 and outfile is not None:
			ka['stdout'] = subprocess.DEVNULL
			ka['stderr'] = subprocess.DEVNULL
		else:
			ka['stdout'] = subprocess.PIPE
			ka['stderr'] = subprocess.DEVNULL if nostderr else subprocess.STDOUT
		if cd:
			ka['cwd'] = td
		logging.debug('Calling program ' + cmda)
		try:
			ans = run_timed.run(cmda, timeout, shell=True, **ka)
		except TimeoutError as e:
			if not noerror:
				raise

		if quiet == 0 and len(ans) > 1 and ans[1] is not None:
			print(ans[1].decode())
		if not noerror and ans[0] is None or ans[0].returncode != 0:
			raise RuntimeError('Command failed, possibly due to program error: ' + cmda)

		if outfile is None:
			if len(ans) > 1 and ans[1] is not None:
				return ans[1].decode()
			else:
				return

		# Read output file
		def readfile(fout):
			fo = pjoin(td, fout)
			logging.debug('Reading temporary output file: ' + fo)
			try:
				t1 = os.stat(fo).st_size
				if sizelimit is not None and t1 > sizelimit:
					raise BufferError('File too large ({}): '.format(t1) + fo)
			except OSError:
				return
			f = open(fo, 'rb')
			ans = f.read()
			f.close()
			return ans

		ans = readfile(outfile) if type(outfile) is str else [
			readfile(x) for x in outfile]
	return ans
