"""
Tests for the repository helpers repo_dataverse.py and repo_zenodo.py.

No test touches the network: ``requests`` and ``pyDataverse`` are replaced in
``sys.modules`` by stubs that record every request and serve canned responses.
The stubs are used whether or not the real packages are installed. The
pyDataverse stub follows pyDataverse 0.3.1, the version utils/conda/meta.yaml
installs; 0.3.2 and later build a different download URL (a bug in the base).
"""

import io
import os
import sys
import types

import pytest


class HTTPError(Exception):
	"""Raised by FakeResponse.raise_for_status, and exported as requests.exceptions.HTTPError."""


class FakeResponse:
	"""Minimal stand-in for a requests/pyDataverse response."""
	def __init__(self, json=None, status=200, content=b''):
		self._json = json
		self.status_code = status
		self.content = content
		self.raw = io.BytesIO(content)

	def json(self):
		return self._json

	def raise_for_status(self):
		if self.status_code >= 400:
			raise HTTPError(self.status_code)

	def __enter__(self):
		return self

	def __exit__(self, *a):
		return False


class FakeHTTPBasicAuth:
	def __init__(self, username, password):
		self.username = username
		self.password = password


@pytest.fixture
def fake_requests(monkeypatch):
	"""
	Stub ``requests`` module. ``routes[(method, url)]`` is the response served;
	``calls`` lists (method, url, keyword arguments), with a file passed as
	``data`` replaced by its bytes.
	"""
	mod = types.ModuleType('requests')
	exc = types.ModuleType('requests.exceptions')
	exc.HTTPError = HTTPError
	auth = types.ModuleType('requests.auth')
	auth.HTTPBasicAuth = FakeHTTPBasicAuth
	mod.exceptions = exc
	mod.auth = auth
	mod.calls = []
	mod.routes = {}

	def method(name):
		def f(url, **ka):
			ka = dict(ka)
			if hasattr(ka.get('data'), 'read'):
				ka['data'] = ka['data'].read()
			mod.calls.append((name, url, ka))
			if (name, url) not in mod.routes:
				raise KeyError('unexpected request {} {}'.format(name, url))
			return mod.routes[(name, url)]
		return f
	mod.get = method('GET')
	mod.put = method('PUT')
	mod.delete = method('DELETE')
	monkeypatch.setitem(sys.modules, 'requests', mod)
	monkeypatch.setitem(sys.modules, 'requests.exceptions', exc)
	monkeypatch.setitem(sys.modules, 'requests.auth', auth)
	return mod


def read_text(path):
	with open(path, 'rb') as f:
		return f.read()


def out_lines(capsys):
	return [x for x in capsys.readouterr().out.splitlines() if x]


######################################################################
# repo_dataverse.py
######################################################################

DV_BASE = 'http://dataverse.test'
DV_DOI = '10.1234/DVN/ABC'


@pytest.fixture
def fake_dataverse(monkeypatch, fake_requests):
	"""
	Stub ``pyDataverse.api`` with NativeApi and DataAccessApi, with the signatures of
	pyDataverse 0.3.1. ``state.files`` is the dataset's file list as the Dataverse API
	returns it, ``state.content`` maps file IDs to bytes, ``state.calls`` records every
	API call. A dataset or upload addressed by a non-persistent ID gets 404, and
	get_datafile records the path it would request: in 0.3.1 ``is_pid=True`` (the
	default) gives ``access/datafile/<id>``, so only that path serves content.
	"""
	state = types.SimpleNamespace(calls=[], files=[], content={}, info_status=200, dataset_status=200,
		upload_status=200, file_status=200)

	class NativeApi:
		def __init__(self, base_url, api_token=None, api_version='v1'):
			state.calls.append(('NativeApi', base_url, api_token))

		def get_info_version(self, auth=False):
			state.calls.append(('get_info_version',))
			return FakeResponse(json={'data': {'version': '5'}}, status=state.info_status)

		def get_dataset(self, identifier, version=':latest', auth=True, is_pid=True):
			state.calls.append(('get_dataset', identifier))
			status = state.dataset_status if is_pid and version == ':latest' else 404
			return FakeResponse(json={'data': {'latestVersion': {'files': state.files}}}, status=status)

		def upload_datafile(self, identifier, filename, json_str=None, is_pid=True):
			state.calls.append(('upload_datafile', identifier, filename, read_text(filename)))
			return FakeResponse(status=state.upload_status if is_pid else 404)

	class DataAccessApi:
		def __init__(self, base_url, api_token=None):
			state.calls.append(('DataAccessApi', base_url, api_token))

		def get_datafile(self, identifier, data_format=None, no_var_header=None, image_thumb=None, is_pid=True,
				auth=False):
			if is_pid and not (data_format or no_var_header or image_thumb):
				path = 'access/datafile/{}'.format(identifier)
				content, status = state.content.get(identifier, b''), state.file_status
			else:
				path = 'access/datafile/:persistentId/?persistentId={}'.format(identifier)
				content, status = b'', 404
			state.calls.append(('get_datafile', path))
			return FakeResponse(content=content, status=status)

	pkg = types.ModuleType('pyDataverse')
	api = types.ModuleType('pyDataverse.api')
	api.NativeApi = NativeApi
	api.DataAccessApi = DataAccessApi
	pkg.api = api
	monkeypatch.setitem(sys.modules, 'pyDataverse', pkg)
	monkeypatch.setitem(sys.modules, 'pyDataverse.api', api)
	return state


def dv_file(label, fid, filename=None):
	return {'label': label, 'dataFile': {'filename': label if filename is None else filename, 'id': fid}}


@pytest.fixture
def dataverse(run_helper, fake_dataverse, fake_requests):
	"""Globals of repo_dataverse.py after one `ls` run against the stubs; state and call logs are then reset."""
	g = run_helper('repo_dataverse.py', 'ls', '--base', DV_BASE, '--doi', DV_DOI)
	fake_dataverse.calls.clear()
	fake_requests.calls.clear()
	return g


class TestDvd:
	def test_init(self, dataverse, fake_dataverse):
		"""The constructor stores url, doi, 'doi:'-prefixed pid and token, and checks the server version."""
		d = dataverse['dvd'](DV_BASE, DV_DOI, token='tok')
		assert (d.url, d.doi, d.pid, d.token) == (DV_BASE, DV_DOI, 'doi:' + DV_DOI, 'tok')
		assert fake_dataverse.calls == [('NativeApi', DV_BASE, 'tok'), ('DataAccessApi', DV_BASE, 'tok'),
			('get_info_version',)]

	def test_init_token_default(self, dataverse, fake_dataverse):
		"""Without a token the APIs are built with token None."""
		d = dataverse['dvd'](DV_BASE, DV_DOI)
		assert d.token is None
		assert fake_dataverse.calls[0] == ('NativeApi', DV_BASE, None)

	def test_test_failure(self, dataverse, fake_dataverse):
		"""test() raises when the server version request fails, so the constructor does too."""
		d = dataverse['dvd'](DV_BASE, DV_DOI)
		fake_dataverse.info_status = 503
		with pytest.raises(HTTPError):
			d.test()
		with pytest.raises(HTTPError):
			dataverse['dvd'](DV_BASE, DV_DOI)

	def test_ls(self, dataverse, fake_dataverse):
		"""ls maps each file label to its data file ID for the dataset's pid."""
		fake_dataverse.files = [dv_file('a.txt', 11), dv_file('b.tar.xz.0', 12)]
		d = dataverse['dvd'](DV_BASE, DV_DOI)
		assert d.ls() == {'a.txt': 11, 'b.tar.xz.0': 12}
		assert ('get_dataset', 'doi:' + DV_DOI) in fake_dataverse.calls

	def test_ls_label_mismatch(self, dataverse, fake_dataverse):
		"""A label that differs from the data file name is rejected."""
		fake_dataverse.files = [dv_file('a.txt', 11, filename='other.txt')]
		with pytest.raises(AssertionError):
			dataverse['dvd'](DV_BASE, DV_DOI).ls()

	def test_ls_duplicate_labels(self, dataverse, fake_dataverse):
		"""Two files with the same label are rejected."""
		fake_dataverse.files = [dv_file('a.txt', 11), dv_file('a.txt', 12)]
		with pytest.raises(AssertionError):
			dataverse['dvd'](DV_BASE, DV_DOI).ls()

	def test_ls_http_error(self, dataverse, fake_dataverse):
		"""A failed dataset request raises."""
		fake_dataverse.dataset_status = 404
		with pytest.raises(HTTPError):
			dataverse['dvd'](DV_BASE, DV_DOI).ls()

	def test_delete_byid(self, dataverse, fake_requests):
		"""delete_byid sends DELETE to the SWORD edit-media URL with the token as basic-auth user."""
		url = DV_BASE + '/dvn/api/data-deposit/v1.1/swordv2/edit-media/file/42'
		fake_requests.routes[('DELETE', url)] = FakeResponse()
		dataverse['dvd'](DV_BASE, DV_DOI, token='tok').delete_byid(42)
		assert len(fake_requests.calls) == 1
		method, u, ka = fake_requests.calls[0]
		assert (method, u) == ('DELETE', url)
		assert (ka['auth'].username, ka['auth'].password) == ('tok', '')

	def test_delete_byid_failure(self, dataverse, fake_requests):
		"""A failed DELETE raises."""
		url = DV_BASE + '/dvn/api/data-deposit/v1.1/swordv2/edit-media/file/42'
		fake_requests.routes[('DELETE', url)] = FakeResponse(status=403)
		with pytest.raises(HTTPError):
			dataverse['dvd'](DV_BASE, DV_DOI, token='tok').delete_byid(42)

	def test_upload_new_file(self, dataverse, fake_dataverse, fake_requests, tmp_path):
		"""A file absent from the dataset is uploaded without any delete."""
		f = tmp_path / 'new.txt'
		f.write_bytes(b'abc')
		fake_dataverse.files = [dv_file('a.txt', 11)]
		dataverse['dvd'](DV_BASE, DV_DOI).upload(str(f))
		assert ('upload_datafile', 'doi:' + DV_DOI, str(f), b'abc') in fake_dataverse.calls
		assert fake_requests.calls == []

	def test_upload_overwrite(self, dataverse, fake_dataverse, fake_requests, tmp_path):
		"""exist='overwrite' deletes the existing file by its ID, then uploads."""
		f = tmp_path / 'a.txt'
		f.write_bytes(b'new')
		fake_dataverse.files = [dv_file('a.txt', 11)]
		url = DV_BASE + '/dvn/api/data-deposit/v1.1/swordv2/edit-media/file/11'
		fake_requests.routes[('DELETE', url)] = FakeResponse()
		dataverse['dvd'](DV_BASE, DV_DOI, token='tok').upload(str(f), exist='overwrite')
		assert [x[:2] for x in fake_requests.calls] == [('DELETE', url)]
		assert fake_dataverse.calls[-1] == ('upload_datafile', 'doi:' + DV_DOI, str(f), b'new')

	def test_upload_skip(self, dataverse, fake_dataverse, fake_requests, tmp_path):
		"""exist='skip' leaves an existing file alone and uploads nothing."""
		f = tmp_path / 'a.txt'
		f.write_bytes(b'new')
		fake_dataverse.files = [dv_file('a.txt', 11)]
		dataverse['dvd'](DV_BASE, DV_DOI).upload(str(f), exist='skip')
		assert not [x for x in fake_dataverse.calls if x[0] == 'upload_datafile']
		assert fake_requests.calls == []

	def test_upload_unknown_exist(self, dataverse, fake_dataverse, tmp_path):
		"""An unknown exist value raises ValueError when the file already exists."""
		f = tmp_path / 'a.txt'
		f.write_bytes(b'new')
		fake_dataverse.files = [dv_file('a.txt', 11)]
		with pytest.raises(ValueError, match='Unknown exist=bad'):
			dataverse['dvd'](DV_BASE, DV_DOI).upload(str(f), exist='bad')

	def test_upload_dst_not_implemented(self, dataverse, tmp_path):
		"""A destination name is not supported."""
		with pytest.raises(NotImplementedError):
			dataverse['dvd'](DV_BASE, DV_DOI).upload(str(tmp_path / 'a.txt'), dst='b.txt')

	def test_upload_failure(self, dataverse, fake_dataverse, tmp_path):
		"""A failed upload raises."""
		f = tmp_path / 'a.txt'
		f.write_bytes(b'new')
		fake_dataverse.upload_status = 500
		with pytest.raises(HTTPError):
			dataverse['dvd'](DV_BASE, DV_DOI).upload(str(f))

	def test_download_byid(self, dataverse, fake_dataverse, tmp_path):
		"""download_byid requests the file by its numeric ID and writes its bytes to the destination."""
		fake_dataverse.content[7] = b'\x00payload'
		dst = tmp_path / 'x.bin'
		dataverse['dvd'](DV_BASE, DV_DOI).download_byid(7, str(dst))
		assert dst.read_bytes() == b'\x00payload'
		assert fake_dataverse.calls[-1] == ('get_datafile', 'access/datafile/7')

	def test_download_byid_failure(self, dataverse, fake_dataverse, tmp_path):
		"""A failed file request raises and writes nothing."""
		fake_dataverse.file_status = 404
		dst = tmp_path / 'x.bin'
		with pytest.raises(HTTPError):
			dataverse['dvd'](DV_BASE, DV_DOI).download_byid(7, str(dst))
		assert not dst.exists()

	def test_download(self, dataverse, fake_dataverse, tmp_path):
		"""download finds the file ID by name and writes to dst."""
		fake_dataverse.files = [dv_file('a.txt', 11), dv_file('b.txt', 12)]
		fake_dataverse.content.update({11: b'A', 12: b'B'})
		dst = tmp_path / 'got.txt'
		dataverse['dvd'](DV_BASE, DV_DOI).download('b.txt', dst=str(dst))
		assert dst.read_bytes() == b'B'

	def test_download_default_dst(self, dataverse, fake_dataverse, tmp_path, monkeypatch):
		"""Without dst the file is written under its own name in the working directory."""
		fake_dataverse.files = [dv_file('a.txt', 11)]
		fake_dataverse.content[11] = b'A'
		monkeypatch.chdir(tmp_path)
		dataverse['dvd'](DV_BASE, DV_DOI).download('a.txt')
		assert (tmp_path / 'a.txt').read_bytes() == b'A'

	def test_download_missing(self, dataverse, fake_dataverse, tmp_path):
		"""A name not in the dataset raises FileNotFoundError."""
		fake_dataverse.files = [dv_file('a.txt', 11)]
		with pytest.raises(FileNotFoundError):
			dataverse['dvd'](DV_BASE, DV_DOI).download('zzz', dst=str(tmp_path / 'z'))


class TestDataverseCommandLine:
	def test_ls(self, run_helper, fake_dataverse, capsys):
		"""`ls` prints the file names sorted, one per line, using --base, --doi and --token."""
		fake_dataverse.files = [dv_file('b.txt', 2), dv_file('a.txt', 1), dv_file('c.txt', 3)]
		run_helper('repo_dataverse.py', 'ls', '--base', DV_BASE, '--doi', DV_DOI, '--token', 'tok')
		assert out_lines(capsys) == ['a.txt', 'b.txt', 'c.txt']
		assert fake_dataverse.calls[0] == ('NativeApi', DV_BASE, 'tok')
		assert ('get_dataset', 'doi:' + DV_DOI) in fake_dataverse.calls

	def test_defaults(self, run_helper, fake_dataverse):
		"""Without options the Harvard Dataverse base URL and the Dictys dataset DOI are used."""
		run_helper('repo_dataverse.py', 'ls')
		assert fake_dataverse.calls[0] == ('NativeApi', 'https://dataverse.harvard.edu', None)
		assert ('get_dataset', 'doi:10.7910/DVN/IJWFWR') in fake_dataverse.calls

	def test_down(self, run_helper, fake_dataverse, tmp_path):
		"""`down` with --dst writes the named file there."""
		fake_dataverse.files = [dv_file('a.txt', 11)]
		fake_dataverse.content[11] = b'hello'
		dst = tmp_path / 'saved.txt'
		run_helper('repo_dataverse.py', 'down', 'a.txt', '--dst', str(dst), '--base', DV_BASE, '--doi', DV_DOI,
			'--token', 'tok')
		assert dst.read_bytes() == b'hello'

	def test_down_default_dst(self, run_helper, fake_dataverse, tmp_path):
		"""`down` without --dst writes into the working directory."""
		fake_dataverse.files = [dv_file('a.txt', 11)]
		fake_dataverse.content[11] = b'hello'
		run_helper('repo_dataverse.py', 'down', 'a.txt', '--base', DV_BASE, cwd=tmp_path)
		assert (tmp_path / 'a.txt').read_bytes() == b'hello'

	@pytest.mark.parametrize('existing,uploads,deletes', [('overwrite', 1, 1), ('skip', 0, 0)])
	def test_up(self, run_helper, fake_dataverse, fake_requests, tmp_path, existing, uploads, deletes):
		"""`up` overwrites or skips an existing file according to --existing."""
		f = tmp_path / 'a.txt'
		f.write_bytes(b'v2')
		fake_dataverse.files = [dv_file('a.txt', 11)]
		fake_requests.routes[('DELETE', DV_BASE + '/dvn/api/data-deposit/v1.1/swordv2/edit-media/file/11')] = FakeResponse()
		run_helper('repo_dataverse.py', 'up', str(f), '--existing', existing, '--base', DV_BASE, '--doi', DV_DOI,
			'--token', 'tok')
		assert len([x for x in fake_dataverse.calls if x[0] == 'upload_datafile']) == uploads
		assert len(fake_requests.calls) == deletes

	def test_up_default_overwrites(self, run_helper, fake_dataverse, fake_requests, tmp_path):
		"""`up` overwrites by default."""
		f = tmp_path / 'a.txt'
		f.write_bytes(b'v2')
		fake_dataverse.files = [dv_file('a.txt', 11)]
		url = DV_BASE + '/dvn/api/data-deposit/v1.1/swordv2/edit-media/file/11'
		fake_requests.routes[('DELETE', url)] = FakeResponse()
		run_helper('repo_dataverse.py', 'up', str(f), '--base', DV_BASE)
		assert [x[:2] for x in fake_requests.calls] == [('DELETE', url)]
		assert fake_dataverse.calls[-1][0] == 'upload_datafile'

	def test_subcommand_required(self, run_helper, fake_dataverse):
		"""A subcommand is required."""
		with pytest.raises(SystemExit):
			run_helper('repo_dataverse.py')

	def test_module_functions(self, dataverse, fake_dataverse, tmp_path, capsys):
		"""down, up and ls take an argparse namespace and act on the dataset it names."""
		import argparse
		fake_dataverse.files = [dv_file('a.txt', 11)]
		fake_dataverse.content[11] = b'x'
		ns = argparse.Namespace(base=DV_BASE, doi=DV_DOI, token=None)
		dataverse['ls'](ns)
		assert out_lines(capsys) == ['a.txt']
		dst = tmp_path / 'd.txt'
		dataverse['down'](argparse.Namespace(src='a.txt', dst=str(dst), **vars(ns)))
		assert dst.read_bytes() == b'x'
		f = tmp_path / 'n.txt'
		f.write_bytes(b'n')
		dataverse['up'](argparse.Namespace(src=str(f), existing='skip', **vars(ns)))
		assert fake_dataverse.calls[-1] == ('upload_datafile', 'doi:' + DV_DOI, str(f), b'n')


######################################################################
# repo_zenodo.py
######################################################################

ZB = 'http://zenodo.test'
ZDOI = '10.5281/zenodo.777'


def zenodo_record(files):
	return FakeResponse(json={'files': [{'key': x, 'links': {'self': y}} for x, y in files.items()]})


def zenodo_deposition(files):
	return FakeResponse(json=[{'filename': x, 'links': {'download': y}} for x, y in files.items()])


@pytest.fixture
def zenodo(run_helper, fake_requests):
	"""Globals of repo_zenodo.py after one `ls` run against the stub; routes and calls are then reset."""
	fake_requests.routes[('GET', ZB + '/api/records/777')] = zenodo_record({})
	g = run_helper('repo_zenodo.py', 'ls', '--base', ZB, '--doi', ZDOI)
	fake_requests.routes.clear()
	fake_requests.calls.clear()
	return g


class TestZenodod:
	def test_init(self, zenodo, fake_requests):
		"""The record ID is the part of the DOI after its last dot; nothing is requested."""
		d = zenodo['zenodod'](ZDOI, url=ZB, token='tok')
		assert (d.url, d.doi, d.pid, d.token) == (ZB, ZDOI, '777', 'tok')
		assert fake_requests.calls == []

	def test_init_defaults(self, zenodo):
		"""The default base URL is zenodo.org and the default token is None."""
		d = zenodo['zenodod'](ZDOI)
		assert (d.url, d.token) == ('https://zenodo.org', None)

	def test_test(self, zenodo, fake_requests):
		"""test() requests the records endpoint and raises on failure."""
		fake_requests.routes[('GET', ZB + '/api/records')] = FakeResponse()
		d = zenodo['zenodod'](ZDOI, url=ZB)
		d.test()
		assert fake_requests.calls == [('GET', ZB + '/api/records', {})]
		fake_requests.routes[('GET', ZB + '/api/records')] = FakeResponse(status=500)
		with pytest.raises(HTTPError):
			d.test()

	@pytest.mark.parametrize('token,params', [(None, {}), ('tok', {'access_token': 'tok'})])
	def test_ls_record(self, zenodo, fake_requests, token, params):
		"""mode='record' maps each file key to its 'self' link, passing the token as access_token."""
		url = ZB + '/api/records/777'
		fake_requests.routes[('GET', url)] = zenodo_record({'a.txt': 'L1', 'b.txt': 'L2'})
		fs = zenodo['zenodod'](ZDOI, url=ZB, token=token).ls(mode='record')
		assert fs == {'a.txt': 'L1', 'b.txt': 'L2'}
		assert fake_requests.calls == [('GET', url, {'params': params})]

	@pytest.mark.parametrize('token,params', [(None, {}), ('tok', {'access_token': 'tok'})])
	def test_ls_deposition(self, zenodo, fake_requests, token, params):
		"""mode='deposition' maps each file name to its download link."""
		url = ZB + '/api/deposit/depositions/777/files'
		fake_requests.routes[('GET', url)] = zenodo_deposition({'a.txt': 'D1'})
		fs = zenodo['zenodod'](ZDOI, url=ZB, token=token).ls(mode='deposition')
		assert fs == {'a.txt': 'D1'}
		assert fake_requests.calls == [('GET', url, {'params': params})]

	def test_ls_all_uses_record_first(self, zenodo, fake_requests):
		"""mode='all' returns the record listing without asking for the deposition when it succeeds."""
		fake_requests.routes[('GET', ZB + '/api/records/777')] = zenodo_record({'a.txt': 'L1'})
		assert zenodo['zenodod'](ZDOI, url=ZB).ls() == {'a.txt': 'L1'}
		assert [x[1] for x in fake_requests.calls] == [ZB + '/api/records/777']

	def test_ls_all_falls_back_to_deposition(self, zenodo, fake_requests):
		"""mode='all' falls back to the deposition listing when the record request fails."""
		fake_requests.routes[('GET', ZB + '/api/records/777')] = FakeResponse(status=404)
		fake_requests.routes[('GET', ZB + '/api/deposit/depositions/777/files')] = zenodo_deposition({'d.txt': 'D'})
		assert zenodo['zenodod'](ZDOI, url=ZB).ls(mode='all') == {'d.txt': 'D'}
		assert len(fake_requests.calls) == 2

	def test_ls_all_fails(self, zenodo, fake_requests):
		"""mode='all' raises RuntimeError when both listings fail."""
		fake_requests.routes[('GET', ZB + '/api/records/777')] = FakeResponse(status=404)
		fake_requests.routes[('GET', ZB + '/api/deposit/depositions/777/files')] = FakeResponse(status=403)
		with pytest.raises(RuntimeError, match='All modes failed'):
			zenodo['zenodod'](ZDOI, url=ZB).ls()

	def test_ls_unknown_mode(self, zenodo):
		"""An unknown mode raises ValueError."""
		with pytest.raises(ValueError, match='Unknown mode bad'):
			zenodo['zenodod'](ZDOI, url=ZB).ls(mode='bad')

	def test_delete_byid(self, zenodo, fake_requests):
		"""Deleting files is not implemented for Zenodo and sends no request."""
		with pytest.raises(NotImplementedError):
			zenodo['zenodod'](ZDOI, url=ZB, token='tok').delete_byid(1)
		assert fake_requests.calls == []

	@pytest.mark.parametrize('token,params', [(None, {}), ('tok', {'access_token': 'tok'})])
	def test_upload_datafile(self, zenodo, fake_requests, tmp_path, token, params):
		"""upload_datafile reads the deposition's bucket and PUTs the file bytes to bucket/dst."""
		f = tmp_path / 'a.txt'
		f.write_bytes(b'data')
		fake_requests.routes[('GET', ZB + '/api/deposit/depositions/777')] = FakeResponse(json={'links': {'bucket': 'http://b.test/bk'}})
		put = FakeResponse(json={'ok': 1})
		fake_requests.routes[('PUT', 'http://b.test/bk/remote.txt')] = put
		ans = zenodo['zenodod'](ZDOI, url=ZB, token=token).upload_datafile(str(f), 'remote.txt')
		assert ans is put
		assert fake_requests.calls == [('GET', ZB + '/api/deposit/depositions/777', {'params': params}),
			('PUT', 'http://b.test/bk/remote.txt', {'data': b'data', 'params': params})]

	@pytest.mark.parametrize('step', ['GET', 'PUT'])
	def test_upload_datafile_failure(self, zenodo, fake_requests, tmp_path, step):
		"""A failed deposition or upload request raises."""
		f = tmp_path / 'a.txt'
		f.write_bytes(b'data')
		fake_requests.routes[('GET', ZB + '/api/deposit/depositions/777')] = FakeResponse(
			json={'links': {'bucket': 'http://b.test/bk'}}, status=401 if step == 'GET' else 200)
		fake_requests.routes[('PUT', 'http://b.test/bk/a.txt')] = FakeResponse(status=500)
		with pytest.raises(HTTPError):
			zenodo['zenodod'](ZDOI, url=ZB).upload_datafile(str(f), 'a.txt')

	def _upload_routes(self, fake_requests, existing):
		fake_requests.routes[('GET', ZB + '/api/records/777')] = zenodo_record({x: 'L' for x in existing})
		fake_requests.routes[('GET', ZB + '/api/deposit/depositions/777')] = FakeResponse(json={'links': {'bucket': 'http://b.test/bk'}})
		fake_requests.routes[('PUT', 'http://b.test/bk/a.txt')] = FakeResponse()

	@pytest.mark.parametrize('existing,exist,nput', [([], 'skip', 1), (['a.txt'], 'overwrite', 1), (['a.txt'], 'skip', 0)])
	def test_upload(self, zenodo, fake_requests, tmp_path, existing, exist, nput):
		"""upload sends the file unless it exists and exist='skip'."""
		f = tmp_path / 'a.txt'
		f.write_bytes(b'v')
		self._upload_routes(fake_requests, existing)
		zenodo['zenodod'](ZDOI, url=ZB).upload(str(f), exist=exist)
		assert len([x for x in fake_requests.calls if x[0] == 'PUT']) == nput

	def test_upload_unknown_exist(self, zenodo, fake_requests, tmp_path):
		"""An unknown exist value raises ValueError when the file exists."""
		f = tmp_path / 'a.txt'
		f.write_bytes(b'v')
		self._upload_routes(fake_requests, ['a.txt'])
		with pytest.raises(ValueError, match='Unknown exist=bad'):
			zenodo['zenodod'](ZDOI, url=ZB).upload(str(f), exist='bad')

	def test_upload_dst_not_implemented(self, zenodo, tmp_path):
		"""A destination name is not supported."""
		with pytest.raises(NotImplementedError):
			zenodo['zenodod'](ZDOI, url=ZB).upload(str(tmp_path / 'a.txt'), dst='b.txt')

	def test_download_byurl(self, zenodo, fake_requests, tmp_path):
		"""download_byurl streams the response body into the destination file."""
		fake_requests.routes[('GET', 'http://f.test/x')] = FakeResponse(content=b'\x01\x02body')
		dst = tmp_path / 'x.bin'
		zenodo['zenodod'].download_byurl('http://f.test/x', str(dst))
		assert dst.read_bytes() == b'\x01\x02body'
		assert fake_requests.calls == [('GET', 'http://f.test/x', {'stream': True})]

	def test_download(self, zenodo, fake_requests, tmp_path):
		"""download looks up the file's link by name and saves it to dst."""
		fake_requests.routes[('GET', ZB + '/api/records/777')] = zenodo_record({'a.txt': 'http://f.test/a', 'b.txt': 'http://f.test/b'})
		fake_requests.routes[('GET', 'http://f.test/b')] = FakeResponse(content=b'B')
		dst = tmp_path / 'got.txt'
		zenodo['zenodod'](ZDOI, url=ZB).download('b.txt', dst=str(dst))
		assert dst.read_bytes() == b'B'

	def test_download_default_dst(self, zenodo, fake_requests, tmp_path, monkeypatch):
		"""Without dst the file is saved under its own name in the working directory."""
		fake_requests.routes[('GET', ZB + '/api/records/777')] = zenodo_record({'a.txt': 'http://f.test/a'})
		fake_requests.routes[('GET', 'http://f.test/a')] = FakeResponse(content=b'A')
		monkeypatch.chdir(tmp_path)
		zenodo['zenodod'](ZDOI, url=ZB).download('a.txt')
		assert (tmp_path / 'a.txt').read_bytes() == b'A'

	def test_download_missing(self, zenodo, fake_requests, tmp_path):
		"""A name not in the listing raises FileNotFoundError."""
		fake_requests.routes[('GET', ZB + '/api/records/777')] = zenodo_record({'a.txt': 'http://f.test/a'})
		with pytest.raises(FileNotFoundError):
			zenodo['zenodod'](ZDOI, url=ZB).download('zzz', dst=str(tmp_path / 'z'))


class TestZenodoCommandLine:
	def test_ls(self, run_helper, fake_requests, capsys):
		"""`ls` prints sorted file names, using --base, --doi and --token."""
		url = ZB + '/api/records/777'
		fake_requests.routes[('GET', url)] = zenodo_record({'b': 'x', 'a': 'y', 'c': 'z'})
		run_helper('repo_zenodo.py', 'ls', '--base', ZB, '--doi', ZDOI, '--token', 'tok')
		assert out_lines(capsys) == ['a', 'b', 'c']
		assert fake_requests.calls == [('GET', url, {'params': {'access_token': 'tok'}})]

	def test_defaults(self, run_helper, fake_requests):
		"""Without options zenodo.org and the Dictys record 6787658 are used."""
		url = 'https://zenodo.org/api/records/6787658'
		fake_requests.routes[('GET', url)] = zenodo_record({})
		run_helper('repo_zenodo.py', 'ls')
		assert fake_requests.calls == [('GET', url, {'params': {}})]

	def test_down(self, run_helper, fake_requests, tmp_path):
		"""`down` with --dst saves the named file there."""
		fake_requests.routes[('GET', ZB + '/api/records/777')] = zenodo_record({'a.txt': 'http://f.test/a'})
		fake_requests.routes[('GET', 'http://f.test/a')] = FakeResponse(content=b'hello')
		dst = tmp_path / 'saved.txt'
		run_helper('repo_zenodo.py', 'down', 'a.txt', '--dst', str(dst), '--base', ZB, '--doi', ZDOI, '--token', 'tok')
		assert dst.read_bytes() == b'hello'

	def test_down_default_dst(self, run_helper, fake_requests, tmp_path):
		"""`down` without --dst saves into the working directory."""
		fake_requests.routes[('GET', ZB + '/api/records/777')] = zenodo_record({'a.txt': 'http://f.test/a'})
		fake_requests.routes[('GET', 'http://f.test/a')] = FakeResponse(content=b'hello')
		run_helper('repo_zenodo.py', 'down', 'a.txt', '--base', ZB, '--doi', ZDOI, cwd=tmp_path)
		assert (tmp_path / 'a.txt').read_bytes() == b'hello'

	@pytest.mark.parametrize('existing,nput', [('overwrite', 1), ('skip', 0)])
	def test_up(self, run_helper, fake_requests, tmp_path, existing, nput):
		"""`up` uploads an existing file again only with --existing overwrite."""
		f = tmp_path / 'a.txt'
		f.write_bytes(b'v2')
		fake_requests.routes[('GET', ZB + '/api/records/777')] = zenodo_record({'a.txt': 'L'})
		fake_requests.routes[('GET', ZB + '/api/deposit/depositions/777')] = FakeResponse(json={'links': {'bucket': 'http://b.test/bk'}})
		fake_requests.routes[('PUT', 'http://b.test/bk/a.txt')] = FakeResponse()
		run_helper('repo_zenodo.py', 'up', str(f), '--existing', existing, '--base', ZB, '--doi', ZDOI, '--token', 'tok')
		puts = [x for x in fake_requests.calls if x[0] == 'PUT']
		assert len(puts) == nput
		if nput:
			assert puts[0][2] == {'data': b'v2', 'params': {'access_token': 'tok'}}

	def test_subcommand_required(self, run_helper, fake_requests):
		"""A subcommand is required."""
		with pytest.raises(SystemExit):
			run_helper('repo_zenodo.py')

	def test_module_functions(self, zenodo, fake_requests, tmp_path, capsys):
		"""down, up and ls take an argparse namespace and act on the record it names."""
		import argparse
		fake_requests.routes[('GET', ZB + '/api/records/777')] = zenodo_record({'a.txt': 'http://f.test/a'})
		fake_requests.routes[('GET', 'http://f.test/a')] = FakeResponse(content=b'x')
		ns = argparse.Namespace(base=ZB, doi=ZDOI, token=None)
		zenodo['ls'](ns)
		assert out_lines(capsys) == ['a.txt']
		dst = tmp_path / 'd.txt'
		zenodo['down'](argparse.Namespace(src='a.txt', dst=str(dst), **vars(ns)))
		assert dst.read_bytes() == b'x'
		f = tmp_path / 'a.txt'
		f.write_bytes(b'n')
		zenodo['up'](argparse.Namespace(src=str(f), existing='skip', **vars(ns)))
		assert not [x for x in fake_requests.calls if x[0] == 'PUT']
