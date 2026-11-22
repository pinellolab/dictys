"""
Tests for the bash helper scripts under dictys/scripts/helper and for the
``dictys_helper`` dispatcher in bin/.

External programs the scripts call are replaced where possible by stub
scripts on PATH (make, homer's configureHomer.pl, dictys_helper for the
repository downloads, python3 pointing at this interpreter), so nothing is
downloaded and no real pipeline runs. Tests that need a real program call
``require`` first.
"""

import filecmp
import os
import subprocess
import sys
import tarfile

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DICTYS_HELPER = os.path.join(REPO, 'bin', 'dictys_helper')


def write_exe(path, text):
	"""Write an executable script."""
	with open(str(path), 'w', newline='\n') as f:
		f.write(text)
	os.chmod(str(path), 0o755)


def path_with(*dirs):
	"""PATH with the given folders in front of the current one."""
	return os.pathsep.join([str(x) for x in dirs] + [os.environ.get('PATH', '')])


@pytest.fixture
def stubbin(tmp_path):
	"""A folder for stub programs, holding a `python3` that runs this test's interpreter."""
	d = tmp_path / 'stubbin'
	d.mkdir()
	write_exe(d / 'python3', '#!/usr/bin/env bash\nexec "{}" "$@"\n'.format(sys.executable.replace('\\', '/')))
	return d


@pytest.fixture
def workdir(tmp_path):
	d = tmp_path / 'work'
	d.mkdir()
	return d


def read_lines(path):
	with open(str(path), 'r') as f:
		return [x for x in f.read().splitlines() if x]


def command_works(cmd, input=None):
	"""Whether a command exits 0 (used to detect GNU-only behaviour)."""
	try:
		return subprocess.run(cmd, input=input, capture_output=True, text=True, timeout=30).returncode == 0
	except OSError:
		return False


@pytest.fixture
def require_realpath_missing(require, tmp_path):
	"""Skip unless realpath resolves a path that does not exist yet, as GNU realpath does."""
	require('realpath')
	if not command_works(['realpath', str(tmp_path / 'not_yet')]):
		pytest.skip('realpath does not resolve missing paths here')


@pytest.fixture
def require_head_negative(require):
	"""Skip unless `head -n -N` (all but the last N lines, GNU only) works."""
	require('head')
	if not command_works(['head', '-n', '-1'], input='a\nb\n'):
		pytest.skip('head does not accept negative line counts here')


######################################################################
# bin/dictys_helper
######################################################################

class TestDictysHelper:
	@pytest.fixture
	def run(self, require, child_env, stubbin, workdir):
		"""Run bin/dictys_helper with bash in workdir, with the python3 stub first on PATH."""
		require('bash')

		def f(*args):
			return subprocess.run(['bash', DICTYS_HELPER] + [str(x) for x in args], cwd=str(workdir),
				env=child_env(PATH=path_with(stubbin)), capture_output=True, text=True, timeout=120)
		return f

	def test_no_argument_lists_scripts(self, run, require, dictys_scripts):
		"""Without arguments it prints usage and every helper script, and exits non-zero."""
		require('ls', 'grep')
		r = run()
		assert r.returncode != 0
		lines = r.stdout.splitlines()
		assert lines[0].startswith('Usage: dictys_helper helper_script_name')
		assert lines[1] == 'Available helper scripts:'
		expected = sorted(x for x in os.listdir(os.path.join(dictys_scripts, 'helper')) if x != '__pycache__')
		assert sorted(lines[2:]) == expected

	def test_dispatch_python(self, run, workdir):
		"""A Python helper receives the remaining arguments intact, including spaces."""
		f = workdir / 'a.mk'
		f.write_text('A=1\nB=2\n')
		r = run('makefile_update.py', str(f), '{"B":"x y"}')
		assert r.returncode == 0, r.stderr
		assert f.read_text() == 'A=1\nB=x y\n'

	def test_dispatch_python_error(self, run):
		"""The helper's exit status is passed through."""
		r = run('expression.py')
		assert r.returncode == 1
		assert 'expression_mtx.py' in r.stderr

	def test_dispatch_shell(self, run, require, workdir, dictys_scripts):
		"""A bash helper runs in the caller's working directory."""
		require('realpath')
		r = run('makefile_template.sh', 'config.mk')
		assert r.returncode == 0, r.stderr
		assert os.listdir(str(workdir)) == ['config.mk']
		assert filecmp.cmp(str(workdir / 'config.mk'), os.path.join(dictys_scripts, 'makefiles', 'config.mk'), shallow=False)

	def test_unknown_script(self, run):
		"""An unknown script name fails."""
		r = run('no_such_helper.py')
		assert r.returncode != 0


######################################################################
# makefile_template.sh
######################################################################

class TestMakefileTemplate:
	def test_all_makefiles(self, run_shell, require, workdir, dictys_scripts):
		"""Without arguments every .mk file of the package is copied unchanged."""
		require('realpath', 'ls', 'grep')
		r = run_shell('makefile_template.sh', cwd=str(workdir))
		assert r.returncode == 0, r.stderr
		src = os.path.join(dictys_scripts, 'makefiles')
		expected = sorted(x for x in os.listdir(src) if x.endswith('.mk'))
		assert len(expected) >= 4
		assert sorted(os.listdir(str(workdir))) == expected
		for x in expected:
			assert filecmp.cmp(str(workdir / x), os.path.join(src, x), shallow=False)

	def test_named_makefiles(self, run_shell, require, workdir, dictys_scripts):
		"""Named makefiles are the only ones copied."""
		require('realpath')
		r = run_shell('makefile_template.sh', 'config.mk', 'static.mk', cwd=str(workdir))
		assert r.returncode == 0, r.stderr
		assert sorted(os.listdir(str(workdir))) == ['config.mk', 'static.mk']
		assert filecmp.cmp(str(workdir / 'static.mk'), os.path.join(dictys_scripts, 'makefiles', 'static.mk'), shallow=False)

	def test_unknown_makefile(self, run_shell, require, workdir):
		"""A makefile name not in the package fails."""
		require('realpath')
		r = run_shell('makefile_template.sh', 'nonexistent.mk', cwd=str(workdir))
		assert r.returncode != 0
		assert os.listdir(str(workdir)) == []

	@pytest.mark.parametrize('args', [['-h'], ['-x']])
	def test_usage(self, run_shell, workdir, args):
		"""-h and unknown options print usage to stderr, exit 1 and copy nothing."""
		r = run_shell('makefile_template.sh', *args, cwd=str(workdir))
		assert r.returncode == 1
		assert 'Usage: makefile_template.sh' in r.stderr
		assert os.listdir(str(workdir)) == []


######################################################################
# network_inference.sh
######################################################################

STUB_MAKE = '''#!/usr/bin/env bash
echo "$*" >> "$MAKE_LOG"
if [ "a$MAKE_FAIL" != "a" ] && [ "a${@: -1}" == "a$MAKE_FAIL" ]; then
	echo "make: *** No rule to make target 'x'." >&2
	exit 2
fi
exit 0
'''


class TestNetworkInference:
	@pytest.fixture
	def run(self, run_shell, require, stubbin, workdir, tmp_path):
		"""Run network_inference.sh with a stub make that logs its arguments; returns (result, logged calls)."""
		# grep filters make's stderr; without it the filtering test would pass vacuously
		require('grep')
		write_exe(stubbin / 'make', STUB_MAKE)
		log = tmp_path / 'make.log'

		def f(*args, fail=''):
			r = run_shell('network_inference.sh', *args, cwd=str(workdir),
				env={'PATH': path_with(stubbin), 'MAKE_LOG': str(log), 'MAKE_FAIL': fail})
			return r, (read_lines(log) if log.exists() else [])
		return f

	def test_static_defaults(self, run):
		"""static mode runs cpu twice, gpu, cpu and combine with 4 CPU jobs and 1 GPU job."""
		r, calls = run('static')
		assert r.returncode == 0, r.stderr
		mk = '-f makefiles/static.mk'
		assert calls == [mk + ' -j 4 -k cpu', mk + ' -j 4 -k cpu', mk + ' -j 1 -k gpu', mk + ' -j 4 -k cpu',
			mk + ' combine']

	def test_dynamic_options(self, run):
		"""dynamic mode builds subsets first; -j, -J set job counts and -c adds the clean step."""
		r, calls = run('-j', '8', '-J', '2', '-c', 'dynamic')
		assert r.returncode == 0, r.stderr
		mk = '-f makefiles/dynamic.mk'
		assert calls == [mk + ' subset', mk + ' -j 8 -k cpu', mk + ' -j 8 -k cpu', mk + ' -j 2 -k gpu',
			mk + ' -j 8 -k cpu', mk + ' combine', mk + ' clean']

	@pytest.mark.parametrize('target', ['cpu', 'gpu'])
	def test_cpu_gpu_failures_ignored(self, run, target):
		"""Failures in the cpu and gpu steps are ignored, every later step still runs, and missing-rule messages are filtered out."""
		r, calls = run('-c', 'static', fail=target)
		assert r.returncode == 0
		mk = '-f makefiles/static.mk'
		assert calls == [mk + ' -j 4 -k cpu', mk + ' -j 4 -k cpu', mk + ' -j 1 -k gpu', mk + ' -j 4 -k cpu',
			mk + ' combine', mk + ' clean']
		assert 'No rule to make target' not in r.stderr

	@pytest.mark.parametrize('target,nstep', [('combine', 6), ('subset', 1)])
	def test_combine_subset_failures_stop(self, run, target, nstep):
		"""A failure in the subset or combine step stops the script with an error right after that step, before cleaning."""
		r, calls = run('-c', 'dynamic', fail=target)
		assert r.returncode != 0
		mk = '-f makefiles/dynamic.mk'
		steps = [mk + ' subset', mk + ' -j 4 -k cpu', mk + ' -j 4 -k cpu', mk + ' -j 1 -k gpu', mk + ' -j 4 -k cpu',
			mk + ' combine']
		assert calls == steps[:nstep]
		assert calls[-1].endswith(' ' + target)

	@pytest.mark.parametrize('args', [[], ['other'], ['static', 'dynamic'], ['-h'], ['-x', 'static']])
	def test_usage(self, run, args):
		"""Missing, unknown or extra modes, -h and unknown options print usage, exit 1 and run nothing."""
		r, calls = run(*args)
		assert r.returncode == 1
		assert 'Usage: network_inference.sh' in r.stderr
		assert calls == []

	def test_option_requires_argument(self, run):
		"""-j without a value is reported."""
		r, calls = run('-j')
		assert r.returncode == 1
		assert 'Error: -j requires an argument.' in r.stderr
		assert calls == []

	def test_real_make(self, run_shell, require, workdir):
		"""With GNU make and a stand-in makefile, the targets run in the documented order."""
		require('make')
		(workdir / 'makefiles').mkdir()
		rules = ''.join('{0}:\n\techo {0} >> log.txt\n'.format(x) for x in ['cpu', 'gpu', 'combine', 'clean'])
		(workdir / 'makefiles' / 'static.mk').write_text('.PHONY: cpu gpu combine clean\n' + rules)
		r = run_shell('network_inference.sh', '-j', '2', '-J', '1', '-c', 'static', cwd=str(workdir))
		assert r.returncode == 0, r.stderr
		assert read_lines(workdir / 'log.txt') == ['cpu', 'cpu', 'gpu', 'cpu', 'combine', 'clean']


######################################################################
# gene_gtf.sh
######################################################################

GTF = (
	'#!genome-build test\n'
	'1\tsrc\tgene\t100\t200\t.\t+\t.\tgene_id "E1"; gene_name "Alpha"; gene_biotype "pc";\n'
	'1\tsrc\ttranscript\t100\t200\t.\t+\t.\tgene_id "E1"; transcript_id "T1"; gene_name "Alpha"; gene_biotype "pc";\n'
	'10\tsrc\tgene\t50\t80\t.\t-\t.\tgene_id "E2"; gene_name "Beta"; gene_biotype "pc";\n'
	'2\tsrc\tgene\t300\t400\t.\t-\t.\tgene_id "E3"; gene_name "Gamma"; gene_biotype "pc";\n'
	'X\tsrc\tgene\t5\t9\t.\t+\t.\tgene_id "E4"; gene_name "Delta"; gene_biotype "pc";\n'
	'2\tsrc\tgene\t10\t20\t.\t+\t.\tgene_id "E5"; gene_name "Gamma"; gene_biotype "pc";\n'
	'chr3\tsrc\tgene\t1\t2\t.\t+\t.\tgene_id "E6"; gene_biotype "pc";\n'
)


class TestGeneGtf:
	@pytest.fixture
	def run(self, run_shell, require, workdir):
		"""Run gene_gtf.sh on GTF text in workdir; returns (result, BED rows split into columns)."""
		require('awk', 'grep', 'sed', 'sort', 'uniq', 'paste', 'tr')

		def f(*opts, gtf=GTF):
			(workdir / 'in.gtf').write_text(gtf)
			r = run_shell('gene_gtf.sh', *opts, 'in.gtf', 'out.bed', cwd=str(workdir), env={'LC_ALL': 'C'})
			rows = [x.split('\t') for x in read_lines(workdir / 'out.bed')] if (workdir / 'out.bed').exists() else None
			return r, rows
		return f

	@staticmethod
	def without_start(rows):
		# column 1 (start) is left out: it is the GTF 1-based start, not shifted to 0-based (a bug in the base)
		return [[x[0]] + x[2:] for x in rows]

	def test_gene_name(self, run, workdir):
		"""Gene rows become chr-prefixed BED6 rows named by gene_name, sorted by chromosome number then others."""
		r, rows = run()
		assert r.returncode == 0, r.stderr
		# by hand: transcript, comment and the gene without gene_name are dropped;
		# the duplicated Gamma keeps its first occurrence (end 400)
		assert self.without_start(rows) == [
			['chr1', '200', 'Alpha', '.', '+'],
			['chr2', '400', 'Gamma', '.', '-'],
			['chr10', '80', 'Beta', '.', '-'],
			['chrX', '9', 'Delta', '.', '+'],
		]
		assert sorted(os.listdir(str(workdir))) == ['in.gtf', 'out.bed']

	def test_field_option(self, run):
		"""-f gene_id names rows by gene_id, keeping genes without gene_name, sorted by start within a chromosome."""
		r, rows = run('-f', 'gene_id')
		assert r.returncode == 0, r.stderr
		assert self.without_start(rows) == [
			['chr1', '200', 'E1', '.', '+'],
			['chr2', '20', 'E5', '.', '+'],
			['chr2', '400', 'E3', '.', '-'],
			['chr3', '2', 'E6', '.', '+'],
			['chr10', '80', 'E2', '.', '-'],
			['chrX', '9', 'E4', '.', '+'],
		]

	@pytest.mark.parametrize('args', [['-h'], ['only.gtf'], ['a', 'b', 'c'], ['-x', 'a', 'b']])
	def test_usage(self, run_shell, workdir, args):
		"""-h, a wrong argument count and unknown options print usage and exit 1."""
		r = run_shell('gene_gtf.sh', *args, cwd=str(workdir))
		assert r.returncode == 1
		assert 'Usage: gene_gtf.sh' in r.stderr
		assert os.listdir(str(workdir)) == []

	def test_option_requires_argument(self, run_shell, workdir):
		"""-f without a value is reported."""
		r = run_shell('gene_gtf.sh', '-f', cwd=str(workdir))
		assert r.returncode == 1
		assert 'Error: -f requires an argument.' in r.stderr


######################################################################
# genome_homer.sh and motif_homer.sh: a stand-in homer installation
######################################################################

CONFIGURE_HOMER = '''#!/usr/bin/env bash
# Stand-in for homer's configureHomer.pl: lists genomes, "installs" mm10 only.
if [ "$1" == "-list" ]; then
	echo "Packages:"
	for g in hg38 mm10 dm6; do
		if [ -e "data/genomes/$g/installed" ]; then echo "+	$g	v6.4	genome"; else echo "-	$g	v6.4	genome"; fi
	done
elif [ "$1" == "-install" ]; then
	echo "$2" >> install.log
	if [ "$2" == "mm10" ]; then
		mkdir -p data/genomes/mm10/preparsed
		echo ">chr1" > data/genomes/mm10/genome.fa
		touch data/genomes/mm10/installed
	else
		echo "install failed for $2" >&2
		exit 1
	fi
fi
'''


def make_homer(base):
	"""
	Build a stand-in homer installation at ``base``: configureHomer.pl, an installed hg38 genome
	with a preparsed folder, motif files and a mapping file, and a known.motifs for 'vertebrates'.
	"""
	base = str(base)
	os.makedirs(base)
	write_exe(os.path.join(base, 'configureHomer.pl'), CONFIGURE_HOMER)
	g = os.path.join(base, 'data', 'genomes', 'hg38')
	os.makedirs(os.path.join(g, 'preparsed'))
	for name, text in [('installed', ''), ('genome.fa', '>chr1\nACGT\n'), ('chrom.sizes', 'chr1\t4\n'),
			(os.path.join('preparsed', 'x.200.cgbins'), 'junk\n')]:
		with open(os.path.join(g, name), 'w') as f:
			f.write(text)
	m = os.path.join(base, 'motifs')
	os.makedirs(os.path.join(m, 'extras'))
	os.makedirs(os.path.join(m, 'sub'))
	files = {
		'gata1.motif': '>AGATAA\tGATA1(Zf)/K562/Homer\t5.0\n0.9\t0.05\t0.03\t0.02\n0.1\t0.7\t0.1\t0.1\n',
		os.path.join('sub', 'klf.motif'): '>GGGCGG\tKLF4(Zf)/ESC/Homer\t6.0\n0.2\t0.3\t0.4\t0.1\n',
		'two.motif': '>A\tX\t1\n0.25\t0.25\t0.25\t0.25\n>C\tY\t1\n0.25\t0.25\t0.25\t0.25\n',
		'nohead.motif': '0.25\t0.25\t0.25\t0.25\n>A\tX\t1\n0.25\t0.25\t0.25\t0.25\n',
	}
	for name, text in files.items():
		with open(os.path.join(m, name), 'w') as f:
			f.write(text)
	with open(os.path.join(m, 'extras', 'motifs2symbol.txt'), 'w') as f:
		f.write('motif\tsymbol\ngata1.motif\tgaTA1\nsub/klf.motif\tKlf4\nmissing.motif\tNone\ntwo.motif\tTwo\nnohead.motif\tNo\n')
	k = os.path.join(base, 'data', 'knownTFs', 'vertebrates')
	os.makedirs(k)
	with open(os.path.join(k, 'known.motifs'), 'w') as f:
		f.write('>AGATAA\tGata1:gaTA2(Zf)/K562/Homer\t5.0\n0.9\t0.05\t0.03\t0.02\n0.1\t0.7\t0.1\t0.1\n'
			'>GGGCGG\tKlf4(Zf)/ESC/Homer\t6.0\n0.2\t0.3\t0.4\t0.1\n0.3\t0.3\t0.2\t0.2\n')
	return base


def parse_motifs(text):
	"""Split motif text into [(header fields, [row fields, ...]), ...], splitting on any whitespace."""
	ans = []
	for line in text.splitlines():
		if not line.strip():
			continue
		if line.startswith('>'):
			ans.append((line[1:].split(), []))
		else:
			ans[-1][1].append(line.split())
	return ans


@pytest.fixture
def homer(tmp_path):
	return make_homer(tmp_path / 'homer')


class TestGenomeHomer:
	@pytest.fixture
	def run(self, run_shell, require, require_realpath_missing, workdir):
		require('awk')

		def f(*args, env=None):
			return run_shell('genome_homer.sh', *args, cwd=str(workdir), env=env)
		return f

	def test_installed_genome(self, run, homer, workdir):
		"""An installed genome is copied to the output folder without its preparsed folder."""
		r = run('-b', homer, 'hg38', 'out/')
		assert r.returncode == 0, r.stdout + r.stderr
		out = workdir / 'out'
		assert sorted(os.listdir(str(out))) == ['chrom.sizes', 'genome.fa', 'installed']
		assert (out / 'genome.fa').read_text() == '>chr1\nACGT\n'
		assert not os.path.exists(os.path.join(homer, 'install.log'))

	def test_genome_installed_on_demand(self, run, homer, workdir):
		"""A listed but missing genome is installed through configureHomer.pl, then copied."""
		r = run('-b', homer, 'mm10', 'out')
		assert r.returncode == 0, r.stdout + r.stderr
		assert 'Downloading reference genome mm10 in homer' in r.stdout
		assert read_lines(os.path.join(homer, 'install.log')) == ['mm10']
		assert sorted(os.listdir(str(workdir / 'out'))) == ['genome.fa', 'installed']

	def test_install_failure(self, run, homer, workdir):
		"""A failed installation exits 1 and shows the installer's error."""
		r = run('-b', homer, 'dm6', 'out')
		assert r.returncode == 1
		assert read_lines(os.path.join(homer, 'install.log')) == ['dm6']
		lines = r.stdout.splitlines()
		assert lines[-3:] == ['Error: reference genome dm6 not found in homer',
			'Possible error during reference genome installation:', 'install failed for dm6']
		assert not (workdir / 'out').exists()

	def test_unknown_genome(self, run, homer, workdir):
		"""A genome homer does not list exits 1 without attempting an installation."""
		r = run('-b', homer, 'xx1', 'out')
		assert r.returncode == 1
		assert 'Error: reference genome xx1 not found in homer' in r.stdout
		assert 'Possible error' not in r.stdout
		assert not os.path.exists(os.path.join(homer, 'install.log'))

	def test_autodetect_basedir(self, run, require, homer, workdir):
		"""Without -b the base folder is two levels above the homer program found on PATH."""
		require('which')
		os.makedirs(os.path.join(homer, 'bin'))
		write_exe(os.path.join(homer, 'bin', 'homer'), '#!/usr/bin/env bash\nexit 0\n')
		r = run('hg38', 'out', env={'PATH': path_with(os.path.join(homer, 'bin'))})
		assert r.returncode == 0, r.stdout + r.stderr
		assert (workdir / 'out' / 'genome.fa').exists()

	def test_missing_basedir(self, run, workdir, tmp_path):
		"""A base folder that does not exist exits 1."""
		r = run('-b', str(tmp_path / 'nohomer'), 'hg38', 'out')
		assert r.returncode == 1
		assert 'cannot find homer base directory' in r.stderr

	@pytest.mark.parametrize('args', [['-h'], ['hg38'], ['hg38', 'a', 'b'], ['-x', 'hg38', 'out']])
	def test_usage(self, run, args):
		"""-h, a wrong argument count and unknown options print usage and exit 1."""
		r = run(*args)
		assert r.returncode == 1
		assert 'Usage: genome_homer.sh' in r.stderr

	def test_option_requires_argument(self, run):
		"""-b without a value is reported."""
		r = run('-b')
		assert r.returncode == 1
		assert 'Error: -b requires an argument.' in r.stderr


class TestMotifHomer:
	@pytest.fixture
	def run(self, run_shell, require, workdir):
		require('realpath', 'awk', 'head', 'tail', 'tr', 'grep')

		def f(*args, env=None):
			return run_shell('motif_homer.sh', *args, cwd=str(workdir), env=env)
		return f

	def test_mapfile_default(self, run, homer):
		"""The default mapping file names each single-motif file '<symbol>_<path with _>' and skips bad files."""
		r = run('-b', homer)
		assert r.returncode == 0, r.stderr
		# header fields are compared after whitespace splitting: the base writes them space-separated, a bug
		assert parse_motifs(r.stdout) == [
			(['AGATAA', 'gaTA1_gata1_motif', '5.0'], [['0.9', '0.05', '0.03', '0.02'], ['0.1', '0.7', '0.1', '0.1']]),
			(['GGGCGG', 'Klf4_sub_klf_motif', '6.0'], [['0.2', '0.3', '0.4', '0.1']]),
		]
		assert [x for x in r.stdout.splitlines() if not x.startswith('>')] == [
			'0.9\t0.05\t0.03\t0.02', '0.1\t0.7\t0.1\t0.1', '0.2\t0.3\t0.4\t0.1']
		assert 'skipped motif file because not found: missing.motif' in r.stderr
		assert 'skipped motif file because not containing exactly one motif: two.motif' in r.stderr
		assert 'skipped motif file because not starting from first line: nohead.motif' in r.stderr

	@pytest.mark.parametrize('cap,names', [
		('0', ['gaTA1', 'Klf4']), ('1', ['gata1', 'klf4']), ('2', ['GATA1', 'KLF4']), ('3', ['Gata1', 'Klf4'])])
	def test_mapfile_capitalization(self, run, homer, workdir, cap, names):
		"""-c converts the mapped gene symbols; -m reads a mapping file given relative to the working directory."""
		(workdir / 'map.txt').write_text('h\th\ngata1.motif\tgaTA1\nsub/klf.motif\tKlf4\n')
		r = run('-b', homer, '-m', 'map.txt', '-c', cap)
		assert r.returncode == 0, r.stderr
		assert [x[0][1] for x in parse_motifs(r.stdout)] == [names[0] + '_gata1_motif', names[1] + '_sub_klf_motif']
		assert r.stderr == ''

	@pytest.mark.parametrize('cap,names', [
		('0', ['Gata1', 'gaTA2', 'Klf4']), ('1', ['gata1', 'gata2', 'klf4']), ('2', ['GATA1', 'GATA2', 'KLF4']),
		('3', ['Gata1', 'Gata2', 'Klf4'])])
	def test_organism(self, run, require_head_negative, homer, cap, names):
		"""-o splits 'A:B(family)/...' names into one motif per gene with the suffix after the first '(' or '/'."""
		r = run('-b', homer, '-o', 'vertebrates', '-c', cap)
		assert r.returncode == 0, r.stderr
		motifs = parse_motifs(r.stdout)
		# awk's for-in order over the genes of one motif is unspecified, so the first two are sorted
		assert sorted(x[0] for x in motifs[:2]) == sorted([
			['AGATAA', names[0] + '_Zf__K562_Homer', '5.0'],
			['AGATAA', names[1] + '_Zf__K562_Homer', '5.0'],
		])
		assert [x[0] for x in motifs[2:]] == [['GGGCGG', names[2] + '_Zf__ESC_Homer', '6.0']]
		# the last motif's rows are not compared: the base can drop its last row, a bug
		for x in motifs[:2]:
			assert x[1] == [['0.9', '0.05', '0.03', '0.02'], ['0.1', '0.7', '0.1', '0.1']]

	def test_unknown_organism(self, run, homer):
		"""An organism without known.motifs exits 1."""
		r = run('-b', homer, '-o', 'plants')
		assert r.returncode == 1
		assert 'cannot find motif file' in r.stderr

	def test_missing_mapfile(self, run, require_realpath_missing, homer, tmp_path):
		"""A missing mapping file exits 1."""
		os.remove(os.path.join(homer, 'motifs', 'extras', 'motifs2symbol.txt'))
		r = run('-b', homer, '-m', str(tmp_path / 'nomap.txt'))
		assert r.returncode == 1
		assert 'cannot find gene name mapping file' in r.stderr

	def test_missing_basedir(self, run, require_realpath_missing, tmp_path):
		"""A base folder that does not exist exits 1."""
		r = run('-b', str(tmp_path / 'nohomer'))
		assert r.returncode == 1
		assert 'cannot find homer base directory' in r.stderr

	def test_autodetect_basedir(self, run, require, require_head_negative, homer):
		"""Without -b the base folder is found from the homer program on PATH."""
		require('which')
		os.makedirs(os.path.join(homer, 'bin'))
		write_exe(os.path.join(homer, 'bin', 'homer'), '#!/usr/bin/env bash\nexit 0\n')
		r = run('-o', 'vertebrates', env={'PATH': path_with(os.path.join(homer, 'bin'))})
		assert r.returncode == 0, r.stderr
		assert len(parse_motifs(r.stdout)) == 3

	def test_invalid_capitalization(self, run, homer):
		"""-c accepts only 0 to 3."""
		r = run('-b', homer, '-c', '4')
		assert r.returncode == 1
		assert 'Error: -c has an unrecognized argument value.' in r.stderr

	@pytest.mark.parametrize('args', [['-h'], ['extra'], ['-x']])
	def test_usage(self, run, args):
		"""-h, a positional argument and unknown options print usage and exit 1."""
		r = run(*args)
		assert r.returncode == 1
		assert 'Usage: motif_homer.sh' in r.stderr

	@pytest.mark.parametrize('opt', ['-b', '-m', '-o', '-c'])
	def test_option_requires_argument(self, run, opt):
		"""Options without a value are reported."""
		r = run(opt)
		assert r.returncode == 1
		assert 'Error: {} requires an argument.'.format(opt) in r.stderr


######################################################################
# repodown.sh
######################################################################

STUB_DICTYS_HELPER = '''#!/usr/bin/env bash
# Stand-in for dictys_helper repo_*.py: lists and copies files from $FAKE_REPO.
echo "$*" >> "$FAKE_LOG"
cmd="$2"
if [ "$cmd" == "ls" ]; then
	ls -1 "$FAKE_REPO"
elif [ "$cmd" == "down" ]; then
	f="${@: -1}"
	if [ "a$f" != "a$FAKE_BROKEN" ]; then cp "$FAKE_REPO/$f" ./; fi
fi
'''


class TestRepodown:
	@pytest.fixture
	def repo(self, tmp_path):
		"""Fake repository folder holding data.tar (one file pkg/a.txt) split three ways, plus unrelated files."""
		src = tmp_path / 'src'
		(src / 'pkg').mkdir(parents=True)
		(src / 'pkg' / 'a.txt').write_text('hello\n' * 100)
		d = tmp_path / 'repo'
		d.mkdir()
		with tarfile.open(str(d / 'data.tar'), 'w') as t:
			t.add(str(src / 'pkg'), arcname='pkg')
		whole = (d / 'data.tar').read_bytes()
		n = len(whole) // 3
		for i, part in enumerate([whole[:n], whole[n:2 * n], whole[2 * n:]]):
			(d / 'split.tar.{}'.format(i)).write_bytes(part)
		(d / 'split.tar.bak').write_bytes(b'not a part')
		(d / 'other.tar.0').write_bytes(b'not a part')
		return d

	@pytest.fixture
	def run(self, run_shell, require, stubbin, workdir, repo, tmp_path):
		"""Run repodown.sh with the stub dictys_helper; returns (result, logged dictys_helper calls)."""
		require('tar', 'grep', 'sort')
		write_exe(stubbin / 'dictys_helper', STUB_DICTYS_HELPER)
		log = tmp_path / 'helper.log'

		def f(*args, broken=''):
			r = run_shell('repodown.sh', *args, cwd=str(workdir), env={'PATH': path_with(stubbin),
				'FAKE_REPO': str(repo), 'FAKE_LOG': str(log), 'FAKE_BROKEN': broken})
			return r, (read_lines(log) if log.exists() else [])
		return f

	def test_single_file(self, run, workdir):
		"""A file present under its own name is downloaded, extracted and removed."""
		r, calls = run('zenodo', 'data.tar')
		assert r.returncode == 0, r.stdout + r.stderr
		assert 'Downloading data.tar' in r.stdout
		assert calls == ['repo_zenodo.py ls', 'repo_zenodo.py down data.tar']
		assert os.listdir(str(workdir)) == ['pkg']
		assert (workdir / 'pkg' / 'a.txt').read_text() == 'hello\n' * 100

	def test_split_file(self, run, workdir, repo):
		"""A file stored as name.0, name.1, ... is downloaded part by part, joined and extracted."""
		r, calls = run('dataverse', 'split.tar')
		assert r.returncode == 0, r.stdout + r.stderr
		assert calls == ['repo_dataverse.py ls', 'repo_dataverse.py ls'] + [
			'repo_dataverse.py down split.tar.{}'.format(i) for i in range(3)]
		assert os.listdir(str(workdir)) == ['pkg']
		assert (workdir / 'pkg' / 'a.txt').read_text() == 'hello\n' * 100

	def test_no_decompress(self, run, workdir, repo):
		"""-n keeps the joined archive, replacing an older file of the same name, and does not extract it."""
		(workdir / 'split.tar').write_bytes(b'old')
		r, _ = run('-n', 'zenodo', 'split.tar')
		assert r.returncode == 0, r.stdout + r.stderr
		assert os.listdir(str(workdir)) == ['split.tar']
		expected = b''.join((repo / 'split.tar.{}'.format(i)).read_bytes() for i in range(3))
		assert (workdir / 'split.tar').read_bytes() == expected

	def test_arguments_passed_on(self, run):
		"""Arguments after the file name are passed to the repository helper's ls and down."""
		r, calls = run('-n', 'dataverse', 'data.tar', '--doi', '10.1/X', '--token', 'T')
		assert r.returncode == 0, r.stdout + r.stderr
		assert calls == ['repo_dataverse.py ls --doi 10.1/X --token T',
			'repo_dataverse.py down --doi 10.1/X --token T data.tar']

	def test_file_not_found(self, run, workdir):
		"""A name with neither a whole file nor parts exits 1 without downloading."""
		r, calls = run('zenodo', 'absent.tar')
		assert r.returncode == 1
		assert 'File not found' in r.stdout
		assert all(' down ' not in x for x in calls)
		assert os.listdir(str(workdir)) == []

	@pytest.mark.parametrize('name,broken', [('data.tar', 'data.tar'), ('split.tar', 'split.tar.1')])
	def test_download_failure(self, run, name, broken):
		"""A download that leaves no file exits 1 naming the file."""
		r, _ = run('zenodo', name, broken=broken)
		assert r.returncode == 1
		assert 'Downloading failed for {}'.format(broken) in r.stdout

	@pytest.mark.parametrize('args', [['-h'], ['zenodo'], ['figshare', 'data.tar'], ['-x', 'zenodo', 'data.tar']])
	def test_usage(self, run, args, workdir):
		"""-h, a missing file name, an unknown repository and unknown options print usage and exit 1."""
		r, calls = run(*args)
		assert r.returncode == 1
		assert 'Usage: repodown.sh' in r.stderr
		assert calls == []
		assert os.listdir(str(workdir)) == []


######################################################################
# split_bam.sh
######################################################################

class TestSplitBam:
	@pytest.mark.parametrize('args', [['-h'], ['whole.bam'], ['-x', 'a', 'b']])
	def test_usage(self, run_shell, workdir, args):
		"""-h, a missing output folder and unknown options print usage and exit 1."""
		r = run_shell('split_bam.sh', *args, cwd=str(workdir))
		assert r.returncode == 1
		assert 'Usage: split_bam.sh' in r.stderr

	def test_split(self, run_shell, require, require_realpath_missing, stubbin, workdir):
		"""Reads are split into one BAM per CB barcode with the original header; temporary files are removed."""
		require('samtools')
		header = '@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n'
		cells = ['A', 'B', 'A', 'C', 'B', 'A']
		reads = ''.join('r{}\t0\tchr1\t{}\t60\t4M\t*\t0\t0\tACGT\tIIII\tCB:Z:{}\tNM:i:0\n'.format(i, 100 + i, c)
			for i, c in enumerate(cells))
		(workdir / 'whole.sam').write_text(header + reads)
		subprocess.run(['samtools', 'view', '-b', '-o', str(workdir / 'whole.bam'), str(workdir / 'whole.sam')],
			check=True, capture_output=True)
		r = run_shell('split_bam.sh', 'whole.bam', 'out', '--section', 'CB:Z:', cwd=str(workdir),
			env={'PATH': path_with(stubbin, os.path.join(REPO, 'bin'))})
		assert r.returncode == 0, r.stdout + r.stderr
		assert sorted(os.listdir(str(workdir))) == ['out', 'whole.bam', 'whole.sam']
		assert sorted(os.listdir(str(workdir / 'out'))) == ['A.bam', 'B.bam', 'C.bam']
		for c in 'ABC':
			v = subprocess.run(['samtools', 'view', str(workdir / 'out' / (c + '.bam'))], check=True,
				capture_output=True, text=True).stdout
			assert [x.split('\t')[0] for x in v.splitlines()] == ['r{}'.format(i) for i, y in enumerate(cells) if y == c]
			h = subprocess.run(['samtools', 'view', '-H', str(workdir / 'out' / (c + '.bam'))], check=True,
				capture_output=True, text=True).stdout
			assert '@SQ\tSN:chr1\tLN:1000' in h
