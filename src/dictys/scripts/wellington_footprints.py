#!/usr/bin/env python3
# Footprinting of peak regions with pyDNase's Wellington algorithm.
# Same command line and 'p value cutoffs' output as pyDNase's own
# wellington_footprints.py, minus what dictys never reads: the null
# randomisations of the FDR cutoff (-fdr, -fdrlimit and -fdriter are accepted
# and ignored), the WIG of all scores and the FDR bed.
#
# Adapted from pyDNase 0.3.0 (https://github.com/jpiper/pyDNase), under the MIT license
# below. Adapted parts: the command line (parse_args), xrange_from_string, the loop over
# regions and the output of footprints in main follow pyDNase's
# scripts/wellington_footprints.py; get_cut_values in make_handler repeats the read loop of
# BAMHandler.__lookupReadsWithoutCache in pyDNase/__init__.py. The changes are the removal
# of the FDR null, the WIG and the FDR bed, and per-region cut counts built with numpy in
# pool workers.
#
# Copyright 2018 Jason Piper
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
import argparse
import multiprocessing as mp
import os
import sys


def parse_args(argv=None):
	parser = argparse.ArgumentParser(description='Footprint regions in a DNase-seq or ATAC-seq experiment with the Wellington algorithm.')
	parser.add_argument('-b', '--bonferroni', action='store_true', default=False)
	parser.add_argument('-sh', '--shoulder-sizes', default='35,36,1', type=str, help='from,to,step')
	parser.add_argument('-fp', '--footprint-sizes', default='11,26,2', type=str, help='from,to,step')
	parser.add_argument('-d', '--one_dimension', action='store_true', default=False)
	parser.add_argument('-fdr', '--FDR_cutoff', default=0.01, type=float, help='Ignored')
	parser.add_argument('-fdriter', '--FDR_iterations', default=100, type=int, help='Ignored')
	parser.add_argument('-fdrlimit', '--FDR_limit', default=-20, type=int, help='Ignored')
	parser.add_argument('-pv', '--pv_cutoffs', default=[-10, -20, -30, -40, -50, -75, -100, -300, -500, -700], type=int, nargs='+')
	parser.add_argument('-dm', '--dont-merge-footprints', action='store_true', default=False)
	parser.add_argument('-o', '--output_prefix', default='', type=str)
	parser.add_argument('-p', default=1, type=int, help='Number of processes, 0 for all cores')
	parser.add_argument('-A', action='store_true', default=False, help='ATAC-seq mode')
	parser.add_argument('regions', help='BED file of the regions to footprint')
	parser.add_argument('reads', help='BAM file of the reads')
	parser.add_argument('outputdir', help='Empty writeable directory for the results')
	return parser.parse_args(argv)


def xrange_from_string(s):
	try:
		a, b, c = map(int, s.split(','))
		ans = list(range(a, b, c))
		assert len(ans) > 0
		return ans
	except Exception:
		raise ValueError('shoulder and footprint sizes must be supplied as from,to,step')


def make_handler(pyDNase, reads_path, atac):
	"""BAMHandler whose per-region cut counts are built with numpy a dict lookup per base (WellingtonC takes lists)."""
	import numpy as np

	class Handler(pyDNase.BAMHandler):
		def get_cut_values(self, vals):
			if isinstance(vals, str):
				return super().get_cut_values(vals)
			chrom, startbp, endbp, flip = vals.chromosome, vals.startbp, vals.endbp, vals.strand
			fw = []
			rv = []
			for r in self.samfile.fetch(chrom, max(startbp, 0), endbp):
				if not r.is_unmapped:
					if r.is_reverse:
						a = int(r.reference_end)
						if a <= endbp + 1:
							rv.append(a)
					else:
						a = int(r.reference_start) - 1
						if a >= startbp:
							fw.append(a)
			out = {}
			for key, pos, off in (('+', fw, self.loffset), ('-', rv, self.roffset)):
				arr = np.zeros(endbp - startbp, dtype=np.int64)
				if pos:
					idx = np.array(pos, dtype=np.int64) - (startbp + off)
					idx = idx[(idx >= 0) & (idx < arr.size)]
					arr += np.bincount(idx, minlength=arr.size)
				out[key] = arr.tolist()
			if flip == '-':
				out['+'], out['-'] = out['-'][::-1], out['+'][::-1]
			return out

	return Handler(reads_path, caching=False, ATAC=atac)


_worker = {}


def init_worker(reads_path, atac, cls_name, shoulders, footprints, fdr_cutoff, fdr_iterations, bonferroni):
	import pyDNase
	from pyDNase import footprinting
	try:
		_worker['reads'] = make_handler(pyDNase, reads_path, atac)
	except Exception as e:
		# An initializer that raises makes Pool respawn workers forever; report the error from calculate instead
		_worker['error'] = e
	_worker['args'] = (getattr(footprinting, cls_name), shoulders, footprints, fdr_cutoff, fdr_iterations, bonferroni)


def calculate(region):
	"""Fetch the reads of one region and score it without the FDR null (worker side)."""
	import pyDNase
	if 'error' in _worker:
		raise _worker['error']
	chrom, start, end, strand = region
	cls, shoulders, footprints, fdr_cutoff, fdr_iterations, bonferroni = _worker['args']
	fp = cls(pyDNase.GenomicInterval(chrom, start, end, strand=strand), _worker['reads'], shoulder_sizes=shoulders, footprint_sizes=footprints, FDR_cutoff=fdr_cutoff, FDR_iterations=fdr_iterations, bonferroni=bonferroni)
	fp.storedscore, fp.storedlengths = fp.calculate()
	fp.reads = None
	return fp


def main(argv=None):
	import pyDNase
	args = parse_args(argv)
	shoulders = xrange_from_string(args.shoulder_sizes)
	footprints = xrange_from_string(args.footprint_sizes)
	assert len([x for x in os.listdir(args.outputdir) if x[0] != '.']) == 0, 'output directory {} is not empty!'.format(args.outputdir)
	prefix = args.output_prefix or os.path.basename(args.reads) + '.' + os.path.basename(args.regions)

	regions = pyDNase.GenomicIntervalSet(args.regions)
	pvdir = os.path.join(args.outputdir, 'p value cutoffs')
	os.makedirs(pvdir)
	ofiles = {x: open(os.path.join(pvdir, prefix + '.WellingtonFootprints.{}.bed'.format(x)), 'a') for x in args.pv_cutoffs}

	def write(fp):
		for x, f in ofiles.items():
			for footprint in fp.footprints(withCutoff=x):
				print(footprint, file=f)

	ordered = [(i.chromosome, i.startbp, i.endbp, i.strand) for y in sorted(regions.intervals.values()) for i in y]
	ncpu = args.p if args.p else mp.cpu_count()
	cls_name = 'wellington1D' if args.one_dimension else 'wellington'
	initargs = (args.reads, args.A, cls_name, shoulders, footprints, args.FDR_cutoff, args.FDR_iterations, args.bonferroni)
	try:
		with mp.Pool(ncpu, init_worker, initargs) as pool:
			for fp in pool.imap(calculate, ordered, chunksize=16):
				write(fp)
	except Exception as e:
		sys.exit('wellington_footprints.py: cannot process {}: {}: {}'.format(args.reads, type(e).__name__, e))
	for f in ofiles.values():
		f.close()


if __name__ == '__main__':
	main()
