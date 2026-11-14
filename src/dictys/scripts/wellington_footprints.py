#!/usr/bin/env python3
# Footprinting of peak regions with pyDNase's Wellington algorithm.
# Same command line and 'p value cutoffs' output as pyDNase's own
# wellington_footprints.py, minus what dictys never reads: the null
# randomisations of the FDR cutoff (-fdr, -fdrlimit and -fdriter are accepted
# and ignored), the WIG of all scores and the FDR bed.
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


def calculate(fp):
	"""Score one region without the FDR null (worker side)."""
	fp.storedscore, fp.storedlengths = fp.calculate()
	return fp


def main(argv=None):
	import pyDNase
	from pyDNase import footprinting
	args = parse_args(argv)
	shoulders = xrange_from_string(args.shoulder_sizes)
	footprints = xrange_from_string(args.footprint_sizes)
	assert len([x for x in os.listdir(args.outputdir) if x[0] != '.']) == 0, 'output directory {} is not empty!'.format(args.outputdir)
	prefix = args.output_prefix or os.path.basename(args.reads) + '.' + os.path.basename(args.regions)

	regions = pyDNase.GenomicIntervalSet(args.regions)
	reads = pyDNase.BAMHandler(args.reads, caching=False, ATAC=args.A)
	pvdir = os.path.join(args.outputdir, 'p value cutoffs')
	os.makedirs(pvdir)
	ofiles = {x: open(os.path.join(pvdir, prefix + '.WellingtonFootprints.{}.bed'.format(x)), 'a') for x in args.pv_cutoffs}

	def write(fp):
		for x, f in ofiles.items():
			for footprint in fp.footprints(withCutoff=x):
				print(footprint, file=f)

	ordered = [x for y in sorted(regions.intervals.values()) for x in y]
	ncpu = args.p if args.p else mp.cpu_count()
	cached = 10 * ncpu
	cls = footprinting.wellington1D if args.one_dimension else footprinting.wellington
	pool = mp.Pool(ncpu)
	for i in ordered:
		fp = cls(i, reads, shoulder_sizes=shoulders, footprint_sizes=footprints, FDR_cutoff=args.FDR_cutoff, FDR_iterations=args.FDR_iterations, bonferroni=args.bonferroni)
		pool.apply_async(calculate, (fp,), callback=write)
		while pool._taskqueue.qsize() > cached:
			pass
	pool.close()
	pool.join()
	for f in ofiles.values():
		f.close()


if __name__ == '__main__':
	main()
