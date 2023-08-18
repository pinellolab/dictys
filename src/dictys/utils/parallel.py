#!/usr/bin/python3
# Lingfei Wang, 2020-2022. All rights reserved.

"""
Module for parallel processing
"""

from typing import Optional
from contextlib import contextmanager

def autocount(n:float)->int:
	"""Return the number of threads based on parameter n. Possible values:
	n=0: 			All CPUs
	n=other int:	n CPUs
	1<n<0:			CPU count*n
	"""
	from multiprocessing import cpu_count
	import logging
	assert n >= 0
	if n == 0:
		n = cpu_count()
	elif n < 1:
		n = max(1, int(n * cpu_count()))
	else:
		n = int(n)
	logging.info('Using {} threads'.format(n))
	return n

def set_num_threads(n:int)->dict[str,Optional[str]]:
	"""Sets the number of threads in numerical calculations of external libraries through environmental variables.
	Parameters:
	n:	Number of threads allowed in external libraries
	Return:
	Dictionary of previous environmental values to allow recovery with function recover_num_threads.
	"""
	import os
	keys = 'OMP_NUM_THREADS,MKL_NUM_THREADS,NUMEXPR_NUM_THREADS,OPENBLAS_NUM_THREADS,OMP_MAX_THREADS,MKL_MAX_THREADS,NUMEXPR_MAX_THREADS,OPENBLAS_MAX_THREADS,VECLIB_MAXIMUM_THREADS'.split(',')

	n = str(n)
	ans = {}
	for xi in keys:
		if xi not in os.environ:
			ans[xi] = None
		else:
			ans[xi] = os.environ[xi]
		os.environ[xi] = n
	return ans

def recover_num_threads(d:dict[str,str])->None:
	"""Recovers the previous environmental variables changed by function set_num_threads.
	Parameters:
	d:	Previous environmental values before setting with set_num_threads. Also return of set_num_threads.
	"""
	import os
	import logging

	t1 = list(filter(lambda x: x not in os.environ, d))
	if len(t1) > 0:
		logging.warning('Environmental variables not found: ' + ','.join(t1))
	t1 = list(filter(lambda x: x in os.environ, d))
	for xi in t1:
		if d[xi] is None:
			del os.environ[xi]
		else:
			os.environ[xi] = d[xi]

@contextmanager
def num_threads(n:float):
	"""
	Context manager for controlling CPU thread count.
	"""
	from threadpoolctl import threadpool_limits
	from joblib import parallel_backend
	n=autocount(n)
	assert isinstance(n,int) and n>=1
	n0=set_num_threads(n)
	try:
		with parallel_backend('threading', n_jobs=n):
			with threadpool_limits(limits=n):
				yield
	finally:
		recover_num_threads(n0)


assert __name__ != "__main__"







































#
