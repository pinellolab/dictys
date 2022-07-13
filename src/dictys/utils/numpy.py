#!/usr/bin/python3
# Lingfei Wang, 2020-2022. All rights reserved.

"""Utility function for numpy"""

from typing import Any,Optional

try:
	import numpy.typing as npt
	NDArray=npt.NDArray
except ModuleNotFoundError:
	NDArray=list
try:
	import numpy.typing as npt
	ArrayLike=npt.ArrayLike
except ModuleNotFoundError:
	ArrayLike=list

def groupby(array:NDArray)->dict[Any,NDArray]:
	"""
	Groups 1-D array to dict(value=np.array(index))

	Parameters
	----------
	array:	numpy.ndarray
		Array to group

	Returns
	-------
	{value:numpy.ndarray}
		Dictionary to convert value to indices
	"""
	from collections import defaultdict
	import numpy as np
	d=defaultdict(list)
	for xi in range(len(array)):
		d[array[xi]].append(xi)
	d=dict(d)
	d={x:np.array(y) for x,y in d.items()}
	return d

def median(v:NDArray,w:Optional[NDArray]=None,**ka)->Any:
	"""
	Weighted median.

	Parameters
	----------
	v:		numpy.ndarray in 1-dimension
		Array of values to find median
	w:		numpy.ndarray in 1-dimension or None
		Weight of each value in the same shape as `v`. If None, use numpy.median.
	ka:		dict
		Keyworad arguments passed to numpy.median if w is None.

	Returns
	-------
	any
		Median value
	"""
	import numpy as np
	if w is None:
		return np.median(v,**ka)
	if len(ka)>0:
		raise TypeError('Weighted median does not accept any other keyword argument.')
	assert v.ndim==1 and v.shape==w.shape
	assert np.isfinite(v).all() and np.isfinite(w).all()
	assert (w>=0).all()
	t1=v.argsort()
	v,w=[x[t1] for x in [v,w]]
	ws=w.cumsum()
	return v[ws.searchsorted(ws[-1]/2)]

def dtype_min(v:NDArray):
	"""
	Find the minimum sized numpy.dtype for integer array.

	Parameters
	----------
	v:
		Array to find smallest size

	Returns
	-------
	numpy.dtype
		Minimum dtype found.
	"""
	import numpy as np
	if np.issubdtype(v.dtype,np.signedinteger):
		typebase='i'
	elif np.issubdtype(v.dtype,np.unsignedinteger):
		typebase='u'
	else:
		assert not np.issubdtype(v.dtype,np.integer)
		raise TypeError('v must be numpy.integer type.')
	itemsize=v.dtype.itemsize
	vrange=[v.min(),v.max()]
	while itemsize>0:
		itemsize_new=itemsize//2
		t1=np.iinfo(np.dtype(f'{typebase}{itemsize_new}'))
		if t1.min<=vrange[0] and t1.max>=vrange[1]:
			itemsize=itemsize_new
		else:
			break
	assert itemsize>0
	return np.dtype(f'{typebase}{itemsize}')


assert __name__ != "__main__"


































#
