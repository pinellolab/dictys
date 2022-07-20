#!/usr/bin/python3
# Michael Hu, 2021, 2022. All rights reserved.

"""
Statistics of networks
"""
import networkx as nx
from dictys.utils.networkx import random_state

def _fruchterman_reingold(
	A,pos, k=None, fixed=None, iterations=50, threshold=1e-4, dim=2, pretendIterations = None, stop = None
):
	"""
	Modified _fruchterman_reingold function from networkx.
	Adds 
	"""
	import numpy as np

	try:
		nnodes, _ = A.shape
	except AttributeError as e:
		msg = "fruchterman_reingold() takes an adjacency matrix as input"
		raise nx.NetworkXError(msg) from e
	
	if pretendIterations is None:
		pretendIterations = iterations
		
	if stop is None:
		stop = iterations

	# make sure positions are of same type as matrix
	pos = pos.astype(A.dtype)

	# optimal distance between nodes
	if k is None:
		k = np.sqrt(1.0 / nnodes)
	# the initial "temperature"  is about .1 of domain area (=1x1)
	# this is the largest step allowed in the dynamics.
	# We need to calculate this in case our fixed positions force our domain
	# to be much bigger than 1x1
	if (iterations != 0):
		t = max(max(pos.T[0]) - min(pos.T[0]), max(pos.T[1]) - min(pos.T[1])) * 0.1 * (float(pretendIterations)/iterations)
	# simple cooling scheme.
#     linearly step down by dt on each iteration so last iteration is size dt.
		dt = t / float(iterations + 1)
	delta = np.zeros((pos.shape[0], pos.shape[0], pos.shape[1]), dtype=A.dtype)
	# the inscrutable (but fast) version
	# this is still O(V^2)
	# could use multilevel methods to speed this up significantly
	for iteration in range(stop):
		# matrix of difference between points
		delta = pos[:, np.newaxis, :] - pos[np.newaxis, :, :]
		# distance between points
		distance = np.linalg.norm(delta, axis=-1)
		# enforce minimum distance of 0.01
		np.clip(distance, 0.01, None, out=distance)
		# displacement "force"
		displacement = np.einsum(
			"ijk,ij->ik", delta, (k * k / distance ** 2 - A * distance / k)
		)
		# update positions
		length = np.linalg.norm(displacement, axis=-1)
		length = np.where(length < 0.01, 0.1, length)
		delta_pos = np.einsum("ij,i->ij", displacement, t / length)
		if fixed is not None:
			# don't change positions of fixed nodes
			delta_pos[fixed] = 0.0
		pos += delta_pos
		# cool temperature
		t -= dt
		err = np.linalg.norm(delta_pos) / nnodes
		if err < threshold*(float(pretendIterations)/iterations):
			break
	return pos

def _sparse_fruchterman_reingold(
	A, pos, k=None, fixed=None, iterations=50, threshold=1e-4, dim=2, pretendIterations = None, stop = None
):
	# Position nodes in adjacency matrix A using Fruchterman-Reingold
	# Entry point for NetworkX graph is fruchterman_reingold_layout()
	# Sparse version
	import numpy as np
	import scipy as sp
	import scipy.sparse  # call as sp.sparse

	try:
		nnodes, _ = A.shape
	except AttributeError as e:
		msg = "fruchterman_reingold() takes an adjacency matrix as input"
		raise nx.NetworkXError(msg) from e
	
	if pretendIterations is None:
		pretendIterations = iterations

	if stop is None:
		stop = iterations

	# make sure we have a LIst of Lists representation
	try:
		A = A.tolil()
	except AttributeError:
		A = (sp.sparse.coo_matrix(A)).tolil()

	# make sure positions are of same type as matrix
	pos = pos.astype(A.dtype)

	# no fixed nodes
	if fixed is None:
		fixed = []

	# optimal distance between nodes
	if k is None:
		k = np.sqrt(1.0 / nnodes)
	# the initial "temperature"  is about .1 of domain area (=1x1)
	# this is the largest step allowed in the dynamics.
	if (iterations != 0):
		t = max(max(pos.T[0]) - min(pos.T[0]), max(pos.T[1]) - min(pos.T[1])) * 0.1 * (float(pretendIterations)/iterations)
	# simple cooling scheme.
	# linearly step down by dt on each iteration so last iteration is size dt.
	dt = t / float(iterations + 1)

	displacement = np.zeros((dim, nnodes))
	for iteration in range(stop):
		displacement *= 0
		# loop over rows
		for i in range(A.shape[0]):
			if i in fixed:
				continue
			# difference between this row's node position and all others
			delta = (pos[i] - pos).T
			# distance between points
			distance = np.sqrt((delta ** 2).sum(axis=0))
			# enforce minimum distance of 0.01
			distance = np.where(distance < 0.01, 0.01, distance)
			# the adjacency matrix row
			Ai = np.asarray(A.getrowview(i).toarray())
			# displacement "force"
			displacement[:, i] += (
				delta * (k * k / distance ** 2 - Ai * distance / k)
			).sum(axis=1)
		# update positions
		length = np.sqrt((displacement ** 2).sum(axis=0))
		length = np.where(length < 0.01, 0.1, length)
		delta_pos = (displacement * t / length).T
		pos += delta_pos
		# cool temperature
		t -= dt
		err = np.linalg.norm(delta_pos) / nnodes
		if err < threshold*(float(pretendIterations)/iterations):
			break
	return pos

@random_state(10)
def fruchterman_reingold_layout(
	G,
	k=None,
	pos=None,
	fixed=None,
	iterations=1000,
	threshold=1e-4,
	weight="weight",
	scale=1,
	center=None,
	dim=2,
	seed=None,
	pretendIterations=50,
	stop = None
):
	"""
	This is an modified version of the Fruchterman-Reingold layout function. (from networkx import fruchterman_reingold_layout)
	This function allows the user to get layout information before it fully converges with the "stop" parameter
	The function also allows for finer steps for each iteration using "iterations" and "pretendIterations"
	pretendIterations (int) = the max number of iterations run in the original function; the ratio of pretendIterations/iterations is the ratio of step sizes between the edited and original functions
	
	stop (int) = the number of iterations run before the returned layout stops prematurely
	Documention below is from the original FR function:
	Position nodes using Fruchterman-Reingold force-directed algorithm.
	The algorithm simulates a force-directed representation of the network
	treating edges as springs holding nodes close, while treating nodes
	as repelling objects, sometimes called an anti-gravity force.
	Simulation continues until the positions are close to an equilibrium.
	There are some hard-coded values: minimal distance between
	nodes (0.01) and "temperature" of 0.1 to ensure nodes don't fly away.
	During the simulation, `k` helps determine the distance between nodes,
	though `scale` and `center` determine the size and place after
	rescaling occurs at the end of the simulation.
	Fixing some nodes doesn't allow them to move in the simulation.
	It also turns off the rescaling feature at the simulation's end.
	In addition, setting `scale` to `None` turns off rescaling.
	Parameters
	----------
	G : NetworkX graph or list of nodes
		A position will be assigned to every node in G.
	k : float (default=None)
		Optimal distance between nodes.  If None the distance is set to
		1/sqrt(n) where n is the number of nodes.  Increase this value
		to move nodes farther apart.
	pos : dict or None  optional (default=None)
		Initial positions for nodes as a dictionary with node as keys
		and values as a coordinate list or tuple.  If None, then use
		random initial positions.
	fixed : list or None  optional (default=None)
		Nodes to keep fixed at initial position.
		ValueError raised if `fixed` specified and `pos` not.
	iterations : int  optional (default=50)
		Maximum number of iterations taken
	threshold: float optional (default = 1e-4)
		Threshold for relative error in node position changes.
		The iteration stops if the error is below this threshold.
	weight : string or None   optional (default='weight')
		The edge attribute that holds the numerical value used for
		the edge weight.  If None, then all edge weights are 1.
	scale : number or None (default: 1)
		Scale factor for positions. Not used unless `fixed is None`.
		If scale is None, no rescaling is performed.
	center : array-like or None
		Coordinate pair around which to center the layout.
		Not used unless `fixed is None`.
	dim : int
		Dimension of layout.
	seed : int, RandomState instance or None  optional (default=None)
		Set the random state for deterministic node layouts.
		If int, `seed` is the seed used by the random number generator,
		if numpy.random.RandomState instance, `seed` is the random
		number generator,
		if None, the random number generator is the RandomState instance used
		by numpy.random.
	Returns
	-------
	pos : dict
		A dictionary of positions keyed by node
	Examples
	--------
	>>> G = nx.path_graph(4)
	>>> pos = nx.spring_layout(G)
	# The same using longer but equivalent function name
	>>> pos = nx.fruchterman_reingold_layout(G)
	"""
	import numpy as np
	from networkx.drawing.layout import _process_params,rescale_layout

	G, center = _process_params(G, center, dim)
	if fixed is not None:
		if pos is None:
			raise ValueError("nodes are fixed without positions given")
		for node in fixed:
			if node not in pos:
				raise ValueError("nodes are fixed without positions given")
		nfixed = {node: i for i, node in enumerate(G)}
		fixed = np.asarray([nfixed[node] for node in fixed])

	if pos is not None:
		# Determine size of existing domain to adjust initial positions
		dom_size = max(coord for pos_tup in pos.values() for coord in pos_tup)
		if dom_size == 0:
			dom_size = 1
		pos_arr = seed.rand(len(G), dim) * dom_size + center

		for i, n in enumerate(G):
			if n in pos:
				pos_arr[i] = np.asarray(pos[n])
	else:
		pos_arr = None
		dom_size = 1

	if len(G) == 0:
		return {}
	if len(G) == 1:
		return {nx.utils.arbitrary_element(G.nodes()): center}
	A = nx.to_numpy_array(G, weight=weight)
	if k is None and fixed is not None:
		# We must adjust k by domain size for layouts not near 1x1
		nnodes, _ = A.shape
		k = dom_size / np.sqrt(nnodes) 
#             what is k and dom_size
	pos = _fruchterman_reingold(
		A, k, pos_arr, fixed, iterations, threshold, dim, seed, pretendIterations, stop
	)
	if fixed is None and scale is not None:
		pos = rescale_layout(pos, scale=scale) + center
	pos = dict(zip(G, pos))
	return pos





































#
