#!/usr/bin/python3
# Lingfei Wang, 2022. All rights reserved.

"""
Dynamic network functions
"""

################################################################
# Split cells into subsets along trajectory
################################################################

def subsets_rna(fi_traj:str,fi_traj_cell_rna:str,fi_coord_rna:str,fo_subsets:str,fo_subset_locs:str,diro_subsets:str,fo_subset_edges:str,ncell:int,noverlap:int,dmax:float)->None:
	"""
	Split RNA cells into overlapping subsets along trajectory.

	Parameters
	----------
	fi_traj:
		Path of input h5 file of trajectory
	fi_traj_cell_rna:
		Path of input h5 file of RNA cell locations on trajectory
	fi_coord_rna:
		Path of input tsv file of RNA cell locations on low dimension. Used for cell names
	fo_subsets:
		Path of output txt file of subset names
	fo_subset_locs:
		Path of output h5 file of subset locations on trajectory
	diro_subsets:
		Path of output folder containing one subfolder for each cell subset. This function adds names_rna.txt to each folder for cell names included.
	fo_subset_edges:
		Path of output tsv file for neighborhood relations between cell subsets
	ncell:
		Number of cells in each subset
	noverlap:
		Expected number of cell overlap between neighboring subsets
	dmax:
		Soft upper bound of pseudotime distance between neighboring subsets
	"""
	#Note: this function is not tested in continuous integration
	import logging
	from os import makedirs,linesep
	from os.path import join as pjoin
	import numpy as np
	import pandas as pd
	from dictys.traj import trajectory,point
	#Load data
	traj=trajectory.from_file(fi_traj)
	points=point.from_file(traj,fi_traj_cell_rna)
	logging.info(f'Reading file {fi_coord_rna}')
	names_cell=np.array(list(pd.read_csv(fi_coord_rna,header=0,index_col=0,sep='\t').index))
	#Construct subsets
	ans=points.subsets(ncell,noverlap,dmax)
	#Formatting output
	n=len(ans[0])
	subsets=point(traj,ans[0],ans[1])
	subset_cells=[linesep.join(names_cell[x])+linesep for x in ans[2].T]
	names_subset=ans[4]
	nodegraph=pd.DataFrame(ans[3],index=names_subset,columns=names_subset)
	#Output subset files
	for xi in range(n):
		makedirs(pjoin(diro_subsets,names_subset[xi]),exist_ok=True)
		fo=pjoin(diro_subsets,names_subset[xi],'names_rna.txt')
		logging.info(f'Writing file {fo}')
		with open(fo,'w') as f:
			f.write(subset_cells[xi])
	#Output single files
	subsets.to_file(fo_subset_locs)
	logging.info(f'Writing file {fo_subset_edges}')
	nodegraph.to_csv(fo_subset_edges,header=True,index=True,sep='\t')
	logging.info(f'Writing file {fo_subsets}')
	with open(fo_subsets,'w') as f:
		f.write(linesep.join(names_subset)+linesep)	

def subsets_atac(fi_traj:str,fi_traj_cell_atac:str,fi_coord_atac:str,fi_subsets:str,fi_subset_locs:str,fo_subset:str,subset_name:str,ncell:int)->None:
	"""
	Split ATAC cells into overlapping subsets along trajectory.

	Parameters
	----------
	fi_traj:
		Path of input h5 file of trajectory
	fi_traj_cell_atac:
		Path of input h5 file of ATAC cell locations on trajectory
	fi_coord_atac:
		Path of input tsv file of ATAC cell locations on low dimension. Used for cell names
	fi_subsets:
		Path of output txt file of subset names
	fi_subset_locs:
		Path of input h5 file of subset locations on trajectory
	fo_subset:
		Path of output txt file for ATAC cell subset names
	subset_name:
		Name of subset to find ATAC cell subset
	ncell:
		Number of cells in subset
	"""
	import logging
	from os import linesep
	import numpy as np
	import pandas as pd
	from dictys.traj import trajectory,point
	from dictys.utils.file import read_txt
	# Load data
	traj=trajectory.from_file(fi_traj)
	points=point.from_file(traj,fi_traj_cell_atac)
	logging.info(f'Reading file {fi_coord_atac}')
	names_cell=np.array(list(pd.read_csv(fi_coord_atac,header=0,index_col=0,sep='\t').index))
	logging.info(f'Reading file {fi_subsets}')
	names_subset=read_txt(fi_subsets,unique=True)
	dict_subset=dict(zip(names_subset,range(len(names_subset))))
	subsets=point.from_file(traj,fi_subset_locs)
	subsets=subsets[[dict_subset[subset_name]]]
	# Find nearest ncell cells for each subset
	points.perturb()
	dist=(subsets-points).ravel()
	t1=np.argpartition(dist,ncell)[:ncell]
	subset_cells=linesep.join(names_cell[t1])+linesep
	# Output subset files
	logging.info(f'Writing file {fo_subset}')
	with open(fo_subset,'w') as f:
		f.write(subset_cells)



















































#
