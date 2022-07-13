# Lingfei Wang, 2022. All rights reserved.
"""
File I/O
"""

def read_txt(fi:str,unique:bool=False)->list[str]:
	"""
	Reading single text file to a list of str. Removing empty lines and heading & trailing spaces.

	Parameters
	----------
	fi:
		Path of input txt file

	Returns
	-------
	List of lines

	"""
	import logging
	logging.info(f'Reading file {fi}')
	with open(fi,'r') as f:
		ans=f.readlines()
	ans=[x.strip() for x in ans]
	ans=list(filter(lambda x:len(x)>0,ans))
	assert not unique or len(ans)==len(set(ans))
	return ans

def read_columns(fi:str,header=0,index_col=0,sep='\t',unique=False,**ka)->list[str]:
	"""
	Reading column names from table file.

	Parameters
	----------
	fi:
		Path of input table file
	others:
		Keyword arguments for pandas.read_csv

	Returns
	-------
	List of column names

	"""
	import logging
	import pandas as pd
	logging.info(f'Reading file {fi}')
	d=pd.read_csv(fi,header=header,index_col=index_col,sep=sep,nrows=1,**ka)
	ans=list(str(x) for x in d.columns)
	assert not unique or len(ans)==len(set(ans))
	return ans

def read_index(fi:str,header=0,index_col=0,sep='\t',unique=False,**ka)->list[str]:
	"""
	Reading row names from table file.

	Parameters
	----------
	fi:
		Path of input table file
	others:
		Keyword arguments for pandas.read_csv

	Returns
	-------
	List of row names

	"""
	import logging
	import pandas as pd
	usecols=[0]
	if index_col not in usecols:
		raise NotImplementedError
	logging.info(f'Reading file {fi}')
	d=pd.read_csv(fi,header=header,index_col=index_col,sep=sep,usecols=usecols,**ka)
	ans=list(str(x) for x in d.index)
	assert not unique or len(ans)==len(set(ans))
	return ans





































assert __name__ != "__main__"
