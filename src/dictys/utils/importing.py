#!/usr/bin/python3
# Lingfei Wang, 2020-2022. All rights reserved.

"""
Pretending to import unimportable modules when not needed
"""
try:
	ModuleNotFoundError
except NameError:
	ModuleNotFoundError = Exception		# pylint: disable=W0622

class dummy_class:
	"""
	Dummy class for unimportable modules
	"""
	def __getattr__(self,key):
		if key=='__mro_entries__':
			return lambda *a,**ka:tuple()
		return self
	def __call__(self,*a,**ka):
		return self


try:
	import torch
except ModuleNotFoundError:		# pylint: disable=W0703
	torch=dummy_class()
try:
	import pyro
except ModuleNotFoundError:		# pylint: disable=W0703
	pyro=dummy_class()
try:
	import matplotlib
	import matplotlib.pyplot as matplotlib_pyplot
	import matplotlib.figure as matplotlib_figure
	matplotlib.pyplot=matplotlib_pyplot
	matplotlib.figure=matplotlib_figure
except ModuleNotFoundError:		# pylint: disable=W0703
	matplotlib=dummy_class()

assert __name__ != "__main__"


































#
