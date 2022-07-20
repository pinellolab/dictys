#!/usr/bin/python3
# Lingfei Wang, 2022. All rights reserved.

"""Utility function for networkx"""

try:
	from networkx.drawing.layout import random_state		# noqa: F401
except ImportError:
	from networkx.drawing.layout import np_random_state as random_state		# noqa: F401

assert __name__ != "__main__"


































#
