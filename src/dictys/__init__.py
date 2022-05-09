#!/usr/bin/python3
# Lingfei Wang, 2020-2022. All rights reserved.

"""
Dictys for dynamic gene regulatory network reconstruction from single-cell multi-omic data
"""

__all__=['dynamic','preproc','chromatin','network','net','plot','traj','utils']

_docstring2argparse_ignore_=['net','plot','traj','utils','curve']

from . import *

assert __name__ != "__main__"
