"""
Command line executor
"""

import logging
from docstring2argparse import docstringrunner

logging.basicConfig(format='%(levelname)s:%(process)d:%(asctime)s:%(pathname)s:%(lineno)d:%(message)s',
	level=logging.WARNING)
# 	level=logging.DEBUG)

if __name__=="__main__":
	docstringrunner(__package__)










#
