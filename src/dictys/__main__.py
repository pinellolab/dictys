"""
Command line executor
"""

from docstring2argparse import docstringrunner
import logging
logging.basicConfig(format='%(levelname)s:%(process)d:%(asctime)s:%(pathname)s:%(lineno)d:%(message)s',level=logging.WARNING)
# level=logging.DEBUG)

if __name__=="__main__":
	docstringrunner(__package__)










#
