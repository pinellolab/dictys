"""
Command line executor
"""

from docstring2argparse import docstringrunner
import logging
logging.basicConfig(
	format=
	'%(levelname)s:%(process)d:%(asctime)s:%(pathname)s:%(lineno)d:%(message)s',
	# level=logging.WARNING)
	level=logging.DEBUG)
docstringrunner(__package__)

# if len(sys.argv) == 1:
# 	p[pkgname].print_help(sys.stderr)
# 	sys.exit(1)
# args=p[pkgname].parse_args()
# args=vars(args)
# logging.basicConfig(
# 	format=
# 	'%(levelname)s:%(process)d:%(asctime)s:%(pathname)s:%(lineno)d:%(message)s',
# 	level=logging.DEBUG if args['verbose'] else logging.WARNING)
# from . import run
# func = getattr(run, args['cmd'])
# func(args)









#
