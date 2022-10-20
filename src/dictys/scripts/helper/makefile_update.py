#!/usr/bin/env python3
# Lingfei Wang, 2022. All rights reserved.
import json
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Updates makefile variable assignments with values provided in json string")
parser.add_argument('makefile_path',type=str,help='Path of makefile to update and rewritten.')
parser.add_argument('json_string',type=str,help='Update to be made in json format: {"variable_name":"new_value"}. Variable names can have "+" suffix to indicate appending to current value.')

args=parser.parse_args()
fmake=args.makefile_path
vs=args.json_string

#Load files
vs=json.loads(vs)
with open(fmake,'r') as f:
	d=f.readlines()

for xi in vs:
	#Find mode
	if xi.endswith('+'):
		vname=xi[:-1]
		mode='append'
	else:
		vname=xi
		mode='replace'

	#Find surgery line
	syms=['=',':=','+=']
	t1=[[x.lstrip().startswith(vname+y) for y in syms] for x in d]
	t1=np.array(t1)
	t2=np.nonzero(t1.any(axis=1))[0]
	if len(t2)==0:
		raise ValueError(f'Variable not found: {vname}')
	t2=t2[-1]
	t3=np.nonzero(t1[t2])[0][0]
	if t3>1:
		raise ValueError(f'Final rule is not assignment: {vname}')

	#Perform surgery
	tl=[len(d[t2]),len(d[t2].lstrip()),len(d[t2].rstrip())]
	if mode=='replace':
		d[t2]=d[t2][:tl[0]-tl[1]]+vname+syms[t3]+vs[xi]+d[t2][tl[2]:]
	elif mode=='append':
		d[t2]=d[t2][:tl[2]]+vs[xi]+d[t2][tl[2]:]
	else:
		raise ValueError(f'Unknown mode {mode}')

#Save file
with open(fmake,'w') as f:
	f.write(''.join(d))
