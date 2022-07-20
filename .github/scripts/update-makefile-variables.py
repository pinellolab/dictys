# Usage: python3 $0 makefile_path json_string
# Updates makefile variable assignments with values provided in json file
# Only edits final assignment. Supports = and := assignments. 
# Json string format {"variable_name":"new_value"}. Variable names can have "+" suffix to indicate addition to current value.
# New makefile will overwrite the input makefile

import json
import sys
import numpy as np

assert len(sys.argv)==3
fmake=sys.argv[1]
vs=sys.argv[2]

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
