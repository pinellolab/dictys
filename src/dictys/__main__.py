import argparse
import docutils.nodes
import sys

def parse_rst(text: str) -> docutils.nodes.document:
	import docutils.parsers.rst
	import docutils.utils
	import docutils.frontend
	parser = docutils.parsers.rst.Parser()
	components = (docutils.parsers.rst.Parser,)
	settings = docutils.frontend.OptionParser(components=components).get_default_values()
	document = docutils.utils.new_document('<rst-doc>', settings=settings)
	parser.parse(text, document)
	return document

def get_docstring(func):
	from inspect import signature,_empty
	from os import linesep
	import numpy as np

	#Annotation type 1
	t2=signature(func)
	a1=[None,[(x,y.annotation,None,y.default) if hasattr(y,'default') and y.default!=_empty else (x,y.annotation,None) for x,y in t2.parameters.items()],t2.return_annotation]

	#Annotation type 2
	#Preprocess before parsing
	t2=linesep.join([x.strip() for x in func.__doc__.split(linesep)]).strip().split(linesep+linesep)
	assert len(t2)>=2
	t3=[parse_rst(x) for x in t2]
	#Postprocess after parsing
	assert all([len(x)==1 for x in t3])
	t3=[x[0] for x in t3]
	assert all([x.tagname in {'section','paragraph'} for x in t3])
	#find sections for Parameters and Returns
	t4=[[len(x)==1 for x in t3],[len(x)==2 and x[0].astext()=='Parameters' for x in t3],[len(x)==2 and x[0].astext()=='Returns' for x in t3]]
	assert all([x^y^z for x,y,z in zip(*t4)])
	t4=[np.nonzero(x)[0] for x in t4[1:]]
	assert all([len(x)<=1 for x in t4])
	t4=[t4[0][0] if len(t4[0])==1 else 0,t4[1] if len(t4[1])==1 else len(t3)]
	a2=[(linesep*2).join([x.astext() for x in t3[:t4[0]]])]
	if t4[0]!=0:
		#Parameters
		t5=linesep.join([t3[t4[0]][1].astext()]+[x.astext() for x in t3[t4[0]+1:t4[1]]]).split(linesep)
		t5=[x.strip() for x in t5]
		t6=[x.split(':') for x in t5]
		assert len(t6)%2==0
		t7=[len(x)==2 and (len(x[1])==0 or x[1].lstrip()!=x[1]) for x in t6]
		assert (np.array(t7)==~(np.arange(len(t7))%2).astype(bool)).all()
		t6=[(t5[2*x].split(':')[0],t5[2*x].split(':')[1].strip(),t5[2*x+1]) for x in range(len(t5)//2)]
		a2.append(t6)
	else:
		a2.append([])
	if t4[1]!=len(t3):
		#Returns
		t5=linesep.join([t3[t4[1]][1].astext()]+[x.astext() for x in t3[t4[1]+1:]]).split(linesep)
		t5=[x.strip() for x in t5]
		t6=[x.split(':') for x in t5]
		assert len(t6)%2==0
		t7=[len(x)==2 and (len(x[1])==0 or x[1].lstrip()!=x[1]) for x in t6]
		assert (np.array(t7)==~(np.arange(len(t7))%2).astype(bool)).all()
		t6=[(t5[2*x].split(':')[0],t5[2*x].split(':')[1].strip(),t5[2*x+1]) for x in range(len(t5)//2)]
		a2.append(t6)
	else:
		a2.append(None)
	#Confirm identical annotations
	assert len(a1[1])==len(a2[1])
	assert all([x[0]==y[0] and x[1].__name__==y[1] for x,y in zip(a1[1],a2[1])])
	assert a1[2]==a2[2]
	#Merge annotations
	a1=[a2[0],[(x[0],x[1],y[2]) if len(x)==3 else (x[0],x[1],y[2],x[3]) for x,y in zip(a1[1],a2[1])],None if a1[2]==None else [(x[0],x[1],y[2]) for x,y in zip(a1[2],a2[2])]]
	return a1

def get_functions(pkgname):
	from inspect import ismodule
	from importlib import import_module
	
	#Find functions
	pkg=import_module(pkgname)
	ans=[]
	mods=[[pkg,[pkgname]]]
	mods_done=[]
	c=0
	while len(mods)>0:
		t1=mods.pop()
		mods_done.append(t1)
		t2=t1[0].__dict__
		ans+=[[t2[x],t1[1]+[x]] for x in t2 if not x.startswith('_') and hasattr(t2[x],'__module__') and t2[x].__module__.startswith('.'.join(t1[1]))]
		mods+=[[t2[x],t1[1]+[x]] for x in t2 if not x.startswith('_') and ismodule(t2[x]) and hasattr(t2[x],'__package__') and t2[x].__package__==pkgname]
		c+=1
		assert c<1E8
	ans={'.'.join(y):x for x,y in ans}
	ans_f={x:get_docstring(y) for x,y in ans.items()}
	ans_m={'.'.join(y):x.__doc__ for x,y in mods_done}
	ans_m={x:y.strip() if y is not None else '' for x,y in ans_m.items()}
	return [ans_f,ans_m]


ka_argparse=dict(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
pkgname=__package__
f,m=get_functions(pkgname)
p=dict()
ps=dict()
p[pkgname]=argparse.ArgumentParser(prog=pkgname,description=m[pkgname],**ka_argparse)
ps[pkgname]=p[pkgname].add_subparsers(help='sub-commands',required=True)

for xi in f:
	t1=xi.split('.')
	for xj in range(1,len(t1)):
		t2='.'.join(t1[:xj])
		if t2 not in p:
			t3='.'.join(t2.split('.')[:-1])
			assert t3 in ps
			p[t2]=ps[t3].add_parser(t2.split('.')[-1],help=m[t2],**ka_argparse)
			ps[t2]=p[t2].add_subparsers(help='sub-commands',required=True)
	t2='.'.join(t1[:-1])
	p[xi]=ps[t2].add_parser(t1[-1],description=f[xi][0],**ka_argparse)
	for xj in f[xi][1]:
		if len(xj)==3:
			p[xi].add_argument(xj[0],type=xj[1],help=xj[2],action='store')
		elif len(xj)==4:
			p[xi].add_argument('--'+xj[0],type=xj[1],help=xj[2],default=xj[3],action='store')
		else:
			assert False

# if len(sys.argv) == 1:
# 	p[pkgname].print_help(sys.stderr)
# 	sys.exit(1)

args=p[pkgname].parse_args()
# args=vars(args)
# logging.basicConfig(
# 	format=
# 	'%(levelname)s:%(process)d:%(asctime)s:%(pathname)s:%(lineno)d:%(message)s',
# 	level=logging.DEBUG if args['verbose'] else logging.WARNING)
# from . import run
# func = getattr(run, args['cmd'])
# func(args)









#
