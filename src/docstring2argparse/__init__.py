# Lingfei Wang, 2022. All rights reserved.
"""
Converts docstring and typing to command-line interface
"""

import argparse

class function_parser_base:
	"""
	Base class for function parser
	"""
	def __call__(self,func):
		return self.parse(func)
	def parse(self,func):
		ans=self._parse(func)
		self.check(ans,func.__name__)
		return ans
	@staticmethod
	def check(ans,name):
		"""
		Format check for self._parse

		Parameters
		----------
		ans:	tuple
			Supposed return of self.parse. Direct return of self._parse.

		Returns
		-------
		tuple
			Same as input `ans` if format is good.

		Raises
		------
		AssertionError
			If incorrectly formatted.
		"""
		from inspect import _empty
		assert len(ans)==4
		assert isinstance(ans[0],str) or ans[0] is None
		assert isinstance(ans[1],str) or ans[1] is None
		assert isinstance(ans[2],list) or ans[2] is None
		assert all(x is None or (isinstance(x,tuple) and len(x)==4 and (isinstance(x[0],str) or x[0] is None) and (isinstance(x[2],str) or x[2] is None) and ((isinstance(x[3],tuple) and len(x[3])==2 and (x[3][0] is None or isinstance(x[3][0],bool) and (x[3][1] is None or x[3][0] is True))) or x[3] is None)) for x in ans[2])
		if ans[3]==_empty:
			raise TypeError(f"{name} function's output type is not annotated.")
		assert ans[3] is None or (isinstance(ans[3],tuple) and len(ans[3])==3 and (isinstance(ans[3][0],str) or ans[3][0] is None) and (isinstance(ans[3][2],str) or ans[3][2] is None)) or all(x is None or (isinstance(x,tuple) and len(x) in {3,4} and (isinstance(x[0],str) or x[0] is None) and (isinstance(x[2],str) or x[2] is None)) for x in ans[3])
		return ans
	def _parse(self,func):
		"""
		Parse function to get documentation for argparse
		
		Parameters
		----------
		func:	function
			Function to parse
		
		Returns
		-------
		doc_short:	str or NoneType
			Short documentation for subcommand introduction. Use None for not parsed.
		doc_long:	str or NoneType
			Long documentation for full description. Use None for not parsed.
		doc_parameters:	list or NoneType
			Documentation for each parameter. Each entry is a tuple (name,type,description,(whether_optional,default_value)) for each parameters. default_value must be None if whether_optional is False. doc_parameters or any of its element can be None for not parsed.
		doc_return:	list or tuple or NoneType
			Documentation for each return. Each entry is a tuple (name,type,description) for each return. doc_return can be such a tuple for single return. doc_return or any of its element can be None for not parsed. Note: for argparse purpose, this return is not checked because it is irrelevant to command I/O.
		"""
		raise NotImplementedError

class function_parser_signature(function_parser_base):
	"""
	Function parser through signature and typing
	"""
	@staticmethod
	def _parse(func):
		from inspect import signature,_empty
		#Signature
		t2=signature(func)
		a1=[None,None,[(x,y.annotation,None,(True,y.default) if hasattr(y,'default') and y.default!=_empty else (False,None)) for x,y in t2.parameters.items()],t2.return_annotation]
		return a1

class function_parser_docstring_numpy(function_parser_base):
	"""
	Function parser through docstring in numpy format
	"""
	@staticmethod
	def _parse_rst(text):
		#Credit: https://stackoverflow.com/questions/12883428/how-to-parse-restructuredtext-in-python
		import docutils.parsers.rst
		import docutils.utils
		import docutils.frontend
		parser = docutils.parsers.rst.Parser()
		components = (docutils.parsers.rst.Parser,)
		settings = docutils.frontend.OptionParser(components=components).get_default_values()
		document = docutils.utils.new_document('<rst-doc>', settings=settings)
		parser.parse(text, document)
		return document
	@classmethod
	def _parse(cls,func):
		from textwrap import dedent
		#Return buffer
		ans=[None,None,None,None]
		
		doc=dedent(func.__doc__)
		# print(doc)
		d=cls._parse_rst(doc)
		#Sections to extract
		secs={'parameters','results'}
		#Sections to ignore
		secs_hide={'yields','receives'}
		#Extract needed sections
		secs={x:list(filter(lambda y:y.tagname=='section' and y['names'][0]==x,d)) for x in secs}		# pylint: disable=W0640
		secs={x:y for x,y in secs.items() if len(y)>0}
		assert all(len(x)==1 for x in secs.values())
		#Remove sections
		t1=list(filter(lambda x:x.tagname=='section' and x['names'][0] in set(secs)|secs_hide,d))
		for xi in t1:
			d.remove(xi)
		#Short description
		t1=list(filter(lambda x:x.tagname=='paragraph',d))
		if len(t1)>0:
			ans[0]=t1[0].astext().strip()
			d.remove(t1[0])
		#Long description
		ans[1]=d.astext().strip()
		for xi in range(2):
			if len(ans[xi])==0:
				ans[xi]=None
				
		#Get parameters & results
		for xi in zip(['parameters','results'],range(2)):
			if xi[0] not in secs:
				continue
			if xi[0]=='results':
				raise NotImplementedError
			d=secs[xi[0]]
			assert len(d)==1
			d=d[0]
			assert [x.tagname for x in d]==['title','definition_list']
			d=d[1]
			assert all(x.tagname=='definition_list_item' for x in d)
			assert all([y.tagname for y in x]==['term', 'definition'] for x in d)
			d=[[y.astext().strip() for y in x] for x in d]
			assert all(len(x[0].split(':'))==2 for x in d)
			d=[[x[0].split(':')[0].strip(),x[0].split(':')[1].strip(),x[1],None] for x in d]
			d=[tuple(x) if len(x[1])>0 else tuple([x[0],None]+x[2:]) for x in d]
			ans[2+xi[1]]=d
		return ans

class function_parser_union(function_parser_base):
	"""
	Function parser that combines parsed results from multiple parsers. Complement results are combined. Conflicts are not allowed.
	"""
	def __init__(self,parsers):
		assert all(isinstance(x,function_parser_base) for x in parsers)
		self.parsers=parsers
	@classmethod
	def union(cls,val):
		v1=list(filter(lambda x: x is not None,val))
		if len(v1)==0:
			#No known values
			return None
		if len(v1)==1:
			#Single known value
			return v1[0]
		#Confirm identical types
		t1=type(v1[0])
		assert all(type(x)==t1 for x in v1[1:])		# pylint: disable=C0123
		if all(x==v1[0] for x in v1[1:]):
			#Identical values
			return v1[0]
		#Identical lengths
		assert all(hasattr(x,'__len__') for x in v1)
		t1=len(v1[0])
		if not all(len(x)==t1 for x in v1[1:]):
			raise ValueError('Union failed with different parameter counts.')
		#Recursive union
		return type(v1[0])([cls.union(x) for x in zip(*v1)])
	def _parse(self,func):
		ans=[x.parse(func) for x in self.parsers]
		return self.union(ans)

def get_functions(pkgname,parser,varname_ignore='_docstring2argparse_ignore_',func_ignore=lambda name,obj:name.startswith('_')):
	"""
	Get all functions to argparse from module
		
	Parameters
	----------
	pkgname:	str
		Name of package/module to argparse
	parser:		docstring2argparse.function_parser_base
		Parser for function documentation
	varname_ignore:	str
		Ignore list variable name in each file. The variable is a list of str. Each str is the name of function to be ignored for argparse.
	func_ignore:	function
		Function takes the name and object to determine if the object should be ignored for argparse or recursive search.
	
	Returns
	-------
	ans_f:	dict
		Dictionary that converts function name (pkgname.submodulename....functionname) to function_parser_base._parse outputs
	ans_m:	dict
		Dictionary that converts submodule name (pkgname.submodulename1....submodulenamen) to module documentation as str
	"""
	from inspect import ismodule
	from importlib import import_module
	#Find functions
	pkg=import_module(pkgname)
	ans=[]
	mods=[[pkg,[pkgname]]]
	mods_done=[]
	c=0
	while len(mods)>0:
		#Recursive search within the module
		t1=mods.pop()
		mods_done.append(t1)
		t2=t1[0].__dict__
		#Restrict to modules & functions
		ans_add=filter(lambda x:hasattr(t2[x],'__module__') and hasattr(t2[x],'__call__') and t2[x].__module__.startswith('.'.join(t1[1])),t2)
		#Withing this module
		mod_add=filter(lambda x:ismodule(t2[x]) and hasattr(t2[x],'__package__') and t2[x].__package__==pkgname,t2)
		if varname_ignore is not None and varname_ignore in t2:
			ans_add,mod_add=[filter(lambda x:x not in t2[varname_ignore],y) for y in [ans_add,mod_add]]
		if func_ignore is not None:
			ans_add,mod_add=[filter(lambda x:not func_ignore(x,t2[x]),y) for y in [ans_add,mod_add]]
		ans+=[[t2[x],t1[1]+[x]] for x in ans_add]
		mods+=[[t2[x],t1[1]+[x]] for x in mod_add]
		c+=1
		assert c<1E8
	ans={'.'.join(y):x for x,y in ans}
	ans_f={x:parser.parse(y) for x,y in ans.items()}
	ans_m={'.'.join(y):x.__doc__ for x,y in mods_done}
	ans_m={x:y.strip() if y is not None else '' for x,y in ans_m.items()}
	return (ans_f,ans_m)

def totype(t):
	"""
	Convert typing.Optional to type.
	"""
	import typing
	if isinstance(t,type):
		return t
	t1=[typing.get_origin(t),typing.get_args(t)]
	if t1[0]==typing.Union and len(t1[1])==2 and sum(x==type(None) for x in t1[1])==1:		# pylint: disable=W0143
		#typing.Optional
		return list(filter(lambda x:x!=type(None),t1[1]))[0]
	raise TypeError('Unknown type: '+repr(t))

def docstringparser(pkgname,parser,ka_argparse=dict(formatter_class=argparse.ArgumentDefaultsHelpFormatter),**ka):
	"""
	Create argparse parser for module
		
	Parameters
	----------
	pkgname:		str
		Name of package/module to argparse
	parser:			function_parser_base
		Parser for function documentation
	ka_argparse:	dict
		Keyword arguments passed to argparse.ArgumentParser
	ka:				dict
		Keyword arguments passed to docstring2argparse.get_functions
	
	Returns
	-------
	argparse.ArgumentParser
		Command line argument parser for module
	"""
	from os import linesep
	#Find functions & modules
	f,m=get_functions(pkgname,parser,**ka)
	#Parsers
	p={}
	#Subparsers
	ps={}
	p[pkgname]=argparse.ArgumentParser(prog=pkgname,description=m[pkgname],**ka_argparse)
	ps[pkgname]=p[pkgname].add_subparsers(help='sub-commands',dest='subcommand_1',required=True)

	for xi in f:
		#Iteratively add parsers for functions
		t1=xi.split('.')
		for xj in range(1,len(t1)):
			#Recursively add subparsers for all parent modules of this function if not present
			t2='.'.join(t1[:xj])
			if t2 in p:
				continue
			t3='.'.join(t2.split('.')[:-1])
			assert t3 in ps
			tka=dict(ka_argparse)
			if m[t2] is not None:
				tka['help']=m[t2]
			p[t2]=ps[t3].add_parser(t2.rsplit('.',maxsplit=1)[-1],**tka)
			ps[t2]=p[t2].add_subparsers(help='sub-commands',dest=f'subcommand_{xj}',required=True)
		#Add subparser for this function
		t2='.'.join(t1[:-1])
		tka=dict(ka_argparse)
		desc=''
		if f[xi][0] is not None:
			tka['help']=f[xi][0]
			desc+=f[xi][0]+linesep*2
		if f[xi][1] is not None:
			desc+=f[xi][1]
		p[xi]=ps[t2].add_parser(t1[-1],description=desc,**tka)
		for xj in f[xi][2]:
			if xj[3][0]:
				#Keyword arguments
				ttype=totype(xj[1])
				if ttype==bool:
					p[xi].add_argument('--'+xj[0],help=xj[2],default=False,action='store_true')
				else:
					p[xi].add_argument('--'+xj[0],type=ttype,help=xj[2],default=xj[3][1],action='store')
			else:
				#Required arguments
				p[xi].add_argument(xj[0],type=totype(xj[1]),help=xj[2],action='store')
	return (p[pkgname],f)

def run_args(pkgname,funcs,args):
	from importlib import import_module
	func=import_module(pkgname)
	fullname=[pkgname]
	lv=1
	while hasattr(args,f'subcommand_{lv}'):
		t1=f'subcommand_{lv}'
		fullname.append(getattr(args,t1))
		func=getattr(func,fullname[-1])
		lv+=1
	fullname='.'.join(fullname)
	a=[getattr(args,x[0]) for x in funcs[fullname][2] if not x[3][0]]
	ka={x[0]:getattr(args,x[0]) for x in funcs[fullname][2] if x[3][0]}
	return func(*a,**ka)

def docstringrunner(pkgname):
	import sys
	parser_func=function_parser_union([function_parser_signature(),function_parser_docstring_numpy()])
	parser_arg,funcs=docstringparser(pkgname,parser_func)
	if len(sys.argv) == 1:
		parser_arg.print_help(sys.stderr)
		sys.exit(1)
	args=parser_arg.parse_args()
	run_args(pkgname,funcs,args)



























#
