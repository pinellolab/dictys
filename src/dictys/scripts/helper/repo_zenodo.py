#!/usr/bin/env python3
# Lingfei Wang, 2022. All rights reserved.
import argparse

class zenodod:
	#Zenodo dataset
	def __init__(self,doi,url='https://zenodo.org',token=None):
		self.url=url
		self.doi=doi
		pid=doi.split('.')[-1]
		self.pid=pid
		self.token=token

	def delete_byid(self,pid):
		raise NotImplementedError
		import requests
		from requests.auth import HTTPBasicAuth

		basic = HTTPBasicAuth(self.token, '')
		url=self.url+f'/dvn/api/data-deposit/v1.1/swordv2/edit-media/file/{pid}'
		resp=requests.delete(url, auth=basic)
		resp.raise_for_status()

	def test(self):
		#Expensive test
		import requests
		resp=requests.get(self.url+'/api/records')
		resp.raise_for_status()

	def ls(self,mode='all'):
		#List files
		import requests
		params={}
		if mode=='record':
			url=self.url+f'/api/records/{self.pid}'
			if self.token is not None:
				params['access_token']=self.token
			resp=requests.get(url,params=params)
			resp.raise_for_status()
			fs={x['key']:x['links']['self'] for x in resp.json()['files']}
		elif mode=='deposition':
			url=self.url+f'/api/deposit/depositions/{self.pid}/files'
			if self.token is not None:
				params['access_token']=self.token
			resp=requests.get(url,params=params)
			resp.raise_for_status()
			fs={x['filename']:x['links']['download'] for x in resp.json()}
		elif mode=='all':
			fs=None
			errs={}
			for mode1 in ['record','deposition']:
				try:
					fs=self.ls(mode=mode1)
					break
				except requests.exceptions.HTTPError as e:
					errs[mode1]=e
			if fs is None:
				from os import linesep
				raise RuntimeError('All modes failed:'+linesep+repr(errs))
		else:
			raise ValueError(f'Unknown mode {mode}')
		return fs

	def upload_datafile(self,src,dst):
		import requests
		params={}
		url=self.url+f'/api/deposit/depositions/{self.pid}'
		if self.token is not None:
			params['access_token']=self.token
		resp=requests.get(url,params=params)
		resp.raise_for_status()
		bucket=resp.json()['links']['bucket']
		
		with open(src, "rb") as fp:
		    resp=requests.put(f'{bucket}/{dst}',data=fp,params=params)
		resp.raise_for_status()
		return resp

	def upload(self,src,dst=None,exist='overwrite'):
		from os.path import basename
		if dst is None:
			dst=basename(src)
		else:
			raise NotImplementedError

		fs=self.ls()
		if dst in fs:
			if exist=='overwrite':
				pass
			elif exist=='skip':
				return
			else:
				raise ValueError(f'Unknown exist={exist}')
		#Upload file
		self.upload_datafile(src,dst)
		
	@staticmethod
	def download_byurl(url,dst):
		import requests
		import shutil
		with requests.get(url, stream=True) as r:
			with open(dst, 'wb') as f:
				shutil.copyfileobj(r.raw, f)

	def download(self,src,dst=None):
		if dst is None:
			dst=src

		#Find file ID
		fs=self.ls()
		if src not in fs:
			raise FileNotFoundError(src)
		#Download file
		self.download_byurl(fs[src],dst)

def down(args):
	d=zenodod(args.doi,token=args.token,url=args.base)
	return d.download(args.src,dst=args.dst)

def up(args):
	d=zenodod(args.doi,token=args.token,url=args.base)
	return d.upload(args.src,exist=args.existing)

def ls(args):
	from os import linesep
	d=zenodod(args.doi,token=args.token,url=args.base)
	print(linesep.join(sorted(d.ls())))

parser=argparse.ArgumentParser(description="File operation with Zenodo repository.")
subparsers=parser.add_subparsers(help='sub-commands',dest='subcommand_1',required=True)
subparser=subparsers.add_parser('up',help='Upload file to Zenodo.',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparser.add_argument('src',type=str,help='File path to upload')
subparser.add_argument('--doi',default='10.5281/zenodo.6787658',type=str,help='Dataset doi')
subparser.add_argument('--token',default=None,type=str,help='Access token')
subparser.add_argument('--base',default='https://zenodo.org',type=str,help='Base URL for Zenodo repository')
subparser.add_argument('--existing',default='overwrite',type=str,help='How to handle existing files in repository. Accepts: overwrite, skip.')
subparser.set_defaults(func=up)

subparser=subparsers.add_parser('down',help='Download file from Zenodo.',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparser.add_argument('src',type=str,help='File name to download')
subparser.add_argument('--dst',default=None,type=str,help='Path to save downloaded file. Defaults to same name in current folder.')
subparser.add_argument('--doi',default='10.5281/zenodo.6787658',type=str,help='Dataset doi')
subparser.add_argument('--token',default=None,type=str,help='Access token')
subparser.add_argument('--base',default='https://zenodo.org',type=str,help='Base URL for Zenodo repository')
subparser.set_defaults(func=down)

subparser=subparsers.add_parser('ls',help='List files from Zenodo.',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparser.add_argument('--doi',default='10.5281/zenodo.6787658',type=str,help='Dataset doi')
subparser.add_argument('--token',default=None,type=str,help='Access token')
subparser.add_argument('--base',default='https://zenodo.org',type=str,help='Base URL for Zenodo repository')
subparser.set_defaults(func=ls)

args=parser.parse_args()
args.func(args)
