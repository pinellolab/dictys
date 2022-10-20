#!/usr/bin/env python3
# Lingfei Wang, 2022. All rights reserved.

import argparse

class dvd:
	#Dataverse dataset
	def __init__(self,url,doi,token=None):
		from pyDataverse.api import NativeApi,DataAccessApi
		self.url=url
		self.doi=doi
		self.pid='doi:'+doi
		self.token=token
		self.api=NativeApi(self.url,self.token)
		self.dapi=DataAccessApi(self.url,self.token)
		self.test()
	def delete_byid(self,pid):
		import requests
		from requests.auth import HTTPBasicAuth

		basic = HTTPBasicAuth(self.token, '')
		url=self.url+f'/dvn/api/data-deposit/v1.1/swordv2/edit-media/file/{pid}'
		resp=requests.delete(url, auth=basic)
		resp.raise_for_status()

	def test(self):
		self.api.get_info_version().raise_for_status()

	def ls(self):
		#List files
		resp=self.api.get_dataset(self.pid)
		resp.raise_for_status()
		fs=resp.json()['data']['latestVersion']['files']
		assert all(x['label']==x['dataFile']['filename'] for x in fs)
		assert len(set(x['label'] for x in fs))==len(fs)
		fs={x['label']:x['dataFile']['id'] for x in fs}
		return fs

	def upload(self,src,dst=None,exist='overwrite'):
		from os.path import basename
		if dst is None:
			dst=basename(src)
		else:
			raise NotImplementedError

		fs=self.ls()
		if dst in fs:
			if exist=='overwrite':
				self.delete_byid(fs[dst])
			elif exist=='skip':
				return
			else:
				raise ValueError(f'Unknown exist={exist}')

		#Upload file
		resp=self.api.upload_datafile(self.pid,src)
		resp.raise_for_status()

	def download_byid(self,id,dst):
		#Download file
		resp=self.dapi.get_datafile(id)
		resp.raise_for_status()
		with open(dst,'wb') as f:
			f.write(resp.content)

	def download(self,src,dst=None):
		if dst is None:
			dst=src

		#Find file ID
		fs=self.ls()
		if src not in fs:
			raise FileNotFoundError(src)
		#Download file
		self.download_byid(fs[src],dst)

def down(args):
	d=dvd(args.base,args.doi,token=args.token)
	return d.download(args.src,dst=args.dst)

def up(args):
	d=dvd(args.base,args.doi,token=args.token)
	return d.upload(args.src,exist=args.existing)

def ls(args):
	from os import linesep
	d=dvd(args.base,args.doi,token=args.token)
	print(linesep.join(sorted(d.ls())))

parser=argparse.ArgumentParser(description="File operation with Dataverse repository.")
subparsers=parser.add_subparsers(help='sub-commands',dest='subcommand_1',required=True)
subparser=subparsers.add_parser('up',help='Upload file to Dataverse.',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparser.add_argument('src',type=str,help='File path to upload')
subparser.add_argument('--doi',default='10.7910/DVN/IJWFWR',type=str,help='Dataset DOI')
subparser.add_argument('--token',default=None,type=str,help='Access token')
subparser.add_argument('--base',default='https://dataverse.harvard.edu',type=str,help='Base URL for dataverse repository')
subparser.add_argument('--existing',default='overwrite',type=str,help='How to handle existing files in repository. Accepts: overwrite, skip.')
subparser.set_defaults(func=up)

subparser=subparsers.add_parser('down',help='Download file from Dataverse.',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparser.add_argument('src',type=str,help='File name to download')
subparser.add_argument('--dst',default=None,type=str,help='Path to save downloaded file. Defaults to same name in current folder.')
subparser.add_argument('--doi',default='10.7910/DVN/IJWFWR',type=str,help='Dataset DOI')
subparser.add_argument('--token',default=None,type=str,help='Access token')
subparser.add_argument('--base',default='https://dataverse.harvard.edu',type=str,help='Base URL for dataverse repository')
subparser.set_defaults(func=down)

subparser=subparsers.add_parser('ls',help='List files from Dataverse.',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparser.add_argument('--doi',default='10.7910/DVN/IJWFWR',type=str,help='Dataset DOI')
subparser.add_argument('--token',default=None,type=str,help='Access token')
subparser.add_argument('--base',default='https://dataverse.harvard.edu',type=str,help='Base URL for dataverse repository')
subparser.set_defaults(func=ls)

args=parser.parse_args()
args.func(args)
