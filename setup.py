#!/usr/bin/python3


pkgname="dictys"
pkgnamefull="Dictys"
version=[0,1,0]
license="BSD-3-Clause"
url="https://github.com/pinellolab/"+pkgname
author="Lingfei Wang, Nikolaos Trasanidis"
author_email="Lingfei.Wang.github@outlook.com, NTRASANIDIS@mgh.harvard.edu"



def pkg_setup():
	from setuptools import setup
	from os import path
	pkgnameu=pkgname[0].upper()+pkgname[1:].lower()
	with open(path.join(path.abspath(path.dirname(__file__)),'README.rst'),encoding='utf-8') as f:
		long_description=f.read()
	setup(name=pkgname,
		version='.'.join(map(str,version)),
		author=author,
		author_email=author_email,
		description=pkgnamefull,
		long_description=long_description,
		long_description_content_type='text/x-rst',
		url=url,
		# download_url=url,
		scripts=['bin/dictys'],
		# include_package_data=True,
		install_requires=['numpy','scipy','pandas','scikit-learn','biothings_client'],
		classifiers=['Development Status :: 2 - Pre-Alpha ',
			'License :: OSI Approved :: BSD License',
			'Environment :: Console',
			# 'Framework :: Pytest',
			'Intended Audience :: Science/Research',
			'Intended Audience :: Developers',
			'Operating System :: OS Independent',
			'Programming Language :: Python :: 3',
			'Topic :: Scientific/Engineering :: Bio-Informatics'],
		license=license,
		packages=[pkgname],
		package_dir={pkgname:path.join('src',pkgname)},
	)

pkg_setup()
