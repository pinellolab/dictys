#!/bin/bash

function usage()
{
	fmt='  %-25s%s\n'
	echo "Usage: $(basename "$0") [-h]" >&2
	echo "Installing Dictys as a conda environment" >&2
	echo "" >&2	
	echo "Parameters" >&2
	printf "$fmt" '-h' 'Display this help.' >&2
	echo "" >&2	
	echo "Accepted environmental variables" >&2
	printf "$fmt" 'CONDAENV_NAME' 'Name of conda environment' >&2
	printf "$fmt" 'PYTHONVERSION_CONDA' 'Python version' >&2
	printf "$fmt" 'CUDAVERSION_CONDA' 'CUDA version for GPU enabled pytorch. Omit or set to null str for CPU powered pytorch.' >&2
	printf "$fmt" '' 'You need a version supported by your driver and by pytorch.' >&2
	printf "$fmt" '' 'See https://stackoverflow.com/questions/9727688/how-to-get-the-cuda-version' >&2
	printf "$fmt" 'COMMIT_VERSION' 'Dictys commit to install. Accepts hash, branch name, and tag. Defaults to master.' >&2
	echo '' >&2
	echo "Advanced options in environmental variables" >&2
	printf "$fmt" 'LOCAL_VERSION' 'Use a local folder to intall Dictys from. If specified, does not download from repo and ignores COMMIT_VERSION.' >&2
	printf "$fmt" 'STEPS' 'Binary keys for which installation steps to perform. Accepts:' >&2
	printf "$fmt" '' '1: conda steps' >&2
	printf "$fmt" '' '2: pip steps' >&2
	printf "$fmt" '' '4: homer update' >&2
	printf "$fmt" '' 'Default: all steps' >&2
	printf "$fmt" 'PIP_OPTIONS' 'Options for pip install' >&2
	exit 1
}

#Parse arguments
while getopts ':h' o; do case "$o" in
	f)	field="$OPTARG";;
	:)	echo "Error: -${OPTARG} requires an argument." >&2;echo >&2;usage;;
	*)	usage;;
	esac
done
shift $((OPTIND-1))

set -ex -o pipefail

if [ "a$CONDAENV_NAME" == "a" ]; then
	CONDAENV_NAME="dictys"
fi
if [ "a$PYTHONVERSION_CONDA" == "a" ]; then
	PYTHONVERSION_CONDA="3.9"
fi
if [ "a$COMMIT_VERSION" == "a" ]; then 
	COMMIT_VERSION="master"
fi
if [ "a$CUDAVERSION_CONDA" == "a" ]; then
	conda_deps="cpuonly"
else
	conda_deps="cudatoolkit=$CUDAVERSION_CONDA"
fi

if [ "a$STEPS" == "a" ] || [ "a$(( STEPS & 1 ))" != "a0" ]; then
	#Install non-pypi dependencies: pytorch, bedtools, homer, samtools, macs2, ffmpeg
	conda create -y -n $CONDAENV_NAME -c bioconda -c conda-forge -c pytorch python=$PYTHONVERSION_CONDA pytorch torchvision torchaudio $conda_deps bedtools homer samtools macs2 ffmpeg
	#You may need "conda activate ..." instead
	. activate $CONDAENV_NAME
fi
if [ "a$STEPS" == "a" ] || [ "a$(( STEPS & 2 ))" != "a0" ]; then
	#Install Dictys
	if [ "a$LOCAL_VERSION" == "a" ]; then
		pip install $PIP_OPTIONS git+https://github.com/pinellolab/dictys.git@$COMMIT_VERSION
	else
		pip install $PIP_OPTIONS "$LOCAL_VERSION"
	fi
	#Correcting matplotlib version due to pyDNase dependency
	pip uninstall -y pyDNase
	pip install -U matplotlib
	pip install $PIP_OPTIONS --no-deps pyDNase
fi
if [ "a$STEPS" == "a" ] || [ "a$(( STEPS & 4 ))" != "a0" ]; then
	#Update homer
	pushd . &> /dev/null
	cd "$(dirname "$(dirname "$(realpath "$(which homer)")")")"
	./configureHomer.pl -update
	chmod u+x configureHomer.pl
	popd &> /dev/null
fi
if [ "a$STEPS" == "a" ] || [ "a$(( STEPS & 1 ))" != "a0" ]; then
	conda deactivate
fi
