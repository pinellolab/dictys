#!/bin/bash

function usage()
{
	fmt='%-25s%s\n'
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

#Install non-pypi dependencies: pytorch, bedtools, homer, samtools, macs2, ffmpeg
conda create -y -n $CONDAENV_NAME -c bioconda -c conda-forge -c pytorch python=$PYTHONVERSION_CONDA pytorch torchvision torchaudio $conda_deps bedtools homer samtools macs2 ffmpeg
#You may need "conda activate ..." instead
. activate $CONDAENV_NAME
#Install pypi dependencies
pip install numpy pandas docutils h5py pyro-ppl==1.6.0 scipy networkx pybedtools pyDNase threadpoolctl joblib matplotlib jupyter adjustText
#Install Dictys
pip install git+https://github.com/pinellolab/dictys.git@$COMMIT_VERSION
#Correcting matplotlib version due to pyDNase dependency
pip uninstall -y pyDNase
pip install -U matplotlib
pip install --no-deps pyDNase
#Update homer
cd "$(dirname "$(dirname "$(realpath "$(which homer)")")")"
./configureHomer.pl -update
chmod u+x configureHomer.pl
conda deactivate
