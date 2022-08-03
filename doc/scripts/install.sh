# Install Dictys
# Accepted environmental variables:
#	CONDAENV_NAME:			Name of conda environment
#	PYTHONVERSION_CONDA:	Python version
#	CUDAVERSION_CONDA:		CUDA version for GPU enabled pytorch. Omit or set to null str for CPU powered pytorch.
#							You need a version supported by your driver and by pytorch.
#							See https://stackoverflow.com/questions/9727688/how-to-get-the-cuda-version
#	COMMIT_VERSION:			Dictys commit to install. Accepts hash, branch name, and tag. Defaults to master.

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
conda deactivate
