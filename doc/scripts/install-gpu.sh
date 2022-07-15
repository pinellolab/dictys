set -ex -o pipefail
if [ "a$PYTHONVERSION_CONDA" == "a" ]; then exit 1; fi
if [ "a$CUDAVERSION_CONDA" == "a" ]; then exit 1; fi

#Install non-pypi dependencies: pytorch, bedtools, homer, samtools, macs2, ffmpeg
#Replace cudatoolkit version with what's supported by your driver. See https://stackoverflow.com/questions/9727688/how-to-get-the-cuda-version.
conda create -y -n dictys_env_name -c bioconda -c conda-forge -c pytorch python=$PYTHONVERSION_CONDA pytorch torchvision torchaudio cudatoolkit=$CUDAVERSION_CONDA bedtools homer samtools macs2 ffmpeg
#You may need "conda activate ..." instead
. activate dictys_env_name
#Install pypi dependencies
pip install numpy pandas docutils h5py pyro-ppl==1.6.0 scipy networkx pybedtools pyDNase threadpoolctl joblib matplotlib jupyter
#Install Dictys
pip install git+https://github.com/pinellolab/dictys.git
#Correcting matplotlib version due to pyDNase dependency
pip uninstall -y pyDNase
pip install -U matplotlib
pip install --no-deps pyDNase
conda deactivate
