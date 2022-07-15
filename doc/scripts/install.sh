#Install non-pypi dependencies: pytorch, bedtools, homer, samtools, macs2, ffmpeg
conda create -y -n dictys_env_name -c bioconda -c conda-forge -c pytorch pytorch torchvision torchaudio cpuonly bedtools homer samtools macs2 ffmpeg
#You may need "conda activate ..." instead
. activate dictys_env_name
#Install pypi dependencies
pip install numpy pandas docutils h5py pyro-ppl==1.6.0 scipy networkx pybedtools pyDNase threadpoolctl joblib matplotlib jupyter
#Correcting matplotlib version due to pyDNase dependency
pip uninstall -y pyDNase
pip install -U matplotlib
pip install --no-deps pyDNase 
#Install Dictys
pip install git+https://github.com/pinellolab/dictys.git
conda deactivate
