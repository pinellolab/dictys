# Installation

## GPU support
Dictys relies on `pytorch` for GPU capacity. To enable GPU support, simply install `pytorch` with GPU support. Specifically, this installs Dictys with GPU capacity from Anaconda:
```
#Install non-pypi dependencies: pytorch, bedtools, homer, samtools, macs2, ffmpeg
#Replace cudatoolkit version with what's supported by your driver. See https://stackoverflow.com/questions/9727688/how-to-get-the-cuda-version.
conda create -y -n dictys_env_name -c bioconda -c conda-forge -c pytorch python=3.9 pytorch torchvision torchaudio cudatoolkit=11.3 bedtools homer samtools macs2 ffmpeg
#You may need "conda activate ..." instead
. activate dictys_env_name
#Install pypi dependencies
pip install numpy pandas docutils h5py pyro-ppl==1.6.0 scipy networkx pybedtools pyDNase threadpoolctl joblib matplotlib jupyter adjustText
#Install Dictys
pip install git+https://github.com/pinellolab/dictys.git
#Correcting matplotlib version due to pyDNase dependency
pip uninstall -y pyDNase
pip install -U matplotlib
pip install --no-deps pyDNase
conda deactivate
```

## Custom installation

Dictys depends on several softwares. There are often multiple methods to install each software if one method fails.
* For `bedtools`, `homer`, `samtools`, `macs2`, and `ffmpeg`, you can install them through `conda`, `module add`, or manually.
* For `pytorch`, you can use any [supported method](https://pytorch.org/get-started/locally/) for your computing platform.
* For other dependencies, you can install with `pip install numpy pandas docutils h5py pyro-ppl==1.6.0 scipy networkx pybedtools pyDNase threadpoolctl joblib matplotlib jupyter adjustText`. This is based on [pyro](https://pyro.ai) 1.6.0.
* For `dictys`, you can install with `pip install git+https://github.com/pinellolab/dictys.git`
* For upgrading `matplotlib`, you need to uninstall `pyDNase` first, upgrade `matplotlib`, and then reinstall `pyDNase` ignoring its dependencies.

You can also use precompiled docker or singularity containers (TBA).

## Issues

Please raise an issue on [github](https://github.com/pinellolab/dictys/issues/new) if none of the above works.
