# Installation

## GPU support
Dictys relies on `pytorch` for GPU capacity. To enable GPU support, simply install `pytorch` with GPU support. Specifically, this installs Dictys with GPU capacity from Anaconda:

.. code-block::

	#Install non-pypi dependencies: pytorch, bedtools, homer, samtools, macs2
	#Replace cudatoolkit version with what's supported by your driver. See https://stackoverflow.com/questions/9727688/how-to-get-the-cuda-version.
	conda create -y -n dictys_env_name -c bioconda -c conda-forge -c pytorch pytorch torchvision torchaudio cudatoolkit=11.1 bedtools homer samtools macs2
	#You may need "conda activate ..." instead
	. activate dictys_env_name
	#Install pypi dependencies
	pip install numpy pandas docutils h5py pyro-ppl==1.6.0 scipy networkx pybedtools pyDNase threadpoolctl joblib matplotlib jupyter
	#Install Dictys
	pip install git+https://github.com/pinellolab/dictys.git
	conda deactivate

## Custom installation

Dictys depends on several softwares. There are often multiple methods to install each software if one method fails.
* For `bedtools`, `homer`, `samtools`, and `macs2`, you can install them through `conda`, `module add`, or manually.
* For `pytorch`, you can use any [supported method](https://pytorch.org/get-started/locally/) for your computing platform.
* For other dependencies, you can install with `pip install numpy pandas docutils h5py pyro-ppl==1.6.0 scipy networkx pybedtools pyDNase threadpoolctl joblib matplotlib jupyter`. This is based on [pyro](https://pyro.ai) 1.6.0.
* For `dictys`, you can install with `pip install git+https://github.com/pinellolab/dictys.git`

You can also use precompiled docker or singularity containers (TBA).

## Issues

Please raise an issue on [github](https://github.com/pinellolab/dictys/issues/new) if none of the above works.

