# Installation

## Custom installation
Dictys depends on several softwares. There are often multiple methods to install each software if one method fails.
* For `bedtools`, `homer`, `samtools`, `macs2`, and `ffmpeg`, you can install them through `conda`, `module add`, or manually.
* For `pytorch`, you can use any [supported method](https://pytorch.org/get-started/locally/) for your computing platform.
* For other dependencies, you can install with `pip install numpy pandas docutils h5py pyro-ppl==1.6.0 scipy networkx pybedtools pyDNase threadpoolctl joblib matplotlib jupyter adjustText`. This is based on [pyro](https://pyro.ai) 1.6.0.
* For `dictys`, you can install with `pip install git+https://github.com/pinellolab/dictys.git`
* For upgrading `matplotlib`, you need to uninstall `pyDNase` first, upgrade `matplotlib`, and then reinstall `pyDNase` ignoring its dependencies.

You can also use precompiled docker or singularity containers (TBA).

Dictys relies on `pytorch` for GPU capacity. To enable GPU support, simply install `pytorch` with GPU support.

## Issues

Please raise an issue on [github](https://github.com/pinellolab/dictys/issues/new) if none of the above works.
