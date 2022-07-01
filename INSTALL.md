# Installation

Dictys depends on several softwares. There are often multiple methods to install each software if one method fails.
* For `bedtools`, `homer`, `samtools`, and `macs2`, you can install them through `conda`, `module add`, or manually.
* For `pytorch`, you can use any [supported method](https://pytorch.org/get-started/locally/) for your computing platform.
* For other dependencies, you can install with `pip install numpy pandas docutils h5py pyro-ppl==1.6.0 scipy networkx pybedtools pyDNase threadpoolctl joblib matplotlib jupyter`.
* For `dictys`, you can install with `pip install git+https://github.com/pinellolab/dictys.git`

You can also use precompile docker or singularity containers (TBA).

Please raise an issue on [github](https://github.com/pinellolab/dictys/issues/new) if none of the above works.
