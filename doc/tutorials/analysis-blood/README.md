# Network analysis for human blood (BMMC)

This is a [Dictys](https://github.com/pinellolab/dictys) tutorial for context specific and dynamic gene regulatory network (GRN) analysis for separate scRNA-seq and scATAC-seq datasets of human blood (BMMC).

These experimental datasets were previously published and annotated in [Granja et al 2019](https://dx.doi.org/10.1038/s41587-019-0332-7) and [Granja et al 2021](https://dx.doi.org/10.1038/s41588-021-00790-6). We aligned them onto the same low dimensional embedding in [our paper](https://doi.org/10.1101/2022.09.14.508036). We inferred the developmental trajectories in the tutorial [trajectory-blood](https://github.com/pinellolab/dictys/tree/master/doc/tutorials/trajectory-blood) (TBA) and performed context specific and dynamic network inference in the tutorial [inference-blood](https://github.com/pinellolab/dictys/tree/master/doc/tutorials/inference-blood) (TBA). The reconstructed networks are provided here as input for this tutorial.

To use this tutorial, [install Dictys](https://github.com/pinellolab/dictys#installation), download this folder, and serve the notebook with jupyter inside Dictys' environment. You can start with `notebooks/static/main.ipynb` or `notebooks/dynamic/main.ipynb` for context specific (i.e. static) or dynamic GRN analyses respective. This also downloads the pre-inferred networks. Then you can try other notebooks with extra details.

If you face any issues or need any assistance, see [FAQ](https://github.com/pinellolab/dictys#faq) and [Issues](https://github.com/pinellolab/dictys#issues).
