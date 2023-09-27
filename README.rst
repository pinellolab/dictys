=========
Dictys
=========
Dictys reconstructs cell-type specific and dynamic gene regulatory networks (GRN) from scRNA-seq and scATAC-seq datasets. Dictys first infers a Transcription Factor (TF) binding network with TF footprinting from single-cell chromatin accessibility. Then Dictys refines the edges with single-cell transcriptome. Dictys addresses traditional challenges in network inference by orienting causality with TF binding information, modeling transcriptional rate to reconstruct cycle-compatible networks, and using probabilistic programming to capture the scRNA-seq process.

Dictys provides network analysis and visualization at global (across all cell types), pairwise (between two cell types) and single GRN levels. Dictys directly quantifies TF regulatory activity from GRN and enables a series of analyses such as cell-type specific TF discovery as regulation markers, differential regulation analysis alongside differential expression, and TF regulatory program illustration through its subnetwork and top activation/repression targets. These GRN-based analyses can capture unique biological insights not available from mean expression.

Dictys infers and analyzes dynamic GRN from scRNA-seq and scATAC-seq datasets along (inferred) trajectories. This avoids artificial cell subsets and potential biases from population imbalance, and allows (pseudo-)time-resolved discovery and investigation of driver TFs and their individual regulations. Dictys provides an integrative network viewer for dynamic GRN visualization of synchronous panels in animation.

Overview
=============

.. image:: https://raw.githubusercontent.com/pinellolab/dictys/master/doc/images/logo.png
   :width: 180

.. image:: https://raw.githubusercontent.com/pinellolab/dictys/master/doc/images/Dictys_overview.png
   :width: 1000


Installation
=============
Dictys has dependencies not in python. The options below automatically install these dependencies. Installation should take ~<10 mins.

Option 1: with Anaconda
-----------------------
First install `Anaconda/Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_. Then, install Dictys and PyTorch **with CPU computation**:

.. code-block::

	conda create -y -n dictys -c conda-forge python=3.9 mamba
	. activate dictys
	mamba install -y -c lingfeiwang -c bioconda -c conda-forge -c pytorch dictys pytorch torchvision torchaudio cpuonly

This will create a conda environment named ``dictys``.

Alternatively, **with GPU computation** for PyTorch (here CUDA 11.7):

.. code-block::

	conda create -y -n dictys -c conda-forge python=3.9 mamba
	. activate dictys
	mamba install -y -c lingfeiwang -c bioconda -c conda-forge -c pytorch -c nvidia dictys pytorch torchvision torchaudio pytorch-cuda=11.7

Or, with earlier versions (here CUDA 11.3, only supported in PyTorch 1):

.. code-block::

	conda create -y -n dictys -c conda-forge python=3.9 mamba
	. activate dictys
	mamba install -y -c lingfeiwang -c bioconda -c conda-forge -c pytorch -c nvidia dictys pytorch==1.12.1 torchvision==0.13.1 torchaudio==0.12.1 cudatoolkit=11.3

Option 2: with `bash script <https://tinyurl.com/dictys>`_
----------------------------------------------------------
First install `Anaconda/Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_. Then, install Dictys and PyTorch **with CPU computation**:

.. code-block::

	wget https://tinyurl.com/dictys -O - | bash

This will create a conda environment named `dictys`.

Alternatively, under a different conda environment name:

.. code-block::

	wget https://tinyurl.com/dictys -O - | CONDAENV_NAME=your_favorite_name bash

Alternatively, **with GPU computation** for PyTorch (here CUDA 11.7):

.. code-block::

	wget https://tinyurl.com/dictys -O - | CUDAVERSION_CONDA=11.7 bash

Option 3: with containers
-------------------------
To pull and run the pre-built docker image for Dictys **with CPU computation**:

.. code-block::

	docker pull lfwa/dictys-cpu
	#Add public ports with '--expose' or '-p' to serve jupyter notebooks and bind mount with '-v' to transfer input/output data
	docker run -it lfwa/dictys-cpu

Inside the container, activate conda environment and serve jupyter notebooks:

.. code-block::

	. activate dictys
	jupyter notebook --allow-root

Then, you can access jupyter notebooks with the exposed or published ports.

Additional notes
----------------
For more advanced installation, see `INSTALL.md <https://github.com/pinellolab/dictys/blob/master/INSTALL.md>`_ and/or edit the `install script <https://tinyurl.com/dictys>`_.

*Note: dynamic network inference is computationally intensive and GPU availability is highly recommended.* Running time depends on the dataset, but it can take weeks or longer without a GPU.

If you need `STREAM <https://github.com/pinellolab/STREAM>`_, `ArchR <https://www.archrproject.com/>`_, or other softwares upstream of Dictys, we recommend to install them in separate environments following their official instructions.

Updating Dictys
----------------
If your minor version **is the latest** (e.g. your installed version is **1.0**.0 and the `latest release <https://github.com/pinellolab/dictys/releases>`_ is **1.0**.9), you can update Dictys to the latest github version with ``pip3 install --no-deps --force-reinstall git+https://github.com/pinellolab/dictys`` inside your Dictys conda environment.

If your minor version **is not the latest** (e.g. your installed version is **1.0**.0 but the `latest release <https://github.com/pinellolab/dictys/releases>`_ is **1.1**.0), you should reinstall Dictys in a new conda environment with any option above.

Tutorials
=========
We provide several tutorials for different data types. Please download each tutorial folder structure before running. Note that these tutorials are not intended to fully replicate the results in the paper due to differences in software versions, computing platforms, various randomness e.g. in `HOMER genome preparsing <http://homer.ucsd.edu/homer/ngs/peakMotifs.html>`_ or `Pytorch algorithms <https://pytorch.org/docs/stable/notes/randomness.html>`_, etc.

1. `short-multiome <https://www.github.com/pinellolab/dictys/blob/master/doc/tutorials/short-multiome>`_: a single-notebook tutorial for the data preparation, inference, and analysis of context specific networks on 10x multiome data for human blood.

2. `full-multiome <https://www.github.com/pinellolab/dictys/blob/master/doc/tutorials/full-multiome>`_: an extended version of the above tutorial with detailed usage.

3. `full-skin <https://www.github.com/pinellolab/dictys/blob/master/doc/tutorials/full-skin>`_: a short tutorial for the inference and analysis of dynamic networks on SHARE-seq data for mouse skin. Contains a simple demonstration to account for covariates.

The network analysis tutorials below use the same reconstructed networks as in the paper and are designed to fully replicate the results.

1. `analysis-blood <https://www.github.com/pinellolab/dictys/blob/master/doc/tutorials/analysis-blood>`_: a simple tutorial for context specific and dynamic network analysis on separate scRNA-seq and scATAC-seq quantifications of human blood as in manuscript.

2. `analysis-skin <https://www.github.com/pinellolab/dictys/blob/master/doc/tutorials/analysis-skin>`_: a simple tutorial for context specific network analysis on SHARE-seq of mouse skin as in manuscript.

Gallery
=======
The figures below are produced with the blood example dataset. You can reproduce them with the `analysis-blood` example. See `Tutorials`_. Each figure is linked to the jupyter notebook that produces it.

Cell-type specific GRN analyses
-------------------------------
`Regulation marker TF discovery <https://nbviewer.org/github/pinellolab/dictys/blob/master/doc/tutorials/analysis-blood/notebooks/static/main.ipynb#Regulation-marker-TF-discovery-with-dot-plot-(global-dotplot.ipynb)>`_

.. image:: https://raw.githubusercontent.com/pinellolab/dictys/master/doc/images/Global_dotplot.png
   :width: 300

`Top activation target heatmap for select TFs <https://nbviewer.org/github/pinellolab/dictys/blob/master/doc/tutorials/analysis-blood/notebooks/static/main.ipynb#Heatmap-of-regulation-strengths-between-select-TFs-and-their-top-targets-in-select-cell-types-(global-heatmap.ipynb)>`_

.. image:: https://raw.githubusercontent.com/pinellolab/dictys/master/doc/images/Global_heatmap.png
   :width: 400

`Differential regulation v.s. differential expression scatter plot; integrative TF rank plot <https://nbviewer.org/github/pinellolab/dictys/blob/master/doc/tutorials/analysis-blood/notebooks/static/main.ipynb#Scatter-plot-and-bar-plot-of-differential-regulation-&-differential-expression-between-two-cell-clusters--(pair-diff.ipynb)>`_

.. image:: https://raw.githubusercontent.com/pinellolab/dictys/master/doc/images/Diff_analysis.png
   :width: 750

`Subnetwork for select TF <https://nbviewer.org/github/pinellolab/dictys/blob/master/doc/tutorials/analysis-blood/notebooks/static/main.ipynb#Draw-target-gene-subnetwork-of-a-TF-(subnet.ipynb)>`_

.. image:: https://raw.githubusercontent.com/pinellolab/dictys/master/doc/images/Subnet.png
   :width: 300
   
Dynamic GRN analysis
--------------------
`Driver TF discovery based on regulatory activity curve <https://nbviewer.org/github/pinellolab/dictys/blob/master/doc/tutorials/analysis-blood/notebooks/dynamic/main.ipynb#TF-discovery-based-on-4-patterns-of-highly-variable-regulatory-activity-over-developmental-trajectory-(discovery.ipynb)>`_

.. image:: https://raw.githubusercontent.com/pinellolab/dictys/master/doc/images/Dynamic_discovery.png
   :width: 1050

`Dynamic GRN animation <https://nbviewer.org/github/pinellolab/dictys/blob/master/doc/tutorials/analysis-blood/notebooks/dynamic/main.ipynb#Animation-visualization-of-dynamic-networks>`_

.. image:: https://raw.githubusercontent.com/pinellolab/dictys/master/doc/images/animation.gif
   :width: 800

FAQ
==========================
* **Can I use fragments files instead of bam files for chromatin accessibility**?

  Dictys uses wellington (pyDNase) which does not accept fragments files. Most journals require raw data for publication so *in theory* you should be able to obtain the bam files or something equivalent. If you really cannot obtain bam files, there are ways that may circumvent this requirement: (i) fork and patch `pyDNase <https://github.com/jpiper/pyDNase>`_ to accept fragments files or (ii) convert fragments files to bams. Note that they are outside the scope of Dictys and we have not tried them. We do not support, endorse, or guarantee the validity of these approaches which are highly experimental in nature.

* **How do I perform network inference faster**?

  1. Get a GPU, such as:
  
     - `Google Colaboratory <https://colab.research.google.com/>`_ offers free GPU access with zero/minimal setup. You can run Dictys on very small datasets for free, or larger datasets with paid membership.
     - Major cloud computing service providers offer GPU access that is orders of magnitude cheaper than a scRNA-seq experiment.
     - High-performance computing cluster with GPU access at institution or other levels. Dedicated computing server. Personal computer with high-end consumer level GPU.
     - People or labs with the above access.
		
  2. Reduce the computational load, such as:
  
     - For context specific networks, choose only cell clusters of your interest. For this, delete the uninterested cell clusters in `data/subsets.txt`.
     - For dynamic networks, use fewer windows. This risks reducing time resolution. Details TBA.
     - Reduce the number of training steps. This risks reducing network quality. Details TBA.
		
  3. Configure properly for a powerful CPU. Details TBA.

* **Why do I see this error:** ``AssertionError: Torch not compiled with CUDA enabled``?
  
  This is because you installed a CPU-only pytorch but tried to run it on GPU. You have several options:
  
  1. To run pytorch on **CPU**, run ``dictys_helper makefile_update.py path/to/config.mk '{"DEVICE": "cpu"}'`` to configure to CPU mode. See `Tutorials`_ to find the right place to run this command.
  2. To run pytorch on **GPU**, reinstall Dictys with the correct options to enable GPU support at `Installation`_.

* **How do I use a large motif database where each motif can map to multiple TFs, such as from SCENIC+?**
  
  You need to first convert the motif database into a ``.motif`` file in `HOMER format <http://homer.ucsd.edu/homer/motif/creatingCustomMotifs.html>`_. Each motif should be named as ``TFNAME_unique-suffix`` where TFNAME should match the gene name in your dataset including capitalization. For multi-TF motifs, merge them as ``TFNAME1,TFNAME2,TFNAME3_unique-suffix`` for best speed, instead of duplicating them under each TF. See ``motifs.motif`` in tutorial inputs to understand file format. **Important**: the `log odds detection threshold <http://homer.ucsd.edu/homer/motif/creatingCustomMotifs.html>`_ column needs to be filled properly.

* **How do I use Dictys on multiple samples?**
  
  We used Dictys on multiple samples in the human blood dataset in `our paper <https://www.nature.com/articles/s41592-023-01971-3>`_ (Figs 2&5). However, we did not need to integrate multiple samples because the published dataset already did that. To check that on your own dataset, please see if samples display unintended separation in the low dimensions. If so, you may want to integrate them properly with any existing software before cell clustering or trajectory inference. These clusters or trajectories are inputs for GRN inference. In addition, you should try both including sample IDs as covariates in Dictys and not including them. We have not comprehensively tested covariate inclusion, so we suggest to choose the option that gives better biology for downstream analysis. See the `Optional: Prepare covariates <https://nbviewer.org/github/pinellolab/dictys/blob/master/doc/tutorials/full-skin/notebooks/main2.ipynb?flush_cache=true#Optional:-Prepare-covariates>`_ section in the full-skin tutorial on how to include covariates.

  To prepare input files for Dictys, please make sure each cell has a unique name across samples in the read count matrix and in the bam file. For read count matrices, you can append sample names to cell names before merging these matrices. For bam files, you can split each of them by cells using the script provided by Dictys in a separate folder for each sample, append sample names to the file names, and then move all the bam files into a single folder.

* **How do I save figures from jupyter notebooks onto the disk?**
  
  You can use ``plt.savefig('output.pdf')`` to save the current figure to disk. See `matplotlib.pyplot.savefig <https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.savefig.html#matplotlib.pyplot.savefig>`_.
  
  Some visualization functions in Dictys return two or more figures, such as ``figs = net.draw_discover(...)``. You can save them separately with ``figs[0].savefig('output1.pdf'); figs[1].savefig('output2.pdf'); ...``. See `matplotlib.figure.savefig <https://matplotlib.org/stable/api/figure_api.html#matplotlib.figure.Figure.savefig>`_ and `issue 15 <https://github.com/pinellolab/dictys/issues/15>`_.

Issues
==========================
Please raise an issue on `github <https://github.com/pinellolab/dictys/issues/new/choose>`_.

References
==========================
`Dictys: dynamic gene regulatory network dissects developmental continuum with single-cell multiomics <https://www.nature.com/articles/s41592-023-01971-3>`_ Nature Methods (2023)
