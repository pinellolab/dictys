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
Installation should take ~<10 mins.

Option 1: with Anaconda
-----------------------
First install `Anaconda/Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_. Then, install Dictys **with CPU computation**:

.. code-block::

	conda create -y -n dictys -c conda-forge python=3.9 mamba
	. activate dictys
	mamba install -y -c lingfeiwang -c bioconda -c conda-forge -c pytorch dictys cpuonly

This will create a conda environment named ``dictys``.

Alternatively, **with GPU computation** (here CUDA 11.3):

.. code-block::

	conda create -y -n dictys -c conda-forge python=3.9 mamba
	. activate dictys
	mamba install -y -c lingfeiwang -c bioconda -c conda-forge -c pytorch dictys cudatoolkit=11.3

Option 2: with `bash script <https://tinyurl.com/dictys>`_
----------------------------------------------------------
First install `Anaconda/Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_. Then, install Dictys **with CPU computation**:

.. code-block::

	wget https://tinyurl.com/dictys -O - | bash

This will create a conda environment named `dictys`.

Alternatively, under a different conda environment name:

.. code-block::

	wget https://tinyurl.com/dictys -O - | CONDAENV_NAME=your_favorite_name bash

Alternatively, **with GPU computation** (here CUDA 11.3):

.. code-block::

	wget https://tinyurl.com/dictys -O - | CUDAVERSION_CONDA=11.3 bash

Option 3: with containers
-------------------------
TBA

Additional notes
----------------
For more advanced installation, see `INSTALL.md <https://github.com/pinellolab/dictys/blob/master/INSTALL.md>`_ and/or edit the `install script <https://tinyurl.com/dictys>`_.

*Note: dynamic network inference is computationally intensive and GPU availability is highly recommended.* Running time depends on the dataset, but it can take weeks or longer without a GPU.

If you need `STREAM <https://github.com/pinellolab/STREAM>`_, `ArchR <https://www.archrproject.com/>`_, or other softwares upstream of Dictys, we recommend to install them in separate environments following their official instructions.

Updating Dictys
----------------
If your minor version **is the latest** (e.g. your installed version is **0.1**.0 and the `latest release <https://github.com/pinellolab/dictys/releases>`_ is **0.1**.9), you can update Dictys to the latest github version with ``pip3 install --no-deps git+https://github.com/pinellolab/dictys`` inside your Dictys conda environment.

If your minor version **is not the latest** (e.g. your installed version is **0.1**.0 but the `latest release <https://github.com/pinellolab/dictys/releases>`_ is **0.2**.0), you should reinstall Dictys in a new conda environment with any option above.

Tutorials
=========
We provide several tutorials for different data types. Please download each tutorial folder structure before running. Note that these tutorials are not intended to replicate the results in the paper due to differences in software versions, computing platforms, various randomness e.g. in `HOMER genome preparsing <http://homer.ucsd.edu/homer/ngs/peakMotifs.html>`_ or `Pytorch algorithms <https://pytorch.org/docs/stable/notes/randomness.html>`_, etc.

1. `short-multiome <https://www.github.com/pinellolab/dictys/blob/master/doc/tutorials/short-multiome>`_: a single-notebook tutorial from data preparation to context specific network analysis on 10x multiome data for human blood.

2. `full-multiome <https://www.github.com/pinellolab/dictys/blob/master/doc/tutorials/full-multiome>`_: an extended version of the above tutorial with detailed usage.

The network analysis tutorials below use the same reconstructed networks as in the paper and are designed to replicate the results.

1. `analysis-blood <https://www.github.com/pinellolab/dictys/blob/master/doc/tutorials/analysis-blood>`_: a simple tutorial for context specific and dynamic network analysis on separate scRNA-seq and scATAC-seq quantifications of human blood as in manuscript.

We are organizing more tutorials for release. For now, you can explore some of them without technical support on `Zenodo <https://zenodo.org/record/6787658>`_ or `Google Colaboratory <https://colab.research.google.com/drive/1XJFpmAKzub-41QyoD6N_OGUgtbaGtU8g?usp=sharing>`_. Note these tutorials are subject to structural change.

Gallery
=======
The figures below are produced with the blood example dataset. You can reproduce them with the `analysis-blood` example. See `Tutorials`_. Each figure is linked to the jupyter notebook that produces it.

Cell-type specific GRN analyses
-------------------------------
`Regulation marker TF discovery <https://www.github.com/pinellolab/dictys/blob/master/doc/notebooks/static/global-dotplot.ipynb>`_

.. image:: https://raw.githubusercontent.com/pinellolab/dictys/master/doc/images/Global_dotplot.png
   :width: 300

`Top activation target heatmap for select TFs <https://www.github.com/pinellolab/dictys/blob/master/doc/notebooks/static/global-heatmap.ipynb>`_

.. image:: https://raw.githubusercontent.com/pinellolab/dictys/master/doc/images/Global_heatmap.png
   :width: 400

`Differential regulation v.s. differential expression scatter plot; integrative TF rank plot <https://www.github.com/pinellolab/dictys/blob/master/doc/notebooks/static/pair-diff.ipynb>`_

.. image:: https://raw.githubusercontent.com/pinellolab/dictys/master/doc/images/Diff_analysis.png
   :width: 750

`Subnetwork for select TF <https://www.github.com/pinellolab/dictys/blob/master/doc/notebooks/static/subnet.ipynb>`_

.. image:: https://raw.githubusercontent.com/pinellolab/dictys/master/doc/images/Subnet.png
   :width: 300
   
Dynamic GRN analysis
--------------------
`Driver TF discovery based on regulatory activity curve <https://www.github.com/pinellolab/dictys/blob/master/doc/notebooks/dynamic/discover.ipynb>`_

.. image:: https://raw.githubusercontent.com/pinellolab/dictys/master/doc/images/Dynamic_discovery.png
   :width: 1050

`Dynamic GRN animation <https://www.github.com/pinellolab/dictys/blob/master/doc/notebooks/dynamic/animation.ipynb>`_

.. image:: https://raw.githubusercontent.com/pinellolab/dictys/master/doc/images/animation.gif
   :width: 800

FAQ
==========================
* **How do I perform network inference faster**?

  1. Get a GPU, such as:
  
     - `Google Colaboratory <https://colab.research.google.com/>`_ offers free GPU access with zero/minimal setup. You can run Dictys on very small datasets for free, or larger datasets with paid membership. See `our tutorial <https://colab.research.google.com/drive/1XJFpmAKzub-41QyoD6N_OGUgtbaGtU8g?usp=sharing>`_.
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
  
  1. To run pytorch on CPU, run ``dictys_helper makefile_update.py path/to/config.mk '{"DEVICE": "cpu"}'`` to configure to CPU mode. See `Tutorials`_ to find the right place to run this command.
  2. To run pytorch on GPU, reinstall Dictys with the correct options to enable GPU support at `Installation`_.

Issues
==========================
Please raise an issue on `github <https://github.com/pinellolab/dictys/issues/new/choose>`_.

References
==========================
`Dictys: dynamic gene regulatory network dissects developmental continuum with single-cell multi-omics <https://www.biorxiv.org/content/10.1101/2022.09.14.508036>`_ bioRxiv (2022)
