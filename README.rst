=========
Dictys
=========
Dictys reconstructs cell-type specific and dynamic gene regulatory networks (GRN) from scRNA-seq and scATAC-seq datasets. Dictys first infers a Transcription Factor (TF) binding network with TF footprinting from single-cell chromatin accessibility. Then Dictys refines the edges with single-cell transcriptome. Dictys addresses traditional challenges in network inference by orienting causality with TF binding information, modeling transcriptional rate to reconstruct cycle-compatible networks, and using probabilistic programming to capture the scRNA-seq process.

Dictys provides network analysis and visualization at global (across all cell types), pairwise (between two cell types) and single GRN levels. Dictys directly quantifies TF regulatory activity from GRN and enables a series of analyses such as cell-type specific TF discovery as regulation markers, differential regulation analysis alongside differential expression, and TF regulatory program illustration through its subnetwork and top activation/repression targets. These GRN-based analyses can capture unique biological insights not available from mean expression.

Dictys infers and analyzes dynamic GRN from scRNA-seq and scATAC-seq datasets along (infered) trajectories. This avoids artificial cell subsets and potential biases from population imbalance, and allows (pseudo-)time-resolved discovery and investigation of driver TFs and their individual regulations. Dictys provides an integrative network viewer for dynamic GRN visualization of synchronous panels in animation.

Overview
=============

.. image:: https://raw.githubusercontent.com/pinellolab/dictys/master/doc/images/logo.png
   :width: 180

.. image:: https://raw.githubusercontent.com/pinellolab/dictys/master/doc/images/Dictys_overview.png
   :width: 1000


Installation
=============
First install `Anaconda <https://www.anaconda.com/>`_. Then, install Dictys **with CPU computation**:

.. code-block::

	wget https://tinyurl.com/dictys -O - | bash

This will create a conda environment `dictys`.

Alternatively, under a different conda environment name:

.. code-block::

	wget https://tinyurl.com/dictys -O - | CONDAENV_NAME=your_favorite_name bash

Alternatively, **with GPU computation** (here with CUDA 11.3):

.. code-block::

	wget https://tinyurl.com/dictys -O - | CUDAVERSION_CONDA=11.3 bash

For more advanced installation, see `INSTALL.md <https://github.com/pinellolab/dictys/blob/master/INSTALL.md>`_ and/or edit `install.sh <https://github.com/pinellolab/dictys/blob/master/doc/scripts/install.sh>`_. *Note: dynamic network inference is computationally intensive and GPU availability is highly recommended.*

If you need `STREAM <https://github.com/pinellolab/STREAM>`_, `ArchR <https://www.archrproject.com/>`_, or other softwares upstream of Dictys, we recommend to install them in separate environments following their official instructions.

Examples
========
Dictys contains two major functions: network inference and network analysis. Network inference contains multiple command-line steps which are summarized in a `make` based pipeline. Network analysis can be achieved purely in jupyter notebooks. You can download the examples on `Zenodo <https://zenodo.org/record/6787658>`_.

Gallery
=======
The figures below are produced with the blood example dataset. You can reproduce them with the `analysis-blood` example. See `Examples`_. Each figure is linked to the jupyter notebook that produces it.

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

Issues
==========================
Please raise an issue on `github <https://github.com/pinellolab/dictys/issues/new>`_.

References
==========================
TBA
