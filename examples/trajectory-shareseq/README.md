This example performs trajectory inference and format it for Dictys pipeline input for joint profiles of single-cell transcriptome and chromatin accessibility. It produces `traj_cell_rna.h5` and `traj_node.h5` files from `coord_rna.tsv.gz` file in the input `data` folder of the dictys pipeline of dynamic GRN reconstruction (examples/pipeline-shareseq).

The first step (`step1.ipynb`) performs trajectory inference (here using STREAM but others are also compatible) in the trajectory inference method's designated environment.

The second step (`step2.ipynb`) reformats the trajectory output for the pipeline with Dictys.

