# Examples
This folder contains examples for three steps:
1. Preparing trajectory input for dynamic GRN inference in `trajectory-*`
2. Cell-type specific and dynamic GRN inference in `inference-*`
3. Cell-type specific and dynamic GRN analysis in `analysis-*`

Four example datasets are included:
1. A small test dataset for cell-type specific GRN inference Erythroid and HSCs derived from dataset 2 in `inference-test`. Dynamic GRN is not availble.
2. Separate profiles of single-cell transcriptome and chromatin accessibility in human blood from [Granja et al 2019](https://doi.org/10.1038/s41587-019-0332-7) in `*-blood`. Because raw reads for human are protected data, we only provide intermediate files for GRN inference. Some steps in GRN inference are not included.
3. Joint profiles of single-cell transcriptome and chromatin accessibility in mouse skin from [Ma et al 2020](https://doi.org/10.1016/j.cell.2020.09.056) in `*-skin`.
4. A 10x multiome dataset from scratch in `*-multiome`.

For each example, the input files are in `*/data` folder and output files are in `*/output`. Below we go through each example one by one.

## Preparing trajectory input for dynamic GRN inference
Trajectory input is needed only for dynamic GRN inference. Generating trajectory input takes two steps: 1) trajectory inference with existing methods in their designated environment; and 2) formatting the trajectory for dynamic GRN inference input.

### Input
* `coord_rna.tsv.gz`: Cells' low-dimensional coordinates from their RNA profiles. If you only need trajectory inference and dynamic GRN reconstruction for a subset of cells, only include these cells in this file.

Only needed for separate profiles of transcriptome and chromatin accessibility of different cells:
* `atac_map.tsv.gz`: Mapping of each cell from scATAC-seq to their most similar cell in scRNA-seq.

### Output
* `traj_node.h5`: Inferred trajectory
* `traj_cell_rna.h5`: Location of each RNA-measured cell on the trajectory

For separate profiles of transcriptome and chromatin accessibility of different cells, these additional output files are produced:
* `traj_cell_atac.h5`: Location of each ATAC-measured cell on the trajectory
* `coord_atac.tsv.gz`: Cells' low-dimensional coordinates from their ATAC profiles

### Running
1. Trajectory inference is performed with [STREAM](TBA) in the examples in `trajectory-*/notebooks/step1.ipynb`. Other trajectory inference methods are also compatible. This step uses cells' low-dimensional coordinates from their RNA profiles (`coord_rna.tsv.gz`) to infer a trajectory consisted of nodes, edges, and cells' locations on the trajectory as intermediate files. Several different types of trajectory inference output are accepted. Here for STREAM, they include edges' terminal nodes (`trajectory-*/tmp/edge.tsv.gz`), the edge each cell is on (`trajectory-*/tmp/branch.tsv.gz`), and each cell's distance to every node (`trajectory-*/tmp/dist.tsv.gz`).

2. Formatting the trajectory involves assembling the above intermediate files into Dictys' class objects and exporting them as the output files. Because Dictys and trajectory inference method can have distinct running environments, this step is separated as `trajectory-*/notebooks/step2.ipynb`.

3. In your custom run, files in `output` folder of trajectory inference should be moved to the `data` folder of dynamic GRN inference. In this example, this is not needed because these files are already included the dynamic GRN inference examples.

## GRN inference
These examples cover context specific and dynamic GRN inference, for joint and separate profiles of single-cell transcriptome and chromatin accessibility.

### Input
* `bams`: Folder that contains one bam file for each (scATAC-seq) cell. File name should be cell name.
* `blacklist.bed` (optional): Bed file of regions to exclude in chromatin accessibility analyses
* `expression.tsv.gz`: Read count matrix of RNA-profiled cells
* `genome`: Folder containing reference genome in HOMER format. Copying it from the original location is recommended because HOMER creates preparsed files in this folder.
* `gene.bed`: Bed file of gene regions and strand information to locate transcription start sites
* `motifs.motif`: Motif file of all motifs in HOMER format. Motifs must be named as TF_... where TF is the TF gene name matching those in `expression.tsv.gz`. [Log odds detection threshold](http://homer.ucsd.edu/homer/motif/creatingCustomMotifs.html) must be valid. Motif file can be obtained from different motif databases, e.g. from [HOCOMOCO](https://hocomoco11.autosome.org/downloads_v11) or aggregated from HOMER.

Only needed for context-specific GRN inference:
* `subsets.txt`: Names of cell subsets/contexts to reconstruct GRNs
* `subsets`: Folder containing each cell subset as a subfolder. Each subfolder contains files `names_rna.txt` and `names_atac.txt` including cell names belonging to this subset from RNA or ATAC sequencing respectively. For joint measurements of RNA and ATAC, these two files should be identical.

Only needed for dynamic GRN inference:
* `coord_rna.tsv.gz`: Cells' low-dimensional coordinates from their RNA profiles.
* `traj_cell_rna.h5`: Location of each RNA-measured cell on the trajectory
* `traj_node.h5`: Trajectory

Only needed for dynamic GRN inference from separate profiles:
* `coord_atac.tsv.gz`: Cells' low-dimensional coordinates from their ATAC profiles
* `traj_cell_atac.h5`: Location of each ATAC-measured cell on the trajectory

### Output
The output file depends on the running mode.
* `static.h5`: Context specific GRNs
* `dynamic.h5`: Dynamic GRNs

### Determining the running environment
The simplest way to run Dictys is within the anaconda environment, whether the anaconda environment itself is inside a docker or singularity container, or directly on the machine. See [installation instructions](../INSTALL.md). Below it is assumed we are inside the anaconda environment for Dictys.

### Configuration
Configuration is needed before running GRN inference.
1. Edit `makefiles/config.mk` to determine number of threads for each process, pytorch device, genome size, and whether the data is a joint quantification of transcriptome and chromatin accessibility.
2. For context specific and dynamic GRN inference, edit their respective makefiles (`makefiles/static.mk` or `makefiles/dynamic.mk`) to change running parameters for each GRN inference step if needed.

### Setting parameters (optional)
If needed, you can change the running parameter for each GRN inference step in the above makefiles. To understand the meaning of the positional parameters or optional parameters, follow the steps below:
* Enter one `inference-*` folder to use
* Run `make -f makefiles/static.mk showcmd` to obtain the command for running Dictys. For example, we get `python3 -m dictys`.
* Find the positional or optional parameter you wish to change or understand (`PARAMS-*` or `KPARAMS-*` in `makefiles/*.mk`). Replace `-` in variable name with space, convert to lower case, and remove the first word. For example, we use `PARAMS-CHROMATIN-MACS2` which is converted to `chromatin macs2`.
* Run the command that combines the two strings above with a `-h` suffix to obtain help. In this example, it is `python3 -m dictys chromatin macs2 -h`. We get
```
usage: dictys chromatin macs2 [-h] [--qcut QCUT] [--nth NTH] [--nmax NMAX] fi_names fi_bam fo_bam fo_bai fo_bed genome_size

Peak calling using macs2. Needs bam files for each cell in a given folder.

positional arguments:
  fi_names     Path of input text file containing one sample/cell name per line for macs2 peak calling
  fi_bam       Path of input folder that contains each cell's bam file by name in fi_names
  fo_bam       Path of output bam file for select samples/cells
  fo_bai       Path of output bai file for select samples/cells
  fo_bed       Path of output bed file of peaks
  genome_size  Genome size input of macs2. Use shortcuts hs or mm for human or mouse.

optional arguments:
  -h, --help   show this help message and exit
  --qcut QCUT  Qvalue cutoff for macs2 (default: 0.05)
  --nth NTH    Number of threads (default: 1)
  --nmax NMAX  Maximum number of peaks to retain, ordered by macs2 score. Use 0 for no limit. (default: 500000)
```
* All file system related positional parameters are accounted for by the makefile pipeline. Only other positional parameters should be specified in `PARAMS-*` in the same order. Most optional parameters can be specified in `KPARAMS-*`.

### Prior checks
This stage make simple prior checks to ensure your inputs are in place and consistent with each other. See `notebooks/check.ipynb`.

### Running GRN inference pipeline
For context specific GRN inference:
```
# Three lines below assume pytorch runs on GPU as configured in makefiles/common.mk
# Run cpu part: TF binding network
# Here 4 is the number of parallel processes to run.
# The nominal total number of CPU cores used is up to 4 * 4 (from number of threads in makefiles/common.mk)
# Change it to fit your machine
make -f makefiles/static.mk -j 4 -k cpu
# Run gpu part: stochastic process network
# Here 2 is the number of parallel processes to run
make -f makefiles/static.mk -j 1 -k gpu
# Run cpu part: postprocessing
make -f makefiles/static.mk -j 4 -k cpu
# If using CPU for pytorch, the commented line below should replace the three lines above
# make -f makefiles/static.mk -j 4 -k cpu
# Save to single file
make -f makefiles/static.mk combine
# Optional cleanup of intermediate files
# make -f makefiles/static.mk clean
```

For dynamic GRN inference:
```
# Run cpu part: Subsetting cells with a moving window
# This step is already done and should be skipped for dynamic GRN inference from the blood example dataset because of protected data
make -f makefiles/dynamic.mk subset
# Three lines below assume pytorch runs on GPU as configured in makefiles/common.mk
# Run cpu part: TF binding network
make -f makefiles/dynamic.mk -j 4 -k cpu
# Run gpu part: stochastic process network
make -f makefiles/dynamic.mk -j 1 -k gpu
# Run cpu part: postprocessing
make -f makefiles/dynamic.mk -j 4 -k cpu
# If using CPU for pytorch, the commented line below should replace the three lines above
# make -f makefiles/dynamic.mk -j 4 -k cpu
# Save to single file
make -f makefiles/dynamic.mk combine
# Optional cleanup of intermediate files
# make -f makefiles/dynamic.mk clean
```

## GRN analysis
TBA

### Input
GRN file copy

3. In your custom run, files in `output` folder of trajectory inference should be moved to the `data` folder of dynamic GRN inference. In this example, this is not needed because these files are already included the dynamic GRN inference examples.

### Output
Most notebooks do not have any explicit output file. However, you can use additional commands to save the computed dataframes and drawn figures to external files.


