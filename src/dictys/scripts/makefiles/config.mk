# Lingfei Wang, 2022. All rights reserved.
#This file contains parameters for whole run and individual steps to be edited for your dataset
#This file should be edited to configure the run
#This file should NOT be directly used for any run with `makefile -f` 

############################################################
# Run environment settings
############################################################
#Which environment to use, corresponding to env_$(ENVMODE).mk file
ENVMODE=none
#Maximum number of CPU threads for each job
#This is only nominated and passed through to other softwares without any guarantee.
NTH=4
#Device name for pyro/pytorch
#Note: cuda devices other than cuda:0 could be incompatible with singularity environment
DEVICE=cpu

############################################################
# Dataset settings
############################################################

#Genome size for Macs2, accept shortcuts like mm & hs
GENOME_MACS2=mm
#Whether dataset is joint profiling of RNA & ATAC of same cell. Separate measurements of two modalities in different cells: 0. Joint measurements: 1.
JOINT=1

#######################################################################################
# Parameters of individual steps shared by cell-type specific and dynamic GRN inference
#######################################################################################

PARAMS-PREPROC-QC_READS:=50 10 0 200 100 0
PARAMS-CHROMATIN-MACS2:=$(GENOME_MACS2)
PARAMS-CHROMATIN-BINLINKING:=20
KPARAMS-PREPROC-SELECTSC_RNA:=
KPARAMS-PREPROC-QC_READS:=
KPARAMS-PREPROC-SELECTSC_ATAC:=
KPARAMS-CHROMATIN-MACS2:=--nth $(NTH)
KPARAMS-CHROMATIN-WELLINGTON:=--nth $(NTH)
KPARAMS-CHROMATIN-HOMER:=--nth $(NTH)
KPARAMS-CHROMATIN-BINDING:=
KPARAMS-CHROMATIN-TSSDIST:=
KPARAMS-CHROMATIN-LINKING:=
KPARAMS-CHROMATIN-BINLINKING:=
KPARAMS-NETWORK-RECONSTRUCT:=--device $(DEVICE) --nth $(NTH)
KPARAMS-NETWORK-INDIRECT:=--nth $(NTH)
KPARAMS-NETWORK-NORMALIZE:=--nth $(NTH)
KPARAMS-NETWORK-TOFILE:=$(KPARAMS-NETWORK-TOFILE_EXTRA)

#######################################################################################
# Parameters of individual steps only in dynamic GRN inference
#######################################################################################

PARAMS-DYNAMIC-SUBSETS_RNA=1000 875 0.05
PARAMS-DYNAMIC-SUBSETS_ATAC=$(firstword $(PARAMS-DYNAMIC-SUBSETS_RNA))
