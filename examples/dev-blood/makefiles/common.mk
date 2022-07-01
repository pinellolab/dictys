
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
DEVICE=cuda:0

############################################################
# Dataset settings
############################################################

#Genome size for Macs2, accept shortcuts like mm & hs
GENOME_MACS2=hs
#Whether dataset is joint profiling of RNA & ATAC of same cell. Separate measurements of two modalities in different cells: 0. Joint measurements: 1.
JOINT=0

############################################################
# Parameters of each step
############################################################

PARAMS_QC_READS:=50 10 0 200 100 0
PARAMS_MACS2:=$(GENOME_MACS2)
PARAMS_BINLINKING:=20
KPARAMS_SELECTSC_RNA:=
KPARAMS_QC_READS:=
KPARAMS_SELECTSC_ATAC:=
KPARAMS_MACS2:=--nth $(NTH)
KPARAMS_WELLINGTON:=--nth $(NTH)
KPARAMS_HOMER:=--nth 1
KPARAMS_BINDING:=
KPARAMS_TSSDIST:=
KPARAMS_LINKING:=
KPARAMS_BINLINKING:=
KPARAMS_RECONSTRUCT:=--device $(DEVICE) --nth $(NTH)
KPARAMS_INDIRECT:=--nth $(NTH)
KPARAMS_NORMALIZE:=--nth $(NTH)

############################################################
# Below should not be changed
############################################################

############################################################
# Run environment settings
############################################################

#Input data folder
DIRI:=data
#Output data folder
ifeq (a$(DIRO),a)
DIRO:=output
endif
#Working input folder
ifeq (a$(DIRTI),a)
ifneq (a$(DIRT),a)
DIRTI:=$(DIRT)
else
$(error Specify working input folder with DIRTI or DIRT)
endif
endif
#Working output folder
ifeq (a$(DIRTO),a)
ifneq (a$(DIRT),a)
DIRTO:=$(DIRT)
else
$(error Specify working output folder with DIRTO or DIRT)
endif
endif
#File for cell subsets
FILE_SUBSET?=$(DIRI)/subsets.txt
#Temporary folder if not set
TMPDIR?=/tmp
#Environment specific makefile
include $(dir $(lastword $(MAKEFILE_LIST)))env_$(ENVMODE).mk
#Test temporary folder size
TMPDIR_SUFFICIENT=$(shell [ $$(df --block-size=1G --output=avail $(FULL_TMPDIR) | tail -n +2) -gt 99 ]; echo $$? )
ifneq ($(TMPDIR_SUFFICIENT),0)
$(warning WARNING: Temporary folder ($(FULL_TMPDIR)) available space >100GB recommended.)
endif

############################################################
# Parameters of each step
############################################################
ifneq (a$(wildcard $(DIRI)/blacklist.bed),a)
KPARAMS_WELLINGTON+=--fi_blacklist $(DIRI)/blacklist.bed
endif

############################################################
# Products by computing mode
############################################################

PRODUCT_NAMES_CPU:=names_rna.txt names_atac0.txt names_atac.txt expression0.tsv.gz expression.tsv.gz names_atac.txt reads.bam reads.bai peaks.bed footprints.bed motifs.bed homer.tsv.gz wellington.tsv.gz binding.tsv.gz tssdist.tsv.gz linking.tsv.gz binlinking.tsv.gz net_nweight.tsv.gz net_iweight.tsv.gz net_inweight.tsv.gz
PRODUCT_NAMES_GPU:=net_weight.tsv.gz net_meanvar.tsv.gz net_covfactor.tsv.gz net_loss.tsv.gz net_stats.tsv.gz
DPRODUCT_NAMES=

SHOW_CPU:=1
SHOW_PYRO:=1

SUBSETS:=$(shell cat $(FILE_SUBSET) 2> /dev/null | tr "\n" " ")
ifeq (a$(SUBSETS),a)
ifneq (a$(MAKECMDGOALS),asubset)
$(error Missing subsets.txt. Run `make subset` first for dynamic network.)
endif
endif

#Determine targets to make based on cpu/gpu mode
ifeq (a$(DEVICE),acpu)
PRODUCT_NAMES_CPU+=$(PRODUCT_NAMES_GPU)
PRODUCT_NAMES_GPU:=
else
ifneq (a$(firstword $(subst :, ,$(DEVICE))),acuda)
$(error This pipeline only supports device cpu and cuda*)
endif
ifeq (a$(MAKECMDGOALS),agpu)
SHOW_CPU:=0
endif
ifeq (a$(MAKECMDGOALS),acpu)
SHOW_PYRO:=0
endif
endif

PRODUCT_CPU=$(foreach c,$(SUBSETS),$(foreach p,$(PRODUCT_NAMES_CPU),$(DIRTO)/$(c)/$(p)))
PRODUCT_GPU=$(foreach c,$(SUBSETS),$(foreach p,$(PRODUCT_NAMES_GPU),$(DIRTO)/$(c)/$(p)))
#Temporary files by called programs
PRODUCT_TMP+=$(foreach c,$(SUBSETS),$(wildcard $(DIRTO)/$(c)/reads.bam.tmp.*.bam))

############################################################
# Recipes
############################################################

.PHONY: cpu gpu subset combine t clean distclean

combine: $(DPRODUCT)

cpu: $(PRODUCT_CPU)

gpu: $(PRODUCT_GPU)

subset: $(PRODUCT_SUBSET)

t:
	@echo $(PRODUCT)
	@echo $(PRODUCT_TMP)
	@echo $(DPRODUCT)

clean:
	$(RM) $(PRODUCT_CPU) $(PRODUCT_GPU) $(PRODUCT_SUBSET) $(PRODUCT_TMP)
	
distclean: clean
	$(RM) $(DPRODUCT)

ifeq ($(SHOW_CPU),1)

$(DIRTO)/%/expression0.tsv.gz: $(DIRI)/expression.tsv.gz $(DIRTI)/%/names_rna.txt
	$(FULL_CMD) $(ENV_MAIN) preproc selects_rna $(KPARAMS_SELECTSC_RNA) $^ $@

$(DIRTO)/%/expression.tsv.gz: $(DIRTI)/%/expression0.tsv.gz
	$(FULL_CMD) $(ENV_MAIN) preproc qc_reads $(KPARAMS_QC_READS) $^ $@ $(PARAMS_QC_READS)

ifeq ($(JOINT),1)
$(DIRTO)/%/names_atac.txt: $(DIRTI)/%/expression.tsv.gz $(DIRTI)/%/names_atac0.txt
	$(FULL_CMD) $(ENV_MAIN) preproc selects_atac $(KPARAMS_SELECTS_ATAC) $^ $@
else
$(DIRTO)/%/names_atac.txt: $(DIRTI)/%/names_atac0.txt
	cp $< $@
endif

$(DIRTO)/%/reads.bam $(DIRTO)/%/reads.bai $(DIRTO)/%/peaks.bed : $(DIRTI)/%/names_atac.txt $(DIRI)/bams
	$(FULL_CMD) $(ENV_MAIN) chromatin macs2 $(KPARAMS_MACS2) $^ $(DIRTO)/$*/reads.bam $(DIRTO)/$*/reads.bai $(DIRTO)/$*/peaks.bed $(PARAMS_MACS2)
	
$(DIRTO)/%/footprints.bed: $(DIRTI)/%/reads.bam $(DIRTI)/%/reads.bai $(DIRTI)/%/peaks.bed
	$(FULL_CMD) $(ENV_MAIN) chromatin wellington $(KPARAMS_WELLINGTON) $^ $@

.NOTPARALLEL: 
$(DIRTO)/%/motifs.bed $(DIRTO)/%/wellington.tsv.gz $(DIRTO)/%/homer.tsv.gz: $(DIRTI)/%/footprints.bed $(DIRI)/motifs.motif $(DIRI)/genome $(DIRTI)/%/expression.tsv.gz
	$(FULL_CMD) $(ENV_MAIN) chromatin homer $(KPARAMS_HOMER) $^ $(DIRTO)/$*/motifs.bed $(DIRTO)/$*/wellington.tsv.gz $(DIRTO)/$*/homer.tsv.gz
	
$(DIRTO)/%/binding.tsv.gz: $(DIRTI)/%/wellington.tsv.gz $(DIRTI)/%/homer.tsv.gz
	$(FULL_CMD) $(ENV_MAIN) chromatin binding $(KPARAMS_BINDING) $^ $@

$(DIRTO)/%/tssdist.tsv.gz: $(DIRTI)/%/expression.tsv.gz $(DIRTI)/%/wellington.tsv.gz $(DIRI)/gff.gff
	$(FULL_CMD) $(ENV_MAIN) chromatin tssdist $(KPARAMS_TSSDIST) $^ $@

$(DIRTO)/%/linking.tsv.gz: $(DIRTI)/%/binding.tsv.gz $(DIRTI)/%/tssdist.tsv.gz
	$(FULL_CMD) $(ENV_MAIN) chromatin linking $(KPARAMS_LINKING) $^ $@

$(DIRTO)/%/binlinking.tsv.gz: $(DIRTI)/%/linking.tsv.gz
	$(FULL_CMD) $(ENV_MAIN) chromatin binlinking $(KPARAMS_BINLINKING) $^ $@ $(PARAMS_BINLINKING)
	
$(DIRTO)/%/net_nweight.tsv.gz: $(DIRTI)/%/net_weight.tsv.gz $(DIRTI)/%/net_meanvar.tsv.gz $(DIRTI)/%/net_covfactor.tsv.gz
	$(FULL_CMD) $(ENV_MAIN) network normalize $(KPARAMS_NORMALIZE) $^ $@

$(DIRTO)/%/net_iweight.tsv.gz: $(DIRTI)/%/net_meanvar.tsv.gz $(DIRTI)/%/net_weight.tsv.gz $(DIRTI)/%/net_covfactor.tsv.gz
	$(FULL_CMD) $(ENV_MAIN) network indirect $(KPARAMS_INDIRECT) --fi_meanvar $^ $@
	
$(DIRTO)/%/net_inweight.tsv.gz: $(DIRTI)/%/net_iweight.tsv.gz $(DIRTI)/%/net_meanvar.tsv.gz $(DIRTI)/%/net_covfactor.tsv.gz
	$(FULL_CMD) $(ENV_MAIN) network normalize $(KPARAMS_NORMALIZE) $^ $@

endif

ifeq ($(SHOW_PYRO),1)

$(DIRTO)/%/net_weight.tsv.gz $(DIRTO)/%/net_meanvar.tsv.gz $(DIRTO)/%/net_covfactor.tsv.gz $(DIRTO)/%/net_loss.tsv.gz $(DIRTO)/%/net_stats.tsv.gz : $(DIRTI)/%/expression.tsv.gz $(DIRTI)/%/binlinking.tsv.gz
	$(FULL_CMD) $(ENV_PYRO) network reconstruct $(KPARAMS_RECONSTRUCT) $^ $(DIRTO)/$*/net_weight.tsv.gz $(DIRTO)/$*/net_meanvar.tsv.gz $(DIRTO)/$*/net_covfactor.tsv.gz $(DIRTO)/$*/net_loss.tsv.gz $(DIRTO)/$*/net_stats.tsv.gz

endif

$(DPRODUCT):
	$(FULL_CMD) $(ENV_MAIN) network tofile $(KPARAMS_TOFILE) $(DIRI) $(DIRTI) $(FILE_SUBSET) $@






































#

