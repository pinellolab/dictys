# Lingfei Wang, 2022. All rights reserved.
#This file contains variables & targets shared between cell-type specific and dynamic GRN inference
#This file should NOT be edited to configure the run
#This file should NOT be directly used for any run with `makefile -f` 

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
KPARAMS-CHROMATIN-WELLINGTON+=--fi_blacklist $(DIRI)/blacklist.bed
endif
ifneq (a$(wildcard $(DIRI)/whitelist.bed),a)
KPARAMS-CHROMATIN-LINKING+=--fi_whitelist $(DIRI)/whitelist.bed
endif
ifneq (a$(wildcard $(DIRI)/covariate.tsv.gz),a)
KPARAMS-NETWORK-RECONSTRUCT+=--fi_cov $(DIRI)/covariate.tsv.gz
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

.PHONY: cpu gpu combine t clean distclean showcmd

combine: $(DPRODUCT)

cpu: $(PRODUCT_CPU)

gpu: $(PRODUCT_GPU)

t:
	@echo $(PRODUCT)
	@echo $(PRODUCT_TMP)
	@echo $(DPRODUCT)

showcmd:
	@echo $(FULL_CMD) $(ENV_MAIN)

clean:
	$(RM) $(PRODUCT_CPU) $(PRODUCT_GPU) $(PRODUCT_SUBSET) $(PRODUCT_TMP)
	
distclean: clean
	$(RM) $(DPRODUCT)

ifeq ($(SHOW_CPU),1)

$(DIRTO)/%/expression0.tsv.gz: $(DIRI)/expression.tsv.gz $(DIRTI)/%/names_rna.txt
	$(FULL_CMD) $(ENV_MAIN) preproc selects_rna $(KPARAMS-PREPROC-SELECTSC_RNA) $^ $@

$(DIRTO)/%/expression.tsv.gz: $(DIRTI)/%/expression0.tsv.gz
	$(FULL_CMD) $(ENV_MAIN) preproc qc_reads $(KPARAMS-PREPROC-QC_READS) $^ $@ $(PARAMS-PREPROC-QC_READS)

ifeq ($(JOINT),1)
$(DIRTO)/%/names_atac.txt: $(DIRTI)/%/expression.tsv.gz $(DIRTI)/%/names_atac0.txt
	$(FULL_CMD) $(ENV_MAIN) preproc selects_atac $(KPARAMS-PREPROC-SELECTS_ATAC) $^ $@
else
$(DIRTO)/%/names_atac.txt: $(DIRTI)/%/names_atac0.txt
	cp $< $@
endif

$(DIRTO)/%/reads.bam $(DIRTO)/%/reads.bai $(DIRTO)/%/peaks.bed: $(DIRTI)/%/names_atac.txt $(DIRI)/bams
	$(FULL_CMD) $(ENV_MAIN) chromatin macs2 $(KPARAMS-CHROMATIN-MACS2) $^ $(DIRTO)/$*/reads.bam $(DIRTO)/$*/reads.bai $(DIRTO)/$*/peaks.bed $(PARAMS-CHROMATIN-MACS2)
	
$(DIRTO)/%/footprints.bed: $(DIRTI)/%/reads.bam $(DIRTI)/%/reads.bai $(DIRTI)/%/peaks.bed
	$(FULL_CMD) $(ENV_MAIN) chromatin wellington $(KPARAMS-CHROMATIN-WELLINGTON) $^ $@

$(DIRTO)/%/motifs.bed $(DIRTO)/%/wellington.tsv.gz $(DIRTO)/%/homer.tsv.gz: $(DIRTI)/%/footprints.bed $(DIRI)/motifs.motif $(DIRI)/genome $(DIRTI)/%/expression.tsv.gz
	$(FULL_CMD) $(ENV_MAIN) chromatin homer $(KPARAMS-CHROMATIN-HOMER) $^ $(DIRTO)/$*/motifs.bed $(DIRTO)/$*/wellington.tsv.gz $(DIRTO)/$*/homer.tsv.gz
	
$(DIRTO)/%/binding.tsv.gz: $(DIRTI)/%/wellington.tsv.gz $(DIRTI)/%/homer.tsv.gz
	$(FULL_CMD) $(ENV_MAIN) chromatin binding $(KPARAMS-CHROMATIN-BINDING) $^ $@

$(DIRTO)/%/tssdist.tsv.gz: $(DIRTI)/%/expression.tsv.gz $(DIRTI)/%/wellington.tsv.gz $(DIRI)/gene.bed
	$(FULL_CMD) $(ENV_MAIN) chromatin tssdist $(KPARAMS-CHROMATIN-TSSDIST) $^ $@

$(DIRTO)/%/linking.tsv.gz: $(DIRTI)/%/binding.tsv.gz $(DIRTI)/%/tssdist.tsv.gz
	$(FULL_CMD) $(ENV_MAIN) chromatin linking $(KPARAMS-CHROMATIN-LINKING) $^ $@

$(DIRTO)/%/binlinking.tsv.gz: $(DIRTI)/%/linking.tsv.gz
	$(FULL_CMD) $(ENV_MAIN) chromatin binlinking $(KPARAMS-CHROMATIN-BINLINKING) $^ $@ $(PARAMS-CHROMATIN-BINLINKING)
	
$(DIRTO)/%/net_nweight.tsv.gz: $(DIRTI)/%/net_weight.tsv.gz $(DIRTI)/%/net_meanvar.tsv.gz $(DIRTI)/%/net_covfactor.tsv.gz
	$(FULL_CMD) $(ENV_MAIN) network normalize $(KPARAMS-NETWORK-NORMALIZE) $^ $@

$(DIRTO)/%/net_iweight.tsv.gz: $(DIRTI)/%/net_meanvar.tsv.gz $(DIRTI)/%/net_weight.tsv.gz $(DIRTI)/%/net_covfactor.tsv.gz
	$(FULL_CMD) $(ENV_MAIN) network indirect $(KPARAMS-NETWORK-INDIRECT) --fi_meanvar $^ $@
	
$(DIRTO)/%/net_inweight.tsv.gz: $(DIRTI)/%/net_iweight.tsv.gz $(DIRTI)/%/net_meanvar.tsv.gz $(DIRTI)/%/net_covfactor.tsv.gz
	$(FULL_CMD) $(ENV_MAIN) network normalize $(KPARAMS-NETWORK-NORMALIZE) $^ $@

endif

ifeq ($(SHOW_PYRO),1)

$(DIRTO)/%/net_weight.tsv.gz $(DIRTO)/%/net_meanvar.tsv.gz $(DIRTO)/%/net_covfactor.tsv.gz $(DIRTO)/%/net_loss.tsv.gz $(DIRTO)/%/net_stats.tsv.gz: $(DIRTI)/%/expression.tsv.gz $(DIRTI)/%/binlinking.tsv.gz
	$(FULL_CMD) $(ENV_PYRO) network reconstruct $(KPARAMS-NETWORK-RECONSTRUCT) $^ $(DIRTO)/$*/net_weight.tsv.gz $(DIRTO)/$*/net_meanvar.tsv.gz $(DIRTO)/$*/net_covfactor.tsv.gz $(DIRTO)/$*/net_loss.tsv.gz $(DIRTO)/$*/net_stats.tsv.gz

endif

$(DPRODUCT):
	mkdir -p $(dir $(DPRODUCT))
	$(FULL_CMD) $(ENV_MAIN) network tofile $(KPARAMS-NETWORK-TOFILE) $(DIRI) $(DIRTI) $(FILE_SUBSET) $@






































#
