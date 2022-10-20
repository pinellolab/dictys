# Lingfei Wang, 2022. All rights reserved.
#This file contains variables & targets only for dynamic GRN inference
#This file should NOT be edited to configure the run
#This file should be directly used for dynamic GRN inference with `makefile -f` 

DPRODUCT=$(DIRO)/dynamic.h5
KPARAMS-NETWORK-TOFILE_EXTRA=--dynamic

#Working folder
DIRT:=tmp_dynamic
FILE_SUBSET=$(DIRTI)/subsets.txt
PRODUCT_SUBSET=$(DIRTO)/subsets.txt $(DIRTO)/subset_locs.h5 $(DIRTO)/subset_edges.tsv.gz

include $(dir $(lastword $(MAKEFILE_LIST)))config.mk
include $(dir $(lastword $(MAKEFILE_LIST)))common.mk

ifeq ($(SHOW_CPU),1)

.PHONY:
subset: $(DIRI)/traj_node.h5 $(DIRI)/traj_cell_rna.h5 $(DIRI)/coord_rna.tsv.gz
	mkdir -p $(DIRTO)
	if ! [ -e $(firstword $(PRODUCT_SUBSET)) ] $(foreach VAR,$(wordlist 2,$(words $(PRODUCT_SUBSET)),$(PRODUCT_SUBSET)),|| ! [ -e $(VAR) ]); then $(FULL_CMD) $(ENV_MAIN) dynamic subsets_rna $(DIRI)/traj_node.h5 $(DIRI)/traj_cell_rna.h5 $(DIRI)/coord_rna.tsv.gz $(DIRTO)/subsets.txt $(DIRTO)/subset_locs.h5 $(DIRTO) $(DIRTO)/subset_edges.tsv.gz $(PARAMS-DYNAMIC-SUBSETS_RNA); fi

ifeq (a$(JOINT),a1)
$(DIRTO)/%/names_atac0.txt: $(DIRTI)/%/names_rna.txt
	cp $< $@
else
$(DIRTO)/%/names_atac0.txt: $(DIRI)/traj_node.h5 $(DIRI)/traj_cell_atac.h5 $(DIRI)/coord_atac.tsv.gz $(DIRTI)/subsets.txt $(DIRTI)/subset_locs.h5
	mkdir -p $(dir $@)
	$(FULL_CMD) $(ENV_MAIN) dynamic subsets_atac $^ $@ $* $(PARAMS-DYNAMIC-SUBSETS_ATAC)
endif

endif






























#
