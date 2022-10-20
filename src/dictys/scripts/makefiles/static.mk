# Lingfei Wang, 2022. All rights reserved.
#This file contains variables & targets only for cell-type specific GRN inference
#This file should NOT be edited to configure the run
#This file should be directly used for cell-type specific GRN inference with `makefile -f` 

#Working folder
DIRT:=tmp_static
DPRODUCT=$(DIRO)/static.h5

include $(dir $(lastword $(MAKEFILE_LIST)))config.mk
include $(dir $(lastword $(MAKEFILE_LIST)))common.mk

ifeq ($(SHOW_CPU),1)

$(DIRTO)/%/names_rna.txt: $(DIRI)/subsets/%/names_rna.txt
	mkdir -p $(dir $@)
	cp $< $@

ifeq (a$(JOINT),a1)
$(DIRTO)/%/names_atac0.txt: $(DIRTI)/%/names_rna.txt
	cp $< $@
else
$(DIRTO)/%/names_atac0.txt: $(DIRI)/subsets/%/names_atac.txt
	mkdir -p $(dir $@)
	cp $< $@
endif

endif





























#
