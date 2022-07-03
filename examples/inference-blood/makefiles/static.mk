############################################################
# Parameters of each step
############################################################

KPARAMS-NETWORK-TOFILE:=

############################################################
# Below should not be changed
############################################################

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

