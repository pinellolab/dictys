############################################################
# Parameters of each step
############################################################

KPARAMS_TOFILE:=

############################################################
# Below should not be changed
############################################################

#Working folder
DIRT:=tmp_static
DPRODUCT=$(DIRO)/static.h5

include Makefile.common

ifeq ($(SHOW_CPU),1)

$(DIRT)/%/names_rna.txt: $(DIRI)/subsets/%/names_rna.txt
	mkdir -p $(dir $@)
	cp $< $@

ifeq (a$(JOINT),a1)
$(DIRT)/%/names_atac0.txt: $(DIRT)/%/names_rna.txt
	cp $< $@
else
$(DIRT)/%/names_atac0.txt: $(DIRI)/subsets/%/names_atac.txt
	mkdir -p $(dir $@)
	cp $< $@
endif

endif





























#

