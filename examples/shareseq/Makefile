#Working folder
DIRT:=tmp_static

include Makefile.common

ifeq ($(SHOW_CPU),1)

$(DIRT)/%/names_rna.txt: $(DIRI)/subsets/%/names_rna.txt $(DIRT)/% .placeholder
	cp $< $@

ifeq (a$(JOINT),a1)
$(DIRT)/%/names_atac0.txt: $(DIRT)/%/names_rna.txt
	cp $< $@
else
$(DIRT)/%/names_atac0.txt: $(DIRI)/subsets/%/names_atac.txt $(DIRT)/%
	cp $< $@
endif

endif
































#

