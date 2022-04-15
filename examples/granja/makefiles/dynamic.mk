
############################################################
# Parameters of each step
############################################################
PARAMS_SUBSETS_RNA=1000 875 0.05
PARAMS_SUBSET_ATAC=$(firstword $(PARAMS_SUBSETS_RNA))
KPARAMS_TOFILE:=

############################################################
# Below should not be changed
############################################################

DPRODUCT=$(DIRO)/dynamic.h5
KPARAMS_TOFILE+=--dynamic

#Working folder
DIRT:=tmp_dynamic
FILE_SUBSET:=$(DIRT)/subsets.txt
PRODUCT_SUBSET:=$(FILE_SUBSET) $(DIRT)/subset_locs.h5 $(DIRT)/subset_edges.tsv.gz

include $(dir $(lastword $(MAKEFILE_LIST)))common.mk

ifeq ($(SHOW_CPU),1)

$(PRODUCT_SUBSET) &: $(DIRI)/traj_node.h5 $(DIRI)/traj_cell_rna.h5 $(DIRI)/coord_rna.tsv.gz
	mkdir -p $(DIRTO)
	$(SINGULARITY_CMD) $(ENV_MAIN) dynamic subsets_rna $(DIRI)/traj_node.h5 $(DIRI)/traj_cell_rna.h5 $(DIRI)/coord_rna.tsv.gz $(DIRTO)/subsets.txt $(DIRTO)/subset_locs.h5 $(DIRTO) $(DIRTO)/subset_edges.tsv.gz $(PARAMS_SUBSETS_RNA)

ifeq (a$(JOINT),a1)
$(DIRTO)/%/names_atac0.txt: $(DIRTI)/%/names_rna.txt
	cp $< $@
else
$(DIRTO)/%/names_atac0.txt: $(DIRI)/traj_node.h5 $(DIRI)/traj_cell_atac.h5 $(DIRI)/coord_atac.tsv.gz $(DIRTI)/subsets.txt $(DIRTI)/subset_locs.h5
	mkdir -p $(dir $@)
	$(SINGULARITY_CMD) $(ENV_MAIN) dynamic subset_atac $^ $@ $* $(PARAMS_SUBSET_ATAC)
endif

endif






























#
