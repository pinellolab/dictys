#Working folder
DIRO:=test/output
DIRTI:=test/tmp_static_expected
DIRTO:=test/tmp_static
include $(dir $(lastword $(MAKEFILE_LIST)))static.mk

$(DIRTO)/%/motifs.bed $(DIRTO)/%/wellington.tsv.gz $(DIRTO)/%/homer.tsv.gz : $(DIRTO)/%/footprints.bed $(DIRI)/motifs.motif $(DIRI)/genome $(DIRTI)/%/expression.tsv.gz
	$(SINGULARITY_CMD) $(ENV_MAIN) chromatin homer $(KPARAMS_HOMER) $^ $(DIRTO)/$*/motifs.bed $(DIRTO)/$*/wellington.tsv.gz $(DIRTO)/$*/homer.tsv.gz



























#

