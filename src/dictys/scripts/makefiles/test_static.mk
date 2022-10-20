# Lingfei Wang, 2022. All rights reserved.
#This file contains file I/O folder redirections for unit testing for cell-type specific GRN inference
#This file should NOT be edited to configure the run
#This file should be directly used for unit testing cell-type specific GRN inference with `makefile -f` 

#Working folder
DIRO:=test/output
DIRTI:=test/tmp_static_expected
DIRTO:=test/tmp_static
include $(dir $(lastword $(MAKEFILE_LIST)))static.mk
























#
