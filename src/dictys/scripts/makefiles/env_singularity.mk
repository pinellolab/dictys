# Lingfei Wang, 2022. All rights reserved.
#This file contains the definition for running environment through singularity
#To use this running environment, set ENVMODE=singularity in config.mk
#This file should be edited to configure the run (top part only)
#This file should NOT be directly used for any run with `makefile -f` 
#This file is still under development

############################################################
# Run environment settings: Singularity
############################################################
#Singularity binary
SINGULARITY_BIN=singularity
#Custom parameters for `singularity exec`
SINGULARITY_PARAMS:=
#Path to singularity image
SINGULARITY_IMAGE=path/to/image.sif
#Custom temporary folder can be specified here, or in environmental variable TMPDIR
SINGULARITY_TMPDIR=$(TMPDIR)

############################################################
# Below should not be changed
############################################################
ifeq (a$(wildcard $(SINGULARITY_IMAGE)),a)
$(error Singularity image file not found: $(SINGULARITY_IMAGE) )
endif

#Conda environments in container
ENV_MAIN=
ENV_PYRO=
#Packages to map from host to container
SINGULARITY_PKGS_MAP=
#Parameters for singularity
SINGULARITY_PARAMS+=-ec --nv
SINGULARITY_PARAMS+=-B $(abspath $(lastword $(MAKEFILE_LIST))):/home/docker/.local/bin/.placeholder:ro
SINGULARITY_PARAMS+=-B $(SINGULARITY_TMPDIR):/tmp:rw -B $(CURDIR):/dictys/work:rw
ifneq ($(words $(foreach p,$(SINGULARITY_PKGS_MAP),$(shell python -c "import $(p) as m;from os.path import dirname,abspath;print(dirname(abspath(m.__file__)))"))),$(words $(SINGULARITY_PKGS_MAP)))
$(error Mapped python packages not found)
endif
SINGULARITY_PARAMS+=$(foreach p,$(SINGULARITY_PKGS_MAP),-B $(shell python -c "import $(p) as m;from os.path import dirname,abspath;print(dirname(abspath(m.__file__)))"):/dictys/lib/$(p):ro)
SINGULARITY_CMD=$(SINGULARITY_BIN) exec $(SINGULARITY_PARAMS) $(SINGULARITY_IMAGE) /usr/bin/env HOME=/home/docker /home/docker/.local/bin/dictys_singularity

FULL_CMD=$(SINGULARITY_CMD)
FULL_TMPDIR=$(SINGULARITY_TMPDIR)





































#
