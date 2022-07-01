############################################################
# Run environment settings: Singularity
############################################################
#Singularity binary
SINGULARITY_BIN=singularity
#Custom parameters for `singularity exec`. Can be used to fix symlinks of data files.
SINGULARITY_PARAMS:=-B /data/pinello/PROJECTS/2020_08_dynet/dynet1/dat/00-raw/bams/granja2019-matched-sepagg:/data/pinello/PROJECTS/2020_08_dynet/dynet1/dat/00-raw/bams/granja2019-matched-sepagg:ro
#Path to singularity image
SINGULARITY_IMAGE=~/dictys3.sif
#Custom temporary folder can be specified here, or in environmental variable TMPDIR
SINGULARITY_TMPDIR=$(TMPDIR)

############################################################
# Below should not be changed
############################################################
#Conda environments in container
ENV_MAIN=main
ENV_PYRO=pyro
#Packages to map from host to container
SINGULARITY_PKGS_MAP=dictys docstring2argparse networkx
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

