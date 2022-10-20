# Lingfei Wang, 2022. All rights reserved.
#This file contains the definition for running environment of no special environment setting
#To use this running environment, set ENVMODE=none in config.mk
#This file should NOT be edited to configure the run
#This file should NOT be directly used for any run with `makefile -f` 

#Conda environments in container
ENV_MAIN=
ENV_PYRO=
FULL_CMD=OPENBLAS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_MAX_THREADS=1 NUMEXPR_MAX_THREADS=1 MKL_MAX_THREADS=1 python3 -m dictys

FULL_TMPDIR=$(TMPDIR)




































#
