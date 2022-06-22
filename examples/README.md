1. Prepare input data in data/ folder.
1. Edit makefiles/common.mk for running parameters, including environment settings, dataset settings and each step's parameters.
2. Edit the correspond makefiles/env_*.mk for the parameters of the corresponding environment.
3. Run below.

Static:

#Run cpu part (TF binding network)
make -f makefiles/static.mk -j 32 -k cpu
#Run gpu part (stochastic process network)
make -f makefiles/static.mk -j 2 -k gpu
#Run cpu part (postprocessing)
make -f makefiles/static.mk -j 32 -k cpu
#Save to single file
make -f makefiles/static.mk combine
#Optional cleanup of intermediate files
make -f makefiles/static.mk clean

Dynamic:

#Run cpu part (subsetting)
make -f makefiles/dynamic.mk subset
#Run cpu part (TF binding network)
make -f makefiles/dynamic.mk -j 32 -k cpu
#Run gpu part (stochastic process network)
make -f makefiles/dynamic.mk -j 2 -k gpu
#Run cpu part (postprocessing)
make -f makefiles/dynamic.mk -j 32 -k cpu
#Save to single file
make -f makefiles/dynamic.mk combine
#Optional cleanup of intermediate files
make -f makefiles/dynamic.mk clean




