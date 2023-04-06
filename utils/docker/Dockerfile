# Lingfei Wang, 2022-2023. All rights reserved.
FROM continuumio/miniconda3
USER root
SHELL ["/bin/bash", "-c"]

#System update
RUN DEBIAN_FRONTEND=noninteractive apt-get update \
	&& DEBIAN_FRONTEND=noninteractive apt-get upgrade -y \
	&& DEBIAN_FRONTEND=noninteractive apt-get install -y curl gawk man pkg-config python3 python3-pip git wget zip unzip xzip \
	&& rm -Rf /var/lib/apt/lists/*

# Install dictys
# Name of conda environment to create
ARG CONDAENV_NAME=dictys
#Commit version to install. If empty, uses local version (./local).
ARG COMMIT_VERSION=master
# Python version
ARG PYTHONVERSION_CONDA=3.9
# CUDA version. When empty, uses CPU instead
ARG CUDAVERSION_CONDA=
COPY local /dictys/local
RUN cd /dictys \
	&& if [ "a${COMMIT_VERSION}" != "a" ]; then wget -O install.sh https://raw.githubusercontent.com/pinellolab/dictys/"${COMMIT_VERSION}"/doc/scripts/install.sh; localpath=""; else cp local/doc/scripts/install.sh ./; localpath="/dictys/local"; fi \
	&& chmod u+x install.sh \
	&& COMMIT_VERSION="${COMMIT_VERSION}" CONDAENV_NAME="${CONDAENV_NAME}" PYTHONVERSION_CONDA="${PYTHONVERSION_CONDA}" CUDAVERSION_CONDA="${CUDAVERSION_CONDA}" LOCAL_PATH="$localpath" ./install.sh \
	&& cd / \
	&& rm -Rf /dictys

#Create entry point
RUN echo '#!/bin/bash' > /usr/bin/run_dictys \
	&& echo ". activate ${CONDAENV_NAME}" >> /usr/bin/run_dictys \
	&& echo 'dictys "$@"' >> /usr/bin/run_dictys \
	&& chmod u+x /usr/bin/run_dictys

#ENTRYPOINT ["/usr/bin/run_dictys"]
	

CMD ["/bin/bash"]
