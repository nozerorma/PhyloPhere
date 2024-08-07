#
#
#  ██████╗ ██╗  ██╗██╗   ██╗██╗      ██████╗ ██████╗ ██╗  ██╗███████╗██████╗ ███████╗
#  ██╔══██╗██║  ██║╚██╗ ██╔╝██║     ██╔═══██╗██╔══██╗██║  ██║██╔════╝██╔══██╗██╔════╝
#  ██████╔╝███████║ ╚████╔╝ ██║     ██║   ██║██████╔╝███████║█████╗  ██████╔╝█████╗  
#  ██╔═══╝ ██╔══██║  ╚██╔╝  ██║     ██║   ██║██╔═══╝ ██╔══██║██╔══╝  ██╔══██╗██╔══╝  
#  ██║     ██║  ██║   ██║   ███████╗╚██████╔╝██║     ██║  ██║███████╗██║  ██║███████╗
#  ╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝ ╚═════╝ ╚═╝     ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚══════╝
#                                                                                    
#                                      
# PHYLOPHERE: A Nextflow pipeline including a complete set
# of phylogenetic comparative tools and analyses for Phenome-Genome studies
#
# Github: https://github.com/nozerorma/PhyloPhere
#
# Author:         Miguel Ramon (miguel.ramon@upf.edu)
#
# File: Dockerfile
#

# Set image origin to micromamba for a stripped down, immutable Debian-Conda environment.
FROM mambaorg/micromamba:1.5.1

LABEL image.author.name "Miguel Ramon"
LABEL image.author.email "miguel.ramon@upf.edu"

# Set root as default user due to some funky behavour when using $MAMBA_USER.
USER root

# Download some basic debian tools + procps (Nextflow dependency).
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive \
    apt-get install -qq -y --no-install-recommends \
        software-properties-common \
        dirmngr \
        ed \
        gpg-agent \
		less \
		locales \
		vim-tiny \
		ca-certificates \
        procps

# Add environment file to the Docker image.
COPY env.yml /tmp/env.yml

# Create micromamba environment with the specifics from env.yml.
RUN micromamba install -y -n base -f /tmp/env.yml && \
    micromamba clean --all --yes

# Set CAASTools as WD and create build directory structure.
# Add all the requirements for CAAStools.
WORKDIR /ct

RUN mkdir -p ./modules ./scripts

ADD  modules/ ./modules/
ADD  scripts/ ./scripts/

# Add ct maintool to wd
ADD ct .

# Make ct and scripts executable
RUN chmod +x ./ct
RUN chmod +x ./scripts/*

# Add ct and conda binary to $PATH
ENV PATH=/ct:$PATH
ENV PATH /opt/conda/envs/base/bin:$PATH

# Set the image to be run inside base environment directly
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
