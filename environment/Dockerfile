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

# Build context: this directory (docker build -t phylophere:local environment/).
# The image holds the software environment only. The pipeline code is not copied into it: Nextflow runs the processes
# with the code of the checkout (baseDir), which the container runtimes mount.

# Set image origin to micromamba for a stripped down, immutable Debian-Conda environment.
FROM mambaorg/micromamba:2.5.0

LABEL image.author.name="Miguel Ramon"
LABEL image.author.email="miguel.ramon@upf.edu"

# Set root as default user due to some funky behavour when using $MAMBA_USER.
USER root

# Reset SHELL to avoid conflicts with micromamba's entrypoint during apt operations
SHELL ["/bin/bash", "-c"]

# Download some basic debian tools + procps (Nextflow dependency) + git (the RERconverge build clones it).
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive \
    apt-get install -qq -y --no-install-recommends \
        dirmngr \
        ed \
        gpg-agent \
        less \
        locales \
        ca-certificates \
        procps \
        git \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

# The environment recipe: the same script and files as a local installation (install_env.sh creates the conda
# environment "phylophere" from phylophere.yml, the environment "vep" from vep.yml, and builds the R packages
# conda does not provide).
WORKDIR /home/mambauser/environment
COPY --chown=$MAMBA_USER:$MAMBA_USER install_env.sh phylophere.yml vep.yml /home/mambauser/environment/
RUN chmod +x install_env.sh \
    && bash install_env.sh /home/mambauser/environment/phylophere.yml \
    && micromamba clean --all --yes

# _entrypoint.sh activates ${ENV_NAME:-base}: without this variable the base environment is activated, its bin comes
# first on the PATH, and python3 and Rscript would not be the ones of the phylophere environment.
ENV ENV_NAME=phylophere
# Without the entrypoint (apptainer exec) the environment is still the first on the PATH. The vep environment is last, so
# `vep` is found and nothing of that environment shadows a program of phylophere.
ENV PATH=/opt/conda/envs/phylophere/bin:$PATH:/opt/conda/envs/vep/bin

# Set the image to be run inside the phylophere environment directly
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
