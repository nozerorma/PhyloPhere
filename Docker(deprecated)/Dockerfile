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
FROM mambaorg/micromamba:2.5.0

LABEL image.author.name="Miguel Ramon"
LABEL image.author.email="miguel.ramon@upf.edu"

# Set root as default user due to some funky behavour when using $MAMBA_USER.
USER root

# Reset SHELL to avoid conflicts with micromamba's entrypoint during apt operations
SHELL ["/bin/bash", "-c"]

# Download some basic debian tools + procps (Nextflow dependency).
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

# Add project files to the Docker image (recursive, preserving directory layout).
# Set Phylophere as WD and create build directory structure.
RUN mkdir -p /home/mambauser/phylophere && chown -R $MAMBA_USER:$MAMBA_USER /home/mambauser/phylophere
WORKDIR /home/mambauser/phylophere

COPY --chown=$MAMBA_USER:$MAMBA_USER phylophere/ /home/mambauser/phylophere/

# Create micromamba environment with the specifics from env.yml.
RUN chmod +x /home/mambauser/phylophere/environment/install_env.sh
RUN bash /home/mambauser/phylophere/environment/install_env.sh /home/mambauser/phylophere/environment/phylophere.yml

# There's a high chance I need to make ct executable as docker does not seem to preserve permissions
RUN chmod +x /home/mambauser/phylophere/subworkflows/CT/local/ct

# Clone hyphy-analyses so that molerate.bf is available inside the container.
# Pinned to master; update the ref when a newer release is required.
RUN git clone --depth 1 --branch master \
        https://github.com/veg/hyphy-analyses.git /home/mambauser/hyphy-analyses && \
    chmod -R a+rX /home/mambauser/hyphy-analyses

# Add ct and conda binary to $PATH
ENV PATH=/home/mambauser/phylophere:$PATH
ENV PATH=/opt/conda/envs/phylophere/bin:$PATH

# Set the image to be run inside base environment directly
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
