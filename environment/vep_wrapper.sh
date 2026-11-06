#!/bin/sh
# vep_wrapper.sh — Run a program of the "vep" conda environment with that environment first on the PATH.
# PhyloPhere | environment/
#
# Installed by the Dockerfile as /usr/local/bin/vep, vep_install, filter_vep, haplo, variant_recoder and vep_convert_cache.
# These programs start with "#!/usr/bin/env perl" and need the perl and the Perl libraries (DBI, Bio::EnsEMBL) of their own
# environment, which the system perl lacks; putting that environment on the global PATH would instead shadow the programs of
# the "phylophere" environment, so it is put first only for the program that is run.

PATH=/opt/conda/envs/vep/bin:$PATH
export PATH
exec "/opt/conda/envs/vep/bin/$(basename "$0")" "$@"
