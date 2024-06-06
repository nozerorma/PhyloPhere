```
██████╗ ██╗  ██╗██╗   ██╗██╗      ██████╗ ██████╗ ██╗  ██╗███████╗██████╗ ███████╗
██╔══██╗██║  ██║╚██╗ ██╔╝██║     ██╔═══██╗██╔══██╗██║  ██║██╔════╝██╔══██╗██╔════╝
██████╔╝███████║ ╚████╔╝ ██║     ██║   ██║██████╔╝███████║█████╗  ██████╔╝█████╗  
██╔═══╝ ██╔══██║  ╚██╔╝  ██║     ██║   ██║██╔═══╝ ██╔══██║██╔══╝  ██╔══██╗██╔══╝  
██║     ██║  ██║   ██║   ███████╗╚██████╔╝██║     ██║  ██║███████╗██║  ██║███████╗
╚═╝     ╚═╝  ╚═╝   ╚═╝   ╚══════╝ ╚═════╝ ╚═╝     ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚══════╝
```
                                                                                                                                                                                                                               
 # PHYLOPHERE 
 ### A Nextflow pipeline including a complete set of phylogenetic comparative tools and analyses for Phenome-Genome studies

 **Github:** https://github.com/nozerorma/caastools/nf-phylophere

 **Author:** Miguel Ramon (miguel.ramon@upf.edu)

# NOTE THAT THIS README IS RATHER OUTDATED, FOR MORE DETAIL PLEASE CONTACT ME.

## Overview

PhyloPhere is a Nextflow pipeline designed to facilitate a comprehensive and reproducible analysis of genome-phenome relationships using various phylogenetic comparative methods (PCMs). The pipeline adheres to Nextflow's guidelines and incorporates the CAAStools toolset, enabling researchers to address diverse contrast scenarios with optimal workload distribution. Additionally, PhyloPhere partially integrates the RERConverge package, offering essential functionalities for RER analysis, including trait string, RER matrix, and the master tree derived from individual gene trees.

## Table of Contents

- [Overview](#overview)
- [Getting Started](#getting-started)
- [Configuration](#configuration)
- [Documentation](#documentation)
- [Contributing](#contributing)
- [License](#license)
- [Usage](#usage)


## Getting Started

To get started with PhyloPhere, follow these steps:

1. Clone the repository: `git clone https://github.com/your-username/phylophere.git`
2. Navigate to the repository: `cd phylophere`
3. Install Nextflow: [Nextflow Installation](https://www.nextflow.io/docs/latest/getstarted.html#installation)
4. Execute the pipeline: `nextflow main.nf`

## Configuration

The configuration files are located in the `conf` directory. The `base.config` file contains the basic configuration parameters, while the `modules.config` file defines module-specific configurations.

## Documentation

- [Authors](./docs/authors.txt)
- [Citations](./docs/CITATIONS.md)
- [Contributors](./docs/contributors.txt)
- [License](./docs/license.txt)
- [Support](./docs/support.txt)

## Contributing

If you would like to contribute to PhyloPhere, please follow the guidelines outlined in [CONTRIBUTING.md](CONTRIBUTING.md).

## License

PhyloPhere is licensed under the [MIT License](LICENSE). See the [LICENSE](LICENSE) file for details.

## Usage

### Docker

#### Constructing the CAAStools Docker Image

To build the CAAStools Docker image, execute the following command:

```bash
$ docker build -t phylophere

```

#### Executing CAAStools with a Mounted Volume
To run PhyloPhere with a mounted volume, use the following command:
```
$ docker run -it --volume <hostdir>:<containerdir> phylophere:latest bash
```

#### Launching CAAStools from an External Image

To launch CAAStools from an external image, use the following command:
```bash

$ docker run -it --volume .:/home docker.io/miralnso/phylophere:latest
```

### Nextflow
Pipeline is optimized to be run in HPC environments using either Docker, Singularity (or Apptainer) or Podman. Most common scenarios will involve the deployment of a Singularity instance (no need for root or daemon capabilities).
To run Phylophere in Singularity:
1. Make necessary modifications to nextflow.config, modules.config, and other relevant scripts according to the needs of your workflow.
2. Execute the pipeline:
```
# CAAS discovery

# If configurations are set in the files, execute as follows:
## Tower enables resource management, must be correctly set in config.
$ nextflow run main.nf -with-singularity -with-tower -profile singularity

# To override specific parameters:
## Note: The default output subdirectory structure is $workDir/results/<timestamp>/<tool>/output.out
$ nextflow run main.nf -with-singularity -with-tower -profile singularity --ct_tool <discovery,resample,bootstrap> --alignment <alignmentsheet_csv/alignment_dir> --output <outdir> --tree <nw_tree> --mode <mode> [...]

# Display help
$ nextflow un main.nf -with-singularity -with-tower -profile singularity --help (--ct_tool <tool>) # tool can be specified for individual use-case help

# RER in continuous mode
run main.nf -with-singularity -with-tower -profile singularity --rer_tool build_trait,continuous
```


It's recommended to initiate the pipeline from main.nf, which establishes the standard environment for tool execution. A help prompt can be shown by using the argument --help. Specific help prompts per tool can be called by using both the ct_tool and help arguments (i.e., --ct_tool discovery --help).

RER objects may also be build using RERConverge toolset by correctly setting the config file and running the *build_trait* option. Analysis for continuous trait correlation to RERs may be performed through *continuous* option in --rer_tool.

3. If running in an HPC environment using Slurm queue manager, configure the *nextflow.config* file accordingly, allocate resources as required in *base.config* and check the SBATCH sample files for an idea of how Nextflow and Slurm work together.

ie. 
```
#!/bin/bash
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
# Github: https://github.com/nozerorma/caastools/nf-phylophere
#
# Author:         Miguel Ramon (miguel.ramon@upf.edu)
#
# File: example-sbatch_discovery-array.sh
#
# DO NOT MODIFY THE QUEUE UNLESS YOU ARE COMPLETELY SURE ABOUT WHAT YOU ARE DOING

# Ensure the Slurm log directory exists
mkdir -p Slurm

# Function to submit the array job
submit_array_job() {
    sbatch --parsable <<'EOF'
#!/bin/bash
#SBATCH --job-name=phylophere-discovery
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=9
#SBATCH -p haswell
#SBATCH --mem=8G
#SBATCH -e Slurm/slurm-%A_%a.err
#SBATCH -o Slurm/slurm-%A_%a.out
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=mail@institution.edu
#SBATCH --array=1-8%8

# Define the directory where trait files are located
DATE_WITH_TIME=$(date "+%Y%m%d-%H%M%S") # add %3N as we want milliseconds too
TRAIT_DIR="/path/goes/here/scratch/user/traitdir"
RESULTS_DIR="/path/goes/here/user/outdir/CT/${DATE_WITH_TIME}"

# Define workdir
WORK_DIR="/path/goes/here/user/PhyloPhere-work"
mkdir -p $WORK_DIR

# Calculate file index
FILE_INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

# Assuming your trait files are named in a way that they can be sorted and indexed
## Modify according to the needs and folder structure, this assumes traitdir/subdir/traits, for traitdir/traits remove **
TRAIT_FILES=($(ls $TRAIT_DIR/**/*.tab))
TRAIT_FILE=${TRAIT_FILES[$FILE_INDEX]}

# Load nextflow module
module load Nextflow

# Run nextflow pipeline
nextflow run main.nf -with-singularity -with-tower -profile singularity \
    -w $WORK_DIR \
    --ct_tool discovery \
    --traitfile $TRAIT_FILE \
    --outdir $RESULTS_DIR \
    --maxbgmiss "0" \
    --maxfgmiss "0" \
    --maxmiss "0" \
    --patterns "1,2,3"
EOF
}

# Function to submit the cleanup job
submit_cleanup_job() {
    local array_job_id=$1
    sbatch --dependency=afterok:$array_job_id <<'EOF'
#!/bin/bash
#SBATCH --job-name=cleanup
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -p haswell
#SBATCH --mem=2G
#SBATCH -e Slurm/cleanup-%j.err
#SBATCH -o Slurm/cleanup-%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mail@institution.edu

rm -rf ${WORK_DIR}/*

EOF
}

# Main script logic
array_job_id=$(submit_array_job)
submit_cleanup_job $array_job_id

```

## Future steps

As of today, this workflow extends support for only continuous trait exploration.  By harmonizing multiple analytical approaches, Phylophere streamlines the study of complex genome-phenome analyses. This not only saves time and effort but also ensures methodological consistency throughout the analytical process. 


Looking ahead, the roadmap includes the integration of the complete RERConverge package and other commonly used PCMs, such as PGLS analysis. The idea is to further expand the pipeline's capabilities by incorporating subsequent steps in our downstream analysis,
specifically the binary trait exploration and a comprehensive functional analysis.

## Disclaimer
Note: The current pipeline is in its early stages and will undergo continuous enhancements. For substantial tasks, it's recommended to utilize CAAStools for CAAS and RERConverge for RER analysis. The Nextflow pipeline has been crafted adhering to the best practices of DSL2.



