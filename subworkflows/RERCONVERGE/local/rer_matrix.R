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
# File: rer_matrix.R
#

# Set up variable to control command line arguments
args <- commandArgs(TRUE)

# Load libraries
library(dplyr)
library(RERconverge)

# Load our traits
traitPath <- args[1]
load(traitPath) # As trait_vector

# Load our trees
treePath <- args[2]
geneTrees <- readRDS(treePath)

# Get residuals from our traitfile
min_sp <- if (length(args) >= 4) as.integer(args[4]) else 10L

# Verify species overlap between trait vector and master tree
master_species <- geneTrees$masterTree$tip.label
n_overlap <- length(intersect(names(trait_vector), master_species))
message(sprintf(
    "INFO rer_matrix: %d/%d trait species found in master tree (useSpecies filter active)",
    n_overlap, length(names(trait_vector))
))
if (n_overlap < min_sp) {
    stop(sprintf(
        "ERROR rer_matrix: only %d species overlap trait vector / master tree (min.sp = %d). Check species name formatting.",
        n_overlap, min_sp
    ))
}
traitRERw <- tryCatch(
    getAllResiduals(geneTrees, useSpecies=names(trait_vector),
        transform = "sqrt", weighted = T, scale = T, min.sp = min_sp),
    error = function(e) {
        message("WARNING: getAllResiduals with weighted=T failed: ", conditionMessage(e))
        message("Retrying with weighted=F ...")
        getAllResiduals(geneTrees, useSpecies=names(trait_vector),
            transform = "sqrt", weighted = F, scale = T, min.sp = min_sp)
    }
)

# Now we save our RERs

## Save to path
saveRDS(traitRERw, args[3])

### DONE ###
