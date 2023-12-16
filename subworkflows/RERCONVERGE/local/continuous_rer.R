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
# File: continuous_rer.R
#

# Set up variable to control command line arguments
args <- commandArgs(TRUE)

# Load libraries
library(dplyr)
library(RERconverge)

# Load traitfile
traitPath <- args[1]
load(traitPath) # As neoplasia_vector

# Load RER matrix
treePath <- args[2]
geneTrees <- readRDS(treePath)

# Convert the trait vector to paths comparable to the paths in the RER matrix.
charpaths <- char2Paths(trait_vector, geneTrees)
char2path_obj <- saveRDS(charpaths, args[3])

# Load RERs
rdsPath <- args[4]
traitRERw <- readRDS(rdsPath)

# Perform continuous RER analysis
res <- correlateWithContinuousPhenotype(traitRERw, charpaths, min.sp = args[6],
    winsorizeRER = args[7], winsorizetrait = args[8])

saveRDS(res, args[5])


# ## Visualize the results
# head(res[order(res$p.adj),])

# ## Check how P-values are sorted
# pdf(args[5])
# p1 <- hist(res$P, breaks=100, main = "P-value distribution in RERs")
# dev.off()

# ## Check how correlation is affected by individual clades
# ### Probably this should be performed in a streamlined way in a different script
# x <- charpaths
# y <- primRERw['A1BG',] # example
# pathNames <- namePathsWSpecies(geneTrees$masterTree) 

# pdf(args[6])
# names(y) <- pathNames
# plot(x,y, cex.axis=1, cex.lab=1, cex.main=1, xlab="Weight Change", 
#     ylab="Evolutionary Rate", main="Gene A1BG Pearson Correlation (Example)",
#     pch=19, cex=1, xlim=c(-2,2))

# text(x, y, labels = names(y), pos=4)
# abline(lm(y~x), col='red',lwd=3)
# dev.off()
