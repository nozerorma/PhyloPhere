#!/usr/bin/env Rscript
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
# File: permulations.R
#
# This script is part of CAASTOOLS.

# A Convergent Amino Acid Substitution identification
# and analysis toolbox
#
# Author:         Fabio Barteri (fabio.barteri@upf.edu)
#
# Contributors:   Alejandro Valenzuela (alejandro.valenzuela@upf.edu)
#                 Xavier Farré (xfarrer@igtp.cat),
#                 David de Juan (david.juan@upf.edu),
#                 Miguel Ramon (miguel.ramon@upf.edu).
#
# SCRIPT NAME: permulations.r
# DESCRIPTION: Permulation script from RERconverge
# DEPENDENCIES: modules in modules/simulation folder
# CALLED BY: simulatte.py


library(tibble)
library(readr)
library(ape)
library(geiger)
library(dplyr)


# Set of RERConverge functions used by CT Resample
simulatevec <- function(namedvec, treewithbranchlengths) {
  rm = ratematrix(treewithbranchlengths, namedvec)
  sims = sim.char(treewithbranchlengths, rm, nsim = 1)
  nam = rownames(sims)
  s = as.data.frame(sims)
  simulatedvec = s[, 1]
  names(simulatedvec) = nam
  vec = simulatedvec
  vec
}
## Simpermvec
simpermvec <- function(namedvec, treewithbranchlengths) {
  vec = simulatevec(namedvec, treewithbranchlengths)
  simsorted = sort(vec)
  realsorted = sort(namedvec)
  l = length(simsorted)
  c = 1
  while (c <= l) {
    simsorted[c] = realsorted[c]
    c = c + 1
  }
  simsorted
}

# Inputs

args = commandArgs(trailingOnly=TRUE)

tree <- args[1]
config.file <- args[2]
number.of.cycles <- args[3]
selection.strategy <- args[4]
phenotypes <- args[5]
outdir <- args[6]
chunk.size <- ifelse(length(args) >= 7, as.integer(args[7]), 500)
include.b0 <- TRUE
if (length(args) >= 8) {
  include.b0 <- tolower(as.character(args[8])) %in% c("1", "true", "t", "yes", "y")
}

# Create output directory if it doesn't exist
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
  write(paste("[INFO]", Sys.time(), "Created output directory:", outdir), stdout())
}

# tree <- "Data/5.Phylogeny/science.abn7829_data_s4.nex.pruned.tree"
# config.file <-"Out/2.CAAS/20240703-def/4.Traitfiles/df4_sp.tab"
# number.of.cycles <- "10"
# selection.strategy <- "phylogeny"
# phenotypes <- "Data/9.CAAS_bootstrap/traitfile.tab"
# outfile <- "Data/9.CAAS_bootstrap/permulated_traits.tab"

# Read the tree object.
#imported.tree <- read_file(tree)

tree.o <- read.tree(tree)
trait <- tree.o$tip.label
l <- length(trait)

# Read the config file

cfg <- read.table(config.file, sep ="\t", header = F)

foreground.df <- subset(cfg, cfg$V2 == "1")
background.df <- subset(cfg, cfg$V2 == "0")

foreground.size <- length(foreground.df$V1)
background.size <- length(background.df$V1)

foreground.species <- foreground.df$V1
background.species <- background.df$V1

# Read the phenotype file

phenotype.df <- read.table(phenotypes, sep = "\t")
foreground.values <- subset(phenotype.df, phenotype.df$V1 %in% foreground.species)$V2
background.values <- subset(phenotype.df, phenotype.df$V1 %in% background.species)$V2

# SEED AND PRUNE
starting.values <- phenotype.df$V2
all.species <- phenotype.df$V1

pruned.tree.o <- drop.tip(tree.o, setdiff(tree.o$tip.label, all.species))

names(starting.values) <- all.species

### SIMULATE
counter = 0
chunk.counter = 1
file.counter = 1

simulated.traits.df <- data.frame(matrix(ncol = 3, nrow = 0))

# Start timing
start.time <- Sys.time()
write(paste("[START]", start.time, "Beginning permulation generation..."), stdout())
write(paste("[INFO] Total cycles:", number.of.cycles, "| Chunk size:", chunk.size), stdout())

if (include.b0) {
  # Write b_0 (original trait configuration) as resample_000.tab
  write(paste("[INFO]", Sys.time(), "Writing b_0 (original trait configuration) to resample_000.tab"), stdout())
  b0.df <- data.frame(matrix(ncol = 3, nrow = 1))
  b0.cycle.tag <- "b_0"
  b0.fg.species.tag <- paste(foreground.species, collapse=",")
  b0.bg.species.tag <- paste(background.species, collapse=",")
  b0.outline <- c(b0.cycle.tag, b0.fg.species.tag, b0.bg.species.tag)
  b0.df[1,] <- b0.outline
  b0.filepath <- file.path(outdir, "resample_000.tab")
  write.table(b0.df, sep="\t", col.names=FALSE, row.names=FALSE, file=b0.filepath, quote=FALSE)
  write(paste("[COMPLETE] b_0 written to:", b0.filepath), stdout())
} else {
  write(paste("[INFO]", Sys.time(), "Skipping b_0 (include_b0=FALSE)"), stdout())
}


for (j in 1:as.integer(number.of.cycles)){
  counter = counter + 1
  cycle.tag = paste("b", as.character(counter), sep = "_")
  permulated_phenotype <- simpermvec(starting.values, pruned.tree.o)
  x <- enframe(permulated_phenotype)


  print("Using strategy: Random")
  # Select potential foreground and background species
  potential.fg.df <- subset(x, value %in% foreground.values)
  potential.bg.df <- subset(x, value %in% background.values)

  potential.fg <- potential.fg.df$name
  potential.bg <- potential.bg.df$name

  fg.species <- sample(potential.fg, foreground.size)
  bg.species <- sample(potential.bg, background.size)

  # Create the output
  fg.species.tag = paste(fg.species, collapse=",")
  bg.species.tag = paste(bg.species, collapse=",")

  outline = c(cycle.tag, fg.species.tag, bg.species.tag)
  simulated.traits.df[nrow(simulated.traits.df) +1,] <- outline
  
  chunk.counter = chunk.counter + 1
  
  # Write chunk to file when chunk size is reached or on final cycle
  if (chunk.counter > chunk.size || counter == as.integer(number.of.cycles)) {
    # Generate filename with zero-padded file number
    filename <- sprintf("resample_%03d.tab", file.counter)
    filepath <- file.path(outdir, filename)
    
    # Write the chunk
    write.table(simulated.traits.df, sep="\t", col.names=FALSE, row.names=FALSE, file=filepath, quote=FALSE)
    
    # Calculate progress and timing
    elapsed <- as.numeric(difftime(Sys.time(), start.time, units="secs"))
    pct.complete <- (counter / as.integer(number.of.cycles)) * 100
    cycles.per.sec <- counter / elapsed
    remaining.cycles <- as.integer(number.of.cycles) - counter
    eta.secs <- remaining.cycles / cycles.per.sec
    eta.mins <- eta.secs / 60
    
    logline <- sprintf("[%s] File %d: %s | Cycles %d-%d | Progress: %.1f%% | Elapsed: %.1f min | ETA: %.1f min",
                      format(Sys.time(), "%H:%M:%S"),
                      file.counter,
                      filename,
                      counter - nrow(simulated.traits.df) + 1,
                      counter,
                      pct.complete,
                      elapsed / 60,
                      eta.mins)
    write(logline, stdout())
    
    # Reset for next chunk
    simulated.traits.df <- data.frame(matrix(ncol = 3, nrow = 0))
    chunk.counter = 1
    file.counter = file.counter + 1
  }
}

# Final summary
end.time <- Sys.time()
total.elapsed <- as.numeric(difftime(end.time, start.time, units="mins"))
write(paste("[COMPLETE]", end.time, "|", number.of.cycles, "cycles in", round(total.elapsed, 2), "minutes"), stdout())
write(paste("[OUTPUT] Generated", file.counter, "permuted trait files + 1 original (b_0) =", file.counter, "total files in:", outdir), stdout())
