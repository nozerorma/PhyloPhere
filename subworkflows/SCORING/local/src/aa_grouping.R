#!/usr/bin/env Rscript
# aa_grouping.R — Amino-acid grouping schemes (US, GS1-GS4) and a residue encoder.
# PhyloPhere | subworkflows/SCORING/local/src/
# Sourced by: nothing in the pipeline (scoring_compute.nf stages it next to scoring_compute.R)
# =============================================================================
# R copy of the five partitions defined in
#   subworkflows/CT_DISAMBIGUATION/local/src/biochem/grouping.py
# (the source of truth; its labels match this table residue by residue). The same
# partitions are defined in subworkflows/CT/local/modules/caas_id.py, after
# Chen & Zou 2025, Mol Ecol Resour 25(1):e70052, doi:10.1111/1755-0998.70052.
#
# The five schemes are independent partitions along different physicochemical axes,
# not a nested refinement. Any cross-scheme descriptor built on them must therefore
# be a set or a profile, never an ordinal level.
#
#   US   Classical CAAS. Strict amino-acid identity: each residue is its own group.
#   GS1  Coarse biochemical recoding, 6 groups. Dayhoff-like partition with no
#        literature source. Basis: broad side-chain family.
#   GS2  Side-chain dipole and volume, 7 groups. Yang 2010 from Shen 2007
#        (doi:10.2174/092986610791760306).
#   GS3  Polarity and volume, 6 groups. Zhang 2000 (doi:10.1007/s002399910007).
#        Basis: polar vs non-polar crossed with small vs large.
#   GS4  Fine-grained biochemical, 12 groups. Textbook functional bins; reactive or
#        structural residues (G, P, K, M, F) are kept apart.
#
# encode_aa_r(aa, scheme) returns the group label. An unknown residue passes through
# as the raw upper-cased residue, as path_scores.py::encode_aa does (it returns
# get_grouping_scheme(raw, scheme) or raw).
# =============================================================================

.AA_GROUPING_SCHEMES <- local({
  .aas <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  US <- setNames(.aas, .aas)  # identity mapping (classical CAAS)
  # GS1: CV // AGPS // NDQE // RHK // ILMFWY // T
  GS1 <- c(
    C="t", V="t",
    A="n", G="n", P="n", S="n",
    N="p", D="p", Q="p", E="p",
    R="b", H="b", K="b",
    I="h", L="h", M="h", F="h", W="h", Y="h",
    T="o"
  )
  # GS2: C // AGV // DE // NQHW // RK // ILFP // YMTS
  GS2 <- c(
    C="c",
    A="s", G="s", V="s",
    D="a", E="a",
    N="n", Q="n", H="n", W="n",
    R="b", K="b",
    I="h", L="h", F="h", P="h",
    Y="x", M="x", T="x", S="x"
  )
  # GS3: C // AGPST // NDQE // RHK // ILMV // FWY
  GS3 <- c(
    C="c",
    A="n", G="n", P="n", S="n", T="n",
    N="s", D="s", Q="s", E="s",
    R="b", H="b", K="b",
    I="l", L="l", M="l", V="l",
    F="g", W="g", Y="g"
  )
  # GS4: C // AILV // ST // NQ // DE // RH // G // P // K // M // F // WY
  GS4 <- c(
    C="c",
    A="h", I="h", L="h", V="h",
    S="o", T="o",
    N="p", Q="p",
    D="a", E="a",
    R="b", H="b",
    G="g",
    P="r",
    K="k",
    M="m",
    F="f",
    W="y", Y="y"
  )
  list(US = US, GS1 = GS1, GS2 = GS2, GS3 = GS3, GS4 = GS4)
})

#' The 5 scheme names, in the fixed profile order used by convergence_schemes.
AA_SCHEME_NAMES <- c("US", "GS1", "GS2", "GS3", "GS4")

#' Encode one amino acid in a grouping scheme.
#'
#' @param aa single-letter residue (case-insensitive); "" / NA -> NA.
#' @param scheme one of US/GS1/GS2/GS3/GS4 (case-insensitive). Unknown scheme
#'   returns the raw residue.
#' @return group label; unknown residue -> raw upper-cased residue (mirrors
#'   path_scores.py::encode_aa).
encode_aa_r <- function(aa, scheme) {
  if (is.null(aa) || length(aa) == 0) return(NA_character_)
  raw <- toupper(trimws(as.character(aa)[1]))
  if (is.na(raw) || !nzchar(raw)) return(NA_character_)
  sc <- toupper(trimws(as.character(scheme)[1]))
  if (is.na(sc) || sc == "US" || !nzchar(sc)) return(raw)
  tbl <- .AA_GROUPING_SCHEMES[[sc]]
  if (is.null(tbl)) return(raw)
  lab <- unname(tbl[raw])
  if (is.null(lab) || is.na(lab)) return(raw)
  lab
}
