#!/usr/bin/env Rscript
# aa_grouping.R — Amino-acid groupings (GS1-GS4) and the groupings a derived side agrees with.
# PhyloPhere | subworkflows/SCORING/local/src/
# Sourced by: scoring_compute.R (scoring_compute.nf stages it next to the script)
# =============================================================================
# The four partitions are independent classifications along different physicochemical axes, not a nested
# refinement (Chen & Zou 2025, Mol Ecol Resour 25(1):e70052, doi:10.1111/1755-0998.70052). They annotate a position;
# they never enter a score.
#
#   GS1  Coarse biochemical recoding, 6 groups. Dayhoff-like partition with no literature source. Basis: broad
#        side-chain family.
#   GS2  Side-chain dipole and volume, 7 groups. Yang 2010 from Shen 2007 (doi:10.2174/092986610791760306).
#   GS3  Polarity and volume, 6 groups. Zhang 2000 (doi:10.1007/s002399910007). Basis: polar vs non-polar crossed
#        with small vs large.
#   GS4  Fine-grained biochemical, 12 groups. Textbook functional bins; reactive or structural residues (G, P, K, M, F)
#        are kept apart.
#
# A position agrees with a grouping when its derived residues share one class and that class is not the class of
# any ancestral residue: the substitution crosses a class boundary and the derived side is uniform in the grouping.
# The residues are read as the VEP annotation reads them, so both consumers agree on what the derived side is.
# =============================================================================

.AA_GROUPINGS <- list(
  # GS1: CV // AGPS // NDQE // RHK // ILMFWY // T
  GS1 = c(
    C="t", V="t",
    A="n", G="n", P="n", S="n",
    N="p", D="p", Q="p", E="p",
    R="b", H="b", K="b",
    I="h", L="h", M="h", F="h", W="h", Y="h",
    T="o"
  ),
  # GS2: C // AGV // DE // NQHW // RK // ILFP // YMTS
  GS2 = c(
    C="c",
    A="s", G="s", V="s",
    D="a", E="a",
    N="n", Q="n", H="n", W="n",
    R="b", K="b",
    I="h", L="h", F="h", P="h",
    Y="x", M="x", T="x", S="x"
  ),
  # GS3: C // AGPST // NDQE // RHK // ILMV // FWY
  GS3 = c(
    C="c",
    A="n", G="n", P="n", S="n", T="n",
    N="s", D="s", Q="s", E="s",
    R="b", H="b", K="b",
    I="l", L="l", M="l", V="l",
    F="g", W="g", Y="g"
  ),
  # GS4: C // AILV // ST // NQ // DE // RH // G // P // K // M // F // WY
  GS4 = c(
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
)

#' Names of the groupings, in the order they are listed in `biochem_groupings`.
AA_GROUPING_NAMES <- names(.AA_GROUPINGS)

#' Residue sets of a position, read as the VEP annotation reads them (VEP/local/src/vep_common.py::anc_der_from_row).
#'
#' The residues are the ones the ASR inferred (`ancestral`, `derived`: comma-separated strings). A position without
#' either falls back to the `caas` pattern "<top residues>/<bottom residues>" read in the direction of its side: on side
#' "top" the bottom residues are ancestral and the top ones derived, on "bottom" the reverse. Residues the derived
#' set shares with the ancestral one are removed from it unless that would empty it.
#'
#' @return list(anc = character vector, der = character vector), upper-case single letters.
.anc_der_sets <- function(ancestral, derived, caas, side) {
  .set <- function(x) {
    v <- toupper(trimws(strsplit(as.character(x), ",", fixed = TRUE)[[1]]))
    unique(v[nzchar(v)])
  }
  anc <- .set(ancestral)
  der <- .set(derived)
  if (!length(anc) && !length(der)) {
    parts <- strsplit(as.character(caas), "/", fixed = TRUE)[[1]]
    .letters <- function(x) if (is.na(x) || !nzchar(x)) character(0) else unique(strsplit(gsub("[^A-Za-z]", "", toupper(x)), "")[[1]])
    top <- .letters(if (length(parts) >= 1) parts[1] else "")
    bot <- .letters(if (length(parts) >= 2) parts[2] else "")
    if (identical(side, "top"))         { anc <- bot; der <- top }
    else if (identical(side, "bottom")) { anc <- top; der <- bot }
    else                                { anc <- character(0); der <- union(top, bot) }
  }
  d <- setdiff(der, anc)
  list(anc = anc, der = if (length(d)) d else der)
}

#' Classes of a residue vector in one grouping; residues the grouping does not define (X, gaps) are left out.
.classes_of <- function(residues, grouping) {
  unique(unname(grouping[residues[residues %in% names(grouping)]]))
}

#' Groupings a derived side agrees with.
#'
#' A position agrees with a grouping when the derived residues share one class and that class is not the class of any
#' ancestral residue. The residue sets are those of `.anc_der_sets`.
#'
#' @param ancestral,derived,caas,side vectors with one element per position (the `ancestral_aa`, `derived_aa`, `caas` and
#'   `side` columns of position_scores.tsv).
#' @return character vector, one element per position: the names of the agreeing groupings joined by "+"
#'   ("GS2+GS4"), or "" when there is none.
derived_groupings <- function(ancestral, derived, caas, side) {
  stopifnot(length(derived) == length(ancestral), length(caas) == length(ancestral), length(side) == length(ancestral))
  vapply(seq_along(derived), function(i) {
    r <- .anc_der_sets(ancestral[i], derived[i], caas[i], side[i])
    hit <- vapply(.AA_GROUPINGS, function(g) {
      d <- .classes_of(r$der, g)
      a <- .classes_of(r$anc, g)
      length(d) == 1L && length(a) > 0L && !(d %in% a)
    }, logical(1))
    paste(names(.AA_GROUPINGS)[hit], collapse = "+")
  }, character(1))
}
