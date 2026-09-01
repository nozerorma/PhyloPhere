#!/usr/bin/env Rscript
# =============================================================================
# PHYLOPHERE: A Nextflow pipeline for Phenome-Genome studies
# File: subworkflows/CT/local/scripts/permulations.R
# =============================================================================
# Generates null trait permulations for CAAS significance calibration.
#
# Strategies (--perm_strategy):
#   FGBG   : plain label shuffle, no phylogenetic structure.
#   BM     : RERconverge permulation — Brownian simulation on the real tree,
#            rank-matched back onto the observed values (simpermvec).
#   lambda : same, but simulated on a tree rescaled by the ML Pagel's lambda of
#            the empirical trait. Rank-matching preserves the empirical marginal
#            exactly, so the ONLY thing the evolutionary model controls is the
#            phylogenetic dependence structure — which is exactly what lambda
#            parameterises. With lambda < 1 the nulls are less clade-clumped and
#            therefore closer to the observed contrast design.
#
# Pool harvesting: the observed design is not a random draw — it is the output of
# an optimiser (TRAIT_ANALYSIS pair selection) that maximises phylogenetic
# independence. Comparing it against unfiltered nulls is not like-for-like, so
# each draw is run through the same selection rule (lean_contrast_selector.R).
# Every accepted permulation carries EXACTLY N_pairs_obs pairs — the null must
# contain the same number of species as the observed design. The tiers relax only
# the independence requirement:
#   Tier 1 : all N_pairs_obs pairs have mod_dunn >= 1
#   Tier 2 : exactly one pair falls below mod_dunn 1 (fallback, only used if
#            Tier 1 cannot fill the pool within max_tries)
# The pool is assembled Tier 1 first, so Tier 2 records only ever appear when
# Tier 1 came up short — which keeps the downstream seeded random subset in
# caas_permulation.nf::SUBSET_RESAMPLE_PERMS Tier-1-dominant by construction.
#
# The species written to resample_*.tab ARE the selector's chosen pairs: FG holds
# the high-extreme member of each pair, BG the low-extreme member, in matched
# pair order.
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(geiger)
})

log_msg <- function(tag, ...) write(paste0("[", tag, "] ", format(Sys.time()), " ", paste0(...)), stdout())

# ── RERconverge permulation primitives ───────────────────────────────────────
simulatevec <- function(namedvec, treewithbranchlengths) {
  rm   <- ratematrix(treewithbranchlengths, namedvec)
  sims <- sim.char(treewithbranchlengths, rm, nsim = 1)
  setNames(as.data.frame(sims)[, 1], rownames(sims))
}

# Rank-matching: the simulated vector supplies the ORDERING, the observed vector
# supplies the VALUES. The permulated marginal distribution is therefore
# identical to the real one by construction.
simpermvec <- function(namedvec, treewithbranchlengths) {
  vec       <- simulatevec(namedvec, treewithbranchlengths)
  simsorted <- sort(vec)
  simsorted[] <- sort(namedvec)
  simsorted
}

# ── CLI ──────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("usage: permulations.R <tree> <config> <cycles> <strategy> <phenotypes> <outdir> ",
       "[chunk_size] [include_b0] [discrete_method] [max_tries] [pheno_col] ",
       "[bottom_quantile] [top_quantile] [n_col] [c_col] [resample_use_n] [trait_type]")
}

arg_or <- function(i, default, cast = as.character) {
  if (length(args) >= i && nzchar(args[i])) cast(args[i]) else default
}

tree.path          <- args[1]
config.file        <- args[2]
number.of.cycles   <- as.integer(args[3])   # size of the harvested pool
selection.strategy <- tolower(args[4])      # "bm" | "fgbg" | "lambda"
phenotypes         <- args[5]
outdir             <- args[6]
chunk.size         <- arg_or(7,  500L,       as.integer)
include.b0         <- tolower(arg_or(8, "true")) %in% c("1", "true", "t", "yes", "y")
discrete_method    <- arg_or(9,  "quintile")
max_tries          <- arg_or(10, 1000000L,   as.integer)
pheno_col_name     <- arg_or(11, "")
bottom_quantile    <- arg_or(12, 0.10,       as.numeric)
top_quantile       <- arg_or(13, 0.90,       as.numeric)
n_col              <- arg_or(14, "")   # denominator column (e.g. adult_necropsy_count)
c_col              <- arg_or(15, "")   # numerator column   (e.g. malignant_count)
resample_use_n     <- tolower(arg_or(16, "true")) %in% c("1", "true", "t", "yes", "y")
trait_type         <- tolower(arg_or(17, "auto"))   # "" / auto / ordinal / continuous

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
  log_msg("INFO", "Created output directory: ", outdir)
}

# ── Locate the lean selector next to this script ─────────────────────────────
# Rscript exposes the script path via `--file=` in the full commandArgs; the
# `sys.frames()$ofile` idiom only works for source()d files and is NULL here.
script_dir <- {
  full <- commandArgs(trailingOnly = FALSE)
  hit  <- grep("^--file=", full, value = TRUE)
  if (length(hit)) dirname(normalizePath(sub("^--file=", "", hit[1]))) else getwd()
}
lean_script <- file.path(script_dir, "lean_contrast_selector.R")

# ── Config: species | binary label | pair (cluster) id ───────────────────────
cfg <- read.table(config.file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
foreground.species <- cfg$V1[cfg$V2 == "1"]
background.species <- cfg$V1[cfg$V2 == "0"]

# N_pairs_obs is ALWAYS read from the V3 cluster ids that the contrast selection
# step wrote into the discovery config. It is deliberately not a parameter: a null
# filtered to a different pair count than the observed design would not be
# calibrating that design.
if (ncol(cfg) < 3) {
  stop("Discovery config '", config.file, "' has no V3 cluster-id column, so the ",
       "observed independent pair count cannot be determined. Permulations must ",
       "match the pair count produced by the contrast selection step.")
}
target_pairs <- suppressWarnings(max(as.integer(cfg$V3), na.rm = TRUE))
if (!is.finite(target_pairs) || target_pairs <= 0L) {
  stop("Could not parse a positive pair count from the config V3 column.")
}
n_fg <- sum(cfg$V2 == "1"); n_bg <- sum(cfg$V2 == "0")
if (n_fg != target_pairs || n_bg != target_pairs) {
  log_msg("WARN", sprintf(
    "Config has %d FG / %d BG species but %d pair ids — permulations will target %d pairs",
    n_fg, n_bg, target_pairs, target_pairs))
}
log_msg("INFO", sprintf("Observed independent pair count from config V3: N_pairs_obs = %d", target_pairs))

# ── Phenotype table ──────────────────────────────────────────────────────────
# Accepts both the headerless two-column traitfile the pipeline emits and a full
# annotated TSV. Sniffed rather than assumed: reading a headerless file with
# header=TRUE silently consumes the first species as a column name.
first_line <- readLines(phenotypes, n = 1L, warn = FALSE)
first_flds <- strsplit(first_line, "\t", fixed = TRUE)[[1]]
has_header <- length(first_flds) < 2 || is.na(suppressWarnings(as.numeric(first_flds[2])))

pheno_raw <- read.delim(phenotypes, sep = "\t", header = has_header,
                        stringsAsFactors = FALSE,
                        na.strings = c("", "NA", "NaN", "nan", "NULL", "null"))
log_msg("INFO", sprintf("Read %d phenotype rows from %s (header=%s)",
                        nrow(pheno_raw), basename(phenotypes), has_header))

sp_col <- if ("species" %in% names(pheno_raw)) "species" else names(pheno_raw)[1]
val_col <- if (nzchar(pheno_col_name) && pheno_col_name %in% names(pheno_raw)) {
  pheno_col_name
} else {
  cand <- setdiff(names(pheno_raw)[vapply(pheno_raw, is.numeric, logical(1))], sp_col)
  if (length(cand)) cand[1] else names(pheno_raw)[2]
}
log_msg("INFO", sprintf("Using species column '%s' and value column '%s'", sp_col, val_col))

phenotype.df <- data.frame(
  species = as.character(pheno_raw[[sp_col]]),
  value   = suppressWarnings(as.numeric(pheno_raw[[val_col]])),
  stringsAsFactors = FALSE
)

# ── Count data → Jeffreys CIs, so the permulation can apply the SAME candidate
# filter the observed contrast selection used. TRAIT_ANALYSIS/3.CI-composition.Rmd
# takes the count path when both n_trait and c_trait exist, and there a pair is a
# valid candidate iff its intervals do not overlap; otherwise it takes the discrete
# path, whose trait_overlap rule ("one top, one bottom") the percentile selection
# already reproduces exactly. Traits with no counts therefore need nothing here.
if (resample_use_n) {
  if (!nzchar(n_col) || !n_col %in% names(pheno_raw)) {
    cand_n <- c("n_trait", "n", "n_population", "sample_size", "total_count", "N")
    found_n <- intersect(cand_n, names(pheno_raw))
    if (length(found_n) > 0) n_col <- found_n[1]
  }
  if (!nzchar(c_col) || !c_col %in% names(pheno_raw)) {
    cand_c <- c("c_trait", "c", "c_cases", "n_cases", "cases", "C")
    found_c <- intersect(cand_c, names(pheno_raw))
    if (length(found_c) > 0) c_col <- found_c[1]
  }
}

use_ci <- resample_use_n && nzchar(n_col) && nzchar(c_col) &&
          n_col %in% names(pheno_raw) && c_col %in% names(pheno_raw)
if (use_ci) {
  if (!requireNamespace("binom", quietly = TRUE)) {
    stop("Count columns were supplied but the 'binom' package is unavailable, so the ",
         "Jeffreys CIs the observed selection used cannot be reproduced.")
  }
  np <- suppressWarnings(as.numeric(pheno_raw[[n_col]]))
  nc <- suppressWarnings(as.numeric(pheno_raw[[c_col]]))
  good <- is.finite(np) & is.finite(nc) & np > 0 & nc >= 0 & nc <= np
  phenotype.df$ci_lb <- NA_real_; phenotype.df$ci_ub <- NA_real_
  if (any(good)) {
    b <- binom::binom.confint(nc[good], np[good], method = "bayes",
                              priors = c(0.5, 0.5), conf.level = 0.95)
    phenotype.df$ci_lb[good] <- b[, 5]; phenotype.df$ci_ub[good] <- b[, 6]
  }
  log_msg("INFO", sprintf("Count columns '%s'/'%s' found (resample_use_n=TRUE): %d/%d species have Jeffreys CIs — ",
                          n_col, c_col, sum(good), length(good)),
          "candidate pairs will be filtered on CI non-overlap (count path)")
} else {
  if (!resample_use_n) {
    log_msg("INFO", "resample_use_n is FALSE in config: using the percentile/top-bottom candidate rule (discrete path)")
  } else {
    log_msg("INFO", "No usable count columns found: using the percentile/top-bottom candidate rule ",
            "(the discrete path's trait_overlap)")
  }
}

keep <- !is.na(phenotype.df$species) & is.finite(phenotype.df$value)
if (use_ci) keep <- keep & is.finite(phenotype.df$ci_lb) & is.finite(phenotype.df$ci_ub)
phenotype.df <- phenotype.df[keep, , drop = FALSE]
if (nrow(phenotype.df) == 0) stop("No valid phenotype rows remain after NA filtering")

# ── Tree ─────────────────────────────────────────────────────────────────────
tree.o <- read.tree(tree.path)
starting.values <- setNames(phenotype.df$value, phenotype.df$species)

pruned.tree <- drop.tip(tree.o, setdiff(tree.o$tip.label, names(starting.values)))
pruned.tree <- multi2di(pruned.tree)
pruned.tree$edge.length[pruned.tree$edge.length <= 0] <- 1e-8
starting.values <- starting.values[pruned.tree$tip.label]
if (any(is.na(starting.values))) stop("Internal error: NA trait values after pruning to tree tips")

missing_fg <- setdiff(foreground.species, pruned.tree$tip.label)
missing_bg <- setdiff(background.species, pruned.tree$tip.label)
if (length(missing_fg) || length(missing_bg)) {
  log_msg("WARN", sprintf("%d FG and %d BG config species absent from the pruned tree",
                          length(missing_fg), length(missing_bg)))
}
log_msg("INFO", sprintf("Pruned tree: %d tips | target_pairs = %d | discrete_method = %s",
                        Ntip(pruned.tree), target_pairs, discrete_method))

# Distances are measured on the REAL tree, always — lambda rescaling belongs to
# the generative model, not to the geometry the Dunn index is evaluated in.
D <- cophenetic(pruned.tree)

# ── Simulation tree (lambda rescaling) ───────────────────────────────────────
simulation_tree <- pruned.tree
lambda_hat <- NA_real_
if (selection.strategy == "lambda") {
  log_msg("INFO", "Fitting ML Pagel's lambda on the empirical trait vector...")
  fitl <- fitContinuous(pruned.tree, starting.values, model = "lambda")
  lambda_hat <- as.numeric(fitl$opt$lambda)
  log_msg("INFO", sprintf("ML Pagel's lambda = %.4f", lambda_hat))
  simulation_tree <- rescale(pruned.tree, "lambda", lambda_hat)
  simulation_tree$edge.length[simulation_tree$edge.length <= 0] <- 1e-8
}

# ── Lean selector ────────────────────────────────────────────────────────────
if (!file.exists(lean_script)) {
  stop("lean_contrast_selector.R not found next to permulations.R (looked in '",
       script_dir, "'). Pair-count and independence filtering cannot be applied.")
}
source(lean_script)
log_msg("INFO", "Lean contrast filtering enabled from ", lean_script)

# ── b_0: the real labeling ───────────────────────────────────────────────────
if (include.b0) {
  write.table(
    data.frame("b_0", paste(foreground.species, collapse = ","),
               paste(background.species, collapse = ",")),
    file = file.path(outdir, "resample_000.tab"),
    sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE
  )
  log_msg("INFO", "Wrote b_0 (original trait configuration) to resample_000.tab")
}

# ── Harvest ──────────────────────────────────────────────────────────────────
# Three quantities govern this step and they are deliberately distinct:
#   max_tries      — draw budget (escalated below if the pool comes up short)
#   pool_size      — accepted permulations to harvest; the position-level null
#   caas_full_perms— drawn downstream from the pool for the FCS null (not seen here)
#
# If the budget is exhausted before the pool is full, the budget is raised by 50%
# and harvesting RESUMES (accepted records are kept, not discarded). After two
# escalations the run fails loudly rather than silently producing a thin null.
start.time <- Sys.time()
log_msg("START", sprintf("Harvesting pool (strategy: %s, pool_size: %d, max_tries: %d)",
                         selection.strategy, number.of.cycles, max_tries))

# Phenotype records in ascending value order — the order simpermvec's rank-matching
# assigns them in. rank_order is the identity on that table, kept explicit so the
# BM/lambda and FGBG branches below index it the same way.
rec_ord    <- order(starting.values)
rec_value  <- unname(starting.values[rec_ord])
n_rec      <- length(rec_value)
rank_order <- seq_len(n_rec)
rec_lb <- rec_ub <- NULL
if (use_ci) {
  ci_map <- setNames(phenotype.df$ci_lb, phenotype.df$species)
  rec_lb <- unname(ci_map[names(starting.values)[rec_ord]])
  ci_map <- setNames(phenotype.df$ci_ub, phenotype.df$species)
  rec_ub <- unname(ci_map[names(starting.values)[rec_ord]])
}

tier1 <- vector("list", number.of.cycles)
tier2 <- vector("list", number.of.cycles)
n1 <- 0L; n2 <- 0L
total_draws <- 0L
reject_reasons <- character(0)

MAX_ESCALATIONS  <- 2L
ESCALATION_FACTOR <- 1.5
budget      <- max_tries
escalations <- 0L

repeat {
  while (n1 < number.of.cycles && total_draws < budget) {
    total_draws <- total_draws + 1L

    # A permulation moves whole PHENOTYPE RECORDS between tips, not bare numbers:
    # when a tip receives another species' prevalence it receives that prevalence's
    # uncertainty too, so its confidence interval travels with it. simpermvec
    # returns tips ordered by their simulated value with the sorted observed values
    # assigned by rank, so rank i carries record i of the value-sorted table.
    ord <- NULL
    pvec <- if (selection.strategy == "fgbg") {
      ord <- sample.int(n_rec)
      setNames(rec_value[ord], names(starting.values))
    } else {
      v <- simpermvec(starting.values, simulation_tree)
      ord <- rank_order                       # ascending value order = record order
      v
    }

    ci_lb_draw <- NULL; ci_ub_draw <- NULL
    if (use_ci) {
      ci_lb_draw <- setNames(rec_lb[ord], names(pvec))
      ci_ub_draw <- setNames(rec_ub[ord], names(pvec))
    }

    e <- evaluate_lean_contrast_selection(
      trait_vec = pvec, D = D, target_pairs = target_pairs,
      discrete_method = discrete_method,
      bottom_quantile = bottom_quantile, top_quantile = top_quantile,
      ci_lb = ci_lb_draw, ci_ub = ci_ub_draw,
      trait_type = trait_type
    )

    if (e$tier == 1L) {
      n1 <- n1 + 1L; tier1[[n1]] <- e
    } else if (e$tier == 2L && n2 < number.of.cycles) {
      n2 <- n2 + 1L; tier2[[n2]] <- e
    } else if (e$tier == 0L) {
      reject_reasons <- c(reject_reasons, e$reason)
    }

    if (total_draws %% 5000 == 0 || n1 == number.of.cycles) {
      el <- as.numeric(difftime(Sys.time(), start.time, units = "secs"))
      log_msg("PROGRESS", sprintf("draws=%d/%d | Tier1=%d/%d | Tier2=%d | %.0f draws/s | %.1f min",
                                  total_draws, budget, n1, number.of.cycles, n2,
                                  total_draws / el, el / 60))
    }
  }

  if (n1 >= number.of.cycles || escalations >= MAX_ESCALATIONS) break

  escalations <- escalations + 1L
  budget <- as.integer(ceiling(budget * ESCALATION_FACTOR))
  log_msg("WARN", sprintf(
    "Pool not filled from Tier 1 (%d/%d) after %d draws — escalating max_tries by %.0f%% to %d (escalation %d of %d)",
    n1, number.of.cycles, total_draws, 100 * (ESCALATION_FACTOR - 1), budget, escalations, MAX_ESCALATIONS))
}

# ── Assemble the pool: Tier 1 first, Tier 2 only to fill a shortfall ─────────
if (n1 >= number.of.cycles) {
  pool <- tier1[seq_len(number.of.cycles)]
  log_msg("INFO", sprintf("Pool filled entirely from Tier 1 (%d/%d) in %d draws%s",
                          number.of.cycles, number.of.cycles, total_draws,
                          if (escalations) sprintf(" after %d escalation(s)", escalations) else ""))
} else if (n1 + n2 >= number.of.cycles) {
  use2 <- number.of.cycles - n1
  pool <- c(tier1[seq_len(n1)], tier2[seq_len(use2)])
  log_msg("WARN", sprintf(
    "Tier 1 exhausted after %d draws and %d escalation(s): topping up with Tier 2 (one pair below Dunn 1). %d Tier-1 + %d Tier-2 = %d/%d records",
    total_draws, escalations, n1, use2, length(pool), number.of.cycles))
} else {
  tab <- sort(table(reject_reasons), decreasing = TRUE)
  top <- seq_len(min(3, length(tab)))
  stop(sprintf(
    paste0("Permulation pool could not be filled: %d Tier-1 + %d Tier-2 = %d of the requested %d ",
           "after %d draws and %d escalation(s) of max_tries (final budget %d).\n",
           "  The downstream FCS null (caas_full_perms) is drawn from this pool, so continuing ",
           "would silently under-power it.\n",
           "  Top rejection reasons: %s\n",
           "  Consider: %sa larger --max_tries, or a smaller --perm_pool_size."),
    n1, n2, n1 + n2, number.of.cycles, total_draws, escalations, budget,
    if (length(tab)) paste(sprintf("%s (%d)", names(tab)[top], as.integer(tab)[top]), collapse = "; ") else "none recorded",
    if (selection.strategy != "lambda") "--perm_strategy lambda (accepts far more designs than BM), " else ""))
}

# ── Write resample chunks + a tier/Dunn manifest ─────────────────────────────
rows      <- vector("list", length(pool))
manifest  <- vector("list", length(pool))
file.counter <- 1L; chunk.start <- 1L

flush_chunk <- function(from, to) {
  fp <- file.path(outdir, sprintf("resample_%03d.tab", file.counter))
  write.table(do.call(rbind, rows[from:to]), file = fp, sep = "\t",
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  file.counter <<- file.counter + 1L
}

for (b in seq_along(pool)) {
  e <- pool[[b]]
  rows[[b]] <- data.frame(
    cycle = paste0("b_", b),
    fg    = paste(e$fg, collapse = ","),
    bg    = paste(e$bg, collapse = ","),
    stringsAsFactors = FALSE
  )
  # fg_values/bg_values are the PERMULATED trait values of the written species,
  # with the quantile cut-offs they were selected against — so the pool can be
  # audited: every fg value >= q_upper, every bg value <= q_lower, in every row.
  # Written at full precision on purpose: a quantile cut-off is itself an
  # observed trait value, so a species sitting exactly ON the threshold is the
  # common case, and rounding the value but not the threshold makes those rows
  # look like violations.
  fmt <- function(x) paste(format(x, digits = 15, trim = TRUE), collapse = ",")
  manifest[[b]] <- data.frame(cycle = paste0("b_", b), tier = e$tier,
                              n_pairs = e$n_pairs, dunn_min = e$dunn_min,
                              n_pairs_below_1 = e$n_below, mode = e$mode,
                              q_lower = e$q_lower, q_upper = e$q_upper,
                              fg_values = fmt(e$fg_values),
                              bg_values = fmt(e$bg_values),
                              stringsAsFactors = FALSE)
  if (b - chunk.start + 1L >= chunk.size || b == length(pool)) {
    flush_chunk(chunk.start, b)
    chunk.start <- b + 1L
  }
}

write.table(do.call(rbind, manifest), file = file.path(outdir, "permulation_manifest.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ── Summary ──────────────────────────────────────────────────────────────────
tiers <- vapply(pool, function(e) e$tier, integer(1))
dunns <- vapply(pool, function(e) e$dunn_min, numeric(1))
elapsed <- as.numeric(difftime(Sys.time(), start.time, units = "mins"))
log_msg("COMPLETE", sprintf(
  "%d records (Tier1=%d, Tier2=%d) from %d draws | acceptance %.2f%% | median overall Dunn %.3f | %.2f min%s",
  length(pool), sum(tiers == 1L), sum(tiers == 2L), total_draws,
  100 * length(pool) / total_draws, median(dunns), elapsed,
  if (is.finite(lambda_hat)) sprintf(" | lambda = %.4f", lambda_hat) else ""))
if (length(reject_reasons)) {
  tab <- sort(table(reject_reasons), decreasing = TRUE)
  for (i in seq_len(min(3, length(tab)))) {
    log_msg("INFO", sprintf("  rejected %6d x %s", as.integer(tab)[i], names(tab)[i]))
  }
}
