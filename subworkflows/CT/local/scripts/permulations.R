#!/usr/bin/env Rscript
# =============================================================================
# PHYLOPHERE: A Nextflow pipeline for Phenome-Genome studies
# File: subworkflows/CT/local/scripts/permulations.R
# =============================================================================
# Generates null trait permulations for CAAS significance calibration.
#
# Strategies (--perm_strategy):
#   OU : (Default) Fits an Ornstein-Uhlenbeck evolutionary model (alpha, sigsq),
#        simulates traits under OU on the tree, and selects candidate pairs
#        using Global Top 1% PSS (or non-overlapping Bayesian CIs for count data).
#   BM : Brownian Motion simulation on the real tree, rank-matched back onto
#        observed values (simpermvec).
#
# Pool harvesting:
# Every accepted permulation carries EXACTLY N_pairs_obs pairs with Dunn validation:
#   Tier 1 : all N_pairs_obs pairs have mod_dunn >= 1 (fully independent)
#   Tier 2 : exactly one pair falls below mod_dunn 1 (fallback)
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(geiger)
})

log_msg <- function(tag, ...) write(paste0("[", tag, "] ", format(Sys.time()), " ", paste0(...)), stdout())

# Evolutionary model fit + BM/OU covariances + AIC model selection come from the
# vendored phyloq engine (pss_core.R), sourced below once script_dir is known.

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
       "[chunk_size] [include_b0] [pss_top_pct] [max_tries] [pheno_col] ",
       "[n_col] [c_col] [resample_use_n] [trait_type]")
}

arg_or <- function(i, default, cast = as.character) {
  if (length(args) >= i && nzchar(args[i])) cast(args[i]) else default
}

tree.path          <- args[1]
config.file        <- args[2]
number.of.cycles   <- as.integer(args[3])   # size of the harvested pool
selection.strategy <- tolower(args[4])      # "ou" | "bm"
phenotypes         <- args[5]
outdir             <- args[6]
chunk.size         <- arg_or(7,  500L,       as.integer)
include.b0         <- tolower(arg_or(8, "true")) %in% c("1", "true", "t", "yes", "y")
pss_top_pct        <- arg_or(9,  0.01,       as.numeric)
max_tries          <- arg_or(10, 1000000L,   as.integer)
pheno_col_name     <- arg_or(11, "")
n_col              <- arg_or(12, "")   # denominator column (e.g. adult_necropsy_count)
c_col              <- arg_or(13, "")   # numerator column   (e.g. malignant_count)
resample_use_n     <- tolower(arg_or(14, "true")) %in% c("1", "true", "t", "yes", "y")
trait_type         <- tolower(arg_or(15, "auto"))
fop_null           <- tolower(arg_or(16, "false")) %in% c("1", "true", "t", "yes", "y")
max_fop            <- arg_or(17, 100L, as.integer)

if (!selection.strategy %in% c("auto", "best_model", "ou", "bm")) {
  log_msg("WARN", sprintf("Unknown strategy '%s', defaulting to 'auto'", selection.strategy))
  selection.strategy <- "auto"
}

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
  log_msg("INFO", "Created output directory: ", outdir)
}

# ── Locate the lean selector next to this script ─────────────────────────────
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

first_line <- readLines(phenotypes, n = 1L, warn = FALSE)
delim_char <- if (grepl(",", first_line) && !grepl("\t", first_line)) "," else "\t"
first_flds <- strsplit(first_line, delim_char, fixed = TRUE)[[1]]
has_header <- length(first_flds) < 2 || is.na(suppressWarnings(as.numeric(first_flds[2])))

pheno_raw <- read.delim(phenotypes, sep = delim_char, header = has_header,
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

# ── Count data → Jeffreys CIs ────────────────────────────────────────────────
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
    stop("Count columns were supplied but the 'binom' package is unavailable.")
  }
  np <- suppressWarnings(as.numeric(pheno_raw[[n_col]]))
  nc <- suppressWarnings(as.numeric(pheno_raw[[c_col]]))
  good <- is.finite(np) & is.finite(nc) & np > 0 & nc >= 0 & nc <= np
  phenotype.df$n_pop <- np   # denominator, for the pair_n ranking tiebreak
  phenotype.df$ci_lb <- NA_real_; phenotype.df$ci_ub <- NA_real_
  if (any(good)) {
    b <- binom::binom.confint(nc[good], np[good], method = "bayes",
                              priors = c(0.5, 0.5), conf.level = 0.95)
    phenotype.df$ci_lb[good] <- b[, 5]; phenotype.df$ci_ub[good] <- b[, 6]
  }
  log_msg("INFO", sprintf("Count columns '%s'/'%s' found: %d/%d species have Jeffreys CIs",
                          n_col, c_col, sum(good), length(good)))
}

keep <- !is.na(phenotype.df$species) & is.finite(phenotype.df$value)
if (use_ci) keep <- keep & is.finite(phenotype.df$ci_lb) & is.finite(phenotype.df$ci_ub)
phenotype.df <- phenotype.df[keep, , drop = FALSE]
if (nrow(phenotype.df) == 0) stop("No valid phenotype rows remain after NA filtering")

# ── Tree ─────────────────────────────────────────────────────────────────────
tree.o <- read.tree(tree.path)
starting.values <- setNames(phenotype.df$value, phenotype.df$species)

pruned.tree <- drop.tip(tree.o, setdiff(tree.o$tip.label, names(starting.values)))
pruned.tree <- multi2di(pruned.tree, random = FALSE)
pruned.tree$edge.length[pruned.tree$edge.length <= 0] <- 1e-8
starting.values <- starting.values[pruned.tree$tip.label]

# ── Ultrametric check (warn-only) ───────────────────────────────────────────
# KEEP IN SYNC with subworkflows/TRAIT_ANALYSIS/local/src/commons.R. Contrast
# independence (Dunn) and the OU/BM PSS assume a time tree; an ML phylogram lets
# a long terminal branch (rate, not time) inflate a species' contrast-pair
# diameter and wrongly fail otherwise-independent pairs. PhyloPhere expects a
# dated species tree — warn rather than rate-smooth an arbitrary phylogram here.
if (!is.ultrametric(pruned.tree, tol = 1e-6)) {
  log_msg("WARN", "permulation tree is NOT ultrametric (phylogram); contrast ",
          "independence assumes a time tree. Supply a dated species tree.")
}

D <- cophenetic(pruned.tree)

# ── Vendored phyloq engine + lean selector ──────────────────────────────────
pss_core_script <- file.path(script_dir, "pss_core.R")
for (.f in c(lean_script, pss_core_script)) {
  if (!file.exists(.f)) stop(basename(.f), " not found next to permulations.R (looked in '", script_dir, "').")
}
source(pss_core_script)   # fit_models / covariances_from_fits / select_model / calculate_pairwise_scores
source(lean_script)       # evaluate_lean_contrast_selection + shared Dunn core
log_msg("INFO", "Lean contrast filtering + phyloq PSS enabled from ", script_dir)

# ── Evolutionary model: phyloq's fit + AIC selection (verbatim), once, on the
#    observed trait. Covariances are held fixed across all permulation draws.
.force_model <- if (selection.strategy %in% c("bm", "ou")) toupper(selection.strategy) else NULL
obs_fits       <- fit_models(pruned.tree, starting.values)
selected_model <- select_model(obs_fits, force_model = .force_model)
obs_cov        <- covariances_from_fits(pruned.tree, obs_fits)
cov_bm         <- obs_cov$BM
cov_ou         <- obs_cov$OU
log_msg("INFO", sprintf("Evolutionary model: %s (AIC BM = %.3f, OU = %.3f; delta = %.3f)",
                        selected_model, fit_aic(obs_fits$BM), fit_aic(obs_fits$OU),
                        fit_aic(obs_fits$BM) - fit_aic(obs_fits$OU)))

# Simulation tree for the rank-match null: OU-rescaled when OU is selected.
if (selected_model == "OU") {
  simulation_tree <- rescale(pruned.tree, "OU", as.numeric(obs_fits$OU$opt$alpha))
  simulation_tree$edge.length[simulation_tree$edge.length <= 0] <- 1e-8
} else {
  simulation_tree <- pruned.tree
}

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
start.time <- Sys.time()
log_msg("START", sprintf("Harvesting pool (strategy: %s, pool_size: %d, max_tries: %d)",
                         selection.strategy, number.of.cycles, max_tries))

rec_ord    <- order(starting.values)
rec_value  <- unname(starting.values[rec_ord])
n_rec      <- length(rec_value)
rank_order <- seq_len(n_rec)
rec_lb <- rec_ub <- rec_n <- NULL
if (use_ci) {
  ci_map <- setNames(phenotype.df$ci_lb, phenotype.df$species)
  rec_lb <- unname(ci_map[names(starting.values)[rec_ord]])
  ci_map <- setNames(phenotype.df$ci_ub, phenotype.df$species)
  rec_ub <- unname(ci_map[names(starting.values)[rec_ord]])
  n_map <- setNames(phenotype.df$n_pop, phenotype.df$species)
  rec_n <- unname(n_map[names(starting.values)[rec_ord]])
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

    sim_v <- simulatevec(starting.values, simulation_tree)
    sim_ord <- order(sim_v)
    
    pvec <- setNames(rec_value[order(sim_ord)], names(sim_v))

    ci_lb_draw <- NULL; ci_ub_draw <- NULL; n_draw <- NULL
    if (use_ci) {
      ci_lb_draw <- setNames(rec_lb[order(sim_ord)], names(sim_v))
      ci_ub_draw <- setNames(rec_ub[order(sim_ord)], names(sim_v))
      n_draw     <- setNames(rec_n[order(sim_ord)], names(sim_v))
    }

    e <- tryCatch(
      evaluate_lean_contrast_selection(
        trait_vec = pvec, D = D, target_pairs = target_pairs,
        tree = pruned.tree, cov_bm = cov_bm, cov_ou = cov_ou,
        selected_model = selected_model,
        ci_lb = ci_lb_draw, ci_ub = ci_ub_draw,
        top_pct = pss_top_pct, n_vec = n_draw,
        ordinal = if (trait_type == "ordinal") TRUE
                  else if (trait_type == "continuous") FALSE
                  else NULL
      ),
      error = function(err) list(tier = 0L, n_pairs = 0L, dunn_min = 0,
                                 n_below = NA_integer_, fg = NULL, bg = NULL,
                                 reason = paste0("error: ", conditionMessage(err)))
    )

    # Retain the permuted vector (+ CI/n draws) on accepted cycles so the FOP
    # mirror can harvest alternative hypotheses around this exact labeling after
    # the pool is assembled.
    if (fop_null && e$tier %in% c(1L, 2L)) {
      e$pvec <- pvec
      e$ci_lb_draw <- ci_lb_draw; e$ci_ub_draw <- ci_ub_draw; e$n_draw <- n_draw
    }

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
    "Tier 1 exhausted after %d draws: topping up with Tier 2. %d Tier-1 + %d Tier-2 = %d/%d records",
    n1, use2, length(pool), number.of.cycles))
} else {
  tab <- sort(table(reject_reasons), decreasing = TRUE)
  top <- seq_len(min(3, length(tab)))
  stop(sprintf(
    paste0("Permulation pool could not be filled: %d Tier-1 + %d Tier-2 = %d of the requested %d ",
           "after %d draws and %d escalation(s) of max_tries (final budget %d).\n",
           "  Top rejection reasons: %s"),
    n1, n2, n1 + n2, number.of.cycles, total_draws, escalations, budget,
    if (length(tab)) paste(sprintf("%s (%d)", names(tab)[top], as.integer(tab)[top]), collapse = "; ") else "none recorded"))
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
  fmt <- function(x) paste(format(x, digits = 15, trim = TRUE), collapse = ",")
  manifest[[b]] <- data.frame(cycle = paste0("b_", b), tier = e$tier,
                              n_pairs = e$n_pairs, dunn_min = e$dunn_min,
                              n_below = e$n_below, mode = e$mode,
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

# ── FOP mirror: per-cycle alternative-hypothesis harvest ─────────────────────
# Mirrors the observed FOP harvest (selection_algorithm.R::fop_pair_sel.f) for
# every accepted permulation cycle, so the null holds the SAME domain-pooled
# statistic scoring_compute.R §2b builds on the observed data. H1 is the cycle's
# already-accepted canonical contrast; H2..Hn are Dunn-independent alternatives
# drawn from the same Voronoi domains, ranked by min-PSS, capped at max_fop.
#   resample_fop.tab         : "<cycle>~H<m>" \t fg_csv \t bg_csv   (fanned discovery input)
#   resample_fop_pairs.tsv   : cycle, hypothesis_id, pair(domain), species1, species2, pss_score
if (fop_null) {
  log_msg("START", sprintf("FOP mirror harvest for %d accepted cycles (max_fop=%d)",
                           length(pool), max_fop))
  fop_rows  <- list()
  pair_rows <- list()
  n_hyp_tot <- 0L
  for (b in seq_along(pool)) {
    e <- pool[[b]]
    cyc <- paste0("b_", b)
    if (is.null(e$pvec) || is.null(e$fg) || is.null(e$bg)) next
    hv <- tryCatch(
      lean_fop_harvest(
        trait_vec = e$pvec, D = D, target_pairs = target_pairs,
        tree = pruned.tree, cov_bm = cov_bm, cov_ou = cov_ou,
        selected_model = selected_model,
        ci_lb = e$ci_lb_draw, ci_ub = e$ci_ub_draw, n_vec = e$n_draw,
        top_pct = pss_top_pct, max_fop = max_fop,
        ordinal = if (trait_type == "ordinal") TRUE
                  else if (trait_type == "continuous") FALSE else NULL,
        canon_pairs = data.frame(species1 = e$fg, species2 = e$bg,
                                 stringsAsFactors = FALSE)),
      error = function(err) { log_msg("WARN", sprintf("FOP harvest %s: %s", cyc, conditionMessage(err))); NULL })
    if (is.null(hv) || length(hv$hypotheses) == 0L) {
      # fall back to H1-only so the cycle still enters the fanned discovery
      fop_rows[[length(fop_rows) + 1L]] <- data.frame(
        cycle = paste0(cyc, "~H1"),
        fg = paste(e$fg, collapse = ","), bg = paste(e$bg, collapse = ","),
        stringsAsFactors = FALSE)
      next
    }
    for (h_id in names(hv$hypotheses)) {
      hd <- hv$hypotheses[[h_id]]
      fop_rows[[length(fop_rows) + 1L]] <- data.frame(
        cycle = paste0(cyc, "~", h_id),
        fg = paste(hd$species1, collapse = ","),
        bg = paste(hd$species2, collapse = ","),
        stringsAsFactors = FALSE)
      pair_rows[[length(pair_rows) + 1L]] <- data.frame(
        cycle = cyc, hypothesis_id = h_id,
        pair = if ("cluster" %in% names(hd)) as.integer(hd$cluster) else seq_len(nrow(hd)),
        species1 = hd$species1, species2 = hd$species2,
        pss_score = if ("pss_score" %in% names(hd)) hd$pss_score else NA_real_,
        stringsAsFactors = FALSE)
    }
    n_hyp_tot <- n_hyp_tot + length(hv$hypotheses)
  }
  if (length(fop_rows) > 0) {
    write.table(do.call(rbind, fop_rows), file = file.path(outdir, "resample_fop.tab"),
                sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
  if (length(pair_rows) > 0) {
    write.table(do.call(rbind, pair_rows), file = file.path(outdir, "resample_fop_pairs.tsv"),
                sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
  log_msg("COMPLETE", sprintf("FOP mirror: %d cycles -> %d hypothesis labelings (mean %.1f/cycle) -> resample_fop.tab",
                              length(pool), n_hyp_tot,
                              if (length(pool)) n_hyp_tot / length(pool) else 0))
}

# ── Summary ──────────────────────────────────────────────────────────────────
tiers <- vapply(pool, function(e) e$tier, integer(1))
dunns <- vapply(pool, function(e) e$dunn_min, numeric(1))
elapsed <- as.numeric(difftime(Sys.time(), start.time, units = "mins"))
log_msg("COMPLETE", sprintf(
  "%d records (Tier1=%d, Tier2=%d) from %d draws | acceptance %.2f%% | median overall Dunn %.3f | %.2f min",
  length(pool), sum(tiers == 1L), sum(tiers == 2L), total_draws,
  100 * length(pool) / total_draws, median(dunns), elapsed))
