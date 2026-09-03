#!/usr/bin/env Rscript
# =============================================================================
# PROBE (not a pipeline step): size the two open unknowns for the FOP
# multi-hypothesis permulation null.
#
#   Unknown 1  precompute grid width  = # distinct candidate pairs the fanned
#              discovery would ever touch across the permulation pool.
#   Unknown 3  per-cycle cost         = median wall time of one null cycle,
#              split PSS vs FOP-harvest vs rest, plus the draws-per-accept ratio.
#
# It replays the null's candidate-selection + a lean FOP harvest on real inputs,
# WITHOUT scoring any ASR path (Unknown 2 already verified: the ASR walk is
# labeling-invariant, so grid width and harvest cost are the whole question).
#
# Mirrors, deliberately and minimally:
#   * lean_contrast_selector.R::evaluate_lean_contrast_selection  (candidate gate
#     + rank), stages reproduced inline so cand_df is observable per draw.
#   * selection_algorithm.R::fop_pair_sel.f  steps 2-4 (Voronoi partition +
#     capped combinatorial harvest + Dunn>=1 filter), reproduced inline lean.
# The shared Dunn core (mod_dunn_lean / overall_dunn_lean / greedy_dunn_select /
# rank_candidates) and the phyloq PSS engine are sourced verbatim, never copied.
#
# Usage:
#   Rscript .probe_fop_null.R \
#     --tree      pruned_tree_file.nwk \
#     --trait     boot_traitfile.tab \
#     --config    traitfile.tab            # observed FG/BG/pair-id config (V3 = pair count)
#     --draws     1000 \
#     --pheno-col malignant_prevalence \
#     --n-col     adult_necropsy_count \
#     --c-col     malignant_count \
#     --strategy  best_model \
#     --pss-top-pct 0.05 \
#     --trait-type  auto \
#     --max-fop   100 \
#     --seed      1998 \
#     --out       probe_fop_null.rds
#
# All flags except --tree/--trait/--config have defaults. When --n-col/--c-col
# resolve, the probe reports the candidate count under BOTH the CI gate (what
# this trait's null actually uses) and the pure PSS top-pct gate, so one run
# covers both regimes.
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(geiger)   # ratematrix / sim.char, for simulatevec
})

# verbatim from permulations.R:32-36 (BM sim on the OU-rescaled tree)
simulatevec <- function(namedvec, treewithbranchlengths) {
  rm   <- ratematrix(treewithbranchlengths, namedvec)
  sims <- sim.char(treewithbranchlengths, rm, nsim = 1)
  setNames(as.data.frame(sims)[, 1], rownames(sims))
}

# ── arg parsing ──────────────────────────────────────────────────────────────
.args <- commandArgs(trailingOnly = TRUE)
.getopt <- function(flag, default = NULL) {
  i <- match(flag, .args)
  if (is.na(i) || i == length(.args)) return(default)
  .args[i + 1L]
}
opt <- list(
  tree        = .getopt("--tree"),
  trait       = .getopt("--trait"),
  config      = .getopt("--config"),
  draws       = as.integer(.getopt("--draws", "1000")),
  pheno_col   = .getopt("--pheno-col", ""),
  n_col       = .getopt("--n-col", ""),
  c_col       = .getopt("--c-col", ""),
  strategy    = tolower(.getopt("--strategy", "best_model")),
  pss_top_pct = as.numeric(.getopt("--pss-top-pct", "0.05")),
  trait_type  = tolower(.getopt("--trait-type", "auto")),
  max_fop     = as.integer(.getopt("--max-fop", "100")),
  seed        = as.integer(.getopt("--seed", "1998")),
  out         = .getopt("--out", "probe_fop_null.rds")
)
stopifnot(!is.null(opt$tree), !is.null(opt$trait), !is.null(opt$config))

msg <- function(...) cat(sprintf("[PROBE %s] ", format(Sys.time(), "%H:%M:%S")),
                         sprintf(...), "\n", sep = "")

# ── locate + source the shared cores (verbatim) ──────────────────────────────
script_dir <- {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd()
}
for (f in c("pss_core.R", "lean_contrast_selector.R")) {
  p <- file.path(script_dir, f)
  if (!file.exists(p)) stop(f, " not found next to the probe (", script_dir, ")")
  source(p)
}
msg("sourced pss_core.R + lean_contrast_selector.R from %s", script_dir)

# ── load inputs, mirroring permulations.R:94-201 ─────────────────────────────
cfg <- read.table(opt$config, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
if (ncol(cfg) < 3) stop("config needs a V3 pair-id column")
target_pairs <- suppressWarnings(max(as.integer(cfg$V3), na.rm = TRUE))
if (!is.finite(target_pairs) || target_pairs <= 0L) stop("bad target_pairs from config V3")
msg("target_pairs (N_pairs_obs) = %d", target_pairs)

first_line <- readLines(opt$trait, n = 1L, warn = FALSE)
delim <- if (grepl(",", first_line) && !grepl("\t", first_line)) "," else "\t"
flds  <- strsplit(first_line, delim, fixed = TRUE)[[1]]
has_header <- length(flds) < 2 || is.na(suppressWarnings(as.numeric(flds[2])))
pheno_raw <- read.delim(opt$trait, sep = delim, header = has_header,
                        stringsAsFactors = FALSE,
                        na.strings = c("", "NA", "NaN", "nan", "NULL", "null"))
sp_col <- if ("species" %in% names(pheno_raw)) "species" else names(pheno_raw)[1]
val_col <- if (nzchar(opt$pheno_col) && opt$pheno_col %in% names(pheno_raw)) {
  opt$pheno_col
} else {
  cand <- setdiff(names(pheno_raw)[vapply(pheno_raw, is.numeric, logical(1))], sp_col)
  if (length(cand)) cand[1] else names(pheno_raw)[2]
}
msg("species col '%s', value col '%s'", sp_col, val_col)

phenotype.df <- data.frame(
  species = as.character(pheno_raw[[sp_col]]),
  value   = suppressWarnings(as.numeric(pheno_raw[[val_col]])),
  stringsAsFactors = FALSE
)

# Jeffreys CIs, mirroring permulations.R:142-174 (auto-detect + explicit override)
n_col <- opt$n_col; c_col <- opt$c_col
if (!nzchar(n_col) || !n_col %in% names(pheno_raw)) {
  hit <- intersect(c("n_trait","n","n_population","sample_size","total_count","N",
                     "adult_necropsy_count","necropsy_count"), names(pheno_raw))
  if (length(hit)) n_col <- hit[1]
}
if (!nzchar(c_col) || !c_col %in% names(pheno_raw)) {
  hit <- intersect(c("c_trait","c","c_cases","n_cases","cases","C",
                     "malignant_count","neoplasia_count"), names(pheno_raw))
  if (length(hit)) c_col <- hit[1]
}
use_ci <- nzchar(n_col) && nzchar(c_col) &&
          n_col %in% names(pheno_raw) && c_col %in% names(pheno_raw) &&
          requireNamespace("binom", quietly = TRUE)
if (use_ci) {
  np <- suppressWarnings(as.numeric(pheno_raw[[n_col]]))
  nc <- suppressWarnings(as.numeric(pheno_raw[[c_col]]))
  good <- is.finite(np) & is.finite(nc) & np > 0 & nc >= 0 & nc <= np
  phenotype.df$n_pop <- np
  phenotype.df$ci_lb <- NA_real_; phenotype.df$ci_ub <- NA_real_
  if (any(good)) {
    b <- binom::binom.confint(nc[good], np[good], method = "bayes",
                              priors = c(0.5, 0.5), conf.level = 0.95)
    phenotype.df$ci_lb[good] <- b[, 5]; phenotype.df$ci_ub[good] <- b[, 6]
  }
  msg("CI gate ON: cols '%s'/'%s', %d/%d species with Jeffreys CIs",
      n_col, c_col, sum(good), length(good))
} else {
  msg("CI gate OFF (no usable count columns) — PSS top-pct gate only")
}

keep <- !is.na(phenotype.df$species) & is.finite(phenotype.df$value)
if (use_ci) keep <- keep & is.finite(phenotype.df$ci_lb) & is.finite(phenotype.df$ci_ub)
phenotype.df <- phenotype.df[keep, , drop = FALSE]

tree.o <- read.tree(opt$tree)
starting.values <- setNames(phenotype.df$value, phenotype.df$species)
pruned.tree <- drop.tip(tree.o, setdiff(tree.o$tip.label, names(starting.values)))
pruned.tree <- multi2di(pruned.tree, random = FALSE)
pruned.tree$edge.length[pruned.tree$edge.length <= 0] <- 1e-8
starting.values <- starting.values[pruned.tree$tip.label]
phenotype.df <- phenotype.df[match(pruned.tree$tip.label, phenotype.df$species), ]
n_sp <- length(starting.values)
msg("analysis species after pruning = %d  (n_choose_2 = %d)", n_sp, choose(n_sp, 2))
if (!is.ultrametric(pruned.tree, tol = 1e-6))
  msg("WARN tree is not ultrametric (phylogram) — Dunn assumes a time tree")

D <- cophenetic(pruned.tree)

# ── observed evolutionary model (phyloq, verbatim), fixed across draws ───────
.force <- if (opt$strategy %in% c("bm", "ou")) toupper(opt$strategy) else NULL
obs_fits       <- fit_models(pruned.tree, starting.values)
selected_model <- select_model(obs_fits, force_model = .force)
obs_cov        <- covariances_from_fits(pruned.tree, obs_fits)
cov_bm <- obs_cov$BM; cov_ou <- obs_cov$OU
msg("model = %s (AIC BM=%.2f OU=%.2f)", selected_model,
    fit_aic(obs_fits$BM), fit_aic(obs_fits$OU))

simulation_tree <- if (selected_model == "OU") {
  st <- rescale(pruned.tree, "OU", as.numeric(obs_fits$OU$opt$alpha))
  st$edge.length[st$edge.length <= 0] <- 1e-8; st
} else pruned.tree

ordinal <- (
  if (opt$trait_type == "ordinal") TRUE
  else if (opt$trait_type == "continuous") FALSE
  else .is_ordinal_vec(starting.values)
)

rec_ord   <- order(starting.values)
rec_value <- unname(starting.values[rec_ord])
rec_lb <- rec_ub <- rec_n <- NULL
if (use_ci) {
  m <- setNames(phenotype.df$ci_lb, phenotype.df$species); rec_lb <- unname(m[names(starting.values)[rec_ord]])
  m <- setNames(phenotype.df$ci_ub, phenotype.df$species); rec_ub <- unname(m[names(starting.values)[rec_ord]])
  m <- setNames(phenotype.df$n_pop, phenotype.df$species); rec_n  <- unname(m[names(starting.values)[rec_ord]])
}

# ── candidate stage: mirror evaluate_lean_contrast_selection stages 1-2 ──────
# Returns the ranked candidate data.frame (cand_df) the FOP harvest draws from,
# plus the mode actually taken. `gate` forces "ci" | "pss" | "auto".
candidate_df <- function(trait_vec, ci_lb = NULL, ci_ub = NULL, n_vec = NULL,
                         gate = "auto") {
  sp <- intersect(names(trait_vec), rownames(D))
  trait_vec <- trait_vec[sp]
  this_ci <- gate != "pss" && !is.null(ci_lb) && !is.null(ci_ub)
  if (this_ci) { lb <- ci_lb[sp]; ub <- ci_ub[sp] }
  this_ord <- gate == "auto" && !this_ci && isTRUE(ordinal)

  tr <- if (length(sp) < length(pruned.tree$tip.label))
          ape::drop.tip(pruned.tree, setdiff(pruned.tree$tip.label, sp)) else pruned.tree
  tsp <- tr$tip.label
  sc <- calculate_pairwise_scores(trait_vec[tsp], tr,
                                  cov_bm[tsp, tsp, drop = FALSE],
                                  cov_ou[tsp, tsp, drop = FALSE], selected_model)
  hi1   <- sc$TraitValue1 >= sc$TraitValue2
  c_hi  <- ifelse(hi1, sc$Species1, sc$Species2)
  c_lo  <- ifelse(hi1, sc$Species2, sc$Species1)
  c_dif <- abs(sc$TraitValue1 - sc$TraitValue2)
  c_dst <- sc$PatristicDistance
  pss   <- sc$FinalScore
  d0 <- c_dif > 0
  c_hi<-c_hi[d0]; c_lo<-c_lo[d0]; c_dif<-c_dif[d0]; c_dst<-c_dst[d0]; pss<-pss[d0]
  if (!length(c_hi)) return(list(df = NULL, mode = "none"))

  if (this_ci) {
    tmask <- lb[c_hi] > ub[c_lo]
  } else if (this_ord) {
    hi <- max(trait_vec, na.rm = TRUE); lo <- min(trait_vec, na.rm = TRUE)
    tmask <- trait_vec[c_hi] >= hi & trait_vec[c_lo] <= lo
  } else {
    tmask <- rep(TRUE, length(c_hi))
  }
  s1 <- which(tmask)
  if (!length(s1)) return(list(df = NULL, mode = if (this_ci) "ci" else if (this_ord) "ordinal" else "pss"))

  if (!this_ci && !this_ord) {
    nkeep  <- min(max(1L, ceiling(length(s1) * opt$pss_top_pct)), length(s1))
    keep_i <- s1[order(pss[s1], decreasing = TRUE)][seq_len(nkeep)]
  } else keep_i <- s1
  keep <- logical(length(c_hi)); keep[keep_i] <- TRUE

  df <- data.frame(species1 = c_hi[keep], species2 = c_lo[keep],
                   distance = c_dst[keep], abs_diff = c_dif[keep],
                   pss_score = pss[keep], stringsAsFactors = FALSE)
  if (!is.null(n_vec)) df$pair_n <- as.numeric(n_vec[df$species1]) + as.numeric(n_vec[df$species2])
  list(df = rank_candidates(df),
       mode = if (this_ci) "ci" else if (this_ord) "ordinal" else "pss")
}

# ── lean FOP harvest: mirror fop_pair_sel.f steps 2-4 ────────────────────────
# canon = greedy_dunn_select(enforce_dunn=TRUE) on cand_df (== pair_sel.f's H1);
# Voronoi-partition species by nearest canonical pair; capped combinatorial draw
# of one in-domain candidate per domain; keep hypotheses with overall Dunn>=1.
# Returns the K per-domain pair sets of every kept hypothesis (for the union).
fop_harvest_lean <- function(cand_df, max_fop = opt$max_fop) {
  if (is.null(cand_df) || nrow(cand_df) < 1) return(list(hyps = list(), K = 0L, n_raw = 0L))
  canon <- greedy_dunn_select(cand_df, D, target = Inf, enforce_dunn = TRUE)
  K <- length(canon$members)
  if (K < 1) return(list(hyps = list(), K = 0L, n_raw = 0L))
  cm <- canon$members
  all_sp <- rownames(D)
  dom <- vapply(all_sp, function(s)
    which.min(vapply(seq_len(K), function(k) min(D[s, cm[[k]][1]], D[s, cm[[k]][2]]), numeric(1))),
    integer(1))
  names(dom) <- all_sp
  pools <- lapply(seq_len(K), function(k) {
    ds <- names(dom)[dom == k]
    cand_df[cand_df$species1 %in% ds & cand_df$species2 %in% ds, , drop = FALSE]
  })
  psz <- vapply(pools, nrow, integer(1))
  if (any(psz == 0L)) return(list(hyps = list(), K = K, n_raw = 0L, degenerate = TRUE))

  ITER_CAP <- as.integer(max_fop) * 20L
  total <- prod(as.numeric(psz))
  idx <- if (total <= ITER_CAP)
    do.call(expand.grid, c(lapply(psz, seq_len), list(KEEP.OUT.ATTRS = FALSE)))
  else
    as.data.frame(lapply(psz, function(n) sample.int(n, ITER_CAP, replace = TRUE)))

  seen <- character(0); hyps <- list(); minpss <- numeric(0)
  for (it in seq_len(nrow(idx))) {
    ci <- as.integer(unlist(idx[it, , drop = FALSE], use.names = FALSE))
    rows <- do.call(rbind, lapply(seq_len(K), function(k) pools[[k]][ci[k], ]))
    sp <- c(rows$species1, rows$species2)
    if (length(unique(sp)) < 2L * K) next
    sig <- paste(sort(sp), collapse = "|")
    if (sig %in% seen) next
    seen <- c(seen, sig)
    mem <- lapply(seq_len(K), function(i) c(rows$species1[i], rows$species2[i]))
    if (overall_dunn_lean(D, mem) >= 1.0) {
      hyps[[length(hyps) + 1L]] <- mem
      minpss <- c(minpss, suppressWarnings(min(rows$pss_score, na.rm = TRUE)))
    }
  }
  n_raw <- length(hyps)
  # Mirror fop_pair_sel.f:288-307: rank Dunn-valid hypotheses by min-PSS desc,
  # keep only the top (max_fop - 1). This is the set the null would actually score.
  if (n_raw > 0L) {
    ord  <- order(-replace(minpss, is.na(minpss), -Inf))
    keep <- head(ord, max(0L, as.integer(max_fop) - 1L))
    hyps <- hyps[keep]
  }
  list(hyps = hyps, K = K, n_raw = n_raw, n_kept = length(hyps))
}

# ── the loop ────────────────────────────────────────────────────────────────
set.seed(opt$seed)
N <- opt$draws
msg("running %d draws (use_ci=%s, ordinal=%s, pss_top_pct=%.3f)",
    N, use_ci, isTRUE(ordinal), opt$pss_top_pct)

pair_sig <- function(p) paste(sort(p), collapse = "~")

rec <- vector("list", N)
union_ci  <- new.env(parent = emptyenv())   # distinct candidate pairs (CI gate)
union_pss <- new.env(parent = emptyenv())   # distinct candidate pairs (PSS gate)
union_hrv <- new.env(parent = emptyenv())   # distinct pairs actually in a kept hypothesis
union_curve <- integer(N)                    # cumulative |union_hrv| after each draw
n_accept <- 0L
t0 <- Sys.time()

for (b in seq_len(N)) {
  sim_v   <- simulatevec(starting.values, simulation_tree)
  sim_ord <- order(sim_v)
  pvec    <- setNames(rec_value[order(sim_ord)], names(sim_v))
  lbd <- ubd <- nd <- NULL
  if (use_ci) {
    lbd <- setNames(rec_lb[order(sim_ord)], names(sim_v))
    ubd <- setNames(rec_ub[order(sim_ord)], names(sim_v))
    nd  <- setNames(rec_n[order(sim_ord)],  names(sim_v))
  }

  # candidate stage — timed as "PSS+gate" (calculate_pairwise_scores dominates)
  t_pss <- system.time({
    prim <- candidate_df(pvec, lbd, ubd, nd, gate = "auto")
  })["elapsed"]
  # secondary gate for the cross-regime number (only when the primary was CI)
  alt_n <- NA_integer_
  if (identical(prim$mode, "ci")) {
    alt <- candidate_df(pvec, NULL, NULL, nd, gate = "pss")
    alt_n <- if (is.null(alt$df)) 0L else nrow(alt$df)
    if (!is.null(alt$df)) for (i in seq_len(nrow(alt$df)))
      assign(pair_sig(c(alt$df$species1[i], alt$df$species2[i])), TRUE, union_pss)
  }

  cand <- prim$df
  n_cand <- if (is.null(cand)) 0L else nrow(cand)
  if (!is.null(cand)) {
    tgt <- if (identical(prim$mode, "ci")) union_ci else union_pss
    for (i in seq_len(nrow(cand)))
      assign(pair_sig(c(cand$species1[i], cand$species2[i])), TRUE, tgt)
  }

  # FOP harvest — timed separately
  t_hrv <- system.time({ H <- fop_harvest_lean(cand) })["elapsed"]
  for (mem in H$hyps) for (p in mem)
    assign(pair_sig(p), TRUE, union_hrv)

  # tier: reuse the real grader for the accept ratio
  e <- tryCatch(
    evaluate_lean_contrast_selection(
      trait_vec = pvec, D = D, target_pairs = target_pairs,
      tree = pruned.tree, cov_bm = cov_bm, cov_ou = cov_ou,
      selected_model = selected_model,
      ci_lb = lbd, ci_ub = ubd, top_pct = opt$pss_top_pct, n_vec = nd,
      ordinal = if (opt$trait_type == "ordinal") TRUE
                else if (opt$trait_type == "continuous") FALSE else NULL),
    error = function(err) list(tier = 0L, reason = conditionMessage(err)))
  if (isTRUE(e$tier == 1L)) n_accept <- n_accept + 1L

  union_curve[b] <- length(ls(union_hrv))
  rec[[b]] <- data.frame(
    draw = b, mode = prim$mode, n_cand = n_cand, n_cand_pss_alt = alt_n,
    fop_K = H$K, fop_hyps = H$n_raw, fop_kept = if (is.null(H$n_kept)) 0L else H$n_kept,
    tier = if (is.null(e$tier)) NA_integer_ else e$tier,
    t_pss = as.numeric(t_pss), t_harvest = as.numeric(t_hrv),
    union_pairs_running = union_curve[b], stringsAsFactors = FALSE)

  if (b %% 50L == 0L || b == N) {
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    msg("draw %d/%d | %.2f s/draw | accepts=%d | union(harvest pairs)=%d",
        b, N, el / b, n_accept, union_curve[b])
  }
}

df <- do.call(rbind, rec)
summ <- list(
  inputs = list(tree = opt$tree, trait = opt$trait, config = opt$config,
                n_species = n_sp, n_choose_2 = choose(n_sp, 2),
                target_pairs = target_pairs, model = selected_model,
                use_ci = use_ci, ordinal = isTRUE(ordinal),
                pss_top_pct = opt$pss_top_pct, draws = N),
  timing = list(
    median_t_pss     = median(df$t_pss),
    median_t_harvest = median(df$t_harvest),
    median_t_cycle   = median(df$t_pss + df$t_harvest),
    p90_t_cycle      = quantile(df$t_pss + df$t_harvest, 0.9, names = FALSE),
    total_wall_s     = as.numeric(difftime(Sys.time(), t0, units = "secs"))),
  acceptance = list(
    tier1_accepts = n_accept, draws = N,
    draws_per_accept = if (n_accept > 0) N / n_accept else NA_real_),
  grid_width = list(
    distinct_candidate_pairs_primary = length(ls(if (use_ci) union_ci else union_pss)),
    distinct_candidate_pairs_pss     = length(ls(union_pss)),
    distinct_harvest_pairs           = length(ls(union_hrv)),
    union_curve                      = union_curve,
    saturation_frac_at_end = {
      uc <- union_curve
      if (length(uc) > 20 && tail(uc, 1) > 0)
        1 - (tail(uc, 1) - uc[max(1, length(uc) - 20)]) / tail(uc, 1)
      else NA_real_
    }),
  per_draw = df)

saveRDS(summ, opt$out)

cat("\n================ PROBE RESULT ================\n")
cat(sprintf("species=%d  n_choose_2=%d  target_pairs=%d  model=%s  gate=%s\n",
            n_sp, choose(n_sp, 2), target_pairs, selected_model,
            if (use_ci) "CI (+PSS alt)" else if (isTRUE(ordinal)) "ordinal" else "PSS"))
cat(sprintf("\nUNKNOWN 1 — grid width\n"))
cat(sprintf("  distinct candidate pairs (primary gate) : %d\n", summ$grid_width$distinct_candidate_pairs_primary))
cat(sprintf("  distinct candidate pairs (PSS %.0f%% gate) : %d\n",
            100 * opt$pss_top_pct, summ$grid_width$distinct_candidate_pairs_pss))
cat(sprintf("  distinct pairs ever in a kept hypothesis : %d   <-- true precompute pair axis\n",
            summ$grid_width$distinct_harvest_pairs))
cat(sprintf("  union curve: %s ... %s  (last 20 draws added %d)\n",
            paste(head(union_curve, 3), collapse = ","),
            paste(tail(union_curve, 3), collapse = ","),
            tail(union_curve, 1) - union_curve[max(1, N - 20)]))
cat(sprintf("\nUNKNOWN 3 — per-cycle cost\n"))
cat(sprintf("  median PSS+gate   : %.3f s\n", summ$timing$median_t_pss))
cat(sprintf("  median FOP harvest: %.3f s\n", summ$timing$median_t_harvest))
cat(sprintf("  median cycle      : %.3f s   (p90 %.3f s)\n",
            summ$timing$median_t_cycle, summ$timing$p90_t_cycle))
cat(sprintf("  draws per tier-1 accept: %.1f   (%d accepts / %d draws)\n",
            summ$acceptance$draws_per_accept, n_accept, N))
cat(sprintf("\nsaved -> %s\n", opt$out))
cat("=============================================\n")
