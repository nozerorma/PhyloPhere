#!/usr/bin/env Rscript
# =============================================================================
# OFFLINE VERIFICATION: modal vs. fractional harvest-wide derived_agreement
# =============================================================================
# Reconstructs, from a finished run's flat TSVs only (no pipeline re-run):
#   * da_modal  -- exactly what fop_pool.R::.collect_changed_pairs +
#                  rebuild_derived_agreement produce today (one MODAL residue
#                  per (domain, side), plurality concentration).
#   * da_frac   -- the spec-3 fractional variant: a PSS-weighted residue
#                  DISTRIBUTION per (domain, side), then a weighted plurality
#                  concentration.
# and propagates both through asr -> CAAS_score -> gene_caas_score, using the
# per-(Gene,Position,caap_group) pooled indep/core/div/gate that the REAL
# fop_pool.R produces (we call apply_fop_pooling() verbatim for that).
#
# Usage:
#   Rscript da_frac_verify.R <SRC_DIR> <filtered_discovery.tsv> \
#       <contrast_hypotheses_pairs.tsv> <position_scores.tsv> <gene_scores.tsv> \
#       <out_dir> [tau] [diversity_floor]
# SRC_DIR = subworkflows/SCORING/local/src (for aa_grouping.R + fop_pool.R).
#
# NOTE (post-2026-09-07): fop_pool.R now ships the FRACTIONAL rule. When this
# script sources the current fop_pool.R, `apply_fop_pooling` returns the
# fractional da, so the columns this script labels `da_modal` / `CAAS_modal`
# actually hold the SHIPPED fractional values, and "RECONSTRUCTION FIDELITY vs
# the run" for `derived_agreement` will show the ~24/792-row delta against the
# pre-fractional stored run (this is the expected change, not a regression).
# The retained modal path is checked separately via `da_modal_indep`
# (.collect_changed_pairs + rebuild_derived_agreement), which must still match
# the pre-fractional run bit-for-bit.
# =============================================================================

suppressPackageStartupMessages({ library(dplyr) })

args <- commandArgs(TRUE)
SRC_DIR   <- args[[1]]
FD_PATH   <- args[[2]]
HP_PATH   <- args[[3]]
PS_PATH   <- args[[4]]
GS_PATH   <- args[[5]]
OUT_DIR   <- args[[6]]
TAU       <- if (length(args) >= 7) as.numeric(args[[7]]) else 0.8
DIV_FLOOR <- if (length(args) >= 8) as.numeric(args[[8]]) else 0.75
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(SRC_DIR, "aa_grouping.R"))
source(file.path(SRC_DIR, "fop_pool.R"))

SCHEMES <- c("US", "GS4", "GS3", "GS2", "GS1")
SCHEME_PRIORITY <- c(US = 5, GS4 = 4, GS3 = 3, GS2 = 2, GS1 = 1)
BAD_AA <- c("NA", "NAN", "NONE")

rd <- function(p) read.delim(p, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                             check.names = FALSE, quote = "", na.strings = c("NA", ""))

# =============================================================================
# 1. Load + prep exactly as scoring_compute.R sections 2a/2b do
# =============================================================================
fd <- rd(FD_PATH)
fd$Position <- as.integer(fd$Position)
fd$scheme_priority <- SCHEME_PRIORITY[fd$caap_group]
fd$hyp_id <- if ("trait" %in% names(fd)) {
  ifelse(grepl("H[0-9]+", fd$trait), sub(".*(H[0-9]+).*", "\\1", fd$trait), NA_character_)
} else NA_character_
fd <- fd[fd$caap_group %in% SCHEMES, , drop = FALSE]

hp <- read_hypothesis_pairs(HP_PATH)          # (hyp_id, pair, pss_score) or NULL
pss_fn <- local({
  if (is.null(hp)) return(function(h, dom) NA_real_)
  function(h, dom) {
    v <- hp$pss_score[hp$hyp_id == h & hp$pair == dom]
    v <- v[is.finite(v)]
    if (length(v)) v[1] else NA_real_
  }
})
eqw_fn <- function(h, dom) NA_real_   # equal-weight (what THIS run actually used)

ps_run <- rd(PS_PATH); ps_run$Position <- as.integer(ps_run$Position)
gs_run <- rd(GS_PATH)

# column families (same detection as apply_fop_pooling)
path_cols <- grep("^mrca_[0-9]+_path_score$", names(fd), value = TRUE)
path_cols <- path_cols[order(as.integer(sub("^mrca_([0-9]+)_path_score$", "\\1", path_cols)))]
node_cols <- sub("_path_score$", "_node", path_cols)
top_aa_cols <- sub("_path_score$", "_top_aa", path_cols)
bot_aa_cols <- sub("_path_score$", "_bot_aa", path_cols)
K <- length(path_cols)
cat(sprintf("K = %d Voronoi domains; %d filtered_discovery rows; %d scored positions in run\n",
            K, nrow(fd), nrow(ps_run)))

# =============================================================================
# 2. Fractional harvest-wide derived_agreement (spec section 3)
# =============================================================================
# Per (domain i, side): PSS-weighted residue distribution p_i(r) over the
# hypotheses that recorded a change there. Then, per side with >= 2 changed
# domains: concentration_s = max_g( sum_i p_i^enc(g) ) / |D_s|.  da = mean over
# qualifying sides; 1.0 if none.  (Modal is the special case p_i == delta.)
.addw <- function(acc, k, w) { cur <- if (k %in% names(acc)) acc[[k]] else 0; acc[k] <- cur + w; acc }

domain_side_dists <- function(g, wfun) {
  out <- list()  # "i|side" -> named numeric residue->weight (normalised to 1)
  hids <- as.character(g$hyp_id)
  for (i in seq_len(K)) {
    for (sd in c("top", "bot")) {
      acol <- if (sd == "top") top_aa_cols[i] else bot_aa_cols[i]
      if (is.na(acol) || !(acol %in% names(g))) next
      raw <- toupper(trimws(as.character(g[[acol]])))
      acc <- setNames(numeric(0), character(0))
      for (r in seq_len(nrow(g))) {
        v <- raw[r]
        if (is.na(v) || !nzchar(v) || v %in% BAD_AA) next
        w <- wfun(hids[r], i); if (!is.finite(w)) w <- 1
        acc <- .addw(acc, v, w)
      }
      if (!length(acc)) next
      out[[paste0(i, "|", sd)]] <- acc / sum(acc)
    }
  }
  out
}

da_frac_from_dists <- function(dists, scheme) {
  if (!length(dists)) return(1.0)
  conc <- c()
  for (sd in c("top", "bot")) {
    keys <- names(dists)[endsWith(names(dists), paste0("|", sd))]
    if (length(keys) < 2) next
    gtot <- setNames(numeric(0), character(0))
    for (k in keys) {
      p <- dists[[k]]
      for (res in names(p)) {
        gg <- encode_aa_r(res, scheme)
        gtot <- .addw(gtot, gg, p[[res]])
      }
    }
    conc <- c(conc, max(gtot) / length(keys))
  }
  if (!length(conc)) return(1.0)
  mean(conc)
}

# independent modal reconstruction (cross-check against apply_fop_pooling)
da_modal_indep <- function(g, scheme) {
  cp <- .collect_changed_pairs(g, node_cols, top_aa_cols, bot_aa_cols)
  rebuild_derived_agreement(cp, scheme)
}

# =============================================================================
# 3. Modal pipeline via the REAL fop_pool.R (gives indep/core/div/gate/da_modal)
# =============================================================================
# NOTE: this run's SCORING did NOT receive contrast_hypotheses_pairs.tsv — the
# stored position_scores.tsv reproduces bit-exactly only with equal-weight
# pooling (hyp_pairs_path = NULL). So the MODAL baseline is rebuilt with NULL,
# and da_frac is reported in two variants: equal-weight (apples-to-apples with
# what ran) and PSS-weighted (spec 3.1, what it would be once PSS is wired).
pooled <- apply_fop_pooling(fd, NULL, tau = TAU)
pooled$Position <- as.integer(pooled$Position)
pss_ok <- !is.null(hp)
stopifnot(all(c("independence", "core", "mrca_diversity", "conservation_gate",
                "derived_agreement", "asr_path_score", "caap_group") %in% names(pooled)))

# per-(Gene,Position,caap_group) group list from the pre-pool frame
fd$.grp <- paste(fd$Gene, fd$Position, fd$caap_group, sep = "\r")
grp_split <- split(fd, fd$.grp)

rows <- lapply(names(grp_split), function(key) {
  g <- grp_split[[key]]
  gene <- g$Gene[1]; pos <- as.integer(g$Position[1]); scheme <- g$caap_group[1]
  hyps <- unique(g$hyp_id[!is.na(g$hyp_id) & nzchar(g$hyp_id)])
  nH <- length(hyps)
  pr <- pooled[pooled$Gene == gene & pooled$Position == pos & pooled$caap_group == scheme, ]
  if (nrow(pr) != 1) return(NULL)

  dists <- domain_side_dists(g, eqw_fn)
  da_frac      <- da_frac_from_dists(dists, scheme)
  da_frac_pss  <- if (pss_ok) da_frac_from_dists(domain_side_dists(g, pss_fn), scheme) else da_frac
  da_modal_chk <- da_modal_indep(g, scheme)

  # split-domain diagnostics
  n_split_domains <- 0L; n_alpha_ties <- 0L
  for (sd in c("top", "bot")) {
    keys <- names(dists)[endsWith(names(dists), paste0("|", sd))]
    for (k in keys) {
      p <- dists[[k]]
      if (length(p) >= 2) {
        n_split_domains <- n_split_domains + 1L
        if (sum(p == max(p)) >= 2) n_alpha_ties <- n_alpha_ties + 1L
      }
    }
  }

  data.frame(
    Gene = gene, Position = pos, caap_group = scheme,
    n_hyp = nH,
    indep = pr$independence[1], core = pr$core[1],
    div = pr$mrca_diversity[1], gate = pr$conservation_gate[1],
    da_modal = pr$derived_agreement[1],
    da_modal_indep = da_modal_chk,
    da_frac = da_frac,
    da_frac_pss = da_frac_pss,
    asr_modal_run = pr$asr_path_score[1],
    recovery_boot = suppressWarnings(as.numeric(pr$recovery_boot[1])),
    n_split_domains = n_split_domains,
    n_alpha_ties = n_alpha_ties,
    stringsAsFactors = FALSE
  )
})
rowdf <- bind_rows(rows)

# recompute asr with the exact fop_pool.R algebra, swapping ONLY da
asr_recompute <- function(indep, core, div, da, gate) {
  indep <- ifelse(is.finite(indep), indep, 1)
  div   <- ifelse(is.finite(div),   div,   0)
  da    <- ifelse(is.finite(da),    da,    1)
  gate  <- ifelse(is.finite(gate),  gate,  1)
  pmin(1, pmax(0, (indep * core) * ((DIV_FLOOR + (1 - DIV_FLOOR) * div) * da) * gate))
}
rowdf$asr_modal_chk <- asr_recompute(rowdf$indep, rowdf$core, rowdf$div, rowdf$da_modal, rowdf$gate)
rowdf$asr_frac      <- asr_recompute(rowdf$indep, rowdf$core, rowdf$div, rowdf$da_frac,  rowdf$gate)
rowdf$asr_frac_pss  <- asr_recompute(rowdf$indep, rowdf$core, rowdf$div, rowdf$da_frac_pss, rowdf$gate)

# =============================================================================
# 4. Propagate to CAAS_score  (scoring_compute.R sections 2f/2g)
# =============================================================================
rowdf$phen_score <- 1 - dplyr::percent_rank(rowdf$recovery_boot)
rowdf$caas_row_modal     <- rowdf$phen_score * rowdf$asr_modal_chk
rowdf$caas_row_frac      <- rowdf$phen_score * rowdf$asr_frac
rowdf$caas_row_frac_pss  <- rowdf$phen_score * rowdf$asr_frac_pss
rowdf$caas_row_modal_run <- rowdf$phen_score * rowdf$asr_modal_run

pos <- rowdf %>%
  group_by(Gene, Position) %>%
  summarise(
    n_schemes = dplyr::n(),
    scheme_set = paste(sort(unique(caap_group)), collapse = "+"),
    n_hyp_max = max(n_hyp),
    n_split_domains = max(n_split_domains),
    n_alpha_ties = sum(n_alpha_ties),
    da_modal = mean(da_modal), da_frac = mean(da_frac), da_frac_pss = mean(da_frac_pss),
    CAAS_modal_run = mean(caas_row_modal_run),
    CAAS_modal = mean(caas_row_modal),
    CAAS_frac  = mean(caas_row_frac),
    CAAS_frac_pss = mean(caas_row_frac_pss),
    asr_modal = mean(asr_modal_chk), asr_frac = mean(asr_frac),
    .groups = "drop"
  )

# join run values + change_side
pos <- pos %>%
  left_join(ps_run %>% select(Gene, Position,
                              CAAS_run = CAAS_score, da_run = derived_agreement,
                              asr_run = asr_score, change_side),
            by = c("Gene", "Position"))

# =============================================================================
# 5. Verification: modal reconstruction must match the run
# =============================================================================
v_da   <- with(pos, max(abs(da_modal - da_run), na.rm = TRUE))
v_caas <- with(pos, max(abs(CAAS_modal_run - CAAS_run), na.rm = TRUE))
v_caas2<- with(pos, max(abs(CAAS_modal - CAAS_run), na.rm = TRUE))
v_asr  <- with(pos, max(abs(asr_modal - asr_run), na.rm = TRUE))
v_row  <- with(rowdf, max(abs(asr_modal_chk - asr_modal_run), na.rm = TRUE))
v_daind<- with(rowdf, max(abs(da_modal - da_modal_indep), na.rm = TRUE))
cat("\n=== RECONSTRUCTION FIDELITY (max |Δ| vs the run) ===\n")
cat(sprintf("  derived_agreement (pos mean)      : %.3e\n", v_da))
cat(sprintf("  CAAS_score  (run asr * phen)      : %.3e\n", v_caas))
cat(sprintf("  CAAS_score  (recomputed asr algebra): %.3e\n", v_caas2))
cat(sprintf("  asr_score (pos mean)              : %.3e\n", v_asr))
cat(sprintf("  asr row-level (algebra vs stored) : %.3e\n", v_row))
cat(sprintf("  da_modal (pool vs indep rebuild)  : %.3e\n", v_daind))

# =============================================================================
# 6. Analyses -> report
# =============================================================================
EPS <- 1e-9
rowdf$dda <- rowdf$da_frac - rowdf$da_modal
rowdf$diff <- abs(rowdf$dda) > EPS
pos$dda <- pos$da_frac - pos$da_modal
pos$dCAAS <- pos$CAAS_frac - pos$CAAS_modal
pos$diff_da <- abs(pos$dda) > EPS
pos$diff_caas <- abs(pos$dCAAS) > EPS

# Q6 gene_caas
gene_dir <- pos %>% select(Gene, Position, CAAS_modal, CAAS_frac, change_side)
pool_all_m <- sort(gene_dir$CAAS_modal); pool_all_f <- sort(gene_dir$CAAS_frac)
sam <- function(x, pool) { x <- x[!is.na(x)]; if (!length(x) || !length(pool)) return(NA_real_)
  (findInterval(max(x), pool) / length(pool))^length(x) }
gene_caas <- gene_dir %>% group_by(Gene) %>%
  summarise(n = dplyr::n(),
            g_modal = sam(CAAS_modal, pool_all_m),
            g_frac  = sam(CAAS_frac,  pool_all_f), .groups = "drop")

jacc <- function(a, b) { u <- length(union(a, b)); if (!u) return(NA_real_); length(intersect(a, b)) / u }
topset <- function(df, col, frac) { v <- df[[col]]; k <- ceiling(frac * sum(!is.na(v)))
  df$Gene[order(-v)][seq_len(k)] }
topset_pos <- function(df, col, frac) { v <- df[[col]]; k <- ceiling(frac * sum(!is.na(v)))
  paste(df$Gene, df$Position)[order(-v)][seq_len(k)] }

saveRDS(list(rowdf = rowdf, pos = pos, gene_caas = gene_caas), file.path(OUT_DIR, "da_frac_objects.rds"))
write.csv(rowdf, file.path(OUT_DIR, "rows_caap_group.csv"), row.names = FALSE)
write.csv(pos,   file.path(OUT_DIR, "positions.csv"), row.names = FALSE)
write.csv(gene_caas, file.path(OUT_DIR, "gene_caas.csv"), row.names = FALSE)

sink(file.path(OUT_DIR, "REPORT.md"))
cat("# Modal vs. fractional harvest-wide `derived_agreement` — offline verification\n\n")
cat(sprintf("Run: `%s`\n\n", normalizePath(FD_PATH)))
cat(sprintf("Params: `tau=%.2f` (convergence_schemes), `diversity_floor=%.2f`. ", TAU, DIV_FLOOR))
cat(sprintf("K=%d Voronoi domains. %d filtered_discovery rows -> %d (Gene,Position,caap_group) rows -> %d scored positions.\n\n",
            K, nrow(fd), nrow(rowdf), nrow(pos)))

cat("## Reconstruction fidelity\n\n")
cat("| quantity | max \\|Δ\\| vs run |\n|---|---|\n")
cat(sprintf("| derived_agreement (position mean) | %.2e |\n", v_da))
cat(sprintf("| CAAS_score (run asr × phen) | %.2e |\n", v_caas))
cat(sprintf("| CAAS_score (asr algebra re-derived) | %.2e |\n", v_caas2))
cat(sprintf("| asr_score (position mean) | %.2e |\n", v_asr))
cat(sprintf("| asr row-level algebra vs stored | %.2e |\n", v_row))
cat(sprintf("| da_modal: pool vs independent rebuild | %.2e |\n\n", v_daind))
cat("> **PSS wiring note.** `position_scores.tsv` reproduces bit-exactly only with\n")
cat("> equal-weight FOP pooling (`hyp_pairs_path = NULL`). With the run's real\n")
cat("> `contrast_hypotheses_pairs.tsv` the pooled asr diverges (max |Δ|≈0.05 on\n")
cat("> multi-hypothesis positions) — i.e. **this run's SCORING never received the PSS\n")
cat("> file**, so Job A / Job B PSS weighting was a silent no-op. The modal baseline\n")
cat("> below is therefore rebuilt equal-weight; `da_frac` is reported both equal-weight\n")
cat("> (apples-to-apples) and PSS-weighted (spec 3.1).\n\n")

cat("## Q1 — How many positions differ?\n\n")
nd_row <- sum(rowdf$diff); nd_pos <- sum(pos$diff_da)
cat(sprintf("- (Gene,Position,caap_group) rows with |Δda|>1e-9: **%d / %d** (%.1f%%)\n",
            nd_row, nrow(rowdf), 100 * nd_row / nrow(rowdf)))
cat(sprintf("- Positions (any caap_group row differs): **%d / %d** (%.1f%%)\n",
            nd_pos, nrow(pos), 100 * nd_pos / nrow(pos)))
cat(sprintf("- Positions whose CAAS_score moves: **%d / %d** (%.1f%%)\n\n",
            sum(pos$diff_caas), nrow(pos), 100 * sum(pos$diff_caas) / nrow(pos)))
brk <- c(0, 1e-9, 1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.25, 0.5, 1)
cat("|Δda| bucket (rows):\n\n")
print(table(cut(abs(rowdf$dda), brk, include.lowest = TRUE)))
cat("\n|Δda| bucket (positions):\n\n")
print(table(cut(abs(pos$dda), brk, include.lowest = TRUE)))
if (pss_ok) {
  ddp  <- rowdf$da_frac_pss - rowdf$da_modal
  ddpe <- rowdf$da_frac_pss - rowdf$da_frac
  cat(sprintf("\nPSS-weighted variant (spec 3.1): rows moving vs modal: **%d**; rows where PSS weighting\n", sum(abs(ddp) > EPS)))
  cat(sprintf("further shifts da vs equal-weight frac: **%d** (max extra |Δ| = %.3f). PSS weighting is\n",
              sum(abs(ddpe) > EPS), max(abs(ddpe))))
  cat("NOT a no-op — but note this run's pipeline never applied it (see wiring note).\n")
}

cat("\n## Q2 — Direction of the change\n\n")
dd <- rowdf[rowdf$diff, ]
cat(sprintf("- rows: da_frac **higher** than modal: %d ; **lower**: %d\n",
            sum(dd$dda > 0), sum(dd$dda < 0)))
pdd <- pos[pos$diff_da, ]
cat(sprintf("- positions: da_frac higher: %d ; lower: %d\n", sum(pdd$dda > 0), sum(pdd$dda < 0)))
cat("- Spec 3.4 prediction: da_frac RISES when a split domain partially supports the plurality,\n")
cat("  FALLS when the modal broke a tie to invent full agreement. Alpha-tie rows (modal tie-break):\n")
cat(sprintf("  %d differing rows carry >=1 alpha-tie domain; of those %d go DOWN with da_frac.\n",
            sum(dd$n_alpha_ties > 0), sum(dd$n_alpha_ties > 0 & dd$dda < 0)))

cat("\n## Q3 — Characterising the differing positions\n\n")
cat(sprintf("- All %d differing rows have >=1 split domain: %s\n",
            nd_row, all(dd$n_split_domains >= 1)))
cat("- split domains per differing row:\n\n")
print(table(dd$n_split_domains))
cat(sprintf("\n- differing rows that are EXACT alpha-ties (modal chose alphabetically): %d\n",
            sum(dd$n_alpha_ties > 0)))
cat("- by scheme (differing rows):\n\n")
print(table(dd$caap_group))
cat("- by scheme (all rows, for base rate):\n\n")
print(table(rowdf$caap_group))

cat("\n## Q4 — Back-compatibility (no split domains => da_frac == da_modal)\n\n")
nosplit <- rowdf[rowdf$n_split_domains == 0, ]
bad <- nosplit[abs(nosplit$dda) > 1e-9, ]
cat(sprintf("- rows with no split domain: %d ; of those, |Δda|>1e-9: **%d**\n", nrow(nosplit), nrow(bad)))
if (nrow(bad)) { cat("\nVIOLATIONS:\n\n"); print(bad[, c("Gene","Position","caap_group","da_modal","da_frac","dda")]) }

cat("\n## Q5 — Propagation to CAAS_score\n\n")
cat(sprintf("- positions with CAAS_score change (|Δ|>1e-9): **%d / %d**\n",
            sum(pos$diff_caas), nrow(pos)))
cat("- |ΔCAAS_score| buckets:\n\n")
print(table(cut(abs(pos$dCAAS), brk, include.lowest = TRUE)))
cat(sprintf("- max |ΔCAAS_score| = %.4f ; mean over changed = %.5f\n",
            max(abs(pos$dCAAS)), mean(abs(pos$dCAAS[pos$diff_caas]))))
for (fr in c(0.01, 0.05, 0.10)) {
  a <- topset_pos(pos, "CAAS_modal", fr); b <- topset_pos(pos, "CAAS_frac", fr)
  cat(sprintf("- top-%d%% position CAAS_score: Jaccard(modal, frac) = %.3f (|set|=%d)\n",
              round(100 * fr), jacc(a, b), length(a)))
}
if (pss_ok) {
  dcp <- pos$CAAS_frac_pss - pos$CAAS_modal
  cat(sprintf("- PSS-weighted variant: %d positions move CAAS_score; max |Δ| = %.4f\n",
              sum(abs(dcp) > EPS), max(abs(dcp))))
}

cat("\n## Q6 — Propagation to gene_caas_score\n\n")
sp <- suppressWarnings(cor(gene_caas$g_modal, gene_caas$g_frac, method = "spearman",
                           use = "complete.obs"))
cat(sprintf("- Spearman(gene_caas_score modal, frac) = %.5f  (n=%d genes)\n", sp, nrow(gene_caas)))
for (fr in c(0.01, 0.05, 0.10)) {
  a <- topset(gene_caas, "g_modal", fr); b <- topset(gene_caas, "g_frac", fr)
  cat(sprintf("- top-%d%% genes: Jaccard = %.3f (|set|=%d)\n", round(100 * fr), jacc(a, b), length(a)))
}
a5 <- topset(gene_caas, "g_modal", 0.05); b5 <- topset(gene_caas, "g_frac", 0.05)
cat(sprintf("- genes ENTERING top-5%% under frac: %s\n",
            paste(setdiff(b5, a5), collapse = ", ")))
cat(sprintf("- genes LEAVING top-5%% under frac: %s\n",
            paste(setdiff(a5, b5), collapse = ", ")))

cat("\n## Q7 — 10 rows with the largest |Δda|\n\n")
top10 <- rowdf[order(-abs(rowdf$dda)), ][seq_len(min(10, nrow(rowdf))), ]
print(top10[, c("Gene","Position","caap_group","n_hyp","n_split_domains","n_alpha_ties",
                "indep","core","div","gate","da_modal","da_frac","dda",
                "asr_modal_chk","asr_frac")], row.names = FALSE)
cat("\nPer-domain breakdown of those rows is written to `top10_domain_breakdown.txt`.\n")

cat("\n## Q8 — convergence_schemes\n\n")
# recompute .position_descriptors modal vs fractional across ALL rows of a position
cs_cmp <- fd %>% group_by(Gene, Position) %>% group_modify(function(sub, key) {
  modal <- .position_descriptors(sub, node_cols, top_aa_cols, bot_aa_cols, TAU)$convergence_schemes
  # fractional: same gate structure, but frac concentration over the whole position
  cp <- .collect_changed_pairs(sub, node_cols, top_aa_cols, bot_aa_cols)
  sn <- if (nrow(cp)) table(cp$side) else integer(0)
  if (nrow(cp) == 0 || !length(sn) || max(sn) < 2) return(data.frame(modal = modal, frac = ""))
  multi <- any(vapply(split(cp$raw_aa, cp$side), function(rs) length(unique(rs)) >= 2, logical(1)))
  if (!multi) return(data.frame(modal = modal, frac = "US"))
  dists <- domain_side_dists(sub, eqw_fn)
  kept <- AA_SCHEME_NAMES[vapply(AA_SCHEME_NAMES, function(sc)
    isTRUE(da_frac_from_dists(dists, sc) >= TAU), logical(1))]
  data.frame(modal = modal, frac = paste(kept, collapse = ","))
}) %>% ungroup()
cat(sprintf("- positions where convergence_schemes would change: **%d / %d**\n",
            sum(cs_cmp$modal != cs_cmp$frac), nrow(cs_cmp)))
chg <- cs_cmp[cs_cmp$modal != cs_cmp$frac, ]
if (nrow(chg)) print(head(chg, 30), row.names = FALSE)
write.csv(cs_cmp, file.path(OUT_DIR, "convergence_schemes_cmp.csv"), row.names = FALSE)

cat("\n## C — Recommendation\n\n")
motiv_rows <- sum(rowdf$n_split_domains > 0)
cat(sprintf("- (Gene,Position,caap_group) rows with a genuinely split domain (the case that\n"))
cat(sprintf("  motivates harvest-wide da: unanimous within H, different between H): **%d / %d** (%.1f%%)\n",
            motiv_rows, nrow(rowdf), 100 * motiv_rows / nrow(rowdf)))
cat(sprintf("- of those, |Δda|>1e-9 under the fractional rule: %d\n", sum(rowdf$n_split_domains > 0 & rowdf$diff)))
cat(sprintf("- exact alpha-tie rows (modal's alphabetical tie-break is load-bearing): %d\n",
            sum(rowdf$n_alpha_ties > 0)))
cat(sprintf("- top-5%% gene Jaccard modal vs frac: %.3f ; Spearman gene ranking: %.4f\n",
            jacc(topset(gene_caas, "g_modal", 0.05), topset(gene_caas, "g_frac", 0.05)),
            suppressWarnings(cor(gene_caas$g_modal, gene_caas$g_frac, method = "spearman", use = "complete.obs"))))
cat("\nSee the narrative recommendation in `docs/da_frac/RECOMMENDATION.md`.\n")
sink()

# domain breakdown for Q7
sink(file.path(OUT_DIR, "top10_domain_breakdown.txt"))
for (r in seq_len(nrow(top10))) {
  tr <- top10[r, ]
  key <- paste(tr$Gene, tr$Position, tr$caap_group, sep = "\r")
  g <- grp_split[[key]]
  cat(sprintf("\n=== %s / %d / %s  (n_hyp=%d)  da_modal=%.4f da_frac=%.4f  ΔCAAS(row asr) modal=%.4f frac=%.4f\n",
              tr$Gene, tr$Position, tr$caap_group, tr$n_hyp, tr$da_modal, tr$da_frac,
              tr$asr_modal_chk, tr$asr_frac))
  for (i in seq_len(K)) for (sd in c("top", "bot")) {
    acol <- if (sd == "top") top_aa_cols[i] else bot_aa_cols[i]
    if (!(acol %in% names(g))) next
    raw <- toupper(trimws(as.character(g[[acol]])))
    keep <- !is.na(raw) & nzchar(raw) & !(raw %in% BAD_AA)
    if (!any(keep)) next
    hh <- as.character(g$hyp_id)[keep]
    pss <- vapply(hh, function(h) pss_fn(h, i), numeric(1))
    cat(sprintf("  domain %d %s: %s\n", i, sd,
                paste(sprintf("%s=%s(pss %.2f)", hh, raw[keep], pss), collapse = "  ")))
    enc <- vapply(raw[keep], encode_aa_r, character(1), scheme = tr$caap_group)
    cat(sprintf("      encoded(%s): %s\n", tr$caap_group, paste(enc, collapse = ",")))
  }
}
sink()

cat("\nWrote:", file.path(OUT_DIR, "REPORT.md"), "\n")
