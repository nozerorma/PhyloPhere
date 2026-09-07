#!/usr/bin/env Rscript
# Tests for fop_pool.R — run: Rscript test_fop_pool.R
suppressPackageStartupMessages(library(dplyr))
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "fop_pool.R"))

ok <- TRUE
check <- function(cond, msg) {
  cat(if (isTRUE(cond)) "PASS  " else "FAIL  ", msg, "\n", sep = "")
  if (!isTRUE(cond)) ok <<- FALSE
}
approx <- function(a, b, tol = 1e-6) is.finite(a) && is.finite(b) && abs(a - b) < tol

# ── _p_at_least_2 matches the Python algebra ────────────────────────────────
check(approx(.p_at_least_2(c(1, 1)), 1),                 "_p_at_least_2(1,1) == 1")
check(approx(.p_at_least_2(c(0.5, 0.5)), 0.25),          "_p_at_least_2(.5,.5) == .25")
check(approx(.p_at_least_2(c(0.9, 0.1)), 0.09),          "_p_at_least_2(.9,.1) == .09")
check(approx(.p_at_least_2(c(0.8, 0.8)), 0.64),          "_p_at_least_2(.8,.8) == .64")
check(approx(.p_at_least_2(c(1)), 0),                    "_p_at_least_2 single == 0")
check(approx(.p_at_least_2(c(0.8, 0.8, 0.8)), 0.896),    "_p_at_least_2(.8x3) == .896")

# ── .wmean ─────────────────────────────────────────────────────────────────
check(approx(.wmean(c(1, 0), c(3, 1)), 0.75),            ".wmean weighted")
check(approx(.wmean(c(1, 0), c(NA, NA)), 0.5),           ".wmean falls back to plain mean")
check(approx(.wmean(c(0.4, NA), c(1, 1)), 0.4),          ".wmean drops non-finite x")

# ── single-hypothesis input is a no-op ─────────────────────────────────────
one <- data.frame(
  Gene = "G", Position = 10L, caap_group = "US", hyp_id = NA_character_,
  asr_path_score = 0.42, independence = 0.9, mrca_diversity = 0.5,
  derived_agreement = 1, conservation_gate = 1, core = 0.6,
  mrca_1_path_score = 0.8, mrca_1_node = "n1",
  mrca_2_path_score = 0.7, mrca_2_node = "n2",
  stringsAsFactors = FALSE
)
res1 <- apply_fop_pooling(one, NULL)
check(nrow(res1) == 1, "single-hyp: one row out")
check(approx(res1$asr_path_score, 0.42), "single-hyp: asr_path_score untouched")
check("n_hypotheses" %in% names(res1) && res1$n_hypotheses == 0, "single-hyp: n_hypotheses == 0")

# ── FOP: two hypotheses, two domains — two-job PSS weighting ───────────────
# H1 rides a strong domain-1 pair (node a, s=0.9, own pss 10) but a WEAK
# domain 2 (own pss 1).  H2's domain-1 pair (node b, s=0.2) has own pss 2.
# Domain 2 is the same pair (node c, s=0.8) for both hypotheses.
#   Job A (per-domain c_i): each pair weighted by its OWN pss in that domain.
#     c_1 = wmean(s=(0.9,0.2), w=(pss(H1,1)=10, pss(H2,1)=2))
#         = (9 + 0.4) / 12 = 9.4/12 = 0.7833333
#     c_2: dedup by node c -> w = max(pss(H1,2)=1, pss(H2,2)=2) = 2, single s
#         -> c_2 = 0.8
#     (a min-of-domains weight would have given H1 weight 1 -> c_1 ~= 0.433;
#      pair-own PSS keeps the strong pair dominant.)
#   Job B (axis pools): per-hypothesis MEAN pss across domains.
#     w(H1) = mean(10, 1) = 5.5 ; w(H2) = mean(2, 2) = 2
#     diversity_pooled = wmean((1, 0.5), w=(5.5, 2)) = (5.5 + 1)/7.5 = 0.8666667
#     (a min weight gave w(H1)=1 -> diversity ~= 0.667; mean-PSS lifts H1.)
fop <- data.frame(
  Gene = "G", Position = 20L, caap_group = "US",
  hyp_id = c("H1", "H2"),
  asr_path_score = c(0.7, 0.2),
  independence = c(1, 1), mrca_diversity = c(1, 0.5),
  derived_agreement = c(1, 1), conservation_gate = c(1, 1),
  core = c(0.72, 0.24),
  mrca_1_path_score = c(0.9, 0.2), mrca_1_node = c("a", "b"),
  mrca_2_path_score = c(0.8, 0.8), mrca_2_node = c("c", "c"),
  stringsAsFactors = FALSE
)
hp <- data.frame(
  hypothesis_id = c("H1", "H1", "H2", "H2"),
  pair = c(1L, 2L, 1L, 2L),
  pss_score = c(10, 1, 2, 2),
  stringsAsFactors = FALSE
)
hp_path <- tempfile(fileext = ".tsv")
write.table(hp, hp_path, sep = "\t", quote = FALSE, row.names = FALSE)

resF <- apply_fop_pooling(fop, hp_path)
c1 <- 9.4 / 12       # 0.7833333
c2 <- 0.8
div <- 6.5 / 7.5     # 0.8666667
check(nrow(resF) == 1, "FOP: collapses to one row")
check(resF$n_hypotheses == 2, "FOP: n_hypotheses == 2")
check(approx(resF$pooled_domain_1, c1), "FOP job A: c_1 pair-own-PSS weighted == 9.4/12")
check(approx(resF$pooled_domain_2, c2), "FOP job A: c_2 (one distinct pair) == 0.8")
check(approx(resF$core, .p_at_least_2(c(c1, c2))), "FOP: pooled core == p>=2(c_1,c_2)")
check(approx(resF$mrca_diversity, div), "FOP job B: pooled diversity == 6.5/7.5 (mean-PSS)")
exp_asr <- 1 * .p_at_least_2(c(c1, c2)) * (0.75 + 0.25 * div) * 1 * 1
check(approx(resF$asr_path_score, exp_asr), sprintf("FOP: pooled asr == %.5f", exp_asr))
check(resF$asr_path_score > min(fop$asr_path_score) &&
      resF$asr_path_score < max(fop$asr_path_score),
      "FOP: pooled asr sits between the hypotheses")
# Job A is pair-PSS-driven: the s=0.9 pair (own pss 10) dominates c_1 even
# though it sits in H1 whose OTHER domain is weak (pss 1). A hypothesis-min
# weight would have sunk it to ~0.43.
check(resF$pooled_domain_1 > 0.75,
      "FOP job A: high-own-PSS pair dominates its c_i despite a weak sibling domain")
# Job B is mean-PSS-driven: H1 (mean pss 5.5) outweighs H2 (2) on the axes,
# so pooled diversity leans toward H1's 1.0, well above the min-weighted 0.667.
check(resF$mrca_diversity > 0.8,
      "FOP job B: mean-PSS weight lifts pooled diversity toward the higher-PSS hypothesis")

# ── FOP without a PSS file -> equal weights ────────────────────────────────
resE <- apply_fop_pooling(fop, NULL)
# c_1 = mean(0.9, 0.2) = 0.55 ; c_2 = 0.8
check(approx(resE$pooled_domain_1, 0.55), "FOP no-PSS: c_1 equal-weight == 0.55")
check(approx(resE$core, .p_at_least_2(c(0.55, 0.8))), "FOP no-PSS: core from equal-weight c")
check(approx(resE$mrca_diversity, 0.75), "FOP no-PSS: diversity equal-weight == 0.75")

# ── DIRECTIONAL core fix ────────────────────────────────────────────────────
# Two hypotheses, same pair in each of 2 domains (dedup is a no-op here so the
# fix is isolated to the direction split). Domain 1's pair changed on TOP
# (0.85); domain 2's DIFFERENT pair changed on BOTTOM (0.80). These are two
# unrelated, opposite-direction single events — not "2 domains agreeing" — so
# the directional core must be 0 (each side has < 2 observations, exactly like
# path_scores.py's own core_top/core_bottom). The pre-fix collapsed pool would
# have scored core = p_at_least_2(.85, .80) =~ 0.68, treating them as if they
# corroborated each other.
fop_dir_opp <- data.frame(
  Gene = "G", Position = 40L, caap_group = "US",
  hyp_id = c("H1", "H2"),
  asr_path_score = c(0.5, 0.5),
  independence = c(1, 1), mrca_diversity = c(0, 0),
  derived_agreement = c(1, 1), conservation_gate = c(1, 1),
  core = c(0.5, 0.5),
  mrca_1_path_score = c(0.85, 0.85), mrca_1_node = c("d1", "d1"),
  mrca_1_top_path_score = c(0.85, 0.85),
  mrca_2_path_score = c(0.80, 0.80), mrca_2_node = c("d2", "d2"),
  mrca_2_bot_path_score = c(0.80, 0.80),
  stringsAsFactors = FALSE
)
resDirOpp <- apply_fop_pooling(fop_dir_opp, NULL)
check(approx(resDirOpp$core, 0),
      "DIRECTIONAL: a top-only domain and a bottom-only domain do NOT corroborate -> core == 0")
check(approx(resDirOpp$pooled_domain_1, 0.85) && approx(resDirOpp$pooled_domain_2, 0.80),
      "DIRECTIONAL: display column pooled_domain_<i> is unaffected (still side-averaged)")

# Same fixture MINUS the directional columns -> old collapsed fallback, which
# does over-credit the opposite-direction pair (documented known approximation).
fop_dir_opp_nocols <- fop_dir_opp[, setdiff(
  names(fop_dir_opp), c("mrca_1_top_path_score", "mrca_2_bot_path_score")
)]
resDirOppFallback <- apply_fop_pooling(fop_dir_opp_nocols, NULL)
check(approx(resDirOppFallback$core, .p_at_least_2(c(0.85, 0.80))),
      "DIRECTIONAL fallback: no side columns -> collapsed pool (core ~= 0.68)")
check(resDirOppFallback$core - resDirOpp$core > 0.3,
      "DIRECTIONAL: the fix materially lowers core vs the pre-fix approximation")

# Contrast: the SAME two scores, but now BOTH on the top side (genuine
# convergence). The directional core must equal the plain p_at_least_2 of the
# two — nothing is lost when domains actually agree in direction.
fop_dir_same <- fop_dir_opp
names(fop_dir_same)[names(fop_dir_same) == "mrca_2_bot_path_score"] <- "mrca_2_top_path_score"
resDirSame <- apply_fop_pooling(fop_dir_same, NULL)
check(approx(resDirSame$core, .p_at_least_2(c(0.85, 0.80))),
      "DIRECTIONAL: two domains agreeing on the SAME side still corroborate normally")

# ── POINT 2: conservation_gate rebuilt from deduplicated conserved pairs ────
# Three hypotheses. H1 and H2 SHARE one conserved pair (same conserved_1_node
# "cn_shared", cons 0.6); H3 carries a hypothesis-unique conserved pair
# ("cn_uniq", cons 0.4). Distinct conserved pairs = {0.6, 0.4} -> the rebuilt
# gate must be 0.5 + 0.5*mean(0.6,0.4) = 0.75, counting the shared pair ONCE —
# NOT the 3-way mean of per-hypothesis gates (0.8, 0.8, 0.7) = 0.76667.
fopc <- data.frame(
  Gene = "G", Position = 30L, caap_group = "US",
  hyp_id = c("H1", "H2", "H3"),
  asr_path_score = c(0.5, 0.5, 0.4),
  independence = c(1, 1, 1), mrca_diversity = c(0, 0, 0),
  derived_agreement = c(1, 1, 1),
  conservation_gate = c(0.8, 0.8, 0.7),   # 0.5 + 0.5*cons per hypothesis
  core = c(0.5, 0.5, 0.4),
  mrca_1_path_score = c(0.7, 0.7, 0.6), mrca_1_node = c("a", "a", "a"),
  mrca_2_path_score = c(0.6, 0.6, 0.5), mrca_2_node = c("b", "b", "b"),
  conserved_1_node = c("cn_shared", "cn_shared", "cn_uniq"),
  conserved_1_cons = c(0.6, 0.6, 0.4),
  stringsAsFactors = FALSE
)
resC <- apply_fop_pooling(fopc, NULL)
check(nrow(resC) == 1, "POINT2: collapses to one row")
check(approx(resC$conservation_gate, 0.75),
      "POINT2: gate == 0.5+0.5*mean(0.6,0.4) == 0.75 (shared pair deduped)")
check(!approx(resC$conservation_gate, mean(c(0.8, 0.8, 0.7))),
      "POINT2: gate is NOT the 3-way per-hypothesis gate mean (0.76667)")

# No conserved_<j>_* columns -> fallback to the per-hypothesis-gate wmean.
fopc_nocols <- fopc[, setdiff(names(fopc), c("conserved_1_node", "conserved_1_cons"))]
resCn <- apply_fop_pooling(fopc_nocols, NULL)
check(approx(resCn$conservation_gate, mean(c(0.8, 0.8, 0.7))),
      "POINT2 fallback: no conserved cols -> per-hypothesis gate mean")

# All-conserved-columns-empty group (FOP structure, novel position) -> gate 1.0.
fopc_empty <- fopc
fopc_empty$conserved_1_node <- NA_character_
fopc_empty$conserved_1_cons <- NA_real_
resCe <- apply_fop_pooling(fopc_empty, NULL)
check(approx(resCe$conservation_gate, 1.0),
      "POINT2: no distinct conserved pairs -> gate == 1.0")

# single-hypothesis with conserved cols present stays a strict no-op
one_c <- one
one_c$conserved_1_node <- "cn"
one_c$conserved_1_cons <- 0.3
res1c <- apply_fop_pooling(one_c, NULL)
check(nrow(res1c) == 1 && approx(res1c$asr_path_score, 0.42),
      "POINT2: single-hyp + conserved cols still a no-op")

# ── POINT 3: harvest-wide, scheme-resolved derived_agreement ───────────────
# BRCA1/96, 3 hypotheses, K=2 domains. Domain 1's pair is I in every hypothesis,
# domain 2's pair is V in every hypothesis (each domain unanimous across the
# harvest -> p_i is a point mass -> the fractional rule reduces to the modal
# plurality exactly). No top-side change. Ancestral "A" throughout.
#   US : bot {I},{V} -> max_g(sum p_i^enc(g))/2 = 1/2 = 0.5  -> harvest-wide da = 0.5
#   GS4: I,V -> both 'h'  -> 2/2 = 1.0        -> harvest-wide da = 1.0
#   GS3: I,V -> both 'l'  -> 1.0
#   GS1/GS2: I,V split    -> 0.5
# convergence_schemes @ tau 0.8 -> exactly {GS3, GS4}.
mk_p3 <- function(grp) data.frame(
  Gene = "BRCA1", Position = 96L, caap_group = grp,
  hyp_id = c("H1", "H2", "H3"),
  asr_path_score = c(0.5, 0.4, 0.45),
  independence = c(1, 1, 1), mrca_diversity = c(0, 0, 0),
  derived_agreement = c(1, 1, 1),          # each hypothesis internally unanimous
  conservation_gate = c(1, 1, 1),
  core = c(0.5, 0.4, 0.45),
  mrca_1_path_score = c(0.7, 0.7, 0.7), mrca_1_node = c("p1", "p3", "p1"),
  mrca_2_path_score = c(0.6, 0.6, 0.6), mrca_2_node = c("p2", "p4", "p2"),
  mrca_1_anc_aa = "A", mrca_2_anc_aa = "A",
  mrca_1_top_aa = "",  mrca_2_top_aa = "",
  mrca_1_bot_aa = c("I", "I", "I"), mrca_2_bot_aa = c("V", "V", "V"),
  stringsAsFactors = FALSE
)
resP3 <- apply_fop_pooling(mk_p3("US"), NULL)
check(nrow(resP3) == 1, "POINT3: collapses to one row")
check(approx(resP3$derived_agreement, 0.5),
      "POINT3: harvest-wide da under US == 0.5 (I/V split, plurality 2/4)")
check(identical(resP3$convergence_schemes, "GS3,GS4"),
      "POINT3: convergence_schemes == GS3,GS4 (da >= 0.8 only there)")

# derived_residues / {top,bottom}_residue_support are produced UPSTREAM now
# (CT_POSTPROC residue_descriptors.py, tested there); apply_fop_pooling only
# carries them through. Absent in this fixture -> stable empty strings.
check(identical(resP3$derived_residues, ""),
      "POINT3: derived_residues empty when the input frame lacks it")
mk_p3_carry <- mk_p3("US"); mk_p3_carry$derived_residues <- "A/IV"
mk_p3_carry$top_residue_support <- ""; mk_p3_carry$bottom_residue_support <- "I:2,V:2"
resP3c <- apply_fop_pooling(mk_p3_carry, NULL)
check(identical(resP3c$derived_residues, "A/IV") &&
      identical(resP3c$bottom_residue_support, "I:2,V:2"),
      "POINT3: apply_fop_pooling carries upstream descriptor columns through")

resP3g <- apply_fop_pooling(mk_p3("GS4"), NULL)
check(approx(resP3g$derived_agreement, 1.0),
      "POINT3: harvest-wide da under GS4 == 1.0 (I,V co-encode to 'h')")
check(identical(resP3g$convergence_schemes, "GS3,GS4"),
      "POINT3: convergence_schemes is scheme-independent (raw residues)")

# Fallback: no *_aa columns -> old per-hypothesis wmean da, empty descriptors.
p3_nocols <- mk_p3("US")[, grep("_aa$", names(mk_p3("US")), invert = TRUE)]
resP3n <- apply_fop_pooling(p3_nocols, NULL)
check(approx(resP3n$derived_agreement, 1.0),
      "POINT3 fallback: no *_aa cols -> per-hypothesis da wmean (all 1.0)")
check(identical(resP3n$convergence_schemes, ""),
      "POINT3 fallback: convergence_schemes empty without raw residues")

# Non-FOP input still gets the descriptor columns for a stable schema.
check(all(c("derived_residues", "top_residue_support", "bottom_residue_support",
            "convergence_schemes") %in% names(res1)),
      "POINT3: single-hyp output carries descriptor columns")
check(identical(res1$convergence_schemes, ""), "POINT3: single-hyp descriptors empty")

# tau argument is honoured.
resP3t <- apply_fop_pooling(mk_p3("US"), NULL, tau = 0.4)
check(identical(resP3t$convergence_schemes, "US,GS1,GS2,GS3,GS4"),
      "POINT3: tau=0.4 admits all five schemes")

# ── convergence_schemes: identity vs chemical vs no-statement ───────────────
# Single exact derived residue across 2 pairs -> "US" ONLY (not all five). The
# GS entries would be tautological: a single residue is 100%-concentrated under
# any partition.
mk_single <- function() data.frame(
  Gene = "G", Position = 5L, caap_group = "US",
  hyp_id = c("H1", "H2", "H3"),
  asr_path_score = c(0.5, 0.5, 0.5), independence = 1, mrca_diversity = 0,
  derived_agreement = 1, conservation_gate = 1, core = 0.5,
  mrca_1_path_score = c(0.7, 0.7, 0.7), mrca_1_node = c("p1", "p3", "p1"),
  mrca_2_path_score = c(0.6, 0.6, 0.6), mrca_2_node = c("p2", "p4", "p2"),
  mrca_1_anc_aa = "L", mrca_2_anc_aa = "L",
  mrca_1_top_aa = c("I", "I", "I"), mrca_2_top_aa = c("I", "I", "I"),
  mrca_1_bot_aa = "", mrca_2_bot_aa = "",
  stringsAsFactors = FALSE
)
check(identical(apply_fop_pooling(mk_single(), NULL)$convergence_schemes, "US"),
      "POINT3: single exact derived residue over 2 pairs -> 'US' only")

# Fewer than 2 changed pairs on any side -> "" (no convergence statement).
mk_onepair <- mk_single()
mk_onepair$mrca_2_top_aa <- ""   # kill the second changed pair
check(identical(apply_fop_pooling(mk_onepair, NULL)$convergence_schemes, ""),
      "POINT3: <2 changed pairs -> convergence_schemes empty")

# convergence_schemes identity gate is MODAL by design (a structural "is there
# any plurality-level disagreement" question): pair 1 is I under H1/H3 but a
# stray V under H2 -> modal I wins -> still identity convergence with pair 2.
mk_hypdisagree <- mk_single()
mk_hypdisagree$mrca_1_top_aa <- c("I", "V", "I")
check(identical(apply_fop_pooling(mk_hypdisagree, NULL)$convergence_schemes, "US"),
      "POINT3: modal identity gate survives one disagreeing hypothesis")

# ── POINT 3 fractional: distribution vs modal ──────────────────────────────────
# Helper: a K-domain, bottom-only harvest group, one residue vector per domain.
mk_frac <- function(dom_residues, grp = "US", pss = NULL) {
  n <- length(dom_residues[[1]])
  K <- length(dom_residues)
  d <- data.frame(
    Gene = "G", Position = 7L, caap_group = grp,
    hyp_id = paste0("H", seq_len(n)),
    asr_path_score = rep(0.5, n), independence = 1, mrca_diversity = 0,
    derived_agreement = 1, conservation_gate = 1, core = 0.5,
    stringsAsFactors = FALSE
  )
  for (i in seq_len(K)) {
    d[[paste0("mrca_", i, "_path_score")]] <- 0.6
    d[[paste0("mrca_", i, "_node")]]       <- paste0("n", i)
    d[[paste0("mrca_", i, "_anc_aa")]]     <- "A"
    d[[paste0("mrca_", i, "_top_aa")]]     <- ""
    d[[paste0("mrca_", i, "_bot_aa")]]     <- dom_residues[[i]]
  }
  d
}

# Partial alignment (spec 3.4): p_N={I:1}, p_C={I:1}, p_S={I:0.4, L:0.6} under US
#   fractional -> (1 + 1 + 0.4) / 3 = 0.80   (modal would tie-break S->L -> 2/3 = 0.667)
resFa <- apply_fop_pooling(
  mk_frac(list(rep("I", 5), rep("I", 5), c("I", "I", "L", "L", "L"))), NULL)
check(approx(resFa$derived_agreement, 0.80),
      "POINT3 frac: partial alignment -> da == 0.80 (modal gives 0.667)")

# Genuine ambiguity (spec 3.4): p_N={I:1}, p_C={I:0.5,L:0.5}, p_S={L:1} under US
#   fractional -> (1 + 0.5 [+0.5 for L] ...) max_g(I:1.5, L:1.5)/3 = 0.50
#   (modal tie-breaks C->I and fabricates {I,I,L} -> 2/3 = 0.667)
resAmb <- apply_fop_pooling(
  mk_frac(list(c("I", "I"), c("I", "L"), c("L", "L"))), NULL)
check(approx(resAmb$derived_agreement, 0.50),
      "POINT3 frac: genuine ambiguity -> da == 0.50 (modal fabricates 0.667)")

# Byte-identity: every domain unanimous across the harvest -> fractional == modal
# for every scheme. Both sides change here.
mk_unanim <- mk_frac(list(rep("I", 3), rep("V", 3)))
mk_unanim$mrca_1_top_aa <- rep("F", 3); mk_unanim$mrca_2_top_aa <- rep("Y", 3)
ncU <- grep("^mrca_\\d+_node$", names(mk_unanim), value = TRUE)
taU <- grep("^mrca_\\d+_top_aa$", names(mk_unanim), value = TRUE)
baU <- grep("^mrca_\\d+_bot_aa$", names(mk_unanim), value = TRUE)
cpU    <- .collect_changed_pairs(mk_unanim, ncU, taU, baU)
distsU <- .collect_changed_pair_dists(mk_unanim, ncU, taU, baU)
for (sc in AA_SCHEME_NAMES) {
  check(approx(rebuild_derived_agreement_frac(distsU, sc),
              rebuild_derived_agreement(cpU, sc)),
        sprintf("POINT3 frac: byte-identical to modal on unanimous domains (%s)", sc))
}

# PSS weighting (spec 3.1): domain 1 = H1:"I"(pss 3) / H2:"L"(pss 1) -> {I:0.75, L:0.25};
# domain 2 = "I" for both. US bottom: max_g(I: 0.75 + 1, L: 0.25) / 2 = 1.75/2 = 0.875.
# Equal weight -> {I:0.5, L:0.5} -> 1.5/2 = 0.75.
mk_pssw <- mk_frac(list(c("I", "L"), c("I", "I")))
hp_pssw <- data.frame(
  hypothesis_id = c("H1", "H2", "H1", "H2"),
  pair = c(1L, 1L, 2L, 2L), pss_score = c(3, 1, 2, 2),
  stringsAsFactors = FALSE)
hp_pssw_path <- tempfile(fileext = ".tsv")
write.table(hp_pssw, hp_pssw_path, sep = "\t", quote = FALSE, row.names = FALSE)
check(approx(apply_fop_pooling(mk_pssw, hp_pssw_path)$derived_agreement, 0.875),
      "POINT3 frac: PSS-weighted da == 0.875 (I gets 0.75 of domain 1)")
check(approx(apply_fop_pooling(mk_pssw, NULL)$derived_agreement, 0.75),
      "POINT3 frac: equal-weight da == 0.75 (same fixture, no PSS file)")

# .position_descriptors picks up the same PSS instrument. domain 1 = I always;
# domain 2 = V,V,V,I (modal V, so cp is multi-residue {I,V} -> GS-list branch).
# I/V split under GS1/GS2 but co-encode under GS3/GS4.
#   equal weight: p_2 = {V:0.75, I:0.25} -> GS1 concentration (1 + 0.25)/2 = 0.625 < tau
#     -> convergence_schemes == "GS3,GS4"
#   weight the lone-I hypothesis (H4) x9: p_2 = {I:0.75, V:0.25} -> GS1 concentration
#     (1 + 0.75)/2 = 0.875 >= tau -> GS1/GS2 also clear.
mk_cs <- mk_frac(list(rep("I", 4), c("V", "V", "V", "I")))
ncC <- grep("^mrca_\\d+_node$", names(mk_cs), value = TRUE)
taC <- grep("^mrca_\\d+_top_aa$", names(mk_cs), value = TRUE)
baC <- grep("^mrca_\\d+_bot_aa$", names(mk_cs), value = TRUE)
cs_eq <- .position_descriptors(mk_cs, ncC, taC, baC, tau = 0.8)$convergence_schemes
cs_wI <- .position_descriptors(mk_cs, ncC, taC, baC, tau = 0.8,
                               pair_pss = function(h, i) if (h == "H4") 9 else 1)$convergence_schemes
check(identical(cs_eq, "GS3,GS4") && identical(cs_wI, "US,GS1,GS2,GS3,GS4"),
      "POINT3 frac: .position_descriptors convergence_schemes responds to pair_pss weighting")

cat("\n", if (ok) "ALL TESTS PASSED" else "SOME TESTS FAILED", "\n", sep = "")
quit(status = if (ok) 0 else 1)
