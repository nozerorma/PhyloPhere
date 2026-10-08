#!/usr/bin/env Rscript
# Calibration of the position p-values against the run's own permulation null.
#
# Each null cycle is taken in turn as the "observed" data and judged against the other N - 1 cycles, with the rules of
# scoring_compute.R: p.emp = (k + 1) / N and the factorized p of .fact_fit / .fact_assign_obs / .fact_p. If the null is
# exchangeable with the observed data, the share of (position, cycle) pairs with p <= alpha cannot exceed alpha, overall
# or within a class of propensity (cycles that score the position). A position the pseudo-observed cycle does not score
# has p = 1. The helpers are read from scoring_compute.R itself, so the check follows the code and not a copy of it.
#
# Usage: Rscript pepc_null_calibration.R <results/<trait>_complete> <out_prefix>
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
res_dir <- args[1]; out_prefix <- args[2]
repo <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), "..", "..", ".."))

# Helpers (TIE_TOL, FACT_*, .fact_*) taken verbatim from the scoring step
src <- readLines(file.path(repo, "subworkflows/SCORING/local/src/scoring_compute.R"))
a <- grep("^TIE_TOL <- ", src); b <- grep("^# ── end of the factorized p helpers", src) - 1L
eval(parse(text = src[a:b]))

nul <- read.delim(gzfile(file.path(res_dir, "caas_permulation/perm_pos_cycle_caas.tsv.gz")), stringsAsFactors = FALSE,
                  colClasses = c(caas_score = "character"))
stopifnot(all(nul$score_aggregation == "us_plus_gs_mean"))
nul$caas_score <- suppressWarnings(as.numeric(nul$caas_score))
nul <- nul[!is.na(nul$caas_score), ]
nul$cyc <- sub("~H.*$", "", nul$cycle)                      # mirror labels collapse to the base cycle
nul$key <- paste(nul$Gene, nul$Position)
# statistic of the cycle: best side of the position
st <- aggregate(caas_score ~ key + cyc, nul, max)
keys <- sort(unique(st$key)); cycs <- unique(st$cyc)
N <- length(cycs)                                           # every cycle of this null scores some position
S <- matrix(NA_real_, length(keys), N, dimnames = list(keys, cycs))
S[cbind(match(st$key, keys), match(st$cyc, cycs))] <- st$caas_score
det <- !is.na(S)
nd_all <- rowSums(det)
stratum <- cut(nd_all, c(0, 5, 20, 100, Inf), labels = c("nd<=5", "6-20", "21-100", ">100"))
cat(sprintf("%s: %d positions scored by the null, N = %d cycles, %d detections\n", basename(res_dir), length(keys), N, sum(det)))

alphas <- c(0.001, 0.01, 0.05)
pe <- pf <- matrix(1, length(keys), N)                       # p of every (position, cycle) pair; 1 where not scored
hit_bh <- matrix(FALSE, N, 4, dimnames = list(NULL, c("emp_0.05", "emp_0.1", "fact_0.05", "fact_0.1")))
long <- which(det, arr.ind = TRUE)                           # (position, cycle) of every detection
for (cy in seq_len(N)) {
  d <- which(det[, cy])                                      # positions this cycle scores
  s <- S[d, cy]
  # empirical p against the other cycles
  k <- rowSums(S[d, , drop = FALSE] >= s - TIE_TOL, na.rm = TRUE) - 1L
  p1 <- ifelse(s <= TIE_TOL, 1, (k + 1) / N)
  # factorized p: classes and pools fitted on the other cycles
  nd_m <- nd_all - det[, cy]                                 # cycles other than c that score each position
  oth <- long[long[, 2] != cy, , drop = FALSE]
  fit <- .fact_fit(S[oth], nd_m[oth[, 1]])
  p2 <- .fact_p(s, nd_m[d], .fact_assign_obs(fit, nd_m[d]), fit, N - 1L)
  pe[d, cy] <- p1; pf[d, cy] <- p2
  # BH over the null universe, pairs not scored by this cycle at p = 1
  for (j in 1:2) {
    q <- p.adjust(if (j == 1) pe[, cy] else pf[, cy], "BH")
    hit_bh[cy, (j - 1) * 2 + 1:2] <- c(any(q < 0.05), any(q < 0.1))
  }
}

# share of (position, cycle) pairs with p <= alpha, overall and by propensity class
rows <- list()
for (cls in c("all", levels(stratum))) {
  m <- if (cls == "all") rep(TRUE, length(keys)) else stratum == cls
  for (al in alphas) for (nm in c("p.emp", "p.emp_fact")) {
    P <- if (nm == "p.emp") pe else pf
    rate <- mean(P[m, , drop = FALSE] <= al + 1e-12)
    rows[[length(rows) + 1]] <- data.frame(stratum = cls, positions = sum(m), alpha = al, p = nm, rate = rate, ratio = rate / al)
  }
}
out <- do.call(rbind, rows)
write.table(out, paste0(out_prefix, ".rates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
bh <- data.frame(adjustment = c("BH of p.emp", "BH of p.emp_fact"),
                 cycles_with_call_0.05 = c(mean(hit_bh[, 1]), mean(hit_bh[, 3])),
                 cycles_with_call_0.1 = c(mean(hit_bh[, 2]), mean(hit_bh[, 4])), N = N)
write.table(bh, paste0(out_prefix, ".bh.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(format(out, digits = 3), row.names = FALSE); print(bh, row.names = FALSE)
