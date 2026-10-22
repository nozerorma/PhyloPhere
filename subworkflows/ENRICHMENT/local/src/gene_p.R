# gene_p.R — per-gene permulation p of the CAAS gene score (WP6, stage GENE_P).
# PhyloPhere | subworkflows/ENRICHMENT/local/src/
#
# S: genes x N matrix of null gene scores (one column per null cycle); obs: observed score per gene (length = nrow(S)).
# Same pseudo-count rule for both sides: p = (1 + m) / (M + 1), with m = number of comparison cycles whose score is at
# least the gene's score (within TIE_TOL) and M = number of comparison cycles.
#   observed: M = number of non-NA null cycles of the gene, m = #{k: S[g,k] >= obs[g] - TIE_TOL}
#   null cycle j, leave-one-out: M = N - 1, m = #{k != j: S[g,k] >= S[g,j] - TIE_TOL}, so p = (1 + m) / N
# The leave-one-out count uses the row sorted once per gene (O(N log N) per gene); it never returns a p of 0, so the
# -log10 p is finite.
TIE_TOL_GENE_P <- 1e-12

gene_p <- function(S, obs, tie_tol = TIE_TOL_GENE_P) {
  S <- as.matrix(S)
  N <- ncol(S)
  stopifnot(N >= 2, length(obs) == nrow(S))
  obs_m <- rowSums(S >= (obs - tie_tol), na.rm = TRUE)            # NA scores count as no comparison
  p_obs <- (1 + obs_m) / (rowSums(!is.na(S)) + 1)
  p_null <- matrix(NA_real_, nrow(S), N, dimnames = dimnames(S))
  for (g in seq_len(nrow(S))) {
    x <- S[g, ]
    ok <- !is.na(x)
    if (sum(ok) < 2) next
    sx <- sort(x[ok])
    n_ok <- length(sx)
    # number of values >= v in the row, for v = each column's score minus the tolerance
    ge_all <- n_ok - findInterval(x[ok] - tie_tol, sx, left.open = TRUE)
    m_loo <- ge_all - 1L                                          # the column itself is always counted once
    p_null[g, ok] <- (1 + m_loo) / (n_ok)                         # M = n_ok - 1 comparisons: (1 + m) / ((n_ok - 1) + 1)
  }
  list(p_obs = p_obs, p_null = p_null, nlog10_obs = -log10(p_obs), nlog10_null = -log10(p_null))
}
