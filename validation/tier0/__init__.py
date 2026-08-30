"""Tier 0 — simulation-based null calibration and planted-positive recovery.

The generative model here is deliberately *not* the pipeline's inference model
(PhyloPhere ASR / FADE default to LG, single profile per site). Tier 0 simulates
under WAG exchangeabilities with **site-specific Dirichlet profiles**
(CAT/C-series-like heterogeneity), so a passing null validates the method rather
than a shared assumption, and background (non-adaptive) convergence actually
arises the way it does in real data.

Modules:
    model.py          PAML rate-matrix parsing, reversible Q construction,
                      per-site profile draws, discrete/continuous gamma rates
    trees.py          tree loading / depth-preserving pruning / pathological
                      trees (star, ladder) / foreground-transition sampling
    simulate.py       evolve a sequence down a tree, inject planted positives
                      (profile-shift and identical-AA convergence), write
                      FASTA + phenotype.tsv + truth.json
    make_replicates.py driver that emits the `null` and `power` replicate sets
"""
