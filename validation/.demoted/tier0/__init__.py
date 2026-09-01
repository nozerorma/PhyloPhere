"""Tier 0 — simulation-based null calibration and planted-positive recovery.

DEMOTED / UNMAINTAINED (D-DIR-01). Moved to validation/.demoted/; package
imports are broken by the move. See validation/.demoted/README.md.


The generative model here is deliberately *not* the pipeline's inference model
(PhyloPhere ASR / FADE default to LG, single profile per site). Tier 0 simulates
under WAG exchangeabilities with **site-specific Dirichlet profiles**
(CAT/C-series-like heterogeneity), so a passing null validates the method rather
than a shared assumption, and background (non-adaptive) convergence actually
arises the way it does in real data.

Modules:
    model.py          PAML rate-matrix parsing, reversible Q, per-site profile
                      draws, gamma rates
    groups.py         GS1-GS4 amino-acid class tables (vendored from the pipeline)
    trees.py          tree loading / depth-preserving pruning / star + ladder
    pheno.py          make_lambda_foreground — BM-under-lambda latent trait ->
                      independent tail pairs; binary / rate / continuous archetype
                      (one per contrast-selection path); lambda in {0, 0.5, 1}
    simulate.py       evolve a sequence down a tree, inject planted positives
                      (identical_aa + grouped_caap), write FASTA + truth
    replicate.py      one (tree, phenotype) draw -> n_genes alignments + support
                      files + truth.json
    build_project.py  the pipeline entrypoint via the GUI generation recipe
    run_replicates.py stage + run the primate_half x {binary, rate, continuous}
                      x {null, power} x lambda {0, 0.5, 1} gate
"""
