# palettes.R — Shared colour palettes and ggplot2 theme for PhyloPhere reports.
# PhyloPhere | subworkflows/TRAIT_ANALYSIS/local/obj/
# Sourced by: 1.Dataset_exploration.Rmd, 2.Phenotype_exploration.Rmd, SCORING_report.Rmd
# =============================================================================

# ── Colour palettes ───────────────────────────────────────────────────────────

# Directional: top (positive phenotype extreme) vs bottom (negative extreme).
PALETTE_DIRECTION <- c(
    top    = "#E63946",   # warm red  — phenotype acceleration
    bottom = "#457B9D",   # steel blue — phenotype deceleration
    all    = "#6D6875"    # muted purple — direction-pooled
)

# Ordinal significance tiers used across all report p-value visualisations.
PALETTE_SIGNIFICANCE <- c(
    "p < 0.01"  = "#E63946",
    "p < 0.05"  = "#F4A261",
    "p < 0.10"  = "#A8DADC",
    "n.s."      = "#CCCCCC"
)

# Sequential for continuous scores (low → high).
PALETTE_SCORE <- c("#F1FAEE", "#A8DADC", "#457B9D", "#1D3557")


# ── ggplot2 theme ─────────────────────────────────────────────────────────────

# Applied via `+ theme_phylo()` in all ggplot calls within pipeline reports.
# Keeps visual language consistent across Rmd outputs without per-plot overrides.
theme_phylo <- function(base_size = 12) {
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
        panel.grid.minor  = ggplot2::element_blank(),
        strip.text        = ggplot2::element_text(face = "bold"),
        legend.position   = "bottom",
        plot.title        = ggplot2::element_text(face = "bold", size = base_size + 2),
        plot.subtitle     = ggplot2::element_text(colour = "grey40"),
    )
}


# ── Helper: significance tier label ──────────────────────────────────────────

sig_tier <- function(pval) {
    dplyr::case_when(
        pval < 0.01 ~ "p < 0.01",
        pval < 0.05 ~ "p < 0.05",
        pval < 0.10 ~ "p < 0.10",
        TRUE        ~ "n.s."
    )
}
