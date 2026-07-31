#!/usr/bin/env python3
# fade_tab.py — FADE module tab (directional selection, HyPhy).
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""Off by default (RUN_FADE=false in the reference scripts)."""

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.modules import FadeConfig
from gui.widgets.common.module_tab import ModuleTabWidget
from gui.widgets.common.specs import FieldSpec, ModuleTabSpec

SPEC = ModuleTabSpec(
    title="FADE",
    blurb=(
        "Directional selection analysis (HyPhy FADE) across the phylogeny, run on "
        "either every gene or just the CAAS hit set (fade_mode)."
    ),
    disclaimer=(
        "Scoring needs this module's output when it's off. Check 'Use precomputed "
        "FADE output' on the Precomputed Run tab — it's auto-derived per phenotype "
        "from one base path, no per-row entry needed."
    ),
    essential_fields=(
        FieldSpec(name="fade_mode", label="Mode", kind="choice", choices=("all", "gene_set")),
    ),
    advanced_fields=(
        FieldSpec(name="fade_postproc_top", label="Gene-set top genes (standalone gene_set mode)", kind="path_file"),
        FieldSpec(name="fade_postproc_bottom", label="Gene-set bottom genes (standalone gene_set mode)", kind="path_file"),
        FieldSpec(name="selection_prep_batch_size", label="Alignment-prep genes per task"),
        FieldSpec(name="fade_batch_size", label="FADE genes per task"),
        FieldSpec(name="fade_bf_threshold", label="Bayes Factor threshold"),
        FieldSpec(name="fade_model", label="Substitution model"),
        FieldSpec(name="lg_dat_path", label="LG substitution matrix path", kind="path_file"),
        FieldSpec(
            name="fade_method",
            label="Inference method",
            kind="choice",
            choices=("Variational-Bayes", "Collapsed-Gibbs", "Metropolis-Hastings"),
        ),
        FieldSpec(name="fade_grid", label="Posterior grid resolution"),
        FieldSpec(name="fade_chains", label="MCMC chains"),
        FieldSpec(name="fade_chain_length", label="MCMC chain length"),
        FieldSpec(name="fade_burn_in", label="MCMC burn-in"),
        FieldSpec(name="fade_samples", label="MCMC samples"),
        FieldSpec(name="fade_concentration", label="Dirichlet concentration prior"),
        FieldSpec(name="fade_min_genes_for_heatmap", label="Min genes for report heatmaps"),
        FieldSpec(name="fade_universe_file", label="FADE tested-gene universe file", kind="path_file"),
    ),
)


class FadeTab(ModuleTabWidget):
    def __init__(self, config: FadeConfig, parent=None):
        super().__init__(SPEC, config, parent)
