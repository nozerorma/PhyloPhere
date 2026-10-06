#!/usr/bin/env python3
# fade_tab.py — FADE module tab (directional selection, HyPhy).
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
FadeTab: field specification of the FADE module (directional selection, HyPhy).

Off by default (FadeConfig.enabled is False in gui/models/modules.py).

Imported by: gui/widgets/main_window.py
"""

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.modules import FadeConfig
from gui.widgets.common.module_tab import ModuleTabWidget
from gui.widgets.common.specs import FieldSpec, ModuleTabSpec, Section

SPEC = ModuleTabSpec(
    title="FADE",
    blurb=(
        "Directional selection analysis (HyPhy FADE) across the phylogeny. "
        "Choose which extreme direction(s) to test, what counts as the "
        "background branch set, and optionally supply your own foreground/"
        "background species assignment instead of the automatic contrast "
        "selection."
    ),
    disclaimer=(
        "Scoring needs this module's output when it's off. Check 'Use precomputed "
        "FADE output' on the Precomputed Run tab — it's auto-derived per phenotype "
        "from one base path, no per-row entry needed."
    ),
    essential_fields=(
        Section("Direction and background"),
        # Rated "default": the direction (top, bottom or both extremes) sets which
        # hypothesis is tested, so a wrong choice answers a different question
        # without any error.
        FieldSpec(
            name="fade_direction",
            label="Foreground direction(s)",
            kind="choice",
            choices=("both", "top", "bottom"),
            importance="default",
        ),
        # Sets the background branch definition the Bayes Factor is computed
        # against, which affects test validity like the choice of a control group.
        FieldSpec(
            name="fade_background_scope",
            label="Background scope",
            kind="choice",
            choices=("all", "opposite"),
            importance="default",
        ),
        # Which branches of a foreground clade are tested: every branch of an
        # all-foreground clade (stem included) or the terminal branches only.
        FieldSpec(
            name="fade_internal_nodes",
            label="Foreground internal branches",
            kind="choice",
            choices=("all_descendants", "none"),
            importance="default",
        ),
        FieldSpec(
            name="fade_species_file",
            label="Custom fg/bg species file (optional)",
            kind="path_file",
            importance="optional",
        ),
        Section("Selection parameters"),
        FieldSpec(
            name="fade_model",
            label="Substitution model",
            kind="choice",
            choices=("LG", "WAG", "JTT", "Dayhoff"),
            importance="default",
        ),
        FieldSpec(name="fade_bf_threshold", label="Bayes Factor threshold", importance="default"),
    ),
    advanced_fields=(
        Section("HyPhy inference and MCMC sampling"),
        FieldSpec(
            name="fade_method",
            label="Inference method",
            kind="choice",
            choices=("Variational-Bayes", "Collapsed-Gibbs", "Metropolis-Hastings"),
            importance="default",
        ),
        FieldSpec(name="fade_grid", label="Posterior grid resolution", importance="default"),
        FieldSpec(name="fade_chains", label="MCMC chains", importance="default"),
        FieldSpec(name="fade_chain_length", label="MCMC chain length", importance="default"),
        FieldSpec(name="fade_burn_in", label="MCMC burn-in", importance="default"),
        FieldSpec(name="fade_samples", label="MCMC samples", importance="default"),
        FieldSpec(name="fade_concentration", label="Dirichlet concentration prior", importance="default"),
        FieldSpec(
            name="lg_dat_path",
            label="Substitution matrix path (optional)",
            kind="path_file",
            importance="optional",
            help=(
                "Path to substitution matrix (.dat) file. When left blank, "
                "automatically resolves to the matching model file in "
                "subworkflows/SELECTION/local/dat/."
            ),
        ),
        Section("Batching and performance"),
        FieldSpec(name="selection_prep_batch_size", label="Alignment-prep genes per task", importance="optional"),
        FieldSpec(name="fade_batch_size", label="FADE genes per task", importance="optional"),
        Section("Reporting and precomputed inputs"),
        FieldSpec(name="fade_min_genes_for_heatmap", label="Min genes for report heatmaps", importance="optional"),
        FieldSpec(name="fade_universe_file", label="FADE tested-gene universe file", kind="path_file", importance="optional"),
    ),
)


class FadeTab(ModuleTabWidget):
    def __init__(self, config: FadeConfig, parent=None):
        super().__init__(SPEC, config, parent)
