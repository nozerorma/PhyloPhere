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
            name="fade_direction", help="Foreground direction: top = high-extreme species, bottom = low-extreme species, both = run both as separate tests.",
            label="Foreground direction(s)",
            kind="choice",
            choices=("both", "top", "bottom"),
            importance="default",
        ),
        # Sets the background branch definition the Bayes Factor is computed
        # against, which affects test validity like the choice of a control group.
        FieldSpec(
            name="fade_background_scope", help="Background branches: all = every non-foreground branch (HyPhy default); opposite = only the other extreme group, middle species excluded.",
            label="Background scope",
            kind="choice",
            choices=("all", "opposite"),
            importance="default",
        ),
        # Which branches of a foreground clade are tested: every branch of an
        # all-foreground clade (stem included) or the terminal branches only.
        FieldSpec(
            name="fade_internal_nodes", help="all_descendants = label every branch of a foreground clade, stem included (HyPhy default); none = foreground terminal branches only.",
            label="Foreground internal branches",
            kind="choice",
            choices=("all_descendants", "none"),
            importance="default",
        ),
        FieldSpec(
            name="fade_species_file", help="Optional species table (species, contrast group, pair) that replaces automatic contrast selection for FADE only.",
            label="Custom fg/bg species file (optional)",
            kind="path_file",
            importance="optional",
        ),
        Section("Selection parameters"),
        FieldSpec(
            name="fade_model", help="Protein substitution model used by FADE; LG is the default, with WAG, JTT and Dayhoff as alternatives.",
            label="Substitution model",
            kind="choice",
            choices=("LG", "WAG", "JTT", "Dayhoff"),
            importance="default",
        ),
        FieldSpec(name="fade_bf_threshold", help="Bayes Factor cutoff for calling a site under directional selection; default 100.", label="Bayes Factor threshold", importance="default"),
    ),
    advanced_fields=(
        Section("HyPhy inference and MCMC sampling"),
        FieldSpec(
            name="fade_method", help="Inference method: Variational-Bayes (default) or the MCMC samplers Collapsed-Gibbs and Metropolis-Hastings.",
            label="Inference method",
            kind="choice",
            choices=("Variational-Bayes", "Collapsed-Gibbs", "Metropolis-Hastings"),
            importance="default",
        ),
        FieldSpec(name="fade_grid", help="Grid resolution for posterior estimation; default 20.", label="Posterior grid resolution", importance="default"),
        FieldSpec(name="fade_chains", help="Number of MCMC chains; used only by the MCMC inference methods, default 5.", label="MCMC chains", importance="default"),
        FieldSpec(name="fade_chain_length", help="Length of each MCMC chain; used only by the MCMC inference methods, default 2000000.", label="MCMC chain length", importance="default"),
        FieldSpec(name="fade_burn_in", help="MCMC steps discarded from the start of each chain; used only by the MCMC inference methods, default 1000000.", label="MCMC burn-in", importance="default"),
        FieldSpec(name="fade_samples", help="Posterior samples kept from the MCMC run; used only by the MCMC inference methods, default 1000.", label="MCMC samples", importance="default"),
        FieldSpec(name="fade_concentration", help="Dirichlet concentration prior of the FADE model; default 0.5.", label="Dirichlet concentration prior", importance="default"),
        Section("Batching and performance"),
        FieldSpec(name="selection_prep_batch_size", help="Genes per alignment-preparation task; default 500. Smaller values give more, shorter Nextflow tasks.", label="Alignment-prep genes per task", importance="optional"),
        FieldSpec(name="fade_batch_size", help="Genes per FADE Nextflow task; default 200. Smaller values give more, shorter tasks.", label="FADE genes per task", importance="optional"),
        Section("Reporting and precomputed inputs"),
        FieldSpec(name="fade_min_genes_for_heatmap", help="Minimum number of genes needed to draw the report heatmaps; default 2.", label="Min genes for report heatmaps", importance="optional"),
        FieldSpec(name="fade_universe_file", help="Gene list defining the FADE tested universe; blank falls back to the cleaned background.", label="FADE tested-gene universe file", kind="path_file", importance="optional"),
    ),
)


class FadeTab(ModuleTabWidget):
    def __init__(self, config: FadeConfig, parent=None):
        super().__init__(SPEC, config, parent)
