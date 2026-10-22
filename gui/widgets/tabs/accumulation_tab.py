#!/usr/bin/env python3
# accumulation_tab.py — Accumulation module tab.
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
AccumulationTab: field specification of the Accumulation module (conf/accumulation.config).

Declares the ModuleTabSpec and binds it to AccumulationConfig through ModuleTabWidget.

Imported by: gui/widgets/main_window.py
"""

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.modules import AccumulationConfig
from gui.widgets.common.module_tab import ModuleTabWidget
from gui.widgets.common.specs import FieldSpec, ModuleTabSpec, Section

SPEC = ModuleTabSpec(
    title="Accumulation",
    blurb=(
        "Tests, gene by gene, whether CAAS positions accumulate on a gene more than its share of "
        "the position universe predicts (exact hypergeometric test, one p-value per gene and "
        "direction, adjusted by Benjamini-Hochberg over the genes with observed CAAS).\n\n"
        "ℹ️ Universe U: the positions CAAStools tested (in genes that survive post-processing) "
        "that score in at least one cycle of the filtered permulation null that SCORING also "
        "uses (null_positions.parquet, auto-wired when the CAAS permulation ran upstream, "
        "or via --caas_null_dir), plus the observed CAAS positions. It is the same "
        "universe as in Scoring and Position enrichment, and a gene's length enters through its "
        "number of positions in U (PoolPositions in the results). A position detected on both "
        "sides counts once in the 'all' direction.\n\n"
        "Positions are treated as independent draws, which understates the variance of genes with "
        "many positions in U: read the results as a ranking, not as a list of significant genes.\n"
    ),
    disclaimer=(
        "Enrichment's accumulation gene lists need this module's output. Check 'Use "
        "precomputed Accumulation output' on the Precomputed Run tab instead."
    ),
    essential_fields=(
        Section("Accumulation parameters"),
        FieldSpec(name="accumulation_fdr", label="FDR threshold", importance="default"),
    ),
)


class AccumulationTab(ModuleTabWidget):
    def __init__(self, config: AccumulationConfig, parent=None):
        super().__init__(SPEC, config, parent)
