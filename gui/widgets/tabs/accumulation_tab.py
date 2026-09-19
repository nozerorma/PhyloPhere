#!/usr/bin/env python3
# accumulation_tab.py — Accumulation module tab.
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.modules import AccumulationConfig
from gui.widgets.common.module_tab import ModuleTabWidget
from gui.widgets.common.specs import FieldSpec, ModuleTabSpec, Section

SPEC = ModuleTabSpec(
    title="Accumulation",
    blurb=(
        "Tests whether CAAS hits accumulate on specific genes more than expected under "
        "a randomized background, using per-site Valdar-variability conservation weighting.\n\n"
        "ℹ️ Randomization type:\n"
        "• naive: uniform random placement across all eligible (unmasked) alignment "
        "positions — the simplest null, ignores conservation entirely.\n"
        "• cons_decile: bins positions into conservation deciles (from the Valdar "
        "variability score) and randomizes within each gene's own decile profile — "
        "controls for genes differing in overall conservation. Default.\n"
        "• permulation: reuses a prior CAAS_PERMULATION null (via --caas_pos_detail_file, "
        "auto-wired when that module ran upstream) — holds the phenotype-tree confound "
        "but is NOT conservation-decile matched. Opt-in."
    ),
    disclaimer=(
        "Enrichment's accumulation gene lists need this module's output. Check 'Use "
        "precomputed Accumulation output' on the Precomputed Run tab instead."
    ),
    essential_fields=(
        Section("Randomization and burden parameters"),
        FieldSpec(
            name="accumulation_randomization_type",
            label="Randomization type",
            kind="choice",
            choices=("cons_decile", "naive", "permulation"),
        ),
        FieldSpec(name="accumulation_n_randomizations", label="Randomizations"),
        FieldSpec(name="accumulation_fdr", label="FDR threshold"),
    ),
    advanced_fields=(
        Section("Entropy and standalone inputs"),
        FieldSpec(
            name="accumulation_entropy_dir",
            label="Entropy directory (optional)",
            kind="path_dir",
            help=(
                "Only used by 'cons_decile' randomization. Leave blank to auto-generate "
                "per-gene Valdar variability files directly from the alignment (a "
                "verbatim port of the ortholog_characterizator scoring, requires "
                "--tax_id to be set — see the Runtime tab). If left blank with no "
                "tax_id available, falls back to raw majority-residue conservation "
                "computed straight from the alignment (no external file needed either "
                "way, just a coarser conservation measure)."
            ),
        ),
    ),
)


class AccumulationTab(ModuleTabWidget):
    def __init__(self, config: AccumulationConfig, parent=None):
        super().__init__(SPEC, config, parent)
