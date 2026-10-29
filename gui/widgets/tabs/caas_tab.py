#!/usr/bin/env python3
# caas_tab.py — CAAS / CT module tab (contrast selection: discovery, resample).
# PhyloPhere | gui/widgets/tabs/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
CaasTab: field specification of the CAAS / contrast-selection module (conf/ct.config).

The discovery and resample sub-steps are two checkboxes on the config
(ct_tool_discovery, ct_tool_resample; see gui/models/modules.py) that jointly build
the comma-separated --ct_tool value (gui/generation/context.py). ModuleTabWidget
has no field kind for that, so CaasTab adds them directly in __init__.

--contrast_selection has no checkbox of its own: run_single.sh.j2 turns it on
whenever CAAS or Disambiguation is enabled, because contrast selection is the only
producer of the foreground/background trait file that both need.

Imported by: gui/widgets/main_window.py
"""

# ── Third-party ───────────────────────────────────────────────────────────────
from PySide6.QtWidgets import QCheckBox, QLabel

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.modules import CaasConfig
from gui.widgets.common.module_tab import ModuleTabWidget
from gui.widgets.common.specs import FieldSpec, ModuleTabSpec, Section

# Texts of the two --ct_tool checkboxes (translated in gui/i18n.py, keyed by the English text).
_CT_TOOL_DISCOVERY_TIP = (
    "Replays the real labeling (b_0) and the permulated labelings over the alignments, in one pass of the "
    "permulation core; the discovery of the real labeling is its b_0 slice. It needs the labelings that "
    "'resample' harvests. With 0 permulation cycles (--caas_full_perms 0) it replays only the real labeling."
)
_CT_TOOL_RESAMPLE_TIP = (
    "Harvests the permulated labelings (the null) from the tree and the trait. It reads no alignments and "
    "can run alone, for example to inspect the null before replaying it."
)
_CT_TOOL_NOTE = (
    "discovery needs the labelings of resample: check both, or reuse a previous harvest from the "
    "Precomputed Run tab."
)

SPEC = ModuleTabSpec(
    title="CAAS / Contrast Selection",
    blurb=(
        "Runs CAAStools' discovery and resample steps to find convergent "
        "amino-acid substitutions (CAAS) associated with the phenotype."
    ),
    disclaimer=(
        "Disambiguation and Accumulation need this module's output. Check Discovery/"
        "Resample on the Precomputed Run tab to feed them precomputed "
        "results instead."
    ),
    essential_fields=(
        Section("Discovery params"),
        FieldSpec(
            name="patterns",
            label="CT patterns",
            kind="multichoice",
            choices=("1", "2", "3", "4"),
            importance="default",
            help="Pattern classes to keep: 1 one-to-one, 2 one foreground vs many background, 3 many vs one, 4 many-to-many. Default 1,2,3.",
        ),
        FieldSpec(
            name="min_divergent_fraction",
            label="Min divergent fraction",
            kind="choice",
            choices=("0.5", "0.75", "1.0"),
            editable=True,
            importance="default",
            help="Fraction (0-1) of contrast pairs that must show divergent residues; 1.0 allows no conserved pair (strict); default 0.5.",
        ),
        FieldSpec(
            name="min_non_overlapping",
            label="Min non-overlapping species",
            kind="choice",
            choices=("2", "3", "4", "5"),
            editable=True,
            importance="default",
            help="Minimum number of non-overlapping (exclusive) species required on FG or BG to call a CAAS (default: 2).",
        ),
        FieldSpec(
            name="strict_conservation_output",
            label="Also write the strict tables (no conserved pair)",
            kind="bool",
            importance="default",
            help=(
                "With a min divergent fraction below 1, also writes the observed tables and null of hypotheses without conserved pairs under <outdir>/strict/."
            ),
        ),
        FieldSpec(
            name="multi_hypothesis",
            label="Multi-hypothesis mode",
            kind="bool",
            importance="default",
            help="Also tests alternative Dunn-independent contrast hypotheses (H1..Hn) for each contrast, not only the canonical one; default on.",
        ),
        FieldSpec(
            name="max_fop",
            label="Max FOP hypotheses (H1..Hn)",
            placeholder="alternative Dunn-independent hypotheses per contrast (observed + null); default 100",
            importance="default",
            help="Cap on alternative Dunn-independent hypotheses (H1..Hn) per contrast, observed and null; default 100.",
        ),
        Section("Resample / permulation params"),
        FieldSpec(
            name="chunk_size",
            label="Resampled groups per output file",
            importance="optional",
            help="Number of permulation cycles written per resample_*.tab output file; default 500.",
        ),
        FieldSpec(
            name="resample_use_n",
            label="Use sample size counts (n/c)",
            kind="bool",
            importance="default",
            help="For count traits with sample sizes (N/C), filters permulations by Jeffreys confidence intervals when those sizes are available; default on.",
        ),
        FieldSpec(
            name="perm_strategy",
            label="Permutation strategy",
            kind="choice",
            choices=("auto", "OU", "BM"),
            importance="default",
            help="Evolutionary model used to simulate the null trait: auto picks BM or OU by AICc; OU and BM force one model; default auto.",
        ),
        FieldSpec(
            name="perm_match_pss",
            label="Match observed PSS profile in the null",
            kind="bool",
            importance="default",
            help="Builds each permulation's canonical pairs so their PSS follows the observed PSS profile; not applied to count traits; default on.",
        ),
        FieldSpec(
            name="perm_match_pss_tol",
            label="PSS matching tolerance",
            placeholder="relative tolerance of each pair; default 0.25",
            importance="optional",
            help="Tolerance of that PSS matching: each null pair's PSS must lie within a factor of (1 + tol) of its observed PSS; default 0.25.",
        ),
        FieldSpec(
            name="caas_full_perms",
            label="Permulations",
            placeholder="accepted permulations to harvest AND replay for the CAAS FCS null; 0 replays only the real labeling",
            importance="default",
            help="Accepted permulations harvested and replayed as the CAAS null; 0 replays only the real labeling, with no null; default 1000.",
        ),
        # Rated "optional" although borderline: it is a draw-budget cap (raised 50%
        # up to twice if the pool falls short) and does not redefine the null, but
        # a very low value can truncate the permulation pool.
        FieldSpec(
            name="max_tries",
            label="Max permulation tries",
            placeholder="draw budget; raised 50% up to twice if the pool falls short",
            importance="optional",
            help="Draw budget for harvesting permulations; raised by 50% up to twice if the pool falls short, then the run fails; default 1000000.",
        ),
    ),
    advanced_fields=(
        Section("Missingness parameters in discovery mode and otherwise"),
        FieldSpec(
            name="caas_config_path",
            label="CAAS config file (auto-derived if empty)",
            kind="path_file",
            placeholder="auto-derived from contrast_selection.nf when empty",
            importance="optional",
            help="CAAStools configuration file; leave empty to derive it automatically from the contrast-selection step.",
        ),
        FieldSpec(
            name="max_bg_gaps_fraction",
            label="Max background gaps fraction",
            importance="default",
            help="Maximum fraction (0-1) of pairs allowed to have a gap in the background; 0.0 allows none; default 0.0.",
        ),
        FieldSpec(
            name="max_fg_gaps_fraction",
            label="Max foreground gaps fraction",
            importance="default",
            help="Maximum fraction (0-1) of pairs allowed to have a gap in the foreground; 0.0 allows none; default 0.0.",
        ),
        FieldSpec(
            name="max_gaps_fraction",
            label="Max any-gaps fraction",
            importance="default",
            help="Maximum fraction (0-1) of pairs allowed to have a gap in either group; 0.0 allows none; default 0.0.",
        ),
        FieldSpec(
            name="max_bg_miss_fraction",
            label="Max background missing fraction",
            importance="default",
            help="Maximum fraction (0-1) of pairs allowed to have missing data in the background; default 0.0.",
        ),
        FieldSpec(
            name="max_fg_miss_fraction",
            label="Max foreground missing fraction",
            importance="default",
            help="Maximum fraction (0-1) of pairs allowed to have missing data in the foreground; default 0.0.",
        ),
        FieldSpec(
            name="max_miss_fraction",
            label="Max any-missing fraction",
            importance="default",
            help="Maximum fraction (0-1) of pairs allowed to have missing data in either group; default 0.0.",
        ),
        FieldSpec(
            name="miss_pair",
            label="Enforce missing pairs",
            kind="bool",
            importance="default",
            help="Requires missing or gapped species to come from matching pairs; applies when the max thresholds are equal; default on.",
        ),
        Section("Batching logic (performance)"),
        FieldSpec(
            name="ct_core_batch_size",
            label="Permulation-null genes per task (replay + ASR replay)",
            importance="optional",
            help="Genes per permulation-core task (replay and ASR replay); 1 = one task per gene; default 20.",
        ),
        FieldSpec(
            name="caas_cycle_block_size",
            label="Permulation-null cycles per block",
            importance="optional",
            help="Permulation cycles per block; each task replays one gene batch against one block, so memory does not grow with caas_full_perms; default 500.",
        ),
        Section("Permulation null"),
        FieldSpec(
            name="caas_export_single",
            label="Export single CAAS null RDS",
            kind="bool",
            importance="optional",
            help="Also writes a single compatibility caas_perms.rds of the null; the blocks and roster remain the canonical outputs; default off.",
        ),
        FieldSpec(
            name="caas_perms_postproc",
            label="Apply post-processing filters to the permulation null",
            kind="bool",
            importance="default",
            help="Applies the observed cluster and gene filters to the permulation null, so the null follows filtered discovery; needs a gene Ensembl file; default on.",
        ),
        FieldSpec(
            name="caas_permulation_enrichment",
            label="Use CAAS permulations for the FCS Wilcoxon null",
            kind="bool",
            importance="default",
            help="Uses the CAAS permulation null as the null of the FCS Wilcoxon enrichment test; default on.",
        ),
        Section("Publishing norms (debug)"),
        FieldSpec(
            name="publish_intermediates",
            label="Publish intermediate files",
            kind="bool",
            importance="optional",
            help="Publishes intermediate discovery, background and resample files for debugging; default off.",
        ),
        Section("Contrast-selection tuning (conf/common.config)"),
        FieldSpec(
            name="pss_top_pct",
            label="PSS top percentile candidate gate",
            importance="default",
            help="For continuous traits, only the top fraction (0-1) of pairs by PSS qualifies as candidate; 0.05 keeps the top 5 percent; default 0.05.",
        ),
        FieldSpec(
            name="pair_center_gate",
            label="Candidate pairs must straddle the trait centre",
            kind="choice",
            choices=("none", "median", "mean"),
            importance="default",
            help="none = no constraint; median or mean = each candidate pair must have one member on each side of that trait center; default none.",
        ),
        FieldSpec(
            name="max_contrasts",
            label="Max contrasts (0 = dynamic)",
            importance="default",
            help="Maximum number of contrasts to select; 0 = dynamic discovery until the Dunn index falls below 1; default 0.",
        ),
        FieldSpec(
            name="min_contrasts",
            label="Minimum foreground contrasts",
            importance="default",
            help="Minimum number of foreground contrasts required to proceed; default 3.",
        ),
    ),
)


class CaasTab(ModuleTabWidget):
    def __init__(self, config: CaasConfig, parent=None):
        super().__init__(SPEC, config, parent)

        # discovery and resample checkboxes jointly build --ct_tool; inserted at
        # the top of the essential fields.
        self.ct_tool_discovery = QCheckBox("discovery")
        self.ct_tool_discovery.setChecked(config.ct_tool_discovery)
        self.ct_tool_discovery.toggled.connect(self._on_ct_tool_discovery)

        self.ct_tool_resample = QCheckBox("resample")
        self.ct_tool_resample.setChecked(config.ct_tool_resample)
        self.ct_tool_resample.toggled.connect(self._on_ct_tool_resample)

        self._ct_tool_label = QLabel("CT tools (--ct_tool)")
        self._ct_tool_note = QLabel(_CT_TOOL_NOTE)
        self._ct_tool_note.setWordWrap(True)
        self._essential_form.insertRow(0, self._ct_tool_label, self.ct_tool_discovery)
        self._essential_form.insertRow(1, "", self.ct_tool_resample)
        self._essential_form.insertRow(2, "", self._ct_tool_note)
        self.ct_tool_discovery.setToolTip(_CT_TOOL_DISCOVERY_TIP)
        self.ct_tool_resample.setToolTip(_CT_TOOL_RESAMPLE_TIP)

    def retranslate(self, lang: str = "en") -> None:
        super().retranslate(lang)
        from gui.i18n import tr
        self._ct_tool_label.setText(tr("CT tools (--ct_tool)", lang))
        self.ct_tool_discovery.setText(tr("discovery", lang))
        self.ct_tool_resample.setText(tr("resample", lang))
        self._ct_tool_note.setText(tr(_CT_TOOL_NOTE, lang))
        self.ct_tool_discovery.setToolTip(tr(_CT_TOOL_DISCOVERY_TIP, lang))
        self.ct_tool_resample.setToolTip(tr(_CT_TOOL_RESAMPLE_TIP, lang))

    def _on_ct_tool_discovery(self, value: bool) -> None:
        self._config.ct_tool_discovery = value
        self.changed.emit()

    def _on_ct_tool_resample(self, value: bool) -> None:
        self._config.ct_tool_resample = value
        self.changed.emit()
