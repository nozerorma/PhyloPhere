"""Render the Tier 1 pipeline entrypoint via the GUI generation recipe.

Same approach as the (demoted) Tier 0 ``build_project.py``: construct the
``ProjectConfig`` the runner GUI uses and render it through
``gui.generation.render`` so the validation runs stay in lockstep with the GUI.

Module set for this stage (MIGUEL: "only fade and caas, not rer"):
    contrast_selection + CAAS/CAAP(discovery,resample,bootstrap,permulation)
    + CT_DISAMBIGUATION(+POSTPROC + ASR robustness) + FADE + SCORING
    OFF: RER, accumulation, enrichment/position-FCS (no bespoke position-GMTs),
         VEP, reporting

Writes ``<repo>/run_phenotype_single.sh`` (the GUI per-run script — must sit at
repo root, it resolves REPO_DIR from its own path) and ``<out>/tier0_env.sh``
(the export preamble the run_*.sh wrappers source, then override per run).
"""

from __future__ import annotations

import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[2]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from gui.generation.render import render_single  # noqa: E402
from gui.models.general import GeneralConfig  # noqa: E402
from gui.models.project import ProjectConfig  # noqa: E402
from gui.models.runtime import PhenotypeRow, RuntimeConfig  # noqa: E402

# reuse the Tier 0 render() helper: identical template-patching for local runs
sys.path.insert(0, str(_REPO / "validation" / ".demoted"))
from tier0.build_project import render  # noqa: E402


def tier1_project(
    *,
    traitname: str,
    trait_type: str = "ordinal",              # "ordinal" | "continuous" | ""
    discrete_method: str = "quintile",        # only used when trait_type != ordinal
    min_divergent_fraction: str = "0.75",
    asr_model: str = "lg",                    # lg | jtt | wag  (ASR-robustness control)
    clade_name: str = "tier1",
    seed: int = 1998,
    perm_pool_size: str = "1000",
    caas_full_perms: str = "200",
    max_tries: str = "200000",
    ct_batch_size: str = "25",
    max_cpus: str = "8",
    max_memory: str = "24.GB",
    max_time: str = "8.h",
) -> ProjectConfig:
    p = ProjectConfig()

    p.general = GeneralConfig(
        project_name="tier1_validation",
        seed=str(seed),
        reporting=False,
        repo_dir=str(_REPO),
        nextflow_plugins_dir=str(Path.home() / ".nextflow" / "plugins"),
    )

    row = PhenotypeRow(trait_class=2, trait=traitname, trait_type=trait_type)
    if trait_type != "ordinal":
        row.discrete_method = discrete_method

    p.runtime = RuntimeConfig(
        resume=False, toy_mode=False, runtime_type="local", batched=False,
        use_tower=False, ali_format="fasta", sp_colname="species",
        clade_name=clade_name, taxon_of_interest="family",
        phenotype_rows=[row],
    )

    p.resources.local_max_cpus = max_cpus
    p.resources.local_max_memory = max_memory
    p.resources.local_max_time = max_time

    # ── CAAS / CAAP ────────────────────────────────────────────────────────
    c = p.modules.caas
    c.enabled = True
    c.ct_tool_discovery = c.ct_tool_resample = c.ct_tool_bootstrap = True
    c.caas_permulation_enrichment = True
    c.caap_mode = True                       # US + GS1..GS4 encodings
    c.patterns = "1,2,3"
    c.perm_strategy = "BM"
    c.perm_pool_size = perm_pool_size
    c.caas_full_perms = caas_full_perms
    c.max_tries = max_tries
    c.min_contrasts = "3"
    c.min_divergent_fraction = min_divergent_fraction
    c.ct_discovery_batch_size = ct_batch_size
    c.ct_bootstrap_batch_size = ct_batch_size

    # ── Disambiguation (+ postproc + ASR robustness) ───────────────────────
    d = p.modules.disambiguation
    d.enabled = True
    d.ct_disambig_asr_mode = "compute"       # no precomputed ASR for these fixtures
    d.ct_disambig_asr_model = asr_model
    d.ct_disambig_convergence_mode = "focal_clade"
    d.asr_robustness = True
    d.gene_filter_mode = "dubious"           # production default; needs gene_ensembl

    # ── FADE (amino-acid in PhyloPhere; runs on the AA alignment) ──────────
    f = p.modules.fade
    f.enabled = True
    f.fade_model = "LG"
    f.fade_method = "Variational-Bayes"

    # ── SCORING ───────────────────────────────────────────────────────────
    p.modules.scoring.enabled = True
    p.modules.scoring.scoring_stress = False

    # ── OFF (this stage) ──────────────────────────────────────────────────
    p.modules.rer.enabled = False
    p.modules.accumulation.enabled = False
    p.modules.vep.enabled = False
    e = p.modules.enrichment
    e.enabled = False
    e.posenrich_enabled = False
    e.scoring_ami = False
    e.scoring_string = False

    return p


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path,
                    default=_REPO / "validation" / "runs" / "tier1" / "_generated")
    ap.add_argument("--trait", default="c4")
    ap.add_argument("--trait-type", default="ordinal")
    ap.add_argument("--mdf", default="0.75")
    ap.add_argument("--asr-model", default="lg")
    ap.add_argument("--clade", default="tier1")
    ap.add_argument("--print", action="store_true")
    a = ap.parse_args()

    proj = tier1_project(
        traitname=a.trait, trait_type=a.trait_type,
        min_divergent_fraction=a.mdf, asr_model=a.asr_model, clade_name=a.clade,
    )
    proj.runtime.alignment_dir = "REPLICATE/align"
    proj.runtime.tree_file = "REPLICATE/tree.nwk"
    proj.runtime.trait_file = "REPLICATE/my_traits.tsv"
    proj.runtime.simple_trait_file = "REPLICATE/my_traits.tsv"
    proj.runtime.results_dir = "REPLICATE/out"
    proj.runtime.work_dir = "REPLICATE/work"
    if a.print:
        print(render_single(proj))
    else:
        paths = render(proj, out_dir=a.out)
        print(f"wrote {paths['single']}")
        print(f"wrote {paths['env']}")
