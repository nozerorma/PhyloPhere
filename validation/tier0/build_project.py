"""Build the Tier 0 pipeline entrypoint via the GUI generation recipe.

Rather than hand-write a `nextflow run` invocation, we construct the same
`ProjectConfig` the runner GUI uses and render it through
`gui.generation.render` — so the Tier 0 runs stay in lockstep with whatever the
GUI produces (`params.json` shape, module gating, `${VAR:-default}` fallbacks).

Module set is DECISIONS.md D-T0-D:
    contrast_selection + CAAS(discovery,resample,bootstrap,permulation)
    + CT_DISAMBIGUATION(+POSTPROC+ASR robustness) + RER + FADE + SCORING
    OFF: accumulation, enrichment/STRING/DOMINO, VEP

Outputs:
    <repo>/run_phenotype_single.sh   the per-run script (must sit at repo root:
                                     it resolves REPO_DIR from its own path;
                                     already covered by the repo .gitignore)
    <out>/tier0_env.sh               the `export` preamble the batch wrapper
                                     would set — sourced by run_replicates.py,
                                     which then overrides the per-replicate paths
"""

from __future__ import annotations

import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[2]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from gui.generation.render import render_batch, render_single  # noqa: E402
from gui.models.general import GeneralConfig  # noqa: E402
from gui.models.project import ProjectConfig  # noqa: E402
from gui.models.runtime import PhenotypeRow, RuntimeConfig  # noqa: E402

TRAITNAME = "sim_trait"


def tier0_project(
    *,
    repo_dir: Path,
    alignment_dir: Path,
    tree_file: Path,
    trait_file: Path,
    results_dir: Path,
    work_dir: Path,
    seed: int = 1998,
    perm_pool_size: str = "1000",
    caas_full_perms: str = "200",
    max_tries: str = "200000",
    rer_perm_batches: str = "5",
    rer_perms_per_batch: str = "100",
    top_quantile: str = "0.90",
    bottom_quantile: str = "0.10",
    min_divergent_fraction: str = "0.5",   # 0.5 = production; 1.0 = strict (no conserved pair)
    ct_batch_size: str = "25",
    max_cpus: str = "8",
    max_memory: str = "24.GB",
    max_time: str = "8.h",
    nextflow_plugins_dir: str = "",
) -> ProjectConfig:
    p = ProjectConfig()

    p.general = GeneralConfig(
        project_name="tier0_validation",
        seed=str(seed),
        # reporting OFF (MIGUEL): the dataset/phenotype-exploration Rmds assume
        # real data (PhyloPic silhouette lookups on taxon names, real taxonomy).
        # contrast_selection runs its own dataset_exploration regardless and
        # emits trait_stats.csv, which is FADE's SELECTION_PREP fallback.
        reporting=False,
        repo_dir=str(repo_dir),
        # the template symlinks this into each run's NXF_HOME; must be a real,
        # writable plugins dir (default: the user's ~/.nextflow/plugins).
        nextflow_plugins_dir=nextflow_plugins_dir or str(Path.home() / ".nextflow" / "plugins"),
    )

    p.runtime = RuntimeConfig(
        resume=False,
        toy_mode=False,
        runtime_type="local",
        batched=False,
        use_tower=False,
        work_dir=str(work_dir),
        results_dir=str(results_dir),
        alignment_dir=str(alignment_dir),
        ali_format="fasta",
        tree_file=str(tree_file),
        trait_file=str(trait_file),
        simple_trait_file=str(trait_file),
        sp_colname="species",
        clade_name="tier0",
        taxon_of_interest="family",
        phenotype_rows=[
            PhenotypeRow(trait_class=2, trait=TRAITNAME, discrete_method="parameterized")
        ],
    )

    p.resources.local_max_cpus = max_cpus
    p.resources.local_max_memory = max_memory
    p.resources.local_max_time = max_time

    # ── CAAS ────────────────────────────────────────────────────────────────
    c = p.modules.caas
    c.enabled = True
    c.ct_tool_discovery = c.ct_tool_resample = c.ct_tool_bootstrap = True
    c.caas_permulation_enrichment = True
    c.caap_mode = True
    c.patterns = "1,2,3"
    c.perm_strategy = "BM"            # D-T0-B
    c.perm_pool_size = perm_pool_size
    c.caas_full_perms = caas_full_perms
    c.max_tries = max_tries
    c.top_quantile = top_quantile
    c.bottom_quantile = bottom_quantile
    c.min_contrasts = "3"
    c.min_divergent_fraction = min_divergent_fraction
    c.ct_discovery_batch_size = ct_batch_size
    c.ct_bootstrap_batch_size = ct_batch_size

    # ── Disambiguation (+ postproc + ASR robustness) ────────────────────────
    d = p.modules.disambiguation
    d.enabled = True
    d.ct_disambig_asr_mode = "compute"      # no precomputed ASR for simulated data
    d.ct_disambig_asr_model = "lg"
    d.ct_disambig_convergence_mode = "focal_clade"
    d.asr_robustness = True
    # "none": the dubious/extreme gene filters need a --gene_ensembl_file mapping
    # synthetic gene ids to genomic coordinates, which does not exist for
    # simulated genes. Tier 0 keeps the full postproc gene pool.
    d.gene_filter_mode = "none"

    # ── RER ────────────────────────────────────────────────────────────────
    r = p.modules.rer
    r.enabled = True
    r.rer_tool_build_trait = r.rer_tool_build_tree = r.rer_tool_build_matrix = True
    r.rer_tool_continuous = True
    r.rer_trait_mode = "auto"       # echo -> binary (2 unique vals); bodysize -> continuous
    # NOT the GUI default "ha_logit" (needs n_trait/c_trait count columns; stop()s
    # without them). MIGUEL (D-T0-F): "auto". echo takes the binary path (no
    # transform); the BM bodysize trait is ~normal and real-valued so auto -> raw
    # (log10 is the natural allometric transform but our BM trait is centred at 0,
    # so log10 is N/A unless we simulate exp(BM) instead — noted, not done).
    r.rer_transform = "auto"
    r.rer_perm_batches = rer_perm_batches
    r.rer_perms_per_batch = rer_perms_per_batch
    r.rer_minsp = "10"

    # ── FADE ───────────────────────────────────────────────────────────────
    f = p.modules.fade
    f.enabled = True
    f.fade_mode = "all"
    f.fade_method = "Variational-Bayes"
    f.fade_model = "LG"
    f.fade_batch_size = "50"
    f.selection_prep_batch_size = "100"

    # ── SCORING on ─────────────────────────────────────────────────────────
    p.modules.scoring.enabled = True
    p.modules.scoring.scoring_stress = False

    # ── OFF (D-T0-D) ───────────────────────────────────────────────────────
    p.modules.accumulation.enabled = False
    p.modules.vep.enabled = False
    e = p.modules.enrichment
    e.enabled = False
    e.posenrich_enabled = False
    e.scoring_ami = False
    e.scoring_string = False

    return p


# --------------------------------------------------------------------------- #
# render
# --------------------------------------------------------------------------- #
_ENV_START = "conda activate phylophere"
_ENV_END = "# --- PHENOTYPE CATALOGUE ---"

_PORTABLE_ACTIVATE = """\
# [tier0] portable env activation (local runs) — see build_project.render()
set +u
if command -v micromamba >/dev/null 2>&1; then
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate phylophere
elif [ -n "${CONDA_EXE:-}" ] || command -v conda >/dev/null 2>&1; then
    __base="$(conda info --base 2>/dev/null || dirname "$(dirname "${CONDA_EXE:-}")")"
    [ -f "$__base/etc/profile.d/conda.sh" ] && . "$__base/etc/profile.d/conda.sh"
    conda activate phylophere
fi
set -u"""


def render(project: ProjectConfig, *, out_dir: Path) -> dict:
    """Write `run_phenotype_single.sh` to the repo root and `tier0_env.sh` (the
    batch wrapper's `export` preamble) to `out_dir`. Returns the paths."""
    repo = Path(project.general.repo_dir)
    single = render_single(project)
    if project.runtime.runtime_type == "local":
        # The template hardcodes `source ~/.bashrc; conda deactivate; conda
        # activate phylophere`, which assumes `conda init` was run. On a
        # micromamba box (conda is a shell function, ~/.bashrc early-returns for
        # non-interactive shells) `conda deactivate` errors under `set -e`.
        # Swap in a portable activation for local Tier 0 runs only.
        marker = "source ~/.bashrc\nconda deactivate\nconda activate phylophere"
        if marker not in single:
            raise RuntimeError("env-activation marker not found in rendered run_single — template changed")
        single = single.replace(marker, _PORTABLE_ACTIVATE)
    single_path = repo / "run_phenotype_single.sh"
    single_path.write_text(single)
    single_path.chmod(0o755)

    batch = render_batch(project, single_runner_filename="run_phenotype_single.sh")
    lines = batch.splitlines()
    try:
        i0 = next(i for i, ln in enumerate(lines) if ln.strip() == _ENV_START) + 1
        i1 = next(i for i, ln in enumerate(lines) if ln.strip() == _ENV_END)
    except StopIteration:  # template changed — fail loud rather than ship half an env
        raise RuntimeError(
            "could not slice the env preamble out of the rendered batch script; "
            "the sbatch_array template markers moved"
        )
    out_dir.mkdir(parents=True, exist_ok=True)
    env_path = out_dir / "tier0_env.sh"
    env_path.write_text(
        "# Generated from gui.generation.render_batch — the export preamble only.\n"
        "# run_replicates.py sources this, then overrides ALI_DIR / TREE_FILE /\n"
        "# SIMPLE_TRAIT_FILE / TOP_QUANTILE / BOTTOM_QUANTILE / SEED / CAAS_OUTBASE\n"
        "# / WORK_BASE per replicate.\n"
        + "\n".join(lines[i0:i1]) + "\n"
    )
    return {"single": single_path, "env": env_path}


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=_REPO / "validation" / "tier0" / "generated")
    ap.add_argument("--print", action="store_true", help="print run_phenotype_single.sh to stdout")
    a = ap.parse_args()
    proj = tier0_project(
        repo_dir=_REPO,
        alignment_dir=Path("REPLICATE/align"),
        tree_file=Path("REPLICATE/tree.nwk"),
        trait_file=Path("REPLICATE/my_traits.tsv"),
        results_dir=Path("REPLICATE/out"),
        work_dir=Path("REPLICATE/work"),
    )
    if a.print:
        print(render_single(proj))
    else:
        paths = render(proj, out_dir=a.out)
        print(f"wrote {paths['single']}")
        print(f"wrote {paths['env']}")
