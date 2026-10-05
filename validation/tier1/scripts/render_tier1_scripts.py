#!/usr/bin/env python3
"""Render the run scripts of a Tier 1 GUI template into a fresh output directory, as the GUI's Generate Scripts does.

The template's results, work and ASR-cache directories are replaced by subdirectories of --outdir, so a rerun never reuses
the work cache or the results of an earlier run. Everything else is the template's: the trait rows, the modules, the
seed and the permulation settings. The scripts are written to --outdir, ready to run.

    python3 validation/tier1/scripts/render_tier1_scripts.py --template gui/templates/tier1_pepc_c4.json \\
        --outdir validation/tier1/output/pepc_unified --evidence-top-n 30
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from gui.generation.render import render_batch, render_single  # noqa: E402
from gui.generation.validate import validate  # noqa: E402
from gui.project_io import load_project  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--template", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--evidence-top-n", type=int, default=0, help="caas_evidence_top_n (0 = off)")
    args = ap.parse_args()

    out = Path(args.outdir).resolve()
    proj = load_project(Path(args.template))
    proj.runtime.results_dir = str(out / "results")
    proj.runtime.work_dir = str(out / "work")
    proj.modules.disambiguation.ct_disambig_asr_cache_dir = str(out / "asr_cache")
    proj.modules.scoring.caas_evidence_top_n = str(args.evidence_top_n)
    errors = validate(proj)
    if errors:
        sys.exit("The template does not validate:\n" + "\n".join(f"  {e}" for e in errors))

    base = proj.runtime.script_base_name.strip() or None
    plural, singular = (base, base) if base else ("phenotypes", "phenotype")
    local = proj.runtime.runtime_type == "local"
    batch_name = f"run_{plural}_local_complete.sh" if local else f"SBATCH_run_{plural}_complete.sh"
    single_name = f"run_{singular}_single_complete.sh"
    out.mkdir(parents=True, exist_ok=True)
    (out / batch_name).write_text(render_batch(proj, postproc_mode="filter", single_runner_filename=single_name))
    (out / single_name).write_text(render_single(proj, postproc_mode="filter"))
    for name in (batch_name, single_name):
        (out / name).chmod(0o755)
    print(f"wrote {out / batch_name} and {out / single_name}")


if __name__ == "__main__":
    main()
