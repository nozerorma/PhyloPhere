"""Stage (and optionally run) the Tier 0 replicate set through the pipeline.

    python -m validation.tier0.run_replicates \
        --out validation/runs/tier0 \
        --archetypes echo,bodysize --sets null,power \
        --n-replicates 20 --n-genes 40 [--trees primate,mammal] [--run]

Per replicate:
  1. prune/scale the tree (validation/tier0/grid.json `trees` block)
  2. `replicate.build_replicate` → <repdir>/{align/, my_traits.tsv, tree.nwk,
     truth.json}
  3. render the GUI project (build_project) once, then write <repdir>/run.sh —
     sources tier0_env.sh, exports the per-replicate ALI_DIR / TREE_FILE /
     SIMPLE_TRAIT_FILE / TOP_QUANTILE / BOTTOM_QUANTILE / SEED / CAAS_OUTBASE /
     WORK_BASE, then calls <repo>/run_phenotype_single.sh 2 sim_trait "" … parameterized

Default is **stage only** — it writes run.sh + run_all.sh + manifest.jsonl and
stops, so the exact pipeline invocations can be reviewed. `--run` executes them.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

_REPO = Path(__file__).resolve().parents[2]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from validation.tier0 import trees as T  # noqa: E402
from validation.tier0.build_project import TRAITNAME, render, tier0_project  # noqa: E402
from validation.tier0.replicate import ReplicateConfig, build_replicate  # noqa: E402

_GRID = Path(__file__).parent / "grid.json"


def _load_tree(spec: dict, rng: np.random.Generator) -> T.PhyloTree:
    if spec.get("builtin") == "star":
        return T.star_tree(spec["k"], spec.get("branch_length", 0.3))
    if spec.get("builtin") == "ladder":
        return T.ladder_tree(spec["k"], spec.get("branch_length", 0.15))
    path = _REPO / spec["path"]
    if not path.exists():
        raise FileNotFoundError(f"{path} — fetch Tier 0 trees (validation/fixtures/tier0/README.md)")
    tree = T.prune_depth_preserving(path, spec["k"], rng)
    if spec.get("scale", 1.0) != 1.0:
        tree = T.scale_branches(tree, spec["scale"])
    return tree


_RUN_SH = """\
#!/bin/bash
# Tier 0 replicate: {rep} | archetype={arch} set={set} tree={tree} seed={seed}
set -Eeuo pipefail
cd "{repo}"
set -a
source "{env_sh}"
ALI_DIR="{align}"
TREE_FILE="{tree_nwk}"
TRAIT_FILE="{traits}"
SIMPLE_TRAIT_FILE="{traits}"
TOP_QUANTILE="{tq}"
BOTTOM_QUANTILE="{bq}"
SEED="{seed}"
CAAS_OUTBASE="{outbase}"
WORK_BASE="{workbase}"
set +a
mkdir -p "$CAAS_OUTBASE" "$WORK_BASE"
exec bash "{repo}/run_phenotype_single.sh" 2 {trait} "" "" "" "" "" parameterized
"""


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="validation.tier0.run_replicates")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--archetypes", default="echo,bodysize")
    ap.add_argument("--sets", default="null,power", help="null = no planted signal; power = planted")
    ap.add_argument("--trees", default="primate", help="comma list of keys in grid.json trees")
    ap.add_argument("--n-replicates", type=int, default=20)
    ap.add_argument("--n-genes", type=int, default=40)
    ap.add_argument("--frac-planted-genes", type=float, default=0.25)
    ap.add_argument("--n-sites", type=int, default=400)
    ap.add_argument("--n-transitions", type=int, default=4)
    ap.add_argument("--seed0", type=int, default=70000)
    ap.add_argument("--run", action="store_true", help="execute the staged run.sh scripts")
    ap.add_argument("--jobs", type=int, default=1, help="parallel replicates when --run")
    a = ap.parse_args(argv)

    grid = json.loads(_GRID.read_text())
    tree_keys = a.trees.split(",")
    archetypes = a.archetypes.split(",")
    sets = a.sets.split(",")
    a.out.mkdir(parents=True, exist_ok=True)

    # render the GUI project once (paths are per-replicate, overridden in run.sh)
    proj = tier0_project(
        repo_dir=_REPO,
        alignment_dir=Path("REPLICATE/align"), tree_file=Path("REPLICATE/tree.nwk"),
        trait_file=Path("REPLICATE/my_traits.tsv"),
        results_dir=Path("REPLICATE/out"), work_dir=Path("REPLICATE/work"),
    )
    paths = render(proj, out_dir=a.out / "generated")
    env_sh = paths["env"].resolve()

    rows: list[dict] = []
    run_scripts: list[Path] = []
    i = 0
    for arch in archetypes:
        for setname in sets:
            for tkey in tree_keys:
                for k in range(a.n_replicates):
                    i += 1
                    seed = a.seed0 + i
                    rng = np.random.default_rng(seed)
                    tree = _load_tree(grid["trees"][tkey], rng)
                    repdir = a.out / f"{arch}_{setname}_{tkey}" / f"rep{k:03d}"
                    rcfg = ReplicateConfig(
                        archetype=arch, n_genes=a.n_genes,
                        frac_planted_genes=a.frac_planted_genes,
                        traitname=TRAITNAME, n_sites=a.n_sites,
                        n_transitions=a.n_transitions,
                        planted=(setname == "power"),
                    )
                    row = build_replicate(tree, rcfg, seed, repdir)
                    row.update(set=setname, tree_key=tkey,
                               dir=str(repdir.relative_to(a.out)))
                    rows.append(row)

                    run_sh = repdir / "run.sh"
                    run_sh.write_text(_RUN_SH.format(
                        rep=row["dir"], arch=arch, set=setname, tree=tkey, seed=seed,
                        repo=_REPO, env_sh=env_sh,
                        align=(repdir / "align").resolve(),
                        tree_nwk=(repdir / "tree.nwk").resolve(),
                        traits=(repdir / "my_traits.tsv").resolve(),
                        tq=row["top_quantile"], bq=row["bottom_quantile"],
                        outbase=(repdir / "out").resolve(),
                        workbase=(repdir / "work").resolve(),
                        trait=TRAITNAME,
                    ))
                    run_sh.chmod(0o755)
                    run_scripts.append(run_sh)

    (a.out / "manifest.jsonl").write_text("".join(json.dumps(r) + "\n" for r in rows))
    run_all = a.out / "run_all.sh"
    run_all.write_text("#!/bin/bash\nset -u\n" + "".join(
        f'echo "=== {s.parent.relative_to(a.out)} ==="; bash "{s.resolve()}" || echo "FAILED: {s}"\n'
        for s in run_scripts
    ))
    run_all.chmod(0o755)

    print(f"staged {len(rows)} replicate(s) under {a.out}")
    print(f"  project script : {paths['single']}")
    print(f"  env preamble   : {paths['env']}")
    print(f"  run all        : {run_all}")

    if a.run:
        print(f"\nexecuting {len(run_scripts)} replicate(s) (jobs={a.jobs}) ...")
        _execute(run_scripts, a.jobs, a.out)
    return 0


def _execute(scripts: list[Path], jobs: int, out: Path) -> None:
    from concurrent.futures import ThreadPoolExecutor

    def one(s: Path) -> tuple[Path, int]:
        log = s.parent / "run.log"
        with log.open("w") as fh:
            rc = subprocess.run(["bash", str(s)], stdout=fh, stderr=subprocess.STDOUT).returncode
        return s, rc

    results = []
    with ThreadPoolExecutor(max_workers=max(1, jobs)) as ex:
        for s, rc in ex.map(one, scripts):
            tag = "ok" if rc == 0 else f"FAIL(rc={rc})"
            print(f"  {tag:12} {s.parent.relative_to(out)}")
            results.append((str(s), rc))
    (out / "run_results.json").write_text(json.dumps(results, indent=2))
    n_fail = sum(1 for _, rc in results if rc)
    print(f"done: {len(results) - n_fail} ok, {n_fail} failed")


if __name__ == "__main__":
    raise SystemExit(main())
