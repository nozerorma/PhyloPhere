"""Emit the Tier 0 ``null`` and ``power`` replicate sets.

    python -m validation.tier0.make_replicates --config validation/tier0/grid.json \
        --out validation/runs/tier0 [--set null|power|all] [--limit N] [--dry-run]

Each replicate directory contains::

    aln.fasta        simulated amino-acid alignment
    phenotype.tsv    species <tab> 1|0
    truth.json       planted sites + mechanism + every parameter + seed
    traits/          the full PhyloPhere trait matrix for this phenotype
                     (harness.trait_matrix — categorical/direct is the one the
                     CAAS/FADE modules consume; the other cells exercise the
                     trait-mode plumbing)

and the set directory gets ``manifest.jsonl`` (one line per replicate) for the
calibration step to iterate over.

The grid config is JSON::

    {
      "trees": {
        "primate":   {"path": "validation/fixtures/tier0/trees/primates_233_subst.tree",
                      "k": 50, "scale": 1.0},
        "primate_x5":{"path": "...primates_233_subst.tree", "k": 50, "scale": 5.0},
        "mammal":    {"path": ".../speciesTree_speciesname_pruned.nh", "k": 50, "scale": 1.0},
        "star":      {"builtin": "star", "k": 12, "branch_length": 0.3},
        "ladder":    {"builtin": "ladder", "k": 16, "branch_length": 0.15}
      },
      "base_sim": {"n_sites": 400, "concentration": 4.0, "gamma_alpha": 0.7, "matrix": "wag"},
      "null":  {"n_replicates": 1000, "transitions": [2, 3, 5, 7]},
      "power": {"n_replicates": 1000, "transitions": [2, 3, 5, 7],
                "n_planted_profile_shift": 20, "n_planted_identical_aa": 20,
                "sweep": {"concentration": [2.0, 4.0, 8.0],
                          "shift_concentration": [1.0, 2.0]}}
    }
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[2]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from validation.harness.phenotype import DatasetSpec, PhenotypeSpec  # noqa: E402
from validation.harness.trait_matrix import emit_trait_matrix  # noqa: E402
from validation.tier0 import trees as T  # noqa: E402
from validation.tier0.simulate import SimConfig, simulate  # noqa: E402


def _load_tree(spec: dict, rng: np.random.Generator) -> T.PhyloTree:
    if spec.get("builtin") == "star":
        return T.star_tree(spec["k"], spec.get("branch_length", 0.3))
    if spec.get("builtin") == "ladder":
        return T.ladder_tree(spec["k"], spec.get("branch_length", 0.15))
    path = _ROOT / spec["path"] if not Path(spec["path"]).is_absolute() else Path(spec["path"])
    if not path.exists():
        raise FileNotFoundError(
            f"{path} — fetch the Tier 0 trees first (see validation/fixtures/tier0/README.md)"
        )
    tree = T.prune_depth_preserving(path, spec["k"], rng)
    if spec.get("scale", 1.0) != 1.0:
        tree = T.scale_branches(tree, spec["scale"])
    return tree


def _phenotype_spec(res, name: str) -> PhenotypeSpec:
    return PhenotypeSpec(
        name=name, kind="categorical",
        labels={t: ("fg" if s else "bg") for t, s in res.phenotype.items()},
        # pairs left empty -> harness order-pairs and records it in the manifest
    )


def _emit_replicate(tree, fg, cfg: SimConfig, seed: int, outdir: Path,
                    tree_name: str, emit_traits: bool) -> dict:
    res = simulate(tree, fg, cfg, seed=seed)
    res.write(outdir)
    # newick for this (possibly pruned/scaled) tree
    _write_newick(tree, outdir / "tree.nwk")
    row = {
        "dir": str(outdir.relative_to(outdir.parents[1])),
        "tree": tree_name, "seed": seed,
        "n_tips": res.truth.n_tips, "n_transitions": res.truth.n_transitions,
        "tree_total_length": round(res.truth.tree_total_length, 4),
        "n_planted": len(res.truth.planted_sites),
        "planted_by_mechanism": _count_mech(res.truth.planted_sites),
        "sim": asdict(cfg),
    }
    if emit_traits:
        ds = DatasetSpec(name=f"tier0_{tree_name}_{seed}", tree=outdir / "tree.nwk",
                         traits=[_phenotype_spec(res, "sim_fg")], provenance="tier0 simulation")
        m = emit_trait_matrix(ds, outdir / "traits")
        row["trait_warnings"] = m["warnings"]
    return row


def _count_mech(planted: dict) -> dict:
    out: dict[str, int] = {}
    for v in planted.values():
        out[v] = out.get(v, 0) + 1
    return out


def _write_newick(tree: T.PhyloTree, path: Path) -> None:
    children: dict[int, list[int]] = {}
    for p, c, bl in tree.edges:
        children.setdefault(p, []).append(c)
    bl_of = {c: bl for _, c, bl in tree.edges}

    def render(node: int) -> str:
        kids = children.get(node, [])
        if not kids:
            return f"{tree.labels[node]}:{bl_of.get(node, 0.0):.8f}"
        inner = ",".join(render(k) for k in kids)
        if node == tree.root:
            return f"({inner});"
        return f"({inner}){tree.labels[node]}:{bl_of.get(node, 0.0):.8f}"

    path.write_text(render(tree.root) + "\n")


def _iter_null(grid: dict, out: Path, limit: int | None):
    conf = grid["null"]
    base = grid["base_sim"]
    n = conf["n_replicates"]
    trans = conf["transitions"]
    tree_names = list(grid["trees"])
    count = 0
    for i in range(n):
        if limit and count >= limit:
            return
        seed = 10_000 + i
        rng = np.random.default_rng(seed)
        tree_name = tree_names[i % len(tree_names)]
        ntr = trans[i % len(trans)]
        tree = _load_tree(grid["trees"][tree_name], rng)
        try:
            fg = T.sample_foreground(tree, ntr, rng, max_clade=3)
        except ValueError:
            fg = T.sample_foreground(tree, max(2, ntr - 1), rng, max_clade=3)
        cfg = SimConfig(**base)  # no planted signal
        yield ("null", tree, fg, cfg, seed, out / "null" / f"rep{i:04d}_{tree_name}", tree_name)
        count += 1


def _iter_power(grid: dict, out: Path, limit: int | None):
    conf = grid["power"]
    base = dict(grid["base_sim"])
    n = conf["n_replicates"]
    trans = conf["transitions"]
    tree_names = list(grid["trees"])
    sweep = conf.get("sweep", {})
    sweep_keys = list(sweep)
    sweep_combos = list(itertools.product(*sweep.values())) or [()]
    count = 0
    for i in range(n):
        if limit and count >= limit:
            return
        seed = 20_000 + i
        rng = np.random.default_rng(seed)
        tree_name = tree_names[i % len(tree_names)]
        ntr = trans[i % len(trans)]
        combo = sweep_combos[i % len(sweep_combos)]
        overrides = dict(zip(sweep_keys, combo))
        tree = _load_tree(grid["trees"][tree_name], rng)
        try:
            fg = T.sample_foreground(tree, ntr, rng, max_clade=3)
        except ValueError:
            fg = T.sample_foreground(tree, max(2, ntr - 1), rng, max_clade=3)
        params = dict(base)
        params.update(overrides)  # sweep overrides base_sim (e.g. concentration)
        params["n_planted_profile_shift"] = conf.get("n_planted_profile_shift", 20)
        params["n_planted_identical_aa"] = conf.get("n_planted_identical_aa", 20)
        cfg = SimConfig(**params)
        tag = "_".join(f"{k}{v}" for k, v in overrides.items())
        yield ("power", tree, fg, cfg, seed,
               out / "power" / f"rep{i:04d}_{tree_name}_{tag}".rstrip("_"), tree_name)
        count += 1


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="validation.tier0.make_replicates")
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--set", choices=["null", "power", "all"], default="all")
    ap.add_argument("--limit", type=int, default=None, help="cap replicates per set (smoke test)")
    ap.add_argument("--no-traits", action="store_true", help="skip trait-matrix emission")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)

    grid = json.loads(args.config.read_text())
    args.out.mkdir(parents=True, exist_ok=True)

    iters = []
    if args.set in ("null", "all"):
        iters.append(("null", _iter_null(grid, args.out, args.limit)))
    if args.set in ("power", "all"):
        iters.append(("power", _iter_power(grid, args.out, args.limit)))

    for set_name, it in iters:
        rows = []
        for _, tree, fg, cfg, seed, repdir, tree_name in it:
            if args.dry_run:
                rows.append({"dir": str(repdir), "tree": tree_name, "seed": seed,
                             "n_tips": len(tree.tips), "n_transitions": fg.n_transitions})
                continue
            repdir.mkdir(parents=True, exist_ok=True)
            rows.append(_emit_replicate(tree, fg, cfg, seed, repdir, tree_name,
                                        emit_traits=not args.no_traits))
        man = args.out / set_name / "manifest.jsonl"
        man.parent.mkdir(parents=True, exist_ok=True)
        man.write_text("".join(json.dumps(r) + "\n" for r in rows))
        print(f"{set_name}: {len(rows)} replicate(s) -> {man}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
