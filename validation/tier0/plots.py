"""Tier 0 gate figures.

    python -m validation.tier0.plots --run validation/runs/tier0_gate2 --out <dir>

Reads a staged run tree (or the compact gate_result/ tree) via
``harness.tier0_adapter`` and writes:

    separation.png        detected-CAAS count per gene, planted (power) vs every
                          gene (null) — the load-bearing gate figure
    precision_at_k.png    per power replicate, top-n_planted genes that are planted
    site_recovery.png     planted-site recall by mechanism x CAAP scheme
    contrast_recovery.png operative-foreground Jaccard per replicate
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

_REPO = Path(__file__).resolve().parents[2]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from validation.harness import tier0_adapter as t0  # noqa: E402

_C = {"planted": "#2f6f4f", "null": "#b04a4a", "binary": "#3a6ea5", "rate": "#c9832b"}


def _archetypes(reps):
    return sorted({r.archetype for r in reps if r.archetype})


def fig_separation(reps, out: Path) -> None:
    arches = _archetypes(reps)
    fig, axes = plt.subplots(1, len(arches), figsize=(5 * len(arches), 4), squeeze=False)
    for ax, a in zip(axes[0], arches):
        planted, null_all = [], []
        for r in reps:
            if r.archetype != a or r.results_dir is None:
                continue
            npos = r.gene_n_positions()
            if r.is_null:
                null_all += list(npos.values())
            else:
                planted += [npos[g] for g in r.planted_genes if g in npos]
        mx = int(max(planted + null_all + [1]))
        bins = range(0, mx + 2)
        ax.hist(null_all, bins=bins, color=_C["null"], alpha=0.75,
                label=f"null-replicate genes (n={len(null_all)})")
        ax.hist(planted, bins=bins, color=_C["planted"], alpha=0.75,
                label=f"planted genes (n={len(planted)})")
        ax.set_title(f"{a}: detected CAAS per gene")
        ax.set_xlabel("n_positions (detected CAAS)")
        ax.set_ylabel("genes")
        ax.set_yscale("log")
        ax.legend(fontsize=8, loc="upper center")
    fig.tight_layout()
    fig.savefig(out / "separation.png", dpi=130)
    plt.close(fig)


def fig_precision_at_k(scored, out: Path) -> None:
    arches = sorted({s.archetype for s in scored})
    fig, ax = plt.subplots(figsize=(7, 4))
    for i, a in enumerate(arches):
        vals = [s.gene_precision_at_k_npos for s in scored if s.archetype == a]
        vals2 = [s.gene_precision_at_k for s in scored if s.archetype == a]
        x = [i - 0.15] * len(vals)
        ax.scatter(x, vals, color=_C.get(a, "#555"), s=40, label=f"{a} (by CAAS count)")
        ax.scatter([i + 0.15] * len(vals2), vals2, color=_C.get(a, "#555"),
                   s=40, marker="x", label=f"{a} (by gene_caas_score)")
    ax.set_xticks(range(len(arches)))
    ax.set_xticklabels(arches)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("precision@k  (top-n_planted genes that are planted)")
    ax.set_title("Gene recovery per power replicate")
    ax.axhline(0.8, ls="--", c="grey", lw=0.8)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "precision_at_k.png", dpi=130)
    plt.close(fig)


def fig_site_recovery(scored, out: Path) -> None:
    arches = sorted({s.archetype for s in scored})
    metrics = [
        ("site_recall", "any scheme"),
        ("identical_aa_us_recall", "identical_aa → US"),
        ("grouped_caap_gs_recall", "grouped_caap → GS"),
        ("grouped_caap_us_leakage", "grouped_caap also US"),
    ]
    fig, ax = plt.subplots(figsize=(8, 4))
    w = 0.8 / len(arches)
    for ai, a in enumerate(arches):
        group = [s for s in scored if s.archetype == a]
        means = [t0._mean([getattr(s, m) for s in group]) for m, _ in metrics]
        xs = [j + ai * w for j in range(len(metrics))]
        ax.bar(xs, means, width=w, color=_C.get(a, "#555"), label=a)
        for x, m in zip(xs, means):
            ax.text(x, m + 0.02, f"{m:.2f}", ha="center", fontsize=8)
    ax.set_xticks([j + w * (len(arches) - 1) / 2 for j in range(len(metrics))])
    ax.set_xticklabels([lbl for _, lbl in metrics], rotation=15, ha="right")
    ax.set_ylim(0, 1.15)
    ax.set_ylabel("mean over power replicates")
    ax.set_title("Planted-site recovery by mechanism × scheme")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "site_recovery.png", dpi=130)
    plt.close(fig)


def fig_contrast_recovery(scored, out: Path) -> None:
    arches = sorted({s.archetype for s in scored})
    fig, ax = plt.subplots(figsize=(7, 4))
    for i, a in enumerate(arches):
        jac = [s.contrast_jaccard for s in scored if s.archetype == a]
        pr = [s.contrast_pairs_recovered / s.contrast_pairs_total
              for s in scored if s.archetype == a and s.contrast_pairs_total]
        ax.scatter([i - 0.12] * len(jac), jac, color=_C.get(a, "#555"), s=40, label=f"{a} Jaccard")
        ax.scatter([i + 0.12] * len(pr), pr, color=_C.get(a, "#555"), s=40, marker="^",
                   label=f"{a} pairs recovered")
    ax.set_xticks(range(len(arches)))
    ax.set_xticklabels(arches)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("operative fg vs planted pairs")
    ax.set_title("Contrast-selection recovery per replicate")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "contrast_recovery.png", dpi=130)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="validation.tier0.plots")
    ap.add_argument("--run", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)
    out = a.out or (a.run / "figures")
    out.mkdir(parents=True, exist_ok=True)

    reps = t0.find_replicates(a.run)
    res = t0.score(a.run)
    scored = res.per_replicate

    fig_separation(reps, out)
    fig_precision_at_k(scored, out)
    fig_site_recovery(scored, out)
    fig_contrast_recovery(scored, out)
    print(f"wrote 4 figures to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
