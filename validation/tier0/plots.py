"""Tier 0 gate figures.

    python -m validation.tier0.plots --run validation/runs/tier0_gate --out <dir>

Writes:
    separation.png          detected-CAAS count per gene, planted vs null, per (archetype, lambda)
    recovery_vs_lambda.png  contrast / site / gene recovery as lambda goes 0 -> 1
    score_separation.png    null-vs-power AUC by quantity, as lambda goes 0 -> 1
    site_recovery.png       planted-site recall by mechanism x CAAP scheme
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

_AC = {"binary": "#3a6ea5", "rate": "#c9832b", "continuous": "#4c9a5a"}


def _lams(reps):
    ls = sorted({round(r.lam, 3) for r in reps if r.results_dir is not None})
    return ls or [-1.0]


def _lfmt(lm) -> str:
    return "n/a" if lm < 0 else f"{lm:g}"


def _arches(reps):
    return sorted({r.archetype for r in reps if r.archetype})


def fig_separation(reps, out: Path) -> None:
    arches, lams = _arches(reps), _lams(reps)
    fig, axes = plt.subplots(len(arches), len(lams),
                             figsize=(4 * len(lams), 3.2 * len(arches)), squeeze=False)
    for ai, a in enumerate(arches):
        for li, lm in enumerate(lams):
            ax = axes[ai][li]
            planted, null_all = [], []
            for r in reps:
                if r.results_dir is None or r.archetype != a:
                    continue
                if lm >= 0 and round(r.lam, 3) != lm:
                    continue
                npos = r.gene_n_positions()
                if r.is_null:
                    null_all += list(npos.values())
                else:
                    planted += [npos[g] for g in r.planted_genes if g in npos]
            mx = int(max(planted + null_all + [1]))
            bins = range(0, mx + 2)
            ax.hist(null_all, bins=bins, color="#b04a4a", alpha=0.75, label="null genes")
            ax.hist(planted, bins=bins, color="#2f6f4f", alpha=0.75, label="planted genes")
            ax.set_yscale("log")
            ax.set_title(f"{a}  lambda={_lfmt(lm)}", fontsize=10)
            if ai == len(arches) - 1:
                ax.set_xlabel("detected CAAS / gene")
            if li == 0:
                ax.set_ylabel("genes")
            if ai == 0 and li == 0:
                ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "separation.png", dpi=130)
    plt.close(fig)


def fig_recovery_vs_lambda(scored, out: Path) -> None:
    arches = sorted({s.archetype for s in scored})
    lams = sorted({round(s.lam, 3) for s in scored})
    series = [
        ("gene_precision_at_k_npos", "gene precision@k"),
        ("site_recall", "site recall"),
        ("site_directional_recall", "site recall (correct side)"),
        ("contrast_anchor_recall", "contrast: true anchors found"),
        ("contrast_pair_exact", "contrast: pairs exact"),
        ("contrast_jaccard", "contrast: lineage Jaccard"),
    ]
    fig, axes = plt.subplots(1, len(arches), figsize=(5.5 * len(arches), 4), squeeze=False)
    for ax, a in zip(axes[0], arches):
        for attr, lbl in series:
            ys = [t0._mean([getattr(s, attr) for s in scored
                            if s.archetype == a and round(s.lam, 3) == lm]) for lm in lams]
            ax.plot(lams, ys, marker="o", label=lbl)
        ax.set_title(a)
        ax.set_xlabel("Pagel's lambda")
        ax.set_ylabel("mean over power replicates")
        ax.set_ylim(0, 1.05)
        ax.set_xticks(lams)
        ax.legend(fontsize=7, loc="lower left")
        # secondary axis: how many independent pairs the latent could support
        ax2 = ax.twinx()
        npos = [t0._mean([s.n_possible_pairs for s in scored
                          if s.archetype == a and round(s.lam, 3) == lm]) for lm in lams]
        ax2.plot(lams, npos, marker="s", ls="--", color="grey", label="n possible pairs")
        ax2.set_ylabel("n independent pairs the trait supports", color="grey")
        ax2.tick_params(axis="y", colors="grey")
    fig.tight_layout()
    fig.savefig(out / "recovery_vs_lambda.png", dpi=130)
    plt.close(fig)


def fig_score_separation(seps, out: Path) -> None:
    arches = sorted({s.archetype for s in seps})
    lams = sorted({round(s.lam, 3) for s in seps})
    quantities = [
        ("auc_npos", "detected-CAAS count"),
        ("auc_pos_neglog_pvalue", "position -log10 p"),
        ("auc_pos_caas_score", "position CAAS_score"),
        ("auc_gene_caas_score", "gene_caas_score"),
    ]
    fig, axes = plt.subplots(1, len(arches), figsize=(5.5 * len(arches), 4), squeeze=False)
    for ax, a in zip(axes[0], arches):
        for attr, lbl in quantities:
            ys = [next((getattr(s, attr) for s in seps
                        if s.archetype == a and round(s.lam, 3) == lm), float("nan"))
                  for lm in lams]
            ax.plot(lams, ys, marker="o", label=lbl)
        ax.axhline(0.5, ls=":", c="grey", lw=0.8)
        ax.set_title(f"{a}: null vs power AUC")
        ax.set_xlabel("Pagel's lambda")
        ax.set_ylabel("AUC(planted vs null)")
        ax.set_ylim(0, 1.05)
        ax.set_xticks(lams)
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "score_separation.png", dpi=130)
    plt.close(fig)


def fig_site_recovery(scored, out: Path) -> None:
    arches = sorted({s.archetype for s in scored})
    metrics = [
        ("site_recall", "any scheme"),
        ("identical_aa_us_recall", "identical_aa -> US"),
        ("grouped_caap_gs_recall", "grouped_caap -> GS"),
        ("grouped_caap_us_leakage", "grouped_caap also US"),
    ]
    fig, ax = plt.subplots(figsize=(8, 4))
    w = 0.8 / len(arches)
    for ai, a in enumerate(arches):
        g = [s for s in scored if s.archetype == a]
        means = [t0._mean([getattr(s, m) for s in g]) for m, _ in metrics]
        xs = [j + ai * w for j in range(len(metrics))]
        ax.bar(xs, means, width=w, color=_AC.get(a, "#555"), label=a)
        for x, m in zip(xs, means):
            if m == m:
                ax.text(x, m + 0.02, f"{m:.2f}", ha="center", fontsize=8)
    ax.set_xticks([j + w * (len(arches) - 1) / 2 for j in range(len(metrics))])
    ax.set_xticklabels([lbl for _, lbl in metrics], rotation=15, ha="right")
    ax.set_ylim(0, 1.15)
    ax.set_ylabel("mean over all power replicates")
    ax.set_title("Planted-site recovery by mechanism x scheme")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "site_recovery.png", dpi=130)
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

    fig_separation(reps, out)
    fig_recovery_vs_lambda(res.per_replicate, out)
    fig_score_separation(res.separation, out)
    fig_site_recovery(res.per_replicate, out)
    print(f"wrote 4 figures to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
