"""Adapter: a Tier 0 ``runs/`` tree -> harness metric inputs.

A staged Tier 0 run (``validation.tier0.run_replicates``) lays out

    <run_root>/<subset>/rep###/
        truth.json                     planted genes / sites / foreground tips
        out*/scoring/gene_scores.tsv
        out*/scoring/position_scores.tsv
        out*/caas_permulation/perm_pos_pval.tsv
        out*/signification/meta_caas/US_meta_caas.tsv
        out*/data_exploration/2.CT/1.Traitfiles/traitfile.tab

This module turns that into:

* ``calibrate`` — pools a chosen p-value column across every **null** context
  (whole null-set replicates + the unplanted genes of power replicates) and
  runs ``metrics.null_calibration`` (the KS uniformity gate).
* ``score`` — per power replicate: planted-gene recovery
  (``metrics.score_metrics`` on ``gene_caas_score``), planted-site
  detection precision/recall (SCORING only emits *detected* positions, so a
  full ROC is not defined), and contrast-selection recovery
  (operative foreground from ``traitfile.tab`` vs ``truth.foreground_tips``).

Indexing note: ``truth["genes"][g]["planted_sites"]`` keys are 0-based
alignment columns and the pipeline's ``Position`` column matches them
directly (verified on the smoke run — 12/12 exact, 0/12 at an offset of 1).
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass, field
from pathlib import Path

from .metrics import NullReport, ScoreReport, null_calibration, score_metrics

# name -> (relative path under the results dir, column)
PCOLS: dict[str, tuple[str, str]] = {
    "meta_caas_boot": ("signification/meta_caas/US_meta_caas.tsv", "pvalue_boot"),
    "meta_caas_hyp": ("signification/meta_caas/US_meta_caas.tsv", "pvalue"),
    "perm_pos_boot": ("caas_permulation/perm_pos_pval.tsv", "null_pvalue_boot"),
    "position_boot": ("scoring/position_scores.tsv", "pvalue_boot"),
    "position_hyp": ("scoring/position_scores.tsv", "pvalue"),
}


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _num(x: str) -> float:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return math.nan
    return v


# --------------------------------------------------------------------------- #
# replicate view
# --------------------------------------------------------------------------- #
@dataclass
class Replicate:
    repdir: Path
    truth: dict
    subset: str

    @property
    def is_null(self) -> bool:
        return not self.truth.get("planted", False)

    @property
    def archetype(self) -> str:
        return self.truth.get("archetype", "")

    @property
    def planted_genes(self) -> set[str]:
        return {g for g, v in self.truth.get("genes", {}).items() if v.get("planted")}

    def planted_sites(self, gene: str) -> set[int]:
        v = self.truth.get("genes", {}).get(gene, {})
        return {int(k) for k in v.get("planted_sites", {})}

    def site_mechanism(self, gene: str) -> dict[int, str]:
        v = self.truth.get("genes", {}).get(gene, {})
        return {int(k): m for k, m in v.get("planted_sites", {}).items()}

    @property
    def foreground_tips(self) -> set[str]:
        return set(self.truth.get("phenotype", {}).get("foreground_tips", []))

    # -- pipeline outputs --------------------------------------------------- #
    @property
    def results_dir(self) -> Path | None:
        hits = sorted(self.repdir.glob("**/scoring/gene_scores.tsv"))
        return hits[0].parent.parent if hits else None

    def _out(self, rel: str) -> Path | None:
        rd = self.results_dir
        if rd is None:
            return None
        p = rd / rel
        return p if p.exists() else None

    def gene_scores(self) -> dict[str, float]:
        p = self._out("scoring/gene_scores.tsv")
        if p is None:
            return {}
        out = {}
        for r in _read_tsv(p):
            out[r["Gene"]] = _num(r.get("gene_caas_score", "nan"))
        return {g: v for g, v in out.items() if not math.isnan(v)}

    def position_rows(self) -> list[dict[str, str]]:
        p = self._out("scoring/position_scores.tsv")
        return _read_tsv(p) if p else []

    def traitfile(self) -> list[tuple[str, str, str]]:
        p = self._out("data_exploration/2.CT/1.Traitfiles/traitfile.tab")
        if p is None:
            return []
        rows = []
        for ln in p.read_text().splitlines():
            f = ln.split("\t")
            if len(f) >= 3:
                rows.append((f[0], f[1], f[2]))
        return rows

    def pvalues(self, pcol: str, *, unplanted_only: bool = False) -> list[float]:
        rel, col = PCOLS[pcol]
        p = self._out(rel)
        if p is None:
            return []
        planted = self.planted_genes
        out = []
        for r in _read_tsv(p):
            if unplanted_only and r.get("Gene") in planted:
                continue
            v = _num(r.get(col, "nan"))
            if not math.isnan(v):
                out.append(v)
        return out


def find_replicates(run_root: Path) -> list[Replicate]:
    reps = []
    for tj in sorted(run_root.glob("**/truth.json")):
        try:
            truth = json.loads(tj.read_text())
        except json.JSONDecodeError:
            continue
        if "genes" not in truth:
            continue
        repdir = tj.parent
        subset = repdir.parent.name
        reps.append(Replicate(repdir=repdir, truth=truth, subset=subset))
    return reps


# --------------------------------------------------------------------------- #
# calibrate
# --------------------------------------------------------------------------- #
@dataclass
class CalibrateResult:
    pcol: str
    n_replicates: int
    n_pvalues: int
    report: NullReport
    verdict: str

    def to_dict(self) -> dict:
        r = self.report
        return {
            "pcol": self.pcol,
            "n_replicates": self.n_replicates,
            "n_pvalues": self.n_pvalues,
            "ks_stat": r.ks_stat,
            "ks_pvalue": r.ks_pvalue,
            "mean": r.mean,
            "frac_exact_zero": r.frac_exact_zero,
            "frac_exact_one": r.frac_exact_one,
            "type1_error": r.type1_error,
            "verdict": self.verdict,
        }


def calibrate(run_root: Path, pcol: str, *, alpha: float = 0.05) -> CalibrateResult:
    if pcol not in PCOLS:
        raise SystemExit(f"--pcol must be one of {sorted(PCOLS)}")
    reps = find_replicates(run_root)
    pooled: list[float] = []
    used = 0
    for rp in reps:
        if rp.is_null:
            vals = rp.pvalues(pcol)
        else:
            vals = rp.pvalues(pcol, unplanted_only=True)
        if vals:
            used += 1
            pooled.extend(vals)
    if not pooled:
        raise SystemExit(
            f"no p-values for --pcol {pcol} under {run_root} "
            f"({len(reps)} replicate(s) found; is the run finished?)"
        )
    rep = null_calibration(pooled)
    return CalibrateResult(pcol, used, len(pooled), rep, rep.verdict(alpha))


# --------------------------------------------------------------------------- #
# score
# --------------------------------------------------------------------------- #
@dataclass
class ReplicateScore:
    subset: str
    rep: str
    archetype: str
    gene: ScoreReport | None
    gene_precision_at_k: float
    gene_planted_ranks: dict[str, int]
    n_genes_scored: int
    site_recall: float
    site_precision: float
    site_recall_by_mechanism: dict[str, float]
    n_sites_planted: int
    n_sites_detected: int
    contrast_jaccard: float
    contrast_fg_precision: float
    contrast_pairs_recovered: int
    contrast_pairs_total: int

    def to_dict(self) -> dict:
        d = dict(self.__dict__)
        d["gene"] = self.gene.as_dict() if self.gene else None
        return d


def _score_one(rp: Replicate) -> ReplicateScore | None:
    planted = rp.planted_genes
    if not planted:
        return None
    gs = rp.gene_scores()
    gene_report = None
    prec_k = math.nan
    ranks: dict[str, int] = {}
    if gs:
        thr = min(gs.values())
        gene_report = score_metrics(gs, planted, thr, higher_is_hit=True)
        ranks = gene_report.positive_ranks
        k = len(planted)
        topk = {g for g, _ in sorted(gs.items(), key=lambda kv: -kv[1])[:k]}
        prec_k = len(topk & planted) / k

    # -- sites: SCORING only lists detected positions --------------------- #
    detected: dict[str, set[int]] = {}
    for r in rp.position_rows():
        detected.setdefault(r["Gene"], set()).add(int(r["Position"]))
    all_planted_sites = 0
    hit = 0
    total_detected = 0
    mech_tot: dict[str, int] = {}
    mech_hit: dict[str, int] = {}
    for g in planted:
        ps = rp.planted_sites(g)
        mech = rp.site_mechanism(g)
        det = detected.get(g, set())
        all_planted_sites += len(ps)
        total_detected += len(det)
        hit += len(ps & det)
        for s, m in mech.items():
            mech_tot[m] = mech_tot.get(m, 0) + 1
            if s in det:
                mech_hit[m] = mech_hit.get(m, 0) + 1
    site_recall = hit / all_planted_sites if all_planted_sites else math.nan
    site_prec = hit / total_detected if total_detected else math.nan
    mech_recall = {m: mech_hit.get(m, 0) / n for m, n in mech_tot.items()}

    # -- contrast selection recovery (D-T0-A) ---------------------------- #
    tf = rp.traitfile()
    fg_truth = rp.foreground_tips
    op_fg = {sp for sp, lab, _ in tf if lab == "1"}
    jac = (
        len(op_fg & fg_truth) / len(op_fg | fg_truth)
        if (op_fg or fg_truth)
        else math.nan
    )
    fg_prec = len(op_fg & fg_truth) / len(op_fg) if op_fg else math.nan
    pairs: dict[str, list[tuple[str, str]]] = {}
    for sp, lab, pid in tf:
        pairs.setdefault(pid, []).append((sp, lab))
    recovered = sum(
        1
        for members in pairs.values()
        if any(lab == "1" and sp in fg_truth for sp, lab in members)
    )

    return ReplicateScore(
        subset=rp.subset,
        rep=rp.repdir.name,
        archetype=rp.archetype,
        gene=gene_report,
        gene_precision_at_k=prec_k,
        gene_planted_ranks=ranks,
        n_genes_scored=len(gs),
        site_recall=site_recall,
        site_precision=site_prec,
        site_recall_by_mechanism=mech_recall,
        n_sites_planted=all_planted_sites,
        n_sites_detected=total_detected,
        contrast_jaccard=jac,
        contrast_fg_precision=fg_prec,
        contrast_pairs_recovered=recovered,
        contrast_pairs_total=len(pairs),
    )


def _mean(xs: list[float]) -> float:
    xs = [x for x in xs if x is not None and not math.isnan(x)]
    return sum(xs) / len(xs) if xs else math.nan


@dataclass
class ScoreResult:
    n_power_replicates: int
    per_replicate: list[ReplicateScore] = field(default_factory=list)
    summary: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "n_power_replicates": self.n_power_replicates,
            "summary": self.summary,
            "per_replicate": [r.to_dict() for r in self.per_replicate],
        }


def score(run_root: Path) -> ScoreResult:
    reps = find_replicates(run_root)
    scored = [s for s in (_score_one(rp) for rp in reps) if s is not None]
    if not scored:
        raise SystemExit(
            f"no power replicates with planted genes under {run_root} "
            f"({len(reps)} replicate(s) found)"
        )
    summary = {
        "gene_roc_auc": _mean([s.gene.roc_auc for s in scored if s.gene]),
        "gene_pr_auc": _mean([s.gene.pr_auc for s in scored if s.gene]),
        "gene_precision_at_k": _mean([s.gene_precision_at_k for s in scored]),
        "site_recall": _mean([s.site_recall for s in scored]),
        "site_precision": _mean([s.site_precision for s in scored]),
        "contrast_jaccard": _mean([s.contrast_jaccard for s in scored]),
        "contrast_fg_precision": _mean([s.contrast_fg_precision for s in scored]),
        "contrast_pairs_recovered_frac": _mean(
            [
                s.contrast_pairs_recovered / s.contrast_pairs_total
                for s in scored
                if s.contrast_pairs_total
            ]
        ),
    }
    mechs: dict[str, list[float]] = {}
    for s in scored:
        for m, v in s.site_recall_by_mechanism.items():
            mechs.setdefault(m, []).append(v)
    summary["site_recall_by_mechanism"] = {m: _mean(v) for m, v in mechs.items()}
    return ScoreResult(n_power_replicates=len(scored), per_replicate=scored, summary=summary)
