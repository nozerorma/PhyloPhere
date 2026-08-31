"""Adapter: a Tier 0 ``runs/`` tree -> gate metrics.

A staged Tier 0 run (``validation.tier0.run_replicates``) lays out

    <run_root>/<archetype>_<set>_<tree>_lam<tag>/rep###/
        truth.json
        out/<trait>_complete/scoring/{gene_scores,position_scores}.tsv
        out/<trait>_complete/scoring/gene_lists/slice_*.tsv
        out/<trait>_complete/data_exploration/2.CT/1.Traitfiles/traitfile.tab

Metrics are reported **per (archetype, lambda) cell** (D-T0-Q). The pipeline is a
prioritisation engine (``project_tier0_scoping_reframe``); its permulation
p-values are foreground-specificity scores, never Uniform(0,1) — no KS gate.

  1. null vs power (occurrence) — AUC of the detected-CAAS count of planted genes
     vs every gene in the matched null replicates. Score-level AUCs (position
     CAAS_score, position -log10 p, gene_caas_score) reported, NOT gated.
  2. prioritisation — precision@k by CAAS count and by gene_caas_score, planted
     rank percentile, slice_global25.
  3. site recovery — detection recall, directional recall (right top/bottom
     side), split by mechanism x CAAP scheme.
  4. contrast recovery — operative pairs (traitfile.tab) vs the TRUE lambda-drawn
     pairs (anchor / partner recall, exact-pair, lineage Jaccard). n_pairs is
     FIXED; n_possible_pairs (how many the latent supported) is reported.
     REPORTED as a lambda-curve, NOT gated.

VERDICT: occurrence separation at every cell + precision@k/site_recall at the
lowest lambda.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass, field
from pathlib import Path

from .metrics import ScoreReport, score_metrics


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _num(x: str) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def _mean(xs) -> float:
    xs = [x for x in xs if x is not None and not (isinstance(x, float) and math.isnan(x))]
    return sum(xs) / len(xs) if xs else math.nan


def _quantile(xs: list[float], q: float) -> float:
    xs = sorted(x for x in xs if not math.isnan(x))
    if not xs:
        return math.nan
    i = q * (len(xs) - 1)
    lo, hi = int(math.floor(i)), int(math.ceil(i))
    return xs[lo] if lo == hi else xs[lo] + (xs[hi] - xs[lo]) * (i - lo)


def _auc(pos: list[float], neg: list[float]) -> float:
    """Mann-Whitney AUC: P(a random positive scores above a random negative)."""
    pos = [x for x in pos if not math.isnan(x)]
    neg = [x for x in neg if not math.isnan(x)]
    if not pos or not neg:
        return math.nan
    wins = 0.0
    for p in pos:
        wins += sum(1.0 for n in neg if n < p) + 0.5 * sum(1.0 for n in neg if n == p)
    return wins / (len(pos) * len(neg))


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
    def lam(self) -> float:
        return float(self.truth.get("phenotype", {}).get("lambda", -1.0))

    @property
    def n_possible_pairs(self) -> float:
        return float(self.truth.get("phenotype", {}).get("notes", {}).get("n_possible_pairs", math.nan))

    @property
    def planted_genes(self) -> set[str]:
        return {g for g, v in self.truth.get("genes", {}).items() if v.get("planted")}

    def gene_direction(self, gene: str) -> str | None:
        return self.truth.get("genes", {}).get(gene, {}).get("direction")

    def planted_sites(self, gene: str) -> set[int]:
        v = self.truth.get("genes", {}).get(gene, {})
        return {int(k) for k in v.get("planted_sites", {})}

    def site_mechanism(self, gene: str) -> dict[int, str]:
        v = self.truth.get("genes", {}).get(gene, {})
        return {int(k): m for k, m in v.get("planted_sites", {}).items()}

    @property
    def foreground_tips(self) -> set[str]:
        return set(self.truth.get("phenotype", {}).get("foreground_tips", []))

    @property
    def partner_tips(self) -> set[str]:
        return set(self.truth.get("phenotype", {}).get("partner_tips", []))

    @property
    def truth_pairs(self) -> list[tuple[str, str]]:
        return [tuple(p) for p in self.truth.get("phenotype", {}).get("pairs", [])]

    # -- pipeline outputs ------------------------------------------------------ #
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

    def gene_scores(self, col: str = "gene_caas_score") -> dict[str, float]:
        p = self._out("scoring/gene_scores.tsv")
        if p is None:
            return {}
        out = {r["Gene"]: _num(r.get(col, "nan")) for r in _read_tsv(p)}
        return {g: v for g, v in out.items() if not math.isnan(v)}

    def gene_n_positions(self) -> dict[str, float]:
        """Detected CAAS positions per gene — an absolute, replicate-comparable
        quantity (unlike gene_caas_score, which is a within-run percentile rank)."""
        p = self._out("scoring/gene_scores.tsv")
        if p is None:
            return {}
        return {r["Gene"]: _num(r.get("n_positions", "nan")) for r in _read_tsv(p)}

    def gene_slice(self, name: str) -> set[str]:
        """slice_global1 / slice_global5 / slice_top5 ... -> the gene set."""
        p = self._out(f"scoring/gene_lists/{name}.tsv")
        if p is None:
            return set()
        return {r["Gene"] for r in _read_tsv(p) if r.get("Gene")}

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


def find_replicates(run_root: Path) -> list[Replicate]:
    reps = []
    for tj in sorted(run_root.glob("**/truth.json")):
        try:
            truth = json.loads(tj.read_text())
        except json.JSONDecodeError:
            continue
        if "genes" not in truth:
            continue
        reps.append(Replicate(repdir=tj.parent, truth=truth, subset=tj.parent.parent.name))
    return reps


# --------------------------------------------------------------------------- #
# per-replicate power score
# --------------------------------------------------------------------------- #
@dataclass
class ReplicateScore:
    subset: str
    rep: str
    archetype: str
    n_genes_scored: int
    gene: ScoreReport | None
    gene_precision_at_k: float          # top-n_planted by gene_caas_score, fraction planted
    gene_precision_at_k_npos: float     # top-n_planted by detected-CAAS count, fraction planted
    gene_planted_rank_pctile_p50: float # median planted-gene rank as a top-fraction (1.0 = rank 1)
    gene_planted_ranks: dict[str, int]
    planted_in_slice_global5: float     # fraction of planted genes in the top-5% slice
    planted_in_slice_global25: float
    lam: float
    n_possible_pairs: float
    site_recall: float
    site_precision: float
    site_recall_by_mechanism: dict[str, float]
    site_directional_recall: float      # planted site detected AND on the right (top/bottom) side
    n_sites_planted: int
    n_sites_detected: int
    identical_aa_us_recall: float
    grouped_caap_gs_recall: float
    grouped_caap_us_leakage: float
    contrast_anchor_recall: float       # true top anchors that appear as fg in traitfile.tab
    contrast_partner_recall: float      # true bottom partners that appear as bg
    contrast_pair_exact: float          # true (a,b) pairs reproduced exactly
    contrast_jaccard: float             # {operative fg+bg} vs {true anchors+partners}
    contrast_n_operative_pairs: int

    def to_dict(self) -> dict:
        d = dict(self.__dict__)
        d["gene"] = self.gene.as_dict() if self.gene else None
        return d


def _score_one(rp: Replicate) -> ReplicateScore | None:
    planted = rp.planted_genes
    if not planted or rp.results_dir is None:
        return None
    gs = rp.gene_scores()
    npos = rp.gene_n_positions()
    gene_report = None
    prec_k = prec_k_npos = rank_pctile = math.nan
    ranks: dict[str, int] = {}
    if gs:
        gene_report = score_metrics(gs, planted, min(gs.values()), higher_is_hit=True)
        ranks = gene_report.positive_ranks
        n = len(gs)
        k = len(planted)
        topk = {g for g, _ in sorted(gs.items(), key=lambda kv: -kv[1])[:k]}
        prec_k = len(topk & planted) / k
        rank_pctile = _quantile([1.0 - (ranks[g] - 1) / max(1, n - 1) for g in planted], 0.50)
    if npos:
        k = len(planted)
        topk_n = {g for g, _ in sorted(npos.items(), key=lambda kv: -kv[1])[:k]}
        prec_k_npos = len(topk_n & planted) / k

    s5, s25 = rp.gene_slice("slice_global5"), rp.gene_slice("slice_global25")
    in_s5 = _mean([1.0 if g in s5 else 0.0 for g in planted]) if s5 else math.nan
    in_s25 = _mean([1.0 if g in s25 else 0.0 for g in planted]) if s25 else math.nan

    # -- sites: SCORING lists only detected positions + scheme_set + change_side  #
    detected: dict[str, dict[int, set[str]]] = {}
    side: dict[str, dict[int, str]] = {}
    for r in rp.position_rows():
        schemes = set((r.get("scheme_set") or "").replace(" ", "").split("+")) - {""}
        g, pos = r["Gene"], int(r["Position"])
        detected.setdefault(g, {})[pos] = schemes
        side.setdefault(g, {})[pos] = (r.get("change_side") or "").strip().lower()

    all_planted = hit = total_detected = dir_hit = 0
    mech_tot: dict[str, int] = {}
    mech_hit: dict[str, int] = {}
    ia_tot = ia_us = gc_tot = gc_gs = gc_leak = 0
    for g in planted:
        ps = rp.planted_sites(g)
        gdir = rp.gene_direction(g)
        det = detected.get(g, {})
        all_planted += len(ps)
        total_detected += len(det)
        hit += len(ps & set(det))
        dir_hit += sum(1 for s in ps if s in det and side.get(g, {}).get(s) == gdir)
        for s, m in rp.site_mechanism(g).items():
            mech_tot[m] = mech_tot.get(m, 0) + 1
            sch = det.get(s)
            if sch is not None:
                mech_hit[m] = mech_hit.get(m, 0) + 1
            if m == "identical_aa":
                ia_tot += 1
                ia_us += bool(sch and "US" in sch)
            elif m == "grouped_caap":
                gc_tot += 1
                if sch and (sch - {"US"}):
                    gc_gs += 1
                    gc_leak += bool("US" in sch)

    # -- contrast selection recovery: operative pairs vs the TRUE (lambda-drawn)
    #    pairs. traitfile.tab = species <tab> 1|0 <tab> pair_id ------------------ #
    tf = rp.traitfile()
    true_anchors = rp.foreground_tips
    true_partners = rp.partner_tips
    true_pairs = {frozenset(p) for p in rp.truth_pairs}
    by_pid: dict[str, dict[str, set[str]]] = {}
    for sp, lab, pid in tf:
        d = by_pid.setdefault(pid, {"1": set(), "0": set()})
        d.setdefault(lab, set()).add(sp)
    op_fg = {sp for d in by_pid.values() for sp in d["1"]}
    op_bg = {sp for d in by_pid.values() for sp in d["0"]}
    op_pairs = {frozenset(list(d["1"]) + list(d["0"]))
                for d in by_pid.values() if d["1"] and d["0"]}
    n_true = max(1, len(true_pairs))
    return ReplicateScore(
        subset=rp.subset, rep=rp.repdir.name, archetype=rp.archetype, lam=rp.lam,
        n_possible_pairs=rp.n_possible_pairs,
        n_genes_scored=len(gs), gene=gene_report,
        gene_precision_at_k=prec_k, gene_precision_at_k_npos=prec_k_npos,
        gene_planted_rank_pctile_p50=rank_pctile,
        gene_planted_ranks=ranks,
        planted_in_slice_global5=in_s5, planted_in_slice_global25=in_s25,
        site_recall=(hit / all_planted if all_planted else math.nan),
        site_precision=(hit / total_detected if total_detected else math.nan),
        site_recall_by_mechanism={m: mech_hit.get(m, 0) / n for m, n in mech_tot.items()},
        site_directional_recall=(dir_hit / all_planted if all_planted else math.nan),
        n_sites_planted=all_planted, n_sites_detected=total_detected,
        identical_aa_us_recall=(ia_us / ia_tot if ia_tot else math.nan),
        grouped_caap_gs_recall=(gc_gs / gc_tot if gc_tot else math.nan),
        grouped_caap_us_leakage=(gc_leak / gc_gs if gc_gs else math.nan),
        contrast_anchor_recall=len(true_anchors & op_fg) / n_true,
        contrast_partner_recall=len(true_partners & op_bg) / n_true,
        contrast_pair_exact=len(true_pairs & op_pairs) / n_true,
        contrast_jaccard=(
            len((true_anchors | true_partners) & (op_fg | op_bg))
            / len((true_anchors | true_partners) | (op_fg | op_bg))
            if (true_anchors or op_fg) else math.nan
        ),
        contrast_n_operative_pairs=len(op_pairs),
    )


# --------------------------------------------------------------------------- #
# score
# --------------------------------------------------------------------------- #
@dataclass
class CellSeparation:
    archetype: str
    lam: float
    n_power_reps: int
    n_null_reps: int
    # occurrence: detected CAAS positions per gene (planted ~n_planted, spurious ~1)
    auc_npos: float
    planted_npos_p50: float
    null_npos_p95: float
    null_npos_max: float
    # score-level (reported, not gated): position CAAS_score, hypergeometric p (as
    # significance), and the size-decorrelated gene_caas_score
    auc_pos_caas_score: float
    auc_pos_neglog_pvalue: float
    auc_gene_caas_score: float
    separated: bool                 # occurrence AUC >= 0.95 AND planted p50 > null p95

    def to_dict(self) -> dict:
        return dict(self.__dict__)


def _separation(archetype: str, lam: float, reps: list[Replicate]) -> CellSeparation | None:
    cell = [r for r in reps if r.archetype == archetype and abs(r.lam - lam) < 1e-6 and r.results_dir]
    power = [r for r in cell if not r.is_null]
    null = [r for r in cell if r.is_null]

    p_npos: list[float] = []
    p_gscore: list[float] = []
    p_pos_caas: list[float] = []
    p_pos_neglp: list[float] = []
    for r in power:
        npos, gsc = r.gene_n_positions(), r.gene_scores()
        planted = r.planted_genes
        psites = {g: r.planted_sites(g) for g in planted}
        p_npos += [npos[g] for g in planted if g in npos]
        p_gscore += [gsc[g] for g in planted if g in gsc]
        for row in r.position_rows():
            if int(row["Position"]) in psites.get(row["Gene"], set()):
                p_pos_caas.append(_num(row.get("CAAS_score", "nan")))
                pv = _num(row.get("pvalue", "nan"))
                p_pos_neglp.append(-math.log10(pv) if pv and pv > 0 else 30.0)

    n_npos: list[float] = []
    n_gscore: list[float] = []
    n_pos_caas: list[float] = []
    n_pos_neglp: list[float] = []
    for r in null:
        n_npos += list(r.gene_n_positions().values())
        n_gscore += list(r.gene_scores().values())
        for row in r.position_rows():
            n_pos_caas.append(_num(row.get("CAAS_score", "nan")))
            pv = _num(row.get("pvalue", "nan"))
            n_pos_neglp.append(-math.log10(pv) if pv and pv > 0 else 30.0)

    if not p_npos or not n_npos:
        return None
    auc_n = _auc(p_npos, n_npos)
    p50 = _quantile(p_npos, 0.50)
    n95 = _quantile(n_npos, 0.95)
    return CellSeparation(
        archetype=archetype, lam=lam,
        n_power_reps=len(power), n_null_reps=len(null),
        auc_npos=auc_n, planted_npos_p50=p50,
        null_npos_p95=n95, null_npos_max=max(n_npos),
        auc_pos_caas_score=_auc(p_pos_caas, n_pos_caas),
        auc_pos_neglog_pvalue=_auc(p_pos_neglp, n_pos_neglp),
        auc_gene_caas_score=_auc(p_gscore, n_gscore),
        separated=(not math.isnan(auc_n) and auc_n >= 0.95 and p50 > n95),
    )


@dataclass
class ScoreResult:
    n_power_replicates: int
    per_replicate: list[ReplicateScore] = field(default_factory=list)
    separation: list[CellSeparation] = field(default_factory=list)
    cells: dict = field(default_factory=dict)          # "archetype|lambda" -> metric means
    verdict: str = ""

    def to_dict(self) -> dict:
        return {
            "verdict": self.verdict,
            "n_power_replicates": self.n_power_replicates,
            "cells": self.cells,
            "separation": [s.to_dict() for s in self.separation],
            "per_replicate": [r.to_dict() for r in self.per_replicate],
        }


def score(run_root: Path) -> ScoreResult:
    reps = find_replicates(run_root)
    scored = [s for s in (_score_one(rp) for rp in reps) if s is not None]
    if not scored:
        raise SystemExit(
            f"no scorable power replicates under {run_root} ({len(reps)} replicate(s) found)"
        )

    cell_keys = sorted({(s.archetype, s.lam) for s in scored})
    seps = [x for x in (_separation(a, lm, reps) for a, lm in cell_keys) if x is not None]

    grp: dict[tuple[str, float], list[ReplicateScore]] = {}
    for s in scored:
        grp.setdefault((s.archetype, s.lam), []).append(s)

    cells: dict = {}
    for (a, lm), g in grp.items():
        mechs: dict[str, list[float]] = {}
        for s in g:
            for m, v in s.site_recall_by_mechanism.items():
                mechs.setdefault(m, []).append(v)
        cells[f"{a}|lambda={lm:g}"] = {
            "n_power_replicates": len(g),
            "n_possible_pairs": _mean([s.n_possible_pairs for s in g]),
            "gene_precision_at_k_npos": _mean([s.gene_precision_at_k_npos for s in g]),
            "gene_precision_at_k_score": _mean([s.gene_precision_at_k for s in g]),
            "gene_planted_rank_pctile_p50": _mean([s.gene_planted_rank_pctile_p50 for s in g]),
            "planted_in_slice_global25": _mean([s.planted_in_slice_global25 for s in g]),
            "site_recall": _mean([s.site_recall for s in g]),
            "site_directional_recall": _mean([s.site_directional_recall for s in g]),
            "site_recall_by_mechanism": {m: _mean(v) for m, v in mechs.items()},
            "identical_aa_us_recall": _mean([s.identical_aa_us_recall for s in g]),
            "grouped_caap_gs_recall": _mean([s.grouped_caap_gs_recall for s in g]),
            "grouped_caap_us_leakage": _mean([s.grouped_caap_us_leakage for s in g]),
            "contrast_anchor_recall": _mean([s.contrast_anchor_recall for s in g]),
            "contrast_partner_recall": _mean([s.contrast_partner_recall for s in g]),
            "contrast_pair_exact": _mean([s.contrast_pair_exact for s in g]),
            "contrast_jaccard": _mean([s.contrast_jaccard for s in g]),
            "contrast_n_operative_pairs": _mean([s.contrast_n_operative_pairs for s in g]),
        }

    # PASS = at EVERY (archetype, lambda) cell the pipeline manufactures no false
    # signal (occurrence separation holds) and finds the planted signal
    # (precision@k on CAAS count, site recall). Contrast recovery is REPORTED as a
    # lambda-curve, not gated — it is expected to degrade as the trait clumps.
    lo_lam = min(lm for _, lm in cell_keys)
    ok = bool(seps) and all(s.separated for s in seps) and all(
        not math.isnan(v["gene_precision_at_k_npos"]) and v["gene_precision_at_k_npos"] >= 0.85
        and not math.isnan(v["site_recall"]) and v["site_recall"] >= 0.8
        for k, v in cells.items() if k.endswith(f"lambda={lo_lam:g}")
    )
    verdict = "PASS" if ok else "REVIEW"

    return ScoreResult(
        n_power_replicates=len(scored),
        per_replicate=scored, separation=seps, cells=cells, verdict=verdict,
    )
