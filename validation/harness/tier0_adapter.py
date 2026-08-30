"""Adapter: a Tier 0 ``runs/`` tree -> gate metrics.

A staged Tier 0 run (``validation.tier0.run_replicates``) lays out

    <run_root>/<archetype>_<set>_<tree>/rep###/
        truth.json
        out/<trait>_complete/scoring/{gene_scores,position_scores}.tsv
        out/<trait>_complete/scoring/gene_lists/slice_*.tsv
        out/<trait>_complete/data_exploration/2.CT/1.Traitfiles/traitfile.tab

The pipeline is a **prioritisation engine** (see project memory
``project_tier0_scoping_reframe``): its permulation p-values are conditional
foreground-specificity scores, never Uniform(0,1). So the gate is NOT a KS
uniformity test. It is:

  1. prioritisation  — planted genes/positions in ``slice_global1`` / ``slice_global5``
  2. null vs power   — pooled per archetype: planted-gene ``gene_caas_score`` in
     power replicates must separate (AUC) from every gene's score in the matched
     null replicates. If a no-signal replicate reaches the same scores, the
     ranking cannot tell signal from noise.
  3. site recovery   — planted-site detection recall, split by mechanism and by
     which CAAP scheme fired (``identical_aa`` -> US, ``grouped_caap`` -> GS only).
  4. contrast recovery — operative foreground from ``traitfile.tab`` vs the
     planted pairs in ``truth``.

``truth["genes"][g]["planted_sites"]`` keys are 0-based columns; the pipeline's
``Position`` column matches them directly (verified: 12/12 exact, 0/12 at +1).
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

    def gene_scores(self) -> dict[str, float]:
        p = self._out("scoring/gene_scores.tsv")
        if p is None:
            return {}
        out = {r["Gene"]: _num(r.get("gene_caas_score", "nan")) for r in _read_tsv(p)}
        return {g: v for g, v in out.items() if not math.isnan(v)}

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
    gene_precision_at_k: float
    gene_planted_ranks: dict[str, int]
    planted_in_slice_global1: float     # fraction of planted genes in the top-1% slice
    planted_in_slice_global5: float
    site_recall: float
    site_precision: float
    site_recall_by_mechanism: dict[str, float]
    n_sites_planted: int
    n_sites_detected: int
    identical_aa_us_recall: float
    grouped_caap_gs_recall: float
    grouped_caap_us_leakage: float
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
    if not planted or rp.results_dir is None:
        return None
    gs = rp.gene_scores()
    gene_report = None
    prec_k = math.nan
    ranks: dict[str, int] = {}
    if gs:
        gene_report = score_metrics(gs, planted, min(gs.values()), higher_is_hit=True)
        ranks = gene_report.positive_ranks
        k = len(planted)
        topk = {g for g, _ in sorted(gs.items(), key=lambda kv: -kv[1])[:k]}
        prec_k = len(topk & planted) / k

    s1, s5 = rp.gene_slice("slice_global1"), rp.gene_slice("slice_global5")
    in_s1 = _mean([1.0 if g in s1 else 0.0 for g in planted]) if s1 else math.nan
    in_s5 = _mean([1.0 if g in s5 else 0.0 for g in planted]) if s5 else math.nan

    # -- sites: SCORING lists only detected positions + their scheme_set -------- #
    detected: dict[str, dict[int, set[str]]] = {}
    for r in rp.position_rows():
        schemes = set((r.get("scheme_set") or "").replace(" ", "").split("+")) - {""}
        detected.setdefault(r["Gene"], {})[int(r["Position"])] = schemes

    all_planted = hit = total_detected = 0
    mech_tot: dict[str, int] = {}
    mech_hit: dict[str, int] = {}
    ia_tot = ia_us = gc_tot = gc_gs = gc_leak = 0
    for g in planted:
        ps = rp.planted_sites(g)
        det = detected.get(g, {})
        all_planted += len(ps)
        total_detected += len(det)
        hit += len(ps & set(det))
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

    # -- contrast selection recovery ---------------------------------------- #
    tf = rp.traitfile()
    fg_truth = rp.foreground_tips
    op_fg = {sp for sp, lab, _ in tf if lab == "1"}
    jac = (len(op_fg & fg_truth) / len(op_fg | fg_truth)) if (op_fg or fg_truth) else math.nan
    fg_prec = (len(op_fg & fg_truth) / len(op_fg)) if op_fg else math.nan
    by_pid: dict[str, list[tuple[str, str]]] = {}
    for sp, lab, pid in tf:
        by_pid.setdefault(pid, []).append((sp, lab))
    recovered = sum(
        1 for members in by_pid.values()
        if any(lab == "1" and sp in fg_truth for sp, lab in members)
    )

    return ReplicateScore(
        subset=rp.subset, rep=rp.repdir.name, archetype=rp.archetype,
        n_genes_scored=len(gs), gene=gene_report,
        gene_precision_at_k=prec_k, gene_planted_ranks=ranks,
        planted_in_slice_global1=in_s1, planted_in_slice_global5=in_s5,
        site_recall=(hit / all_planted if all_planted else math.nan),
        site_precision=(hit / total_detected if total_detected else math.nan),
        site_recall_by_mechanism={m: mech_hit.get(m, 0) / n for m, n in mech_tot.items()},
        n_sites_planted=all_planted, n_sites_detected=total_detected,
        identical_aa_us_recall=(ia_us / ia_tot if ia_tot else math.nan),
        grouped_caap_gs_recall=(gc_gs / gc_tot if gc_tot else math.nan),
        grouped_caap_us_leakage=(gc_leak / gc_gs if gc_gs else math.nan),
        contrast_jaccard=jac, contrast_fg_precision=fg_prec,
        contrast_pairs_recovered=recovered, contrast_pairs_total=len(by_pid),
    )


# --------------------------------------------------------------------------- #
# score
# --------------------------------------------------------------------------- #
@dataclass
class ArchetypeSeparation:
    archetype: str
    n_power_reps: int
    n_null_reps: int
    n_planted_scores: int
    n_null_scores: int
    auc_planted_vs_null: float      # 1.0 = perfect separation
    planted_score_p50: float
    null_score_p95: float
    null_score_max: float
    separated: bool                 # planted p50 > null p95 AND auc >= 0.9

    def to_dict(self) -> dict:
        return dict(self.__dict__)


def _separation(archetype: str, reps: list[Replicate]) -> ArchetypeSeparation | None:
    power = [r for r in reps if r.archetype == archetype and not r.is_null and r.results_dir]
    null = [r for r in reps if r.archetype == archetype and r.is_null and r.results_dir]
    planted_scores: list[float] = []
    for r in power:
        gs = r.gene_scores()
        planted_scores += [gs[g] for g in r.planted_genes if g in gs]
    null_scores: list[float] = []
    for r in null:
        null_scores += list(r.gene_scores().values())
    if not planted_scores or not null_scores:
        return None
    auc = _auc(planted_scores, null_scores)
    p50 = _quantile(planted_scores, 0.50)
    n95 = _quantile(null_scores, 0.95)
    return ArchetypeSeparation(
        archetype=archetype,
        n_power_reps=len(power), n_null_reps=len(null),
        n_planted_scores=len(planted_scores), n_null_scores=len(null_scores),
        auc_planted_vs_null=auc,
        planted_score_p50=p50, null_score_p95=n95,
        null_score_max=max(null_scores),
        separated=(not math.isnan(auc) and auc >= 0.9 and p50 > n95),
    )


@dataclass
class ScoreResult:
    n_power_replicates: int
    per_replicate: list[ReplicateScore] = field(default_factory=list)
    separation: list[ArchetypeSeparation] = field(default_factory=list)
    summary: dict = field(default_factory=dict)
    verdict: str = ""

    def to_dict(self) -> dict:
        return {
            "verdict": self.verdict,
            "n_power_replicates": self.n_power_replicates,
            "summary": self.summary,
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

    archetypes = sorted({s.archetype for s in scored})
    seps = [s for s in (_separation(a, reps) for a in archetypes) if s is not None]

    by_arch: dict[str, list[ReplicateScore]] = {}
    for s in scored:
        by_arch.setdefault(s.archetype, []).append(s)

    summary: dict = {}
    for a, group in by_arch.items():
        mechs: dict[str, list[float]] = {}
        for s in group:
            for m, v in s.site_recall_by_mechanism.items():
                mechs.setdefault(m, []).append(v)
        summary[a] = {
            "n_power_replicates": len(group),
            "gene_precision_at_k": _mean([s.gene_precision_at_k for s in group]),
            "gene_roc_auc": _mean([s.gene.roc_auc for s in group if s.gene]),
            "planted_in_slice_global1": _mean([s.planted_in_slice_global1 for s in group]),
            "planted_in_slice_global5": _mean([s.planted_in_slice_global5 for s in group]),
            "site_recall": _mean([s.site_recall for s in group]),
            "site_recall_by_mechanism": {m: _mean(v) for m, v in mechs.items()},
            "identical_aa_us_recall": _mean([s.identical_aa_us_recall for s in group]),
            "grouped_caap_gs_recall": _mean([s.grouped_caap_gs_recall for s in group]),
            "grouped_caap_us_leakage": _mean([s.grouped_caap_us_leakage for s in group]),
            "contrast_jaccard": _mean([s.contrast_jaccard for s in group]),
            "contrast_pairs_recovered_frac": _mean(
                [s.contrast_pairs_recovered / s.contrast_pairs_total
                 for s in group if s.contrast_pairs_total]
            ),
        }

    ok = bool(seps) and all(s.separated for s in seps) and all(
        (not math.isnan(v["planted_in_slice_global5"]) and v["planted_in_slice_global5"] >= 0.8)
        for v in summary.values()
    )
    verdict = "PASS" if ok else "REVIEW"

    return ScoreResult(
        n_power_replicates=len(scored),
        per_replicate=scored, separation=seps, summary=summary, verdict=verdict,
    )
