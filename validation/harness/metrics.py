"""Scoring and null-calibration metrics for the validation tiers.

Two families:

* ``score_metrics`` — given predicted scores and a known-positive set, at a
  threshold: precision / recall / F1 / MCC, plus threshold-free PR-AUC, ROC-AUC,
  and the rank of every known positive.
* ``null_calibration`` — given the permulation / bootstrap p-values from a
  no-signal run: KS distance to Uniform(0,1), observed type-I error at nominal
  alphas, and QQ points. This is the Tier 0 gate and the Tier 2 sanity check.

Pure stdlib + numpy (no scipy): the KS test uses the asymptotic
Kolmogorov distribution, which is fine for n >= ~50 null replicates.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np


# --------------------------------------------------------------------------- #
# score metrics
# --------------------------------------------------------------------------- #
@dataclass
class ScoreReport:
    n_pred_pos: int
    n_truth: int
    tp: int
    fp: int
    fn: int
    precision: float
    recall: float
    f1: float
    mcc: float
    pr_auc: float
    roc_auc: float
    positive_ranks: dict[str, int] = field(default_factory=dict)  # 1 = top
    positive_percentiles: dict[str, float] = field(default_factory=dict)

    def as_dict(self) -> dict:
        return self.__dict__


def score_metrics(
    scores: dict[str, float],
    truth: set[str],
    threshold: float,
    *,
    higher_is_hit: bool = True,
    n_negatives: int | None = None,
) -> ScoreReport:
    """
    scores        : id -> score (e.g. gene -> composite score, position -> -log10 p)
    truth         : ids that are true positives (must be a subset of scores' keys,
                    missing ones are counted as rank = len+1 / not recovered)
    threshold     : call a hit if (score >= threshold) when higher_is_hit else (<=)
    n_negatives   : universe size for ROC when the negative set is larger than
                    what appears in ``scores`` (rare; default uses scores' keys)
    """
    keys = list(scores)
    s = np.array([scores[k] for k in keys], dtype=float)
    if not higher_is_hit:
        s = -s
        threshold = -threshold
    y = np.array([1 if k in truth else 0 for k in keys], dtype=int)

    called = s >= threshold
    tp = int(np.sum(called & (y == 1)))
    fp = int(np.sum(called & (y == 0)))
    fn = int(np.sum(~called & (y == 1)))
    tn = int(np.sum(~called & (y == 0)))

    precision = tp / (tp + fp) if (tp + fp) else float("nan")
    recall = tp / (tp + fn) if (tp + fn) else float("nan")
    f1 = (2 * precision * recall / (precision + recall)
          if precision and recall and not math.isnan(precision) and not math.isnan(recall)
          else float("nan"))
    denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = ((tp * tn - fp * fn) / denom) if denom else float("nan")

    order = np.argsort(-s, kind="stable")  # best first
    ranked_keys = [keys[i] for i in order]
    rank_of = {k: r + 1 for r, k in enumerate(ranked_keys)}
    n = len(keys)
    pos_ranks = {k: rank_of.get(k, n + 1) for k in truth}
    pos_pct = {k: 1.0 - (rank_of.get(k, n + 1) - 1) / max(1, n - 1) for k in truth}

    return ScoreReport(
        n_pred_pos=tp + fp,
        n_truth=len(truth),
        tp=tp, fp=fp, fn=fn,
        precision=precision, recall=recall, f1=f1, mcc=mcc,
        pr_auc=_pr_auc(s, y),
        roc_auc=_roc_auc(s, y, n_negatives),
        positive_ranks=pos_ranks,
        positive_percentiles=pos_pct,
    )


def _roc_auc(s: np.ndarray, y: np.ndarray, n_negatives: int | None) -> float:
    pos = s[y == 1]
    neg = s[y == 0]
    if len(pos) == 0 or (len(neg) == 0 and not n_negatives):
        return float("nan")
    # Mann–Whitney U statistic / (n_pos * n_neg), ties = 0.5
    wins = 0.0
    for p in pos:
        wins += np.sum(neg < p) + 0.5 * np.sum(neg == p)
    nn = n_negatives if n_negatives else len(neg)
    return wins / (len(pos) * nn)


def _pr_auc(s: np.ndarray, y: np.ndarray) -> float:
    if np.sum(y) == 0:
        return float("nan")
    order = np.argsort(-s, kind="stable")
    y = y[order]
    tp = np.cumsum(y)
    fp = np.cumsum(1 - y)
    precision = tp / np.maximum(tp + fp, 1)
    recall = tp / np.sum(y)
    # step integration over recall
    rec = np.concatenate([[0.0], recall])
    prec = np.concatenate([[1.0], precision])
    return float(np.sum(np.diff(rec) * prec[1:]))


# --------------------------------------------------------------------------- #
# null calibration
# --------------------------------------------------------------------------- #
@dataclass
class NullReport:
    n: int
    ks_stat: float
    ks_pvalue: float          # H0: p-values ~ Uniform(0,1); small => miscalibrated
    type1_error: dict[float, float]  # nominal alpha -> observed fraction p <= alpha
    mean: float
    frac_exact_zero: float
    frac_exact_one: float
    qq: list[tuple[float, float]] = field(default_factory=list)  # (expected, observed)

    def verdict(self, alpha: float = 0.05) -> str:
        if self.ks_pvalue < 0.05:
            return "MISCALIBRATED (KS rejects uniformity)"
        obs = self.type1_error.get(alpha)
        if obs is not None and obs > 2 * alpha:
            return f"INFLATED type-I error at alpha={alpha}: {obs:.3f}"
        if obs is not None and obs < 0.3 * alpha:
            return f"CONSERVATIVE at alpha={alpha}: {obs:.3f}"
        return "OK"


def null_calibration(pvalues, alphas=(0.10, 0.05, 0.01), n_qq: int = 50) -> NullReport:
    p = np.asarray([x for x in pvalues if x is not None and not np.isnan(x)], dtype=float)
    n = len(p)
    if n == 0:
        raise ValueError("no finite p-values")
    p_sorted = np.sort(p)
    cdf_emp = np.arange(1, n + 1) / n
    d_plus = np.max(cdf_emp - p_sorted)
    d_minus = np.max(p_sorted - (np.arange(0, n) / n))
    ks = float(max(d_plus, d_minus))

    t1 = {a: float(np.mean(p <= a)) for a in alphas}
    qi = np.linspace(0, 1, n_qq)
    qq = [(float(q), float(np.quantile(p, q))) for q in qi]

    return NullReport(
        n=n,
        ks_stat=ks,
        ks_pvalue=_kolmogorov_sf(ks * math.sqrt(n)),
        type1_error=t1,
        mean=float(np.mean(p)),
        frac_exact_zero=float(np.mean(p == 0.0)),
        frac_exact_one=float(np.mean(p == 1.0)),
        qq=qq,
    )


def _kolmogorov_sf(x: float, terms: int = 100) -> float:
    """P(sqrt(n)*D > x) — asymptotic Kolmogorov distribution survival function."""
    if x <= 0:
        return 1.0
    s = sum((-1) ** (k - 1) * math.exp(-2 * (k * x) ** 2) for k in range(1, terms + 1))
    return max(0.0, min(1.0, 2 * s))
