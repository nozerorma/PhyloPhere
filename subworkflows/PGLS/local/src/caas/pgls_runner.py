#!/usr/bin/env python3
# ================================ src/caas/pgls_runner.py ================================

from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, Set

import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize_scalar
from scipy.stats import chi2, hypergeom

from stats.stats_utils import bh_fdr, bh_fdr_per_group

RESULT_COLUMNS = [
    "Gene",
    "Position",
    "trait",
    "CAAP_Group",
    "change_side",
    "is_caas_position",
    "n_used",
    "n_top_aa",
    "n_bottom_aa",
    "beta_top",
    "lambda_full",
    "lambda_null",
    "lrt_stat",
    "p_pgls_two_sided",
    "p_pgls_pos",       # one-tailed: H1 = beta_top > 0  (β positive; top-group AA enriched in high-trait species)
    "p_pgls_neg",       # one-tailed: H1 = beta_top < 0  (β negative; top-group AA enriched in low-trait species)
    "q_p_pgls_two_sided",
    "sig_p_pgls_two_sided",
    "q_p_pgls_pos",
    "sig_p_pgls_pos",
    "q_p_pgls_neg",
    "sig_p_pgls_neg",
    "p_pgls_directional",
    "q_p_pgls_directional",
    "sig_p_pgls_directional",
    "pgls_decile_score",
    "pgls_weighted_decile",
    "pgls_char_score",
]

EXCESS_COLUMNS = [
    "trait",
    "Gene",
    "CAAP_Group",
    "N_testable",
    "K_sig_testable",
    "n_testable_caas",
    "k_sig_caas",
    "rate_testable",
    "rate_caas",
    "excess_ratio",
    "p_hyper_excess",
    "q_hyper_excess",
    "sig_hyper_excess",
]

DIAG_COLUMNS = [
    "Gene",
    "Position",
    "trait",
    "mapped_tips",
    "n_top_aa",
    "n_bottom_aa",
    "trait_model",
    "trait_theta",
    "diag_type",
    "model",
    "k",
    "theta",
    "ll",
    "aicc",
    "lrt_df",
    "lrt_stat",
    "lrt_p",
    "is_best_aicc",
]

SCHEME_WEIGHTS = {
    "US": 0.5,
    "GS4": 0.2,
    "GS3": 0.1,
    "GS2": 0.1,
    "GS1": 0.1,
    "GS0": 0.0,
}


def _decile_score_from_pvalues(values: pd.Series) -> pd.Series:
    """Map p-values to the same 0..1 decile scale used in scoring."""
    vals = pd.to_numeric(values, errors="coerce")
    out = pd.Series(np.nan, index=vals.index, dtype=float)
    ok = vals.notna()
    n = int(ok.sum())
    if n == 0:
        return out
    if n == 1:
        out.loc[ok] = 1.0
        return out
    ranks = vals.loc[ok].rank(method="first", ascending=True)
    ntile = np.ceil(ranks * 10.0 / n).astype(int).clip(lower=1, upper=10)
    out.loc[ok] = 1.0 - (ntile - 1) / 9.0
    return out


def _symmetrize(V: np.ndarray) -> np.ndarray:
    V = np.asarray(V, float)
    return 0.5 * (V + V.T)


def _ensure_psd(V: np.ndarray, eps: float = 1e-8) -> np.ndarray:
    V = _symmetrize(np.asarray(V, float))
    if V.ndim != 2 or V.shape[0] != V.shape[1]:
        raise ValueError("VCV matrix must be square")
    if V.size == 0:
        return V
    try:
        np.linalg.cholesky(V)
        return V
    except np.linalg.LinAlgError:
        pass
    eigvals = np.linalg.eigvalsh(V)
    min_eig = float(np.min(eigvals))
    if min_eig < eps:
        V = V + np.eye(V.shape[0]) * (eps - min_eig)
    try:
        np.linalg.cholesky(V)
        return V
    except np.linalg.LinAlgError:
        jitter = max(eps, 1e-6)
        for _ in range(8):
            candidate = V + np.eye(V.shape[0]) * jitter
            try:
                np.linalg.cholesky(candidate)
                return candidate
            except np.linalg.LinAlgError:
                jitter *= 10.0
    raise np.linalg.LinAlgError("Unable to regularize covariance matrix")


def _gls_fit(y, X, V, ridge=1e-8):
    import numpy.linalg as la

    y = np.asarray(y, float)
    X = np.asarray(X, float)
    V = np.asarray(V, float)
    n = len(y)
    A = _ensure_psd(V, eps=max(ridge, 1e-10))
    c, lower = cho_factor(A, check_finite=False)
    ViX = cho_solve((c, lower), X, check_finite=False)
    Viy = cho_solve((c, lower), y, check_finite=False)
    XtViX = X.T @ ViX
    XtViy = X.T @ Viy
    try:
        beta = la.solve(XtViX, XtViy)
    except la.LinAlgError:
        beta = la.pinv(XtViX) @ XtViy
    resid = y - X @ beta
    sse = float(resid.T @ cho_solve((c, lower), resid, check_finite=False))
    diagL = np.diag(c)
    logdetA = 2.0 * np.sum(np.log(np.abs(diagL)))
    s2_ml = sse / n
    ll = -0.5 * (n * np.log(2.0 * np.pi) + logdetA + n * np.log(s2_ml) + sse / s2_ml)
    return dict(beta=beta, ll=ll)


def _lrt(y, A, X_full, X_null):
    fit_full = _gls_fit(y, X_full, A)
    fit_null = _gls_fit(y, X_null, A)
    df = X_full.shape[1] - X_null.shape[1]
    stat = max(0.0, 2.0 * (fit_full["ll"] - fit_null["ll"]))
    p = chi2.sf(stat, df)
    return stat, df, p, fit_full


def _one_tailed_top(two_sided_p: float, beta_top: float) -> float:
    """One-tailed p-value for H1: beta_top > 0 (amino acid enriched in high-trait / top species)."""
    if not np.isfinite(two_sided_p) or not np.isfinite(beta_top):
        return np.nan
    if beta_top <= 0:
        return 1.0
    return min(1.0, max(0.0, two_sided_p / 2.0))


def _one_tailed_bottom(two_sided_p: float, beta_top: float) -> float:
    """One-tailed p-value for H1: beta_top < 0 (amino acid enriched in low-trait / bottom species)."""
    if not np.isfinite(two_sided_p) or not np.isfinite(beta_top):
        return np.nan
    if beta_top >= 0:
        return 1.0
    return min(1.0, max(0.0, two_sided_p / 2.0))


def _precompute_tree_caches(tree):
    root = tree.root
    terms = list(tree.get_terminals())
    name2node = {str(t.name): t for t in terms}
    paths = {str(t.name): [root, *tree.get_path(t)] for t in terms}

    def _edge_len(parent, child):
        bl = getattr(child, "branch_length", None)
        return float(bl) if bl is not None else 0.0

    depth = {root: 0.0}
    for tip, path in paths.items():
        acc = 0.0
        prev = root
        for node in path[1:]:
            acc += _edge_len(prev, node)
            if node not in depth:
                depth[node] = acc
            prev = node
        depth[name2node[tip]] = acc
    H = max(depth.get(t, 0.0) for t in terms) if terms else 0.0
    height = {node: (H - d) for node, d in depth.items()}
    parent = {}
    for path in paths.values():
        for p, c in zip(path[:-1], path[1:]):
            parent[c] = p
    return dict(root=root, name2node=name2node, paths=paths, depth=depth, height=height, parent=parent)


def _mrca_depth_nm(name_i, name_j, paths, depth):
    pi = paths[name_i]
    pj = paths[name_j]
    m = min(len(pi), len(pj))
    k = 0
    while k < m and pi[k] is pj[k]:
        k += 1
    mrca = pi[k - 1]
    return depth[mrca]


def _vcv_bm_from_cache(tip_names: List[str], cache) -> np.ndarray:
    n = len(tip_names)
    V = np.zeros((n, n), dtype=float)
    depth = cache["depth"]
    paths = cache["paths"]
    name2node = cache["name2node"]
    di = [depth[name2node[nm]] for nm in tip_names]
    for i in range(n):
        V[i, i] = di[i]
        for j in range(i + 1, n):
            dm = _mrca_depth_nm(tip_names[i], tip_names[j], paths, depth)
            V[i, j] = V[j, i] = dm
    return V


def _patristic_D_from_cache(tip_names: List[str], cache) -> np.ndarray:
    n = len(tip_names)
    D = np.zeros((n, n), dtype=float)
    depth = cache["depth"]
    paths = cache["paths"]
    name2node = cache["name2node"]
    di = [depth[name2node[nm]] for nm in tip_names]
    for i in range(n):
        for j in range(i + 1, n):
            dm = _mrca_depth_nm(tip_names[i], tip_names[j], paths, depth)
            dij = di[i] + di[j] - 2.0 * dm
            D[i, j] = D[j, i] = dij
    return D


def _apply_lambda(V: np.ndarray, lam: float) -> np.ndarray:
    lam = float(np.clip(lam, 0.0, 1.0))
    D = np.diag(np.diag(V))
    return D + lam * (V - D)


def _vcv_ou_from_cache(tip_names: List[str], cache, alpha: float) -> np.ndarray:
    alpha = float(max(alpha, 1e-12))
    D = _patristic_D_from_cache(tip_names, cache)
    C = np.exp(-alpha * D)
    Vbm = _vcv_bm_from_cache(tip_names, cache)
    diag = np.diag(Vbm)
    S = np.sqrt(np.maximum(diag, 1e-12))
    V = C * (S[:, None] * S[None, :])
    np.fill_diagonal(V, diag)
    return V


def _transformed_depths(cache, model: str, theta: float) -> Dict[object, float]:
    depth = cache["depth"]
    height = cache["height"]
    parent = cache["parent"]
    root = cache["root"]

    def elen(p, c):
        bl = max(0.0, depth[c] - depth[p])
        if model == "kappa":
            return max(0.0, float(bl) ** float(theta))
        if model == "delta":
            return max(0.0, float(height[p]) ** float(theta) - float(height[c]) ** float(theta))
        if model == "eb":
            return max(0.0, float(bl) * np.exp(float(theta) * float(height[p])))
        return bl

    dprime = {root: 0.0}
    changed = True
    while changed:
        changed = False
        for c, p in parent.items():
            if p in dprime and c not in dprime:
                dprime[c] = dprime[p] + elen(p, c)
                changed = True
    return dprime


def _vcv_from_model(tip_names: List[str], cache, model: str, theta: float | None) -> np.ndarray:
    if model == "bm" or theta is None:
        V = _vcv_bm_from_cache(tip_names, cache)
    elif model == "lambda":
        V = _apply_lambda(_vcv_bm_from_cache(tip_names, cache), float(theta))
    elif model == "ou":
        V = _vcv_ou_from_cache(tip_names, cache, float(theta))
    else:
        dprime = _transformed_depths(cache, model, float(theta))
        n = len(tip_names)
        V = np.zeros((n, n), dtype=float)
        paths = cache["paths"]
        name2node = cache["name2node"]
        di = [dprime[name2node[nm]] for nm in tip_names]
        for i in range(n):
            V[i, i] = di[i]
            for j in range(i + 1, n):
                pi = paths[tip_names[i]]
                pj = paths[tip_names[j]]
                m = min(len(pi), len(pj))
                k = 0
                while k < m and pi[k] is pj[k]:
                    k += 1
                mrca = pi[k - 1]
                V[i, j] = V[j, i] = dprime[mrca]
    return _ensure_psd(V)


def _theta_bounds(model: str, bounds: str) -> tuple[float, float]:
    if bounds == "phylolm":
        return (0.01, 0.99) if model == "lambda" else (1e-3, 5.0) if model in {"delta", "kappa"} else (1e-3, 3.0)
    return (0.0, 1.0) if model == "lambda" else (1e-4, 10.0)


def _fit_site_model(y: np.ndarray,
                    X: np.ndarray,
                    tips: List[str],
                    cache,
                    model: str,
                    theta: float | None,
                    bounds: str) -> dict:
    if model == "lambda" and theta is None:
        lo, hi = _theta_bounds("lambda", bounds)

        def obj(th):
            try:
                V = _vcv_from_model(tips, cache, "lambda", th)
                return -_gls_fit(y, X, V)["ll"]
            except Exception:
                return np.inf

        opt = minimize_scalar(obj, bounds=(lo, hi), method="bounded", options={"xatol": 1e-3})
        if not opt.success or not np.isfinite(opt.x):
            raise RuntimeError("Failed to fit lambda for site-level PGLS")
        theta_hat = float(opt.x)
    else:
        theta_hat = theta

    V = _vcv_from_model(tips, cache, model, theta_hat)
    fit = _gls_fit(y, X, V)
    fit["theta"] = theta_hat
    fit["model"] = model
    return fit


def _aicc(ll: float, k: int, n: int) -> float:
    denom = max(n - k - 1, 1)
    return -2.0 * ll + 2.0 * k + (2.0 * k * (k + 1)) / denom


def _site_lrt(y: np.ndarray,
              X_full: np.ndarray,
              X_null: np.ndarray,
              tips: List[str],
              cache,
              model: str,
              theta: float | None,
              bounds: str) -> dict:
    n = len(y)
    full = _fit_site_model(y, X_full, tips, cache, model, theta, bounds)
    null = _fit_site_model(y, X_null, tips, cache, model, theta, bounds)
    df = X_full.shape[1] - X_null.shape[1]
    stat = max(0.0, 2.0 * (full["ll"] - null["ll"]))
    # chi2.sf can underflow to exact 0.0 for large statistics; floor at the
    # smallest representable positive float so downstream log/BH-FDR stay valid.
    p_two = max(float(chi2.sf(stat, df)), np.finfo(float).tiny)
    beta_top = float(full["beta"][1]) if len(full["beta"]) > 1 else np.nan
    # k_full = regression params (X_full cols) + 1 for sigma^2, + 1 if theta fitted
    k_full = X_full.shape[1] + 1 + (1 if (model != "bm" and theta is None) else 0)
    return {
        "lrt_stat": float(stat),
        "lrt_df": int(df),
        "p_two_sided": p_two,
        "p_one_tailed_top":    float(_one_tailed_top(p_two, beta_top)),
        "p_one_tailed_bottom": float(_one_tailed_bottom(p_two, beta_top)),
        "beta_top": beta_top,
        "lambda_full": full.get("theta", np.nan),
        "lambda_null": null.get("theta", np.nan),
        # per-site diagnostics
        "diag_model":  full.get("model", model),
        "diag_theta":  full.get("theta", np.nan),
        "diag_ll":     float(full["ll"]),
        "diag_k":      k_full,
        "diag_aicc":   _aicc(float(full["ll"]), k_full, n),
    }


def _directional_summary(change_side: str,
                         p_top: float,
                         p_bottom: float,
                         q_top: float,
                         q_bottom: float,
                         sig_top: bool,
                         sig_bottom: bool) -> tuple[float, float, bool]:
    side = str(change_side or "none").strip().lower()
    if side == "top":
        return p_top, q_top, bool(sig_top)
    if side == "bottom":
        # x = (observed_dir == "top"), so beta_top > 0 for a real bottom-CAAS signal.
        # The correct one-sided test is p_top (H1: beta_top > 0), same as for top changes.
        return p_top, q_top, bool(sig_top)
    if side == "both":
        p_vals = [v for v in [p_top, p_bottom] if np.isfinite(v)]
        q_vals = [v for v in [q_top, q_bottom] if np.isfinite(v)]
        return (
            float(min(p_vals)) if p_vals else np.nan,
            float(min(q_vals)) if q_vals else np.nan,
            bool(sig_top) or bool(sig_bottom),
        )
    return np.nan, np.nan, False


def _coalesce_change_side(series: pd.Series) -> str:
    vals = {str(v).strip().lower() for v in series.dropna().astype(str)}
    if "both" in vals or ({"top", "bottom"}.issubset(vals)):
        return "both"
    if "top" in vals:
        return "top"
    if "bottom" in vals:
        return "bottom"
    return "none"


def compute_excess_by_gene_trait(res_df: pd.DataFrame, fdr_alpha: float = 0.05) -> pd.DataFrame:
    if res_df.empty:
        return pd.DataFrame(columns=EXCESS_COLUMNS)

    tested = res_df[res_df["n_used"].fillna(0).astype(float) > 0].copy()
    if tested.empty:
        return pd.DataFrame(columns=EXCESS_COLUMNS)

    group_cols = ["trait", "Gene", "CAAP_Group"]
    rows = []
    for key, grp in tested.groupby(group_cols, dropna=False):
        trait, gene, caap_group = key
        N = int(len(grp))
        K = int((grp["sig_p_pgls_two_sided"] == True).sum())
        n = int((grp["is_caas_position"] == True).sum())
        k = int(((grp["is_caas_position"] == True) & (grp["sig_p_pgls_two_sided"] == True)).sum())

        rate_testable = (K / N) if N > 0 else np.nan
        rate_caas = (k / n) if n > 0 else np.nan
        excess_ratio = (rate_caas / rate_testable) if np.isfinite(rate_caas) and np.isfinite(rate_testable) and rate_testable > 0 else np.nan
        p_hyper = float(hypergeom.sf(k - 1, N, K, n)) if (N > 0 and n > 0) else np.nan

        rows.append({
            "trait": trait,
            "Gene": gene,
            "CAAP_Group": caap_group,
            "N_testable": N,
            "K_sig_testable": K,
            "n_testable_caas": n,
            "k_sig_caas": k,
            "rate_testable": rate_testable,
            "rate_caas": rate_caas,
            "excess_ratio": excess_ratio,
            "p_hyper_excess": p_hyper,
        })

    out = pd.DataFrame(rows)
    out["q_hyper_excess"] = np.nan
    out["sig_hyper_excess"] = False
    for trait, grp in out.groupby("trait"):
        idx = grp.index
        qvals = bh_fdr(np.asarray(out.loc[idx, "p_hyper_excess"].values, float))
        out.loc[idx, "q_hyper_excess"] = qvals
        out.loc[idx, "sig_hyper_excess"] = qvals <= float(fdr_alpha)
    return out.reindex(columns=EXCESS_COLUMNS)


def _first_tip_match(
    species: str,
    tip_set: Set[str],
    sp2tax: Dict[str, str] | None,
    taxid2names: Dict[str, Set[str]] | None = None,
) -> str | None:
    s = str(species).strip()
    if s in tip_set:
        return s
    if sp2tax is not None:
        tid = str(sp2tax.get(s, "")).strip()
        if tid and tid in tip_set:
            return tid
        # Tree tips are binomial names, not taxon IDs.  Walk every known name
        # for this taxon (scientific name + all synonyms) and return the first
        # that matches a tree tip.
        if tid and taxid2names:
            for name in taxid2names.get(tid, set()):
                if name in tip_set:
                    return name
    return None


def _trait_model_select(trait: str,
                        y_map: Dict[str, float],
                        tree,
                        tips_all: List[str],
                        candidates: List[str],
                        bounds: str) -> dict:
    y = np.array([y_map[t] for t in tips_all if t in y_map], float)
    tips = [t for t in tips_all if t in y_map]
    n = len(tips)
    if n < 3:
        res = dict(model="bm", theta=None, ll=np.nan, aicc=np.inf, converged=False, k=1)
        return {"best": res, "fits": [res]}

    cache = _precompute_tree_caches(tree)
    fits = []

    def _aicc(ll, k, n_obs):
        return -2.0 * ll + 2.0 * k + (2.0 * k * (k + 1)) / max(n_obs - k - 1, 1)

    if "bm" in candidates:
        V = _ensure_psd(_vcv_bm_from_cache(tips, cache))
        fit = _gls_fit(y, np.ones((n, 1)), V)
        fits.append(dict(model="bm", theta=None, ll=fit["ll"], aicc=_aicc(fit["ll"], 1, n), converged=True, k=1))

    def _bounds(model):
        if bounds == "phylolm":
            return (0.01, 0.99) if model == "lambda" else (1e-3, 5.0) if model in {"delta", "kappa"} else (1e-3, 3.0)
        return (0.0, 1.0) if model == "lambda" else (1e-4, 10.0)

    for model in [c for c in candidates if c != "bm"]:
        lo, hi = _bounds(model)

        def obj(th):
            try:
                V = _ensure_psd(_vcv_from_model(tips, cache, model, th))
                return -_gls_fit(y, np.ones((n, 1)), V)["ll"]
            except Exception:
                return np.inf

        res = minimize_scalar(obj, bounds=(lo, hi), method="bounded", options={"xatol": 1e-3})
        if res.success and np.isfinite(res.x):
            theta = float(res.x)
            V = _ensure_psd(_vcv_from_model(tips, cache, model, theta))
            fit = _gls_fit(y, np.ones((n, 1)), V)
            fits.append(dict(model=model, theta=theta, ll=fit["ll"], aicc=_aicc(fit["ll"], 2, n), converged=True, k=2))
        else:
            fits.append(dict(model=model, theta=np.nan, ll=-np.inf, aicc=np.inf, converged=False, k=2))

    best = min(fits, key=lambda d: d["aicc"])
    logging.info(f"[TRAIT trait={trait}] best={best['model']} theta={best.get('theta')} AICc={best['aicc']:.3f}")
    return {"best": best, "fits": fits}


def _lrt_bm_vs_alt(bm_fit: dict, alt_fit: dict) -> tuple[float, int, float]:
    if not (bm_fit.get("converged") and alt_fit.get("converged")):
        return (np.nan, 0, np.nan)
    ll0, k0 = bm_fit["ll"], bm_fit["k"]
    ll1, k1 = alt_fit["ll"], alt_fit["k"]
    if not (np.isfinite(ll0) and np.isfinite(ll1)) or k1 <= k0:
        return (np.nan, 0, np.nan)
    stat = max(0.0, 2.0 * (ll1 - ll0))
    df = int(k1 - k0)
    if df == 1:
        p = 0.5 * (1.0 - chi2.cdf(stat, 1))
    else:
        p = 1.0 - chi2.cdf(stat, df)
    return (float(stat), df, float(p))


def run_caas_pgls(valid_df: pd.DataFrame,
                  tree,
                  tip_set: Set[str],
                  sp2tax: Dict[str, str] | None,
                  taxid2names: Dict[str, Set[str]] | None = None,
                  n_map: Dict[str, float] | None = None,
                  min_per_class: int = 4,
                  select_cfg: dict = None,
                  threads: int = 1,
                  fdr_alpha: float = 0.05) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    tips_all = [str(t.name) for t in tree.get_terminals()]
    traits = sorted(valid_df["trait"].unique())

    best_for_trait: Dict[str, dict] = {}
    trait_diag_rows: List[dict] = []
    for trait in traits:
        sub = (
            valid_df[valid_df["trait"] == trait][["species", "trait_value"]]
            .drop_duplicates("species")
            .dropna(subset=["trait_value"])
        )
        sub["tip"] = sub["species"].map(lambda sp: _first_tip_match(sp, tip_set, sp2tax, taxid2names))
        sub = sub.dropna(subset=["tip"])
        y_map = {row["tip"]: float(row["trait_value"]) for _, row in sub.iterrows()}

        if select_cfg["mode"] == "fixed":
            fixed_model = select_cfg["fixed_model"]
            best_for_trait[trait] = {
                "best": dict(
                    model=fixed_model,
                    theta=select_cfg["fixed_theta"],
                    ll=np.nan,
                    aicc=np.inf,
                    converged=True,
                    k=1 if fixed_model == "bm" else 2,
                ),
                "fits": [],
            }
        else:
            best_for_trait[trait] = _trait_model_select(
                trait, y_map, tree, tips_all, candidates=select_cfg["candidates"], bounds=select_cfg["bounds"]
            )

        fits = best_for_trait[trait]["fits"]
        if fits:
            bm = next((f for f in fits if f["model"] == "bm"), None)
            for fit in fits:
                stat, df, p = (
                    (np.nan, 0, np.nan)
                    if bm is None
                    else _lrt_bm_vs_alt(bm, fit) if fit["model"] != "bm" else (np.nan, 0, np.nan)
                )
                trait_diag_rows.append({
                    "Gene": "__TRAIT__",
                    "Position": -1,
                    "trait": trait,
                    "mapped_tips": np.nan,
                    "n_top_aa": np.nan,
                    "n_bottom_aa": np.nan,
                    "trait_model": fit["model"],
                    "trait_theta": fit.get("theta", np.nan),
                    "diag_type": "trait_model",
                    "model": fit["model"],
                    "k": int(fit.get("k", np.nan)),
                    "theta": fit.get("theta", np.nan),
                    "ll": fit.get("ll", np.nan),
                    "aicc": fit.get("aicc", np.nan),
                    "lrt_df": df,
                    "lrt_stat": stat,
                    "lrt_p": p,
                    "is_best_aicc": bool(fit is min(fits, key=lambda d: d["aicc"])),
                })

    cache = _precompute_tree_caches(tree)
    groups = list(valid_df.groupby(["Gene", "Position", "trait", "CAAP_Group"], dropna=False))
    total = len(groups)
    logging.info(f"[PGLS] evaluating {total} sites across {len(traits)} trait(s) with threads={threads}")

    def _eval_one(key, df_site):
        gene, pos, trait, caap_group = key
        change_side = "none"
        if "change_side" in df_site.columns:
            change_side = _coalesce_change_side(df_site["change_side"])
        is_caas_position = bool(df_site.get("is_caas_position", pd.Series([False])).astype(bool).any())
        base = {
            "Gene": gene,
            "Position": int(pos),
            "trait": trait,
            "CAAP_Group": caap_group,
            "change_side": change_side,
            "is_caas_position": is_caas_position,
        }
        diag_row = {"Gene": gene, "Position": int(pos), "trait": trait}
        sub = df_site[df_site["observed_dir"].isin(["top", "bottom"])].dropna(subset=["trait_value"]).copy()
        sub["tip"] = sub["species"].map(lambda sp: _first_tip_match(sp, tip_set, sp2tax, taxid2names))
        sub = sub.dropna(subset=["tip"]).drop_duplicates(subset=["tip"])
        n_top = int((sub["observed_dir"] == "top").sum())
        n_bottom = int((sub["observed_dir"] == "bottom").sum())
        mapped_n = int(len(sub))
        best = best_for_trait[trait]["best"]
        diag_row.update({
            "mapped_tips": mapped_n,
            "n_top_aa": n_top,
            "n_bottom_aa": n_bottom,
            "trait_model": best.get("model", "bm"),
            "trait_theta": best.get("theta", None),
        })

        if mapped_n == 0 or n_top < min_per_class or n_bottom < min_per_class:
            row = {
                **base,
                "n_used": 0,
                "n_top_aa": n_top,
                "n_bottom_aa": n_bottom,
                "beta_top": np.nan,
                "lambda_full": np.nan,
                "lambda_null": np.nan,
                "lrt_stat": np.nan,
                "p_pgls_two_sided": np.nan,
                "p_pgls_pos": np.nan,
                "p_pgls_neg": np.nan,
            }
            return row, diag_row

        tips_sub = list(sub["tip"])
        y = sub["trait_value"].astype(float).to_numpy()
        x = (sub["observed_dir"] == "top").astype(int).to_numpy()

        # Add log(n) as a nuisance covariate when a sample-count map is provided.
        # Both full and null models include it so the AA test is net of sampling effort.
        if n_map:
            log_n = np.array([n_map.get(sp, np.nan) for sp in sub["species"]], float)
            valid_n = np.isfinite(log_n)
            if valid_n.any():
                log_n = np.where(valid_n, log_n, np.nanmedian(log_n))
                X_full = np.column_stack([np.ones_like(x, float), log_n, x])
                X_null = np.column_stack([np.ones_like(x, float), log_n])
            else:
                X_full = np.column_stack([np.ones_like(x, float), x])
                X_null = np.ones((len(x), 1), float)
        else:
            X_full = np.column_stack([np.ones_like(x, float), x])
            X_null = np.ones((len(x), 1), float)

        row = {**base, "n_used": int(len(x)), "n_top_aa": n_top, "n_bottom_aa": n_bottom}
        try:
            site_fit = _site_lrt(
                y,
                X_full,
                X_null,
                tips_sub,
                cache,
                best["model"],
                best.get("theta"),
                select_cfg["bounds"],
            )
            row["beta_top"] = site_fit["beta_top"]
            row["lambda_full"] = site_fit["lambda_full"]
            row["lambda_null"] = site_fit["lambda_null"]
            row["lrt_stat"] = site_fit["lrt_stat"]
            row["p_pgls_two_sided"] = site_fit["p_two_sided"]
            row["p_pgls_pos"] = site_fit["p_one_tailed_top"]
            row["p_pgls_neg"] = site_fit["p_one_tailed_bottom"]
            diag_row.update({
                "diag_type":   "site_lrt",
                "model":       site_fit["diag_model"],
                "k":           site_fit["diag_k"],
                "theta":       site_fit["diag_theta"],
                "ll":          site_fit["diag_ll"],
                "aicc":        site_fit["diag_aicc"],
                "lrt_df":      site_fit["lrt_df"],
                "lrt_stat":    site_fit["lrt_stat"],
                "lrt_p":       site_fit["p_two_sided"],
                "is_best_aicc": True,
            })
        except Exception as e:
            logging.error(f"[PGLS] LRT failed for {gene}:{pos} trait={trait} model={best}: {e}")
            row["beta_top"] = np.nan
            row["lambda_full"] = np.nan
            row["lambda_null"] = np.nan
            row["lrt_stat"] = np.nan
            row["p_pgls_two_sided"] = np.nan
            row["p_pgls_pos"] = np.nan
            row["p_pgls_neg"] = np.nan

        return row, diag_row

    results: List[dict] = []
    diags: List[dict] = []

    with ThreadPoolExecutor(max_workers=max(1, int(threads))) as ex:
        futs = {ex.submit(_eval_one, key, group): key for key, group in groups}
        for i, fut in enumerate(as_completed(futs), start=1):
            row, drow = fut.result()
            results.append(row)
            diags.append(drow)
            if i % 200 == 0 or i == total:
                logging.info(f"[PGLS] progress {i}/{total} ({100.0 * i / total:0.1f}%)")

    diags.extend(trait_diag_rows)

    res_df = pd.DataFrame(results, columns=[
        "Gene", "Position", "trait", "CAAP_Group", "change_side", "is_caas_position",
        "n_used", "n_top_aa", "n_bottom_aa", "beta_top", "lambda_full", "lambda_null", "lrt_stat",
        "p_pgls_two_sided", "p_pgls_pos", "p_pgls_neg"
    ])
    diag_df = pd.DataFrame(diags, columns=DIAG_COLUMNS)
    res_df = bh_fdr_per_group(res_df, "p_pgls_two_sided", group_col="trait", alpha=float(fdr_alpha))
    res_df = bh_fdr_per_group(res_df, "p_pgls_pos", group_col="trait", alpha=float(fdr_alpha))
    res_df = bh_fdr_per_group(res_df, "p_pgls_neg", group_col="trait", alpha=float(fdr_alpha))
    dir_vals = res_df.apply(
        lambda r: _directional_summary(
            r.get("change_side", "none"),
            r.get("p_pgls_pos", np.nan),
            r.get("p_pgls_neg", np.nan),
            r.get("q_p_pgls_pos", np.nan),
            r.get("q_p_pgls_neg", np.nan),
            bool(r.get("sig_p_pgls_pos", False)),
            bool(r.get("sig_p_pgls_neg", False)),
        ),
        axis=1,
        result_type="expand",
    )
    dir_vals.columns = ["p_pgls_directional", "q_p_pgls_directional", "sig_p_pgls_directional"]
    res_df = pd.concat([res_df, dir_vals], axis=1)

    # Scoring-ready proxy: deciles of directional p-values, weighted by CAAP group,
    # then summed per Gene×Position×trait exactly like the biochemical score.
    res_df["pgls_decile_score"] = np.nan
    tested_mask = pd.to_numeric(res_df["n_used"], errors="coerce").fillna(0) > 0
    if tested_mask.any():
        for trait, grp in res_df.loc[tested_mask].groupby("trait"):
            idx = grp.index
            res_df.loc[idx, "pgls_decile_score"] = _decile_score_from_pvalues(grp["p_pgls_directional"])

    res_df["pgls_weighted_decile"] = (
        res_df["CAAP_Group"].map(SCHEME_WEIGHTS).fillna(0.0) *
        pd.to_numeric(res_df["pgls_decile_score"], errors="coerce")
    )
    char_df = (
        res_df.loc[tested_mask]
        .groupby(["Gene", "Position", "trait"], dropna=False)["pgls_weighted_decile"]
        .sum(min_count=1)
        .reset_index()
        .rename(columns={"pgls_weighted_decile": "pgls_char_score"})
    )
    if not char_df.empty:
        res_df = res_df.merge(char_df, on=["Gene", "Position", "trait"], how="left")
    else:
        res_df["pgls_char_score"] = np.nan

    excess_df = compute_excess_by_gene_trait(res_df, fdr_alpha=float(fdr_alpha))
    res_df = res_df.reindex(columns=RESULT_COLUMNS)
    return res_df, diag_df, excess_df
