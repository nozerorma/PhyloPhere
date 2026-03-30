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
from scipy.stats import chi2

from stats.stats_utils import bh_fdr_per_group

RESULT_COLUMNS = [
    "Gene",
    "Position",
    "trait",
    "n_used",
    "n_top_aa",
    "n_bottom_aa",
    "beta_top",
    "lambda_full",
    "lambda_null",
    "lrt_stat",
    "p_pgls_top",       # one-tailed: H1 = beta_top > 0  (amino acid enriched in high-trait species)
    "p_pgls_bottom",    # one-tailed: H1 = beta_top < 0  (amino acid enriched in low-trait species)
    "q_p_pgls_top",
    "sig_p_pgls_top",
    "q_p_pgls_bottom",
    "sig_p_pgls_bottom",
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


def _site_lrt(y: np.ndarray,
              X_full: np.ndarray,
              X_null: np.ndarray,
              tips: List[str],
              cache,
              model: str,
              theta: float | None,
              bounds: str) -> dict:
    full = _fit_site_model(y, X_full, tips, cache, model, theta, bounds)
    null = _fit_site_model(y, X_null, tips, cache, model, theta, bounds)
    df = X_full.shape[1] - X_null.shape[1]
    stat = max(0.0, 2.0 * (full["ll"] - null["ll"]))
    p_two = chi2.sf(stat, df)
    beta_top = float(full["beta"][1]) if len(full["beta"]) > 1 else np.nan
    return {
        "lrt_stat": float(stat),
        "lrt_df": int(df),
        "p_two_sided": float(p_two),
        "p_one_tailed_top":    float(_one_tailed_top(p_two, beta_top)),
        "p_one_tailed_bottom": float(_one_tailed_bottom(p_two, beta_top)),
        "beta_top": beta_top,
        "lambda_full": full.get("theta", np.nan),
        "lambda_null": null.get("theta", np.nan),
    }


def _first_tip_match(species: str, tip_set: Set[str], sp2tax: Dict[str, str] | None) -> str | None:
    s = str(species).strip()
    if s in tip_set:
        return s
    if sp2tax is not None:
        tid = str(sp2tax.get(s, s)).strip()
        if tid in tip_set:
            return tid
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
                  min_per_class: int,
                  select_cfg: dict,
                  threads: int = 1,
                  fdr_alpha: float = 0.05) -> tuple[pd.DataFrame, pd.DataFrame]:
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
        sub["tip"] = sub["species"].map(lambda sp: _first_tip_match(sp, tip_set, sp2tax))
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
    groups = list(valid_df.groupby(["Gene", "Position", "trait"]))
    total = len(groups)
    logging.info(f"[PGLS] evaluating {total} sites across {len(traits)} trait(s) with threads={threads}")

    def _eval_one(key, df_site):
        gene, pos, trait = key
        base = {"Gene": gene, "Position": int(pos), "trait": trait}
        diag_row = {"Gene": gene, "Position": int(pos), "trait": trait}
        sub = df_site[df_site["observed_dir"].isin(["top", "bottom"])].dropna(subset=["trait_value"]).copy()
        sub["tip"] = sub["species"].map(lambda sp: _first_tip_match(sp, tip_set, sp2tax))
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
                "p_pgls_top": np.nan,
                "p_pgls_bottom": np.nan,
            }
            return row, diag_row

        tips_sub = list(sub["tip"])
        y = sub["trait_value"].astype(float).to_numpy()
        x = (sub["observed_dir"] == "top").astype(int).to_numpy()
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
            row["p_pgls_top"]    = site_fit["p_one_tailed_top"]
            row["p_pgls_bottom"] = site_fit["p_one_tailed_bottom"]
        except Exception as e:
            logging.error(f"[PGLS] LRT failed for {gene}:{pos} trait={trait} model={best}: {e}")
            row["beta_top"] = np.nan
            row["lambda_full"] = np.nan
            row["lambda_null"] = np.nan
            row["lrt_stat"] = np.nan
            row["p_pgls_top"]    = np.nan
            row["p_pgls_bottom"] = np.nan

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

    res_df = pd.DataFrame(results, columns=RESULT_COLUMNS[:-4])
    diag_df = pd.DataFrame(diags, columns=DIAG_COLUMNS)
    res_df = bh_fdr_per_group(res_df, "p_pgls_top",    group_col="trait", alpha=float(fdr_alpha))
    res_df = bh_fdr_per_group(res_df, "p_pgls_bottom", group_col="trait", alpha=float(fdr_alpha))
    res_df = res_df.reindex(columns=RESULT_COLUMNS)
    return res_df, diag_df
