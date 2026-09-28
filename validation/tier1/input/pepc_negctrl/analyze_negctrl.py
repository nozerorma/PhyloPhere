#!/usr/bin/env python3
"""Calibration summary for the PEPC negative-control runs.

For each nc?? result dir under output/pepc_negctrl/results/:
  * family p-values: observed positions at their p.emp, positions detected by
    >= 1 null cycle but not by the observed data at p = 1;
  * BH calls, recomputed over the family (one test per position) so that an
    observed position the null never re-detects enters at k = 0, i.e.
    p = 1/(N + 1). scoring_compute.R leaves such positions NA when fewer than
    half of the observed positions are re-detected (its coordinate-mismatch
    guard), which a candidate set of one or two positions can trigger;
    `n_guard_na` counts them;
  * permutation-FDR (SAM-style) calls on the pooled max-over-sides score;
  * observed detection count against the per-cycle null detection counts.
Pooled over controls, P(p <= alpha) is compared with alpha.

Usage: analyze_negctrl.py <results_dir> <overlap_tsv> <out_tsv>
"""
import sys, glob, os
import numpy as np
import pandas as pd

res_dir, overlap_tsv, out_tsv = sys.argv[1:4]
KW = dict(sep="\t", keep_default_na=False, na_values=["", "NA"])
ALPHAS = [0.01, 0.05, 0.1]


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1.0)
    return out


def sam_q(obs_scores, null_scores, n_cycles):
    o = np.asarray(obs_scores, float)
    if not len(o):
        return o
    ns = np.sort(null_scores)
    order = np.argsort(-o)
    s = o[order]
    exp_null = (len(ns) - np.searchsorted(ns, s, side="left")) / n_cycles
    n_obs = np.array([(o >= x).sum() for x in s])
    f = np.minimum(1.0, exp_null / n_obs)
    q = np.minimum.accumulate(f[::-1])[::-1]
    out = np.empty(len(o))
    out[order] = q
    return out


rows, fam_all, obs_all = [], [], []
for d in sorted(glob.glob(os.path.join(res_dir, "nc*_complete"))):
    trait = os.path.basename(d).replace("_complete", "")
    ps_f = os.path.join(d, "scoring", "position_scores.tsv")
    cyc_f = os.path.join(d, "caas_permulation", "perm_pos_cycle_caas.tsv.gz")
    disc_f = os.path.join(d, "caastools", "discovery.tab")
    if not (os.path.exists(ps_f) and os.path.exists(cyc_f)):
        n_disc = (sum(1 for _ in open(disc_f)) - 1) if os.path.exists(disc_f) else None
        # Zero observed discoveries: the pipeline stops before the null, and
        # the run makes no calls by construction.
        status = "no_discovery" if n_disc == 0 else "missing"
        rows.append(dict(trait=trait, status=status, n_obs=0 if n_disc == 0 else np.nan,
                         n_bh_lt_05=0 if n_disc == 0 else np.nan,
                         n_bh_lt_10=0 if n_disc == 0 else np.nan,
                         n_sam_lt_10=0 if n_disc == 0 else np.nan))
        continue
    ps = pd.read_csv(ps_f, **KW)
    obs = ps.groupby("Position").agg(s=("CAAS_score", "max"), p=("p.emp", "first"))
    cyc = pd.read_csv(cyc_f, sep="\t")
    cyc = cyc[cyc.n_schemes > 0].copy()
    cyc["m"] = cyc.caas_sum / cyc.n_schemes
    pooled = cyc.groupby(["cycle", "Position"]).m.max().reset_index()
    n_cyc = 1000  # caas_full_perms in run_negctrl_local.sh; zero-detection cycles are absent from cyc
    null_pos = set(pooled.Position)
    guard_na = obs.p.isna() & ~obs.index.isin(null_pos)
    obs.loc[guard_na, "p"] = 1.0 / (n_cyc + 1)
    fam = pd.concat([obs.p, pd.Series(1.0, index=sorted(null_pos - set(obs.index)))])
    fam_adj = pd.Series(bh(fam.values), index=fam.index)
    obs["padj"] = fam_adj.loc[obs.index].values
    q = sam_q(obs.s.values, pooled.m.values, n_cyc)
    per_cycle = pooled.groupby("cycle").size()
    per_cycle = per_cycle.reindex([f"b_{i}" for i in range(1, n_cyc + 1)], fill_value=0)
    rows.append(dict(
        trait=trait, status="ok", n_obs=len(obs), n_family=len(fam),
        null_det_median=float(per_cycle.median()),
        p_null_det_ge_obs=float((per_cycle >= len(obs)).mean()),
        min_p_emp=float(obs.p.min()) if len(obs) else np.nan,
        n_p_le_05=int((obs.p <= 0.05).sum()),
        min_p_emp_adj=float(obs.padj.min()) if len(obs) else np.nan,
        n_guard_na=int(guard_na.sum()),
        n_bh_lt_05=int((obs.padj < 0.05).sum()),
        n_bh_lt_10=int((obs.padj < 0.1).sum()),
        min_sam_q=float(q.min()) if len(q) else np.nan,
        n_sam_lt_10=int((q < 0.1).sum()),
    ))
    fam_all.append(fam.values)
    obs_all.append(obs.p.values)

tab = pd.DataFrame(rows)
ov = pd.read_csv(overlap_tsv, sep="\t")
tab = tab.merge(ov[["trait", "phi", "fg_shared_with_obs"]], on="trait", how="left")
tab.to_csv(out_tsv, sep="\t", index=False)
pd.set_option("display.width", 200)
print(tab.to_string(index=False))

ok = tab[tab.status == "ok"]
done = tab[tab.status.isin(["ok", "no_discovery"])]
print(f"\nstatus: {tab.status.value_counts().to_dict()}")
if len(ok):
    fam_p = np.concatenate(fam_all)
    obs_p = np.concatenate(obs_all)
    print(f"controls with a null (pooled P below): {len(ok)}; with zero discoveries: {(tab.status == 'no_discovery').sum()}")
    print("P(p <= alpha), pooled over controls with a null:")
    for a in ALPHAS:
        print(f"  alpha={a:<5} family (null-universe, valid if <= alpha): {np.mean(fam_p <= a):.4f}"
              f"   observed-detected only: {np.mean(obs_p <= a):.4f}")
    print(f"runs with >= 1 BH call (p.emp_adj < 0.05): {(done.n_bh_lt_05 > 0).sum()}/{len(done)}")
    print(f"runs with >= 1 BH call (p.emp_adj < 0.1): {(done.n_bh_lt_10 > 0).sum()}/{len(done)}")
    print(f"runs with >= 1 permutation-FDR call (q < 0.1): {(done.n_sam_lt_10 > 0).sum()}/{len(done)}")
    print(f"observed detections above null median: {(ok.n_obs > ok.null_det_median).sum()}/{len(ok)}; "
          f"median P(null >= obs) = {ok.p_null_det_ge_obs.median():.3f}")
