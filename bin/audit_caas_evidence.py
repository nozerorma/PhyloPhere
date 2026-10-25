#!/usr/bin/env python3
"""Read-only CAAS p-value audit and held-out calibration from canonical blocks.

This is validation tooling, not a production evidence selector. Undetected test
positions have p=1. All calibration classes and tails use reference blocks only.
The fixed evaluation domain is explicitly the supplied full-run positional family;
results conditional on this domain do not establish validity on unobserved sites.
Run large datasets on Slurm, never a cluster login node.
"""
from __future__ import annotations
import argparse
import csv
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np

TOL = 1e-12
ALPHAS = np.array([0.0001, 0.001, 0.005, 0.01, 0.025, 0.05])


def bh(p):
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    ok = np.flatnonzero(~np.isnan(p))
    order = ok[np.argsort(p[ok], kind="stable")]
    q = p[order] * len(order) / np.arange(1, len(order) + 1)
    out[order] = np.minimum(1, np.minimum.accumulate(q[::-1])[::-1])
    return out


def propensity_breaks(nd, classes=20):
    values, counts = np.unique(nd[nd > 0], return_counts=True)
    if not len(values):
        return np.array([], float)
    weight = values.astype(float) * counts
    cls = np.minimum(classes, np.ceil(classes * np.cumsum(weight) / weight.sum() - 1e-9))
    last = np.flatnonzero(np.r_[cls[1:] != cls[:-1], True])
    return values[last][:-1].astype(float)


def assign_class(breaks, nd):
    return np.searchsorted(breaks, np.maximum(nd, 1), side="left")


def factorized_p(score, nd, breaks, pools, cycles):
    """Historical removed formula, including observed +1 class assignment and empty pool."""
    score, nd = np.broadcast_arrays(np.asarray(score, float), np.asarray(nd, int))
    cls = assign_class(breaks, nd + 1)
    out = (nd + 1) / (cycles + 1)
    out = out.astype(float, copy=True)
    for k, pool in enumerate(pools):
        mask = cls == k
        if not np.any(mask):
            continue
        count = len(pool) - np.searchsorted(pool, score[mask] - TOL, side="left")
        out[mask] *= (count + 1) / (len(pool) + 1) if len(pool) else 1
    out = np.minimum(1, out)
    out[score <= TOL] = 1
    out[np.isnan(score)] = np.nan
    return out


class PositionReference:
    """Sorted sparse detection scores, including each position's non-detection mass."""
    def __init__(self, ids, values, n_positions, cycles):
        self.cycles = cycles
        self.nd = np.bincount(ids, minlength=n_positions)
        order = np.lexsort((values, ids))
        self.sorted = np.empty(len(ids), dtype=[("id", "i4"), ("value", "f8")])
        self.sorted["id"] = ids[order]
        self.sorted["value"] = values[order]
        self.end = np.cumsum(self.nd)
        self.breaks = propensity_breaks(self.nd)
        cls = assign_class(self.breaks, self.nd[ids])
        self.pools = [np.sort(values[cls == k]) for k in range(len(self.breaks) + 1)]
        sums = np.bincount(ids, weights=values, minlength=n_positions)
        self.mean = np.divide(sums, self.nd, out=np.zeros(n_positions), where=self.nd > 0)

    def direct(self, ids, score):
        query = np.empty(len(ids), dtype=self.sorted.dtype)
        query["id"], query["value"] = ids, score - TOL
        at = np.searchsorted(self.sorted, query, side="left")
        out = (self.end[ids] - at + 1) / (self.cycles + 1)
        out[score <= TOL] = 1
        return out

    def factorized(self, ids, score):
        return factorized_p(score, self.nd[ids], self.breaks, self.pools, self.cycles)


def load_module(path):
    spec = importlib.util.spec_from_file_location("audit_null_position_stats", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def write_csv(path, rows):
    rows = list(rows)
    if not rows:
        return
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)


def calibration(ids, values, cycles, block_ids, roster, n_positions, per_block, out):
    """Two complementary block folds; test cycles sampled before looking at scores."""
    rng = np.random.default_rng(1998)
    summaries, cycle_rows, group_rows = [], [], []
    for fold in (0, 1):
        reference_blocks = np.flatnonzero(np.arange(len(block_ids)) % 2 == fold)
        test_blocks = np.flatnonzero(np.arange(len(block_ids)) % 2 != fold)
        reference_cycles = np.concatenate([block_ids[i] for i in reference_blocks])
        test_cycles = np.concatenate([rng.choice(block_ids[i], min(per_block, len(block_ids[i])), replace=False) for i in test_blocks])
        train = np.isin(cycles, reference_cycles)
        ref = PositionReference(ids[train], values[train], n_positions, len(reference_cycles))
        del train
        prop_group = np.digitize(ref.nd, [1, 2, 6, 21, 101, 501, 2501])
        cls = assign_class(ref.breaks, ref.nd + 1)
        quartile = np.zeros(n_positions, dtype=int)
        for c in np.unique(cls):
            at = np.flatnonzero((cls == c) & (ref.nd > 0))
            if len(at):
                qs = np.quantile(ref.mean[at], [.25, .5, .75])
                quartile[at] = np.searchsorted(qs, ref.mean[at], side="right") + 1
        groups = [("all", np.zeros(n_positions, int)) , ("reference_detections", prop_group), ("class_mean_quartile", cls * 5 + quartile)]
        sizes = {name: np.bincount(g) for name, g in groups}
        tallies = {name: np.zeros((3, len(ALPHAS), len(sizes[name]), len(test_cycles)), np.int64) for name, _ in groups}
        selected = np.isin(cycles, test_cycles)
        ti, tv, tc = ids[selected], values[selected], cycles[selected]
        del selected
        order = np.argsort(tc, kind="stable")
        ti, tv, tc = ti[order], tv[order], tc[order]
        single_pool = np.sort(ref.sorted["value"])
        for j, cy in enumerate(test_cycles):
            lo, hi = np.searchsorted(tc, cy, side="left"), np.searchsorted(tc, cy, side="right")
            ix, v = ti[lo:hi], tv[lo:hi]
            direct, fact = ref.direct(ix, v), ref.factorized(ix, v)
            single = factorized_p(v, ref.nd[ix], np.array([]), [single_pool], ref.cycles)
            for method, p in enumerate((direct, fact, single)):
                label = ("p.emp", "p.emp_fact", "single_pool_control")[method]
                full = np.ones(n_positions); full[ix] = p
                q = bh(full)
                for k, alpha in enumerate(ALPHAS):
                    hit = p <= alpha
                    for name, group in groups:
                        tallies[name][method, k, :, j] = np.bincount(group[ix[hit]], minlength=len(sizes[name]))
                    cycle_rows.append(dict(fold=fold, cycle=roster[cy], method=label, alpha=alpha, rejected=int(hit.sum()), family=n_positions, rate=float(hit.sum()/n_positions), detected=len(ix), bh_005=int(np.sum(q < .05))))
        for name, group in groups:
            for g, size in enumerate(sizes[name]):
                if not size: continue
                for method, label in enumerate(("p.emp", "p.emp_fact", "single_pool_control")):
                    for k, alpha in enumerate(ALPHAS):
                        counts = tallies[name][method, k, g]
                        rates = counts / size
                        group_rows.append(dict(fold=fold, grouping=name, group=g, method=label, alpha=alpha, positions=int(size), test_cycles=len(test_cycles), rejected=int(counts.sum()), pairs=int(size*len(test_cycles)), rate=float(rates.mean()), sd_across_cycles=float(rates.std(ddof=1)), cycle_q025=float(np.quantile(rates,.025)),cycle_q975=float(np.quantile(rates,.975))))
        summaries.append(dict(fold=fold, reference_cycles=len(reference_cycles), evaluation_cycles=[roster[c] for c in test_cycles], reference_blocks=reference_blocks.tolist(), reference_detections=int(len(ref.sorted)), class_breaks=ref.breaks.tolist(), classes=len(ref.pools)))
        print("calibration fold", fold, "done", flush=True)
    write_csv(out / "calibration_cycles.tsv", cycle_rows)
    write_csv(out / "calibration_strata.tsv", group_rows)
    (out / "calibration_design.json").write_text(json.dumps(summaries, indent=2))


def main():
    import pandas as pd
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo", required=True, type=Path)
    ap.add_argument("--null-dir", required=True, type=Path)
    ap.add_argument("--postproc", required=True, type=Path)
    ap.add_argument("--tested", required=True, type=Path)
    ap.add_argument("--universe", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    ap.add_argument("--published-positions", type=Path)
    ap.add_argument("--calibration-per-block", type=int, default=10)
    a = ap.parse_args(); out = a.output_dir; out.mkdir(parents=True, exist_ok=True)
    nps = load_module(a.repo / "subworkflows/SCORING/local/src/null_position_stats.py")
    og, op, ov = nps.observed_best(a.postproc)
    roster = a.null_dir.joinpath("cycle_roster.txt").read_text().splitlines()
    assert len(roster) == len(set(roster)) and "b_0" not in roster
    counts_meta = json.loads((a.null_dir / "counts.json").read_text())
    assert counts_meta["cycles"] == len(roster)
    assert counts_meta["score_rule"] == nps.SCORE_RULE
    N = len(roster); cmap = {c:i for i,c in enumerate(roster)}
    null = pq.read_table(a.null_dir / "null_positions.parquet")
    classes, breaks, sizes, metrics = nps.build_universe(null, og, op, a.tested, a.universe, roster)
    nps.write_universe(out, classes, breaks, sizes, metrics, N, N)
    vocab = nps.GeneCodes(set(pc.unique(null["gene"]).to_pylist()) | set(og))
    ck = vocab.keys(classes["gene"], classes["position"].to_numpy())
    ok = vocab.keys(og, op)
    family_keys = np.union1d(ck, ok)
    obs_id = np.searchsorted(family_keys, ok)
    m = len(family_keys)
    chunks_i, chunks_v, chunks_c, block_ids = [], [], [], []
    for block in sorted((a.null_dir / "blocks").glob("blk_*")):
        cy = block.joinpath("cycles.txt").read_text().splitlines()
        block_ids.append(np.array([cmap[c] for c in cy], int))
        t = pq.read_table(block / "positions.parquet")
        key = vocab.keys(t["gene"], t["position"].to_numpy())
        idx = np.searchsorted(family_keys, key)
        keep = (idx < m) & (family_keys[np.minimum(idx, m-1)] == key)
        idx = idx[keep].astype(np.int32)
        val = t["all"].to_numpy()[keep]
        mapped_cycles = pc.index_in(t["cycle"], value_set=pa.array(roster)).to_numpy()[keep]
        assert np.all(np.isfinite(mapped_cycles))
        cycodes = mapped_cycles.astype(np.int16)
        assert np.all(np.isin(cycodes, block_ids[-1]))
        assert np.all(np.isfinite(val)) and np.all((val > 0) & (val <= 1+TOL))
        assert len(np.unique(idx.astype(np.int64)*N+cycodes)) == len(idx), "duplicate position/cycle"
        chunks_i.append(idx); chunks_v.append(val); chunks_c.append(cycodes)
        print(block.name, len(idx), "family detections", flush=True)
    ids, values, cycles = np.concatenate(chunks_i), np.concatenate(chunks_v), np.concatenate(chunks_c)
    del chunks_i, chunks_v, chunks_c, null
    assert len(np.concatenate(block_ids)) == N and len(np.unique(np.concatenate(block_ids))) == N
    nd = np.bincount(ids, minlength=m)
    assert np.array_equal(nd[np.searchsorted(family_keys,ck)],classes["nd"].to_numpy())
    threshold = np.full(m, np.inf); threshold[obs_id] = ov-TOL
    k = np.bincount(ids[values >= threshold[ids]], minlength=m)[obs_id]
    p = (k+1)/(N+1)
    p[ov <= TOL] = 1
    all_cls = assign_class(breaks, nd[ids])
    pools = [np.sort(values[all_cls == c]) for c in range(len(breaks)+1)]
    assert [len(p) for p in pools] == [sizes[c+1] for c in range(len(pools))]
    fact = factorized_p(ov, nd[obs_id], breaks, pools, N)
    guard = len(og)>=10 and np.mean(nd[obs_id]>0)<.5
    if guard: p[nd[obs_id]==0]=np.nan;fact[nd[obs_id]==0]=np.nan
    pe=np.ones(m);pf=np.ones(m);pe[obs_id]=p;pf[obs_id]=fact
    qe,qf=bh(pe)[obs_id],bh(pf)[obs_id]
    table=pd.DataFrame(dict(Gene=og,Position=op,score=ov,nd=nd[obs_id],k_emp=k,cls_null=assign_class(breaks,nd[obs_id])+1,cls_observed=assign_class(breaks,nd[obs_id]+1)+1,**{"p.emp":p,"p.emp_fact":fact,"p.adj_bh":qe,"p.adj_bh_fact":qf}))
    table.to_csv(out/'position_evidence.tsv',sep='\t',index=False)
    pp=pd.read_csv(a.postproc,sep='\t',keep_default_na=False,na_values=[''],float_precision='round_trip')
    side=pp[['Gene','Position','side','asr_path_score']].drop_duplicates(['Gene','Position','side'])
    side=side[side['side'].isin(['top','bottom']) & (pd.to_numeric(side['asr_path_score'],errors='coerce')>0)]
    side.merge(table,on=['Gene','Position'],validate='many_to_one').to_csv(out/'position_evidence_by_side.tsv',sep='\t',index=False)
    published_error=None
    if a.published_positions:
        published=pd.read_csv(a.published_positions,sep='\t').drop_duplicates(['Gene','Position'])
        merged=table.merge(published,on=['Gene','Position'],suffixes=('_audit','_published'),validate='one_to_one')
        assert len(merged)==len(table)==len(published)
        published_error={c:float(np.nanmax(np.abs(merged[c+'_audit']-merged[c+'_published']))) for c in ['p.emp','p.emp_fact','p.adj_bh','p.adj_bh_fact']}
        assert max(published_error.values())<1e-12,published_error
    b0=pq.read_table(a.null_dir/'b0/positions.parquet').to_pandas()
    b0['Gene']=b0['gene'];b0['Position']=b0['position']
    cmp=table.merge(b0,on=['Gene','Position'],how='outer',indicator=True)
    b0_comparison=dict(observed=len(table),b0=len(b0),shared=int((cmp['_merge']=='both').sum()),observed_only=int((cmp['_merge']=='left_only').sum()),b0_only=int((cmp['_merge']=='right_only').sum()),max_shared_score_difference=float(np.nanmax(np.abs(cmp['score']-cmp['all']))))
    summary=dict(N=N,family=m,observed=len(table),side_rows=len(side),observed_genes=len(set(og)),null_unseen=int(np.sum(nd[obs_id]==0)),guard=bool(guard),metrics=metrics,classes=len(pools),class_assignment_shifts=int(np.sum(table.cls_null!=table.cls_observed)),published_max_error=published_error,b0_comparison=b0_comparison,source_paths={k:str(v) for k,v in vars(a).items() if isinstance(v,Path)})
    (out/'summary.json').write_text(json.dumps(summary,indent=2))
    # Free full-fit pools before constructing held-out reference distributions.
    del pools, all_cls, pe, pf
    if a.calibration_per_block:calibration(ids,values,cycles,block_ids,roster,m,a.calibration_per_block,out)
    print(json.dumps(summary,indent=2),flush=True)

if __name__ == '__main__':
    main()
