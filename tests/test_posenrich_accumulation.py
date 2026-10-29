"""Level 1: share of detections per gene calibrated with the null cycles, against a dense oracle."""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip('pyarrow')
pytest.importorskip('scipy')

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'subworkflows/POSENRICH/local/src'
sys.path.insert(0, str(SRC))
spec = importlib.util.spec_from_file_location('accumulation_under_test', SRC / 'accumulation_test.py')
acc = importlib.util.module_from_spec(spec)
sys.modules['accumulation_under_test'] = acc
spec.loader.exec_module(acc)


class Data:
    """A synthetic universe, a null of exchangeable cycles with a variable total of detections, and an observed vector."""

    def __init__(self, seed=0, n_genes=7, n_cycles=48, empty_cycle=None):
        rng = np.random.default_rng(seed)
        sizes = rng.integers(5, 15, n_genes)
        self.gene = np.concatenate([[f'G{g}'] * n for g, n in enumerate(sizes)])
        self.position = np.concatenate([np.arange(n) for n in sizes])
        self.gene_idx = np.repeat(np.arange(n_genes), sizes)
        P = len(self.gene)
        p0 = rng.beta(1.2, 5.0, P)
        scale = rng.uniform(0.3, 1.5, n_cycles)             # the total of a cycle varies
        self.null = rng.random((n_cycles, P)) < np.clip(p0[None, :] * scale[:, None], 0, 1)
        if empty_cycle is not None:
            self.null[empty_cycle] = False
        self.observed = rng.random(P) < np.clip(p0 * 1.4, 0, 1)
        self.top = self.observed & (rng.random(P) < 0.6)
        self.bottom = self.observed & ~self.top | (self.observed & (rng.random(P) < 0.1))
        self.cycles = [f'b_{i + 1}' for i in range(n_cycles)]


def write_inputs(tmp_path, d, n_blocks):
    tmp_path.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({'gene': d.gene, 'position': d.position, 'nd': 1, 'observed_positive': d.observed.astype(int),
                  'observed_top': d.top.astype(int), 'observed_bottom': d.bottom.astype(int)}) \
        .sample(frac=1, random_state=2).to_csv(tmp_path / 'universe.tsv.gz', sep='\t', index=False)
    blocks = []
    for k, cyc in enumerate(np.array_split(np.arange(len(d.cycles)), n_blocks)):
        block = tmp_path / f'blk_{k:04d}'
        block.mkdir()
        c_idx, p_idx = np.nonzero(d.null[cyc])
        pd.DataFrame({'gene': d.gene[p_idx], 'cycle': [d.cycles[cyc[i]] for i in c_idx],
                      'position': d.position[p_idx].astype(np.int32), 'all': 0.5}).to_parquet(block / 'positions.parquet', index=False)
        (block / 'cycles.txt').write_text('\n'.join(d.cycles[i] for i in cyc) + '\n')
        blocks.append(str(block))
    return blocks


def oracle(d):
    """Dense share of every gene in the observed data and in every cycle with a detection."""
    n_genes = d.gene_idx.max() + 1
    count = lambda m: np.array([m[..., d.gene_idx == g].sum(-1) for g in range(n_genes)]).T
    obs_c = count(d.observed)
    obs_share = obs_c / obs_c.sum()
    null_c = count(d.null)
    K = null_c.sum(axis=1)
    share = null_c[K > 0] / K[K > 0, None]
    ge = (share >= obs_share - acc.TIE_TOL * (1 + np.abs(obs_share))).sum(axis=0)
    return obs_c, (1 + ge) / (share.shape[0] + 1), share.shape[0], share.mean(axis=0)


def run(tmp_path, d, n_blocks=1):
    blocks = write_inputs(tmp_path, d, n_blocks)
    out = tmp_path / 'out'
    acc.main(['--universe', str(tmp_path / 'universe.tsv.gz'), '--blocks', *blocks, '--output-dir', str(out)])
    return pd.read_csv(out / 'accumulation_results.tsv', sep='\t').set_index('gene'), \
        json.loads((out / 'accumulation_manifest.json').read_text())


@pytest.mark.parametrize('seed,n_blocks', [(0, 1), (1, 3), (2, 5)])
def test_matches_the_dense_oracle_whatever_the_block_partition(tmp_path, seed, n_blocks):
    d = Data(seed=seed)
    res, manifest = run(tmp_path, d, n_blocks)
    obs_c, p, n, mean = oracle(d)
    tested = obs_c > 0
    np.testing.assert_array_equal(res['positives'].to_numpy(), obs_c)
    np.testing.assert_allclose(res['p'].to_numpy()[tested], p[tested])
    assert res['p'][~tested].isna().all() and res['q'][~tested].isna().all()
    np.testing.assert_allclose(res['null_share_mean'].to_numpy(), mean)
    assert manifest['counts']['n_null_cycles'] == n == res['n_null_cycles'].iloc[0]


def test_descriptive_side_counts_and_universe_sizes(tmp_path):
    d = Data(seed=3)
    res, manifest = run(tmp_path, d)
    assert res['positives_top'].sum() == d.top.sum() and res['positives_bottom'].sum() == d.bottom.sum()
    assert (res['positives_top'] <= res['positives']).all()
    assert res['n_positions'].sum() == len(d.gene) == manifest['counts']['n_universe']
    assert res['share'].sum() == pytest.approx(1.0)


def test_cycles_without_detections_are_left_out(tmp_path):
    d = Data(seed=4, empty_cycle=7)
    res, manifest = run(tmp_path, d)
    assert manifest['counts']['n_null_cycles'] == len(d.cycles) - 1
    assert manifest['counts']['n_null_cycles_without_positives'] == 1
    obs_c, p, _, _ = oracle(d)
    np.testing.assert_allclose(res['p'].to_numpy()[obs_c > 0], p[obs_c > 0])


def test_q_is_benjamini_hochberg_over_the_tested_genes(tmp_path):
    d = Data(seed=5)
    res, _ = run(tmp_path, d)
    tested = res['p'].notna()
    p = res.loc[tested, 'p'].to_numpy()
    order = np.argsort(p, kind='stable')
    ranked = np.minimum.accumulate((p[order] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    expected = np.empty(len(p))
    expected[order] = np.minimum(ranked, 1.0)
    np.testing.assert_allclose(res.loc[tested, 'q'].to_numpy(), expected)
