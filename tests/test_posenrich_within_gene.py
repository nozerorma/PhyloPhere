"""Level 2 engine: dense oracle, independence from the block partition, hold-out false-positive rate, exclusions, CLI."""
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
spec = importlib.util.spec_from_file_location('posenrich_within_gene_under_test', SRC / 'posenrich_within_gene.py')
wg = importlib.util.module_from_spec(spec)
sys.modules['posenrich_within_gene_under_test'] = wg
spec.loader.exec_module(wg)


class Data:
    """A synthetic universe with a null of exchangeable cycles and an observed vector."""

    def __init__(self, seed=0, n_genes=8, n_cycles=60, n_terms=6, per_gene=(6, 14)):
        rng = np.random.default_rng(seed)
        sizes = rng.integers(per_gene[0], per_gene[1], n_genes)
        self.gene = np.concatenate([[f'G{g}'] * n for g, n in enumerate(sizes)])
        self.position = np.concatenate([np.arange(n) for n in sizes])
        self.gene_idx = np.repeat(np.arange(n_genes), sizes)
        P = len(self.gene)
        p0 = rng.beta(1.2, 6.0, P)                       # detectability heterogeneous across positions
        gene_active = rng.uniform(0.3, 0.9, n_genes)

        def draw(n):
            active = rng.random((n, n_genes)) < gene_active
            hit = (rng.random((n, P)) < p0) & active[:, self.gene_idx]
            return np.where(hit, rng.uniform(0.05, 1.0, (n, P)), 0.0)

        self.null = draw(n_cycles)
        self.observed = draw(1)[0]
        self.members = {}
        for t in range(n_terms):
            mask = rng.random(P) < rng.uniform(0.1, 0.5)
            self.members[f'T{t}'] = mask
        self.members['ALLG0'] = self.gene_idx == 0       # covers every position of gene 0
        self.cycles = [f'b_{i + 1}' for i in range(n_cycles)]


def write_inputs(tmp_path, d, partitions):
    """Write universe, observed, members, collections and the blocks of the given cycle partition."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({'gene': d.gene, 'position': d.position, 'nd': 1, 'observed_positive': 0}).sample(frac=1, random_state=3) \
        .to_csv(tmp_path / 'universe.tsv.gz', sep='\t', index=False)
    pos = d.observed > 0
    pd.DataFrame({'Gene': d.gene[pos], 'Position': d.position[pos], 'CAAS_score': d.observed[pos], 'side': 'top'}) \
        .to_csv(tmp_path / 'observed.tsv', sep='\t', index=False)
    registry = tmp_path / 'registry'
    registry.mkdir()
    rows = []
    for term, mask in d.members.items():
        for g, p in zip(d.gene[mask], d.position[mask]):
            rows.append(('colA' if term != 'ALLG0' else 'colB', term, term, g, p, 'r', 'direct_caas', ''))
    pd.DataFrame(rows, columns=['collection_id', 'term_id', 'term_label', 'Gene', 'Position', 'source_record_id',
                                'mapping_status', 'detail']).to_csv(registry / 'sites.tsv', sep='\t', index=False)
    pd.DataFrame(columns=['collection_id', 'term_id', 'term_label', 'Gene', 'start', 'end', 'source_record_id',
                          'mapping_status', 'detail']).to_csv(registry / 'regions.tsv', sep='\t', index=False)
    pd.DataFrame({'collection_id': ['colA', 'colB'], 'informative_within_protein': ['1', '1']}) \
        .to_csv(registry / 'collections.tsv', sep='\t', index=False)
    blocks = []
    for k, cyc in enumerate(partitions):
        block = tmp_path / f'blk_{k:04d}'
        block.mkdir()
        c_idx, p_idx = np.nonzero(d.null[cyc, :])
        frame = pd.DataFrame({'gene': d.gene[p_idx], 'cycle': [d.cycles[cyc[i]] for i in c_idx],
                              'position': d.position[p_idx].astype(np.int32), 'all': d.null[cyc][c_idx, p_idx]})
        # rows outside U (another gene, an unknown position) must be ignored
        extra = pd.DataFrame({'gene': ['X', d.gene[0]], 'cycle': [d.cycles[cyc[0]]] * 2,
                              'position': np.array([0, 999], dtype=np.int32), 'all': [0.9, 0.9]})
        frame = pd.concat([frame, extra], ignore_index=True)
        frame.to_parquet(block / 'positions.parquet', index=False)
        (block / 'cycles.txt').write_text('\n'.join(d.cycles[i] for i in cyc) + '\n')
        blocks.append(str(block))
    return blocks


def partitions(n_cycles, n_blocks):
    return [list(a) for a in np.array_split(np.arange(n_cycles), n_blocks)]


def engine(tmp_path, d, n_blocks=1, min_nonzero=3, **kwargs):
    blocks = write_inputs(tmp_path, d, partitions(len(d.cycles), n_blocks))
    U = wg.load_universe(tmp_path / 'universe.tsv.gz')
    members = wg.load_members(tmp_path / 'registry', ['colA', 'colB'], U)
    design = wg.build_design(members, U)
    w = wg.load_observed_weights(tmp_path / 'observed.tsv', U)
    out = wg.run(U, design, blocks, w, min_nonzero=min_nonzero, **kwargs)
    return U, design, w, out


def oracle(d, min_nonzero):
    """Dense statistics of every term, computed gene by gene with no sparse algebra."""
    res = {}
    for term, mask in d.members.items():
        pairs = []
        for g in np.unique(d.gene_idx[mask]):
            idx = d.gene_idx == g
            n, m = idx.sum(), (mask & idx).sum()
            if m == n:
                continue

            def r(W, idx=idx, n=n, m=m):      # bound now: the loop variables change before r is called
                return W[..., mask & idx].sum(-1) - W[..., idx].sum(-1) * m / n
            null_r = r(d.null)
            sd = null_r.std(ddof=1)
            if (np.abs(null_r) > 1e-12).sum() < min_nonzero or sd <= 0:
                continue
            pairs.append((r, sd, null_r))
        if not pairs:
            res[term] = None
            continue
        out = {}
        n_cyc = len(d.null)
        for name, scale in (('normalized', True), ('simple', False)):
            t_obs = sum(r(d.observed) / (sd if scale else 1.0) for r, sd, _ in pairs)
            t_null = np.zeros(n_cyc)
            for r, sd, null_r in pairs:
                if not scale:
                    t_null += null_r
                    continue
                s1, s2 = null_r.sum(), (null_r ** 2).sum()
                var = ((s2 - null_r ** 2) - (s1 - null_r) ** 2 / (n_cyc - 1)) / (n_cyc - 2)   # without the cycle itself
                t_null += np.where(var > 0, null_r / np.sqrt(np.where(var > 0, var, 1.0)), 0.0)
            ge = int((t_null >= t_obs - wg.TIE_TOL * (1 + abs(t_obs))).sum())
            out[name] = (float(t_obs), (1 + ge) / (n_cyc + 1), ge)
        res[term] = out
    return res


@pytest.mark.parametrize('seed', [0, 1, 2])
def test_engine_matches_dense_oracle(tmp_path, seed):
    d = Data(seed=seed)
    U, design, w, out = engine(tmp_path, d)
    results, by_gene = wg.results_tables(U, design, out, w)
    expected = oracle(d, 3)
    assert any(v is not None for v in expected.values())
    for term, exp in expected.items():
        row = results[results.term_id == term].iloc[0]
        if exp is None:
            assert row.status == 'not_identifiable'
            continue
        assert row.stat_normalized == pytest.approx(exp['normalized'][0], rel=1e-9, abs=1e-9)
        assert row.stat_simple == pytest.approx(exp['simple'][0], rel=1e-9, abs=1e-9)
        assert row.p_normalized == pytest.approx(exp['normalized'][1])
        assert row.p_simple == pytest.approx(exp['simple'][1])


def test_result_does_not_depend_on_the_block_partition(tmp_path):
    d = Data(seed=5, n_cycles=42)
    ref = None
    for n_blocks in (1, 7, 42):
        U, design, w, out = engine(tmp_path / f'p{n_blocks}', d, n_blocks=n_blocks)
        assert out['n_null'] == 42
        if ref is None:
            ref = out
            continue
        for stat in wg.STATS:
            np.testing.assert_array_equal(out['ge'][stat], ref['ge'][stat])
            np.testing.assert_allclose(out['T'][stat], ref['T'][stat], rtol=1e-9, atol=1e-9)
        np.testing.assert_array_equal(out['eligible'], ref['eligible'])
        np.testing.assert_allclose(out['fit'].sd, ref['fit'].sd, rtol=1e-9, atol=1e-12)


def test_exclusions_have_reasons(tmp_path):
    d = Data(seed=3)
    U, design, w, out = engine(tmp_path, d, min_nonzero=10_000)       # no pair reaches the stability rule
    results, by_gene = wg.results_tables(U, design, out, w)
    assert (results.status == 'not_identifiable').all() and results.p_normalized.isna().all()
    allg0 = by_gene[(by_gene.term_id == 'ALLG0') & (by_gene.gene == 'G0')].iloc[0]
    assert allg0.excluded_reason == 'no_nonmembers' and allg0.eligible == 0
    assert set(by_gene.excluded_reason) <= {'no_nonmembers', 'unstable_sd'}
    U, design, w, out = engine(tmp_path / 'ok', d, min_nonzero=3)
    results, by_gene = wg.results_tables(U, design, out, w)
    assert by_gene[(by_gene.term_id == 'ALLG0') & (by_gene.gene == 'G0')].excluded_reason.iloc[0] == 'no_nonmembers'
    assert (by_gene.loc[by_gene.eligible == 1, 'excluded_reason'] == '').all()


def test_effect_sizes_use_only_eligible_genes(tmp_path):
    d = Data(seed=7)
    U, design, w, out = engine(tmp_path, d)
    results, by_gene = wg.results_tables(U, design, out, w)
    for term in ('T0', 'T1'):
        sub = by_gene[(by_gene.term_id == term) & (by_gene.eligible == 1)]
        row = results[results.term_id == term].iloc[0]
        if sub.empty:
            continue
        freq_in = sub.positives_members.sum() / sub.n_members.sum()
        out_n = (sub.n_positions_gene - sub.n_members).sum()
        freq_out = (sub.positives_gene - sub.positives_members).sum() / out_n
        assert row.frequency_members == pytest.approx(freq_in) and row.frequency_rest == pytest.approx(freq_out)


def test_holdout_false_positive_rate(tmp_path):
    d = Data(seed=11, n_genes=30, n_cycles=420, n_terms=40, per_gene=(8, 20))
    U, design, w, out = engine(tmp_path, d, min_nonzero=10, holdout=120, holdout_seed=4)
    assert out['n_null'] == 300 and len(out['held']) == 120
    table = wg.holdout_table(design, out)
    for stat in wg.STATS:
        sub = table[(table.statistic == stat) & (table.collection_id == 'colA')]
        fpr = {a: float(sub[sub.alpha == a].fpr.iloc[0]) for a in (0.01, 0.05)}
        assert sub.n_pairs.iloc[0] >= 1000
        assert 0.025 <= fpr[0.05] <= 0.085, (stat, fpr)
        assert fpr[0.01] <= 0.03, (stat, fpr)


def test_duplicate_null_rows_are_rejected(tmp_path):
    d = Data(seed=2, n_cycles=6)
    blocks = write_inputs(tmp_path, d, partitions(6, 1))
    frame = pd.read_parquet(Path(blocks[0]) / 'positions.parquet')
    pd.concat([frame, frame.iloc[[0]]]).to_parquet(Path(blocks[0]) / 'positions.parquet', index=False)
    U = wg.load_universe(tmp_path / 'universe.tsv.gz')
    with pytest.raises(ValueError, match='more than one row'):
        wg.read_block(blocks[0], U)


def test_command_line_writes_tables_and_manifest(tmp_path):
    d = Data(seed=1)
    blocks = write_inputs(tmp_path, d, partitions(len(d.cycles), 3))
    out = tmp_path / 'out'
    wg.main(['--universe', str(tmp_path / 'universe.tsv.gz'), '--observed', str(tmp_path / 'observed.tsv'),
             '--registry-dir', str(tmp_path / 'registry'),
             '--blocks', *blocks, '--output-dir', str(out), '--min-nonzero-cycles', '3'])
    results = pd.read_csv(out / 'within_gene_results.tsv', sep='\t')
    assert {'p_normalized', 'q_normalized', 'frequency_difference', 'intensity_difference'} <= set(results.columns)
    assert len(pd.read_csv(out / 'within_gene_by_gene.tsv', sep='\t')) > 0
    manifest = json.loads((out / 'within_gene_manifest.json').read_text())
    assert manifest['counts']['n_null_cycles'] == len(d.cycles) and manifest['counts']['n_blocks'] == 3
    wg.main(['--universe', str(tmp_path / 'universe.tsv.gz'), '--observed', str(tmp_path / 'observed.tsv'),
             '--registry-dir', str(tmp_path / 'registry'), '--collections', 'colA', '--blocks', *blocks,
             '--output-dir', str(tmp_path / 'hold'), '--min-nonzero-cycles', '3', '--holdout', '10'])
    assert (tmp_path / 'hold' / 'holdout_fpr.tsv').is_file()
