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
        frame['top'] = frame['all']
        frame['bottom'] = np.nan
        frame.to_parquet(block / 'positions.parquet', index=False)
        (block / 'cycles.txt').write_text('\n'.join(d.cycles[i] for i in cyc) + '\n')
        blocks.append(str(block))
    return blocks


def partitions(n_cycles, n_blocks):
    return [list(a) for a in np.array_split(np.arange(n_cycles), n_blocks)]


def engine(tmp_path, d, n_blocks=1, min_genes=2, **kwargs):
    blocks = write_inputs(tmp_path, d, partitions(len(d.cycles), n_blocks))
    U = wg.load_universe(tmp_path / 'universe.tsv.gz')
    members = wg.load_members(tmp_path / 'registry', ['colA', 'colB'], U)
    design = wg.build_design(members, U)
    score = wg.load_observed_scores(tmp_path / 'observed.tsv', U)
    out = wg.run(U, design, blocks, score, min_genes=min_genes, **kwargs)
    return U, design, score, out


def oracle(d, min_genes):
    """Dense statistic of every term, computed gene by gene with no sparse algebra: positives of the term minus the
    positives of its genes spread evenly."""
    res = {}
    for term, mask in d.members.items():
        genes = np.unique(d.gene_idx[mask])
        if len(genes) < min_genes:
            res[term] = None
            continue
        pairs = []
        for g in genes:
            idx = d.gene_idx == g
            n, m = idx.sum(), (mask & idx).sum()
            if m == n:
                continue

            def r(W, idx=idx, mask=mask, n=n, m=m):      # bound now: the loop variables change before r is called
                B = W > 0
                return B[..., mask & idx].sum(-1) - B[..., idx].sum(-1) * m / n
            pairs.append(r)
        if not pairs:
            res[term] = None
            continue
        t_obs = sum(r(d.observed) for r in pairs)
        t_null = sum(r(d.null) for r in pairs)
        ge = int((t_null >= t_obs - wg.TIE_TOL * (1 + abs(t_obs))).sum())
        res[term] = (float(t_obs), (1 + ge) / (len(d.null) + 1), ge)
    return res


@pytest.mark.parametrize('seed', [0, 1, 2])
def test_engine_matches_dense_oracle(tmp_path, seed):
    d = Data(seed=seed)
    U, design, score, out = engine(tmp_path, d)
    results, by_gene = wg.results_tables(U, design, out, score)
    expected = oracle(d, 2)
    assert any(v is not None for v in expected.values())
    for term, exp in expected.items():
        row = results[results.term_id == term].iloc[0]
        if exp is None:
            assert row.status != 'tested' and np.isnan(row.p)
            continue
        assert row.status == 'tested'
        assert row.statistic == pytest.approx(exp[0], rel=1e-9, abs=1e-9)
        assert row.p == pytest.approx(exp[1])
        # the statistic is the positives of the term minus the positives expected from the eligible genes
        assert row.positives_members_eligible - row.expected_positives_members == pytest.approx(row.statistic)


def test_result_does_not_depend_on_the_block_partition(tmp_path):
    d = Data(seed=5, n_cycles=42)
    ref = None
    for n_blocks in (1, 7, 42):
        U, design, score, out = engine(tmp_path / f'p{n_blocks}', d, n_blocks=n_blocks)
        assert out['n_null'] == 42
        if ref is None:
            ref = out
            continue
        np.testing.assert_array_equal(out['ge'], ref['ge'])
        np.testing.assert_allclose(out['T'], ref['T'], rtol=1e-9, atol=1e-9)
        np.testing.assert_array_equal(out['eligible'], ref['eligible'])


def test_exclusions_are_fixed_by_the_annotation(tmp_path):
    d = Data(seed=3)
    U, design, score, out = engine(tmp_path, d)
    results, by_gene = wg.results_tables(U, design, out, score)
    # ALLG0 covers every position of one gene: below the minimum of genes, no p, whatever the observed positives
    row = results[results.term_id == 'ALLG0'].iloc[0]
    assert row.status == 'not_identifiable' and np.isnan(row.p) and row.n_genes_with_members_in_U == 1
    allg0 = by_gene[(by_gene.term_id == 'ALLG0') & (by_gene.gene == 'G0')].iloc[0]
    assert allg0.excluded_reason == 'no_nonmembers' and allg0.eligible == 0
    assert set(by_gene.excluded_reason) <= {'no_nonmembers', 'term_below_min_genes', ''}
    assert (by_gene.loc[by_gene.eligible == 1, 'excluded_reason'] == '').all()
    # the same annotation with the observed positives removed excludes the same pairs
    U, design, score0, out0 = engine(tmp_path / 'none', d)
    _, by_gene0 = wg.results_tables(U, design, out0, np.zeros_like(score))
    np.testing.assert_array_equal(by_gene.eligible.to_numpy(), by_gene0.eligible.to_numpy())
    # min_genes = 1 tests the one-gene term only when it has non-members
    U, design, score1, out1 = engine(tmp_path / 'one', d, min_genes=1)
    results1, _ = wg.results_tables(U, design, out1, score1)
    assert results1[results1.term_id == 'ALLG0'].iloc[0].status == 'not_identifiable'


def test_effect_sizes_use_only_eligible_genes(tmp_path):
    d = Data(seed=7)
    U, design, score, out = engine(tmp_path, d)
    results, by_gene = wg.results_tables(U, design, out, score)
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
    U, design, score, out = engine(tmp_path, d, holdout=120, holdout_seed=4)
    assert out['n_null'] == 300 and len(out['held']) == 120
    table = wg.holdout_table(design, out)
    sub = table[table.collection_id == 'colA']
    fpr = {a: float(sub[sub.alpha == a].fpr.iloc[0]) for a in (0.01, 0.05)}
    assert sub.n_pairs.iloc[0] >= 1000
    assert 0.025 <= fpr[0.05] <= 0.085, fpr
    assert fpr[0.01] <= 0.03, fpr


def test_duplicate_null_rows_are_rejected(tmp_path):
    d = Data(seed=2, n_cycles=6)
    blocks = write_inputs(tmp_path, d, partitions(6, 1))
    frame = pd.read_parquet(Path(blocks[0]) / 'positions.parquet')
    pd.concat([frame, frame.iloc[[0]]]).to_parquet(Path(blocks[0]) / 'positions.parquet', index=False)
    U = wg.load_universe(tmp_path / 'universe.tsv.gz')
    with pytest.raises(ValueError, match='more than one row'):
        wg.read_block(blocks[0], U)


def test_cli_requires_verified_universe(tmp_path):
    d = Data(seed=1)
    blocks = write_inputs(tmp_path, d, partitions(len(d.cycles),3))
    with pytest.raises((ValueError, FileNotFoundError)):
        wg.main(['--universe', str(tmp_path/'universe.tsv.gz'), '--observed',str(tmp_path/'observed.tsv'),
                 '--registry-dir',str(tmp_path/'registry'), '--blocks',*blocks,'--output-dir',str(tmp_path/'out')])


def test_source_coverage_changes_centering_and_eligible_genes():
    U=wg.Universe(['G','H'],wg.pair_keys([0,0,0,1,1],[0,1,2,0,1]),np.array([3,2]))
    members=pd.DataFrame({'collection_id':['C','C'],'term_id':['T','T'],'uidx':[0,3]})
    # H has no nonmembers in source coverage; it must not satisfy min_genes=2.
    design=wg.build_design(members,U,coverage={'C':np.array([0,1,3])})
    term_ok,pairs=wg.eligibility(design,2)
    assert not term_ok.any() and not pairs.any()
    term_ok,pairs=wg.eligibility(design,1)
    assert term_ok.all() and pairs.tolist()==[True,False]
    positive=np.array([1.,0.,1.,1.,0.])
    import scipy.sparse as sp
    T=wg.term_statistic(sp.csr_matrix(positive.reshape(1,-1)),design,np.flatnonzero(pairs),sp.csr_matrix([[1.]]),20000)
    assert T[0,0]==pytest.approx(.5)  # G: 1 - 1*(1/2); its uncovered detection at 2 contributes nothing.


def test_missing_null_is_not_a_negative_result():
    U=wg.Universe(['G'],wg.pair_keys([0,0],[0,1]),np.array([2]))
    members=pd.DataFrame({'collection_id':['C'],'term_id':['T'],'uidx':[0]})
    design=wg.build_design(members,U)
    with pytest.raises(ValueError,match='Null blocks'):
        wg.run(U,design,[],np.array([1.,0.]),min_genes=1)
