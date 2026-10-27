"""Universe U: definition, gene list as input, tested-set rule, manifest, and agreement with SCORING and Accumulation."""
import importlib.util
import json
import math
import subprocess
import sys
from argparse import Namespace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import hypergeom

pytest.importorskip('pyarrow')

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'subworkflows/ENRICHMENT/local/src'
UNIVERSE = SRC / 'posenrich_universe.py'
SCORING_SRC = ROOT / 'subworkflows/SCORING/local/src/null_position_stats.py'
ACCUMULATION = ROOT / 'subworkflows/CT_ACCUMULATION/local/accumulation_test.py'


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


uni = load(UNIVERSE, 'posenrich_universe_under_test')

# G1 tests 0..9, G2 0..5, G3 0..4. The null scores G1 {1,2,3,12,15}, G2 {0,1}, G3 {0}; G3 is outside the gene list.
BACKGROUND = 'G1\t' + ','.join(map(str, range(10))) + '\nG2\t' + ','.join(map(str, range(6))) + \
             '\nG3\t' + ','.join(map(str, range(5))) + '\n'
NULL = pd.DataFrame({'gene': ['G1'] * 5 + ['G2'] * 2 + ['G3'],
                     'position': [1, 2, 3, 12, 15, 0, 1, 0],
                     'nd': [5, 2, 1, 9, 2, 3, 4, 7],
                     'first_rank': [0] * 8})
# Observed rows: G1:1 on both sides, G1:5 with a non-positive score, G1:12 (null-scored, untested), G1:20 (neither),
# G2:0, and G3:0 (gene outside the list).
OBSERVED = ('Gene\tPosition\tside\tasr_path_score\n'
            'G1\t1\ttop\t0.9\nG1\t1\tbottom\t0.8\nG1\t5\ttop\t0.0\nG1\t12\ttop\t0.4\nG1\t20\tbottom\t0.5\n'
            'G2\t0\ttop\t0.7\nG3\t0\ttop\t0.6\n')


@pytest.fixture
def files(tmp_path):
    (tmp_path / 'genes.txt').write_text('G1\nG2\n')
    (tmp_path / 'background.output').write_text(BACKGROUND)
    NULL.to_parquet(tmp_path / 'null_positions.parquet', index=False)
    (tmp_path / 'observed.tsv').write_text(OBSERVED)
    return tmp_path


def build(files, genes=None, background='background.output'):
    genes = genes or uni.load_universe_genes(files / 'genes.txt')
    return uni.build_universe(genes, uni.load_null_positive(files / 'null_positions.parquet'),
                              uni.load_observed_positive(files / 'observed.tsv'),
                              files / background if background else None)


def as_pairs(universe):
    table = universe.table()
    return {(g, int(p)): int(nd) for g, p, nd in zip(table.gene, table.position, table.nd)}


def test_universe_is_tested_and_null_positive_plus_observed_positives(files):
    u = build(files)
    # G1:15 is null-positive but untested; G1:12 and G1:20 are observed positives and stay; G3 is outside the list.
    assert as_pairs(u) == {('G1', 1): 5, ('G1', 2): 2, ('G1', 3): 1, ('G1', 12): 9, ('G1', 20): 0,
                           ('G2', 0): 3, ('G2', 1): 4}
    assert u.metrics['tested_set'] == 'used'
    assert u.metrics['n_universe'] == 7
    assert u.metrics['n_observed_positive'] == 4          # G1:1, G1:12, G1:20, G2:0 (G1:5 is not positive; G3 is dropped)
    assert u.metrics['n_observed_outside_tested'] == 2    # G1:12 and G1:20
    assert u.metrics['n_observed_not_scored_by_null'] == 1  # G1:20


def test_observed_sides_and_both_count_once(files):
    u = build(files)
    names = lambda keys: {(u.genes[k >> uni.POS_BITS], int(k & 0xFFFFFFFF)) for k in keys}
    assert names(u.top) == {('G1', 1), ('G1', 12), ('G2', 0)}
    assert names(u.bottom) == {('G1', 1), ('G1', 20)}
    assert u.observed.size == 4


def test_gene_list_is_an_input(files):
    u = build(files, genes=['G1', 'G2', 'G3'])
    assert ('G3', 0) in as_pairs(u) and u.metrics['n_universe'] == 8


def test_incompatible_tested_set_is_ignored(files):
    (files / 'far.output').write_text('G1\t100,101\nG2\t100\nG3\t100\n')
    u = build(files, background='far.output')
    assert u.metrics['tested_set'] == 'ignored_incompatible'
    assert ('G1', 15) in as_pairs(u)                       # every null-positive position of the list


def test_without_background_every_null_positive_position_of_the_list(files):
    u = build(files, background=None)
    assert u.metrics['tested_set'] == 'not_given' and ('G1', 15) in as_pairs(u)


def test_positive_rule_is_the_core_rule():
    core = load(ROOT / 'subworkflows/CT_DISAMBIGUATION/local/src/core/scores.py', 'core_scores_under_test')
    for value in (0, -1.5, 0.1, 3, float('nan'), float('inf'), None, '0.5', 'x'):
        try:
            expected = core.position_score(value) is not None
        except (TypeError, ValueError):
            expected = False
        assert uni.is_positive_score(value) == (expected and math.isfinite(float(value)) if expected else False)


def test_position_scores_table_is_accepted(tmp_path):
    (tmp_path / 'p.tsv').write_text('Gene\tPosition\tCAAS_score\tside\nG1\t1\t0.5\ttop\nG1\t2\t0.2\tbottom\n')
    df = uni.load_observed_positive(tmp_path / 'p.tsv')
    assert len(df) == 2


def test_command_writes_table_and_manifest(files, tmp_path):
    out = tmp_path / 'out'
    subprocess.run([sys.executable, str(UNIVERSE), '--background-genes', str(files / 'genes.txt'),
                    '--background-positions', str(files / 'background.output'),
                    '--null-positions', str(files / 'null_positions.parquet'),
                    '--observed', str(files / 'observed.tsv'), '--output-dir', str(out)], check=True)
    table = pd.read_csv(out / 'universe.tsv.gz', sep='\t')
    assert len(table) == 7 and table['observed_positive'].sum() == 4
    manifest = json.loads((out / 'universe_manifest.json').read_text())
    assert manifest['gene_list_sha256'] == uni.sha256_file(files / 'genes.txt')
    assert manifest['metrics']['n_universe'] == 7


def test_agrees_with_scoring_universe(files):
    scoring = load(SCORING_SRC, 'scoring_null_stats_under_test')
    import pyarrow as pa
    obs = uni.load_observed_positive(files / 'observed.tsv')
    obs = obs[obs.gene.isin(['G1', 'G2'])].drop_duplicates(['gene', 'position'])  # SCORING keeps observed rows of genes outside the list
    classes, _, _, metrics = scoring.build_universe(
        pa.Table.from_pandas(NULL, preserve_index=False), list(obs.gene), obs.position.to_numpy(),
        files / 'background.output', files / 'genes.txt', ['b_1', 'b_2'])
    scoring_keys = set(zip(classes['gene'].to_pylist(), classes['position'].to_pylist()))
    u = build(files)
    in_null = set(zip(NULL.gene, NULL.position))
    assert {k for k in as_pairs(u) if k in in_null} == scoring_keys
    assert metrics['n_null_u'] == len(scoring_keys)


def test_accumulation_uses_the_shared_universe(files, tmp_path):
    acc = load(ACCUMULATION, 'accumulation_under_test')
    for side, K, actual_g1 in (('both', 4, 3), ('top', 3, 2), ('bottom', 2, 2)):
        out = tmp_path / f'{side}.csv'
        acc.run(Namespace(caas_csv=str(files / 'observed.tsv'), background_positions=str(files / 'background.output'),
                          bg_caas=str(files / 'genes.txt'), null_positions=str(files / 'null_positions.parquet'),
                          change_side=side, output=str(out)))
        res = pd.read_csv(out).set_index('Gene')
        assert res.loc['G1', 'PoolPositions'] == 5 and res.loc['G2', 'PoolPositions'] == 2
        assert res.loc['G1', 'ActualCount'] == actual_g1
        expected = hypergeom.sf(actual_g1 - 1, 7, 5, K)
        assert res.loc['G1', 'PValue'] == pytest.approx(expected)
        assert res['ActualCount'].sum() == K
