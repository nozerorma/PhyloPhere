"""The gene list of U is recorded by SCORING and checked by the registry (digest of the list as a set)."""
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

pytest.importorskip('pyarrow')

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'subworkflows/POSENRICH/local/src'
REGISTRY = SRC / 'posenrich_registry.py'
SCORING_SRC = ROOT / 'subworkflows/SCORING/local/src/null_position_stats.py'


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


sys.path.insert(0, str(SRC))
uni = load(SRC / 'posenrich_universe.py', 'posenrich_universe_identity_test')
scoring = load(SCORING_SRC, 'scoring_null_stats_identity_test')
registry_fixture = load(ROOT / 'tests/test_posenrich_registry.py', 'registry_fixture_identity_test')


def test_digest_does_not_depend_on_order_header_or_repeats(tmp_path):
    a = tmp_path / 'a.txt'
    b = tmp_path / 'b.txt'
    a.write_text('G2\nG1\nG3\n')
    b.write_text('Gene\nG1\nG3\nG2\nG1\n')
    n_a, digest_a = scoring.gene_list_identity(a)
    n_b, digest_b = scoring.gene_list_identity(b)
    assert (n_a, digest_a) == (n_b, digest_b) and n_a == 3
    assert digest_a == uni.gene_list_digest(['G1', 'G2', 'G3'])
    assert digest_a != uni.gene_list_digest(['G1', 'G2'])


def test_scoring_universe_records_the_gene_list(tmp_path):
    null = pd.DataFrame({'gene': ['G1', 'G1', 'G2'], 'position': [1, 2, 0], 'nd': [3, 1, 2], 'first_rank': [0, 1, 0]})
    null.to_parquet(tmp_path / 'null_positions.parquet', index=False)
    (tmp_path / 'filtered_discovery.tsv').write_text('Gene\tPosition\tside\tasr_path_score\nG1\t1\ttop\t0.5\n')
    (tmp_path / 'roster.txt').write_text('b_1\nb_2\n')
    (tmp_path / 'background.output').write_text('G1\t0,1,2\nG2\t0\n')
    (tmp_path / 'genes.txt').write_text('G1\nG2\n')
    out = tmp_path / 'universe'
    rc = scoring.main(['universe', '--null-positions', str(tmp_path / 'null_positions.parquet'),
                       '--postproc', str(tmp_path / 'filtered_discovery.tsv'), '--tested', str(tmp_path / 'background.output'),
                       '--universe', str(tmp_path / 'genes.txt'), '--roster', str(tmp_path / 'roster.txt'),
                       '--output-dir', str(out)])
    assert rc == 0
    base = pd.read_csv(out / 'null_universe_base.tsv', sep='\t')
    row = base[base.metric == 'gene_list_sha256'].iloc[0]
    assert int(row.value) == 2 and row.note == uni.gene_list_digest(['G1', 'G2'])


def run_registry(tmp_path, null_universe_text):
    tmp_path.mkdir(parents=True, exist_ok=True)
    files, maps = registry_fixture.make_fixture(tmp_path)
    Path(files['null_universe']).write_text(null_universe_text)
    out = tmp_path / 'registry'
    result = subprocess.run([sys.executable, str(REGISTRY), '--map_dir', maps, '--cleaned_background', files['background'],
                             '--position_scores', files['scores'], '--postproc', files['postproc'],
                             '--null_universe', files['null_universe'], '--output_dir', str(out)],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return json.loads((out / 'manifest.json').read_text())['universe']['gene_list_check'], result.stderr


def test_registry_reports_match_mismatch_and_absence(tmp_path):
    base = registry_fixture.NULL_UNIVERSE_TSV
    # the fixture's gene list is G1 only
    check, err = run_registry(tmp_path / 'match', base + f'gene_list_sha256\tNA\t1\t{uni.gene_list_digest(["G1"])}\n')
    assert check['status'] == 'match' and 'WARNING: the registry gene list' not in err
    check, err = run_registry(tmp_path / 'mismatch', base + f'gene_list_sha256\tNA\t2\t{uni.gene_list_digest(["G1", "G2"])}\n')
    assert check['status'] == 'mismatch' and check['scoring_n_genes'] == 2 and 'WARNING: the registry gene list' in err
    check, err = run_registry(tmp_path / 'absent', base)
    assert check['status'] == 'not_recorded' and check['scoring_digest'] is None
