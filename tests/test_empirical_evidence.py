"""Scientific and project contracts for replacing score percentiles by empirical evidence."""
import copy
import os
from pathlib import Path
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_shared_evidence_oracles(tmp_path):
    r = shutil.which('Rscript')
    if not r:
        pytest.skip('Rscript required')
    result = subprocess.run([r, 'tests/test_scoring_evidence.R', str(ROOT), str(tmp_path)],
                            cwd=ROOT, env=dict(os.environ, TMPDIR=str(tmp_path)),
                            capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'strict boundary' in result.stdout and 'batch invariance' in result.stdout


def test_percentile_project_migration_never_reinterprets_fraction(caplog):
    from gui.models.serialization import from_dict, to_dict
    data = {'modules': {'scoring': {'scoring_gene_top_pct': '0.25', 'scoring_position_top_pct': '0.10'}}}
    project = from_dict(copy.deepcopy(data))
    scoring = to_dict(project)['modules']['scoring']
    assert scoring['evidence_p_threshold'] == '0.05'
    assert scoring['evidence_fdr_threshold'] == '0.05'
    assert not any(k in scoring for k in data['modules']['scoring'])
    assert 'not a p-value threshold' in caplog.text
    explicit = from_dict({'modules': {'scoring': {'evidence_p_threshold': '0.02'}, 'rer': {'rer_pval_threshold': '0.01', 'rer_pval_column': 'p.adj'}}})
    assert explicit.modules.scoring.evidence_p_threshold == '0.02'
    assert not hasattr(explicit.modules.rer, 'rer_pval_column')
