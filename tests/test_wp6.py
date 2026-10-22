"""Directed WP6 contracts; TMPDIR and pytest basetemp must be on the validation disk."""
import os
from pathlib import Path
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]
R = os.environ.get('RSCRIPT', shutil.which('Rscript') or '')

@pytest.mark.parametrize('script', ['test_wp6_blocks.R', 'test_wp6_rer_batches.R', 'test_fcs_auc_contract.R', 'test_wp6_fcs_cli.R', 'test_wp6_rer_cli.R'])
def test_r_contract(script, tmp_path):
    if not R:
        pytest.skip('Rscript required')
    result = subprocess.run([R, str(ROOT / 'tests' / script)], cwd=ROOT,
                            env=dict(os.environ, TMPDIR=str(tmp_path)), capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'PASS' in result.stdout


def test_new_parameters_roundtrip():
    from gui.models.modules import EnrichmentConfig, RerConfig
    assert EnrichmentConfig().fcs_block_batch_genes == '1000'
    assert RerConfig().rer_batch_seed == '1998'


def test_new_parameters_reach_rendered_scripts():
    from gui.generation.render import render_batch, render_single
    from gui.project_io import load_project
    project = load_project(ROOT / 'gui/templates/tier2_marine.json')
    project.modules.rer.rer_batch_seed = '73'
    project.modules.enrichment.fcs_block_batch_genes = '17'
    batch, single = render_batch(project), render_single(project)
    assert 'export RER_BATCH_SEED="73"' in batch
    assert 'export FCS_BLOCK_BATCH_GENES="17"' in batch
    assert '"rer_batch_seed": "${RER_BATCH_SEED:-1998}"' in single
    assert '"fcs_block_batch_genes": "${FCS_BLOCK_BATCH_GENES:-1000}"' in single
