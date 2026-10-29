"""Exercise generated reuse branches with nounset, without launching Nextflow."""
import os
from pathlib import Path
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('strict', [False, True])
@pytest.mark.parametrize('has_roster', [False, True])
def test_reuse_branch_resolves_canonical_null_under_nounset(tmp_path, strict, has_roster):
    from gui.project_io import load_project
    from gui.generation.render import render_single

    project = load_project(ROOT / 'gui/templates/tier2_marine.json')
    script = render_single(project, postproc_mode='filter',
                           reuse_exploratory=True, strict_pass=strict)
    branch = script[script.index('if [ "${POSTPROC_MODE:-}"'):script.index('# A Nextflow home')]
    source = tmp_path / 'marine_exploratory'
    reused = source / 'strict' if strict else source
    null = reused / 'caas_permulation'
    null.mkdir(parents=True)
    (null / 'cycle_summary.tsv').write_text('fixture\n')
    if has_roster:
        (null / 'cycle_roster.txt').write_text('b_1\n')
    env = dict(os.environ, CAAS_OUTBASE=str(tmp_path), TRAIT='marine', TAG='',
               POSTPROC_MODE='filter', RUN_FADE='false', RUN_RER='false')
    env.pop('CAAS_CYCLE_ROSTER_FILE', None)
    check = '\nprintf "%s\\n" "$CAAS_NULL_DIR" "$RESULTS_BASE" "$RUN_CAAS" "$RUN_DISAMBIGUATION"\n'
    run = subprocess.run(['bash', '-c', 'set -Eeuo pipefail\n' + branch + check],
                         env=env, capture_output=True, text=True)
    assert run.returncode == 0, run.stdout + run.stderr
    assert run.stdout.splitlines() == [str(null), str(tmp_path / ('marine_strict' if strict else 'marine_complete')),
                                       'false', 'false']
    assert ('cycle_roster.txt' in run.stderr) is (not has_roster)
    assert 'null_positions.parquet' not in run.stderr


@pytest.mark.parametrize('skipped', [False, True])
def test_runner_skip_flag_reflects_current_invocation(tmp_path, skipped):
    from gui.project_io import load_project
    from gui.generation.render import render_single

    script = render_single(load_project(ROOT / 'gui/templates/tier2_marine.json'))
    start = script.index('# The flag describes this invocation')
    end = script.index('# remove_work:', start)
    invocation = script[start:end]
    flag = tmp_path / 'low_contrasts.skip'
    flag.write_text('stale\n')
    stub = 'nextflow() { ' + ('printf current > "$RESULTS_BASE/low_contrasts.skip";' if skipped else ':;') + ' }\n'
    run = subprocess.run(['bash', '-c', 'set -Eeuo pipefail\n' + stub + invocation],
                         env=dict(os.environ, RESULTS_BASE=str(tmp_path), REPO_DIR=str(ROOT)),
                         capture_output=True, text=True)
    assert run.returncode == 0, run.stdout + run.stderr
    assert flag.exists() is skipped
    if skipped:
        assert flag.read_text() == 'current'
