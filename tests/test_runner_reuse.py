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


def test_strict_pass_precomputed_caas_tables_come_from_strict_dir(tmp_path):
    """Precomputed discovery in the strict pass must read <complete>/strict, never the complete tables."""
    from gui.project_io import load_project
    from gui.generation.render import render_single

    project = load_project(ROOT / 'gui/templates/tier2_marine.json')
    script = render_single(project, postproc_mode='filter', reuse_exploratory=True, strict_pass=True)
    block = script[script.index('PRECOMP_OUTDIR="${PRECOMP_BASE_PATH:-}/${TRAIT}"'):script.index('if [ "${PRECOMP_USE_RER:-false}" = true ]')]
    complete = tmp_path / 'complete'
    strict = complete / 'strict'
    (strict / 'meta_caas/meta_caas').mkdir(parents=True)
    (strict / 'caastools').mkdir()
    (strict / 'caas_permulation').mkdir()
    (strict / 'caas_permulation/cycle_roster.txt').write_text('b_1\n')
    (strict / 'caastools/discovery.tab').write_text('x\n')
    (strict / 'meta_caas/meta_caas/meta_caas.tsv').write_text('x\n')
    env = dict(os.environ, PRECOMP_BASE_PATH=str(tmp_path), TRAIT='nope', PASS_SOURCE_OUTDIR=str(complete),
               CAAS_REUSE_OUTDIR=str(strict), FADE_SOURCE_OUTDIR='', RER_SOURCE_OUTDIR='',
               PRECOMP_USE_DISCOVERY='true', RUN_FADE='false', RUN_RER='false')
    check = '\nprintf "%s\\n" "$DISCOVERY_FROM" "$META_CAAS_FROM" "$BACKGROUND_INPUT" "$CAAS_NULL_DIR"\n'
    run = subprocess.run(['bash', '-c', 'set -Eeuo pipefail\n' + block + check],
                         env=env, capture_output=True, text=True)
    assert run.returncode == 0, run.stdout + run.stderr
    assert all(line.startswith(str(strict) + '/') for line in run.stdout.splitlines()), run.stdout
