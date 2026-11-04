"""Precomputed RER and FADE sources are chosen by the files they are read from, in the template and its Python mirror."""
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from gui.models.precomputed import PrecomputedConfig, derive_paths  # noqa: E402

TEMPLATE = (ROOT / 'gui/generation/templates/run_single.sh.j2').read_text()
TRAIT = 'trait'


def chunk(start, end):
    """Shell text between two literal markers of the template, markers included at the start only."""
    i = TEMPLATE.index(start)
    return TEMPLATE[i:TEMPLATE.index(end, i)]


SELECT = chunk('FADE_SOURCE_OUTDIR=""\nRER_SOURCE_OUTDIR=""', 'PRECOMP_OUTDIR="$PASS_SOURCE_OUTDIR"')
RER_BLOCK = chunk('if [ "${PRECOMP_USE_RER:-false}" = true ]; then', 'if [ "${PRECOMP_USE_FADE:-false}" = true ]; then')


def touch(root, *parts):
    p = Path(root, *parts)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text('x')
    return p


def raw_run(root):
    """Outdir of a run that computed RER."""
    touch(root, 'rerconverge/rer_results', f'{TRAIT}.continuous.output')
    touch(root, 'rerconverge/rer_results', f'{TRAIT}.continuous.perms.rds/blk.rds')
    touch(root, 'rerconverge/rer_results', 'rerconverge_summary_trait.tsv')
    return root


def validated_run(root, mode='continuous'):
    """Outdir of a run that reused a precomputed RER: only the validated result and its null are published."""
    touch(root, 'rerconverge/rer_results/rer_inputs', 'validated_rer.output')
    touch(root, 'rerconverge/rer_results/normalized_null', 'blk_000001.rds')
    touch(root, 'rerconverge/rer_results', 'rerconverge_summary_trait.tsv')
    Path(root, 'rerconverge/rer_results/rer_run_metadata.json').write_text(f'{{"analysis_mode":"{mode}"}}')
    return root


def bash(script, env):
    return subprocess.run(['bash', '-c', 'set -u\n' + script], env=env, capture_output=True, text=True)


def select(normal_dirs, precomp):
    arr = ' '.join(f'"{d}"' for d in normal_dirs)
    out = bash(f'NORMAL_RUN_OUTDIRS=({arr}); TRAIT={TRAIT}; PRECOMP_OUTDIR="{precomp}"\n{SELECT}\n'
               'echo "RER=$RER_SOURCE_OUTDIR"; echo "FADE=$FADE_SOURCE_OUTDIR"', {})
    assert out.returncode == 0, out.stderr
    return dict(l.split('=', 1) for l in out.stdout.splitlines())


def rer_vars(source):
    out = bash(f'TRAIT={TRAIT}; PRECOMP_USE_RER=true; RER_SOURCE_OUTDIR="{source}"; PRECOMP_OUTDIR=/none\n{RER_BLOCK}\n'
               'echo "CONT=${RER_CONTINUOUS_FILE:-}"; echo "BIN=${RER_BINARY_FILE:-}"; echo "PERMS=${RER_PERMS_FILE:-}"; '
               'echo "SUMMARY=${SCORING_RER_INPUT:-}"', {})
    return out, dict(l.split('=', 1) for l in out.stdout.splitlines() if '=' in l)


def test_selection_skips_a_run_with_an_empty_rer_directory(tmp_path):
    # The complete run has the directory (report and tables) but neither the raw nor the validated result.
    complete = tmp_path / 'complete'
    touch(complete, 'rerconverge/rer_results', '5.RERconverge_report.html')
    source = raw_run(tmp_path / 'def2')
    assert select([complete], source)['RER'] == str(source)


def test_selection_prefers_the_normal_run_and_accepts_the_validated_result(tmp_path):
    complete = validated_run(tmp_path / 'complete')
    source = raw_run(tmp_path / 'def2')
    assert select([complete], source)['RER'] == str(complete)


def test_selection_finds_fade_tables_only_where_the_summary_exists(tmp_path):
    complete = tmp_path / 'complete'
    touch(complete, 'selection/fade/top', '6.FADE_report_top.html')
    source = tmp_path / 'def2'
    touch(source, 'selection/fade/top', 'fade_summary_top.tsv')
    assert select([complete], source)['FADE'] == str(source)


def test_selection_is_empty_when_no_run_has_the_files(tmp_path):
    assert select([tmp_path / 'a'], tmp_path / 'b') == {'RER': '', 'FADE': ''}


def test_block_takes_the_raw_result_and_its_perms(tmp_path):
    out, v = rer_vars(raw_run(tmp_path / 'run'))
    assert out.returncode == 0, out.stderr
    base = tmp_path / 'run/rerconverge/rer_results'
    assert v['CONT'] == str(base / f'{TRAIT}.continuous.output')
    assert v['PERMS'] == str(base / f'{TRAIT}.continuous.perms.rds')
    assert v['SUMMARY'] == str(base / 'rerconverge_summary_trait.tsv')


def test_block_takes_the_validated_result_and_its_null(tmp_path):
    out, v = rer_vars(validated_run(tmp_path / 'run'))
    assert out.returncode == 0, out.stderr
    base = tmp_path / 'run/rerconverge/rer_results'
    assert v['CONT'] == str(base / 'rer_inputs/validated_rer.output') and v['BIN'] == ''
    assert v['PERMS'] == str(base / 'normalized_null')


def test_block_routes_a_validated_binary_result(tmp_path):
    out, v = rer_vars(validated_run(tmp_path / 'run', mode='binary'))
    assert out.returncode == 0, out.stderr
    assert v['BIN'].endswith('rer_inputs/validated_rer.output') and v['CONT'] == ''


def test_block_leaves_perms_empty_without_a_null(tmp_path):
    run = validated_run(tmp_path / 'run')
    (run / 'rerconverge/rer_results/normalized_null/blk_000001.rds').unlink()
    (run / 'rerconverge/rer_results/normalized_null').rmdir()
    out, v = rer_vars(run)
    assert out.returncode == 0, out.stderr
    assert v['PERMS'] == ''


def test_block_stops_when_nothing_is_there(tmp_path):
    out, _ = rer_vars(tmp_path / 'empty')
    assert out.returncode == 1 and 'precomputed RER' in out.stderr


def test_python_mirror_follows_the_same_rule(tmp_path):
    cfg = PrecomputedConfig(base_path=str(tmp_path), use_rer=True)
    validated_run(tmp_path / TRAIT)
    paths = {label: p for label, p, _ in derive_paths(cfg, TRAIT)}
    assert paths['rer_continuous_file / rer_binary_file'].endswith('rer_inputs/validated_rer.output')
    assert paths['rer_perms_file / scoring_rer_perms_input'].endswith('normalized_null')
    raw_run(tmp_path / 'other')
    labels = [label for label, _, _ in derive_paths(cfg, 'other')]
    assert 'rer_continuous_file' in labels
