"""Contrast gate schema/metadata tests and real Nextflow resume/publication checks."""
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'subworkflows/CT/local/scripts/check_contrasts.py'
spec = importlib.util.spec_from_file_location('contrast_gate', SCRIPT)
gate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)


def make_design(project, n=3):
    text = ''.join(f'F{i}\t1\t{i}\nB{i}\t0\t{i}\n' for i in range(1, n + 1))
    (project / 'traits.tab').write_text(text)
    directory = project / 'traits_dir'
    directory.mkdir(exist_ok=True)
    (directory / 'traitfile_H1.tab').write_text(text)
    (directory / 'traitfile_H2.tab').write_text(text.replace('F', 'X').replace('B', 'Y'))
    rows = ['hypothesis_id\tpair\tspecies1\tspecies2\tpss_score']
    for hyp, fg, bg in [('H1', 'F', 'B'), ('H2', 'X', 'Y')]:
        rows.extend(f'{hyp}\t{i}\t{fg}{i}\t{bg}{i}\t1' for i in range(1, n + 1))
    (directory / 'contrast_hypotheses_pairs.tsv').write_text('\n'.join(rows) + '\n')
    (directory / 'contrast_hypotheses_summary.tsv').write_text(
        f'hypothesis_id\tnum_pairs\nH1\t{n}\nH2\t{n}\n')
    (project / 'perm.tab').write_text('species\tscore\nF1\t1\nB1\t0\n')
    return dict(traitfile='traits.tab', permulation_traitfile='perm.tab',
                trait_dir='traits_dir', min_contrasts=3, traitname='score', multi_hypothesis=False)


@pytest.mark.parametrize('multi', [False, True])
def test_metadata_exactly_matches_active_hypotheses(tmp_path, monkeypatch, multi):
    args = make_design(tmp_path)
    args['multi_hypothesis'] = multi
    monkeypatch.chdir(tmp_path)
    gate.gate(args)
    output = tmp_path / 'traitfiles_ok_dir'
    expected = {'H1', 'H2'} if multi else {'H1'}
    assert {p.stem.removeprefix('traitfile_') for p in output.glob('traitfile_H*.tab')} == expected
    for name in ['contrast_hypotheses_pairs.tsv', 'contrast_hypotheses_summary.tsv']:
        with (output / name).open() as fh:
            assert {row['hypothesis_id'] for row in gate.csv.DictReader(fh, delimiter='\t')} == expected


@pytest.mark.parametrize('minimum', [0, -1, 'oops', 1.5, True])
def test_invalid_thresholds_fail(tmp_path, monkeypatch, minimum):
    args = make_design(tmp_path); args['min_contrasts'] = minimum
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError, match='positive integer'):
        gate.gate(args)
    assert not (tmp_path / 'traitfile_ok.tab').exists()


@pytest.mark.parametrize('bad', ['F1\t1\t1\n', 'F1\t1\t1\nF1\t0\t1\n',
                               'F1\t1\t2\nB1\t0\t2\n'])
def test_invalid_pairs_fail(tmp_path, bad):
    path = tmp_path / 'bad.tab'; path.write_text(bad)
    with pytest.raises(ValueError):
        gate.read_design(path)


def test_metadata_species_mismatch_fails(tmp_path, monkeypatch):
    args = make_design(tmp_path)
    path = tmp_path / 'traits_dir/contrast_hypotheses_pairs.tsv'
    path.write_text(path.read_text().replace('F1', 'WRONG'))
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError, match='mismatch'):
        gate.gate(args)


def test_empty_design_reaches_the_gate_as_a_valid_skip(tmp_path, monkeypatch):
    args = make_design(tmp_path, 0)
    monkeypatch.chdir(tmp_path)
    gate.gate(args)
    status = json.loads((tmp_path / 'gate_status.json').read_text())
    assert status['skipped'] and status['n_contrasts'] == 0
    assert (tmp_path / 'low_contrasts.skip').exists()
    assert not (tmp_path / 'traitfile_ok.tab').exists()


@pytest.fixture
def project(tmp_path):
    if not shutil.which('nextflow'):
        pytest.skip('Run in the phylophere environment')
    scripts = tmp_path / 'subworkflows/CT/local/scripts'
    scripts.mkdir(parents=True); shutil.copy2(SCRIPT, scripts)
    (tmp_path / 'main.nf').write_text(
        f"include {{ GATE_CONTRASTS }} from '{ROOT}/subworkflows/CT/ct_check_min_contrasts.nf'\n"
        "workflow {\n"
        "  result = GATE_CONTRASTS(Channel.value(file('traits.tab')),Channel.value(file('perm.tab')),Channel.value(file('traits_dir')))\n"
        "  result.skip_flag.view { f -> assert file(\"${params.outdir}/low_contrasts.skip\").exists(); 'skip-ready' }\n"
        "  result.traitfile_out.view { f -> assert !file(\"${params.outdir}/low_contrasts.skip\").exists(); 'pass-ready' }\n"
        "}\n")
    (tmp_path / 'audit.config').write_text(
        "process.executor='local'\nparams { outdir='out'; min_contrasts=3; traitname='score'; multi_hypothesis=false }\n")
    return tmp_path


def run_project(project, resume=False):
    command = ['nextflow', '-log', 'nextflow.log', 'run', 'main.nf', '-c', 'audit.config', '-ansi-log', 'false']
    if resume:
        command += ['-resume']
    result = subprocess.run(command, cwd=project, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr
    return result.stdout + result.stderr


def test_sentinel_reconciles_before_emission_including_resume(project):
    make_design(project, 2)
    assert 'skip-ready' in run_project(project)
    sentinel = project / 'out/low_contrasts.skip'
    assert sentinel.exists()
    sentinel.unlink()
    assert 'skip-ready' in run_project(project, resume=True)
    assert sentinel.exists()  # cached task must restore it before consumers run
    make_design(project, 3)
    assert 'pass-ready' in run_project(project)
    assert not sentinel.exists()
    sentinel.write_text('stale')
    assert 'pass-ready' in run_project(project, resume=True)
    assert not sentinel.exists()
    output = next(d for d in (project / 'work').glob('*/*/traitfiles_ok_dir'))
    assert not (output / 'traitfile_H2.tab').exists()
    assert '\nH2\t' not in (output / 'contrast_hypotheses_pairs.tsv').read_text()

# Reuse only the fixture construction from REPORTING; these tests exercise the
# CONTRAST_SELECTION composition and its newly wrapped gate.
from test_reporting_workflow import reporting_project


@pytest.mark.parametrize(('reporting', 'prune'), [(False, False), (True, False), (False, True), (True, True)])
def test_contrast_workflow_composes_all_reporting_pruning_branches(reporting_project, reporting, prune):
    project = reporting_project
    local = project / 'subworkflows/TRAIT_ANALYSIS/local'
    bodies = {
        '3.CI-composition.Rmd': '''
stopifnot(as.character(params$pss_top_pct)=='0')
stopifnot(file.exists(file.path(params$output_dir,'dataset_marker.txt')))
dir.create(file.path(params$output_dir,'1.Data-exploration/5.CI_overlaps'),recursive=TRUE,showWarnings=FALSE)
writeLines('ci',file.path(params$output_dir,'ci_marker.txt'))
saveRDS(list(),file.path(params$output_dir,'1.Data-exploration/5.CI_overlaps/observed_selection.rds'))
writeLines('F1\\t1\\t1\\nB1\\t0\\t2',file.path(params$output_dir,'1.Data-exploration/5.CI_overlaps/candidate_species.tab'))
''',
        '4.Independent_contrasts.Rmd': '''
stopifnot(as.character(params$pss_top_pct)=='0')
stopifnot(file.exists(file.path(params$output_dir,'ci_marker.txt')))
td<-file.path(params$output_dir,'2.CT/1.Traitfiles');dir.create(td,recursive=TRUE)
writeLines(c('F1\\t1\\t1','B1\\t0\\t1','F2\\t1\\t2','B2\\t0\\t2','F3\\t1\\t3','B3\\t0\\t3'),file.path(td,'traitfile.tab'))
file.copy(file.path(td,'traitfile.tab'),file.path(td,'traitfile_H1.tab'))
saveRDS(list(),file.path(td,'observed_selection.rds'))
pd<-file.path(params$output_dir,'2.CT/2.Permulation_traitfiles');dir.create(pd,recursive=TRUE)
writeLines('species\\tscore\\nF1\\t1',file.path(pd,'permulation_traitfile.tab'))
writeLines('contrast',file.path(params$output_dir,'contrast_marker.txt'))
''',
    }
    for name, body in bodies.items():
        source = (ROOT / 'subworkflows/TRAIT_ANALYSIS/local' / name).read_text()
        header = source[:source.index('\n---', 4) + 4]
        (local / name).write_text(header + '\n```{r}\n' + body + '\n```\n')
    scripts = project / 'subworkflows/CT/local/scripts'
    scripts.mkdir(parents=True)
    for name in ['check_contrasts.py', 'lean_contrast_selector.R', 'pss_core.R']:
        shutil.copy2(SCRIPT.parent / name, scripts)
    (project / 'main.nf').write_text(
        f"include {{ CONTRAST_SELECTION }} from '{ROOT}/workflows/contrast_selection.nf'\n"
        'workflow { CONTRAST_SELECTION(Channel.value(file(params.trait_file)),Channel.value(file(params.tree_file))) }\n')
    with (project / 'audit.config').open('a') as fh:
        fh.write('\nparams.min_contrasts=3\nparams.multi_hypothesis=false\nparams.max_fop=3\n')
    command = ['nextflow', '-log', 'contrast.log', 'run', 'main.nf', '-c', 'audit.config',
               '-ansi-log', 'false', '--reporting', str(reporting).lower(),
               '--prune_data', str(prune).lower(), '--pss_top_pct', '0']
    result = subprocess.run(command, cwd=project, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr
    tasks = {name: [] for name in ['dataset', 'phenotype', 'ci', 'contrast']}
    for script in (project / 'work').glob('*/*/.command.sh'):
        text = script.read_text()
        for stage, report in [('dataset', '1.Dataset_exploration.Rmd'), ('phenotype', '2.Phenotype_exploration.Rmd'),
                              ('ci', '3.CI-composition.Rmd'), ('contrast', '4.Independent_contrasts.Rmd')]:
            if report in text:
                tasks[stage].append(script.parent)
    assert len(tasks['dataset']) == len(tasks['ci']) == len(tasks['contrast']) == 1
    assert len(tasks['phenotype']) == int(reporting)
    assert not (tasks['dataset'][0] / 'data_exploration/ci_marker.txt').exists()
    assert not (tasks['ci'][0] / 'data_exploration/contrast_marker.txt').exists()
    assert (tasks['contrast'][0] / 'data_exploration/contrast_marker.txt').exists()
    assert not (project / 'out/low_contrasts.skip').exists()
