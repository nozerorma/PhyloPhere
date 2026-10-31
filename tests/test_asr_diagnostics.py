"""Report 9 contract, full renders, safe archives and real Nextflow resume."""
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile

import pytest

ROOT = Path(__file__).resolve().parents[1]
LOCAL = Path('subworkflows/ASR_DIAGNOSTICS/local')
STYLE = Path('subworkflows/TRAIT_ANALYSIS/local/src/report_style.R')
spec = importlib.util.spec_from_file_location('asr_extract', ROOT / LOCAL / 'src/extract_diagnostics.py')
extractor = importlib.util.module_from_spec(spec)
spec.loader.exec_module(extractor)


def archive(path, entries):
    with tarfile.open(path, 'w:gz') as out:
        for name, kind in entries:
            item = tarfile.TarInfo(name)
            if kind == 'file':
                item.size = 4
                out.addfile(item, io.BytesIO(b'data'))
            else:
                item.type = tarfile.SYMTYPE if kind == 'symlink' else tarfile.LNKTYPE
                item.linkname = '/outside'
                out.addfile(item)
    return path


@pytest.mark.parametrize('prefix', ['', 'ct_disambiguation/', './out/nested/ct_disambiguation/'])
def test_archive_layouts(tmp_path, prefix):
    path = archive(tmp_path / "input's archive.tgz", [(prefix + name, 'file') for name in
        ['caas_convergence_master.csv', 'observed_scoring.json']])
    extractor.extract_contract(path, tmp_path / 'selected')
    assert {p.name for p in (tmp_path / 'selected').iterdir()} == {
        'caas_convergence_master.csv', 'observed_scoring.json'}


@pytest.mark.parametrize('extra', [
    [('../outside', 'file')], [('/outside', 'file')], [('link', 'symlink')], [('link', 'hardlink')],
    [('caas_convergence_master.csv', 'file')],
    [('second/caas_convergence_master.csv', 'file'), ('second/observed_scoring.json', 'file')],
])
def test_unsafe_or_ambiguous_archives(tmp_path, extra):
    path = archive(tmp_path / 'bad.tgz', [('caas_convergence_master.csv', 'file'),
                                         ('observed_scoring.json', 'file'), *extra])
    with pytest.raises(ValueError):
        extractor.extract_contract(path, tmp_path / 'selected')
    assert not (tmp_path / 'selected').exists()


def test_incomplete_archive(tmp_path):
    path = archive(tmp_path / 'bad.tgz', [('caas_convergence_master.csv', 'file')])
    with pytest.raises(ValueError, match='exactly one'):
        extractor.extract_contract(path, tmp_path / 'selected')


def test_contract_and_full_renders(tmp_path):
    result = subprocess.run(['Rscript', str(ROOT / 'tests/test_asr_diagnostics.R'), str(tmp_path / 'renders')],
                            cwd=ROOT, capture_output=True, text=True, timeout=180)
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.count('PASS full render:') == 5


@pytest.fixture
def project(tmp_path):
    assert shutil.which('nextflow') and shutil.which('Rscript'), 'Use conda run -n phylophere'
    project = tmp_path / 'project with spaces'
    for rel in [Path('workflows/asr_diagnostics.nf'), Path('subworkflows/ASR_DIAGNOSTICS/asr_diagnostics.nf'), STYLE]:
        dest = project / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(ROOT / rel, dest)
    shutil.copytree(ROOT / LOCAL, project / LOCAL)
    data = project / "observed's data"
    data.mkdir()
    (data / 'caas_convergence_master.csv').write_text('gene,msa_pos,side,asr_path_score\n001.2,0,top,0.8\n')
    genes = {'001.2': {'positions': 1, 'sha256': hashlib.sha256(b'0\n').hexdigest()}}
    batch = dict(version=1, complete=True, score_rule='caas_pos_species_pss_v4', genes=genes,
                 excluded_genes=[], posterior_threshold=0.2, model='lg', tree_sha256='tree',
                 taxid_sha256=None, weights_sha256=None, weighted=False)
    (data / 'observed_scoring.json').write_text(json.dumps(dict(version=1, complete=True,
        genes=genes, excluded_genes=[], batches=[batch])))
    (project / 'main.nf').write_text("""
include { ASR_DIAGNOSTICS } from './workflows/asr_diagnostics.nf'
workflow { ASR_DIAGNOSTICS(params.live ? Channel.value(file(params.disambiguation_dir)) : null) }
""")
    (project / 'audit.config').write_text("""
process.executor = 'local'
params {
  disambiguation_dir = "observed's data"
  outdir = 'out'; live = false
  use_singularity = false; use_apptainer = false
  ct_disambig_posterior_threshold = 0.1
}
""")
    return project


def run(project, *extra):
    result = subprocess.run(['nextflow', '-log', 'nextflow.log', 'run', 'main.nf', '-c', 'audit.config',
        '-w', 'work', '-ansi-log', 'false', *extra], cwd=project, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr


def report_tasks(project):
    return [p.parent for p in (project / 'work').glob('*/*/.command.sh')
            if "rmarkdown::render('9.ASR_diagnostics.Rmd'" in p.read_text()]


def test_render_paths_and_resume(project):
    run(project, '--live', 'true')
    tasks = report_tasks(project)
    assert len(tasks) == 1
    assert (project / 'out/html_reports/9.ASR_diagnostics.html').exists()
    assert (project / 'out/asr_diagnostics/tsv/asr_path_summary_statistics.tsv').exists()
    payload = json.loads((tasks[0] / 'report_params.json').read_text())
    assert "observed's data" in payload['disambig_dir']
    run(project, '--live', 'true', '-resume')
    assert len(report_tasks(project)) == 1
    for n, rel in enumerate([LOCAL / '9.ASR_diagnostics.Rmd', LOCAL / 'src/diagnostics_io.R', STYLE], start=2):
        source = project / rel
        source.write_text(source.read_text() + '\n<!-- changed -->\n' if source.suffix == '.Rmd'
                          else source.read_text() + '\n# changed\n')
        run(project, '--live', 'true', '-resume')
        assert len(report_tasks(project)) == n
    # Detect actual content changes inside directory inputs, even at identical size/mtime.
    source = project / "observed's data/caas_convergence_master.csv"
    stat = source.stat()
    source.write_text(source.read_text().replace('0.8', '0.7'))
    os.utime(source, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    run(project, '--live', 'true', '-resume')
    assert len(report_tasks(project)) == 5
    summary = (project / 'out/asr_diagnostics/tsv/asr_path_summary_statistics.tsv').read_text()
    assert '0.7' in summary


def test_nested_archive_workflow(project):
    path = project / "observed's archive.tgz"
    with tarfile.open(path, 'w:gz') as out:
        out.add(project / "observed's data", arcname='run/nested/ct_disambiguation')
    run(project, '--disambiguation_dir', path.name)
    assert (project / 'out/html_reports/9.ASR_diagnostics.html').exists()
    assert len(report_tasks(project)) == 1


def test_regeneration_stages_new_dependency(tmp_path):
    from gui.generation.report_registry import detect_reports, build_script
    out = tmp_path / 'out'
    (out / 'html_reports').mkdir(parents=True)
    (out / 'html_reports/9.ASR_diagnostics.html').touch()
    (out / 'ct_disambiguation').mkdir()
    reports = {r.id: r for r in detect_reports(str(out), ROOT)}
    script = build_script(reports['asr_diagnostics'], ROOT, use_singularity=False)
    assert 'report_style.R' in script
    result = subprocess.run(['bash','-n'], input=script, text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
