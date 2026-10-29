"""Report 13 regeneration recovers the published component directories and invents no setting."""
from pathlib import Path
from test_report_regeneration import _make_outdir

ROOT = Path(__file__).resolve().parents[1]
RMD = ROOT / 'subworkflows/POSENRICH/local/13.POSENRICH_report.Rmd'


def _script(tmp_path, drop=()):
    from gui.generation.report_registry import detect_reports, build_script
    out = _make_outdir(tmp_path)
    for name in drop:
        for path in sorted((out / 'posenrich' / name).rglob('*'), reverse=True):
            path.rmdir() if path.is_dir() else path.unlink()
        (out / 'posenrich' / name).rmdir()
    report = next(r for r in detect_reports(str(out), ROOT) if r.id == 'posenrich')
    return report, build_script(report, ROOT, use_singularity=False)


def test_regeneration_passes_every_component_directory(tmp_path):
    report, script = _script(tmp_path)
    assert report.runnable and report.traitname == 'trait'
    for param in ('universe_dir', 'registry_dir', 'characterization_dir', 'accumulation_dir', 'within_gene_dir',
                  'proximity_dir', 'traitname'):
        assert param in script
    assert "13.POSENRICH_report_trait.html" in script and 'report_presentation.R' in script
    for invented in ('fdr_threshold', 'ranking', 'aggregation', 'n_permutations', 'seed', 'radios'):
        assert f'{invented} ' not in script.split('rmarkdown::render')[1]


def test_regeneration_without_levels_still_runs(tmp_path):
    # Level 1, level 2 and proximity are optional: the report states them as not available.
    report, script = _script(tmp_path, drop=('accumulation/top', 'accumulation/bottom', 'accumulation/all', 'accumulation',
                                             'within_gene', 'proximity'))
    assert report.runnable
    assert 'within_gene_dir      = NULL' in script and 'proximity_dir        = NULL' in script


def test_report_template_reads_the_published_manifests():
    template = RMD.read_text()
    assert 'wg_manifest$parameters' in template and 'universe_manifest' in template
