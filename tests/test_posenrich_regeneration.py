"""Modular report regeneration recovers saved component provenance."""
from pathlib import Path
from test_report_regeneration import _make_outdir

ROOT = Path(__file__).resolve().parents[1]


def test_modular_regeneration_does_not_invent_settings(tmp_path):
    from gui.generation.report_registry import detect_reports, build_script
    out = _make_outdir(tmp_path)
    report = next(r for r in detect_reports(str(out), ROOT) if r.id == 'posenrich_modular')
    script = build_script(report, ROOT, use_singularity=False)
    for invented in ('ranking              =', 'aggregation          =',
                     'n_permutations       =', 'seed                 =', 'radios               ='):
        assert invented not in script
    template = (ROOT / 'subworkflows/ENRICHMENT/local/posenrich_modular_report.Rmd').read_text()
    assert 'loc_manifest$parameters' in template and 'prox_manifest$parameters' in template
