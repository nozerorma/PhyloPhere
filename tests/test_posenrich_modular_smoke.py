"""Smoke test dirigido: el subworkflow POSENRICH_MODULAR corre sobre fixtures (Nextflow local)."""
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tests'))
from test_posenrich_registry import make_fixture  # noqa: E402

NEXTFLOW = shutil.which('nextflow')


SMOKE_MAIN = """\
nextflow.enable.dsl=2
// Main script generado por tests/test_posenrich_modular_smoke.py: vive en la raíz del
// repo para que baseDir de los procesos apunte a la raíz. No es un archivo versionado.
include { POSENRICH_MODULAR } from './subworkflows/ENRICHMENT/posenrich_modular.nf'

params.position_registry_memory = '500 MB'
params.position_module_memory = '500 MB'
params.fixture_dir = null
params.outdir = 'results'
params.traitname = 'smoke'
params.posenrich_ranking = 'caas_score'
params.posenrich_aggregation = 'mean'
params.posenrich_n_permutations = 199
params.posenrich_seed = 7
params.posenrich_radios = '0,5,10'

workflow {
    assert params.fixture_dir : "POSENRICH_MODULAR smoke: --fixture_dir is required"
    def fx = params.fixture_dir
    POSENRICH_MODULAR(
        Channel.value(file("${fx}/position_scores.tsv")),
        Channel.value(file("${fx}/filtered_discovery.tsv")),
        Channel.value(file("${fx}/null_universe_summary.tsv")),
        Channel.value(file("${fx}/maps")),
        Channel.value(file("${fx}/cleaned_background.txt")),
        Channel.value(file('NO_ANNOT_ENSEMBL')),
        Channel.value(file('NO_ANNOT_DOMAIN')),
        Channel.value(file('NO_ANNOT_UCR')),
        Channel.value(file('NO_ANNOT_FUBAR')),
        Channel.value(file('NO_FADE_TOP')),
        Channel.value(file('NO_FADE_BOTTOM')),
        Channel.value(file("${fx}/custom.tsv"))
    )
}
"""


@pytest.mark.skipif(NEXTFLOW is None, reason='nextflow no disponible en PATH')
def test_posenrich_modular_nextflow_smoke(tmp_path):
    files, maps = make_fixture(tmp_path)
    fixture = tmp_path / 'fixture'
    fixture.mkdir()
    shutil.copy(files['scores'], fixture / 'position_scores.tsv')
    shutil.copy(files['postproc'], fixture / 'filtered_discovery.tsv')
    shutil.copy(files['null_universe'], fixture / 'null_universe_summary.tsv')
    shutil.copy(files['background'], fixture / 'cleaned_background.txt')
    shutil.copy(files['custom'], fixture / 'custom.tsv')
    shutil.copytree(maps, fixture / 'maps')

    outdir = tmp_path / 'out'
    work = tmp_path / 'work'
    # Main script temporal en la raíz del repo (baseDir de los procesos = raíz).
    smoke_main = ROOT / '.smoke_posenrich_modular.nf'
    smoke_main.write_text(SMOKE_MAIN)
    try:
        cmd = [NEXTFLOW, 'run', str(smoke_main),
               '--fixture_dir', str(fixture), '--outdir', str(outdir), '-w', str(work),
               '-ansi-log', 'false']
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(ROOT),
                                timeout=280)
    finally:
        smoke_main.unlink(missing_ok=True)
    log = result.stdout + result.stderr
    assert result.returncode == 0, log[-4000:]
    assert 'POSENRICH_MODULAR' in log

    pos = outdir / 'posenrich'
    for rel in ('registry/positions.tsv', 'registry/position_evidence.tsv', 'registry/members.tsv',
                'registry/collections.tsv', 'registry/coverage.tsv', 'registry/pending_mappings.tsv',
                'registry/manifest.json',
                'characterization/position_characterization.tsv', 'characterization/collection_summary.tsv',
                'characterization/ranking_summary.tsv', 'characterization/universe_summary.tsv',
                'localization/localization_results.tsv', 'localization/localization_by_gene.tsv',
                'localization/localization_manifest.json',
                'global/global_results.tsv', 'global/global_by_gene.tsv',
                'proximity/proximity_protein_distances.tsv', 'proximity/proximity_protein_summary.tsv',
                'proximity/proximity_genomic_distances.tsv', 'proximity/proximity_genomic_summary.tsv'):
        assert (pos / rel).exists(), rel
    html = list(pos.glob('POSENRICH_modular_report_*.html')) + \
        list((outdir / 'html_reports').glob('POSENRICH_modular_report_*.html'))
    assert html, 'informe modular no publicado'
    # El GMT export del registro es compatible con build_position_gmt.py.
    assert (pos / 'registry/gmts/pfam_domains.gmt').exists()


@pytest.mark.skipif(NEXTFLOW is None, reason='nextflow no disponible en PATH')
@pytest.mark.parametrize('modular, legacy', [(True, None), (False, None), (False, 'plain'), (False, 'gzip'), (False, 'bad_checksum'), (False, 'wrong_alignment')])
def test_enrichment_routes_canonical_characterization(tmp_path, modular, legacy):
    files, maps = make_fixture(tmp_path)
    main = ROOT / '.smoke_posenrich_enrichment.nf'
    out = tmp_path / 'out'
    work = tmp_path / 'work'
    script = """
nextflow.enable.dsl=2
include { ENRICHMENT } from './workflows/enrichment.nf'
params.traitname = 'route'
params.fcs_enabled = false
params.scoring_ami = false
params.scoring_string = false
params.posenrich = false
workflow {
    def missing_id = 0
    def absent = { missing_id++; Channel.value(file("NO_FILE_${missing_id}")) }
    ENRICHMENT(
        absent(), absent(), absent(), absent(), absent(), absent(),
        Channel.value(file(params.fixture_background)),
        absent(), absent(), absent(), absent(), absent(),
        Channel.value(file(params.fixture_scores)),
        Channel.value(file(params.fixture_discovery)),
        Channel.value(file(params.fixture_null)),
        absent(), absent(), absent(), absent(), absent(), absent(),
        absent(), absent(), absent(), absent(), absent(), absent(), absent()
    )
}
"""
    main.write_text(script)
    try:
        cmd = [NEXTFLOW, 'run', str(main), '--outdir', str(out), '-w', str(work),
               '--position_annotation_memory', '500 MB', '--position_registry_memory', '500 MB', '--position_module_memory', '500 MB', '--fcs_enabled', 'false', '--scoring_ami', 'false', '--scoring_string', 'false',
               '--posenrich', 'false', '--posenrich_modular', str(modular).lower(),
               '--fixture_background', files['background'], '--fixture_scores', files['scores'],
               '--fixture_discovery', files['postproc'], '--fixture_null', files['null_universe'],
               '--caas_map_dir', maps, '--posenrich_custom_file', files['custom'],
               '-ansi-log', 'false']
        if legacy:
            legacy_path = tmp_path / ('legacy.tsv.gz' if legacy == 'gzip' else 'legacy.tsv')
            text = 'gene\tpfam_id\ttarget_name\tdescription\tclan_acc\tclan_name\tali_start\tali_end\n'
            if legacy in ('bad_checksum', 'wrong_alignment'):
                import json
                text = text.rstrip('\n')+'\tcoordinate_space\n'
                legacy_path.write_text(text)
                meta = {'table_sha256': '0'*64} if legacy == 'bad_checksum' else {}
                if legacy == 'wrong_alignment':
                    original = tmp_path/'original_alignments';original.mkdir()
                    actual = tmp_path/'actual_alignments';actual.mkdir()
                    meta['alignment_dir'] = str(original)
                    cmd += ['--alignment', str(actual)]
                Path(str(legacy_path)+'.manifest.json').write_text(json.dumps(meta))
            elif legacy == 'gzip':
                import gzip
                with gzip.open(legacy_path, 'wt') as handle:
                    handle.write(text)
            else:
                legacy_path.write_text(text)
            cmd += ['--posenrich_domains', 'true', '--domain_variability_file', str(legacy_path)]
        result = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, timeout=280)
    finally:
        main.unlink(missing_ok=True)
    log = result.stdout + result.stderr
    if legacy:
        assert result.returncode != 0
        expected = 'checksum mismatch' if legacy == 'bad_checksum' else 'belongs to alignment' if legacy == 'wrong_alignment' else 'has no coordinate_space'
        assert expected in log
        assert 'Submitted process' not in log
        return
    assert result.returncode == 0, log[-5000:]
    assert ('POSITION_ANNOTATIONS' in log) is (not modular)
    assert list((out / 'html_reports').glob('15.Comparison_report_*.html'))
    compare_tasks = [p.parent for p in work.rglob('.command.sh') if '15.Comparison_report.Rmd' in p.read_text()]
    assert len(compare_tasks) == 1
    staged = compare_tasks[0] / 'position_characterization.tsv'
    assert staged.is_file()
    if modular:
        assert staged.resolve().parent.name == 'characterization'
        assert staged.read_bytes() == (out / 'posenrich/characterization/position_characterization.tsv').read_bytes()
        import csv
        with (out / 'posenrich/localization/localization_results.tsv').open() as f:
            rows = list(csv.DictReader(f, delimiter='\t'))
        assert {r['rank_direction'] for r in rows} == {'global', 'top', 'bottom'}
        assert all(r['p_mc'] == 'NA' and r['B_valid'] == '0' for r in rows)
