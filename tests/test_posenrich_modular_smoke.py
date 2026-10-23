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
