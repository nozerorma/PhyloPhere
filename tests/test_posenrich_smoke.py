"""Smoke tests dirigidos de Nextflow: el workflow POSENRICH corre sobre fixtures y el informe 15
de ENRICHMENT acepta la caracterización de POSENRICH como entrada opcional."""
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tests'))
from test_posenrich_registry import make_fixture, run_registry  # noqa: E402

NEXTFLOW = shutil.which('nextflow')
CHARACTERIZE = ROOT / 'subworkflows/POSENRICH/local/src/posenrich_characterize.py'


# Los main scripts viven en la raíz del repo para que baseDir de los procesos apunte a la raíz.
# No son archivos versionados.
POSENRICH_MAIN = """\
nextflow.enable.dsl=2
include { POSENRICH } from './workflows/posenrich.nf'

params.position_registry_memory = '500 MB'
params.position_module_memory = '500 MB'
params.fixture_dir = null
params.outdir = 'results'
params.traitname = 'smoke'
params.posenrich_radios = '0,5,10'
params.null_source = 'live'     // live: null positions and blocks go in as channels; prior: --caas_null_dir holds them
params.level1 = 'true'         // level 1 (Accumulation) on the universe table

workflow {
    assert params.fixture_dir : "POSENRICH smoke: --fixture_dir is required"
    def fx = params.fixture_dir
    def live = params.null_source == 'live'
    POSENRICH(
        params.level1.toString() == 'true',
        Channel.value(file("${fx}/filtered_discovery.tsv")),
        Channel.value(file("${fx}/cleaned_background.txt")),
        Channel.value(file("${fx}/background.output")),
        live ? Channel.value(file("${fx}/caas_permulation/null_positions.parquet")) : Channel.empty(),
        live ? Channel.value(files("${fx}/caas_permulation/blocks/*", type: 'dir').sort()) : Channel.empty(),
        Channel.value(file("${fx}/position_scores.tsv")),
        Channel.value(file("${fx}/filtered_discovery.tsv")),
        Channel.value(file("${fx}/null_universe_summary.tsv")),
        Channel.empty(), Channel.empty()
    )
}
"""

ENRICHMENT_MAIN = """\
nextflow.enable.dsl=2
include { ENRICHMENT } from './workflows/enrichment.nf'
params.traitname = 'route'
workflow {
    def missing_id = 0
    def absent = { missing_id++; Channel.value(file("NO_FILE_${missing_id}")) }
    ENRICHMENT(
        absent(), absent(), absent(), absent(), absent(),
        Channel.value(file(params.fixture_background)),
        absent(), absent(), absent(), absent(), absent(),
        Channel.value(file(params.fixture_scores)),
        absent(), absent(), absent(), absent(), absent(),
        absent(), absent(), absent(), absent(), absent(), absent(), absent(),
        params.fixture_char ? Channel.value(file(params.fixture_char)) : Channel.empty()
    )
}
"""


def stage_fixture(tmp_path):
    files, maps = make_fixture(tmp_path)
    fixture = tmp_path / 'fixture'
    fixture.mkdir()
    shutil.copy(files['scores'], fixture / 'position_scores.tsv')
    shutil.copy(files['postproc'], fixture / 'filtered_discovery.tsv')
    shutil.copy(files['null_universe'], fixture / 'null_universe_summary.tsv')
    shutil.copy(files['background'], fixture / 'cleaned_background.txt')
    shutil.copy(files['custom'], fixture / 'custom.tsv')
    shutil.copytree(maps, fixture / 'maps')
    stage_null(fixture)
    return files, maps, fixture


N_CYCLES = 40


def stage_null(fixture):
    """Tested positions of G1 and a permulation null of two blocks of exchangeable cycles (caas_permulation layout)."""
    (fixture / 'background.output').write_text('G1\t0,1,2,3\n')
    null = fixture / 'caas_permulation'
    rng = np.random.default_rng(5)
    rows = [(f'b_{c}', p, float(rng.uniform(0.05, 1))) for c in range(1, N_CYCLES + 1) for p in range(4) if rng.random() < 0.5]
    frame = pd.DataFrame(rows, columns=['cycle', 'position', 'all']).astype({'position': np.int32})
    frame.insert(0, 'gene', 'G1')
    nd = frame.groupby('position').size().rename('nd').reset_index().assign(gene='G1')
    null.mkdir()
    nd.to_parquet(null / 'null_positions.parquet', index=False)
    for k, half in enumerate((range(1, N_CYCLES // 2 + 1), range(N_CYCLES // 2 + 1, N_CYCLES + 1))):
        block = null / 'blocks' / f'blk_{k}'
        block.mkdir(parents=True)
        names = [f'b_{c}' for c in half]
        frame[frame['cycle'].isin(names)].to_parquet(block / 'positions.parquet', index=False)
        (block / 'cycles.txt').write_text('\n'.join(names) + '\n')


def run_nextflow(main_name, script, tmp_path, args):
    main = ROOT / main_name
    main.write_text(script)
    try:
        cmd = [NEXTFLOW, 'run', str(main), '--outdir', str(tmp_path / 'out'), '-w', str(tmp_path / 'work'),
               '-ansi-log', 'false', *args]
        return subprocess.run(cmd, capture_output=True, text=True, cwd=str(ROOT), timeout=280)
    finally:
        main.unlink(missing_ok=True)


@pytest.mark.skipif(NEXTFLOW is None, reason='nextflow no disponible en PATH')
@pytest.mark.parametrize('null_source,level1', [('live', 'true'), ('prior', 'false')])
def test_posenrich_nextflow_smoke(tmp_path, null_source, level1):
    files, maps, fixture = stage_fixture(tmp_path)
    extra = ['--caas_null_dir', str(fixture / 'caas_permulation')] if null_source == 'prior' else []
    result = run_nextflow('.smoke_posenrich.nf', POSENRICH_MAIN, tmp_path, [
        '--fixture_dir', str(fixture), '--caas_map_dir', str(fixture / 'maps'), '--scoring', 'true',
        '--null_source', null_source, '--level1', level1, '--traitname', 'smoke', '--posenrich_min_genes', '1', *extra,
        '--posenrich_domains', 'false', '--posenrich_ucr', 'false',
        '--posenrich_custom_file', str(fixture / 'custom.tsv')])
    log = result.stdout + result.stderr
    assert result.returncode == 0, log[-4000:]
    assert 'POSITION_REGISTRY' in log and 'POSENRICH_WITHIN_GENE' in log
    pos = tmp_path / 'out' / 'posenrich'
    for rel in ('registry/positions.tsv', 'registry/position_evidence.tsv', 'registry/regions.tsv', 'registry/sites.tsv',
                'registry/collections.tsv', 'registry/coverage.tsv', 'registry/pending_mappings.tsv',
                'registry/manifest.json',
                'characterization/position_characterization.tsv', 'characterization/collection_summary.tsv',
                'characterization/ranking_summary.tsv', 'characterization/universe_summary.tsv',
                'universe/universe.tsv.gz', 'universe/universe_manifest.json',
                'within_gene/within_gene_results.tsv', 'within_gene/within_gene_by_gene.tsv',
                'within_gene/within_gene_manifest.json',
                'proximity/proximity_protein_distances.tsv', 'proximity/proximity_protein_summary.tsv',
                'proximity/proximity_genomic_distances.tsv', 'proximity/proximity_genomic_summary.tsv',
                'proximity/proximity_protein_hits.tsv', 'proximity/proximity_protein_by_gene.tsv'):
        assert (pos / rel).exists(), rel
    # Report 13 is published in posenrich/ and html_reports/, and states level 1 when it did not run.
    for folder in (pos, tmp_path / 'out' / 'html_reports'):
        reports = list(folder.glob('13.POSENRICH_report_smoke.html'))
        assert reports, f'POSENRICH report not published in {folder}'
    text = (pos / '13.POSENRICH_report_smoke.html').read_text()
    assert 'Level 2: concentration within genes' in text
    assert ('Level 1: not available' in text) == (level1 == 'false')
    if level1 == 'true':
        # Level 1 reads the universe table of level 2: the same positions.
        acc = pd.read_csv(pos / 'accumulation/accumulation_results.tsv', sep='\t')
        universe = pd.read_csv(pos / 'universe/universe.tsv.gz', sep='\t')
        assert acc['n_positions'].sum() == len(universe)
        assert acc['positives'].sum() == universe['observed_positive'].sum()
        assert (pos / 'accumulation/accumulation_manifest.json').exists()
    # El motor usa los 40 ciclos de los dos bloques y la lista de genes del universo coincide con la del registro.
    manifest = json.loads((pos / 'within_gene/within_gene_manifest.json').read_text())
    assert manifest['counts']['n_null_cycles'] == N_CYCLES and manifest['counts']['n_blocks'] == 2
    assert manifest['gene_list_check']['status'] == 'match'
    assert manifest['parameters']['min_genes'] == 1
    # El registro no escribe ningún GMT.
    assert not (pos / 'registry/gmts').exists()


LEGACY_CASES = {
    'bad_checksum': 'checksum mismatch',
    'wrong_alignment': 'belongs to alignment',
    'plain': 'has no coordinate_space',
    'gzip': 'has no coordinate_space',
    'ucr_plain': 'has no coordinate_space',
    'ucr_gzip': 'has no coordinate_space',
}


@pytest.mark.skipif(NEXTFLOW is None, reason='nextflow no disponible en PATH')
@pytest.mark.parametrize('legacy', sorted(LEGACY_CASES))
def test_posenrich_rejects_unverified_references(tmp_path, legacy):
    files, maps, fixture = stage_fixture(tmp_path)
    legacy_path = tmp_path / ('legacy.tsv.gz' if legacy in ('gzip', 'ucr_gzip') else 'legacy.tsv')
    text = 'gene\tpfam_id\ttarget_name\tdescription\tclan_acc\tclan_name\tali_start\tali_end\n'
    if legacy.startswith('ucr_'):
        text = 'gene\tposition\tregion_type\tmethod\tvariability\n'
    args = ['--fixture_dir', str(fixture), '--caas_map_dir', str(fixture / 'maps'), '--scoring', 'true']
    if legacy in ('bad_checksum', 'wrong_alignment'):
        legacy_path.write_text(text.rstrip('\n') + '\tcoordinate_space\n')
        meta = {'table_sha256': '0' * 64} if legacy == 'bad_checksum' else {}
        if legacy == 'wrong_alignment':
            original = tmp_path / 'original_alignments'
            original.mkdir()
            actual = tmp_path / 'actual_alignments'
            actual.mkdir()
            meta['alignment_dir'] = str(original)
            args += ['--alignment', str(actual)]
        Path(str(legacy_path) + '.manifest.json').write_text(json.dumps(meta))
    elif legacy in ('gzip', 'ucr_gzip'):
        import gzip
        with gzip.open(legacy_path, 'wt') as handle:
            handle.write(text)
    else:
        legacy_path.write_text(text)
    if legacy.startswith('ucr_'):
        args += ['--posenrich_domains', 'false', '--posenrich_ucr', 'true', '--ucr_positions_file', str(legacy_path)]
    else:
        args += ['--posenrich_ucr', 'false', '--posenrich_domains', 'true', '--domain_variability_file', str(legacy_path)]
    result = run_nextflow('.smoke_posenrich_refs.nf', POSENRICH_MAIN, tmp_path, args)
    log = result.stdout + result.stderr
    assert result.returncode != 0
    assert LEGACY_CASES[legacy] in log
    assert 'Submitted process' not in log


@pytest.mark.skipif(NEXTFLOW is None, reason='nextflow no disponible en PATH')
@pytest.mark.parametrize('with_char', [True, False])
def test_enrichment_reads_optional_characterization(tmp_path, with_char):
    files, maps = make_fixture(tmp_path)
    args = ['--fcs_enabled', 'false', '--scoring_ami', 'false', '--scoring_string', 'false',
            '--fixture_background', files['background'], '--fixture_scores', files['scores'],
            '--fixture_char', '']
    if with_char:
        registry = run_registry(tmp_path, files, maps)
        char_dir = tmp_path / 'characterization'
        done = subprocess.run([sys.executable, str(CHARACTERIZE), '--registry_dir', str(registry),
                               '--output_dir', str(char_dir)], capture_output=True, text=True)
        assert done.returncode == 0, done.stderr
        args[-1] = str(char_dir / 'position_characterization.tsv')
    else:
        args = args[:-2]
    result = run_nextflow('.smoke_enrichment.nf', ENRICHMENT_MAIN, tmp_path, args)
    log = result.stdout + result.stderr
    assert result.returncode == 0, log[-5000:]
    assert list((tmp_path / 'out' / 'html_reports').glob('14.Comparison_report_*.html'))
    compare_tasks = [p.parent for p in (tmp_path / 'work').rglob('.command.sh') if '14.Comparison_report.Rmd' in p.read_text()]
    assert len(compare_tasks) == 1
    staged = compare_tasks[0] / 'position_characterization.tsv'
    if with_char:
        assert staged.is_file()
        assert staged.read_bytes() == (tmp_path / 'characterization/position_characterization.tsv').read_bytes()
    else:
        # El sentinel NO_POSITION_CHAR queda como enlace sin destino: el informe lo trata como ausente.
        assert staged.is_symlink() and not staged.exists()
        assert staged.resolve().name == 'NO_POSITION_CHAR'
