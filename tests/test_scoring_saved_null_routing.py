"""Exercise orchestration and SCORING channels together, with cheap task bodies."""
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which('nextflow')


def _mock_process(name, n_inputs, outputs, body=None):
    return ('process ' + name + ' {\ninput:\n'
            + ''.join(f'val x{i}\n' for i in range(n_inputs))
            + 'output:\n'
            + ''.join(f"path '{o}', emit: {o}\n" for o in outputs)
            + 'script:\n"""\n'
            + (body or 'touch ' + ' '.join(outputs)) + '\n"""\n}\n')


@pytest.mark.skipif(NEXTFLOW is None, reason='Nextflow required')
@pytest.mark.parametrize('mode', ['saved', 'live_discovery', 'live_disambiguation'])
def test_saved_null_reaches_caas_seed_preparation(tmp_path, mode):
    # Keep the production workflow and main's replay decision. Replace only the
    # expensive task bodies; the final seed preparation runs the real Python CLI.
    source = (ROOT / 'workflows/scoring.nf').read_text()
    source = re.sub(r'^include .*$', '', source, flags=re.M)
    source = source.replace('#!/usr/bin/env nextflow', '')
    evidence_outputs = ['position_scores', 'gene_scores', 'fcs_stats', 'gene_evidence',
                        'metadata', 'gene_lists', 'position_lists', 'calibration']
    body = '''mkdir -p gene_lists
    for rank in global top bottom; do
        printf 'Gene\\tranking\\tnominal\\n' > gene_lists/active_\\${rank}.tsv
    done
    echo "${x4}" > null_source
    touch ''' + ' '.join(o for o in evidence_outputs if o != 'gene_lists')
    mocks = _mock_process('SCORING_NULL_POSITIONS', 7, ['counters', 'fit', 'base'])
    mocks += _mock_process('SCORING_COMPUTE', 11, ['gene_scores', 'position_scores', 'fcs_stats',
                          'fcs_stats_rer', 'fcs_stats_fade',
                          'null_universe_summary', 'gene_correlations'])
    mocks += _mock_process('SCORING_EVIDENCE', 6, evidence_outputs + ['null_source'],
                          body + '\n')
    # null_source is a diagnostic mock output, not part of the production contract.
    mocks += _mock_process('SCORING_REPORT', 17, ['report'])
    (tmp_path / 'scoring.nf').write_text(mocks + source)
    main_source = (ROOT / 'main.nf').read_text()
    decision = re.search(r'def run_caas_permulation\s*= (.+)', main_source).group(1)
    fixture = tmp_path / 'prior'
    (fixture / 'blocks/blk_1').mkdir(parents=True)
    for name in ['null_positions.parquet', 'counts.json', 'cycle_roster.txt']:
        (fixture / name).touch()
    for name in ['postproc.tsv', 'background.txt']:
        (tmp_path / name).touch()
    script = '''nextflow.enable.dsl=2
include { SCORING } from './scoring.nf'
def toBool(v) { v == true || v == 'true' }
process PREPARE_SEEDS {
    input:
    path seeds
    output:
    path 'prepared'
    script:
    """
    ''' + sys.executable + ' ' + str(ROOT / 'subworkflows/ENRICHMENT/local/src/prepare_ami_lists.py') + ''' --method caas --sources ${seeds} --output-dir prepared
    """
}
workflow {
    def ran_discovery = params.mode == 'live_discovery'
    def run_ct_disambiguation = params.mode == 'live_disambiguation'
    def run_caas_permulation = ''' + decision + '''
    println "REPLAY=${run_caas_permulation}"
    // The failed recovery had no caas_config: implicit replay outputs were empty.
    def live_positions = !run_caas_permulation ? null :
        (ran_discovery || run_ct_disambiguation ? Channel.value(file("${params.caas_null_dir}/null_positions.parquet")) : Channel.empty())
    def live_blocks = !run_caas_permulation ? null :
        (ran_discovery || run_ct_disambiguation ? Channel.value([file("${params.caas_null_dir}/blocks/blk_1")]) : Channel.empty())
    SCORING(Channel.value(file(params.scoring_postproc_input)), null, null, null, null, null, null,
        Channel.value(file(params.scoring_background_input)), null, null,
        live_positions, live_blocks, null, null, null, null, null, null, null, null, null)
    PREPARE_SEEDS(SCORING.out.gene_lists)
    PREPARE_SEEDS.out.view { "SEEDS_READY=" + it.resolve('canonical_lists.json').text }
}
'''
    (tmp_path / 'main.nf').write_text(script)
    params = dict(mode=mode, traitname='regression', outdir=str(tmp_path / 'out'),
                  enrichment=True, caas_permulation_enrichment=True, caas_null_dir=str(fixture),
                  scoring_postproc_input=str(tmp_path / 'postproc.tsv'),
                  scoring_background_input=str(tmp_path / 'background.txt'))
    (tmp_path / 'params.json').write_text(json.dumps(params))
    result = subprocess.run([NEXTFLOW, '-log', str(tmp_path / 'run.log'), 'run',
                             str(tmp_path / 'main.nf'), '-c', '/dev/null',
                             '-params-file', str(tmp_path / 'params.json'),
                             '-w', str(tmp_path / 'work'), '-ansi-log', 'false'],
                            cwd=tmp_path, capture_output=True, text=True, timeout=90)
    log = result.stdout + result.stderr
    assert result.returncode == 0, log
    assert f"REPLAY={'false' if mode == 'saved' else 'true'}" in log, log
    assert 'SEEDS_READY=' in log, log
    assert all(f'"{rank}"' in log for rank in ('global', 'top', 'bottom')), log
