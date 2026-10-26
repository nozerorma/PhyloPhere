"""Calibration must emit one task input even when no live null blocks arrive."""
from pathlib import Path
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which('nextflow')


@pytest.mark.skipif(NEXTFLOW is None, reason='Nextflow required')
@pytest.mark.parametrize('mode', ['live', 'saved', 'absent', 'empty_live'])
def test_calibration_blocks_reach_evidence_task(tmp_path, mode):
    # Execute the production channel expression, including its empty-channel behavior.
    source = (ROOT / 'workflows/scoring.nf').read_text()
    expression = source.split('def calibration_blocks = ', 1)[1].split('\n        def evidence_out', 1)[0]
    prior = tmp_path / 'prior_null'
    prior.mkdir()
    live = tmp_path / 'blk_0001'
    live.mkdir()
    (live / 'marker').write_text('live')
    (prior / 'marker').write_text('saved')
    script = tmp_path / 'main.nf'
    script.write_text('''
nextflow.enable.dsl=2
params.mode = null
params.caas_null_dir = null
params.fixture = null
process EVIDENCE_INPUT {
    input:
    path blocks, stageAs: 'live/*'
    output:
    path 'received.txt'
    script:
    """
    cat live/*/marker > received.txt 2>/dev/null || touch received.txt
    """
}
workflow {
    def null_blocks_ch = params.mode == 'live'
        ? Channel.value([file("${params.fixture}/blk_0001")])
        : (params.mode == 'empty_live' ? Channel.empty() : null)
    def calibration_blocks = ''' + expression + '''
    EVIDENCE_INPUT(calibration_blocks)
    EVIDENCE_INPUT.out.view { "EVIDENCE_RESULT=" + it.text.trim() }
}
''')
    cmd = [NEXTFLOW, '-log', str(tmp_path / 'nextflow.log'), 'run', str(script),
           '-c', '/dev/null', '-w', str(tmp_path / 'work'), '-ansi-log', 'false',
           '--mode', mode, '--fixture', str(tmp_path)]
    if mode in ('live', 'saved'):
        cmd += ['--caas_null_dir', str(prior)]
    result = subprocess.run(cmd, cwd=tmp_path, text=True, capture_output=True, timeout=90)
    log = result.stdout + result.stderr
    assert result.returncode == 0, log
    assert 'EVIDENCE_INPUT' in log, log
    expected = {'live': 'live', 'saved': 'saved', 'absent': '', 'empty_live': ''}[mode]
    assert f'EVIDENCE_RESULT={expected}\n' in log, log
