"""Retry policy of conf/resources.config, run through the real configuration.

A process that scales its resources with task.attempt must be retried; a process that declares
errorStrategy 'terminate' in its own definition must not be.
"""
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which('nextflow')

# exit 137 on the first attempt; the retry prints the escalated request
RETRIED = '''process TREE_CLEANUP {
    label 'process_low'
    input: val x
    output: stdout
    script:
    """
    echo "attempt=${task.attempt} mem=${task.memory}"
    if [ ${task.attempt} -eq 1 ]; then exit 137; fi
    """
}
workflow { TREE_CLEANUP(1).view() }
'''

# keeps the 'terminate' of its own definition: the first failure ends the run
TERMINATED = '''process RER_TRAIT {
    label 'process_low'
    errorStrategy 'terminate'
    input: val x
    output: stdout
    script:
    """
    exit 137
    """
}
workflow { RER_TRAIT(1).view() }
'''


def run(tmp_path, source):
    (tmp_path / 'main.nf').write_text(source)
    return subprocess.run(
        [NEXTFLOW, '-log', str(tmp_path / 'run.log'), 'run', 'main.nf', '-c', str(ROOT / 'nextflow.config'),
         '-profile', 'local', '-w', str(tmp_path / 'work'), '-ansi-log', 'false'],
        cwd=tmp_path, capture_output=True, text=True, timeout=120)


@pytest.mark.skipif(NEXTFLOW is None, reason='nextflow not on PATH')
def test_process_that_scales_with_attempt_is_retried(tmp_path):
    r = run(tmp_path, RETRIED)
    log = r.stdout + r.stderr
    assert r.returncode == 0, log
    assert 'Execution is retried (1)' in log
    assert 'attempt=2 mem=1000 MB' in log


@pytest.mark.skipif(NEXTFLOW is None, reason='nextflow not on PATH')
def test_process_with_own_terminate_is_not_retried(tmp_path):
    r = run(tmp_path, TERMINATED)
    log = r.stdout + r.stderr
    assert r.returncode != 0, log
    assert 'is retried' not in log
