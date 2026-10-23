"""Runs the WP9 scientific-control oracles without changing productive statistics."""
import os
from pathlib import Path
import shutil
import subprocess
import pytest

ROOT=Path(__file__).resolve().parents[1]

def test_wp9_independent_oracles(tmp_path):
    r=shutil.which('Rscript')
    if not r:pytest.skip('Rscript required')
    result=subprocess.run([r,'tests/test_wp9.R'],cwd=ROOT,env=dict(os.environ,TMPDIR=str(tmp_path)),capture_output=True,text=True,timeout=120)
    assert result.returncode==0,result.stdout+result.stderr
    assert 'PASS WP9' in result.stdout


def test_wp9_existing_launcher_uses_new_snapshot_and_preserves_bounds(tmp_path):
    repo=tmp_path/'checkout'
    (repo/'adhoc_scripts').mkdir(parents=True)
    (repo/'adhoc_scripts/wp9_gene_stat_control.R').write_text('# frozen validation entry\n')
    subprocess.run(['git','init','-q',str(repo)],check=True)
    subprocess.run(['git','add','.'],cwd=repo,check=True)
    subprocess.run(['git','-c','user.name=WP9 Test','-c','user.email=wp9@example.invalid','commit','-qm','fixture'],cwd=repo,check=True)
    capture=tmp_path/'args.txt'
    fake=tmp_path/'fake_rscript'
    fake.write_text('#!/bin/bash\nprintf "%s\\n" "$@" > "$WP9_CAPTURE"\n')
    fake.chmod(0o755)
    out=tmp_path/'new_run'
    env=dict(os.environ,WP9_CONTROL='true',TRAIT='neoplasia',MODE='strict',WP9_OUT=str(out),
             PHYLOPHERE_ROOT=str(repo),RSCRIPT=str(fake),WP9_CAPTURE=str(capture),RANKINGS='global',
             STORE=str(tmp_path/'saved_store'),UNIVERSE=str(tmp_path/'universe'),MAX_G='0',EXPECTED_CYCLES='10000')
    launcher=ROOT/'adhoc_scripts/fcs_calibrated_ranking_control.sbatch'
    run=subprocess.run(['bash',str(launcher)],env=env,capture_output=True,text=True)
    assert run.returncode==0,run.stdout+run.stderr
    args=capture.read_text().splitlines()
    assert Path(args[0]).read_text()=='# frozen validation entry\n'
    assert args[args.index('--max-genes')+1]=='0'
    assert args[args.index('--expected-cycles')+1]=='10000'
    assert args[args.index('--variants')+1]=='max,fmax,beta025,beta05,beta075,size_adj_max'
    assert (Path(str(out)+'_code')/'commit.txt').exists()
    repeat=subprocess.run(['bash',str(launcher)],env=env,capture_output=True,text=True)
    assert repeat.returncode==2 and 'already exists' in repeat.stderr
    env['SLURM_ARRAY_TASK_ID']='2'
    array=subprocess.run(['bash',str(launcher)],env=env,capture_output=True,text=True)
    assert array.returncode==2 and 'do not use an array' in array.stderr
