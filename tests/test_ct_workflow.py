"""Actual local Nextflow checks of CT dispatch, staging, defaults and reuse."""
from pathlib import Path
import json
import shutil
import subprocess

import pytest

ROOT=Path(__file__).resolve().parents[1]
pytestmark=pytest.mark.skipif(not shutil.which('nextflow') or not shutil.which('Rscript'),reason='Run in phylophere')


@pytest.fixture
def project(tmp_path):
    (tmp_path/'subworkflows').symlink_to(ROOT/'subworkflows',target_is_directory=True)
    (tmp_path/'curated_tree.nwk').write_text('(((A:1,B:1):1,(C:1,D:1):1):1,(E:1,F:1):2);\n')
    (tmp_path/'cfg.tab').write_text('B\t1\t1\nA\t0\t1\nD\t1\t2\nC\t0\t2\nF\t1\t3\nE\t0\t3\n')
    (tmp_path/'values.csv').write_text('species,sample_size,score rate\nA,1,.1\nB,2,.9\nC,3,.2\nD,4,.8\nE,5,.3\nF,6,.7\n')
    (tmp_path/'main.nf').write_text(
        f"include {{ CT }} from '{ROOT}/workflows/ct.nf'\n"
        "workflow { CT(null,null,Channel.value(file('curated_tree.nwk')),null) }\n")
    (tmp_path/'audit.config').write_text(
        f"includeConfig '{ROOT}/conf/common.config'\nincludeConfig '{ROOT}/conf/ct.config'\n"
        "process.executor='local'\nprocess.cpus=2\nparams { ct_tool='resample'; enrichment=false; alignment=null; tree='MISSING_RAW_TREE.nwk'; caas_config='cfg.tab'; my_traits='values.csv'; traitname='score rate'; outdir='out'; use_singularity=false; use_apptainer=false; caas_full_perms=3; max_tries=20; pss_top_pct=0; perm_match_pss=false; max_fop=3; publish_intermediates=true }\n")
    return tmp_path


def run(project,*args):
    return subprocess.run(['nextflow','-log','review.log','run','main.nf','-c','audit.config','-ansi-log','false',*args],cwd=project,text=True,capture_output=True,timeout=120)


def test_resample_only_needs_no_alignments_uses_curated_tree_and_actual_auto_default(project):
    r=run(project)
    assert r.returncode==0,r.stdout+r.stderr
    out=next((project/'out/resample').glob('*.resampled.output'))
    provenance=json.loads((out/'permulation_provenance.json').read_text())
    assert provenance['strategy']=='auto' and provenance['selected_cycles']==3
    assert (project/'out/caastools/resample.tab').stat().st_size>0
    scripts=list((project/'work').glob('*/*/.command.sh'))
    resample=next(p for p in scripts if 'permulations.R' in p.read_text())
    assert 'traitvalues/values.csv' in resample.read_text()
    assert (resample.parent/'nw_tree.nwk').resolve()==(project/'curated_tree.nwk').resolve()
    assert not any('prepare_labelings.py' in p.read_text() for p in scripts)


def test_zero_null_needs_no_phenotypes_and_publishes_empty_concat(project):
    r=run(project,'--caas_full_perms','0','--my_traits','')
    assert r.returncode==0,r.stdout+r.stderr
    assert (project/'out/caastools/resample.tab').read_text()==''
    provenance=json.loads(next((project/'out/resample').glob('*/permulation_provenance.json')).read_text())
    assert provenance['selected_cycles']==0


@pytest.mark.parametrize('tool',['resampl','true'])
def test_unknown_or_boolean_tool_fails(project,tool):
    r=run(project,'--ct_tool',tool)
    assert r.returncode!=0 and ('ct_tool' in r.stdout+r.stderr)


def test_short_reused_null_fails_in_actual_subset_process(project):
    (project/'resample').mkdir();(project/'resample/resample_001.tab').write_text('b_1\tB,D,F\tA,C,E\n')
    (project/'main.nf').write_text(
        f"include {{ SUBSET_RESAMPLE_PERMS }} from '{ROOT}/subworkflows/CT/caas_permulation.nf'\n"
        "workflow { SUBSET_RESAMPLE_PERMS(Channel.value(file('resample')),3,2,Channel.value(file('cfg.tab'))) }\n")
    r=run(project)
    assert r.returncode!=0 and 'only 1 available' in r.stdout+r.stderr


def test_normalized_upstream_traits_use_species_even_with_custom_raw_header(project):
    (project/'main.nf').write_text(
        f"include {{ CT }} from '{ROOT}/workflows/ct.nf'\n"
        "workflow { CT(Channel.value(file('cfg.tab')),Channel.value(file('values.csv')),Channel.value(file('curated_tree.nwk')),null) }\n")
    r=run(project,'--sp_colname','sample id')
    assert r.returncode==0,r.stdout+r.stderr


def test_custom_standalone_species_header_is_preserved(project):
    p=project/'values.csv';p.write_text(p.read_text().replace('species,','sample id,'))
    r=run(project,'--sp_colname','sample id')
    assert r.returncode==0,r.stdout+r.stderr


def test_saved_null_does_not_force_alignment_replay(project):
    r=run(project,'--caas_full_perms','0','--enrichment','true','--caas_null_dir','external_saved_null')
    assert r.returncode==0,r.stdout+r.stderr
    assert not any('prepare_labelings.py' in p.read_text() for p in (project/'work').glob('*/*/.command.sh'))


def test_external_observed_scoring_requires_replay_inputs(project):
    r=run(project,'--ct_disambiguation','true','--discovery_from','external_discovery.tab')
    assert r.returncode!=0 and 'nonempty alignment directory' in r.stdout+r.stderr
