"""Serial Nextflow integration: real rendering, audited coverage and content caching."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import pytest
ROOT=Path(__file__).resolve().parents[1]
pytestmark=pytest.mark.skipif(not shutil.which('nextflow'),reason='Run in phylophere')


@pytest.fixture
def project(tmp_path):
    for name in ['FADE','SELECTION']:
        shutil.copytree(ROOT/'subworkflows'/name,tmp_path/'subworkflows'/name,ignore=shutil.ignore_patterns('__pycache__'))
    style=tmp_path/'subworkflows/TRAIT_ANALYSIS/local/src';style.mkdir(parents=True)
    shutil.copy2(ROOT/'subworkflows/TRAIT_ANALYSIS/local/src/report_style.R',style/'report_style.R')
    (tmp_path/'workflows').mkdir();shutil.copy2(ROOT/'workflows/fade.nf',tmp_path/'workflows/fade.nf')
    (tmp_path/'alignments').mkdir()
    (tmp_path/'alignments/001.fa').write_text('>A\nACDEFGHIKLMNPQRSTVWY\n>X.A\nACDEYGHIKLMNPQRSTVWY\n>B\nVCDEFGHIKLMNPQRSTAWY\n>C\nACDEFGHIKLMNPQKSTVWY\n')
    (tmp_path/'alignments/NA.fa').write_text('>X.A\nACDEFGHIKLMNPQRSTVWY\n>C\nACDEYGHIKLMNPQRSTVWY\n')
    (tmp_path/'tree.nwk').write_text('((A:0.1,X.A:0.1):0.1,(B:0.1,C:0.1):0.1);')
    (tmp_path/'top.txt').write_text('A\n');(tmp_path/'bottom.txt').write_text('B\n')
    (tmp_path/'main.nf').write_text('''include { FADE } from './workflows/fade'
workflow {
    def inputs=Channel.fromPath('alignments/*.fa').map { f -> tuple(f.baseName,f) }
    def fan=inputs.multiMap { g,f -> top:tuple(g,f);bottom:tuple(g,f) }
    FADE(fan.top,fan.bottom,Channel.value(file('tree.nwk')),Channel.value(file('top.txt')),Channel.value(file('bottom.txt')))
}
''')
    (tmp_path/'audit.config').write_text('''process.executor='local'
process.cpus=1
trace.overwrite=true
params {
    outdir='out'; fade_direction='both'; fade_background_scope='all'; fade_internal_nodes='all_descendants'
    fade_model='LG'; fade_method='Variational-Bayes'; fade_grid=5; fade_concentration=0.5; fade_bf_threshold=125
    fade_batch_size=2; use_singularity=false;use_apptainer=false; traitname="cancer's trait"
}
''')
    (tmp_path/'bin').mkdir()
    fake=tmp_path/'bin/hyphy'
    fake.write_text('''#!/usr/bin/env python
import json,os,sys
args=dict(zip(sys.argv[2::2],sys.argv[3::2]))
if os.environ.get('FADE_TEST_FAIL'):sys.exit(42)
n=len(next(line.strip() for line in open(args['--alignment']) if not line.startswith('>')))
aas='ACDEFGHIKLMNPQRSTVWY'
j={'runtime':'test-double','input':{'number of sites':n},'MLE':{'headers':[['rate',''],['bias',''],['P',''],['BayesFactor[bias>0]','']],'content':{a:{'0':[[1,1,.1,1] for _ in range(n)]} for a in aas}}}
j['MLE']['content']['A']['0'][0][3]=float(os.environ.get('FADE_TEST_BF','150'))
if os.environ.get('FADE_TEST_PARTIAL'):del j['MLE']['content']['Y']
open(args['--output'],'w').write(json.dumps(j))
''');fake.chmod(0o755)
    return tmp_path


def run(project,*args,success=True,real=False,extra=None):
    env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1'
    if not real:env['PATH']=str(project/'bin')+os.pathsep+env['PATH']
    env.update(extra or {})
    result=subprocess.run(['nextflow','-log','review.log','run','main.nf','-c','audit.config','-ansi-log','false','-with-trace','trace.tsv',*args],cwd=project,env=env,text=True,capture_output=True,timeout=180)
    assert (result.returncode==0)==success,result.stdout+result.stderr
    return result.stdout+result.stderr


def test_actual_workflow(project):
    run(project)
    for direction in ['top','bottom']:
        out=project/'out/selection/fade'/direction
        metadata=json.loads((out/'fade_run_metadata.json').read_text())
        assert (metadata['expected'],metadata['tested'],metadata['excluded'])==(2,1,1)
        assert (out/'gene_lists/background.txt').read_text()=='001\n'
        assert (out/'fade_inputs/settings.json').is_file()
        html=(out/f'6.FADE_report_{direction}.html').read_text()
        assert "cancer" in html and 'LG' in html and '125' in html
    assert 'X.A{Foreground}' not in next((project/'work').glob('*/*/001.top_fg.nwk')).read_text()


@pytest.mark.parametrize('batch,scope,direction',[('1','all','both'),('2','opposite','both'),('2','all','top')])
def test_modes_and_direction_gating(project,batch,scope,direction):
    run(project,'--fade_batch_size',batch,'--fade_background_scope',scope,'--fade_direction',direction)
    for d in (['top','bottom'] if direction=='both' else ['top']):
        out=project/'out/selection/fade'/d
        assert (out/'gene_lists/background.txt').read_text()=='001\n'
        assert json.loads((out/'fade_run_metadata.json').read_text())['settings']['background_scope']==scope
    if direction=='top':assert not (project/'out/selection/fade/bottom').exists()


@pytest.mark.parametrize('empty',['no_inputs','excluded'])
def test_empty_results_have_complete_reports_and_lists(project,empty):
    (project/'alignments/001.fa').unlink()
    if empty=='no_inputs':(project/'alignments/NA.fa').unlink()
    run(project)
    for d in ['top','bottom']:
        out=project/'out/selection/fade'/d
        assert (out/'gene_lists/background.txt').read_text()==''
        assert (out/'gene_lists'/f'fade_{d}_significant.txt').read_text()==''
        meta=json.loads((out/'fade_run_metadata.json').read_text())
        assert meta['complete'] and meta['tested']==0
        assert meta['expected']==(0 if empty=='no_inputs' else 1)
        assert (out/f'6.FADE_report_{d}.html').stat().st_size>1000


@pytest.mark.parametrize('batch,flag',[('1','FADE_TEST_FAIL'),('2','FADE_TEST_FAIL'),('2','FADE_TEST_PARTIAL')])
def test_technical_failures_cannot_publish_complete_tables(project,batch,flag):
    log=run(project,'--fade_batch_size',batch,success=False,extra={flag:'1'})
    assert 'FADE technical failures' in log
    assert not list((project/'out').rglob('fade_summary_*.tsv'))
    audits=list((project/'work').glob('*/*/*fade.audit.json'))
    assert any(json.loads(p.read_text())['status']=='technical_failure' for p in audits)
    assert list((project/'work').glob('*/*/*.hyphy.log'))


@pytest.mark.parametrize('args,message',[
    (['--fade_model','Dayhoff'],'Unsupported FADE model'),
    (['--fade_batch_size','0'],'positive integer'),
    (['--fade_bf_threshold','0'],'Invalid FADE bf_threshold'),
    (['--fade_internal_nodes','typo'],'Invalid FADE internal-node mode'),
])
def test_invalid_settings_fail_early(project,args,message):
    assert message in run(project,*args,success=False)


def stat_edit(path,old,new):
    st=path.stat();body=path.read_text();assert old in body and len(old)==len(new)
    path.write_text(body.replace(old,new));os.utime(path,ns=(st.st_atime_ns,st.st_mtime_ns))


def test_resume_content_helpers_and_regeneration(project):
    run(project,'--fade_direction','top')
    run(project,'--fade_direction','top','-resume')
    assert 'CACHED' in (project/'trace.tsv').read_text()
    stat_edit(project/'alignments/001.fa','ACDEFGHIKLMNPQRSTVWY','VCDEFGHIKLMNPQRSTVWY')
    run(project,'--fade_direction','top','-resume')
    assert 'COMPLETED' in (project/'trace.tsv').read_text()
    audit=json.loads((project/'out/selection/fade/audit/001.top.fade.audit.json').read_text())
    import hashlib
    assert audit['input_sha256']['alignment']==hashlib.sha256((project/'alignments/001.fa').read_bytes()).hexdigest()
    stat_edit(project/'subworkflows/FADE/local/src/fade_json_to_tables.py',"'.17g'","'.16g'")
    run(project,'--fade_direction','top','-resume')
    assert 'COMPLETED' in (project/'trace.tsv').read_text()
    # A report-only code edit must preserve the inference cache.
    stat_edit(project/'subworkflows/FADE/local/6.FADE_report.Rmd','Target amino-acid frequency','Target amino-acid FREQUENCY')
    run(project,'--fade_direction','top','-resume')
    trace=(project/'trace.tsv').read_text()
    assert any('FADE_BATCHED' in line and 'CACHED' in line for line in trace.splitlines())
    assert any('FADE_RENDER' in line and 'COMPLETED' in line for line in trace.splitlines())
    sys_path=str(ROOT)
    import sys
    if sys_path not in sys.path:sys.path.insert(0,sys_path)
    from gui.generation.report_registry import detect_reports,build_script
    reports=detect_reports(project/'out',ROOT)
    report=next(r for r in reports if r.spec.id=='fade_top')
    script=project/'regenerate.sh';script.write_text(build_script(report,ROOT,False))
    subprocess.run(['bash','-n',str(script)],check=True)
    subprocess.run(['bash',str(script)],check=True,capture_output=True,text=True,timeout=90)
    html=(project/'out/html_reports/6.FADE_report_top.html').read_text()
    assert '125' in html and 'LG' in html


@pytest.mark.skipif(not shutil.which('hyphy'),reason='HyPhy unavailable')
def test_real_hyphy_and_shared_baseline(project):
    run(project,real=True)
    audit=json.loads((project/'out/selection/fade/audit/001.bottom.fade.audit.json').read_text())
    assert audit['status']=='complete' and audit['baseline_reused']
    assert audit['n_sites_tested']==20


def test_infinite_bf_is_preserved_in_tables_lists_and_report(project):
    run(project,'--fade_direction','top',extra={'FADE_TEST_BF':'inf'})
    out=project/'out/selection/fade/top'
    assert 'Inf' in (out/'fade_summary_top.tsv').read_text()
    assert (out/'gene_lists/fade_top_significant.txt').read_text()=='001\n'
    assert 'infinite' in (out/'6.FADE_report_top.html').read_text()


def test_legacy_regeneration_does_not_invent_original_settings(project):
    run(project,'--fade_direction','top')
    out=project/'out/selection/fade/top'
    shutil.rmtree(out/'fade_inputs');(out/'fade_run_metadata.json').unlink()
    import sys
    if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
    from gui.generation.report_registry import detect_reports,build_script
    report=next(r for r in detect_reports(project/'out',ROOT) if r.id=='fade_top')
    script=project/'legacy.sh';script.write_text(build_script(report,ROOT,False))
    subprocess.run(['bash',str(script)],check=True,capture_output=True,text=True,timeout=90)
    html=(project/'out/html_reports/6.FADE_report_top.html').read_text()
    assert 'original cutoff (not recorded)' in html
    assert 'Original settings and completion audit were not recorded' in html


@pytest.mark.parametrize('batch',['1','2'])
def test_dotted_and_quoted_gene_names(project,batch):
    gene="a.filtered.fa.v2 'q'"
    (project/'alignments/001.fa').rename(project/'alignments'/f'{gene}.fa')
    run(project,'--fade_batch_size',batch,'--fade_direction','top')
    assert (project/'out/selection/fade/top/gene_lists/background.txt').read_text()==gene+'\n'


@pytest.mark.parametrize('body',['',' A\nB\n','A\nA\n'])
def test_invalid_or_overlapping_species_pools_fail(project,body):
    (project/'top.txt').write_text(body)
    run(project,success=False)


def test_report_rejects_inputs_changed_since_completion(project):
    run(project,'--fade_direction','top')
    saved=project/'out/selection/fade/top/fade_inputs/gene_counts.tsv'
    saved.write_text(saved.read_text().replace('150','151'))
    import sys
    if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
    from gui.generation.report_registry import detect_reports,build_script
    report=next(r for r in detect_reports(project/'out',ROOT) if r.id=='fade_top')
    script=project/'bad-report.sh';script.write_text(build_script(report,ROOT,False))
    p=subprocess.run(['bash',str(script)],capture_output=True,text=True,timeout=90)
    assert p.returncode!=0


def test_repository_resources_bind_to_render_and_gene_list_tasks(project):
    cfg=project/'audit.config'
    cfg.write_text(f"includeConfig '{ROOT}/conf/resources.config'\n"+cfg.read_text()+
                   "\nprocess.resourceLimits=[cpus:1,memory:2.GB,time:1.h]\n")
    run(project,'--fade_direction','top')
    render=next(p for p in (project/'work').glob('*/*/.command.sh') if 'render_fade.R' in p.read_text())
    assert 'VROOM_THREADS=1' in (render.parent/'.command.run').read_text()
    assert (project/'out/selection/fade/top/gene_lists/background.txt').read_text()=='001\n'


def test_all_mode_reuses_prepared_alignment(project):
    run(project,'--fade_direction','top')
    audit=next((project/'work').glob('*/*/001.top.fade.audit.json'))
    assert not (audit.parent/'001.top.fa').exists()
    assert (audit.parent/'fastas/001.fa').is_file()
