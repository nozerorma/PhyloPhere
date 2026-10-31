"""Real serial SCORING workflow, canonical null, and cache content regressions."""
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from test_scoring import ROOT, SRC, observed, null


def command(project,*args,ok=True):
    r=subprocess.run(['nextflow','-log','run.log','run','main.nf','-params-file','params.json','-c','audit.config',
                      '-ansi-log','false','-with-trace','trace.tsv',*args],cwd=project,text=True,capture_output=True,timeout=150)
    if ok:assert r.returncode==0,r.stdout+r.stderr+(project/'run.log').read_text()[-5000:]
    else:assert r.returncode!=0
    return r


def trace(project):
    with (project/'trace.tsv').open() as f:return {x['name'].split(':')[-1].split(' (')[0]:x['status'] for x in csv.DictReader(f,delimiter='\t')}


@pytest.fixture
def project(tmp_path):
    for name in ['SCORING','ENRICHMENT','CT_DISAMBIGUATION']:
        shutil.copytree(ROOT/'subworkflows'/name,tmp_path/'subworkflows'/name,ignore=shutil.ignore_patterns('__pycache__'))
    style=tmp_path/'subworkflows/TRAIT_ANALYSIS/local/src';style.mkdir(parents=True)
    shutil.copy2(ROOT/'subworkflows/TRAIT_ANALYSIS/local/src/report_style.R',style/'report_style.R')
    (tmp_path/'workflows').mkdir();shutil.copy2(ROOT/'workflows/scoring.nf',tmp_path/'workflows/scoring.nf')
    genes=observed(tmp_path/'observed.tsv')
    (tmp_path/'genes.txt').write_text('\n'.join(genes)+'\n')
    (tmp_path/'tested').write_text(''.join(f'{g}\t0,1,2\n' for g in genes))
    prov=dict(version=1,scientific=dict(model='LG',tree_sha256='x'*64,posterior_threshold=.8,taxid_sha256=None,weights_sha256=None),
              cluster=dict(enabled=True,remove=True,minlen=3,maxcaas=.7,min_positions=None),gene_filter=dict(mode='none'),strict=False,maps={},lengths={})
    (tmp_path/'observed.tsv.scoring.json').write_text(json.dumps(dict(prov,observed_sha256=null.file_digest(tmp_path/'observed.tsv'))))
    prior=tmp_path/'prior';block=prior/'blocks/blk_000001';block.mkdir(parents=True)
    cycles=[f'b_{i}' for i in range(1,20)]
    (block/'cycles.txt').write_text('\n'.join(cycles)+'\n');(prior/'cycle_roster.txt').write_text('\n'.join(cycles)+'\n')
    pq.write_table(pa.table(dict(gene=[genes[0]],cycle=['b_1'],position=[0],all=[.2])),block/'positions.parquet')
    pq.write_table(pa.table(dict(gene=[genes[0]],position=[0],nd=[1],first_rank=[0])),prior/'null_positions.parquet')
    counts=dict(cycles=19,cycles_present=['b_1'],scored_positions=1,units=1,removed_units=0,score_rule='caas_pos_species_pss_v4',scoring_provenance=dict(prov,universe=__import__('posenrich_fixtures').universe_provenance(tmp_path/'tested')))
    (block/'counts.json').write_text(json.dumps(counts));(prior/'counts.json').write_text(json.dumps(dict(counts,cycles_present=1)))
    code="g<-readLines('genes.txt');cy<-readLines('prior/cycle_roster.txt');m<-matrix(0,length(g),length(cy),dimnames=list(g,cy));m[1,1]<-.2;saveRDS(list(genes=g,cycles=cy,gene_stat='max',n_units=1,stat=list(global=m,top=m,bottom=m*0)),'prior/blocks/blk_000001/matrices.rds')"
    r=subprocess.run(['Rscript','-e',code],cwd=tmp_path,text=True,capture_output=True);assert r.returncode==0,r.stderr
    import hashlib
    (block/'scoring_artifacts.json').write_text(json.dumps(dict(version=1,files={n:hashlib.md5((block/n).read_bytes()).hexdigest() for n in ['matrices.rds','positions.parquet','cycles.txt','counts.json']})))
    (prior/'cycle_summary.tsv').write_text('cycle\tdirection\tn\tn_genes\n'+''.join(f'{c}\t{d}\t{int(c=="b_1" and d!="bottom")}\t{int(c=="b_1" and d!="bottom")}\n' for c in cycles for d in ['global','top','bottom']))
    (prior/'cycle_hist.tsv').write_text('cycle\tdirection\tbin\tn\nb_1\tglobal\t20\t1\nb_1\ttop\t20\t1\n')
    params=dict(traitname='audit',outdir='out',gene_stat='max',caas_null_dir='prior',posenrich_background_file='tested',scoring_background_input='genes.txt',
                use_singularity=False,use_apptainer=False,evidence_p_threshold=.1,evidence_fdr_threshold=.2,caas_permulation_enrichment=False)
    (tmp_path/'params.json').write_text(json.dumps(params))
    (tmp_path/'audit.config').write_text("process.executor='local'\nprocess.cpus=1\nexecutor.queueSize=1\ntrace.overwrite=true\nprocess.beforeScript='export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VROOM_THREADS=1'\n")
    (tmp_path/'main.nf').write_text('''include { SCORING } from './workflows/scoring'
workflow {
 SCORING(Channel.value(file('observed.tsv')),null,null,null,null,null,null,Channel.value(file('genes.txt')),
         null,null,null,null,null,null,null,null,null,null,null,Channel.value(file('tested')),null)
}
''')
    return tmp_path


def test_canonical_null_evidence_and_resume(project):
    command(project)
    out=project/'out/scoring'
    with (out/'gene_evidence.tsv').open() as fh:rows=list(csv.DictReader(fh,delimiter='\t',quoting=csv.QUOTE_NONE))
    assert len(rows)==15
    assert all(float(x['p.gen'])==.05 for x in rows if x['ranking']=='global')
    assert (out/'10.Scoring_report_audit.html').exists()
    command(project,'-resume');assert all(s=='CACHED' for s in trace(project).values())
    # Content-only change at identical size and mtime must invalidate compute.
    p=project/'observed.tsv';st=p.stat();p.write_text(p.read_text().replace('0.8','0.9'))
    os.utime(p,ns=(st.st_atime_ns,st.st_mtime_ns))
    contract=project/'observed.tsv.scoring.json';m=json.loads(contract.read_text());m['observed_sha256']=null.file_digest(p);contract.write_text(json.dumps(m))
    command(project,'-resume');assert trace(project)['SCORING_COMPUTE']!='CACHED'
    helper=project/'subworkflows/SCORING/local/src/scoring_io.R';helper.write_text(helper.read_text()+'\n# cache regression\n')
    command(project,'-resume');assert trace(project)['SCORING_COMPUTE']!='CACHED'
    report=project/'subworkflows/SCORING/local/10.Scoring_report.Rmd';report.write_text(report.read_text()+'\n<!-- report cache regression -->\n')
    command(project,'-resume');t=trace(project);assert t['SCORING_COMPUTE']=='CACHED' and t['SCORING_REPORT']!='CACHED'


def test_wrong_observed_null_settings_fail(project):
    p=project/'observed.tsv.scoring.json';m=json.loads(p.read_text());m['scientific']['model']='WAG';p.write_text(json.dumps(m))
    r=command(project,ok=False)
    assert 'scientific mismatch' in r.stdout+r.stderr+(project/'run.log').read_text()

@pytest.mark.parametrize('mode',['absent','silent','empty_observed'])
def test_missing_and_valid_empty_evidence(project,mode):
    if mode=='absent':
        p=project/'params.json';cfg=json.loads(p.read_text());cfg['caas_null_dir']='';p.write_text(json.dumps(cfg))
    elif mode=='silent':
        prior=project/'prior';block=prior/'blocks/blk_000001'
        pq.write_table(pa.table(dict(gene=pa.array([],pa.string()),cycle=pa.array([],pa.string()),position=pa.array([],pa.int32()),all=pa.array([],pa.float64()))),block/'positions.parquet')
        pq.write_table(pa.table(dict(gene=pa.array([],pa.string()),position=pa.array([],pa.int32()),nd=pa.array([],pa.int32()),first_rank=pa.array([],pa.int32()))),prior/'null_positions.parquet')
        for p in [block/'counts.json',prior/'counts.json']:
            cfg=json.loads(p.read_text());cfg['scored_positions']=0;cfg['cycles_present']=[] if p.parent==block else 0;p.write_text(json.dumps(cfg))
        r=subprocess.run(['Rscript','-e',"p<-'prior/blocks/blk_000001/matrices.rds';b<-readRDS(p);b$stat<-lapply(b$stat,function(m)m*0);saveRDS(b,p)"],cwd=project,text=True,capture_output=True);assert r.returncode==0,r.stderr
        summary=prior/'cycle_summary.tsv';summary.write_text(summary.read_text().replace('\t1\t1','\t0\t0'))
        (prior/'cycle_hist.tsv').write_text('cycle\tdirection\tbin\tn\n')
        import hashlib
        (block/'scoring_artifacts.json').write_text(json.dumps(dict(version=1,files={n:hashlib.md5((block/n).read_bytes()).hexdigest() for n in ['matrices.rds','positions.parquet','cycles.txt','counts.json']})))
    else:
        p=project/'observed.tsv';gs=(project/'genes.txt').read_text().splitlines()
        with p.open('w') as fh:
            w=csv.writer(fh,delimiter='\t',lineterminator='\n');w.writerow(['Gene','Position','side','asr_path_score'])
            for g in gs:w.writerow([g,0,'none',0])
        p=project/'observed.tsv.scoring.json';cfg=json.loads(p.read_text());cfg['observed_sha256']=null.file_digest(project/'observed.tsv');p.write_text(json.dumps(cfg))
    command(project)
    out=project/'out/scoring'
    with (out/'gene_evidence.tsv').open() as fh:rows=list(csv.DictReader(fh,delimiter='\t',quoting=csv.QUOTE_NONE))
    if mode=='absent':assert all(x['p.gen']=='NA' for x in rows)
    elif mode=='empty_observed':assert all(float(x['p.gen'])==1 for x in rows)
    else:assert all(float(x['p.gen'])==.05 for x in rows if x['ranking']=='global')
    assert (out/'10.Scoring_report_audit.html').exists()
    if mode!='absent':assert (out/'scoring_inputs/validated_inputs.json').exists()


@pytest.mark.parametrize('artifact',['aggregate','positions','matrix'])
def test_corrupt_canonical_null_rejected(project,artifact):
    prior=project/'prior';block=prior/'blocks/blk_000001'
    if artifact=='aggregate':
        p=prior/'null_positions.parquet';t=pq.read_table(p).to_pydict();t['nd']=[2];pq.write_table(pa.table(t),p)
    elif artifact=='positions':
        p=block/'positions.parquet';t=pq.read_table(p);pq.write_table(pa.concat_tables([t,t]),p)
    else:
        code="p<-'prior/blocks/blk_000001/matrices.rds';b<-readRDS(p);b$stat$global[1,1]<-.3;b$stat$top[1,1]<-.3;saveRDS(b,p)"
        r=subprocess.run(['Rscript','-e',code],cwd=project,text=True,capture_output=True);assert r.returncode==0,r.stderr
    r=command(project,ok=False)
    assert 'SCORING_EVIDENCE' not in {k for k,v in trace(project).items() if v=='COMPLETED'}
