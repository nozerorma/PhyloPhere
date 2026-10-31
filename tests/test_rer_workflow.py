"""Serial Nextflow integration and resume checks in a copied project."""
import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest

from test_rer import ROOT, rexpr


@pytest.fixture
def project(tmp_path):
    for name in ['RERCONVERGE','ENRICHMENT']:
        shutil.copytree(ROOT/'subworkflows'/name,tmp_path/'subworkflows'/name)
    style=tmp_path/'subworkflows/TRAIT_ANALYSIS/local/src';style.mkdir(parents=True)
    shutil.copy2(ROOT/'subworkflows/TRAIT_ANALYSIS/local/src/report_style.R',style/'report_style.R')
    (tmp_path/'workflows').mkdir()
    shutil.copy2(ROOT/'workflows/rerconverge.nf',tmp_path/'workflows/rerconverge.nf')
    (tmp_path/'main.nf').write_text('''include { RER_MAIN } from './workflows/rerconverge'
workflow {
    def traits = params.multiple ? Channel.of(file('traits.tsv'),file('traits.tsv')) : Channel.value(file('traits.tsv'))
    RER_MAIN(traits,Channel.empty(),Channel.empty(),Channel.empty())
    RER_MAIN.out.summary_tsv.map { it.name }.collectFile(name:'summary_files.txt',storeDir:params.outdir,newLine:true)
    RER_MAIN.out.perms.map { it.name }.collectFile(name:'null_files.txt',storeDir:params.outdir,newLine:true)
}
''')
    config=dict(traitname='trait',sp_colname='species',my_traits='traits.tsv',gene_trees='trees.nwk',
                rer_tool='',rer_trait_mode='auto',rer_transform='none',trait_out='',trees_out='',matrix_out='',
                rer_continuous_file='',rer_binary_file='',rer_perms_file='',rer_gene_scores='',
                rer_perm_batches=0,rer_perms_per_batch=2,rer_batch_seed=0,rer_minsp=3,rer_min_pos=2,
                winsorizeRER=0,winsorizeTrait=0,rer_binary_clade='all',rer_top_n_labels=3,
                evidence_p_threshold=.05,evidence_fdr_threshold=.05,use_singularity=False,use_apptainer=False,
                outdir='out',enrichment=False,scoring=False,multiple=False)
    (tmp_path/'params.json').write_text(json.dumps(config))
    (tmp_path/'audit.config').write_text("process.executor='local'\nprocess.cpus=1\nexecutor.queueSize=1\ntrace.overwrite=true\nprocess.beforeScript='export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1'\n")
    (tmp_path/'traits.tsv').write_text('species\ttrait\n'+''.join(f'{x}\t{i+1}\n' for i,x in enumerate('ABCDEFGHIJKLMNOP')))
    # At least 20 genes with complete coverage allow native weighted estimation.
    rexpr(tmp_path,"suppressPackageStartupMessages(library(ape));set.seed(61);t<-rtree(16,rooted=FALSE,tip.label=LETTERS[1:16]);ls<-vapply(1:24,function(i){t$edge.length<-runif(nrow(t$edge),.2,2);paste0('g.',i,'\\t',write.tree(t))},'');writeLines(ls,'trees.nwk');r<-data.frame(Rho=c(.4,-.5),N=c(8,8),P=c(.02,.03),p.adj=c(.03,.03),p.perm=c(.02,.1),row.names=c('001','NA'));saveRDS(r,'saved.rds')")
    return tmp_path


def run(project,*args,success=True):
    r=subprocess.run(['nextflow','-log','audit.log','run','main.nf','-params-file','params.json','-c','audit.config',
                      '-ansi-log','false','-with-trace','trace.tsv',*args],cwd=project,text=True,capture_output=True,timeout=240)
    if success: assert r.returncode==0,r.stdout+r.stderr+(project/'audit.log').read_text()[-4000:]
    else: assert r.returncode!=0,r.stdout+r.stderr
    return r.stdout+r.stderr


def trace(project):
    import csv
    with (project/'trace.tsv').open() as f:
        return list(csv.DictReader(f,delimiter='\t'))


def test_precomputed_and_stage_specific_resume(project):
    run(project,'--rer_continuous_file','saved.rds')
    results=project/'out/rerconverge/rer_results'
    assert (results/'rerconverge_summary_trait.tsv').exists()
    assert (results/'rer_inputs/checksums.json').exists()
    assert (project/'out/rerconverge/gene_lists/background.txt').read_text().splitlines()==['001','NA']
    run(project,'--rer_continuous_file','saved.rds','-resume')
    assert all(x['status']=='CACHED' for x in trace(project))
    rmd=project/'subworkflows/RERCONVERGE/local/5.RERconverge_report.Rmd'
    old=rmd.stat();text=rmd.read_text();rmd.write_text(text.replace('Relative Evolutionary Rate Convergence Analysis','Relative Evolutionary Rate Association Analysis'))
    # Equal-size content change with restored mtime must invalidate only the render.
    os.utime(rmd,ns=(old.st_atime_ns,old.st_mtime_ns))
    run(project,'--rer_continuous_file','saved.rds','-resume')
    statuses={x['name'].split(':')[-1].split(' ')[0]:x['status'] for x in trace(project)}
    assert statuses['RER_REPORT']=='COMPLETED' and statuses['RER_TABLES']=='CACHED' and statuses['RER_GENE_LISTS']=='CACHED'
    rexpr(project,"r<-readRDS('saved.rds');r$Rho[1]<-.7;saveRDS(r,'saved.rds')")
    run(project,'--rer_continuous_file','saved.rds','-resume')
    assert '0.7' in (results/'rerconverge_summary_trait.tsv').read_text()


@pytest.mark.parametrize('args,message',[
    (['--rer_continuous_file','saved.rds','--rer_binary_file','saved.rds'],'exactly one precomputed'),
    (['--rer_tool','continuous','--rer_continuous_file','saved.rds'],'new and precomputed'),
    (['--rer_tool','unknown'],'Unknown RER step'),
    (['--rer_tool','build_trait','--multiple','true'],'expected one trait table'),
    (['--rer_trait_mode','wrong'],'Invalid rer_trait_mode'),
])
def test_invalid_workflow_inputs(project,args,message):
    assert message in run(project,*args,success=False)


def test_real_continuous_inference_and_null(project):
    args=['--rer_tool','build_trait,build_tree,build_matrix,continuous','--rer_trait_mode','continuous',
          '--rer_perm_batches','1','--winsorizeRER','1','--winsorizeTrait','2']
    run(project,*args)
    results=project/'out/rerconverge/rer_results'
    meta=json.loads((results/'rer_run_metadata.json').read_text())
    assert meta['null_status']=='verified' and meta['n_perms']==2
    assert meta['inference']['transformation']['effective']=='none'
    assert meta['inference']['matrix_contract']['use.weights'] is True
    assert (results/'trait.continuous.perms.rds/rer_null_manifest.json').exists()
    rexpr(project,"source('subworkflows/RERCONVERGE/local/rer_permulations.R');r<-readRDS('out/rerconverge/rer_results/trait.continuous.output');m<-readRDS('out/rerconverge/rer_objects/trait.RERmatrix.output');p<-readRDS('out/rerconverge/rer_results/trait.char2path.output');ref<-rer_score_paths(m,p,'continuous',3,winR=1,winT=2);stopifnot(sum(is.finite(r$Rho))>0,isTRUE(all.equal(r$Rho,ref$Rho)),isTRUE(all.equal(r$P,ref$P)))")
    run(project,*args,'-resume')
    assert all(x['status']=='CACHED' for x in trace(project))
    rmd=project/'subworkflows/RERCONVERGE/local/5.RERconverge_report.Rmd'
    stat=rmd.stat();rmd.write_text(rmd.read_text().replace('Relative Evolutionary Rate Convergence Analysis','Relative Evolutionary Rate Association Analysis'))
    os.utime(rmd,ns=(stat.st_atime_ns,stat.st_mtime_ns))
    run(project,*args,'-resume')
    statuses={x['name'].split(':')[-1].split(' ')[0]:x['status'] for x in trace(project)}
    assert statuses['RER_CONT']=='CACHED' and statuses['RER_MATRIX']=='CACHED' and statuses['RER_REPORT']=='COMPLETED'
    traits=project/'traits.tsv';stat=traits.stat();original=traits.read_text()
    traits.write_text(original.replace('A\t1\n','A\t2\n'))
    assert traits.stat().st_size==stat.st_size
    os.utime(traits,ns=(stat.st_atime_ns,stat.st_mtime_ns))
    run(project,*args,'-resume')
    statuses={x['name'].split(':')[-1].split(' ')[0]:x['status'] for x in trace(project)}
    assert statuses['RER_TRAIT']=='COMPLETED' and statuses['RER_CONT']=='COMPLETED'


def test_real_binary_inference(project):
    (project/'traits.tsv').write_text('species\ttrait\n'+''.join(f'{x}\t{int(i<4)}\n' for i,x in enumerate('ABCDEFGHIJKLMNOP')))
    run(project,'--rer_tool','build_trait,build_tree,build_matrix,binary','--rer_perm_batches','1','--rer_binary_clade','terminal','--winsorizeRER','1')
    meta=json.loads((project/'out/rerconverge/rer_results/rer_run_metadata.json').read_text())
    assert meta['analysis_mode']=='binary' and meta['null_status']=='verified'
    rexpr(project,"source('subworkflows/RERCONVERGE/local/rer_permulations.R');r<-readRDS('out/rerconverge/rer_results/trait.binary.output');m<-readRDS('out/rerconverge/rer_objects/trait.RERmatrix.output');p<-readRDS('out/rerconverge/rer_results/trait.fg_paths.output');ref<-rer_score_paths(m,p,'binary',3,min.pos=2,winR=1);stopifnot(sum(is.finite(r$Rho))>0,isTRUE(all.equal(r$Rho,ref$Rho)),isTRUE(all.equal(r$P,ref$P)))")
