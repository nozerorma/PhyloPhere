"""Real R harvest checks: strict policy, configured columns, matching and RNG."""
import csv
import json
from pathlib import Path
import shutil
import subprocess

import pytest

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/'subworkflows/CT/local/scripts/permulations.R'
pytestmark=pytest.mark.skipif(not shutil.which('Rscript'),reason='Run in phylophere')


@pytest.fixture
def inputs(tmp_path):
    (tmp_path/'tree.nwk').write_text('(((A:1,B:1):1,(C:1,D:1):1):1,(E:1,F:1):2);\n')
    (tmp_path/'cfg.tab').write_text('B\t1\t1\nA\t0\t1\nD\t1\t2\nC\t0\t2\nF\t1\t3\nE\t0\t3\n')
    values=[.1,.9,.2,.8,.3,.7]
    with (tmp_path/'values.csv').open('w') as out:
        out.write('family,sample id,sample_size,score rate,n,c\n')
        out.writelines(f'F,{s},{i+1},{v},100,{int(v*100)}\n' for i,(s,v) in enumerate(zip('ABCDEF',values)))
    return tmp_path


def harvest(p,name='out',*,cycles=3,workers=1,fop=False,match=False,relaxed=False,
            column='score rate',n='',c='',use_n=False,mode='continuous',strategy='BM',max_fop=3):
    args=['Rscript',str(SCRIPT),str(p/'tree.nwk'),str(p/'cfg.tab'),str(cycles),strategy,
          str(p/'values.csv'),str(p/name),'2','0','20',column,n,c,str(use_n).lower(),
          mode,str(fop).lower(),str(max_fop),str(workers),'1998',str(match).lower(),'.25','none',
          str(relaxed).lower(),'sample id']
    return subprocess.run(args,text=True,capture_output=True,timeout=90)


def test_auto_default_and_spaced_headers_select_the_requested_column(inputs):
    r=harvest(inputs,cycles=1,strategy='auto')
    assert r.returncode==0,r.stdout+r.stderr
    assert "value column 'score rate'" in r.stdout
    assert 'sample_size\'' not in r.stdout


@pytest.mark.parametrize('kwargs,pattern',[
    ({'column':'typo'},'Phenotype column not found'),
    ({'n':'typo','c':'c','use_n':True},'column not found'),
    ({'strategy':'typo'},'Unknown perm_strategy'),
    ({'cycles':-1},'caas_full_perms'),
    ({'max_fop':0},'max_fop'),
])
def test_bad_requested_inputs_fail_instead_of_falling_back(inputs,kwargs,pattern):
    r=harvest(inputs,**kwargs)
    assert r.returncode!=0 and pattern in r.stdout+r.stderr


def test_fractional_counts_are_rejected(inputs):
    path=inputs/'values.csv';path.write_text(path.read_text().replace(',100,10\n',',100,10.5\n'))
    r=harvest(inputs,n='n',c='c',use_n=True)
    assert r.returncode!=0 and 'Invalid integer counts' in r.stderr+r.stdout


def test_serial_parallel_fop_are_identical_and_preserve_the_marginal(inputs):
    for workers in [1,2]:
        r=harvest(inputs,name=f'out{workers}',workers=workers,fop=True)
        assert r.returncode==0,r.stdout+r.stderr
    for name in ['resample_001.tab','resample_002.tab','fop_labelings.tab','fop_pairs.tsv','permulation_manifest.tsv']:
        assert (inputs/'out1'/name).read_bytes()==(inputs/'out2'/name).read_bytes()
    rows=list(csv.DictReader((inputs/'out1/permulation_manifest.tsv').open(),delimiter='\t'))
    assert len(rows)==3 and {r['tier'] for r in rows}=={'1'}
    assert all(sorted(map(float,(r['fg_values']+','+r['bg_values']).split(',')))==[.1,.2,.3,.7,.8,.9] for r in rows)
    assert not (inputs/'out1/.fop_candidates').exists()
    prepared=inputs/'prepared';prepared.mkdir()
    r=subprocess.run(['python',str(SCRIPT.parent/'prepare_labelings.py'),'--source',str(inputs/'out1'),
                      '--config',str(inputs/'cfg.tab'),'--n','3','--block-size','2','--multi','false'],
                     cwd=prepared,text=True,capture_output=True)
    assert r.returncode==0,r.stdout+r.stderr
    assert json.loads((prepared/'ct_preparation.json').read_text())['harvest_provenance_available']



def test_requested_pss_matching_is_strict_and_relaxation_is_recorded(inputs):
    # B/D cannot be the ordinal extreme pair, so its requested profile is unavailable.
    (inputs/'cfg.tab').write_text('B\t1\t1\nD\t0\t1\n')
    strict=harvest(inputs,name='strict',cycles=1,match=True,mode='ordinal')
    assert strict.returncode!=0 and 'Strict null' in strict.stdout+strict.stderr
    relaxed=harvest(inputs,name='relaxed',cycles=1,match=True,mode='ordinal',relaxed=True)
    assert relaxed.returncode==0,relaxed.stdout+relaxed.stderr
    provenance=json.loads((inputs/'relaxed/permulation_provenance.json').read_text())
    assert provenance['allow_relaxed_null'] and provenance['relaxed_reasons']
    assert not provenance['match_pss_applied']


@pytest.mark.parametrize('counts',[False,True])
def test_observed_bundle_is_reused_and_changed_settings_are_rejected(inputs,counts):
    # Use the same preparation API as the observed composition report.
    code='''
source('subworkflows/TRAIT_ANALYSIS/local/src/selection_algorithm.R')
source('subworkflows/TRAIT_ANALYSIS/local/src/contrast_inputs.R')
p<-commandArgs(TRUE)[1];counts<-commandArgs(TRUE)[2]=='true';d<-read.csv(file.path(p,'values.csv'),check.names=FALSE)
d$species<-d[['sample id']];tree<-ape::read.tree(file.path(p,'tree.nwk'))
ci<-data.frame(species=d$species,trait=d[['score rate']],n_population=NA_real_,n_cases=NA_real_,trait_ci_lb=NA_real_,trait_ci_ub=NA_real_)
if(counts) {
 intervals<-binom::binom.confint(d$c,d$n,methods='bayes',prior.shape1=.5,prior.shape2=.5)
 ci$n_population<-d$n;ci$n_cases<-d$c;ci$trait_ci_lb<-intervals$lower;ci$trait_ci_ub<-intervals$upper
}
b<-prepare_observed_contrasts(ci,tree,0,FALSE,'none','BM');b$trait_ci<-ci
b$provenance<-contrast_provenance(d,tree,'score rate',if(counts)'n' else '',if(counts)'c' else '',0,FALSE,'none','BM',1998)
saveRDS(b,file.path(p,'observed_selection.rds'))
'''
    script=inputs/'bundle.R';script.write_text(code)
    r=subprocess.run(['Rscript',str(script),str(inputs),str(counts).lower()],cwd=ROOT,text=True,capture_output=True,timeout=60)
    assert r.returncode==0,r.stderr
    r=harvest(inputs,cycles=1,use_n=counts,n='n' if counts else '',c='c' if counts else '')
    assert r.returncode==0,r.stdout+r.stderr
    assert json.loads((inputs/'out/permulation_provenance.json').read_text())['observed_model_reused']
    r=harvest(inputs,name='changed',cycles=1,strategy='OU',use_n=counts,n='n' if counts else '',c='c' if counts else '')
    assert r.returncode!=0 and 'bundle does not match' in r.stdout+r.stderr


def test_strict_multi_hypothesis_harvest_matches_real_observed_bundle(inputs):
    code='''
source('subworkflows/TRAIT_ANALYSIS/local/src/selection_algorithm.R')
source('subworkflows/TRAIT_ANALYSIS/local/src/contrast_inputs.R')
p<-commandArgs(TRUE)[1];set.seed(1998)
tree<-ape::read.tree(text='(((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1,((I:1,J:1):1,(K:1,L:1):1):1);')
ape::write.tree(tree,file.path(p,'tree.nwk'))
d<-data.frame('sample id'=LETTERS[1:12],'score rate'=rep(c(.1,.9,.2,.8),3),check.names=FALSE)
write.csv(d,file.path(p,'values.csv'),row.names=FALSE);d$species<-d[['sample id']]
ci<-data.frame(species=d$species,trait=d[['score rate']],n_population=NA_real_,n_cases=NA_real_,trait_ci_lb=NA_real_,trait_ci_ub=NA_real_)
b<-prepare_observed_contrasts(ci,tree,0,FALSE,'none','BM');b$trait_ci<-ci
b$provenance<-contrast_provenance(d,tree,'score rate','','',0,FALSE,'none','BM',1998)
saveRDS(b,file.path(p,'observed_selection.rds'))
f<-fop_pair_sel.f(b$ctx,ordinal=FALSE,top_pct=0,max_contrasts=3,max_fop=3,seed=1998,candidate_result=b$candidates)
stopifnot(length(f$hypotheses)==3)
rows<-list()
for(h in names(f$hypotheses)) {
 hd<-f$hypotheses[[h]]
 write.table(gen_traitfiles(hd),file.path(p,paste0('traitfile_',h,'.tab')),sep='\\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
 rows[[h]]<-data.frame(hypothesis_id=h,pair=hd$cluster,species1=hd$species1,species2=hd$species2,pss_score=hd$pss_score)
}
file.copy(file.path(p,'traitfile_H1.tab'),file.path(p,'cfg.tab'),overwrite=TRUE)
write.table(do.call(rbind,rows),file.path(p,'contrast_hypotheses_pairs.tsv'),sep='\\t',row.names=FALSE,quote=FALSE)
'''
    script=inputs/'observed_multi.R';script.write_text(code)
    r=subprocess.run(['Rscript',str(script),str(inputs)],cwd=ROOT,text=True,capture_output=True,timeout=60)
    assert r.returncode==0,r.stdout+r.stderr
    r=harvest(inputs,name='multi',cycles=2,fop=True,match=True)
    assert r.returncode==0,r.stdout+r.stderr
    labels=(inputs/'multi/fop_labelings.tab').read_text().splitlines()
    assert len(labels)==6 and {row.split('\t')[0] for row in labels}=={
        f'b_{cycle}~H{h}' for cycle in [1,2] for h in [1,2,3]}
    provenance=json.loads((inputs/'multi/permulation_provenance.json').read_text())
    assert provenance['observed_model_reused'] and not provenance['relaxed_reasons']
    assert provenance['match_pss_applied']


def test_existing_harvest_artifacts_are_not_mixed_or_overwritten(inputs):
    r=harvest(inputs,cycles=1)
    assert r.returncode==0,r.stdout+r.stderr
    original=(inputs/'out/resample_001.tab').read_bytes()
    r=harvest(inputs,cycles=2)
    assert r.returncode!=0 and 'fresh output directory' in r.stdout+r.stderr
    assert (inputs/'out/resample_001.tab').read_bytes()==original
