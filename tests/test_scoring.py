"""SCORING import contracts and empirical universes, without inference engines."""
import csv
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'subworkflows/SCORING/local/src'
sys.path.insert(0,str(SRC))
from scoring_input import read_rows
spec=importlib.util.spec_from_file_location('scoring_null',SRC/'null_position_stats.py')
null=importlib.util.module_from_spec(spec);spec.loader.exec_module(null)


def run(args,wd):
    return subprocess.run(args,cwd=wd,text=True,capture_output=True,timeout=90)


def observed(path,genes=('001','NA','a.b','a"b',"a'b")):
    with path.open('w') as fh:
        w=csv.writer(fh,delimiter='\t',lineterminator='\n');w.writerow(['Gene','Position','side','asr_path_score'])
        for g in genes:w.writerow([g,0,'top',.8])
        w.writerow([genes[0],1,'none',0])
    return genes


@pytest.mark.parametrize('score,side,pos', [('NaN','top','0'),('inf','top','0'),('1.1','top','0'),('0','top','0'),('.1','none','0'),('.1','wrong','0'),('.1','top','-1'),('.1','top','0.5')])
def test_invalid_observed_rejected(tmp_path,score,side,pos):
    p=tmp_path/'p.tsv';p.write_text(f'Gene\tPosition\tside\tasr_path_score\ng\t{pos}\t{side}\t{score}\n')
    with pytest.raises(ValueError):read_rows(p)


def test_literal_compute_evidence_and_duplicate_sources(tmp_path):
    for f in SRC.glob('*'):
        if f.is_file():shutil.copy(f,tmp_path/f.name)
    p=tmp_path/'input.tsv';genes=observed(p)
    args=[sys.executable,str(SRC/'observed_core_scores.py'),'--input',str(p),'--positions-out','core_pos.tsv','--genes-out','core_gene.tsv']
    r=run(args,tmp_path);assert r.returncode==0,r.stderr
    # Raw TSV from current FADE/RER producers: embedded quotes are literal.
    (tmp_path/'fade.tsv').write_text('gene\tmax_bf\n'+''.join(f'{g}\t101\n' for g in genes))
    (tmp_path/'rer.tsv').write_text('gene\tRho\tP\tp.perm\tp.perm.adj\n'+''.join(f'{g}\t.5\t.001\t.03\t.03\n' for g in genes))
    compute=['Rscript','scoring_compute.R','--postproc','input.tsv','--core_positions','core_pos.tsv','--core_genes','core_gene.tsv','--fade_top','fade.tsv','--rer','rer.tsv','--nominal_threshold','.02','--fdr_threshold','.02']
    r=run(compute,tmp_path);assert r.returncode==0,r.stdout+r.stderr
    (tmp_path/'universe.txt').write_text('\n'.join(genes)+'\n');(tmp_path/'blocks').mkdir()
    r=run(['Rscript',str(SRC/'scoring_evidence.R'),'--genes','gene_scores.tsv','--positions','position_scores.tsv','--stats','fcs_stats.tsv','--universe','universe.txt','--blocks','blocks'],tmp_path)
    assert r.returncode==0,r.stdout+r.stderr
    for name,col in [('gene_scores.tsv','Gene'),('position_scores.tsv','Gene'),('fcs_stats.tsv','gene')]:
        with (tmp_path/name).open() as fh:rows=list(csv.DictReader(fh,delimiter='\t',quoting=csv.QUOTE_NONE))
        assert {x[col] for x in rows}==set(genes)
    with (tmp_path/'fcs_stats.tsv').open() as fh:stats=list(csv.DictReader(fh,delimiter='\t'))
    assert all(x['flag_rer_acc']=='FALSE' for x in stats)
    (tmp_path/'fade.tsv').write_text('gene\tmax_bf\n001\t101\n001\t102\n')
    r=run(compute,tmp_path);assert r.returncode!=0 and 'Duplicate FADE' in r.stderr


def test_fixed_family_and_strict_membership(tmp_path):
    (tmp_path/'bg').write_text('001\t0,1,2\nNA\t0,1\n')
    (tmp_path/'genes').write_text('001\nNA\n')
    t=pa.table(dict(gene=['001'],position=[0],nd=[1],first_rank=[0]))
    classes,_,_,metrics=null.build_universe(t,['001'],np.array([2]),tmp_path/'bg',tmp_path/'genes',['b_1','b_2'])
    assert metrics['family']==5 and classes.num_rows==1
    with pytest.raises(ValueError,match='outside declared'):
        null.build_universe(t,['001'],np.array([99]),tmp_path/'bg',tmp_path/'genes',['b_1'])
    with pytest.raises(ValueError,match='explicit tested'):
        null.build_universe(t,[],np.array([],dtype=int),None,tmp_path/'genes',['b_1'])


def test_count_duplicate_cycle_and_reordered_reduce(tmp_path):
    p=tmp_path/'p.parquet';classes=pa.table(dict(gene=['g'],position=[0],nd=[2],cls=[1]))
    pq.write_table(pa.table(dict(gene=['g','g'],cycle=['b_1','b_1'],position=[0,0],all=[.8,.8])),p)
    with pytest.raises(ValueError,match='Duplicate'):
        null.count_block(p,classes,np.array([]),['g'],np.array([0]),np.array([.8]))
    pq.write_table(pa.table(dict(gene=['h'],position=[0],k_emp=[1])),p)
    with pytest.raises(ValueError,match='keys/order'):
        null.reduce_counts([p],classes,['g'],np.array([0]),np.array([.8]),1)


def test_silent_null_fixed_bh_and_counter_bounds(tmp_path):
    code=f'''source('{SRC}/null_position_p.R')
obs<-data.frame(Gene=letters[1:10],Position=0,.obs=.8)
ctr<-data.frame(Gene=letters[1:10],Position=0,k_emp=0,nd=0)
x<-null_position_p(obs,ctr,NULL,0,99,n_family=100)
stopifnot(all(x$pos$p.emp==.01),all(x$pos$p.adj_bh==.1),x$n_family==100)
y<-try(null_position_p(obs,transform(ctr,k_emp=1),NULL,0,99,n_family=100),silent=TRUE)
stopifnot(inherits(y,'try-error'))
'''
    r=run(['Rscript','-e',code],tmp_path);assert r.returncode==0,r.stderr


def test_roundtrip_precision_and_literal_ids(tmp_path):
    code=f'''source('{SRC}/scoring_io.R')
x<-data.frame(Gene=c('001','NA','a"b'),p=c(.05-1e-16,.05+1e-16,1/19),stringsAsFactors=FALSE)
scoring_write(x,'precision.tsv');y<-scoring_read('precision.tsv')
stopifnot(identical(x$Gene,y$Gene),identical(x$p,y$p))
'''
    r=run(['Rscript','-e',code],tmp_path);assert r.returncode==0,r.stderr


def test_numeric_missing_is_distinct_from_invalid(tmp_path):
    code=f'''source('{SRC}/scoring_io.R')
stopifnot(is.na(scoring_number(NA_real_,'p',0,1,missing=TRUE)))
for(x in list(NaN,'oops',Inf,-.1,1.1,TRUE,FALSE)) stopifnot(inherits(try(scoring_number(x,'p',0,1,missing=TRUE),silent=TRUE),'try-error'))
stopifnot(scoring_number(Inf,'BF',0,Inf,positive_inf=TRUE)==Inf)
'''
    r=run(['Rscript','-e',code],tmp_path);assert r.returncode==0,r.stderr
