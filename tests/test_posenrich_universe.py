"""Fixed tested universe, literal IDs, scientific provenance and empirical count families."""
import json
import sys
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from posenrich_fixtures import canonical, digest, seal_block
ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'subworkflows/POSENRICH/local/src'
sys.path.insert(0,str(SRC))
import posenrich_universe as u
from posenrich_contract import observed_table, verify_universe
from accumulation_test import run as accumulation
from argparse import Namespace


def inputs(tmp_path, silent=False):
    genes=['001','NA','G.v1','"Q"','empty']
    (tmp_path/'genes').write_text('\n'.join(genes)+'\n')
    (tmp_path/'tested').write_text(''.join(f'{g}\t'+('' if g=='empty' else '0,1,2')+'\n' for g in genes))
    obs=tmp_path/'observed.tsv'
    pd.DataFrame([['001',0,'top',.8],['001',0,'bottom',.4],['NA',1,'top',.5],['"Q"',2,'none',0]],columns=['Gene','Position','side','asr_path_score']).to_csv(obs,sep='\t',index=False)
    block=canonical(tmp_path/'null',obs,genes,[] if silent else [('001','b_1',0,.9,np.nan),('NA','b_2',0,np.nan,.3),('001','b_2',99,.5,np.nan)],tested=tmp_path/'tested')
    args=[sys.executable,str(SRC/'posenrich_universe.py'),'--background-genes',str(tmp_path/'genes'),'--background-positions',str(tmp_path/'tested'),'--null-positions',str(tmp_path/'null/null_positions.parquet'),'--observed',str(obs),'--observed-contract',str(obs)+'.scoring.json','--roster',str(tmp_path/'null/cycle_roster.txt'),'--counts',str(tmp_path/'null/counts.json'),'--blocks',str(block),'--output-dir',str(tmp_path/'universe')]
    return genes,block,args


def test_fixed_family_keeps_undetected_and_silent_positions(tmp_path):
    genes,block,args=inputs(tmp_path)
    r=subprocess.run(args,capture_output=True,text=True);assert r.returncode==0,r.stderr
    t=pd.read_csv(tmp_path/'universe/universe.tsv.gz',sep='\t',dtype={'gene':str},keep_default_na=False)
    assert len(t)==12 and set(t.gene)==set(genes)-{'empty'}
    assert t.observed_positive.sum()==2 and t.nd.sum()==2
    m=verify_universe(tmp_path/'universe/universe.tsv.gz',blocks=[block])
    assert m['genes']==genes and m['null']['N']==5
    output=tmp_path/'acc/accumulation_all_results.csv'
    accumulation(Namespace(universe=str(tmp_path/'universe/universe.tsv.gz'),universe_manifest=str(tmp_path/'universe/universe_manifest.json'),blocks=[block],change_side='both',output=str(output)))
    a=pd.read_csv(output,keep_default_na=False).set_index('Gene')
    assert float(a.loc['001','EmpiricalP'])==pytest.approx(2/6)
    assert float(a.loc['NA','EmpiricalP'])==pytest.approx(2/6)
    assert float(a.loc['G.v1','EmpiricalP'])==1
    assert a.loc['empty','status']=='no_tested_positions' and a.loc['empty','PValue']=='NA'
    assert float(a.loc['001','EmpiricalQ'])==pytest.approx(2/3)


def test_entirely_silent_null_is_valid(tmp_path):
    _,_,args=inputs(tmp_path,True)
    r=subprocess.run(args,capture_output=True,text=True);assert r.returncode==0,r.stderr
    t=pd.read_csv(tmp_path/'universe/universe.tsv.gz',sep='\t')
    assert len(t)==12 and t.nd.sum()==0


@pytest.mark.parametrize('kind',['duplicate','bad_score','bad_side','fractional','quotes'])
def test_observed_contract(tmp_path,kind):
    p=tmp_path/'s.tsv';text='Gene\tPosition\tside\tCAAS_score\n001\t0\ttop\t.5\nNA\t1\tbottom\t.4\n"Q"\t2\ttop\t.8\n'
    if kind=='duplicate':text+='001\t0\ttop\t.5\n'
    elif kind=='bad_score':text=text.replace('.5','inf')
    elif kind=='bad_side':text=text.replace('top','none')
    elif kind=='fractional':text=text.replace('001\t0','001\t0.5')
    p.write_text(text)
    if kind=='quotes':assert observed_table(p).Gene.tolist()==['001','NA','"Q"']
    else:
        with pytest.raises(ValueError):observed_table(p)


def test_observed_outside_tested_rejected_and_no_fallback(tmp_path):
    genes,_,_=inputs(tmp_path)
    null=u.load_null_positive(tmp_path/'null/null_positions.parquet')
    obs=u.load_observed_positive(tmp_path/'observed.tsv');obs.loc[0,'position']=88
    with pytest.raises(ValueError,match='outside declared'):u.build_universe(genes,null,obs,tmp_path/'tested')
    with pytest.raises(ValueError,match='Explicit tested'):u.build_universe(genes,null,obs,None)


@pytest.mark.parametrize('kind',['checksum','scientific','aggregate','missing_roster','duplicate_tested'])
def test_incompatible_inputs_fail(tmp_path,kind):
    _,block,args=inputs(tmp_path)
    if kind=='checksum':(block/'positions.parquet').write_bytes(b'corrupt')
    elif kind=='scientific':
        p=tmp_path/'observed.tsv.scoring.json';m=json.loads(p.read_text());m['scientific']['model']='WAG';p.write_text(json.dumps(m))
    elif kind=='aggregate':
        p=tmp_path/'null/null_positions.parquet';t=pd.read_parquet(p);t.loc[0,'nd']=4;t.to_parquet(p,index=False)
    elif kind=='missing_roster':(block/'cycles.txt').unlink()
    else:(tmp_path/'tested').write_text('001\t0,0\n')
    r=subprocess.run(args,capture_output=True,text=True);assert r.returncode!=0
