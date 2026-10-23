"""WP7 integrity and explicit compatibility export contracts."""
import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]
R = shutil.which('Rscript')
CLI = ROOT / 'subworkflows/SCORING/local/src/scoring_caas_perms.R'
PREPARE = ROOT / 'subworkflows/ENRICHMENT/local/src/fcs_prepare_blocks.R'

def run(tmp, *args):
    if not R:
        pytest.skip('Rscript required')
    return subprocess.run([R, *map(str,args)], cwd=tmp, env=dict(os.environ,TMPDIR=str(tmp)), capture_output=True,text=True)

def fixture(tmp):
    result = run(tmp, '-e', '''
      dir.create('blk_0001');dir.create('blk_0002')
      for(i in 1:2) {
        g<-c('G1','G2');cy<-paste0('b_',i)
        m<-matrix(c(0.5,0),2,1,dimnames=list(g,cy))
        b<-list(genes=g,cycles=cy,gene_stat='max',n_units=1,
                stat=list(global=m,top=m,bottom=m))
        saveRDS(b,sprintf('blk_%04d/matrices.rds',i))
        writeLines(cy,sprintf('blk_%04d/cycles.txt',i))
      }
      writeLines(c('b_1','b_2'),'roster.txt')
    ''')
    assert result.returncode==0,result.stderr
    return sorted(tmp.glob('blk_*/matrices.rds'))

@pytest.mark.parametrize('damage', ['missing','duplicate_roster','bad_gene','bad_order','gene_stat','overlap','missing_rank'])
def test_corruption_is_rejected_before_export(tmp_path,damage):
    paths=fixture(tmp_path)
    code={
        'missing': "unlink('blk_0002/matrices.rds')",
        'duplicate_roster': "writeLines(c('b_1','b_2','b_2'),'roster.txt')",
        'bad_gene': "b$genes<-c('G1','G1')",
        'bad_order': "rownames(b$stat$global)<-rev(b$genes)",
        'gene_stat': "b$gene_stat<-'fmax'",
        'overlap': "b$cycles<-'b_1'",
        'missing_rank': "b$stat$bottom<-NULL",
    }[damage]
    if damage not in ('missing','duplicate_roster'):
        code="b<-readRDS('blk_0002/matrices.rds');"+code+";saveRDS(b,'blk_0002/matrices.rds')"
    assert run(tmp_path,'-e',code).returncode==0
    result=run(tmp_path,CLI,'--blocks',*paths,'--cycles','roster.txt','--output','export.rds')
    assert result.returncode!=0
    assert not (tmp_path/'export.rds').exists()

def test_validate_only_preserves_empty_cycles_without_global_rds(tmp_path):
    paths=fixture(tmp_path)
    run(tmp_path,'-e',"b<-readRDS('blk_0002/matrices.rds');b$n_units<-0;b$stat<-lapply(b$stat,function(m) m*0);saveRDS(b,'blk_0002/matrices.rds')")
    result=run(tmp_path,CLI,'--blocks',*paths,'--cycles','roster.txt','--validate-only','--metadata','null_metadata.json')
    assert result.returncode==0,result.stderr
    assert not (tmp_path/'caas_perms.rds').exists()
    meta=json.loads((tmp_path/'null_metadata.json').read_text())
    assert meta['cycles']==['b_1','b_2']
    assert meta['gene_stat']=='max'
    result=run(tmp_path,CLI,'--blocks',*paths,'--cycles','roster.txt','--output','export.rds')
    assert result.returncode==0,result.stderr
    result=run(tmp_path,'-e',"x<-readRDS('export.rds');stopifnot(identical(colnames(x$caas_corStat_byrank$global),c('b_1','b_2')),all(x$caas_corStat_byrank$global[,2]==0))")
    assert result.returncode==0,result.stderr

def test_saved_store_missing_block_is_rejected(tmp_path):
    paths=fixture(tmp_path)
    store=tmp_path/'live/store'
    (store/'blocks').mkdir(parents=True)
    for p in paths: shutil.move(p.parent,store/'blocks'/p.parent.name)
    shutil.copy(tmp_path/'roster.txt',store/'cycle_roster.txt')
    shutil.rmtree(store/'blocks/blk_0002')
    result=run(tmp_path,PREPARE,'caas','max',7,'enabled')
    assert result.returncode!=0 and 'missing from blocks' in result.stderr

def test_export_parameter_roundtrips_and_renders():
    from gui.generation.render import render_single, render_batch
    from gui.project_io import load_project
    from gui.models.modules import CaasConfig
    assert CaasConfig().caas_export_single is False
    p=load_project(ROOT/'gui/templates/tier2_marine.json')
    assert '"caas_export_single": false' in render_single(p)
    p.modules.caas.caas_export_single=True
    assert '"caas_export_single": true' in render_single(p)
    assert render_batch(p)

@pytest.mark.parametrize('damage', ['roster_extra','universe_duplicate','parts_stat','empty_cycle_id'])
def test_additional_integrity_guards(tmp_path,damage):
    paths=fixture(tmp_path)
    extra=[]
    if damage=='roster_extra':
        (tmp_path/'roster.txt').write_text('b_1\nb_2\nb_3\n')
    elif damage=='universe_duplicate':
        (tmp_path/'universe.txt').write_text('G1\nG1\n')
        extra=['--universe','universe.txt']
    else:
        code="b<-readRDS('blk_0002/matrices.rds');"
        if damage=='parts_stat':
            code+="b$parts<-list(max=lapply(b$stat,function(m) m*0),fmax=b$stat,n=b$stat);"
        else:
            code+="b$cycles<-'';"
        code+="saveRDS(b,'blk_0002/matrices.rds')"
        assert run(tmp_path,'-e',code).returncode==0
    result=run(tmp_path,CLI,'--blocks',*paths,'--cycles','roster.txt',*extra,'--validate-only')
    assert result.returncode!=0

def test_saved_metadata_content_mismatch_is_rejected(tmp_path):
    paths=fixture(tmp_path)
    result=run(tmp_path,CLI,'--blocks',*paths,'--cycles','roster.txt','--validate-only','--metadata','null_metadata.json')
    assert result.returncode==0,result.stderr
    store=tmp_path/'live/store'
    (store/'blocks').mkdir(parents=True)
    for p in paths: shutil.move(p.parent,store/'blocks'/p.parent.name)
    shutil.copy(tmp_path/'roster.txt',store/'cycle_roster.txt')
    meta=json.loads((tmp_path/'null_metadata.json').read_text());meta['genes']=['wrong']
    (store/'null_metadata.json').write_text(json.dumps(meta))
    result=run(tmp_path,PREPARE,'caas','max',7,'enabled')
    assert result.returncode!=0 and 'metadata disagrees' in result.stderr

def test_reducer_refuses_a_missing_whole_block_against_replay_input(tmp_path):
    # Even a cycle with no detections must have a block; expected cycles come from replay input.
    from test_wp5_null_blocks_cli import _tile, _merge, _labelings, LOCAL
    tile=_tile(tmp_path/'tile',{})
    result,block=_merge(tmp_path,[tile],_labelings(tmp_path/'one.tab',['b_1']),'block')
    assert result.returncode==0,result.stderr
    all_labels=_labelings(tmp_path/'all.tab',['b_0','b_1','b_2'])
    result=subprocess.run([os.sys.executable,str(LOCAL/'reduce_null_blocks.py'),'--blocks',str(block),
        '--labelings',str(all_labels),'--output-dir',str(tmp_path/'out')],capture_output=True,text=True)
    assert result.returncode!=0 and 'missing or unexpected block cycles' in result.stderr
    assert not (tmp_path/'out/cycle_roster.txt').exists()

def test_missing_entire_saved_blocks_directory_is_an_error(tmp_path):
    (tmp_path/'live/store').mkdir(parents=True)
    (tmp_path/'live/store/cycle_roster.txt').write_text('b_1\n')
    result=run(tmp_path,PREPARE,'caas','max',7,'enabled')
    assert result.returncode!=0 and 'Missing canonical' in result.stderr

def test_valid_saved_metadata_is_accepted(tmp_path):
    paths=fixture(tmp_path)
    result=run(tmp_path,CLI,'--blocks',*paths,'--cycles','roster.txt','--validate-only','--metadata','null_metadata.json')
    assert result.returncode==0,result.stderr
    store=tmp_path/'live/store'
    (store/'blocks').mkdir(parents=True)
    for p in paths: shutil.move(p.parent,store/'blocks'/p.parent.name)
    shutil.copy(tmp_path/'roster.txt',store/'cycle_roster.txt')
    shutil.copy(tmp_path/'null_metadata.json',store/'null_metadata.json')
    result=run(tmp_path,PREPARE,'caas','max',7,'enabled')
    assert result.returncode==0,result.stderr
    assert len(list((tmp_path/'blocks').glob('*.rds')))==2
