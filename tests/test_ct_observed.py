"""Observed contracts and shared observed/null scientific invariants."""
import csv
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
LOCAL = ROOT / 'subworkflows/CT_DISAMBIGUATION/local'
sys.path.insert(0, str(LOCAL))
from src.core.labelings import read_trait_pairs, observed_pss
from src.core.discovery import DiscoveryRows
from src.core.observed import is_strict_row, observed_entries
from src.utils.io_utils import find_gene_alignment, read_alignment, read_alignment_manifest
from src.convergence.path_scores import compute_domain_scores
from src.convergence.fop_pool import pool_domains
from src.convergence.disambiguate_single import _derive_convergent_call
from src.asr import asr_single as asr
from src.asr.cache import prepare_config, validate_entry


def put(p, text):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)
    return p


@pytest.fixture
def inputs(tmp_path):
    tmp_path=tmp_path/"input's directory"
    put(tmp_path/'design.tab','A\t1\t1\nB\t0\t1\nC\t1\t2\nD\t0\t2\nE\t1\t3\nF\t0\t3\nG\t1\t4\nH\t0\t4\n')
    put(tmp_path/'tree.nwk','(((A:.1,B:.1):.1,(C:.1,D:.1):.1):.1,((E:.1,F:.1):.1,(G:.1,H:.1):.1):.1);\n')
    base='ACDEFGHIKLMNPQRSTVWY'*4
    seqs=''.join(f'>{s}\n'+('D' if s in 'ACEG' else 'A')+base[1:]+'\n' for s in 'ABCDEFGH')
    put(tmp_path/'alignments/AC010616.1.fa',seqs)
    put(tmp_path/'discovery.tab','gene\tposition\ttrait\tcaas\tpattern\nAC010616.1\t0\tH1\tDDDD/AAAA\t1\n')
    return tmp_path


@pytest.mark.parametrize('text',[
 'A\t1\t2\nB\t0\t2\n', 'A\t1\t1\nA\t0\t1\n',
 'A\t1\t1\nC\t1\t1\nB\t0\t1\n', 'A\t1\t1\n', 'A\t2\t1\nB\t0\t1\n'])
def test_malformed_designs_fail(tmp_path,text):
    with pytest.raises(ValueError): read_trait_pairs(put(tmp_path/'design.tab',text))


@pytest.mark.parametrize('row', ['H1\t1.5\tA\tB\t1', 'H1\t1\tA\tB\tnan',
                               'H1\t1\tA\tB\t-1', 'H1\t1\tB\tA\t1', 'H2\t1\tA\tB\t1'])
def test_invalid_supplied_weights_fail(tmp_path,row):
    p=put(tmp_path/'weights.tsv','hypothesis_id\tpair\tspecies1\tspecies2\tpss_score\n'+row+'\n')
    with pytest.raises(ValueError): observed_pss(p,{1:[('A','B')]})


@pytest.mark.parametrize('extra', [dict(position='-1'), dict(trait='H2'), dict(pattern='9'),
                                  dict(shared_residue_count='1'),dict(is_conserved_meta='TRUE')])
def test_discovery_validation_precedes_scoring(tmp_path,extra):
    row=dict(gene='G.1',position='0',trait='H1',caas='DD/AA',pattern='1');row.update(extra)
    p=tmp_path/'discovery.tab'
    with p.open('w') as fh:
        writer=csv.DictWriter(fh,fieldnames=list(row),delimiter='\t');writer.writeheader();writer.writerow(row)
    with pytest.raises(ValueError): DiscoveryRows(p,tmp_path,{1:[('A','B')]})


def test_duplicate_hypothesis_rows_fail_and_grouping_is_disk_backed(inputs):
    p=inputs/'discovery.tab';p.write_text(p.read_text()+p.read_text().splitlines()[1]+'\n')
    with pytest.raises(ValueError,match='Duplicate discovery'): DiscoveryRows(p,inputs)
    p.write_text('\n'.join(p.read_text().splitlines()[:2])+'\n')
    rows=DiscoveryRows(p,inputs)
    assert not isinstance(rows,dict) and list(rows)==['AC010616.1']
    rows.close()


@pytest.mark.parametrize('aa',['x','X','?','-','.','B','Z','J','*'])
def test_missing_residues_never_score_as_derived_states(aa):
    root=SimpleNamespace(node_id=1,parent=None)
    pairs=[dict(pair_id=d,node_id=1,focal_state='A',top_tip_mode=aa,bottom_tip_mode='A') for d in [1,2]]
    sides=compute_domain_scores(pairs,{1:{'A':1.0}},{1:root})
    assert not sides['top']['domain_scores'] and not sides['top']['domain_der']


def test_strict_overlap_and_asymmetric_pair_annotations():
    assert not is_strict_row(dict(caas='DDA/AEE',is_conserved_meta='TRUE',conserved_pair='0:',shared_residue_count='1'))
    assert not is_strict_row(dict(caas='DA/AD',conserved_pair='0:'))
    assert is_strict_row(dict(caas='DD/AA',conserved_pair='0:'))
    entry=observed_entries('G',[dict(position='0',trait='H1',pattern='4',caas='DA/AD',conserved_pair='0:')])[0]
    assert entry.is_conserved_meta
    details=[dict(top_tip_mode='W',bottom_tip_mode=None),dict(top_tip_mode='Y',bottom_tip_mode='A')]
    assert _derive_convergent_call([SimpleNamespace(caas='WWY/AA',pair_details=details)],'caas')=='WY/A'


def test_nonfinite_weights_are_rejected_even_in_direct_pool_api():
    with pytest.raises(ValueError,match='PSS'):
        pool_domains([dict(hyp='H1',sides={})],{('H1',1):('A','B',float('nan'))})


def test_full_ids_manifests_formats_and_ambiguity(inputs):
    assert find_gene_alignment(inputs/'alignments','AC010616.1').name=='AC010616.1.fa'
    from src.data.loaders import load_ensembl_genes
    assert load_ensembl_genes(put(inputs/'genes.tsv','gene\nAC010616.1\nAC010616.2\n'))=={'AC010616.1','AC010616.2'}
    clustal=put(inputs/'ali/G.aln','CLUSTAL W\n\nA AAAA\nB AAAA\n')
    assert read_alignment(clustal).get_alignment_length()==4
    m=put(inputs/'manifest.tsv','custom.id\tAC010616.1.fa\n')
    assert read_alignment_manifest(m,inputs/'alignments')['custom.id'].is_file()
    shutil.copy(inputs/'alignments/AC010616.1.fa',inputs/'alignments/AC010616.1.phy')
    from src.utils.io_utils import _scan_alignment_dir
    _scan_alignment_dir.cache_clear()
    with pytest.raises(ValueError,match='Ambiguous'):find_gene_alignment(inputs/'alignments','AC010616.1')


def test_cache_keys_change_with_scientific_inputs_and_legacy_is_not_read(inputs,monkeypatch):
    a=asr.load_alignment_and_mappings(inputs/'alignments/AC010616.1.fa',None,'G')
    t=asr.load_and_match_tree(inputs/'tree.nwk',a,None)
    config=asr.SingleGeneASRConfig(inputs/'alignments/AC010616.1.fa',inputs/'tree.nwk',output_dir=inputs/'cache',model='poisson')
    prepared=prepare_config(config,a,t)
    from dataclasses import replace
    assert prepare_config(replace(config,model='wag'),a,t).output_dir!=prepared.output_dir
    t.tree.root.branch_length=.234
    assert prepare_config(config,a,t).output_dir!=prepared.output_dir
    t=asr.load_and_match_tree(inputs/'tree.nwk',a,None)
    legacy=put(inputs/'cache/asr_G/rst','legacy');put(legacy.parent/'rst1','legacy');put(legacy.parent/'tree_paml.nwk','legacy')
    with pytest.raises(FileNotFoundError):asr.load_precomputed_asr('G',prepared,a)
    calls=[]
    def fake(gene,cfg,*args):
        calls.append(1)
        folder=cfg.output_dir/f'asr_{gene}';folder.mkdir()
        for name in ['rst','rst1','tree_paml.nwk']:put(folder/name,'generated')
        return SimpleNamespace(posteriors_node={1:{1:{'A':1.}}})
    monkeypatch.setattr(asr,'_run_asr_pipeline',fake)
    monkeypatch.setattr(asr,'parse_paml_rst_node_level',lambda *a,**k:({1:{1:{'A':1.}}},{}))
    monkeypatch.setattr(asr,'_write_node_id_map',lambda *a:None)
    asr.run_asr_pipeline('G',prepared,True,a,t)
    asr.run_asr_pipeline('G',prepared,True,a,t)
    assert len(calls)==1 and legacy.read_text()=='legacy'
    validate_entry(prepared.output_dir/'asr_G',prepared.provenance)
    put(prepared.output_dir/'asr_G/rst','corrupt')
    with pytest.raises(ValueError):validate_entry(prepared.output_dir/'asr_G',prepared.provenance)


def cli(script,args,cwd):
    process=subprocess.Popen([sys.executable,str(LOCAL/script),*map(str,args)],cwd=cwd,text=True,
                             stdout=subprocess.PIPE,stderr=subprocess.PIPE,start_new_session=True)
    try:out,err=process.communicate(timeout=90)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid,signal.SIGTERM);process.communicate();raise
    return process.returncode,out+err


def scoring_args(p,out):
    return ['--alignment-dir',p/'alignments','--tree',p/'tree.nwk','--discovery',p/'discovery.tab',
            '--design',p/'design.tab','--output-dir',p/out,'--asr-cache-dir',p/'cache',
            '--workers','1','--asr-model','poisson','--posterior-threshold','0.1']


def contract_args(p,shards,out):
    return ['--b0-dirs',p/shards,'--design',p/'design.tab','--discovery-file',p/'discovery.tab','--output-dir',p/out]


def test_worker_failures_do_not_produce_completed_contracts(inputs):
    shutil.rmtree(inputs/'alignments');(inputs/'alignments').mkdir()
    code,output=cli('observed_b0_main.py',scoring_args(inputs,'shards'),inputs)
    assert code!=0 and 'observed genes failed' in output
    assert not (inputs/'shards/observed_scoring.json').exists()
    code,output=cli('contract_main.py',contract_args(inputs,'shards','out'),inputs)
    assert code!=0 and 'completion record' in output
    assert not (inputs/'out/ct_disambiguation').exists()


def test_empty_discovery_needs_no_asr_or_alignment(inputs):
    p=inputs/'discovery.tab';p.write_text(p.read_text().splitlines()[0]+'\n')
    shutil.rmtree(inputs/'alignments')
    code,output=cli('observed_b0_main.py',scoring_args(inputs,'shards'),inputs)
    assert code==0,output
    code,output=cli('contract_main.py',contract_args(inputs,'shards','out'),inputs)
    assert code==0,output
    assert len((inputs/'out/ct_disambiguation/caas_convergence_master.csv').read_text().splitlines())==1
    assert 'shared_residue_count' in (inputs/'out/meta_caas/meta_caas.tsv').read_text()
    assert not (inputs/'cache').exists()


@pytest.mark.skipif(not shutil.which('codeml'),reason='PAML codeml unavailable')
@pytest.mark.parametrize('model',['poisson','lg'])
def test_real_paml_scoring_cache_reuse_and_contract_coverage(inputs,model):
    for name in ['shards','second']:
        args=scoring_args(inputs,name);args[args.index('--asr-model')+1]=model
        code,output=cli('observed_b0_main.py',args,inputs)
        assert code==0,output
    import gzip
    first=gzip.open(inputs/'shards/AC010616.1.master.csv.gz','rt').read()
    second=gzip.open(inputs/'second/AC010616.1.master.csv.gz','rt').read()
    assert first==second
    rows=list(csv.DictReader(first.splitlines()))
    assert rows and all(r['gene']=='AC010616.1' and math.isfinite(float(r['asr_path_score'])) for r in rows)
    assert len(list((inputs/'cache').rglob('provenance.json')))==1
    from src.core.driver import load_gene_context, score_labelings, pool_labelings
    ctx=load_gene_context('AC010616.1',str(inputs/'alignments'),str(inputs/'tree.nwk'),None,model,str(inputs/'cache'),.1)
    design=read_trait_pairs(inputs/'design.tab')
    labeling=(['A','C','E','G'],['B','D','F','H'])
    entry=observed_entries('AC010616.1',[dict(position='0',trait='H1',caas='DDDD/AAAA',pattern='1')])
    reduced=pool_labelings(score_labelings(ctx,'AC010616.1',['b_0'],{'b_0':labeling},{'b_0':entry},.1),None,{'b_0':labeling})[0][1]
    assert {r.side:r.asr_path_score for r in reduced}==pytest.approx({r['side']:float(r['asr_path_score']) for r in rows})
    code,output=cli('contract_main.py',contract_args(inputs,'shards','out'),inputs)
    assert code==0,output
    (inputs/'second/AC010616.1.master.csv.gz').unlink()
    code,output=cli('contract_main.py',contract_args(inputs,'second','incomplete'),inputs)
    assert code!=0 and 'coverage differs' in output


def test_incomplete_asr_is_not_published(inputs,monkeypatch):
    a=asr.load_alignment_and_mappings(inputs/'alignments/AC010616.1.fa',None,'G')
    t=asr.load_and_match_tree(inputs/'tree.nwk',a,None)
    cfg=prepare_config(asr.SingleGeneASRConfig(inputs/'alignments/AC010616.1.fa',inputs/'tree.nwk',
                                              output_dir=inputs/'cache',model='poisson'),a,t)
    class FakeReconstructor:
        def __init__(self,*args):pass
        def reconstruct_gene(self,gene,alignment,tree,output_dir):
            d=output_dir/f'asr_{gene}'
            for name in ['rst','rst1','tree_paml.nwk']:put(d/name,'partial')
            return d/'rst'
    monkeypatch.setattr(asr,'ASRReconstructor',FakeReconstructor)
    monkeypatch.setattr(asr,'parse_paml_rst_node_level',lambda *a,**k:({1:{1:{'A':1.}}},{}))
    with pytest.raises(RuntimeError,match='diagnostic artifacts'):
        asr.run_asr_pipeline('G',cfg,False,a,t)
    assert not (cfg.output_dir/'asr_G/provenance.json').exists()
    assert list(cfg.output_dir.glob('.failed-*/failure.json'))


def test_old_tiles_cannot_be_relabelled_as_new_scores(tmp_path):
    from src.utils.null_tables import write_tables,concat_tables,read_units
    from src.core.scores import SCORE_RULE
    write_tables(tmp_path/'tile',*concat_tables([]))
    assert json.loads((tmp_path/'tile/meta.json').read_text())['score_rule']==SCORE_RULE
    assert read_units(tmp_path/'tile').num_rows==0
    (tmp_path/'tile/meta.json').write_text(json.dumps({'score_rule':'caas_pos_species_pss'}))
    with pytest.raises(ValueError,match='regenerate observed and null'):
        read_units(tmp_path/'tile')


def test_discovery_merge_rejects_mixed_headers(tmp_path):
    from src.core.contract import write_discovery
    a=put(tmp_path/'a.tsv','gene\tposition\tcaas\tpattern\nG\t0\tDD/AA\t1\n')
    b=put(tmp_path/'b.tsv','gene\tposition\tcaas\tpattern\tshared_residue_count\nH\t0\tDD/AA\t1\t0\n')
    with pytest.raises(ValueError,match='schemas differ'):write_discovery([a,b],tmp_path/'merged.tsv')


def test_shell_command_prefix_preserves_quoted_paths(tmp_path):
    import shlex
    entry=put(tmp_path/"entry point's script",'#!/bin/bash\nexec "$@"\n');entry.chmod(0o755)
    program=put(tmp_path/'fake ct','#!/bin/bash\nexit 0\n');program.chmod(0o755)
    manifest=put(tmp_path/'manifest.tsv','G\tG.fa\n')
    code=subprocess.run(['bash',str(ROOT/'subworkflows/CT/local/scripts/run_ct_perm_replay_batch.sh'),
                         '--batch-id','quoted','--manifest',str(manifest),'--caas-config','unused',
                         '--resampled-path','unused','--workers','1','--ali-format','fasta',
                         '--ct-bin',shlex.join([str(entry),str(program)])],cwd=tmp_path,capture_output=True,text=True,timeout=15)
    assert code.returncode==0,code.stdout+code.stderr


def test_saved_null_rule_is_checked_before_loading_data(tmp_path):
    counts=put(tmp_path/'counts.json',json.dumps({'score_rule':'caas_pos_species_pss','cycles':1,'cycles_present':1}))
    r=subprocess.run([sys.executable,str(ROOT/'subworkflows/SCORING/local/src/null_position_stats.py'),
                      'universe','--null-positions','MISSING','--postproc','MISSING','--roster','MISSING',
                      '--counts-json',str(counts),'--output-dir',str(tmp_path/'out')],text=True,capture_output=True)
    assert r.returncode!=0 and 'Regenerate observed and null artifacts together' in r.stderr
