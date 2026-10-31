"""CT_POSTPROC quality universe, strict contracts, null symmetry and real workflow."""
import csv
import hashlib
import importlib.util
import json
import logging
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
from types import SimpleNamespace

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'subworkflows/CT_POSTPROC/local/src'
sys.path.insert(0, str(ROOT / 'subworkflows/CT_DISAMBIGUATION/local'))
sys.path.insert(0, str(SRC))
from postproc_io import load_observed, coverage, parameter_token
from residue_descriptors import add_species_tally, tally_position
from cleanup_background import load_global_background_genes
from src.core.postproc import ctrain, load_gene_lengths
from src.core.columns import index_files, gene_columns
from src.utils.null_tables import CycleReducer
spec = importlib.util.spec_from_file_location('ct_cluster', SRC / 'filter_caas_clusters-param.py')
cluster = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cluster)


def observed(directory, rows=None):
    directory.mkdir(parents=True, exist_ok=True)
    rows = rows if rows is not None else [
        ['001.2',0,'top',.8],['001.2',0,'bottom',.6],
        ['001.2',1,'top',.7],['001.2',2,'none',0],['other',4,'top',.9]]
    df = pd.DataFrame(rows, columns=['gene','msa_pos','side','asr_path_score'])
    master = directory / 'caas_convergence_master.csv'
    df.to_csv(master,index=False)
    genes = coverage(df)
    batch = dict(version=1,complete=True,score_rule='caas_pos_species_pss_v4',genes=genes,
                 excluded_genes=[],model='lg',posterior_threshold=.2,tree_sha256='tree',
                 weighted=False,taxid_sha256=None,weights_sha256=None)
    stamp = directory / 'observed_scoring.json'
    stamp.write_text(json.dumps(dict(version=1,complete=True,genes=genes,excluded_genes=[],batches=[batch])))
    return master, stamp


@pytest.mark.parametrize('change', ['fractional','negative','nan','duplicate','conflict','wrong_side','positive_none','zero_top','coverage','incomplete','legacy','batch','settings'])
def test_invalid_observed(tmp_path,change):
    master, stamp = observed(tmp_path)
    df=pd.read_csv(master,keep_default_na=False,dtype=str).astype(object)
    rec=json.loads(stamp.read_text())
    if change in ['fractional','negative']: df.loc[0,'msa_pos']=1.9 if change=='fractional' else -1
    elif change=='nan': df.loc[0,'asr_path_score']='nan'
    elif change=='duplicate': df=pd.concat([df,df.iloc[[0]]])
    elif change=='conflict': df.loc[0,'side']='none'; df.loc[0,'asr_path_score']=0
    elif change=='wrong_side': df.loc[0,'side']='both'
    elif change=='positive_none': df.loc[3,'asr_path_score']=.2
    elif change=='zero_top': df.loc[0,'asr_path_score']=0
    elif change=='coverage': rec['genes']={}
    elif change=='incomplete': rec['complete']=False
    elif change=='legacy': rec['batches'][0]['score_rule']='old'
    elif change=='batch': rec['batches']*=2
    elif change=='settings': rec['batches'][0]['posterior_threshold']=float('inf')
    df.to_csv(master,index=False);stamp.write_text(json.dumps(rec))
    with pytest.raises((ValueError,TypeError)): load_observed(master,stamp)


def test_preserve_identifiers_scores_and_diagnostic_universe(tmp_path):
    master,stamp=observed(tmp_path)
    result=subprocess.run([sys.executable,str(SRC/'prepare_postproc_input.py'),'--input',str(master),
        '--completion',str(stamp),'--output','prepared.tsv','--removed-output','removed.tsv'],cwd=tmp_path,capture_output=True,text=True)
    assert result.returncode==0,result.stderr
    prepared=pd.read_csv(tmp_path/'prepared.tsv',sep='\t',keep_default_na=False,dtype={'Gene':str})
    assert prepared.Gene.iloc[0]=='001.2' and len(prepared)==5
    flags=cluster.filterCAAS(tmp_path/'prepared.tsv',.7,3,logging.getLogger('test'))
    flags=pd.read_csv(flags,sep='\t',dtype={'Gene':str})
    assert len(flags)==4 and sum(flags.clustering_flag=='Discarded')==3
    universe=json.loads((tmp_path/'quality_universe.json').read_text())
    assert universe['genes']['001.2']['positions']==3
    assert parameter_token(.701)!=parameter_token(.709)


def test_separate_minima_and_null_symmetry():
    assert ctrain([0,2],.6,3)==[]
    assert ctrain([0,2],.6,3,min_positions=2)==[0,2]
    assert ctrain([0,1],.6,3,min_positions=2)==[]
    # Every qualifying interval needs the minimum count, even inside a larger gene.
    assert ctrain([0,2,100],.6,3)==[0,2]  # Legacy minimum was checked per gene.
    assert ctrain([0,2,100],.6,3,min_positions=3)==[]
    rows=[SimpleNamespace(position=p,side=s,asr_path_score=v) for p,s,v in [(0,'top',.8),(1,'bottom',.7),(2,'none',0)]]
    result=CycleReducer(True,True,3,.7)("b_1",rows)
    assert result.n_pos==3 and result.clust and result.positions==[]
    rows=[SimpleNamespace(position=p,side='top',asr_path_score=.8) for p in [0,2]]
    assert CycleReducer(True,True,3,.6,clust_min_positions=2)('b_1',rows).positions==[]


def test_train_matches_exhaustive_rule():
    rng=random.Random(42)
    for _ in range(200):
        ps=sorted(rng.sample(range(30),rng.randrange(1,12)))
        span=rng.randrange(1,8); count=rng.randrange(1,6); density=rng.choice([0,.3,.6,.7,1])
        expected=set()
        for i in range(len(ps)):
            for j in range(i,len(ps)):
                width=ps[j]-ps[i]+1
                if j-i+1>=count and width>=span and (j-i+1)/width>=density: expected.update(ps[i:j+1])
        assert ctrain(ps,density,span,min_positions=count)==sorted(expected)



def test_legacy_train_matches_previous_implementation():
    rng=random.Random(77)
    for _ in range(200):
        ps=sorted(rng.sample(range(30),rng.randrange(1,12)))
        minimum=rng.randrange(1,10); density=rng.choice([0,.3,.6,.7,1])
        expected=set()
        if len(ps)>=minimum:
            for r in range(len(ps)):
                for l in range(r+1):
                    span=ps[r]-ps[l]+1
                    if span>=minimum and (r-l+1)/span>=density: expected.update(ps[l:r+1])
        assert ctrain(ps,density,minimum)==sorted(expected)


def test_dotted_maps_lengths_and_ambiguity(tmp_path):
    maps=tmp_path/'maps';(maps/'nested').mkdir(parents=True)
    path=maps/'nested/001.2.human.map.tsv'
    path.write_text('ori_codon_col\tstatus\tprot_ali_col\n1\tselected\t1\n2\tremoved\t\n3\tselected\t2\n')
    lengths=tmp_path/'lengths.tsv';lengths.write_text('gene\tlength\n001.2\t99\nmissing\t\n')
    assert gene_columns(index_files(maps,'.map.tsv',['001.2']),'001.2','.map.tsv')=={0:1,1:3}
    # The null merge resolves with annotation genes; observed resolves with candidate genes.
    assert load_gene_lengths(str(lengths),str(maps))=={'001.2':3.0}
    assert load_gene_lengths(str(lengths),str(maps),genes=['001.2'])=={'001.2':3.0}
    (maps/'001.2.other.map.tsv').write_text(path.read_text())
    with pytest.raises(ValueError,match='several'):load_gene_lengths(str(lengths),str(maps))


def test_cluster_coverage_is_required(tmp_path):
    from postproc_io import validate_cluster
    discovery=pd.DataFrame({'Gene':['a','a'],'Position':[0,1]})
    incomplete=pd.DataFrame({'Gene':['a'],'Position':[0],'clustering_flag':['Discarded']})
    with pytest.raises(ValueError,match='coverage'):validate_cluster(discovery,incomplete)


def pairs(path):
    path.write_text('hypothesis_id\tpair\tspecies1\tspecies2\tpss_score\nH1\t1\ta\tb\t1\n')
    return path


def test_nested_dotted_alignment_and_ambiguous_residues(tmp_path):
    ali=tmp_path/'ali';(ali/'nested').mkdir(parents=True)
    (ali/'nested/001.2.fasta').write_text('>a\nAX-\n>b\nAXA\n')
    df=pd.DataFrame({'Gene':['001.2']*3,'Position':[0,1,2],'participating_hypotheses':['H1']*3})
    out=add_species_tally(df,ali,pairs(tmp_path/'pairs.tsv'))
    assert out.n_conserved_pairs.tolist()==['1/1','','']
    assert out.n_top_species.tolist()==['1','0','0']
    assert out.n_bottom_species.tolist()==['1','0','1']
    df.loc[0,'participating_hypotheses']='H2'
    with pytest.raises(ValueError,match='Unknown'):add_species_tally(df,ali,tmp_path/'pairs.tsv')
    df.loc[0,'participating_hypotheses']='H1'
    (ali/'001.2.fasta').write_text('>a\nA\n>b\nA\n')
    with pytest.raises(ValueError,match='Ambiguous'):add_species_tally(df,ali,tmp_path/'pairs.tsv')


def test_literal_background_and_empty(tmp_path):
    bg=tmp_path/'bg';bg.write_text('NA\n001.2\n')
    assert load_global_background_genes(bg)==['001.2','NA']
    bg.write_text('');assert load_global_background_genes(bg)==[]


def run_cmd(args,cwd):
    result=subprocess.run(args,cwd=cwd,capture_output=True,text=True,timeout=180)
    assert result.returncode==0,result.stdout+result.stderr
    return result


@pytest.fixture
def project(tmp_path):
    project=tmp_path/"project's files";project.mkdir()
    for directory in ['subworkflows/CT_POSTPROC','subworkflows/CT_DISAMBIGUATION/local',
                      'subworkflows/ASR_DIAGNOSTICS','subworkflows/CT/local/modules']:
        shutil.copytree(ROOT/directory,project/directory,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    for rel in ['workflows/ct_postproc.nf','workflows/asr_diagnostics.nf','conf/ct_postproc.config',
                'subworkflows/TRAIT_ANALYSIS/local/src/report_style.R']:
        dest=project/rel;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy(ROOT/rel,dest)
    master,stamp=observed(project/'observed')
    (project/'lengths.tsv').write_text('gene\tlength\tchr\tstart\tend\n001.2\t100\t1\t1\t300\nother\t80\t1\t400\t640\n')
    (project/'background.txt').write_text('001.2\nother\nNA\n')
    (project/'main.nf').write_text("""
include { CT_POSTPROC } from './workflows/ct_postproc.nf'
workflow { CT_POSTPROC(params.live ? Channel.value(file(params.disambiguation_input)) : null, null, null, null) }
""")
    (project/'audit.config').write_text("""
includeConfig 'conf/ct_postproc.config'
process.executor='local'
process.cpus=1
params {
  disambiguation_input='observed/caas_convergence_master.csv'
  background_input='background.txt'; gene_ensembl_file='lengths.tsv'
  outdir='out'; live=false; asr_diagnostics=false; alignment=''; tax_id=''; ali_format='auto'
  use_singularity=false;use_apptainer=false;seed=1998
}
""")
    return project


def run_nf(project,*extra):
    return run_cmd(['nextflow','-log','nextflow.log','run','main.nf','-c','audit.config','-w','work','-ansi-log','false',*extra],project)


def tasks(project):
    return [p for p in (project/'work').glob('*/*/.command.sh') if "rmarkdown::render('8.Characterization_report.Rmd'" in p.read_text()]


def test_workflow_filter_exploration_empty_and_resume(project):
    run_nf(project,'--live','true')
    out=project/'out'
    assert (out/'html_reports/8.Characterization_report.html').exists()
    kept=pd.read_csv(out/'postproc/gene_filtering/filtered_discovery.tsv',sep='\t',dtype={'Gene':str})
    assert kept.Gene.tolist()==['other']
    assert 'NA' in (out/'postproc/cleaned_backgrounds/cleaned_background_main.txt').read_text().splitlines()
    assert len(tasks(project))==1
    run_nf(project,'--live','true','-resume');assert len(tasks(project))==1
    # Rmd/helper/style changes invalidate report cache.
    for i,rel in enumerate(['subworkflows/CT_POSTPROC/local/8.Characterization_report.Rmd',
                           'subworkflows/CT_POSTPROC/local/src/report_io.R',
                           'subworkflows/TRAIT_ANALYSIS/local/src/report_style.R'],2):
        p=project/rel;p.write_text(p.read_text()+ ('\n<!-- audit -->\n' if p.suffix=='.Rmd' else '\n# audit\n'))
        run_nf(project,'--live','true','-resume');assert len(tasks(project))==i
    p=project/'observed/caas_convergence_master.csv';stat=p.stat();p.write_text(p.read_text().replace('0.9','0.5'))
    os.utime(p,ns=(stat.st_atime_ns,stat.st_mtime_ns))
    run_nf(project,'-resume');assert len(tasks(project))==5
    run_nf(project,'--caas_postproc_mode','exploratory','--minlen_values','2,3','--maxcaas_values','0.701,0.709','--remove_caas_clusters','false','-resume')
    settings=json.loads((out/'postproc/postproc_inputs/postproc_settings.json').read_text())
    assert settings['remove_caas_clusters'] is False
    kept=pd.read_csv(out/'postproc/gene_filtering/filtered_discovery.tsv',sep='\t');assert len(kept)==5
    summary=pd.read_csv(out/'postproc/summary_statistics/filter_summary.tsv',sep='\t');assert len(summary)==5
    assert max(summary.DiscardedCount)<=4
    run_nf(project,'--filter_minlen','3','--filter_min_positions','2','--filter_maxcaas','0.6','-resume')
    settings=json.loads((out/'postproc/postproc_inputs/postproc_settings.json').read_text())
    assert settings['selected_filter_min_positions']==2 and settings['legacy_gene_count_gate'] is False
    observed(project/'observed',rows=[])
    run_nf(project,'--gene_filter_mode','stratified','-resume')
    assert pd.read_csv(out/'postproc/gene_filtering/filtered_discovery.tsv',sep='\t').empty
    assert (out/'html_reports/8.Characterization_report.html').exists()


def test_directory_fingerprints_and_tally_cache(project):
    ali=project/'alignments';(ali/'nested').mkdir(parents=True)
    (ali/'nested/001.2.fasta').write_text('>a\nAAA\n>b\nAAA\n')
    (ali/'other.fasta').write_text('>a\nAAAAA\n>b\nAAAAA\n')
    maps=project/'maps';(maps/'nested').mkdir(parents=True)
    map_file=maps/'nested/001.2.ref.map.tsv'
    map_file.write_text('ori_codon_col\tstatus\tprot_ali_col\n1\tselected\t1\n2\tselected\t2\n3\tselected\t3\n4\tremoved\t\n')
    hp=pairs(project/'pairs.tsv')
    (project/'main.nf').write_text("""
include { CT_POSTPROC } from './workflows/ct_postproc.nf'
workflow { CT_POSTPROC(null,null,null,Channel.value(file('pairs.tsv'))) }
""")
    run_nf(project,'--alignment','alignments','--caas_map_dir','maps')
    initial=(project/'out/postproc/preprocessed/postproc_disambiguation_input.tsv').read_text()
    assert 'A:1' in initial
    source=ali/'nested/001.2.fasta';stat=source.stat()
    source.write_text(source.read_text().replace('AAA','CAA',1));os.utime(source,ns=(stat.st_atime_ns,stat.st_mtime_ns))
    run_nf(project,'--alignment','alignments','--caas_map_dir','maps','-resume')
    after=(project/'out/postproc/preprocessed/postproc_disambiguation_input.tsv').read_text()
    assert initial!=after and 'C:1' in after
    stat=map_file.stat()
    map_file.write_text(map_file.read_text().replace('3\tselected\t3\n4\tremoved\t','3\tremoved\t\n4\tselected\t3'))
    assert map_file.stat().st_size==stat.st_size
    os.utime(map_file,ns=(stat.st_atime_ns,stat.st_mtime_ns))
    run_nf(project,'--alignment','alignments','--caas_map_dir','maps','-resume')
    flags=pd.read_csv(project/'out/postproc/filter_selected/postproc_disambiguation_input.filtered.minlen3.maxcaas70.tsv',sep='\t')
    # Span changed from 3 to 4, density .75 still qualifies; task must rerun despite unchanged output.
    cmds=[p for p in (project/'work').glob('*/*/.command.sh') if 'filter_caas_clusters-param.py' in p.read_text()]
    assert len(cmds)==3
    # Completion contents are hashed separately from the master.
    stamp=project/'observed/observed_scoring.json';stat=stamp.stat()
    stamp.write_text(stamp.read_text().replace('"posterior_threshold": 0.2','"posterior_threshold": 0.3'))
    os.utime(stamp,ns=(stat.st_atime_ns,stat.st_mtime_ns))
    run_nf(project,'--alignment','alignments','--caas_map_dir','maps','-resume')
    cmds=[p for p in (project/'work').glob('*/*/.command.sh') if 'prepare_postproc_input.py' in p.read_text()]
    assert len(cmds)==3


def test_report_regeneration_uses_saved_contract(project):
    # A complete published snapshot is prepared by a real workflow.
    run_nf(project)
    from gui.generation.report_registry import detect_reports, build_script
    report={r.id:r for r in detect_reports(str(project/'out'),ROOT)}['postproc']
    assert report.runnable,report.missing_required
    script=build_script(report,ROOT,use_singularity=False)
    assert 'postproc_settings.json' in script and 'report_style.R' in script
    syntax=subprocess.run(['bash','-n'],input=script,text=True,capture_output=True)
    assert syntax.returncode==0,syntax.stderr
    result=subprocess.run(['bash'],input=script,text=True,capture_output=True,timeout=120)
    assert result.returncode==0,result.stdout+result.stderr
    assert (project/'out/html_reports/8.Characterization_report.html.bak').exists()
    (project/'out/postproc/postproc_inputs/postproc_settings.json').unlink()
    report={r.id:r for r in detect_reports(str(project/'out'),ROOT)}['postproc']
    assert not report.runnable


@pytest.mark.parametrize('kind',['missing_background','ambiguous_background','incomplete_observed'])
def test_workflow_rejects_invalid_inputs(project,kind):
    extra=[]
    if kind.endswith('background'):
        bg=project/'bgdir';bg.mkdir()
        if kind=='ambiguous_background':
            (bg/'a.background_genes').write_text('001.2\n')
            (bg/'b.background_genes').write_text('other\n')
        extra=['--background_input','bgdir']
    else:
        stamp=project/'observed/observed_scoring.json'
        rec=json.loads(stamp.read_text());rec['complete']=False;stamp.write_text(json.dumps(rec))
    result=subprocess.run(['nextflow','-log','nextflow.log','run','main.nf','-c','audit.config',
        '-w','work','-ansi-log','false',*extra],cwd=project,capture_output=True,text=True,timeout=120)
    assert result.returncode!=0
    assert not (project/'out/html_reports/8.Characterization_report.html').exists()
    if kind.endswith('background'):assert 'exactly one background_genes' in result.stdout+result.stderr
    else:assert 'Invalid observed completion' in result.stdout+result.stderr


def test_optional_taxon_mapping_and_invalid_supplied_design(tmp_path):
    ali=tmp_path/'ali';ali.mkdir();(ali/'g.fasta').write_text('>11\nA\n>22\nA\n')
    mapping=tmp_path/'taxid.tsv';mapping.write_text('species\ttax_id\na\t11\nb\t22\n')
    df=pd.DataFrame({'Gene':['g'],'Position':[0],'participating_hypotheses':['H1']})
    design=pairs(tmp_path/'pairs.tsv')
    result=add_species_tally(df,ali,design,taxid_mapping=mapping)
    assert result.n_conserved_pairs.tolist()==['1/1']
    design.write_text(design.read_text().replace('H1\t1','H1\t2'))
    with pytest.raises(ValueError,match='contiguous'):add_species_tally(df,ali,design,taxid_mapping=mapping)


def test_full_pipeline_help_compiles_modified_null_connection(tmp_path):
    # Compile the production imports without scheduling scientific tasks.
    result=subprocess.run(['nextflow','-log','help.log','run',str(ROOT/'main.nf'),
                   '--help','--ct_postproc','true','-ansi-log','false'],cwd=tmp_path,capture_output=True,text=True,timeout=120)
    # HELP deliberately calls System.exit(1); its exit status is not a compile failure.
    assert result.returncode==1
    assert 'filter_min_positions' in result.stdout+result.stderr
    assert 'ERROR' not in result.stdout+result.stderr




def test_length_source_units_are_explicit(tmp_path):
    annotation=tmp_path/'genes.tsv'
    annotation.write_text('gene\tlength\tstart\tend\ngenomic\t100\t101\t200\nunknown\t50\t\t\n')
    sources={}
    assert load_gene_lengths(str(annotation),sources=sources)=={'genomic':100,'unknown':50}
    assert sources=={'genomic':'genomic_bp','unknown':'annotation_units'}
