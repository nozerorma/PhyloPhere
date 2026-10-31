"""FADE scientific and data contracts, independent of workflow scheduling."""
import copy
import json
from pathlib import Path
import subprocess
import sys
import pytest
ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'subworkflows/FADE/local/src'
sys.path.insert(0,str(SRC));sys.path.insert(0,str(ROOT/'subworkflows/SELECTION/local/src'))
from annotate_tree_fg import annotate
from fade_json_to_tables import AA_LETTERS, GENES_COLUMNS, SITES_COLUMNS, reduce_json, write_tables
from fade_merge_tables import validate


def result(n=2,bf=150):
    j={'input':{'number of sites':n},'MLE':{'headers':[['rate',''],['bias',''],['Prob[bias>0]',''],['BayesFactor[bias>0]','']],
       'content':{aa:{'0':[[1,1,.1,1] for _ in range(n)]} for aa in AA_LETTERS}}}
    j['MLE']['content']['A']['0'][0][3]=bf
    return j


def test_exact_annotation_and_collision_free_internal_labels(tmp_path):
    (tmp_path/'tree').write_text('((A:1,X.A:1):1,(IntNode1:1,B:1):1);')
    (tmp_path/'fa').write_text('>A\nACDE\n>X.A\nACDE\n>IntNode1\nACDE\n>B\nACDE\n')
    (tmp_path/'fg').write_text('A\n')
    output=tmp_path/'annotated'
    a=annotate(tmp_path/'tree',tmp_path/'fa',tmp_path/'fg',output)
    text=output.read_text()
    assert 'A{Foreground}:' in text and 'X.A{Foreground}' not in text
    assert 'IntNode1{Foreground}' not in text and '_IntNode1' in text
    assert a['foreground_internal_branches']==0


@pytest.mark.parametrize('mode,count',[('none',2),('all_descendants',3)])
def test_internal_clade_stem_is_labelled(tmp_path,mode,count):
    (tmp_path/'tree').write_text('((A:1,B:1):1,(C:1,D:1):1);')
    (tmp_path/'fa').write_text(''.join(f'>{s}\nACDE\n' for s in 'ABCD'))
    (tmp_path/'fg').write_text('A\nB\n')
    out=tmp_path/'annotated';annotate(tmp_path/'tree',tmp_path/'fa',tmp_path/'fg',out,internal_mode=mode)
    assert out.read_text().count('{Foreground}')==count


def test_infinite_bf_and_zero_based_sites_are_consistent():
    row,sites=reduce_json(result(bf=float('inf')),'NA','top',100,2)
    assert row[3:7]==[1,float('inf'),'A',1]
    assert sites==[['top','NA',1,0,float('inf'),'A']]


@pytest.mark.parametrize('kind',['negative','nan','partial','unequal','partitions','header','coverage','shortrow'])
def test_incomplete_or_invalid_json_fails(kind):
    j=result()
    if kind=='negative':j['MLE']['content']['A']['0'][0][3]=-1
    if kind=='nan':j['MLE']['content']['A']['0'][0][3]=float('nan')
    if kind=='partial':del j['MLE']['content']['Y']
    if kind=='unequal':j['MLE']['content']['Y']['0'].pop()
    if kind=='partitions':j['MLE']['content']['A']['1']=[[1,1,.1,20]]
    if kind=='header':j['MLE']['headers'][3][0]='wrong'
    if kind=='coverage':j['input']['number of sites']=3
    if kind=='shortrow':j['MLE']['content']['A']['0'][0].pop()
    with pytest.raises(ValueError):reduce_json(j,'g','top',100,2)


def validated_rows():
    row,sites=reduce_json(result(),'001','top',100,2)
    genes=[dict(zip(GENES_COLUMNS,map(str,row)))];sites=[dict(zip(SITES_COLUMNS,map(str,r))) for r in sites]
    audit=[dict(gene='001',direction='top',status='complete',alignment_columns=2,n_significant_sites=1,retained_taxa=['A','B'],foreground_taxa=['A'])]
    return genes,sites,audit


@pytest.mark.parametrize('kind',['duplicate_gene','duplicate_site','missing_audit','missing_gene','orphan','counts','coordinate','top'])
def test_merge_rejects_inconsistent_tables(kind):
    g,s,a=validated_rows()
    if kind=='duplicate_gene':g+=copy.deepcopy(g)
    if kind=='duplicate_site':s+=copy.deepcopy(s)
    if kind=='missing_audit':a=[]
    if kind=='missing_gene':g=[]
    if kind=='orphan':s[0]['gene']='orphan'
    if kind=='counts':g[0]['n_sig_site_aa']='2'
    if kind=='coordinate':s[0]['position']='0.5'
    if kind=='top':s[0]['max_bf']='149'
    with pytest.raises(ValueError):validate(g,s,a,[['001','top']],100)


def test_literal_gene_lists_and_valid_empty(tmp_path):
    from fade_merge_tables import write,SUMMARY
    for genes in [[],[dict(gene=g,n_sites_tested=20,n_sig_site_aa=1,max_bf=150,top_target_aa='A',top_site=1) for g in ['001','NA','g.v2']]]:
        write(tmp_path/'summary',SUMMARY,genes)
        subprocess.run([sys.executable,'-B',str(SRC/'fade_gene_lists.py'),'--summary','summary','--direction','top'],cwd=tmp_path,check=True)
        assert (tmp_path/'background.txt').read_text().splitlines()==[g['gene'] for g in genes]
        assert (tmp_path/'fade_top_significant.txt').read_text().splitlines()==[g['gene'] for g in genes]


def test_valid_and_invalid_biological_exclusion_contract():
    audit=dict(gene='g',direction='top',status='excluded',alignment_columns=20,
               retained_taxa=['B','C'],foreground_taxa=[],reason='no_foreground')
    validate([],[],[audit],[['g','top']],100)
    audit['reason']='technical_failure'
    with pytest.raises(ValueError):validate([],[],[audit],[['g','top']],100)


def test_opposite_background_preserves_columns_and_root_requirement(tmp_path):
    (tmp_path/'tree').write_text('((A:1,B:1):1,(C:1,D:1):1);')
    (tmp_path/'fa').write_text(''.join(f'>{s}\nAC-DE\n' for s in 'ABCD'))
    (tmp_path/'fg').write_text('A\n');(tmp_path/'bg').write_text('C\n')
    result=annotate(tmp_path/'tree',tmp_path/'fa',tmp_path/'fg',tmp_path/'annotated',
                    bg_path=tmp_path/'bg',fasta_out=tmp_path/'filtered')
    assert result['retained_taxa']==['A','C'] and result['alignment_columns']==5
    assert (tmp_path/'filtered').read_text()=='>A\nAC-DE\n>C\nAC-DE\n'
    (tmp_path/'tree').write_text('[&U] ((A:1,B:1):1,(C:1,D:1):1);')
    with pytest.raises(ValueError,match='rooted'):annotate(tmp_path/'tree',tmp_path/'fa',tmp_path/'fg',tmp_path/'annotated')
