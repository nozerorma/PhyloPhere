"""RER scripts run serially in isolated directories, using the installed RERconverge."""
import json
from pathlib import Path
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]
LOCAL = ROOT / 'subworkflows/RERCONVERGE/local'


def rrun(wd, script, *args, success=True):
    result = subprocess.run(['Rscript', str(script), *map(str,args)], cwd=wd,
                            text=True, capture_output=True, timeout=120)
    if success:
        assert result.returncode == 0, result.stdout + result.stderr
    else:
        assert result.returncode != 0, result.stdout + result.stderr
    return result.stdout + result.stderr


def rexpr(wd, expr):
    return rrun(wd, '-e', expr)


@pytest.fixture
def inputs(tmp_path):
    (tmp_path/'traits.tsv').write_text('species\ttrait\tn\tc\nA\t2\t10\t2\nB\t4\t10\t4\nC\t4\t10\t4\nD\t2\t10\t2\n')
    (tmp_path/'trees.nwk').write_text('g.v1\t((A:1,B:1):1,C:1,D:1);\n')
    return tmp_path


def build(wd, mode='auto'):
    return rrun(wd, LOCAL/'build_rer_trait.R', 'traits.tsv','species','trait','trait.RData','type.txt','n','c','none',mode)


def test_r_contracts():
    rrun(ROOT, ROOT/'tests/test_rer.R')
    rrun(ROOT, ROOT/'tests/test_wp6_rer_batches.R')


def test_trait_mode_and_counts(inputs):
    build(inputs,'continuous')
    rexpr(inputs,"load('trait.RData');stopifnot(identical(unname(trait_vector),c(2,4,4,2)),length(n_vector)==4)")
    build(inputs,'binary')
    rexpr(inputs,"load('trait.RData');stopifnot(identical(unname(trait_vector),c(0,1,1,0)),length(n_vector)==4,identical(unname(raw_trait_vector),c(2,4,4,2)))")
    rrun(inputs,LOCAL/'rer_trait_mode.R','trait.RData','auto','mode.txt')
    assert (inputs/'mode.txt').read_text().strip()=='binary'
    assert 'conflicts' in rrun(inputs,LOCAL/'rer_trait_mode.R','trait.RData','continuous','mode.txt',success=False)


@pytest.mark.parametrize('text,message',[
    ('species\ttrait\nA\t1\nA\t2\n','duplicate'),
    ('species\ttrait\nA A\t1\nA_A\t2\n','duplicate'),
    ('species\ttrait\nA\tfoo\nB\t2\n','finite'),
    ('species\ttrait\nA\tInf\nB\t2\n','finite'),
    ('species\ttrait\tn\tc\nA\t1\t10\t11\nB\t2\t10\t1\n','Counts require'),
])
def test_invalid_traits(inputs,text,message):
    (inputs/'traits.tsv').write_text(text)
    log=rrun(inputs,LOCAL/'build_rer_trait.R','traits.tsv','species','trait','trait.RData','type.txt',
             'n' if '\tn\tc' in text else '', 'c' if '\tn\tc' in text else '', 'none','auto',success=False)
    assert message in log


def tree(wd,success=True):
    return rrun(wd,LOCAL/'rer_master_tree.R','trees.nwk','traits.tsv','species','pruned.txt','trees.rds','NO_FILE','trait',success=success)


def test_single_tree_and_duplicate_ids(inputs):
    tree(inputs)
    assert (inputs/'pruned.txt').read_text().startswith('g.v1\t')
    (inputs/'trees.nwk').write_text((inputs/'trees.nwk').read_text()*2)
    assert 'duplicate' in tree(inputs,success=False)


def test_normalized_species_and_exclusion_audit(inputs):
    (inputs/'traits.tsv').write_text('species\ttrait\nA a\t1\nB\t2\nC\t3\nD\t4\n')
    (inputs/'trees.nwk').write_text('001\t((A_a:1,B:1):1,C:1,D:1);\nNA\t(X:1,Y:1,Z:1);\n')
    tree(inputs)
    assert (inputs/'pruned.txt').read_text().startswith('001\t')
    assert 'NA\texcluded\tfewer_than_three_trait_taxa' in (inputs/'pruned.txt.audit.tsv').read_text()


def result(wd,empty=False):
    code="r<-data.frame(Rho=c(.5,-.6,NA),N=c(5,6,0),P=c(.01,0,NA),p.adj=c(.02,0,NA),p.perm=c(.01,.1,NA),row.names=c('001','NA',\"g'quote\"));"
    rexpr(wd,code+("r<-r[FALSE,,drop=FALSE];" if empty else '')+"saveRDS(r,'result.rds')")


def tables(wd,null='NO_FILE',success=True):
    return rrun(wd,LOCAL/'rer_tables.R','result.rds',null,'trait',.05,.05,success=success)


@pytest.mark.parametrize('empty',[False,True])
def test_tables_and_literal_gene_lists(tmp_path,empty):
    result(tmp_path,empty)
    tables(tmp_path)
    rrun(tmp_path,LOCAL/'rer_gene_lists.R','rerconverge_summary_trait.tsv','NO_FILE',.05)
    assert (tmp_path/'background.txt').read_text().splitlines()==([] if empty else ['001','NA',"g'quote"])
    assert (tmp_path/'rer_significant.txt').read_text().splitlines()==([] if empty else ['001'])
    meta=json.loads((tmp_path/'rer_run_metadata.json').read_text())
    assert meta['null_status']=='unavailable' and meta['inference_provenance']=='legacy_unknown'


def test_duplicate_annotation_rejected(tmp_path):
    result(tmp_path); tables(tmp_path)
    (tmp_path/'gs.tsv').write_text('gene\tscore_top\n001\t2\n001\t3\n')
    assert 'duplicate' in rrun(tmp_path,LOCAL/'rer_gene_lists.R','rerconverge_summary_trait.tsv','gs.tsv',.05,success=False)


def test_manifest_and_tampering(tmp_path):
    rexpr(tmp_path,f"source('{LOCAL}/rer_batches.R');r<-data.frame(Rho=c(.2,-.3),N=c(5,5),P=c(.2,.3),p.adj=c(.4,.4),row.names=c('001','NA'));r<-rer_write_batches(r,2,2,0,function(n)list(corRho=matrix(c(.1,.2,.3,.4),2,dimnames=list(rownames(r),NULL)),corStat=matrix(1:4,2,dimnames=list(rownames(r),NULL))),'null');saveRDS(r,'result.rds')")
    tables(tmp_path,'null')
    meta=json.loads((tmp_path/'rer_run_metadata.json').read_text()); assert meta['null_status']=='verified'
    block=tmp_path/'null/blk_000001.rds';block.write_bytes(block.read_bytes()+b'changed')
    assert 'checksum' in tables(tmp_path,'null',success=False)


@pytest.mark.parametrize('empty',[False,True])
def test_report_standalone(tmp_path,empty):
    result(tmp_path,empty); tables(tmp_path)
    for name in ['5.RERconverge_report.Rmd','rer_contract.R']:
        shutil.copy2(LOCAL/name,tmp_path/name)
    shutil.copy2(ROOT/'subworkflows/TRAIT_ANALYSIS/local/src/report_style.R',tmp_path/'report_style.R')
    rrun(tmp_path,LOCAL/'render_rer.R','validated_rer.output','trait',5)
    html=(tmp_path/'5.RERconverge_report.html').read_text()
    assert 'original inference settings and provenance are unknown' in html
    if empty: assert 'Valid empty result' in html


def test_all_unavailable_and_no_null(tmp_path):
    rexpr(tmp_path,"r<-data.frame(Rho=NA,N=0,P=NA,p.adj=NA,row.names='001');saveRDS(r,'result.rds')")
    tables(tmp_path)
    rrun(tmp_path,LOCAL/'rer_gene_lists.R','rerconverge_summary_trait.tsv','NO_FILE',.05)
    assert (tmp_path/'background.txt').read_text()=='001\n'
    assert not (tmp_path/'rer_significant.txt').read_text()


def test_literal_tree_gene_quotes(inputs):
    (inputs/'trees.nwk').write_text("g'quote\t((A:1,B:1):1,C:1,D:1);\n")
    tree(inputs)
    rexpr(inputs,"t<-readRDS('trees.rds');stopifnot(identical(names(t$trees),\"g'quote\"))")


def test_taxid_ambiguity_rejected(inputs):
    (inputs/'tax.tsv').write_text('species\ttax_id\nA\t1\nB\t1\n')
    assert 'ambiguous' in rrun(inputs,LOCAL/'rer_master_tree.R','trees.nwk','traits.tsv','species','p.txt','t.rds','tax.tsv','trait',success=False)


def test_negative_branch_lengths(inputs):
    (inputs/'trees.nwk').write_text('g\t((A:-1,B:1):1,C:1,D:1);\n')
    assert 'nonnegative' in tree(inputs,success=False)


def test_missing_null_blocks_rejected(tmp_path):
    rexpr(tmp_path,f"source('{LOCAL}/rer_batches.R');r<-data.frame(Rho=c(.2,-.3),N=c(5,5),P=c(.2,.3),p.adj=c(.4,.4),row.names=c('001','NA'));r<-rer_write_batches(r,2,2,0,function(n)list(corRho=matrix(c(.1,.2,.3,.4),2,dimnames=list(rownames(r),NULL)),corStat=matrix(1:4,2,dimnames=list(rownames(r),NULL))),'null');saveRDS(r,'result.rds')")
    (tmp_path/'null/blk_000002.rds').unlink()
    assert 'Incomplete' in tables(tmp_path,'null',success=False)


def test_optional_annotations_keep_identifiers_and_order(tmp_path):
    result(tmp_path); tables(tmp_path)
    (tmp_path/'gs.tsv').write_text('gene\tscore_top\tflag_fade\nNA\t7\tTRUE\n001\t3\tFALSE\n')
    rrun(tmp_path,LOCAL/'rer_gene_lists.R','rerconverge_summary_trait.tsv','gs.tsv',.05)
    rexpr(tmp_path,f"source('{LOCAL}/rer_contract.R');x<-rer_table('fcs_stats.tsv');stopifnot(identical(x$gene,c('001','NA',\"g'quote\")),identical(x$score_top,c('3','7','NA')))")


def test_regeneration_snapshot_checksums(tmp_path):
    result(tmp_path);tables(tmp_path)
    rexpr(tmp_path,f"source('{LOCAL}/rer_contract.R');rer_verify_snapshot('rer_inputs/validated_rer.output')")
    tsv=tmp_path/'rer_inputs/rerconverge_summary_trait.tsv';tsv.write_text(tsv.read_text()+'changed')
    assert 'checksum mismatch' in rrun(tmp_path,'-e',f"source('{LOCAL}/rer_contract.R');rer_verify_snapshot('rer_inputs/validated_rer.output')",success=False)


def test_gui_regeneration_uses_snapshot_and_binary_fallback(tmp_path):
    from gui.generation.report_registry import detect_reports, build_script

    result(tmp_path)
    rrun(tmp_path,LOCAL/'rer_tables.R','result.rds','NO_FILE','trait',.02,.03)
    out=tmp_path/'out'
    html=out/'html_reports/5.RERconverge_report.html';html.parent.mkdir(parents=True);html.write_text('fixture')
    results=out/'rerconverge/rer_results';results.mkdir(parents=True)
    shutil.copy2(tmp_path/'result.rds',results/'trait.binary.output')
    reports={r.id:r for r in detect_reports(str(out),ROOT)}
    assert reports['rer'].runnable
    assert str(reports['rer'].slots[0].path).endswith('trait.binary.output')
    shutil.copytree(tmp_path/'rer_inputs',results/'rer_inputs')
    reports={r.id:r for r in detect_reports(str(out),ROOT)}
    assert str(reports['rer'].slots[0].path).endswith('rer_inputs/validated_rer.output')
    script=build_script(reports['rer'],ROOT,use_singularity=False)
    run=subprocess.run(['bash'],input=script,cwd=tmp_path,text=True,capture_output=True,timeout=120)
    assert run.returncode==0,run.stdout+run.stderr
    content=html.read_text()
    assert '0.02' in content and 'Trait: trait' in content
