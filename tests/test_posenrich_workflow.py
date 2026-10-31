"""SCORING/POSENRICH integration on the same canonical null, literal IDs and fixed families."""
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
import os
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from test_scoring_workflow import project as scoring_project
from posenrich_fixtures import annotation, seal_block, digest

ROOT=Path(__file__).resolve().parents[1]


def run(p, *args, ok=True):
    r=subprocess.run(['nextflow','-log','posenrich.log','run','main.nf','-params-file','params.json','-c','audit.config',
                      '-ansi-log','false','-with-trace','trace.tsv',*args],cwd=p,capture_output=True,text=True,timeout=240)
    log=r.stdout+r.stderr+(p/'posenrich.log').read_text()[-6000:]
    assert (r.returncode==0)==ok,log
    return log


def rows(p, quoted=False):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t',quoting=csv.QUOTE_MINIMAL if quoted else csv.QUOTE_NONE))


def trace(p):
    return {r['name'].split(':')[-1].split(' (')[0]:r['status'] for r in rows(p/'trace.tsv')}


@pytest.fixture
def integrated(scoring_project):
    p=scoring_project
    shutil.copytree(ROOT/'subworkflows/POSENRICH',p/'subworkflows/POSENRICH',ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copytree(ROOT/'subworkflows/VEP',p/'subworkflows/VEP',ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copy2(ROOT/'workflows/posenrich.nf',p/'workflows/posenrich.nf')
    block=p/'prior/blocks/blk_000001'
    table=pq.read_table(block/'positions.parquet')
    table=table.append_column('top',table['all']).append_column('bottom',pa.array([float('nan')],pa.float64()))
    pq.write_table(table,block/'positions.parquet');seal_block(block)
    genes=(p/'genes.txt').read_text().splitlines()
    (p/'custom.tsv').write_text(''.join(f'{g}\t0\tT\tpredeclared marker\n' for g in genes))
    annotation(p/'custom.tsv',genes,tested=p/'tested',background=p/'genes.txt')
    params=json.loads((p/'params.json').read_text());params.update(scoring=True,posenrich_domains=False,posenrich_ucr=False,
        accumulation_fdr=.2,posenrich_custom_file='custom.tsv',posenrich_min_genes=2,posenrich_holdout=0,posenrich_radios='0,5,10',
        position_registry_memory='1 GB',position_module_memory='1 GB',posenrich_within_gene_memory='1 GB')
    (p/'params.json').write_text(json.dumps(params))
    (p/'main.nf').write_text('''include { SCORING } from './workflows/scoring'
include { POSENRICH } from './workflows/posenrich'
workflow {
 def bg = Channel.value(file('genes.txt'))
 def obs = Channel.value(file('observed.tsv'))
 def tested = Channel.value(file('tested'))
 def sc = SCORING(obs,null,null,null,null,null,null,bg,null,null,null,null,null,null,null,null,null,null,null,tested,null)
 POSENRICH(true,obs,bg,tested,Channel.value(file('prior/null_positions.parquet')),
           Channel.value([file('prior/blocks/blk_000001')]),sc.position_scores,sc.postproc,sc.null_universe_summary,
           Channel.empty(),Channel.empty())
}
''')
    return p


def test_joint_families_and_cache(integrated):
    p=integrated
    run(p)
    scores=rows(p/'out/scoring/position_scores.tsv')
    genes=set((p/'genes.txt').read_text().splitlines())
    assert {r['Gene'] for r in scores}==genes
    for r in scores: assert float(r['p.emp'])==pytest.approx(.05)
    positional_gene=rows(p/'out/scoring/gene_evidence.tsv')
    assert all(float(r['p.gen'])==pytest.approx(.05) for r in positional_gene if r['ranking']=='global')
    universe=p/'out/posenrich/universe/universe.tsv.gz'
    import pandas as pd
    u=pd.read_csv(universe,sep='\t',dtype={'gene':str},keep_default_na=False)
    assert len(u)==15 and u.observed_positive.sum()==5
    with (p/'out/posenrich/accumulation/all/accumulation_all_results.csv').open() as f:
        acc=list(csv.DictReader(f))
    assert {r['Gene'] for r in acc}==genes
    assert all(int(r['PoolPositions'])==3 for r in acc)
    assert float(next(r for r in acc if r['Gene']=='001')['EmpiricalP'])==pytest.approx(.1)
    assert all(float(r['EmpiricalP'])==pytest.approx(.05) for r in acc if r['Gene']!='001')
    for side in ('top','bottom','all'):
        assert (p/f'out/posenrich/accumulation/{side}/accumulation_{side}_results.csv.manifest.json').is_file()
    wg=rows(p/'out/posenrich/within_gene/within_gene_results.tsv',True)
    assert wg[0]['status']=='tested' and float(wg[0]['p'])==pytest.approx(.05)
    assert {r['Gene'] for r in rows(p/'out/posenrich/characterization/position_characterization.tsv',True)}==genes
    assert (p/'out/posenrich/13.POSENRICH_report_audit.html').is_file()
    run(p,'-resume');assert all(s=='CACHED' for s in trace(p).values())
    # Regenerate from the published nested accumulation directories; no inference is rerun.
    rmd=p/'subworkflows/POSENRICH/local/13.POSENRICH_report.Rmd'
    out=p/'out/posenrich'
    params=','.join(f"{key}_dir='{out/name}'" for key,name in [('universe','universe'),('registry','registry'),('characterization','characterization'),('accumulation','accumulation'),('within_gene','within_gene'),('proximity','proximity')])
    regenerated=subprocess.run(['Rscript','-e',f"rmarkdown::render('{rmd}',params=list({params},traitname='audit'),output_file='{p/'regenerated.html'}',quiet=TRUE)"],capture_output=True,text=True,cwd=p)
    assert regenerated.returncode==0,regenerated.stdout+regenerated.stderr
    html=(p/'regenerated.html').read_text()
    assert 'Empirical expected' in html and 'hypergeometric' in html
    settings=json.loads((out/'report_settings.json').read_text())
    assert settings['fdr_threshold']==.2
    assert 'Displayed BH threshold:' in html and '0.2 (producer p/q are preserved)' in ' '.join(html.split())
    assert all(s=='CACHED' for s in trace(p).values())
    # Add never-detected tested positions in the same genes. Raw maxima and counts remain fixed;
    # positional BH, uniform opportunity and centered within-gene residuals must use the enlarged U.
    before={ (r['Gene'],r['Position'],r['side']):r for r in scores }
    before_stat=float(wg[0]['statistic'])
    (p/'tested').write_text(''.join(f'{g}\t0,1,2,3,4,5\n' for g in genes))
    annotation(p/'custom.tsv',genes,tested=p/'tested',background=p/'genes.txt')
    run(p,'-resume')
    after=rows(p/'out/scoring/position_scores.tsv')
    for r in after:
        old=before[(r['Gene'],r['Position'],r['side'])]
        assert r['p.emp']==old['p.emp']
        assert float(r['p.adj_bh'])>float(old['p.adj_bh'])
    newer=rows(p/'out/scoring/gene_evidence.tsv')
    assert [(r['gene'],r['ranking'],r['p.gen'],r['p.gen.adj']) for r in newer]==[(r['gene'],r['ranking'],r['p.gen'],r['p.gen.adj']) for r in positional_gene]
    with (p/'out/posenrich/accumulation/all/accumulation_all_results.csv').open() as f:
        acc2={r['Gene']:r for r in csv.DictReader(f)}
    for r in acc:
        assert acc2[r['Gene']]['EmpiricalP']==r['EmpiricalP']
        assert float(acc2[r['Gene']]['PValue']) != float(r['PValue'])
    wg2=rows(p/'out/posenrich/within_gene/within_gene_results.tsv',True)
    assert float(wg2[0]['statistic'])!=before_stat
    # Same-size/mtime content change in an annotation must invalidate the registry.
    f=p/'custom.tsv';st=f.stat();meta=p/'custom.tsv.manifest.json';meta_st=meta.stat();meta_size=meta.stat().st_size;f.write_text(f.read_text().replace('predeclared','predeclaret'))
    os.utime(f,ns=(st.st_atime_ns,st.st_mtime_ns));annotation(f,genes,tested=p/'tested',background=p/'genes.txt')
    assert meta.stat().st_size==meta_size
    os.utime(meta,ns=(meta_st.st_atime_ns,meta_st.st_mtime_ns))
    run(p,'-resume');assert trace(p)['POSITION_REGISTRY']!='CACHED'


def test_reference_and_null_corruption_fail(integrated):
    p=integrated
    meta=p/'observed.tsv.scoring.json';m=json.loads(meta.read_text());m['scientific']['model']='WAG';meta.write_text(json.dumps(m))
    log=run(p,ok=False)
    assert 'scientific mismatch' in log


def test_unverified_annotation_remains_descriptive(integrated):
    p=integrated
    (p/'custom.tsv.manifest.json').unlink()
    run(p)
    out=p/'out/posenrich'
    assert rows(out/'within_gene/within_gene_results.tsv',True)==[]
    m=json.loads((out/'within_gene/within_gene_manifest.json').read_text())
    assert m['state']=='no_annotations_in_universe' and m['counts']['n_null_cycles']==19
    r=json.loads((out/'registry/manifest.json').read_text())
    assert r['annotations']['custom_marker_file']['status']=='unverified_source'
    assert rows(out/'registry/sites.tsv',True)
    assert (out/'13.POSENRICH_report_audit.html').is_file()
