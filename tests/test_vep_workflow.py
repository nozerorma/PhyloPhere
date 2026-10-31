"""Serial real-Nextflow integration and content-cache validation for VEP."""
import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest
from test_vep import ROOT, fixture, database, cache_fixture, vc

pytestmark=pytest.mark.skipif(not shutil.which('nextflow'),reason='Run in phylophere')


@pytest.fixture
def project(tmp_path):
    shutil.copytree(ROOT/'subworkflows/VEP',tmp_path/'subworkflows/VEP',ignore=shutil.ignore_patterns('__pycache__'))
    (tmp_path/'workflows').mkdir();shutil.copy2(ROOT/'workflows/vep.nf',tmp_path/'workflows/vep.nf')
    fixture(tmp_path);database(tmp_path);database(tmp_path,'cosmic');cache_fixture(tmp_path/'cache')
    (tmp_path/'main.nf').write_text('''include { VEP } from './workflows/vep'
workflow {
    def scores=params.empty_channel ? Channel.empty() : (params.two_inputs ? Channel.of(file('scores.tsv'),file('scores.tsv')) : Channel.value(file('scores.tsv')))
    VEP(scores)
}
''')
    (tmp_path/'audit.config').write_text('''process.executor='local'
process.cpus=1
process.maxForks=1
executor.queueSize=1
trace.overwrite=true
params {
 outdir='out'; caas_map_dir='maps';gene_ensembl_file='proteins.tsv';vep_primateai_db='primateai.tsv';cosmic_db='cosmic.tsv'
 vep_ensembl=true;vep_mode='offline';vep_cache_dir='cache';vep_species='homo_sapiens';vep_assembly='GRCh38'
 vep_db_host='ensembldb.ensembl.org';vep_db_port=3306;two_inputs=false;empty_channel=false
}
''')
    (tmp_path/'bin').mkdir()
    fake=tmp_path/'bin/vep'
    fake.write_text('''#!/usr/bin/env python3
import os,sys
if '--help' in sys.argv:
 print('ensembl-vep : 115.0');sys.exit(0)
if os.environ.get('VEP_FAIL'):sys.exit(42)
args=sys.argv
out=args[args.index('--output_file')+1]
with open(out,'w') as f:
 f.write('## fake VEP 115\\n#Uploaded_variation\\tLocation\\tAllele\\tGene\\tFeature\\tFeature_type\\tConsequence\\tProtein_position\\tAmino_acids\\tExtra\\n')
 for line in open('candidate_variants.vcf'):
  if line.startswith('#'):continue
  chrom,pos,vid,ref,alt,*_=line.rstrip().split('\\t')
  f.write(f'{vid}\\t1:{pos}\\t{alt}\\tENSG1\\tENST1.2\\tTranscript\\tmissense_variant\\t1\\tA/T\\tCANONICAL=YES;ENSP=ENSP1.2\\n')
if os.environ.get('VEP_WARN'):open('vep_warnings.txt','w').write('reference mismatch')
''');fake.chmod(0o755)
    return tmp_path


def run(project,*args,ok=True,extra=None):
    env=os.environ.copy();env['PATH']=str(project/'bin')+os.pathsep+env['PATH'];env.update(extra or {})
    r=subprocess.run(['nextflow','-log','review.log','run','main.nf','-c','audit.config','-ansi-log','false','-with-trace','trace.tsv',*args],
                     cwd=project,env=env,text=True,capture_output=True,timeout=150)
    assert (r.returncode==0)==ok,r.stdout+r.stderr
    return r.stdout+r.stderr


def trace(project):
    return {r['name'].split(':')[-1].split(' (')[0]:r['status'] for r in vc.table(project/'trace.tsv')}


def test_workflow_serial_sources_and_snapshots(project):
    run(project)
    for name in ['primateai_mapped.tsv','cosmic_scores.tsv','ensembl_vep_mapped.tsv']:
        out=project/'out/vep'/name
        assert len(list(vc.table(out)))==2
        meta=vc.read_json(str(out)+'.manifest.json')
        assert meta['complete'] and meta['output_sha256']==vc.digest(out)
    assert (project/'out/vep/vep_inputs/validated_inputs.json').is_file()
    assert (project/'out/vep/ensembl_raw/candidate_variants.vcf').is_file()


def test_resume_hashes_content_sources_and_resource_bytes(project):
    run(project);run(project,'-resume')
    assert set(trace(project).values())=={'CACHED'}
    src=project/'subworkflows/VEP/local/src/map_to_cosmic.py'
    st=src.stat();src.write_text(src.read_text().replace('Census','CENSUS'));os.utime(src,ns=(st.st_atime_ns,st.st_mtime_ns))
    run(project,'-resume')
    tr=trace(project)
    assert tr['COSMIC_MAP']=='COMPLETED' and tr['PRIMATEAI_MAP']=='CACHED' and tr['VEP_PREPARE']=='CACHED'
    # Database content changes with preserved size and timestamp must invalidate its annotation.
    db=project/'primateai.tsv';st=db.stat();db.write_text(db.read_text().replace('0.8','0.7'));os.utime(db,ns=(st.st_atime_ns,st.st_mtime_ns))
    meta=Path(str(db)+'.vep.json');m=vc.read_json(meta);m['sha256']=vc.digest(db);meta.write_text(json.dumps(m))
    run(project,'-resume')
    assert trace(project)['PRIMATEAI_MAP']=='COMPLETED'
    assert next(vc.table(project/'out/vep/primateai_mapped.tsv'))['score_PAI3D']=='0.7'
    # A changed cache file without an updated manifest must fail, even when stat is identical.
    f=project/'cache/genome.fa';st=f.stat();f.write_text(f.read_text().replace('GCT','ACT'));os.utime(f,ns=(st.st_atime_ns,st.st_mtime_ns))
    run(project,'-resume',ok=False)
    assert any('checksum' in p.read_text() for p in (project/'work').glob('*/*/ensembl_failure.json'))


@pytest.mark.parametrize('flag',['VEP_FAIL','VEP_WARN'])
def test_technical_failure_stops_without_success_table(project,flag):
    run(project,ok=False,extra={flag:'1'})
    assert not (project/'out/vep/ensembl_vep_mapped.tsv').exists()
    assert list((project/'work').glob('*/*/ensembl_failure.json'))


@pytest.mark.parametrize('arg',['--empty_channel','--two_inputs'])
def test_input_cardinality_is_checked(project,arg):
    log=run(project,arg,'true',ok=False)
    assert 'exactly one' in log


def test_empty_scoring_publishes_headers_and_audits(project):
    scores=project/'scores.tsv';scores.write_text(scores.read_text().splitlines()[0]+'\n')
    run(project)
    for name in ['primateai_mapped.tsv','cosmic_scores.tsv','ensembl_vep_mapped.tsv']:
        assert list(vc.table(project/'out/vep'/name))==[]
        assert vc.read_json(str(project/'out/vep'/name)+'.manifest.json')['complete']


def test_missing_database_is_an_error_not_skip(project):
    (project/'cosmic.tsv').unlink()
    run(project,ok=False)


def test_database_mode_is_not_reused(project):
    run(project,'--vep_mode','database');run(project,'--vep_mode','database','-resume')
    tr=trace(project)
    assert tr['ENSEMBL_VEP_ANNOTATE']=='COMPLETED' and tr['PRIMATEAI_MAP']=='CACHED'


def test_legacy_reference_excluded_with_coverage(project):
    (project/'maps/001.map.tsv.vep.json').unlink()
    run(project)
    rows=list(vc.table(project/'out/vep/cosmic_scores.tsv'))
    assert {r['status'] for r in rows}=={'unverified_reference'}


@pytest.mark.parametrize('source',['primateai','ensembl','disabled'])
def test_independent_sources_and_optional_protein_table(project,source):
    cfg={'vep_primateai_db':'primateai.tsv' if source=='primateai' else '',
         'cosmic_db':'','vep_ensembl':source=='ensembl',
         'gene_ensembl_file':'proteins.tsv' if source=='ensembl' else ''}
    cfg['vep_mode']='auto'  # A prepared manifest deterministically resolves to offline.
    (project/'case.json').write_text(json.dumps(cfg))
    run(project,'-params-file','case.json')
    expected={'primateai':['primateai_mapped.tsv'],'ensembl':['ensembl_vep_mapped.tsv'],'disabled':[]}[source]
    assert sorted(p.name for p in (project/'out/vep').glob('*.tsv'))==expected
