"""Serial local Nextflow checks, including resume after stat-preserving edits."""
import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]
pytestmark = pytest.mark.skipif(not shutil.which('nextflow'), reason='Run in phylophere')


@pytest.fixture
def project(tmp_path):
    # Copy sources so cache tests never edit the shared checkout.
    selection=tmp_path/'subworkflows/SELECTION'
    shutil.copytree(ROOT/'subworkflows/SELECTION', selection, ignore=shutil.ignore_patterns('__pycache__'))
    (tmp_path/'subworkflows/CT').mkdir()
    shutil.copy2(ROOT/'subworkflows/CT/ct_alignment_files.nf',tmp_path/'subworkflows/CT/ct_alignment_files.nf')
    (tmp_path/'alignments').mkdir()
    (tmp_path/'alignments/g.v1.fa').write_text('>A\nACDE\n>B\nXXXX\n>C\nAAAA\n')
    (tmp_path/'alignments/h.phy').write_text('2 4\nA  AAAA\nB  CCCC\n')
    (tmp_path/'alignments/empty.fa').write_text('>Z\nAAAA\n')
    (tmp_path/'tree.nwk').write_text('(A:1,B:1);\n')
    (tmp_path/'pool.tab').write_text('A\t1\t9\nB\t0\t3\n')
    (tmp_path/'main.nf').write_text('''include { listAlignmentFiles; sampleAlignmentFiles } from './subworkflows/CT/ct_alignment_files'
include { SELECTION_PREP } from './subworkflows/SELECTION/selection_prep'
workflow {
    def species = params.multiple_species ? Channel.of(file('pool.tab'), file('pool.tab')) : Channel.value(file('pool.tab'))
    def trees = params.multiple_trees ? Channel.of(file('tree.nwk'), file('tree.nwk')) : Channel.empty()
    def discovery = params.discovery ? Channel.value(file(params.discovery)) : Channel.empty()
    if (params.toy_mode && !params.discovery) {
        file('ct_sample.txt').text = sampleAlignmentFiles(listAlignmentFiles(params.alignment), params.toy_n, params.seed).collect { it.baseName }.join('\\n')
    }
    def prep = SELECTION_PREP(species, trees, Channel.empty(), Channel.empty(), discovery)
    prep.prep_audit.map { it.name }.collectFile(name:'audit_files.txt', storeDir:params.outdir, newLine:true)
    prep.filtered_fasta_top_ch.map { gid, f -> "${gid}\\t${f.text.replace('\\n', '|')}" }.collectFile(name:'top.tsv', storeDir:params.outdir, newLine:true)
    prep.filtered_fasta_bottom_ch.map { gid, f -> "${gid}\\t${f.text.replace('\\n', '|')}" }.collectFile(name:'bottom.tsv', storeDir:params.outdir, newLine:true)
}
''')
    (tmp_path/'audit.config').write_text("""process.executor='local'
process.cpus=1
trace.overwrite=true
params {
    alignment='alignments'; tree='tree.nwk'; outdir='out'; selection_prep_batch_size=2
    use_singularity=false; use_apptainer=false; toy_mode=false; toy_n=2; seed=0
    multiple_species=false; multiple_trees=false; discovery=''
}
""")
    return tmp_path


def run(project, *args, success=True):
    result=subprocess.run(['nextflow','-log','review.log','run','main.nf','-c','audit.config',
        '-ansi-log','false','-with-trace','trace.tsv',*args],cwd=project,text=True,capture_output=True,timeout=120)
    if success: assert result.returncode == 0, result.stdout+result.stderr
    else: assert result.returncode != 0, result.stdout+result.stderr
    return result.stdout+result.stderr


def rows(path):
    return set(path.read_text().splitlines())


def test_single_batch_parity_and_last_singleton_batch(project):
    run(project)
    batch=rows(project/'out/top.tsv')
    assert batch==rows(project/'out/bottom.tsv')
    assert {r.split('\t')[0] for r in batch}=={'g.v1','h'}
    audits={p.name:json.loads(p.read_text()) for p in (project/'out/selection/prep_audit').glob('*.json')}
    assert len(audits)==3 and audits['empty.prep.json']['status']=='excluded_no_usable_tree_taxa'
    run(project,'--selection_prep_batch_size','1','--outdir','single')
    assert batch==rows(project/'single/top.tsv')==rows(project/'single/bottom.tsv')
    assert {p.name:json.loads(p.read_text()) for p in (project/'single/selection/prep_audit').glob('*.json')}==audits


@pytest.mark.parametrize('args, message', [
    (['--selection_prep_batch_size','0'],'positive integer'),
    (['--selection_prep_batch_size','1.5'],'positive integer'),
    (['--multiple_species','true'],'expected one species pool'),
    (['--multiple_trees','true'],'expected one tree'),
    (['--alignment','missing'],'alignment directory not found'),
    (['--toy_mode','true','--toy_n','0'],'positive integer'),
])
def test_invalid_inputs_fail(project,args,message):
    assert message in run(project,*args,success=False)


@pytest.mark.parametrize('batch', ['1','2'])
def test_corrupt_gene_aborts_actual_task(project,batch):
    (project/'alignments/bad.phy').write_text('broken')
    log=run(project,'--selection_prep_batch_size',batch,success=False)
    assert 'Alignment preparation failed' in log and 'bad' in log


def test_duplicate_gene_ids_fail(project):
    (project/'alignments/h.fa').write_text('>A\nAAAA\n')
    assert 'unique gene IDs' in run(project,success=False)


def test_toy_discovery_uses_named_column_and_preserves_dots(project):
    (project/'discovery.tsv').write_text('position\tgene\n0\tg.v1\n2\tg.v1\n')
    run(project,'--toy_mode','true','--discovery','discovery.tsv')
    assert {r.split('\t')[0] for r in rows(project/'out/top.tsv')}=={'g.v1'}
    (project/'discovery.tsv').write_text('position\tgene\n')
    run(project,'--toy_mode','true','--discovery','discovery.tsv','--outdir','empty')
    assert not list((project/'empty/selection/prep_audit').glob('*.json'))


def test_toy_missing_alignment_fails(project):
    (project/'discovery.tsv').write_text('gene\nmissing.v2\n')
    assert 'no supported alignment' in run(project,'--toy_mode','true','--discovery','discovery.tsv',success=False)


def stat_preserving_edit(path, old, new):
    before=path.stat(); body=path.read_text(); assert len(old)==len(new) and old in body
    path.write_text(body.replace(old,new)); os.utime(path, ns=(before.st_atime_ns,before.st_mtime_ns))


@pytest.mark.parametrize('batch', ['1', '2'])
def test_resume_checks_content_and_helpers(project, batch):
    config=project/'audit.config'
    config.write_text(config.read_text().replace('selection_prep_batch_size=2', 'selection_prep_batch_size='+batch))
    run(project)
    run(project,'-resume')
    assert 'CACHED' in (project/'trace.tsv').read_text()
    stat_preserving_edit(project/'alignments/g.v1.fa','ACDE','VCDE')
    run(project,'-resume')
    assert any('VCDE' in r for r in rows(project/'out/top.tsv'))
    assert 'COMPLETED' in (project/'trace.tsv').read_text()
    stat_preserving_edit(project/'pool.tab','A\t1\t9\nB\t0\t3','A\t0\t9\nB\t1\t3')
    run(project,'-resume')
    assert (project/'out/selection/species_sets/top_species.txt').read_text()=='B\n'
    stat_preserving_edit(project/'tree.nwk','B:1','C:1')
    run(project,'-resume')
    audits=list((project/'work').glob('*/*/g.v1.prep.json'))
    assert any('C' in json.loads(p.read_text())['retained_taxa'] for p in audits)
    helper=project/'subworkflows/SELECTION/local/src/prepare_alignments.py'
    stat_preserving_edit(helper,'excluded_no_usable_tree_taxa','excluded_no_usable_tree_TAXA')
    run(project,'-resume')
    assert json.loads((project/'out/selection/prep_audit/empty.prep.json').read_text())['status']=='excluded_no_usable_tree_TAXA'
    assert 'COMPLETED' in (project/'trace.tsv').read_text()


def test_seed_zero_matches_ct_sample_and_is_reproducible(project):
    run(project,'--toy_mode','true')
    expected=set((project/'ct_sample.txt').read_text().splitlines())
    assert {p.name.removesuffix('.prep.json') for p in (project/'out/selection/prep_audit').glob('*.json')}==expected
    run(project,'--toy_mode','true','--outdir','repeat','-resume')
    assert {p.name.removesuffix('.prep.json') for p in (project/'repeat/selection/prep_audit').glob('*.json')}==expected


def test_all_genes_excluded_still_emit_complete_audit(project):
    for path in (project/'alignments').iterdir():
        path.write_text('>Z\nAAAA\n')
    run(project)
    assert not (project/'out/top.tsv').exists()
    assert not (project/'out/bottom.tsv').exists()
    assert len(rows(project/'out/audit_files.txt'))==3
    assert all(json.loads(p.read_text())['status']=='excluded_no_usable_tree_taxa'
               for p in (project/'out/selection/prep_audit').glob('*.json'))


@pytest.mark.parametrize('batch', ['1', '2'])
def test_full_identifiers_survive_suffix_and_shell_quoting(project, batch):
    gid="alpha.filtered.fa.v2 'quoted'"
    (project/'alignments'/f'{gid}.fa').write_text('>A\nACDE\n>B\nCCCC\n')
    run(project,'--selection_prep_batch_size',batch)
    assert gid in {r.split('\t')[0] for r in rows(project/'out/top.tsv')}
    assert (project/'out/selection/prep_audit'/f'{gid}.prep.json').is_file()
