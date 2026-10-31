"""Sequential local Nextflow contract checks; one CPU per task."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

ROOT=Path(__file__).resolve().parents[1]
LOCAL=ROOT/'subworkflows/CT_DISAMBIGUATION/local'
pytestmark=pytest.mark.skipif(not shutil.which('nextflow'),reason='Run in phylophere')


@pytest.fixture
def project(tmp_path):
    p=tmp_path/"project's directory";p.mkdir()
    shutil.copytree(LOCAL,p/'subworkflows/CT_DISAMBIGUATION/local',ignore=shutil.ignore_patterns('__pycache__'))
    for name in ['subworkflows/CT_DISAMBIGUATION/ct_observed.nf','subworkflows/CT/caas_permulation.nf',
                 'subworkflows/CT/local/modules/design_io.py','workflows/ct_observed.nf']:
        target=p/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy(ROOT/name,target)
    (p/"design's file.tab").write_text('A\t1\t1\nB\t0\t1\n')
    (p/"discovery's file.tab").write_text('gene\tposition\ttrait\tcaas\tpattern\n')
    (p/'tree.nwk').write_text('(A:1,B:1);\n')
    (p/'alignments').mkdir();(p/'alignments/G.fa').write_text('>A\nAAAA\n>B\nAAAA\n')
    (p/'audit.config').write_text('''
process.executor='local'
process.cpus=1
process.maxForks=1
params {
 caas_config="design's file.tab"; discovery_from="discovery's file.tab"; tree='tree.nwk'
 alignment='alignments'; outdir='out'; use_singularity=false; use_apptainer=false
 ct_disambig_asr_model='poisson'; ct_disambig_posterior_threshold=0.1
 ct_disambig_asr_cache_dir="cache's folder"; ct_disambig_max_tasks_per_child=2
 ct_disambig_hypotheses_pairs=''; scoring_hypotheses_pairs=''; gene_ensembl_file='';ali_format='auto';live=false
}
''')
    (p/'main.nf').write_text('''
include { CT_OBSERVED } from './workflows/ct_observed.nf'
workflow {
 CT_OBSERVED(Channel.value(file(params.discovery_from)),
             params.live ? Channel.value(file(params.caas_config)) : null,
             params.live ? Channel.value(file(params.tree)) : null,null,Channel.value(file('NO_FILE')))
}
''')
    return p


def run(p,*args,success=True):
    r=subprocess.run(['nextflow','-log','review.log','run','main.nf','-c','audit.config','-w','work','-ansi-log','false',*args],
                     cwd=p,text=True,capture_output=True,timeout=120)
    assert (r.returncode==0)==success,r.stdout+r.stderr
    return r.stdout+r.stderr


def tasks(p):
    return [f for f in (p/'work').glob('*/*/.command.sh') if 'run_observed.py' in f.read_text()]


@pytest.mark.parametrize('live',[False,True])
def test_empty_outputs_stage_quoted_paths_and_use_no_asr(project,live):
    run(project,'--live',str(live).lower())
    assert len((project/'out/ct_disambiguation/caas_convergence_master.csv').read_text().splitlines())==1
    assert 'shared_residue_count' in (project/'out/meta_caas/meta_caas/meta_caas.tsv').read_text()
    args=json.loads((tasks(project)[0].parent/'observed_args.json').read_text())
    assert "inputs/design/design's file.tab" in args['scoring']
    assert not (project/"cache's folder").exists()


def test_explicit_missing_weights_fail_instead_of_falling_back(project):
    output=run(project,'--ct_disambig_hypotheses_pairs','MISSING.tsv',success=False)
    assert 'Requested hypothesis weights file not found' in output


def test_source_and_alignment_changes_invalidate_resume(project):
    nested=project/'alignments/nested';nested.mkdir()
    alignment=nested/'G.fa';(project/'alignments/G.fa').rename(alignment)
    design=project/'selected design';design.mkdir()
    trait=design/'traitfile_H1.tab';shutil.copy(project/"design's file.tab",trait)
    def execute(*args):
        return run(project,'--caas_config','selected design',*args)
    execute()
    execute('-resume');assert len(tasks(project))==1
    src=project/'subworkflows/CT_DISAMBIGUATION/local/observed_b0_main.py'
    src.write_text(src.read_text()+'\n# source change\n')
    execute('-resume');assert len(tasks(project))==2
    before=alignment.stat()
    alignment.write_text(alignment.read_text().replace('AAAA','AAAD'))
    os.utime(alignment,ns=(before.st_atime_ns,before.st_mtime_ns))
    execute('-resume');assert len(tasks(project))==3
    before=trait.stat();trait.write_text('C\t1\t1\nD\t0\t1\n')
    os.utime(trait,ns=(before.st_atime_ns,before.st_mtime_ns))
    execute('-resume');assert len(tasks(project))==4


def test_weights_beside_selected_design_are_used(project):
    d=project/'selected design';d.mkdir()
    (d/'traitfile_H1.tab').write_text('A\t1\t1\nB\t0\t1\n')
    (d/'traitfile_H2.tab').write_text('C\t1\t1\nD\t0\t1\n')
    (d/'contrast_hypotheses_pairs.tsv').write_text('hypothesis_id\tpair\tspecies1\tspecies2\tpss_score\nH1\t1\tA\tB\t1\nH2\t1\tC\tD\t2\n')
    run(project,'--caas_config','selected design')
    ledger=json.loads((project/'out/ct_disambiguation/observed_scoring.json').read_text())
    assert ledger['batches'][0]['weights_sha256'] and ledger['batches'][0]['weighted']


def test_live_core_assembler_accepts_validated_empty_slice(project):
    b0=project/'b0';b0.mkdir()
    (b0/'G.b0.discovery.tsv').write_text('gene\tposition\ttrait\tcaas\tpattern\n')
    (b0/'G.b0.background').write_text('G\tNULL\n')
    r=subprocess.run([sys.executable,str(LOCAL/'observed_b0_main.py'),'--b0-dir',str(b0),
        '--design',str(project/"design's file.tab"),'--alignment-dir',str(project/'alignments'),
        '--tree',str(project/'tree.nwk'),'--output-dir',str(b0),'--asr-cache-dir',str(project/'cache'),'--workers','1'],capture_output=True,text=True,timeout=30)
    assert r.returncode==0,r.stdout+r.stderr
    (project/'main.nf').write_text('''
include { CAAS_CORE_OBSERVED } from './subworkflows/CT/caas_permulation.nf'
workflow { CAAS_CORE_OBSERVED(Channel.value([file('b0')]),Channel.value(file(params.caas_config))) }
''')
    run(project)
    assert (project/'out/caastools/background_genes.output').read_text()==''
    assert len((project/'out/ct_disambiguation/caas_convergence_master.csv').read_text().splitlines())==1


@pytest.mark.skipif(not shutil.which('codeml'),reason='PAML codeml unavailable')
def test_real_live_tile_preserves_dotted_gene_and_overlap_strictness(project):
    shutil.copytree(ROOT/'subworkflows/CT/local',project/'subworkflows/CT/local',dirs_exist_ok=True,
                    ignore=shutil.ignore_patterns('__pycache__'))
    cfg=project/"design's file.tab"
    cfg.write_text(''.join(f'{a}\t1\t{i}\n{b}\t0\t{i}\n' for i,(a,b) in enumerate(zip('ACEG','BDFH'),1)))
    (project/'tree.nwk').write_text('(((A:.1,B:.1):.1,(C:.1,D:.1):.1):.1,((E:.1,F:.1):.1,(G:.1,H:.1):.1):.1);\n')
    base='ACDEFGHIKLMNPQRSTVWY'*4
    states=dict(zip('ACEG','DDAA'));states.update(zip('BDFH','AEEE'))
    alignment=project/'alignments/AC010616.1.fa'
    alignment.write_text(''.join(f'>{s}\n{states[s]}{base[1:]}\n' for s in 'ABCDEFGH'))
    block=project/'b0';block.mkdir()
    (block/'resample_perms.tab').write_text('b_0\tA,C,E,G\tB,D,F,H\n')
    with (project/'audit.config').open('a') as fh:
        fh.write("""
params {
 caas_perms_postproc=false; caas_map_dir=''; strict_conservation_output=true; min_divergent_fraction=0.5
 min_non_overlapping=2; patterns='1,2,3,4';miss_pair=false;ali_format='fasta'
 max_bg_gaps_fraction=0;max_fg_gaps_fraction=0;max_gaps_fraction=0
 max_bg_miss_fraction=0;max_fg_miss_fraction=0;max_miss_fraction=0
 max_gaps_per_position=0.5
}
""")
    (project/'main.nf').write_text('''
include { CAAS_CORE_BATCHED } from './subworkflows/CT/caas_permulation.nf'
workflow {
 CAAS_CORE_BATCHED(Channel.value(tuple('batch_1',1,'AC010616.1\\tAC010616.1.fa\\n',[file('alignments/AC010616.1.fa')],file('b0'))),
                   file(params.caas_config),file(params.tree),file('NO_FILE'),file('NO_FILE'))
}
''')
    run(project)
    command=next(f for f in (project/'work').glob('*/*/.command.sh') if 'disambiguation_perms_main.py' in f.read_text())
    task=command.parent
    assert (task/'b0_observed/AC010616.1.master.csv.gz').exists()
    assert not list((task/'b0_observed_strict').glob('*.master.csv.gz'))
    ledger=json.loads((task/'b0_observed_strict/observed_scoring.json').read_text())
    assert ledger['complete'] and ledger['genes']=={}
    import pyarrow.parquet as pq
    assert pq.read_table(task/'null_tables_strict/scores.parquet').num_rows==0
    assert len(list((task/"cache's folder").rglob('provenance.json')))==1
