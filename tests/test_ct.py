"""CT paired-design, residue, subset, CLI and batch-execution regressions."""
import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
CT = ROOT / 'subworkflows/CT/local'
sys.path.insert(0, str(CT))
from modules.design_io import read_labelings, read_design
from modules.perm_replay_io import simtrait_revive
from modules.perm_replay import _reconstruct_hit, _conserved_fields, VectorizedPermReplay

spec = importlib.util.spec_from_file_location('ct_prepare', CT/'scripts/prepare_labelings.py')
prep = importlib.util.module_from_spec(spec); spec.loader.exec_module(prep)


@pytest.fixture
def inputs(tmp_path):
    design = tmp_path/'design'; design.mkdir()
    (design/'traitfile_H1.tab').write_text('A\t1\t1\nB\t0\t1\nC\t1\t2\nD\t0\t2\n')
    (design/'contrast_hypotheses_pairs.tsv').write_text(
        'hypothesis_id\tpair\tspecies1\tspecies2\tpss_score\nH1\t1\tA\tB\t1\nH1\t2\tC\tD\t1\n')
    harvest = tmp_path/'harvest';harvest.mkdir()
    (harvest/'fop_labelings.tab').write_text(
        'b_2~H1\tA,C\tB,D\nb_2~H2\tE,C\tF,D\nb_1~H1\tA,C\tB,D\nb_1~H2\tE,C\tF,D\n')
    with (harvest/'fop_pairs.tsv').open('w') as out:
        out.write('cycle\thypothesis_id\tpair\tspecies1\tspecies2\tpss_score\n')
        for cycle in [2,1]:
            for hyp,fg,bg in [('H1','A','B'),('H2','E','F')]:
                out.write(f'b_{cycle}\t{hyp}\t1\t{fg}\t{bg}\t1\nb_{cycle}\t{hyp}\t2\tC\tD\t1\n')
    return design,harvest


def test_inactive_hypotheses_are_removed_with_matching_weights(inputs,tmp_path,monkeypatch):
    design,harvest=inputs;monkeypatch.chdir(tmp_path)
    prep.prepare(harvest,design,2,1,False)
    labels=read_labelings([tmp_path/'resample_perms.tab'])
    assert list(labels)==['b_0~H1','b_1~H1','b_2~H1']
    for block in (tmp_path/'blocks').iterdir():
        assert 'H2' not in (block/'resample_perms.tab').read_text()
        assert '\tH2\t' not in (block/'fop_pairs.tsv').read_text()


@pytest.mark.parametrize('n', [3,5])
def test_short_harvest_fails(inputs,tmp_path,monkeypatch,n):
    design,harvest=inputs;monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError,match='only 2 available'):
        prep.prepare(harvest,design,n,1,False)
    assert not (tmp_path/'resample_perms.tab').exists()


def test_zero_cycles_is_observed_only(inputs,tmp_path,monkeypatch):
    design,harvest=inputs;monkeypatch.chdir(tmp_path)
    prep.prepare(harvest,design,0,2,False)
    assert set(read_labelings([tmp_path/'resample_perms.tab']))=={'b_0~H1'}
    assert [p.name for p in (tmp_path/'blocks').iterdir()]==['b0']


def test_missing_and_wrong_pair_weights_fail(inputs,tmp_path,monkeypatch):
    design,harvest=inputs;monkeypatch.chdir(tmp_path)
    path=harvest/'fop_pairs.tsv';text=path.read_text();path.write_text(text.replace('\tA\tB\t','\tWRONG\tB\t'))
    with pytest.raises(ValueError,match='inconsistent'):
        prep.prepare(harvest,design,2,2,False)
    path.unlink()
    with pytest.raises(ValueError,match='requires matching'):
        prep.prepare(harvest,design,2,2,False)


def test_relaxed_missing_weights_are_explicit(inputs,tmp_path,monkeypatch):
    design,harvest=inputs;monkeypatch.chdir(tmp_path)
    (harvest/'fop_pairs.tsv').unlink()
    prep.prepare(harvest,design,2,2,False,True)
    assert json.loads((tmp_path/'ct_preparation.json').read_text())['relaxed_reasons']


def test_relaxed_reuse_requires_opt_in_and_keeps_generation_provenance(inputs,tmp_path,monkeypatch):
    design,harvest=inputs;monkeypatch.chdir(tmp_path)
    (harvest/'permulation_provenance.json').write_text(json.dumps({'relaxed_reasons':['Tier-2 draws used']}))
    with pytest.raises(ValueError,match='contains relaxed'):
        prep.prepare(harvest,design,2,2,False)
    prep.prepare(harvest,design,2,2,False,True)
    assert 'Tier-2 draws used' in json.loads((tmp_path/'ct_preparation.json').read_text())['relaxed_reasons']


@pytest.mark.parametrize('text', ['b_1\tA,C\tB\n','bad\n','b_1\tA,A\tB,C\n',
    'b_1\tA\tB\nb_1\tC\tD\n','b_1\tA\tA\n','b_x\tA\tB\n'])
def test_labeling_reader_rejects_malformed_designs(tmp_path,text):
    p=tmp_path/'bad.tab';p.write_text(text)
    with pytest.raises(ValueError):simtrait_revive(p)


@pytest.mark.parametrize('text',['A\t1\t1\n','A\t1\t2\nB\t0\t2\n','A\t1\t1\nA\t0\t1\n'])
def test_observed_design_validation(tmp_path,text):
    p=tmp_path/'bad.tab';p.write_text(text)
    with pytest.raises(ValueError):read_design(p)


@pytest.fixture
def cfg(tmp_path):
    p=tmp_path/'labelings.tab';p.write_text('b_0\tF1,F2,F3,F4\tB1,B2,B3,B4\n')
    return simtrait_revive(p)


def position(residues):
    return {sp:aa+'@0' for sp,aa in zip(['F1','F2','F3','F4','B1','B2','B3','B4'],residues)}


def test_pair_identity_survives_asymmetric_gaps(cfg):
    hit=_reconstruct_hit(cfg,position('-AAACCAC'),'b_0',1)
    assert hit['is_match'] and hit['conserved_pairs']=='1:3'


def test_shared_type_overlap_is_distinct_from_conserved_pairs(cfg):
    hit=_reconstruct_hit(cfg,position('AAACCDEF'),'b_0',1)
    assert hit['is_match'] and hit['conserved_pairs']=='0:'
    assert _conserved_fields(hit['conserved_pairs'],hit['shared_residue_count'])==['TRUE','0:','1']


@pytest.mark.parametrize('residues',['aaaacccc','AAAACCCC'])
def test_residue_case_agrees_between_kernel_and_export(cfg,residues):
    pos=position(residues)
    counts,hits=VectorizedPermReplay(cfg,list(pos)).count([pos],'G',0,0,0,0,0,0,1,'1,2,3',collect_hits=True)
    assert counts['G@0']==1 and _reconstruct_hit(cfg,pos,'b_0',1)['is_match']


def test_unknown_symbol_rejected(cfg):
    with pytest.raises(ValueError,match='Unsupported'):
        _reconstruct_hit(cfg,position('AAAA!CCC'),'b_0',1)


def test_cli_preserves_dotted_gene_and_emits_separate_overlap_count(tmp_path):
    config=tmp_path/'cfg.tab';config.write_text(''.join(f'F{i}\t1\t{i}\nB{i}\t0\t{i}\n' for i in range(1,5)))
    labels=tmp_path/'labs.tab';labels.write_text('b_0\tF1,F2,F3,F4\tB1,B2,B3,B4\n')
    ali=tmp_path/'G.v1.fa';ali.write_text(''.join(f'>{s}\n{v.split("@")[0]}\n' for s,v in position('aaaacccc').items()))
    out=tmp_path/'discovery.tab'
    cmd=[sys.executable,str(CT/'ct'),'perm-replay','-a',str(ali),'-t',str(config),'-s',str(labels),'--fmt','fasta','--max_conserved','1','--export_b0_discovery',str(out)]
    r=subprocess.run(cmd,text=True,capture_output=True)
    assert r.returncode==0,r.stderr
    with out.open() as fh:row=next(csv.DictReader(fh,delimiter='\t'))
    assert row['gene']=='G.v1' and row['shared_residue_count']=='0'


def test_finished_worker_failure_is_collected():
    source=(CT/'scripts/run_ct_perm_replay_batch.sh').read_text()
    fn=source[source.index('wait_for_all() {'):source.index('\n# Each of')]
    cmd='set -euo pipefail\nbatch_id=probe\nterminate_children() { :; }\n'+fn+'\n(exit 42) &\nworker_pids=("$!")\nsleep 0.1\nwait_for_all\n'
    assert subprocess.run(['bash','-c',cmd],capture_output=True).returncode==42


def test_batch_accepts_container_prefix_as_one_argument(tmp_path):
    # A local entrypoint stand-in checks the same two-part argv without a container.
    fake=tmp_path/'entrypoint';fake.write_text('#!/bin/bash\nexec "$@"\n');fake.chmod(0o755)
    program=tmp_path/'fake-ct';program.write_text('#!/bin/bash\nexit 0\n');program.chmod(0o755)
    manifest=tmp_path/'manifest.tsv';manifest.write_text('G\tG.fa\n')
    cmd=['bash',str(CT/'scripts/run_ct_perm_replay_batch.sh'),'--batch-id','test','--manifest',str(manifest),'--caas-config','unused','--resampled-path','unused','--workers','1','--ali-format','fasta','--ct-bin',f'{fake} {program}']
    r=subprocess.run(cmd,cwd=tmp_path,text=True,capture_output=True,timeout=10)
    assert r.returncode==0,r.stderr


def test_pair_metadata_export_is_accepted_by_the_scoring_boundary(cfg):
    # Read the newly extended row with the existing consumer, in a fresh process
    # to avoid mixing the pipeline's independent src packages in pytest imports.
    code=f"""
import sys
sys.path.insert(0,{str(ROOT/'subworkflows/CT_DISAMBIGUATION/local')!r})
from src.utils.perm_stream import entry_from_fields
entry=entry_from_fields(['b_1','G','0','AAA/CCAC','2','TRUE','1:3','1'])
assert entry.conserved_pair=='3' and entry.is_conserved_meta
"""
    r=subprocess.run([sys.executable,'-c',code],text=True,capture_output=True)
    assert r.returncode==0,r.stderr
