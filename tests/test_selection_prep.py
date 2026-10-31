"""Selection preparation: format parity, exclusions and fail-closed contracts."""
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'subworkflows/SELECTION/local/src'
sys.path.insert(0, str(SRC))
from alignment_io import read_alignment, tree_taxa
from extract_extreme_species import load_species_sets
from prepare_alignments import prepare


@pytest.mark.parametrize('body', [
    'A\t1\nA\t0\nB\t0\n', 'A\t1\nB\t0\nC\t2\n',
    'A\t1\nB\n', '\t1\nB\t0\n', 'A\t1\nA\t1\nB\t0\n',
    'A\t1\nB\t0\nC\t1\textra\textra\n', 'A A\t1\nB\t0\n',
    'A\t1\n', '',
    'species,global_label\nA,high_extreme\nB,low_extreme\nA,low_extreme\n',
    'species,global_label\nA,high_extreme\nB,low_extreme\nC,typo\n',
])
def test_invalid_species_pools_fail(tmp_path, body):
    p = tmp_path/'species'; p.write_text(body)
    with pytest.raises(ValueError):
        load_species_sets(p)


@pytest.mark.parametrize('body', [
    'A\t1\t91\nB\t0\t2\nC\t1\t91\n',
    'species\tcontrast_group\tpair\nA\t1\tx\nB\t0\tx\nC\t1\ty\n',
    'species,global_label\nA,high_extreme\nB,low_extreme\nC,high_extreme\nD,intermediate\nE,normal\n',
])
def test_global_pool_is_not_a_paired_design(tmp_path, body):
    p = tmp_path/'species'; p.write_text(body)
    assert load_species_sets(p) == ({'A', 'C'}, {'B'})


@pytest.mark.parametrize('body', [
    '>A\nAAAA\n>A\nAAAA\n', '>A\nAAAA\n>B\nAAA\n',
    '>A\n\n>B\nAAAA\n', '>\nAAAA\n', '>A\nAA!A\n',
    '', 'not PHYLIP\n',
])
def test_malformed_alignments_fail(tmp_path, body):
    p = tmp_path/'input.fa'; p.write_text(body)
    with pytest.raises(ValueError):
        read_alignment(p)


def test_formats_produce_identical_sequences_and_exclusions(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    fa = tmp_path/'x.fa'; fa.write_text('>A\nACDE\n>B\nX?-X\n>C\nAAAA\n')
    phy = tmp_path/'x.phy'; phy.write_text('3 4\nA  ACDE\nB  X?-X\nC  AAAA\n')
    a = prepare('fa', fa, {'A','B'}); b = prepare('phy', phy, {'A','B'})
    assert (tmp_path/'fa.filtered.fa').read_text() == (tmp_path/'phy.filtered.fa').read_text()
    for key in ['status','retained_taxa','excluded_uninformative','excluded_off_tree','alignment_columns']:
        assert a[key] == b[key]
    assert a['retained_taxa'] == ['A']
    assert a['excluded_uninformative'] == ['B'] and a['excluded_off_tree'] == ['C']


@pytest.mark.parametrize('body', ['>A\nXXXX\n', '>C\nAAAA\n'])
def test_biological_empty_is_audited_and_removes_stale_output(tmp_path, monkeypatch, body):
    monkeypatch.chdir(tmp_path)
    p = tmp_path/'input'; p.write_text(body)
    (tmp_path/'g.filtered.fa').write_text('stale')
    result = prepare('g', p, {'A','B'})
    assert result['status'] == 'excluded_no_usable_tree_taxa'
    assert not (tmp_path/'g.filtered.fa').exists()
    assert json.loads((tmp_path/'g.prep.json').read_text()) == result


@pytest.mark.parametrize('tree', ['();', '(A:1,A:1);', '((:1,B:1):1,C:1);', 'broken'])
def test_invalid_trees_fail(tmp_path, tree):
    p=tmp_path/'tree.nwk';p.write_text(tree)
    with pytest.raises(Exception): tree_taxa(p)


def test_batch_failure_is_nonzero_and_records_every_gene(tmp_path):
    (tmp_path/'alignments').mkdir()
    (tmp_path/'alignments/good.fa').write_text('>A\nAAAA\n>B\nCCCC\n')
    (tmp_path/'alignments/bad.phy').write_text('broken')
    (tmp_path/'tree.nwk').write_text('(A:1,B:1);')
    (tmp_path/'manifest').write_text('good\tgood.fa\nbad\tbad.phy\n')
    result=subprocess.run(['bash',str(SRC/'run_prep_alignments_batch.sh'),'--batch-id','test',
        '--manifest','manifest','--workers','2','--runner-mode','local','--tree','tree.nwk',
        '--local-dir',str(SRC)],cwd=tmp_path,text=True,capture_output=True)
    assert result.returncode != 0 and 'bad' in result.stderr
    assert json.loads((tmp_path/'bad.prep.json').read_text())['status']=='technical_failure'
    assert json.loads((tmp_path/'good.prep.json').read_text())['status']=='prepared'


def test_literal_species_ids_and_documented_header(tmp_path):
    p=tmp_path/'pool.tsv';p.write_text('species\tgroup\nNA\t1\nB\t0\n')
    assert load_species_sets(p)==({'NA'},{'B'})


def test_standalone_helpers_share_the_preparation_rules(tmp_path):
    raw=tmp_path/'raw.phy';raw.write_text('3 4\nA  ACDE\nB  XXXX\nC  AAAA\n')
    tree=tmp_path/'tree.nwk';tree.write_text('(A:1,B:1);')
    converted=tmp_path/'converted.fa'
    subprocess.run([sys.executable,'-B',str(SRC/'phylip_to_fasta.py'),str(raw),str(converted)],check=True)
    output=tmp_path/'filtered.fa'
    subprocess.run([sys.executable,'-B',str(SRC/'filter_fasta_to_tree.py'),
        '--tree',str(tree),'--fasta',str(converted),'--output',str(output)],check=True)
    assert output.read_text()=='>A\nACDE\n'
