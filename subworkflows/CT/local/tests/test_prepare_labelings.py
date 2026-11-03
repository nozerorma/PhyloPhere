"""The null replay must keep each cycle's own complete FOP design."""
import csv
import json
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts' / 'prepare_labelings.py'


def test_variable_hypotheses_and_strict_policy(tmp_path):
    design, source, run = (tmp_path / x for x in ('design', 'source', 'run'))
    for folder in (design, source, run):
        folder.mkdir()
    for hyp in ('H1', 'H2'):
        (design / f'traitfile_{hyp}.tab').write_text('a\t1\t1\nb\t0\t1\n')
    (design / 'contrast_hypotheses_pairs.tsv').write_text(
        'hypothesis_id\tpair\tspecies1\tspecies2\tpss_score\n'
        'H1\t1\ta\tb\t1\nH2\t1\ta\tb\t1\n')
    (source / 'fop_labelings.tab').write_text(
        'b_1~H1\ta\tb\nb_2~H1\ta\tb\nb_2~H2\ta\tb\nb_2~H3\ta\tb\n')
    (source / 'fop_pairs.tsv').write_text(
        'cycle\thypothesis_id\tpair\tspecies1\tspecies2\tpss_score\n'
        'b_1\tH1\t1\ta\tb\t1\nb_2\tH1\t1\ta\tb\t1\n'
        'b_2\tH2\t1\ta\tb\t2\nb_2\tH3\t1\ta\tb\t3\n')
    (source / 'permulation_provenance.json').write_text(json.dumps({'match_fop_applied': False}))
    base = [sys.executable, str(SCRIPT), '--source', str(source), '--config', str(design),
            '--n', '2', '--block-size', '2', '--multi', 'true']
    proc = subprocess.run(base + ['--match-fop', 'false'], cwd=run, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    with (run / 'blocks/blk_0001/fop_pairs.tsv').open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    assert [(r['cycle'], r['hypothesis_id'], r['pss_score']) for r in rows] == [
        ('b_1', 'H1', '1.0'), ('b_2', 'H1', '1.0'), ('b_2', 'H2', '2.0'), ('b_2', 'H3', '3.0')]
    assert json.loads((run / 'ct_preparation.json').read_text())['match_fop'] is False
    proc = subprocess.run(base + ['--match-fop', 'true'], cwd=run, capture_output=True, text=True)
    assert proc.returncode != 0
    assert 'does not match active hypotheses' in proc.stderr
    # Legacy harvests cannot silently claim that their hypothesis counts were unrestricted.
    (source / 'permulation_provenance.json').unlink()
    proc = subprocess.run(base + ['--match-fop', 'false'], cwd=run, capture_output=True, text=True)
    assert proc.returncode != 0
    assert 'require harvest provenance' in proc.stderr
    proc = subprocess.run(base + ['--match-fop', 'false', '--allow-relaxed', 'true'],
                          cwd=run, capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    assert 'unverified free hypothesis-count policy' in json.loads(
        (run / 'ct_preparation.json').read_text())['relaxed_reasons'][0]
    # Missing pair weights and missing hypothesis IDs must not silently change pooling.
    (source / 'permulation_provenance.json').write_text(json.dumps({'match_fop_applied': False}))
    pair_file = source / 'fop_pairs.tsv'
    pair_file.write_text(pair_file.read_text().replace('b_2\tH3\t1\ta\tb\t3\n', ''))
    proc = subprocess.run(base + ['--match-fop', 'false'], cwd=run, capture_output=True, text=True)
    assert proc.returncode != 0
    assert 'does not cover active labelings' in proc.stderr
    (source / 'fop_labelings.tab').write_text('b_1~H1\ta\tb\nb_2~H1\ta\tb\nb_2~H3\ta\tb\n')
    proc = subprocess.run(base + ['--match-fop', 'false'], cwd=run, capture_output=True, text=True)
    assert proc.returncode != 0
    assert 'noncontiguous hypothesis IDs' in proc.stderr
