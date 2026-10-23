"""Global association: AUC vs the compatible background, composition, LOO, explicit no-calibration."""
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'subworkflows/ENRICHMENT/local/src'
REGISTRY = SRC / 'posenrich_registry.py'
GLOBAL = SRC / 'posenrich_global.py'

sys.path.insert(0, str(ROOT / 'tests'))
from test_posenrich_registry import read_tsv, write  # noqa: E402
from test_posenrich_localization import make_localization_fixture  # noqa: E402


def run_global(registry_dir, out_dir, extra=()):
    result = subprocess.run([sys.executable, str(GLOBAL), '--registry_dir', str(registry_dir),
                             '--output_dir', str(out_dir), *extra], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return out_dir


def test_global_auc_composition_loo_and_no_calibration(tmp_path):
    registry_dir = make_localization_fixture(tmp_path)
    out = run_global(registry_dir, tmp_path / 'global')
    rows = {(r['collection_id'], r['term_id']): r for r in read_tsv(out / 'global_results.tsv')}

    term_a = rows[('custom', 'TERM_A')]
    # Members: G1{0,1} = 0.9, 0.8; G2{3,4} = 0.25, 0.2. Background: other 8 positions.
    # Exceedances 8 + 8 + 4 + 3, plus half credit for the 0.2-0.2 tie: 23.5/32.
    assert term_a['auc'] == repr(23.5 / 32)
    assert term_a['n_members'] == '4'
    assert term_a['n_background'] == '8'
    assert term_a['n_genes'] == '2'
    assert term_a['status'] == 'descriptive'
    assert term_a['p'] == 'NA'
    assert term_a['calibration_status'] == 'not_available'
    # Leave-one-gene-out: without G1 the AUC is G2's 0.25; without G2 it is G1's 1.0.
    assert term_a['loo_auc_min'] == repr(0.25)
    assert term_a['loo_auc_max'] == repr(1.0)

    by_gene = {(r['term_id'], r['Gene']): r for r in read_tsv(out / 'global_by_gene.tsv')}
    assert by_gene[('TERM_A', 'G1')]['auc_without_gene'] == repr(0.25)
    assert by_gene[('TERM_A', 'G2')]['auc_without_gene'] == repr(1.0)

    # A whole-protein term still gets a descriptive global AUC (association is a
    # different question from within-protein contrast).
    term_full = rows[('custom', 'TERM_FULL')]
    assert term_full['status'] == 'descriptive'
    assert term_full['auc'] != 'NA'

    manifest = json.loads((out / 'global_manifest.json').read_text())
    assert manifest['calibration']['p_localization_reused'] is False
    assert any('no contraste dentro de proteina' in l for l in manifest['limits'])


def test_global_members_outside_background_are_audited(tmp_path):
    # Without MAP, custom members are accepted but have no score: audited, not dropped.
    bg = write(tmp_path / 'bg.txt', 'G1\n')
    scores = write(tmp_path / 'position_scores.tsv',
                   'Gene\tPosition\tCAAS_score\tside\tp.emp\tp.adj_bh\nG1\t0\t0.9\ttop\t0.01\t0.2\n')
    custom = write(tmp_path / 'custom.tsv', 'G1\t0\tTERM\tmember\nG1\t9\tTERM\tno score\n')
    registry_dir = tmp_path / 'registry'
    result = subprocess.run([sys.executable, str(REGISTRY), '--cleaned_background', bg,
                             '--position_scores', scores, '--custom_marker_file', custom,
                             '--output_dir', str(registry_dir)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    out = run_global(registry_dir, tmp_path / 'global')
    rows = {(r['collection_id'], r['term_id']): r for r in read_tsv(out / 'global_results.tsv')}
    term = rows[('custom', 'TERM')]
    assert term['n_members_total'] == '2'
    assert term['n_members'] == '1'  # only the scored position enters the contrast
    assert term['status'] == 'not_identifiable'  # one member, no background


def test_global_minus_log10_pemp_sensitivity(tmp_path):
    registry_dir = make_localization_fixture(tmp_path)
    out = run_global(registry_dir, tmp_path / 'global', extra=['--ranking', 'minus_log10_pemp'])
    rows = {(r['collection_id'], r['term_id']): r for r in read_tsv(out / 'global_results.tsv')}
    term_a = rows[('custom', 'TERM_A')]
    # p.emp 0.01 -> 2.0 for G1 members {0,1}; p.emp 0.2 -> 0.699 for the other 10
    # positions (8 background + G2 members {3,4}). Average ranks: 2.0 -> 11.5,
    # 0.699 -> 5.5. Member rank sum = 2*11.5 + 2*5.5 = 34; AUC = (34 - 10) / 32 = 0.75.
    assert float(term_a['auc']) == 0.75
