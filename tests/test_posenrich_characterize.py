"""Characterization component: annotation layers, collection/ranking summaries, explicit states."""
import csv
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'subworkflows/ENRICHMENT/local/src'
REGISTRY = SRC / 'posenrich_registry.py'
CHARACTERIZE = SRC / 'posenrich_characterize.py'

sys.path.insert(0, str(ROOT / 'tests'))
from test_posenrich_registry import make_fixture, run_registry, read_tsv, write  # noqa: E402


def run_characterize(registry_dir, out_dir):
    result = subprocess.run([sys.executable, str(CHARACTERIZE),
                             '--registry_dir', str(registry_dir), '--output_dir', str(out_dir)],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return out_dir


def test_characterization_layers_and_summaries(tmp_path):
    files, maps = make_fixture(tmp_path)
    registry_dir, _ = run_registry(tmp_path, files, maps)
    out = run_characterize(registry_dir, tmp_path / 'char')

    rows = {(r['Gene'], int(r['Position'])): r for r in read_tsv(out / 'position_characterization.tsv')}
    # Annotated or positively scored positions only.
    assert set(rows) == {('G1', 0), ('G1', 2), ('G1', 3)}
    assert rows[('G1', 3)]['CAAS_score'] == 'NA'
    g10, g12 = rows[('G1', 0)], rows[('G1', 2)]
    assert g10['pfam_domain'] == '' and g10['ucr_region'] == 'flank_up'
    assert g10['fubar_selection'] == 'positive' and g10['fade_top_sig'] == '1'
    assert g10['CAAS_score'] == '0.9' and g10['sides'] == 'bottom;top'
    assert g10['p.emp'] == '0.01' and g10['caas_nominal'] == 'TRUE'
    assert g12['pfam_domain'] == 'D1;D2' and g12['pfam_clan'] == 'C1'
    assert g12['ucr_region'] == 'core' and g12['fubar_selection'] == 'negative'
    assert g12['fade_bottom_sig'] == '1' and g12['fade_top_sig'] == '0'
    assert g12['evaluation_status'] == 'detected_positive'

    coll = {(r['collection_id'], r['term_id']): r for r in read_tsv(out / 'collection_summary.tsv')}
    assert coll[('ucr_core', 'global')]['n_members'] == '1'
    assert coll[('ucr_core', 'global')]['n_members_positive'] == '1'
    assert coll[('ucr_core', 'global')]['prop_positive_of_evaluated'] == '1.0'
    assert coll[('custom', 'PP2A')]['n_members_nominal'] == '1'
    assert coll[('custom', 'PP2A')]['n_members_adjusted'] == '0'
    # Coverage collection: every FUBAR-tested position, not only hits.
    assert coll[('fubar_tested', 'global')]['n_members'] == '2'

    ranking = {r['ranking']: r for r in read_tsv(out / 'ranking_summary.tsv')}
    assert ranking['global']['n_positions'] == '3'
    assert ranking['global']['n_scored'] == '2'
    assert ranking['global']['n_nominal'] == '1'
    assert ranking['global']['n_adjusted'] == '0'
    assert ranking['top']['n_positions'] == '2'
    assert ranking['bottom']['n_positions'] == '1'
    assert ranking['global']['N_cycles'] == '100'
    assert 'q25' in ranking['global']  # distribution without selection percentiles

    universe = read_tsv(out / 'universe_summary.tsv')
    states = {r['evaluation_status']: r['n_positions'] for r in universe}
    assert states['detected_positive'] == '2'
    assert states['evaluated_no_detection'] == '1'
    assert states['map_only_not_in_discovery'] == '1'
    assert states['coverage:n_positions_annotable'] == '3'

    manifest = json.loads((out / 'characterization_manifest.json').read_text())
    assert manifest['states']['characterization'] == 'executed'
    assert manifest['states']['localization'] == 'not_executed'
    assert 'sin inferencia de enriquecimiento' in manifest['limits']


def test_characterization_without_scores_or_annotations(tmp_path):
    # Registry with only a background: no positions, no evidence, no members.
    bg = write(tmp_path / 'bg.txt', 'G1\nG2\n')
    registry_dir = tmp_path / 'registry'
    result = subprocess.run([sys.executable, str(REGISTRY), '--cleaned_background', bg,
                             '--output_dir', str(registry_dir)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    out = run_characterize(registry_dir, tmp_path / 'char')
    assert (out / 'position_characterization.tsv').read_text().count('\n') == 1
    assert (out / 'collection_summary.tsv').read_text().count('\n') == 1
    ranking = read_tsv(out / 'ranking_summary.tsv')
    assert all(r['n_positions'] == '0' and r['n_scored'] == '0' for r in ranking)
    assert all(r['N_cycles'] == 'NA' for r in ranking)
    universe = read_tsv(out / 'universe_summary.tsv')
    assert universe[0]['evaluation_status'] == 'NA'
    manifest = json.loads((out / 'characterization_manifest.json').read_text())
    assert manifest['states']['characterization'] == 'executed'


def test_characterization_does_not_invent_flags(tmp_path):
    # Evidence without canonical flags: counts stay zero, nothing is derived silently.
    files, maps = make_fixture(tmp_path)
    scores = write(tmp_path / 'scores_noflags.tsv',
                   'Gene\tPosition\tCAAS_score\tside\tp.emp\tp.adj_bh\nG1\t0\t0.9\ttop\t0.01\t0.2\n')
    registry_dir = tmp_path / 'registry'
    result = subprocess.run([sys.executable, str(REGISTRY), '--map_dir', maps,
                             '--cleaned_background', files['background'],
                             '--position_scores', scores, '--output_dir', str(registry_dir)],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    out = run_characterize(registry_dir, tmp_path / 'char')
    rows = read_tsv(out / 'position_characterization.tsv')
    assert rows[0]['caas_nominal'] == 'NA' and rows[0]['caas_adjusted'] == 'NA'
    ranking = {r['ranking']: r for r in read_tsv(out / 'ranking_summary.tsv')}
    assert ranking['global']['n_nominal'] == '0' and ranking['global']['n_adjusted'] == '0'
