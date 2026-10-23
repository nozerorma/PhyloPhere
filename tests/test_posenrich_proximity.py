"""Proximity: protein distances in reference residues, genomic in nt, separate contracts."""
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'subworkflows/ENRICHMENT/local/src'
REGISTRY = SRC / 'posenrich_registry.py'
PROXIMITY = SRC / 'posenrich_proximity.py'

sys.path.insert(0, str(ROOT / 'tests'))
from test_posenrich_registry import MAP_TSV, read_tsv, write  # noqa: E402


def make_proximity_fixture(tmp_path):
    """G1: cols 0 (residue 10, chr1:1000000), 1 (no residue, no nt), 2 (residue 12,
    chr1:1000005). A custom SITE at col 2; a Pfam region covering col 2."""
    maps = tmp_path / 'maps'
    maps.mkdir()
    write(maps / 'G1.map.tsv', MAP_TSV)
    bg = write(tmp_path / 'bg.txt', 'G1\n')
    scores = write(tmp_path / 'position_scores.tsv',
                   'Gene\tPosition\tCAAS_score\tside\tp.emp\tp.adj_bh\n'
                   'G1\t0\t0.9\ttop\t0.01\t0.2\n'
                   'G1\t1\t0.8\ttop\t0.05\t0.4\n'
                   'G1\t2\t0.7\ttop\t0.2\t0.9\n')
    domains = write(tmp_path / 'domains.tsv',
                    'gene\tpfam_id\ttarget_name\tdescription\tclan_acc\tclan_name\tali_start\tali_end\tcoordinate_space\n'
                    'G1\tPF1\tD1\tdesc\tNA\tNA\t12\t12\treference_protein_1based\n')
    custom = write(tmp_path / 'custom.tsv', 'G1\t2\tSITE\tpoint site\n')
    out = tmp_path / 'registry'
    result = subprocess.run([sys.executable, str(REGISTRY), '--map_dir', str(maps),
                             '--cleaned_background', bg, '--position_scores', scores,
                             '--domain_variability_file', domains,
                             '--custom_marker_file', custom, '--output_dir', str(out)],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return out


def run_proximity(registry_dir, out_dir, extra=()):
    result = subprocess.run([sys.executable, str(PROXIMITY), '--registry_dir', str(registry_dir),
                             '--output_dir', str(out_dir), *extra], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return out_dir


def test_protein_distances_and_regions(tmp_path):
    registry_dir = make_proximity_fixture(tmp_path)
    out = run_proximity(registry_dir, tmp_path / 'prox')
    rows = {(r['collection_id'], r['term_id'], int(r['Position'])): r
            for r in read_tsv(out / 'proximity_protein_distances.tsv')}

    # Point site at col 2 (residue 12): col 0 (residue 10) is 2 residues away.
    site = rows[('custom', 'SITE', 0)]
    assert site['distance_residues'] == '2' and site['in_region'] == '0'
    assert site['status'] == 'measured'
    # On the site itself: 0.
    assert rows[('custom', 'SITE', 2)]['distance_residues'] == '0'
    # Position without reference residue: coverage lack, not infinite distance.
    nores = rows[('custom', 'SITE', 1)]
    assert nores['status'] == 'no_reference_residue'
    assert nores['distance_residues'] == 'NA'
    # Region (Pfam domain covering col 2): member position is inside (0, in_region=1);
    # col 0 is 2 residues from the region.
    region_member = rows[('pfam_domains', 'PF1', 2)]
    assert region_member['distance_residues'] == '0' and region_member['in_region'] == '1'
    region_near = rows[('pfam_domains', 'PF1', 0)]
    assert region_near['distance_residues'] == '2' and region_near['in_region'] == '0'

    summary = {(r['collection_id'], r['term_id']): r
               for r in read_tsv(out / 'proximity_protein_summary.tsv')}
    site_sum = summary[('custom', 'SITE')]
    assert site_sum['n_positions_evaluated'] == '3'
    assert site_sum['n_positions_measured'] == '2'
    assert site_sum['coverage_proportion'] == repr(2 / 3)
    assert site_sum['prop_within_0'] == repr(1 / 2)
    assert site_sum['prop_within_5'] == '1.0'
    assert site_sum['status'] == 'descriptive'
    assert site_sum['control_status'] == 'not_validated'


def test_genomic_distances_separate_contract(tmp_path):
    registry_dir = make_proximity_fixture(tmp_path)
    out = run_proximity(registry_dir, tmp_path / 'prox')
    rows = {(r['collection_id'], r['term_id'], int(r['Position'])): r
            for r in read_tsv(out / 'proximity_genomic_distances.tsv')}
    # Site at chr1:1000005: col 0 (chr1:1000000) is 5 nt away; col 2 coincides.
    assert rows[('custom', 'SITE', 0)]['distance_nt'] == '5'
    assert rows[('custom', 'SITE', 2)]['distance_nt'] == '0'
    # Position without genomic coordinate: coverage lack.
    assert rows[('custom', 'SITE', 1)]['status'] == 'no_genomic_coordinate'
    summary = {(r['collection_id'], r['term_id']): r
               for r in read_tsv(out / 'proximity_genomic_summary.tsv')}
    assert summary[('custom', 'SITE')]['median_distance'] == '2.5'
    assert summary[('custom', 'SITE')]['control_status'] == 'not_validated'


def test_proximity_without_map_reports_coverage_not_infinite(tmp_path):
    bg = write(tmp_path / 'bg.txt', 'G1\n')
    scores = write(tmp_path / 'position_scores.tsv',
                   'Gene\tPosition\tCAAS_score\tside\tp.emp\tp.adj_bh\nG1\t0\t0.9\ttop\t0.01\t0.2\n')
    custom = write(tmp_path / 'custom.tsv', 'G1\t0\tSITE\tpoint site\n')
    registry_dir = tmp_path / 'registry'
    result = subprocess.run([sys.executable, str(REGISTRY), '--cleaned_background', bg,
                             '--position_scores', scores, '--custom_marker_file', custom,
                             '--output_dir', str(registry_dir)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    out = run_proximity(registry_dir, tmp_path / 'prox')
    rows = read_tsv(out / 'proximity_protein_distances.tsv')
    assert rows[0]['status'] == 'no_reference_residue'
    assert rows[0]['distance_residues'] == 'NA'
    genomic = read_tsv(out / 'proximity_genomic_distances.tsv')
    assert genomic[0]['status'] == 'no_genomic_coordinate'
    manifest = json.loads((out / 'proximity_manifest.json').read_text())
    assert manifest['control']['status'] == 'not_validated'
    assert 'proximidad genómica y proteica no se mezclan' in manifest['limits']
