"""Proximity: protein distances in reference residues, genomic in nt, separate contracts."""
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'subworkflows/POSENRICH/local/src'
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
    from posenrich_fixtures import reference, annotation
    identity=reference(maps/'G1.map.tsv','G1')
    bg = write(tmp_path / 'bg.txt', 'G1\n')
    scores = write(tmp_path / 'position_scores.tsv',
                   'Gene\tPosition\tCAAS_score\tside\tp.emp\tp.adj_bh\n'
                   'G1\t0\t0.9\ttop\t0.01\t0.2\n'
                   'G1\t1\t0.8\ttop\t0.05\t0.4\n'
                   'G1\t2\t0.7\ttop\t0.2\t0.9\n')
    domains = write(tmp_path / 'domains.tsv',
                    'gene\tpfam_id\ttarget_name\tdescription\tclan_acc\tclan_name\tali_start\tali_end\tcoordinate_space\n'
                    'G1\tPF1\tD1\tdesc\tNA\tNA\t12\t12\treference_protein_1based\n')
    annotation(domains, ['G1'], identities={'G1':identity})
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
    assert 'control_status' not in site_sum

    # Overlap by identity of the position; near = another position within the largest radius, with its closest member.
    hits = {(r['collection_id'], r['term_id'], int(r['Position'])): r for r in read_tsv(out / 'proximity_protein_hits.tsv')}
    assert hits[('custom', 'SITE', 2)]['relation'] == 'overlap' and hits[('custom', 'SITE', 2)]['distance_residues'] == '0'
    near = hits[('custom', 'SITE', 0)]
    assert near['relation'] == 'near' and near['distance_residues'] == '2'
    assert near['nearest_site_position'] == '2' and near['nearest_site_residue'] == '12'
    assert ('custom', 'SITE', 1) not in hits                      # no reference residue: no distance, no hit
    by_gene = {(r['collection_id'], r['term_id'], r['Gene']): r for r in read_tsv(out / 'proximity_protein_by_gene.tsv')}
    row = next(r for k, r in by_gene.items() if k[:2] == ('custom', 'SITE'))
    assert (row['n_detections_measured'], row['n_overlaps'], row['n_near'], row['min_distance_residues']) == ('2', '1', '1', '0')


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
    assert 'control_status' not in summary[('custom', 'SITE')]


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
    assert 'control' not in manifest
    assert 'misma proteina' in manifest['contracts']['protein']
    assert 'mismo gen y cromosoma' in manifest['contracts']['genomic']
    assert 'proximidad genómica y proteica no se mezclan' in manifest['limits']


def test_best_side_scores_and_empty_collection_filter(tmp_path):
    registry_dir = make_proximity_fixture(tmp_path)
    evidence = registry_dir / 'position_evidence.tsv'
    lines = evidence.read_text().splitlines()
    columns = lines[0].split('\t')
    duplicate = lines[1].split('\t')
    duplicate[columns.index('side')] = 'bottom'
    duplicate[columns.index('CAAS_score')] = '0.99'
    evidence.write_text('\n'.join(lines + ['\t'.join(duplicate)]) + '\n')
    from posenrich_contract import finish_manifest
    meta=json.loads((registry_dir/'manifest.json').read_text())
    finish_manifest(registry_dir,'manifest.json',meta)
    out = run_proximity(registry_dir, tmp_path / 'prox')
    rows = read_tsv(out / 'proximity_protein_distances.tsv')
    queries = [r for r in rows if r['collection_id'] == 'custom' and r['Position'] == '0']
    assert len(queries) == 1 and queries[0]['CAAS_score'] == '0.99'
    empty = run_proximity(registry_dir, tmp_path / 'empty', ['--collections', 'not_present'])
    assert read_tsv(empty / 'proximity_protein_distances.tsv') == []
    assert read_tsv(empty / 'proximity_genomic_summary.tsv') == []
    assert json.loads((empty / 'proximity_manifest.json').read_text())['states']['proximity'] == 'no_terms'


def test_streamed_summaries_match_exact_distance_statistics(tmp_path):
    import csv
    import importlib.util
    from contextlib import ExitStack

    spec = importlib.util.spec_from_file_location('proximity', PROXIMITY)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # Unsorted distances, mixed missingness, odd/even medians, and term reset.
    with ExitStack() as stack:
        sink = module.DistanceOutput(stack, str(tmp_path), 'protein', 'hg38_aa_pos',
                                     'distance_residues', [0, 5, 10])
        for term, distances in [('a', [9, 0, 2, None]), ('b', [10, 2]), ('c', [None])]:
            for pos, distance in enumerate(distances):
                sink.append(dict(collection_id='custom', term_id=term, Gene='G1', Position=pos,
                                 hg38_aa_pos=pos, distance_residues='NA' if distance is None else distance,
                                 in_region=0, CAAS_score=.9,
                                 status='measured' if distance is not None else 'no_reference_residue'))
            sink.finish_term('custom', term)
            assert len(sink.distances) == 0
    with (tmp_path / 'proximity_protein_summary.tsv').open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    assert [r['median_distance'] for r in rows] == ['2.0', '6.0', 'NA']
    assert rows[0]['mean_distance'] == repr(11 / 3)
    assert rows[0]['coverage_proportion'] == '0.75'
    assert rows[0]['prop_within_5'] == repr(2 / 3)
    assert rows[1]['n_positions_evaluated'] == '2'
    # Nearest-site lookup must respect separated intervals, not fill their gap.
    for value in range(31):
        assert module.nearest_distance(value, [3, 4, 20, 21]) == min(
            abs(value - site) for site in [3, 4, 20, 21])
