"""Common position registry: evaluation states, members, coverage, manifest and GMT export."""
import csv
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'subworkflows/ENRICHMENT/local/src'
REGISTRY = SRC / 'posenrich_registry.py'
BUILDER = SRC / 'build_position_gmt.py'

MAP_TSV = ('ori_codon_col\tstatus\tprot_ali_col\thg38_aa_pos\thg38_nt_coord\n'
           '1\tselected\t1\t10\tchr1:1000000\n'
           '2\tremoved\tNA\t11\tchr1:1000002\n'
           '3\tselected\t2\tNA\tNA\n'
           '4\tselected\t3\t12\tchr1:1000005\n')

SCORES_TSV = ('Gene\tPosition\tCAAS_score\tside\tp.emp\tp.adj_bh\tp.emp_fact\tp.adj_bh_fact\tnd\tnull_unseen\tcaas_nominal\tcaas_adjusted\n'
              'G1\t0\t0.9\ttop\t0.01\t0.2\t0.001\t0.05\t11\tFALSE\tTRUE\tFALSE\n'
              'G1\t0\t0.8\tbottom\t0.01\t0.2\t0.001\t0.05\t11\tFALSE\tTRUE\tFALSE\n'
              'G1\t2\t0.5\ttop\t0.2\t0.9\t0.01\t0.2\t3\tFALSE\tFALSE\tFALSE\n'
              'OUT\t5\t0.7\ttop\t0.04\t0.3\t0.004\t0.1\t2\tFALSE\tTRUE\tFALSE\n')

DISCOVERY_TSV = ('Gene\tPosition\tcaas\tconvergence_type\ttag_support\tcaas_support\tparticipating_hypotheses\tn_hypotheses\tn_species\tside\tasr_path_score\n'
                 'G1\t0\tF/L\tconvergent\tCAAS_x:1\tFFF/FLL:1\tH1\t1\t2.0\ttop\t0.9\n'
                 'G1\t0\tF/L\tconvergent\tCAAS_x:1\tFFF/FLL:1\tH1\t1\t2.0\tbottom\t0.8\n'
                 'G1\t2\tF/L\tconvergent\tCAAS_x:1\tFFF/FLL:1\tH1\t1\t2.0\ttop\t0.5\n'
                 'G1\t3\tF/L\tconvergent\tCAAS_x:1\tFFF/FLL:1\tH1\t1\t2.0\ttop\t0.0\n'
                 'OUT\t7\tF/L\tconvergent\tCAAS_x:1\tFFF/FLL:1\tH1\t1\t2.0\ttop\t0.6\n')

NULL_UNIVERSE_TSV = ('metric\tcycles\tvalue\tnote\n'
                     'N_cycles\tNA\t100\tcycle_roster.txt\n'
                     'null_positive_positions_in_universe\tNA\t50\t\n'
                     'observed_positive_positions\tNA\t3\t\n'
                     'observed_unseen_positions\tNA\t1\t\n'
                     'universe_size\tNA\t52\t\n')

DOMAINS_TSV = ('gene\tpfam_id\ttarget_name\tdescription\tclan_acc\tclan_name\tali_start\tali_end\tcoordinate_space\n'
               'G1\tPF1\tD1\tdesc\tCL1\tC1\t11\t12\treference_protein_1based\n'
               'G1\tPF2\tD2\tdesc\tNA\tNA\t12\t12\treference_protein_1based\n')

UCR_TSV = ('gene\tposition\tregion_type\tmethod\tvariability\n'
           'G1\t4\tcore\tsliding\t0.1\n'
           'G1\t1\tflank_up\tsliding\t0.2\n'
           'G1\t99\tcore\tsliding\t0.3\n'
           'G1\t4\tcore\tabsolute\t0.1\n')

FUBAR_TSV = ('gene\thg38_aa_pos\tis_pos_hit\tis_neg_hit\n'
             'G1\t10\t1\t0\n'
             'G1\t12\t0\t1\n'
             'G1\t77\t0\t0\n')

FADE_TOP_CSV = 'gene,position,max_bf,target_aa\nG1,0,100,S\nG1,5,100,T\n'
FADE_BOTTOM_CSV = 'gene,position,max_bf,target_aa\nG1,2,150,T\n'
CUSTOM_TSV = 'G1\t0\tPP2A\tvalidated mapping fixture\n'


def write(path, content):
    path.write_text(content)
    return str(path)


def make_fixture(tmp_path):
    maps = tmp_path / 'maps'
    maps.mkdir()
    write(maps / 'G1.map.tsv', MAP_TSV)
    files = {
        'background': write(tmp_path / 'bg.txt', 'G1\n'),
        'scores': write(tmp_path / 'position_scores.tsv', SCORES_TSV),
        'postproc': write(tmp_path / 'filtered_discovery.tsv', DISCOVERY_TSV),
        'null_universe': write(tmp_path / 'null_universe_summary.tsv', NULL_UNIVERSE_TSV),
        'domains': write(tmp_path / 'domains.tsv', DOMAINS_TSV),
        'ucr': write(tmp_path / 'ucr.tsv', UCR_TSV),
        'fubar': write(tmp_path / 'fubar.tsv', FUBAR_TSV),
        'fade_top': write(tmp_path / 'fade_top.csv', FADE_TOP_CSV),
        'fade_bottom': write(tmp_path / 'fade_bottom.csv', FADE_BOTTOM_CSV),
        'custom': write(tmp_path / 'custom.tsv', CUSTOM_TSV),
    }
    return files, str(maps)


def run_registry(tmp_path, files, maps, extra=()):
    out = tmp_path / 'registry'
    gmt = tmp_path / 'gmt_export'
    cmd = [sys.executable, str(REGISTRY),
           '--map_dir', maps, '--cleaned_background', files['background'],
           '--position_scores', files['scores'], '--postproc', files['postproc'],
           '--null_universe', files['null_universe'],
           '--domain_variability_file', files['domains'], '--ucr_positions_file', files['ucr'],
           '--fubar_sites_file', files['fubar'], '--fade_sites_top_file', files['fade_top'],
           '--fade_sites_bottom_file', files['fade_bottom'], '--custom_marker_file', files['custom'],
           '--output_dir', str(out), '--gmt_dir', str(gmt), *extra]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return out, gmt


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def test_registry_states_members_coverage_manifest(tmp_path):
    files, maps = make_fixture(tmp_path)
    out, gmt = run_registry(tmp_path, files, maps)

    positions = {(r['Gene'], int(r['Position'])): r for r in read_tsv(out / 'positions.tsv')}
    assert set(positions) == {('G1', 0), ('G1', 1), ('G1', 2), ('G1', 3)}
    assert positions[('G1', 0)]['evaluation_status'] == 'detected_positive'
    assert positions[('G1', 2)]['evaluation_status'] == 'detected_positive'
    # In discovery but without positive score: evaluated without detection, not a zero score.
    assert positions[('G1', 3)]['evaluation_status'] == 'evaluated_no_detection'
    assert positions[('G1', 3)]['exclusion_cause'] == 'no_map'
    # MAP column without CT detection: evaluation unknown, not a negative.
    assert positions[('G1', 1)]['evaluation_status'] == 'map_only_not_in_discovery'
    # Coordinates come from explicit MAP headers, never from row counting.
    assert positions[('G1', 0)]['ori_codon_col'] == '1'
    assert positions[('G1', 0)]['hg38_aa_pos'] == '10'
    assert positions[('G1', 2)]['hg38_aa_pos'] == '12'
    assert positions[('G1', 2)]['hg38_nt_coord'] == 'chr1:1000005'

    evidence = read_tsv(out / 'position_evidence.tsv')
    assert len(evidence) == 3  # the OUT row is outside the universe
    g10 = [r for r in evidence if r['Position'] == '0']
    assert {r['side'] for r in g10} == {'top', 'bottom'}
    # p.emp is the pooled max-over-sides p: one shared test, not two directional ones.
    assert len({r['p.emp'] for r in g10}) == 1
    assert all(r['p_shared'] == 'TRUE' for r in evidence)
    assert all(r['N_cycles'] == '100' for r in evidence)
    assert all(r['evidence_status'] == 'evaluated' for r in evidence)
    # Directional scores stay separate from the shared p.
    assert {r['CAAS_score'] for r in g10} == {'0.9', '0.8'}

    members = {(r['collection_id'], r['term_id'], r['Gene'], int(r['Position']))
               for r in read_tsv(out / 'members.tsv')}
    assert ('pfam_domains', 'PF1', 'G1', 2) in members
    assert ('pfam_domains', 'PF2', 'G1', 2) in members
    assert ('pfam_clans', 'CL1', 'G1', 2) in members
    assert ('ucr_core', 'global', 'G1', 2) in members
    assert ('ucr_flank', 'global', 'G1', 0) in members
    assert ('fubar_positive', 'global', 'G1', 0) in members
    assert ('fubar_purifying', 'global', 'G1', 2) in members
    assert ('fade_top_sig', 'global', 'G1', 0) in members
    assert ('fade_bottom_sig', 'global', 'G1', 2) in members
    assert ('genomic_windows_1mb', 'chr1_0M', 'G1', 0) in members
    assert ('genomic_windows_1mb', 'chr1_1M', 'G1', 2) in members
    assert ('custom', 'PP2A', 'G1', 0) in members

    collections = {(r['collection_id'], r['term_id']): r for r in read_tsv(out / 'collections.tsv')}
    assert collections[('genomic_windows_1mb', 'chr1_1M')]['informative_within_protein'] == '0'
    assert collections[('pfam_domains', 'PF1')]['informative_within_protein'] == '1'
    assert collections[('genomic_windows_1mb', 'chr1_1M')]['analyses'] == 'proximity_genomic'
    assert collections[('custom', 'PP2A')]['status'] == 'included'
    assert collections[('pfam_domains', 'PF1')]['source_sha256'] != 'NA'

    coverage = read_tsv(out / 'coverage.tsv')
    uni = {r['metric']: r for r in coverage if r['scope'] == 'universe'}
    assert uni['n_genes']['value'] == '1'
    assert uni['n_positions_annotable']['value'] == '3'
    assert uni['n_positions_evaluated_discovery']['value'] == '3'
    assert uni['N_cycles']['value'] == '100'
    gene1 = [r for r in coverage if r['scope'] == 'gene' and r['id'] == 'G1']
    assert {r['metric']: r['value'] for r in gene1}['n_map_positions'] == '3'
    col = [r for r in coverage if r['scope'] == 'collection' and r['id'] == 'ucr_core']
    assert {r['metric']: r for r in col}['coverage_proportion']['denominator'] == '3'

    pending = read_tsv(out / 'pending_mappings.tsv')
    reasons = {(r['source'], r['gene'], r['raw_position'], r['reason']) for r in pending}
    assert ('ucr', 'G1', '99', 'raw_position_not_in_map') in reasons
    assert ('fubar', 'G1', '77', 'reference_residue_not_in_map') in reasons
    assert ('fade', 'G1', '5', 'position_outside_map') in reasons
    assert all(r['note'] for r in pending)

    manifest = json.loads((out / 'manifest.json').read_text())
    assert manifest['contract_version'] == 1
    assert manifest['positions']['by_status']['detected_positive'] == 2
    assert manifest['positions']['by_status']['evaluated_no_detection'] == 1
    assert manifest['evidence']['p_available'] == 3
    assert manifest['evidence']['null_cycles'] == 100
    assert manifest['audited_losses']['score_rows_outside_universe'] == 1
    assert manifest['audited_losses']['discovery_rows_outside_universe'] == 1
    by_id = {c['collection_id']: c for c in manifest['collections']}
    assert by_id['custom']['status'] == 'included'
    assert by_id['ucr_core']['n_members'] == 1
    assert manifest['sources']['position_scores']['sha256'] != 'NA'

    # The GMT export is byte-identical to build_position_gmt.py over the same inputs.
    ref = tmp_path / 'gmt_ref'
    result = subprocess.run([sys.executable, str(BUILDER), '--map_dir', maps,
                             '--cleaned_background', files['background'],
                             '--domain_variability_file', files['domains'],
                             '--ucr_positions_file', files['ucr'],
                             '--fubar_sites_file', files['fubar'],
                             '--fade_sites_top_file', files['fade_top'],
                             '--fade_sites_bottom_file', files['fade_bottom'],
                             '--custom_marker_file', files['custom'],
                             '--output_dir', str(ref)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    for name in ('pfam_domains.gmt', 'pfam_clans.gmt', 'genomic_locations.gmt',
                 'custom_features.gmt', 'characterization_layers.tsv'):
        assert (gmt / name).read_text() == (ref / name).read_text(), name


def test_registry_without_map_scores_or_null(tmp_path):
    bg = write(tmp_path / 'bg.txt', 'G1\nG2\n')
    scores = write(tmp_path / 'position_scores.tsv',
                   'Gene\tPosition\tCAAS_score\tside\tp.emp\tp.adj_bh\nG2\t4\t0.3\ttop\tNA\tNA\n')
    out = tmp_path / 'out'
    result = subprocess.run([sys.executable, str(REGISTRY), '--cleaned_background', bg,
                             '--position_scores', scores, '--output_dir', str(out)],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    positions = read_tsv(out / 'positions.tsv')
    assert len(positions) == 1
    assert positions[0]['evaluation_status'] == 'detected_positive'
    assert positions[0]['exclusion_cause'] == 'no_map'
    assert positions[0]['hg38_nt_coord'] == 'NA'
    evidence = read_tsv(out / 'position_evidence.tsv')
    assert evidence[0]['evidence_status'] == 'p_unavailable'
    assert evidence[0]['N_cycles'] == 'NA'
    manifest = json.loads((out / 'manifest.json').read_text())
    assert manifest['evidence']['null_cycles'] == 'NA'
    assert all(c['status'] == 'missing' for c in manifest['collections'])
    assert (out / 'members.tsv').read_text().count('\n') == 1
    assert (out / 'pending_mappings.tsv').read_text().count('\n') == 1


def test_registry_custom_pending_rows_and_malformed_file(tmp_path):
    maps = tmp_path / 'maps'
    maps.mkdir()
    write(maps / 'G1.map.tsv', MAP_TSV)
    bg = write(tmp_path / 'bg.txt', 'G1\n')
    # Invalid rows become pending mappings, never negatives; the collection stays empty.
    custom = write(tmp_path / 'custom.tsv', 'G1\t9\tBAD\toutside map\nG1\t-1\tNEG\tnegative\n')
    out = tmp_path / 'out'
    result = subprocess.run([sys.executable, str(REGISTRY), '--map_dir', str(maps),
                             '--cleaned_background', bg, '--custom_marker_file', custom,
                             '--output_dir', str(out)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    pending = read_tsv(out / 'pending_mappings.tsv')
    assert {(r['reason']) for r in pending} == {'position_outside_map', 'negative_position'}
    manifest = json.loads((out / 'manifest.json').read_text())
    custom = [c for c in manifest['collections'] if c['collection_id'] == 'custom'][0]
    assert custom['status'] == 'empty'
    # A malformed row rejects the whole file: the collection is incompatible, not silent.
    bad = write(tmp_path / 'custom_bad.tsv', 'G1\tX\tBROKEN\tnot a position\n')
    out2 = tmp_path / 'out2'
    result = subprocess.run([sys.executable, str(REGISTRY), '--map_dir', str(maps),
                             '--cleaned_background', bg, '--custom_marker_file', bad,
                             '--output_dir', str(out2)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    collections = read_tsv(out2 / 'collections.tsv')
    assert all(r['status'] == 'incompatible' for r in collections if r['collection_id'] == 'custom')
    manifest = json.loads((out2 / 'manifest.json').read_text())
    custom = [c for c in manifest['collections'] if c['collection_id'] == 'custom'][0]
    assert custom['status'] == 'incompatible' and 'note' in custom


def test_registry_sentinels_and_empty_inputs(tmp_path):
    out = tmp_path / 'out'
    result = subprocess.run([sys.executable, str(REGISTRY),
                             '--map_dir', 'NO_ANNOT_MAP', '--cleaned_background', 'NO_BACKGROUND',
                             '--position_scores', 'NO_POS_SCORES', '--output_dir', str(out)],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    manifest = json.loads((out / 'manifest.json').read_text())
    assert manifest['universe']['n_genes'] == 0
    assert manifest['positions']['total'] == 0
