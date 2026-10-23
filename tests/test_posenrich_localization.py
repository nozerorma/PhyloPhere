"""Within-protein localization: manual AUC oracle, conditional control, identifiability states."""
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'subworkflows/ENRICHMENT/local/src'
REGISTRY = SRC / 'posenrich_registry.py'
LOCALIZATION = SRC / 'posenrich_localization.py'

sys.path.insert(0, str(ROOT / 'tests'))
from test_posenrich_registry import read_tsv, write  # noqa: E402


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


loc = load_module(LOCALIZATION, 'posenrich_localization')


def make_localization_fixture(tmp_path, custom_text=None, scores_text=None):
    """Two genes x six positions; TERM_A members hold the high scores of G1 and the
    low scores of G2; TERM_FULL covers a whole protein; TERM_ONE has one member."""
    maps = tmp_path / 'maps'
    maps.mkdir(parents=True, exist_ok=True)
    for gene, chrom in (('G1', 'chr1'), ('G2', 'chr2')):
        lines = ['ori_codon_col\tstatus\tprot_ali_col\thg38_aa_pos\thg38_nt_coord']
        for i in range(6):
            lines.append(f'{i + 1}\tselected\t{i + 1}\t{i + 1}\t{chrom}:{1000 * (i + 1)}')
        write(maps / f'{gene}.map.tsv', '\n'.join(lines) + '\n')
    bg = write(tmp_path / 'bg.txt', 'G1\nG2\n')
    if scores_text is None:
        rows = ['Gene\tPosition\tCAAS_score\tside\tp.emp\tp.adj_bh\tnd\tnull_unseen\tcaas_nominal\tcaas_adjusted']
        g1_scores = [0.9, 0.8, 0.7, 0.2, 0.1, 0.05]
        g2_scores = [0.5, 0.4, 0.3, 0.25, 0.2, 0.1]
        for gene, scores in (('G1', g1_scores), ('G2', g2_scores)):
            for pos, score in enumerate(scores):
                p = '0.01' if (gene == 'G1' and pos < 2) else '0.2'
                rows.append(f'{gene}\t{pos}\t{score}\ttop\t{p}\t0.5\t2\tFALSE\tFALSE\tFALSE')
        scores_text = '\n'.join(rows) + '\n'
    scores = write(tmp_path / 'position_scores.tsv', scores_text)
    if custom_text is None:
        custom_text = ('G1\t0\tTERM_A\tlocalized high\nG1\t1\tTERM_A\tlocalized high\n'
                       'G2\t3\tTERM_A\tlocalized low\nG2\t4\tTERM_A\tlocalized low\n'
                       + ''.join(f'G1\t{i}\tTERM_FULL\twhole protein\n' for i in range(6))
                       + 'G2\t5\tTERM_ONE\tsingle member\n')
    custom = write(tmp_path / 'custom.tsv', custom_text)
    out = tmp_path / 'registry'
    result = subprocess.run([sys.executable, str(REGISTRY), '--map_dir', str(maps),
                             '--cleaned_background', bg, '--position_scores', scores,
                             '--custom_marker_file', custom, '--output_dir', str(out)],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return out


def run_localization(registry_dir, out_dir, extra=()):
    result = subprocess.run([sys.executable, str(LOCALIZATION), '--registry_dir', str(registry_dir),
                             '--output_dir', str(out_dir), '--n_permutations', '199',
                             '--seed', '7', *extra], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return out_dir


def test_auc_oracle_with_ties():
    import numpy as np
    scores = np.array([0.9, 0.8, 0.8, 0.1])
    ranks = loc.average_ranks(scores)
    assert list(ranks) == [4.0, 2.5, 2.5, 1.0]
    mask = np.array([1, 1, 0, 0], dtype=float)
    # (6.5 - 3) / 4 = 0.875: half credit to the tie at 0.8
    assert loc.auc_from_ranks(ranks, mask) == 0.875
    # all ties -> 0.5
    ranks = loc.average_ranks(np.array([0.5, 0.5, 0.5, 0.5]))
    assert loc.auc_from_ranks(ranks, np.array([1, 1, 0, 0], dtype=float)) == 0.5
    # one group only -> None (gene not informative)
    assert loc.auc_from_ranks(ranks, np.array([1, 1, 1, 1], dtype=float)) is None


def test_localization_effect_aggregation_and_control(tmp_path):
    registry_dir = make_localization_fixture(tmp_path)
    out = run_localization(registry_dir, tmp_path / 'loc')
    rows = {(r['collection_id'], r['term_id']): r for r in read_tsv(out / 'localization_results.tsv')}

    term_a = rows[('custom', 'TERM_A')]
    # G1: members {0,1} hold the top scores -> AUC 1.0; G2: members {3,4} hold low
    # scores -> 2/8 = 0.25. Equal-weight mean = 0.625; pair weights tie -> same.
    assert term_a['status'] == 'executed'
    assert term_a['auc'] == repr(0.625)
    assert term_a['auc_mean'] == repr(0.625)
    assert term_a['auc_pair_weighted'] == repr(0.625)  # pair weights tie here
    assert term_a['auc_minus_half'] == repr(0.125)
    assert term_a['n_genes_informative'] == '2'
    assert term_a['n_genes_evaluated'] == '2'
    assert term_a['n_members'] == '4'
    assert term_a['n_background_positions'] == '8'
    assert term_a['B_valid'] == '199'
    assert 0 < float(term_a['p_mc']) <= 1
    assert float(term_a['auc_control_mean']) != float(term_a['auc'])

    by_gene = {(r['term_id'], r['Gene']): r for r in read_tsv(out / 'localization_by_gene.tsv')}
    assert by_gene[('TERM_A', 'G1')]['auc'] == '1.0'
    assert by_gene[('TERM_A', 'G2')]['auc'] == '0.25'

    # A term covering a whole protein has no within-protein contrast.
    term_full = rows[('custom', 'TERM_FULL')]
    assert term_full['status'] == 'not_identifiable'
    assert term_full['coverage_status'] == 'no_within_protein_contrast'
    # A single low-scoring member: informative gene, AUC 0, p_mc = 1 (minimum).
    term_one = rows[('custom', 'TERM_ONE')]
    assert term_one['auc'] == '0.0'
    assert term_one['p_mc'] == '1.0'
    # Genomic windows are marked non-informative for this analysis and stay out.
    assert all(r['collection_id'] != 'genomic_windows_1mb' for r in rows.values())


def test_localization_all_ties_gives_p_one(tmp_path):
    scores_rows = ['Gene\tPosition\tCAAS_score\tside\tp.emp\tp.adj_bh']
    for pos in range(6):
        scores_rows.append(f'G1\t{pos}\t0.5\ttop\t0.2\t0.5')
    custom = ''.join(f'G1\t{i}\tFLAT\tequal scores\n' for i in (0, 1))
    registry_dir = make_localization_fixture(tmp_path, custom_text=custom,
                                            scores_text='\n'.join(scores_rows) + '\n')
    out = run_localization(registry_dir, tmp_path / 'loc')
    rows = {(r['collection_id'], r['term_id']): r for r in read_tsv(out / 'localization_results.tsv')}
    flat = rows[('custom', 'FLAT')]
    assert flat['auc'] == '0.5'
    assert flat['p_mc'] == '1.0'  # every permutation ties the observed AUC


def test_localization_reproducible_and_order_inert(tmp_path):
    registry_dir = make_localization_fixture(tmp_path)
    out1 = run_localization(registry_dir, tmp_path / 'loc1')
    out2 = run_localization(registry_dir, tmp_path / 'loc2')
    assert (out1 / 'localization_results.tsv').read_text() == \
        (out2 / 'localization_results.tsv').read_text()
    # Same members in a different file order: identical results.
    custom = ('G2\t4\tTERM_A\tlocalized low\nG1\t1\tTERM_A\tlocalized high\n'
              'G2\t3\tTERM_A\tlocalized low\nG1\t0\tTERM_A\tlocalized high\n')
    registry_dir_b = make_localization_fixture(tmp_path / 'b', custom_text=custom)
    out3 = run_localization(registry_dir_b, tmp_path / 'loc3')
    rows1 = {(r['collection_id'], r['term_id']): r for r in read_tsv(out1 / 'localization_results.tsv')}
    rows3 = {(r['collection_id'], r['term_id']): r for r in read_tsv(out3 / 'localization_results.tsv')}
    assert rows1[('custom', 'TERM_A')]['auc'] == rows3[('custom', 'TERM_A')]['auc']
    assert rows1[('custom', 'TERM_A')]['p_mc'] == rows3[('custom', 'TERM_A')]['p_mc']


def test_localization_minus_log10_pemp_sensitivity(tmp_path):
    # p.emp shared between sides; positions without p are not evaluable.
    rows = ['Gene\tPosition\tCAAS_score\tside\tp.emp\tp.adj_bh']
    g1 = [0.9, 0.8, 0.7, 0.2, 0.1, 0.05]
    pemp = ['0.01', '0.01', '0.1', '0.2', '0.5', 'NA']
    for pos, (score, p) in enumerate(zip(g1, pemp)):
        rows.append(f'G1\t{pos}\t{score}\ttop\t{p}\t0.5')
        rows.append(f'G1\t{pos}\t{score}\tbottom\t{p}\t0.5')  # shared p, directional score
    scores_text = '\n'.join(rows) + '\n'
    custom = 'G1\t0\tTERM_A\thigh\nG1\t1\tTERM_A\thigh\n'
    registry_dir = make_localization_fixture(tmp_path, custom_text=custom, scores_text=scores_text)
    out = run_localization(registry_dir, tmp_path / 'loc', extra=['--ranking', 'minus_log10_pemp'])
    results = {(r['collection_id'], r['term_id']): r for r in read_tsv(out / 'localization_results.tsv')}
    term_a = results[('custom', 'TERM_A')]
    # Members hold p.emp = 0.01 -> score 2.0, the maximum; position 5 (p NA) is not evaluable.
    assert term_a['auc'] == '1.0'
    by_gene = read_tsv(out / 'localization_by_gene.tsv')
    assert by_gene[0]['n_positions_evaluable'] == '5'


def test_control_calibration_under_no_effect():
    # Synthetic null: iid scores unrelated to membership, heterogeneous genes.
    # The conditional control must not manufacture significance: p_mc stays
    # approximately uniform (loose bound; 30 replicas, B=199).
    import numpy as np
    rng = np.random.default_rng(0)
    ps = []
    for rep in range(30):
        genes = {}
        for g in range(40):
            scores = rng.random(20)
            members = set(rng.choice(20, size=5, replace=False).tolist())
            genes[f'G{g}'] = {'positions': {i: float(scores[i]) for i in range(20)},
                              'member_positions': members, 'n_members_evaluable': 5}
        term_data = {'genes': genes, 'n_members_total': 200, 'description': 'synthetic'}
        summary, per_gene, n_inf = loc.run_term(term_data, 'mean', 199,
                                                np.random.default_rng(100 + rep))
        assert n_inf == 40
        ps.append(summary['p_mc'])
    frac = sum(p < 0.05 for p in ps) / len(ps)
    assert frac <= 0.2, f'control descalibrado: {frac} de p_mc < 0.05 bajo nulo sintético'
    assert all(0 < p <= 1 for p in ps)


def test_localization_no_hit_selection_before_auc(tmp_path):
    # Members with the LOWEST scores: the background is every evaluable position,
    # not a pre-selected hit list.
    rows = ['Gene\tPosition\tCAAS_score\tside\tp.emp\tp.adj_bh']
    for pos, score in enumerate([0.9, 0.8, 0.7, 0.2, 0.1, 0.05]):
        rows.append(f'G1\t{pos}\t{score}\ttop\t0.2\t0.5')
    scores_text = '\n'.join(rows) + '\n'
    custom = 'G1\t4\tLOW\tlow member\nG1\t5\tLOW\tlow member\n'
    registry_dir = make_localization_fixture(tmp_path, custom_text=custom, scores_text=scores_text)
    out = run_localization(registry_dir, tmp_path / 'loc')
    results = {(r['collection_id'], r['term_id']): r for r in read_tsv(out / 'localization_results.tsv')}
    assert results[('custom', 'LOW')]['auc'] == '0.0'
    assert results[('custom', 'LOW')]['p_mc'] == '1.0'
