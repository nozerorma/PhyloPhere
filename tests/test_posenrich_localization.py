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
    assert term_a['status'] == 'descriptive_pending_control'
    assert term_a['auc'] == repr(0.625)
    assert term_a['auc_mean'] == repr(0.625)
    assert term_a['auc_pair_weighted'] == repr(0.625)  # pair weights tie here
    assert term_a['auc_minus_half'] == repr(0.125)
    assert term_a['n_genes_informative'] == '2'
    assert term_a['n_genes_evaluated'] == '2'
    assert term_a['n_members'] == '4'
    assert term_a['n_background_positions'] == '8'
    assert term_a['B_valid'] == '0'
    assert term_a['p_mc'] == 'NA' and term_a['q'] == 'NA'
    assert term_a['auc_control_mean'] == 'NA'
    assert term_a['ranking_population'] == 'canonical_positive_detections_only'

    by_gene = {(r['term_id'], r['Gene']): r for r in read_tsv(out / 'localization_by_gene.tsv')}
    assert by_gene[('TERM_A', 'G1')]['auc'] == '1.0'
    assert by_gene[('TERM_A', 'G2')]['auc'] == '0.25'

    # A term covering a whole protein has no within-protein contrast.
    term_full = rows[('custom', 'TERM_FULL')]
    assert term_full['status'] == 'not_identifiable'
    assert term_full['coverage_status'] == 'no_within_protein_contrast'
    # A single low-scoring member: informative gene, AUC 0; production inference remains withheld.
    term_one = rows[('custom', 'TERM_ONE')]
    assert term_one['auc'] == '0.0'
    assert term_one['p_mc'] == 'NA'
    # Genomic windows are marked non-informative for this analysis and stay out.
    assert all(r['collection_id'] != 'genomic_windows_1mb' for r in rows.values())


def test_localization_all_ties_is_descriptive(tmp_path):
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
    assert flat['p_mc'] == 'NA'  # production inference is withheld


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
    assert results[('custom', 'LOW')]['p_mc'] == 'NA'


def test_prototype_term_stream_independent_of_other_terms():
    import numpy as np
    data = {'genes': {'G': {'positions': {0: .9, 1: .8, 2: .3, 3: .1},
                          'member_positions': {0, 2}}}}
    def evaluate(term):
        return loc.run_term(data, 'mean', 199, loc.term_rng(7, 'caas_score', 'top', 'custom', term))[0]
    before = evaluate('TARGET')
    evaluate('OTHER')
    assert evaluate('TARGET') == before
    assert before['B_valid'] == 199


def test_directional_views_and_no_production_inference(tmp_path):
    scores = ('Gene\tPosition\tCAAS_score\tside\tp.emp\tp.adj_bh\n'
              'G1\t0\t0.9\ttop\t0.1\t0.5\nG1\t1\t0.1\ttop\t0.2\t0.5\n'
              'G1\t0\t0.1\tbottom\t0.1\t0.5\nG1\t1\t0.9\tbottom\t0.2\t0.5\n')
    reg = make_localization_fixture(tmp_path, custom_text='G1\t0\tTERM\tmember\n', scores_text=scores)
    out = run_localization(reg, tmp_path / 'loc', extra=['--scope', 'all'])
    rows = {r['rank_direction']: r for r in read_tsv(out / 'localization_results.tsv')}
    assert {k: float(v['auc']) for k, v in rows.items()} == {'global': .5, 'top': 1., 'bottom': 0.}
    assert all(r['p_mc'] == 'NA' and r['q'] == 'NA' and r['B_valid'] == '0' for r in rows.values())
    manifest = json.loads((out / 'localization_manifest.json').read_text())
    assert manifest['control']['inference'] == 'withheld'


def test_fubar_uses_tested_coverage_not_unannotated_sites():
    import pandas as pd
    members = pd.DataFrame([
        {'collection_id': 'fubar_tested', 'term_id': 'global', 'Gene': 'G', 'Position': '0'},
        {'collection_id': 'fubar_tested', 'term_id': 'global', 'Gene': 'G', 'Position': '1'},
        {'collection_id': 'fubar_positive', 'term_id': 'global', 'Gene': 'G', 'Position': '0'},
    ])
    collections = pd.DataFrame([{'collection_id': 'fubar_positive', 'term_id': 'global',
        'informative_within_protein': '1', 'status': 'included', 'description': 'test'}])
    scores = {('G', 0): .9, ('G', 1): .5, ('G', 2): 1.}
    term = loc.build_gene_table(members, collections, scores, None)[('fubar_positive', 'global')]
    effect, genes, _ = loc.run_term(term, 'mean', 0, loc.term_rng(1, 'caas_score', 'global', 'fubar_positive', 'global'))
    assert effect['auc'] == 1. and genes[0]['n_background'] == 1
    assert term['coverage_status'] == 'fubar_tested_positions'


def test_prototype_all_ties_p_one():
    data = {'genes': {'G': {'positions': dict.fromkeys(range(4), .5), 'member_positions': {0, 1}}}}
    summary, _, _ = loc.run_term(data, 'mean', 99, loc.term_rng(7, 'caas_score', 'top', 'custom', 'flat'))
    assert summary['auc'] == .5 and summary['p_mc'] == 1.


def test_streamed_terms_preserve_aggregation_and_outside_member_audits():
    import pandas as pd
    members = [
        dict(collection_id='custom', term_id='T', Gene='G1', Position='0'),
        dict(collection_id='custom', term_id='T', Gene='G2', Position='0'),
        dict(collection_id='custom', term_id='T', Gene='G2', Position='1'),
        dict(collection_id='custom', term_id='T', Gene='outside', Position='99'),
        dict(collection_id='fubar_positive', term_id='F', Gene='G1', Position='0'),
        dict(collection_id='fubar_positive', term_id='F', Gene='G2', Position='0'),
    ] + [dict(collection_id='fubar_tested', term_id='tested', Gene=g, Position=str(p))
         for g in ('G1', 'G2') for p in (0, 1)]
    collections = pd.DataFrame([
        dict(collection_id=cid, term_id=term, informative_within_protein='1',
             status='included', description='test')
        for cid, term in [('custom', 'T'), ('fubar_positive', 'F')]])
    scores = {('G1', 0): .9, ('G1', 1): .1, ('G1', 2): .5,
              ('G2', 0): .2, ('G2', 1): .3, ('G2', 2): .9, ('G2', 3): .8,
              ('G3', 0): .4, ('G3', 1): .2}
    terms = dict(loc.iter_gene_tables(iter(members), collections, scores, None))
    custom = terms[('custom', 'T')]
    assert custom['n_members_total'] == 4 and len(custom['genes']) == 3
    for aggregation, expected in [('mean', .5), ('pair_weighted', 1 / 3)]:
        result, genes, n = loc.run_term(custom, aggregation, 0, None)
        assert result['auc'] == expected and n == 2
        assert [g['n_background'] for g in genes] == [2, 2]
    fubar = terms[('fubar_positive', 'F')]
    result, genes, n = loc.run_term(fubar, 'mean', 0, None)
    assert result['auc'] == .5 and n == 2
    assert all(g['n_positions_evaluable'] == 2 for g in genes)
    assert fubar['coverage_status'] == 'fubar_tested_positions'


def test_many_terms_do_not_retain_per_term_score_copies():
    import pandas as pd
    import tracemalloc
    # The input is fixed; adding terms must not multiply resident score tables.
    scores = {(f'G{g:02}', p): float(400 - p) for g in range(20) for p in range(400)}
    collections = pd.DataFrame([
        dict(collection_id='custom', term_id=f'T{t:03}', informative_within_protein='1',
             status='included', description='test') for t in range(250)])
    members = [dict(collection_id='custom', term_id=f'T{t:03}', Gene='G00', Position='0')
               for t in range(250)]
    tracemalloc.start()
    try:
        count = 0
        for _, term in loc.iter_gene_tables(iter(members), collections, scores, None):
            result, _, n = loc.run_term(term, 'mean', 0, None)
            assert result['auc'] == 1. and n == 1
            count += 1
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    assert count == 250
    assert peak < 8 * 1024 * 1024, f'term iteration allocated {peak} bytes'


def test_protein_composition_is_not_within_protein_localization(tmp_path):
    """The user's primary question must not inherit a global protein-composition effect."""
    import math
    scores = ['Gene\tPosition\tCAAS_score\tside\tp.emp\tp.adj_bh']
    for gene, value in [('G1', .9), ('G2', .1)]:
        scores.extend(f'{gene}\t{pos}\t{value}\ttop\t.2\t.5' for pos in range(6))
    custom = ''.join(f'G1\t{pos}\tCOMPOSITION\tsites\n' for pos in range(5)) + 'G2\t0\tCOMPOSITION\tsites\n'
    registry = make_localization_fixture(tmp_path, custom_text=custom, scores_text='\n'.join(scores)+'\n')
    local_out = run_localization(registry, tmp_path/'local')
    local = read_tsv(local_out/'localization_results.tsv')[0]
    assert local['n_genes_informative'] == '2' and float(local['auc']) == .5
    global_out = tmp_path/'global'
    result = subprocess.run([sys.executable, str(SRC/'posenrich_global.py'), '--registry_dir', str(registry),
                             '--output_dir', str(global_out)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    global_row = read_tsv(global_out/'global_results.tsv')[0]
    assert math.isclose(float(global_row['auc']), 5/6)
    assert global_row['p'] == local['p_mc'] == 'NA'
