"""Positional universe downstream of CT: strict universe, null reducer, background.output contract and provenance digest."""
import importlib.util
import random
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tests'))
sys.path.insert(0, str(ROOT / 'subworkflows/CT/local'))
LOCAL = ROOT / 'subworkflows/CT_DISAMBIGUATION/local'
sys.path.insert(0, str(LOCAL))

import test_positional_universe as tu
from src.core import contract, universe as cu
from src.utils.null_tables import CycleReducer
from src.utils.perm_stream import SlimRecord


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_strict_universe_is_the_zero_conserved_universe(tmp_path):
    """--export_b0_universe_strict is the same test with max_conserved 0: a subset of the tolerant universe, equal to the
    oracle evaluated with max_conserved 0."""
    from types import SimpleNamespace
    from modules.perm_replay import run_perm_replay_on_alignment
    from modules.perm_replay_io import labelings_to_multicfg
    seen_strict_smaller = 0
    for seed in range(40):
        rng = random.Random(seed)
        hyps, cols, present, lim, _ = tu._case(rng)
        lim['max_conserved'] = max(1, lim['max_conserved'])
        cfg = labelings_to_multicfg({f'b_0~H{i + 1}': (tuple(f for f, _ in h), tuple(b for _, b in h)) for i, h in enumerate(hyps)})
        sliced = SimpleNamespace(genename='G', species=sorted(present),
                                 d=[{s: f'{c[s]}@{p}' for s in sorted(present)} for p, c in cols])
        a, b = tmp_path / f'a{seed}', tmp_path / f'b{seed}'
        g, m = lim['gaps'], lim['miss']
        run_perm_replay_on_alignment(None, cfg, sliced, g[0], g[1], g[2], m[0], m[1], m[2], lim['patterns'],
                                     miss_pair=lim['miss_pair'], max_conserved=lim['max_conserved'],
                                     min_non_overlapping=lim['min_non'], export_b0_universe=str(a),
                                     export_b0_universe_strict=str(b))
        read = lambda p: set() if (t := p.read_text().rstrip('\n').split('\t')[1]) == 'NULL' else set(t.split(','))
        strict, tolerant = read(b), read(a)
        want = tu._oracle(hyps, cols, present, tu._lim_for_oracle(dict(lim, max_conserved=0)))
        assert strict == want and strict <= tolerant
        seen_strict_smaller += strict < tolerant
    assert seen_strict_smaller > 0


def test_reducer_scores_only_universe_rows_but_judges_all_detections():
    recs = [SlimRecord(p, s, 'top', 1, 1) for p, s in ((1, 0.9), (2, 0.5), (3, 0.4), (50, 0.7))]
    full = CycleReducer(False, False)('b_1', recs)
    cut = CycleReducer(False, False, universe=[2, 3, 99])('b_1', recs)
    assert (full.n_pos, full.max_all, full.n_rows) == (4, 0.9, 4)
    assert (cut.n_pos, cut.max_all, cut.n_rows, cut.positions) == (4, 0.5, 2, [2, 3])
    outside = CycleReducer(False, False, universe=[99])('b_1', recs)
    assert outside is not None and outside.n_pos == 4 and outside.n_rows == 0 and outside.positions == []


def test_background_output_is_the_universe_and_keeps_the_gene_family(tmp_path):
    (tmp_path / 'g1.b0.background').write_text('G1\t1,2,3\n')
    (tmp_path / 'g2.b0.background').write_text('G2\t4,5\n')
    (tmp_path / 'g3.b0.background').write_text('G3\tNULL\n')
    (tmp_path / 'g1.b0.universe').write_text('G1\t2\n')
    (tmp_path / 'g2.b0.universe').write_text('G2\tNULL\n')
    (tmp_path / 'g3.b0.universe').write_text('G3\tNULL\n')
    out, genes = tmp_path / 'background.output', tmp_path / 'background_genes.output'
    contract.write_background(sorted(tmp_path.glob('*.b0.background')), sorted(tmp_path.glob('*.b0.universe')), out, genes)
    assert out.read_text() == 'G1\t2\nG2\tNULL\nG3\tNULL\n'
    assert genes.read_text() == 'G1\n'          # G2 has an empty universe: it is not in the gene family


def test_universe_outside_the_tested_positions_or_genes_is_rejected(tmp_path):
    (tmp_path / 'a.b0.background').write_text('G1\t1,2\n')
    (tmp_path / 'a.b0.universe').write_text('G1\t3\n')
    with pytest.raises(ValueError, match='not tested'):
        contract.write_background([tmp_path / 'a.b0.background'], [tmp_path / 'a.b0.universe'], tmp_path / 'o', tmp_path / 'g')
    (tmp_path / 'a.b0.universe').write_text('G9\t1\n')
    with pytest.raises(ValueError, match='different genes'):
        contract.write_background([tmp_path / 'a.b0.background'], [tmp_path / 'a.b0.universe'], tmp_path / 'o', tmp_path / 'g')


def test_provenance_digest_agrees_between_tiles_and_scoring():
    stats = _load(ROOT / 'subworkflows/SCORING/local/src/null_position_stats.py', 'nps')
    for gene, pos in (('G1', [5, 1, 9]), ('001', []), ('a.b', [0])):
        assert stats.universe_digest(gene, pos) == cu.digest(gene, pos)
    assert cu.digest('G', [1, 2]) != cu.digest('G', [1, 3])
