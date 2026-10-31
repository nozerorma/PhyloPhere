"""Positional universe: the positions b_0 examined that pass the CAAS test under some orientation of the observed pairs.

The oracle is written per position, per hypothesis and per orientation straight from the CT limits and
`caas_id.check_pattern`; the code under test replays the orientations with the vectorized kernel. Thresholds, pair
count, hypotheses, patterns, gaps, missing species and miss_pair vary at random, nothing is fixed to one design."""
import itertools
import random
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
CT = ROOT / 'subworkflows/CT/local'
sys.path.insert(0, str(CT))
from modules.caas_id import check_pattern
from modules.alimport import GAP_SYMBOLS, STANDARD_AAS
from modules.perm_replay import run_perm_replay_on_alignment
from modules.perm_replay_io import labelings_to_multicfg
from modules.universe import orientation_labelings

AMINO = sorted(STANDARD_AAS)[:5]          # few states: conserved columns and convergent columns are both frequent
SYMBOLS = AMINO + ['-', 'X']              # '-' and the ambiguity code are both gaps for the kernel


def _limit(value):
    return None if value == 'NO' else int(value)


def _side_ok(count, limit):
    return limit is None or count <= limit


def _pair_ok(pairs, fg, bg, flagged, gfg, gbg, gall, enabled):
    """miss_pair: with equal side limits, flagged fg and bg members must come from the same pairs (both sides non-empty)."""
    equal = (gfg is not None and gbg is not None and gfg == gbg) or (gfg is None and gbg is None and gall is not None)
    if not (enabled and equal):
        return True
    pf = {k for k, (f, _b) in enumerate(pairs) if f in flagged}
    pb = {k for k, (_f, b) in enumerate(pairs) if b in flagged}
    return not (pf and pb and pf != pb)


def _valid(pairs, column, present, lim):
    """The limits of a labeling at a column: gaps, missing species and miss_pair. pairs = [(fg, bg), ...]."""
    fg = [f for f, _ in pairs]
    bg = [b for _, b in pairs]
    gapped = {s for s in fg + bg if s in present and column[s] in GAP_SYMBOLS}
    missing = {s for s in fg + bg if s not in present}
    for flagged, (a, b, c) in ((gapped, lim['gaps']), (missing, lim['miss'])):
        nf, nb = sum(s in flagged for s in fg), sum(s in flagged for s in bg)
        if not (_side_ok(nf, a) and _side_ok(nb, b) and _side_ok(nf + nb, c)):
            return False
        if not _pair_ok(pairs, fg, bg, flagged, a, b, c, lim['miss_pair']):
            return False
    return True


def _hit(pairs, column, present, lim):
    if not _valid(pairs, column, present, lim):
        return False
    fg = ''.join(column[f] for f, _ in pairs if f in present and column[f] in STANDARD_AAS)
    bg = ''.join(column[b] for _, b in pairs if b in present and column[b] in STANDARD_AAS)
    ok, pattern, *_ = check_pattern(fg, bg, max_conserved=lim['max_conserved'], min_non_overlapping=lim['min_non'])
    return ok and pattern in lim['patterns']


def _examined(pairs, column, present, lim):
    fg = [f for f, _ in pairs]
    bg = [b for _, b in pairs]
    has = lambda side: any(s in present and column[s] != '-' for s in side)
    return has(fg) and has(bg) and _valid(pairs, column, present, lim)


def _oracle(hypotheses, columns, present, lim):
    out = set()
    for pos, column in columns:
        if not any(_examined(h, column, present, lim) for h in hypotheses):
            continue
        for h in hypotheses:
            k = len(h)
            if any(_hit([(b, f) if (mask >> j) & 1 else (f, b) for j, (f, b) in enumerate(h)], column, present, lim)
                   for mask in range(1 << k)):
                out.add(str(pos))
                break
    return out


def _case(rng):
    k = rng.choice([2, 3, 4])
    n_hyp = rng.choice([1, 2, 3])
    species = [f'S{i:02d}' for i in range(rng.randint(2 * k, 2 * k + 6))]
    hypotheses = []
    for _ in range(n_hyp):
        chosen = rng.sample(species, 2 * k)
        hypotheses.append([(chosen[2 * j], chosen[2 * j + 1]) for j in range(k)])
    lim = dict(
        gaps=tuple(rng.choice(['NO', '0', '1', '2']) for _ in range(3)),
        miss=tuple(rng.choice(['NO', '0', '1', '2']) for _ in range(3)),
        miss_pair=rng.random() < 0.6,
        max_conserved=rng.randint(0, k - 1),
        min_non=rng.choice([1, 2, 2, 3]),
        patterns=rng.choice(['1,2,3', '1,2,3,4', '1', '2,3', '1,3']),
    )
    present = {s for s in species if rng.random() < 0.85}
    columns = []
    for pos in range(rng.randint(40, 80)):
        weights = rng.choice([[8, 8, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1], [20, 1, 1, 1, 1, 3, 1]])
        columns.append((pos, {s: rng.choices(SYMBOLS, weights)[0] for s in present}))
    return hypotheses, columns, present, lim, species


def _lim_for_oracle(lim):
    """The limits as (fg, bg, all) tuples of ints/None, the order the kernel takes them in."""
    conv = lambda t: tuple(_limit(x) for x in t)
    return dict(lim, gaps=conv(lim['gaps']), miss=conv(lim['miss']),
                patterns=set(lim['patterns'].split(',')))


@pytest.mark.parametrize('seed', range(60))
def test_universe_matches_the_per_orientation_oracle(seed, tmp_path):
    rng = random.Random(seed)
    hypotheses, columns, present, lim, species = _case(rng)
    got = _run_simple(tmp_path, hypotheses, columns, present, lim)
    want = _oracle(hypotheses, columns, present, _lim_for_oracle(lim))
    assert got == want, (seed, sorted(got ^ want)[:10])


def test_orientations_keep_the_pairs_and_cover_every_swap():
    out = orientation_labelings({'b_0~H1': (('A', 'C', 'E'), ('B', 'D', 'F'))})
    assert len(out) == 8
    for fg, bg in out.values():
        assert {frozenset(p) for p in zip(fg, bg)} == {frozenset('AB'), frozenset('CD'), frozenset('EF')}
    assert out['b_0~H1#0'] == (('A', 'C', 'E'), ('B', 'D', 'F'))
    assert out['b_0~H1#7'] == (('B', 'D', 'F'), ('A', 'C', 'E'))
    assert len({v for v in out.values()}) == 8


def test_pair_bound_is_explicit():
    with pytest.raises(ValueError, match='exceed the bound'):
        orientation_labelings({'b_0': (tuple('ABCDEFGHIJKLM'), tuple('abcdefghijklm'))})


def test_observed_detections_are_inside_the_universe(tmp_path):
    """Mask 0 is the observed orientation, so a column the observed labeling detects is always eligible."""
    rng = random.Random(7)
    for _ in range(25):
        hypotheses, columns, present, lim, _ = _case(rng)
        lim = dict(lim, gaps=('NO',) * 3, miss=('NO',) * 3, miss_pair=False)
        L = _lim_for_oracle(lim)
        got = _run_simple(tmp_path, hypotheses, columns, present, lim)
        for pos, column in columns:
            if any(_hit(h, column, present, L) for h in hypotheses):
                assert str(pos) in got


def _run_simple(tmp_path, hypotheses, columns, present, lim):
    labelings = {f'b_0~H{i + 1}': (tuple(f for f, _ in h), tuple(b for _, b in h)) for i, h in enumerate(hypotheses)}
    sliced = SimpleNamespace(genename='G', species=sorted(present),
                             d=[{s: f'{col[s]}@{pos}' for s in sorted(present)} for pos, col in columns])
    out = tmp_path / 'u.tsv'
    g, m = lim['gaps'], lim['miss']
    run_perm_replay_on_alignment(None, labelings_to_multicfg(labelings), sliced, g[0], g[1], g[2], m[0], m[1], m[2],
                                 lim['patterns'], miss_pair=lim['miss_pair'], max_conserved=lim['max_conserved'],
                                 min_non_overlapping=lim['min_non'], export_b0_universe=str(out))
    text = out.read_text().rstrip('\n').split('\t')[1]
    return set() if text == 'NULL' else set(text.split(','))


def test_dunn_index_ignores_which_pair_member_is_foreground(tmp_path):
    """The universe re-orients pairs without re-selecting them: the original Dunn functions of the contrast selector give
    the same value for all 2**K orientations of a hypothesis (they read each pair as an unordered set of species)."""
    import shutil
    rscript = shutil.which('Rscript')
    if rscript is None:
        pytest.skip('Rscript is not available')
    script = tmp_path / 'dunn.R'
    script.write_text(f'''
suppressPackageStartupMessages(library(ape))
source("{CT}/scripts/lean_contrast_selector.R")
set.seed(11)
worst <- 0; n <- 0
for (rep in 1:30) {{
  K <- sample(2:5, 1)
  tree <- rtree(2 * K + sample(0:4, 1)); D <- cophenetic.phylo(tree)
  sp <- sample(tree$tip.label, 2 * K)
  pairs <- lapply(seq_len(K), function(j) c(sp[2 * j - 1], sp[2 * j]))
  base <- overall_dunn_lean(D, pairs)
  for (mask in 0:(2^K - 1)) {{
    swapped <- lapply(seq_len(K), function(j) if (bitwAnd(bitwShiftR(mask, j - 1), 1L) == 1L) rev(pairs[[j]]) else pairs[[j]])
    worst <- max(worst, abs(overall_dunn_lean(D, swapped) - base)); n <- n + 1
  }}
}}
cat(n, worst, "\\n")
''')
    run = subprocess.run([rscript, str(script)], capture_output=True, text=True)
    assert run.returncode == 0, run.stderr[-800:]
    n, worst = run.stdout.split()[-2:]
    assert int(n) > 100 and float(worst) == 0.0


# ── Downstream: the null scores only the universe, the contract writes it, SCORING recognises it ──────────────────────

LOCAL = ROOT / 'subworkflows/CT_DISAMBIGUATION/local'
sys.path.insert(0, str(LOCAL))
sys.path.insert(0, str(ROOT / 'subworkflows/SCORING/local/src'))


def _records(rows):
    from src.utils.perm_stream import SlimRecord
    return [SlimRecord(p, s, side, 1, 1) for p, s, side in rows]


def test_reducer_takes_the_gene_maximum_over_the_universe_only():
    from src.utils.null_tables import CycleReducer
    rows = [(10, 0.9, 'top'), (20, 0.4, 'top'), (30, 0.7, 'bottom'), (40, 0.0, 'top')]
    free = CycleReducer(False, False)('b_1', _records(rows))
    inside = CycleReducer(False, False, universe=[20, 40])('b_1', _records(rows))
    assert free.max_all == 0.9 and free.n_all == 3
    assert inside.max_all == 0.4 and inside.n_all == 1 and inside.n_top == 1 and inside.n_bottom == 0
    assert inside.positions == [20] and inside.n_rows == 2          # scored 20 and the unscored 40
    assert inside.n_pos == free.n_pos == 4                           # detections still judge gene removal and trains


def test_a_cycle_detected_only_outside_the_universe_is_a_silent_unit_not_a_missing_one():
    from src.utils.null_tables import CycleReducer
    unit = CycleReducer(False, False, universe=[99])('b_1', _records([(10, 0.9, 'top')]))
    assert unit is not None and unit.n_all == 0 and unit.n_rows == 0 and unit.n_pos == 1
    assert unit.max_all != unit.max_all                              # NaN: no score in the universe


def test_background_is_the_universe_and_the_gene_family_is_the_genes_with_universe_positions(tmp_path):
    from src.core.contract import write_background
    b = tmp_path / 'a.b0.background'; u = tmp_path / 'a.b0.universe'
    b.write_text('G1\t1,2,3,4\nG2\t5,6\nG3\tNULL\n')
    u.write_text('G1\t2,4\nG2\tNULL\nG3\tNULL\n')
    write_background([b], [u], tmp_path / 'background.output', tmp_path / 'genes.output')
    assert (tmp_path / 'background.output').read_text() == 'G1\t2,4\nG2\tNULL\nG3\tNULL\n'
    # G2 was tested but has no eligible position and G3 was never tested: neither is in the gene family
    assert (tmp_path / 'genes.output').read_text() == 'G1\n'


@pytest.mark.parametrize('tested,universe,message', [
    ('G1\t1,2\n', 'G1\t3\n', 'not tested'),
    ('G1\t1,2\n', 'G1\t1\nG2\t1\n', 'different genes'),
    ('G1\t1,2\n', 'G1\t1\nG1\t2\n', 'Duplicate gene'),
])
def test_background_rejects_inconsistent_universe(tmp_path, tested, universe, message):
    from src.core.contract import write_background
    (tmp_path / 'a.b0.background').write_text(tested); (tmp_path / 'a.b0.universe').write_text(universe)
    with pytest.raises(ValueError, match=message):
        write_background([tmp_path / 'a.b0.background'], [tmp_path / 'a.b0.universe'],
                         tmp_path / 'b.output', tmp_path / 'g.output')


def test_universe_identity_agrees_between_tile_provenance_and_scoring():
    from src.core.universe import digest, parse_line
    import null_position_stats as scoring
    gene, positions = parse_line('G.1\t7,3,12\n')
    assert positions == (3, 7, 12)
    assert digest(gene, positions) == scoring.universe_digest(gene, [12, 3, 7])
    assert digest('G', ()) == scoring.universe_digest('G', [])
    assert digest('G', (1,)) != digest('G', (2,))


# ── SCORING rejects a null scored over another universe (validate_null) ───────────────────────────────────────────────

def _null_with_universe(tmp_path, tested_text, recorded):
    """A minimal canonical null block whose provenance records `recorded` ({gene: positions} as the tiles saw them, digested
    with the tile-side function; None = a null made before the universe) and the call validate_null makes against
    `tested_text` (the content of background.output)."""
    import json
    import pyarrow.parquet as pq
    sys.path.insert(0, str(ROOT / 'tests'))
    from posenrich_fixtures import canonical, seal_block
    from src.core.universe import digest
    import null_position_stats as scoring
    observed = tmp_path / 'observed.tsv'; observed.write_text('Gene\tPosition\nG1\t1\n')
    block = canonical(tmp_path / 'null', observed, ['G1', 'G2'], rows=[('G1', 'b_1', 1, .5, .2), ('G2', 'b_2', 4, .3, .6)])
    counts = json.loads((block / 'counts.json').read_text())
    prov = counts['scoring_provenance']
    if recorded is None:
        prov.pop('universe', None)
    else:
        prov['universe'] = {g: digest(g, p) for g, p in recorded.items()}
    (block / 'counts.json').write_text(json.dumps(counts)); seal_block(block)
    background = tmp_path / 'background.output'; background.write_text(tested_text)
    root = tmp_path / 'null'
    contract = json.loads(Path(str(observed) + '.scoring.json').read_text())
    return lambda: scoring.validate_null(
        [block], (root / 'cycle_roster.txt').read_text().split(), pq.read_table(root / 'null_positions.parquet'),
        json.loads((root / 'counts.json').read_text()), None, contract, scoring.read_tested(background))


def test_validate_null_accepts_the_universe_the_tiles_used(tmp_path):
    _null_with_universe(tmp_path, 'G1\t1,2\nG2\t4\n', {'G1': (1, 2), 'G2': (4,)})()


@pytest.mark.parametrize('tested,recorded,message', [
    ('G1\t1,2\nG2\t4\n', None, 'without positional-universe provenance'),
    ('G1\t1,2,3\nG2\t4\n', {'G1': (1, 2), 'G2': (4,)}, 'differs from the tested positions'),
    ('G1\t1,2\n', {'G1': (1, 2), 'G2': (4,)}, 'differs from the tested positions'),
])
def test_validate_null_rejects_a_null_scored_over_another_universe(tmp_path, tested, recorded, message):
    with pytest.raises(ValueError, match=message):
        _null_with_universe(tmp_path, tested, recorded)()
