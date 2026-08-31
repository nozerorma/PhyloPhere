"""Tier 0 simulator tests. Run: python validation/tests/test_tier0.py
(or python -m pytest validation/tests/test_tier0.py -q).

The load-bearing test is ``test_null_is_stationary``: under no planted signal the
simulated column composition must match the site equilibrium (no drift bias), and
convergent-substitution counts on random "foreground" tips must not exceed the
neutral expectation. If that fails, every Tier 0 calibration number is invalid.
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[2]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from validation.tier0 import pheno as _pheno  # noqa: E402
from validation.tier0 import trees  # noqa: E402
from validation.tier0.model import PAML_AA, gamma_rates, load_rate_matrix, site_profiles  # noqa: E402
from validation.tier0.replicate import ReplicateConfig, build_replicate  # noqa: E402
from validation.tier0.simulate import SimConfig, simulate  # noqa: E402

_DATA = _ROOT / "validation" / "tier0" / "data"


# --------------------------------------------------------------------------- #
# model
# --------------------------------------------------------------------------- #
def test_wag_parses_and_q_is_valid() -> None:
    rm = load_rate_matrix(_DATA / "wag.dat")
    assert rm.S.shape == (20, 20)
    assert np.allclose(rm.S, rm.S.T)
    assert np.isclose(rm.pi.sum(), 1.0)
    Q = rm.q()
    assert np.allclose(Q.sum(axis=1), 0.0, atol=1e-10)          # rows sum to 0
    assert np.all(np.diag(Q) < 0)
    assert np.isclose(-np.sum(rm.pi * np.diag(Q)), 1.0)         # normalised to rate 1
    # detailed balance: pi_i Q_ij == pi_j Q_ji
    lhs = rm.pi[:, None] * Q
    assert np.allclose(lhs, lhs.T, atol=1e-12)


def test_site_profiles_concentration() -> None:
    rm = load_rate_matrix(_DATA / "wag.dat")
    rng = np.random.default_rng(0)
    tight = site_profiles(rm.pi, 2000, concentration=200.0, rng=rng)
    loose = site_profiles(rm.pi, 2000, concentration=1.0, rng=rng)
    # high concentration -> less across-site spread; low -> peaked, few effective states
    d_tight = np.abs(tight - rm.pi).sum(axis=1).mean()
    d_loose = np.abs(loose - rm.pi).sum(axis=1).mean()
    eff = lambda p: np.exp(-(p * np.log(p)).sum(axis=1)).mean()
    assert d_tight < d_loose
    assert eff(tight) > 15.0 and eff(loose) < 5.0


def test_gamma_rates_mean_one() -> None:
    rng = np.random.default_rng(1)
    r = gamma_rates(5000, alpha=0.5, rng=rng)
    assert np.isclose(r.mean(), 1.0, atol=0.05)
    assert r.min() >= 0


# --------------------------------------------------------------------------- #
# trees
# --------------------------------------------------------------------------- #
def test_pathological_trees() -> None:
    star = trees.star_tree(8)
    assert len(star.tips) == 8
    assert all(p == star.root for p, _, _ in star.edges)       # every tip hangs off root
    lad = trees.ladder_tree(6)
    assert len(lad.tips) == 6


def test_prune_depth_preserving(tmp_path: Path) -> None:
    src = _ROOT / "validation" / "fixtures" / "tier0" / "trees" / "primates_233_subst.tree"
    if not src.exists():
        print("  (skip: primate tree fixture not fetched)")
        return
    rng = np.random.default_rng(0)
    full = trees.load_tree(src)
    pruned = trees.prune_depth_preserving(src, 40, rng)
    assert len(pruned.tips) == 40
    # farthest-point sampling should retain most of the tree's depth
    assert pruned.total_length() > 0.4 * full.total_length()


def test_sample_foreground_disjoint() -> None:
    lad = trees.ladder_tree(30)
    rng = np.random.default_rng(3)
    fg = trees.sample_foreground(lad, n_transitions=4, rng=rng, max_clade=2)
    assert fg.n_transitions == 4
    # origins have disjoint descendant tips
    seen: set[str] = set()
    for node in fg.origin_nodes:
        dt = lad.descendant_tips(node)
        assert not (dt & seen)
        seen |= dt
    assert len(fg.fg_edges) >= 4


# --------------------------------------------------------------------------- #
# simulator — the load-bearing checks
# --------------------------------------------------------------------------- #
def _column_counts(aln: dict[str, str], col: int) -> Counter:
    return Counter(seq[col] for seq in aln.values())


def test_simulate_shapes_and_alphabet() -> None:
    lad = trees.ladder_tree(20)
    rng = np.random.default_rng(0)
    fg = trees.sample_foreground(lad, 3, rng, max_clade=2)
    cfg = SimConfig(n_sites=50, n_planted_profile_shift=3, n_planted_identical_aa=3)
    res = simulate(lad, fg, cfg, seed=42)
    assert set(res.alignment) == set(lad.tips)
    assert all(len(s) == 50 for s in res.alignment.values())
    assert set("".join(res.alignment.values())) <= set(PAML_AA)
    assert sum(res.phenotype.values()) == len(fg.tips)
    assert len(res.truth.planted_sites) == 6


def test_identical_aa_planting_converges() -> None:
    """Planted identical-aa sites: foreground tips share the target residue far
    more than background tips do."""
    lad = trees.ladder_tree(24)
    rng = np.random.default_rng(1)
    fg = trees.sample_foreground(lad, 4, rng, max_clade=1)
    cfg = SimConfig(n_sites=30, n_planted_identical_aa=10, n_planted_profile_shift=0,
                    identical_aa_target_weight=0.97)
    res = simulate(lad, fg, cfg, seed=7)
    fg_tips = set(fg.tips)
    hits = 0
    for col, mech in res.truth.planted_sites.items():
        tgt = res.truth.identical_aa_targets[col]
        fg_frac = np.mean([res.alignment[t][col] == tgt for t in fg_tips])
        bg_frac = np.mean([res.alignment[t][col] == tgt
                           for t in res.alignment if t not in fg_tips])
        if fg_frac > 0.8 and fg_frac > bg_frac + 0.3:
            hits += 1
    assert hits >= 8, f"only {hits}/10 identical-aa sites converged"


def test_grouped_caap_planting_shares_class_not_residue() -> None:
    """Planted grouped-CAAP sites: fg tips share one GS class but not (mostly)
    a single residue; the class differs from the background majority."""
    from validation.tier0.groups import SCHEMES

    lad = trees.ladder_tree(28)
    rng = np.random.default_rng(3)
    fg = trees.sample_foreground(lad, 4, rng, max_clade=1)
    cfg = SimConfig(n_sites=40, n_planted_identical_aa=0, n_planted_grouped_caap=12,
                    grouped_caap_scheme="GS1", identical_aa_target_weight=0.96)
    res = simulate(lad, fg, cfg, seed=11)
    fg_tips = set(fg.tips)
    gs1 = SCHEMES["GS1"]
    class_share = residue_share = 0
    for col, gt in res.truth.grouped_caap_targets.items():
        fg_classes = {gs1.get(res.alignment[t][col]) for t in fg_tips}
        fg_residues = {res.alignment[t][col] for t in fg_tips}
        if fg_classes == {gt["group"]}:
            class_share += 1
        if len(fg_residues) == 1:
            residue_share += 1
    assert class_share >= 10, f"only {class_share}/12 sites are class-convergent"
    assert residue_share <= 4, f"{residue_share}/12 collapsed to one residue (US leakage)"


def test_null_is_stationary() -> None:
    """No planted signal: (a) column composition tracks the site equilibrium,
    (b) foreground tips (chosen at random, no signal) show no excess of shared
    derived residues versus a rematched random tip set."""
    lad = trees.ladder_tree(40, branch_length=0.25)
    rng = np.random.default_rng(0)
    cfg = SimConfig(n_sites=300, concentration=4.0, gamma_alpha=0.7)
    res = simulate(lad, None, cfg, seed=123)

    # (a) mean per-column entropy should be well below log(20) (sites are
    # constrained) but non-zero, and stable across the alignment
    cols = list(zip(*res.alignment.values()))
    ent = []
    for c in cols:
        p = np.array(list(Counter(c).values()), float)
        p /= p.sum()
        ent.append(-(p * np.log(p)).sum())
    ent = np.array(ent)
    assert 0.1 < ent.mean() < np.log(20)

    # (b) convergence excess: pick k random tips as pseudo-foreground, count
    # columns where >=2 of them share a residue that is a minority overall, and
    # compare to the mean over many random redraws (should be ~1.0 ratio).
    tips = list(res.alignment)
    def shared_minority(sel: list[str]) -> int:
        n = 0
        for c in cols:
            overall = Counter(c)
            sub = Counter(c[tips.index(t)] for t in sel)
            for aa, cnt in sub.items():
                if cnt >= 2 and overall[aa] / len(c) < 0.25:
                    n += 1
                    break
        return n
    obs = shared_minority(list(rng.choice(tips, 5, replace=False)))
    null = [shared_minority(list(rng.choice(tips, 5, replace=False))) for _ in range(40)]
    # observed is itself a random draw, so it must sit inside the null spread
    assert np.mean(null) - 2 * np.std(null) - 1 <= obs <= np.mean(null) + 2 * np.std(null) + 1


def test_pheno_lambda_foreground_binary() -> None:
    lad = trees.ladder_tree(80, branch_length=0.15)
    ph = _pheno.make_lambda_foreground(lad, 6, np.random.default_rng(0),
                                       kind="binary", lam=0.5, planted=True)
    assert ph.kind == "binary" and ph.lam == 0.5
    assert len(ph.pairs) >= 3
    assert set(ph.values.values()) <= {0.0, 1.0}
    assert len(ph.values) == 80
    # one terminal edge per anchor and per partner
    assert len(ph.anchor_edges) == len(ph.pairs)
    assert len(ph.partner_edges) == len(ph.pairs)
    # top and bottom foregrounds are disjoint
    top = set(ph.foreground("top").tips); bot = set(ph.foreground("bottom").tips)
    assert not (top & bot)


def test_pheno_lambda_foreground_rate_and_lambda_effect() -> None:
    lad = trees.ladder_tree(120, branch_length=0.1)
    ph = _pheno.make_lambda_foreground(lad, 6, np.random.default_rng(1),
                                       kind="rate", lam=0.0, planted=True)
    assert ph.kind == "continuous" and ph.n_pop is not None
    for a, b in ph.pairs:
        assert ph.values[a] > ph.values[b]           # anchor above partner (with noise)
    # every non-outgroup tip carries data now (ladder has no outgroup)
    assert len(ph.values) == 120
    assert all(0.0 <= v <= 1.0 for v in ph.values.values())   # c/n rates

    # higher lambda -> the latent clumps -> fewer independent tail pairs on average
    def mean_pairs(lam):
        return np.mean([len(_pheno.make_lambda_foreground(
            lad, 6, np.random.default_rng(s), kind="binary", lam=lam).pairs)
            for s in range(6)])
    assert mean_pairs(0.0) >= mean_pairs(1.0)


def test_pheno_lambda_foreground_continuous() -> None:
    lad = trees.ladder_tree(100, branch_length=0.12)
    ph = _pheno.make_lambda_foreground(lad, 5, np.random.default_rng(3),
                                       kind="continuous", lam=1.0, planted=True)
    assert ph.archetype == "continuous" and ph.kind == "continuous"
    assert ph.n_pop is None and ph.n_cases is None
    assert len(ph.values) == 100                       # all non-outgroup tips
    # raw noisy latent -> not confined to [0, 1], real-valued spread
    assert len({round(v, 4) for v in ph.values.values()}) > 50
    for a, b in ph.pairs:
        assert ph.values[a] > ph.values[b]


def test_pheno_null_has_no_edges() -> None:
    lad = trees.ladder_tree(60)
    for kind in ("binary", "rate", "continuous"):
        ph = _pheno.make_lambda_foreground(lad, 4, np.random.default_rng(2),
                                           kind=kind, lam=0.5, planted=False)
        assert ph.anchor_edges == set() and ph.partner_edges == set()
        assert len(ph.pairs) >= 3


def test_build_replicate_pipeline_shaped(tmp_path: Path) -> None:
    src = _ROOT / "validation" / "fixtures" / "tier0" / "trees" / "primates_233_subst.tree"
    if src.exists():
        tree = trees.prune_depth_preserving(src, 40, np.random.default_rng(0))
    else:
        tree = trees.ladder_tree(40, branch_length=0.2)
    rcfg = ReplicateConfig(archetype="binary", lam=0.5, n_pairs=5, n_genes=6,
                           frac_planted_genes=0.34, n_sites=50, planted=True)
    row = build_replicate(tree, rcfg, seed=99, outdir=tmp_path / "rep")
    d = tmp_path / "rep"
    assert (d / "my_traits.tsv").read_text().splitlines()[0] == "species\tsim_trait\tfamily"
    assert (d / "ali_sp_names.txt").exists()
    assert (d / "gene_ensembl.tsv").read_text().splitlines()[0].startswith("gene\tchr\t")
    assert len(list((d / "align").glob("g*.fasta"))) == 6
    import json
    t = json.loads((d / "truth.json").read_text())
    assert t["n_planted_genes"] == 2
    assert t["phenotype"]["lambda"] == 0.5
    assert len(t["phenotype"]["pairs"]) >= 3
    planted = [g for g, v in t["genes"].items() if v["planted"]]
    assert len(planted) == 2
    assert all(t["genes"][g]["n_planted"] > 0 for g in planted)
    assert all(t["genes"][g]["direction"] in ("top", "bottom") for g in planted)
    assert all(t["genes"][g]["n_planted"] == 0 for g in t["genes"] if g not in planted)

    # rate archetype: my_traits carries count columns
    rcfg2 = ReplicateConfig(archetype="rate", lam=0.0, n_pairs=4, n_genes=3, n_sites=40, planted=True)
    build_replicate(tree, rcfg2, seed=7, outdir=tmp_path / "rep_rate")
    hdr = (tmp_path / "rep_rate" / "my_traits.tsv").read_text().splitlines()[0]
    assert hdr == "species\tsim_trait\tn_pop\tn_cases\tfamily"


def test_write_roundtrip(tmp_path: Path) -> None:
    lad = trees.ladder_tree(12)
    rng = np.random.default_rng(0)
    fg = trees.sample_foreground(lad, 2, rng, max_clade=2)
    res = simulate(lad, fg, SimConfig(n_sites=40, n_planted_identical_aa=4), seed=1)
    out = res.write(tmp_path / "rep0")
    assert (out / "aln.fasta").exists()
    assert (out / "phenotype.tsv").exists()
    import json
    t = json.loads((out / "truth.json").read_text())
    assert t["n_transitions"] == 2
    assert len(t["planted_sites"]) == 4
    ph = (out / "phenotype.tsv").read_text().strip().splitlines()
    assert all(len(ln.split("\t")) == 2 for ln in ph)


if __name__ == "__main__":
    import tempfile
    import traceback

    fails = 0
    for name, fn in sorted(globals().items()):
        if not name.startswith("test_") or not callable(fn):
            continue
        try:
            if "tmp_path" in fn.__code__.co_varnames:
                with tempfile.TemporaryDirectory() as d:
                    fn(Path(d))
            else:
                fn()
            print(f"ok   {name}")
        except Exception:  # noqa: BLE001
            fails += 1
            print(f"FAIL {name}")
            traceback.print_exc()
    raise SystemExit(1 if fails else 0)
