"""Smoke tests for the validation harness. Run: python -m pytest validation/tests -q
(or plain `python validation/tests/test_harness.py`)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

_ROOT = Path(__file__).resolve().parents[2]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from validation.harness.phenotype import DatasetSpec, PhenotypeSpec  # noqa: E402
from validation.harness.trait_matrix import emit_trait_matrix, percentilize  # noqa: E402
from validation.harness.metrics import null_calibration, score_metrics  # noqa: E402
from validation.harness import truthset  # noqa: E402


def _toy_dataset(tmp: Path) -> DatasetSpec:
    tips = [f"sp{i}" for i in range(12)]
    tree = "(" + ",".join(f"{t}:0.1" for t in tips) + ");"
    tp = tmp / "tree.nwk"
    tp.write_text(tree + "\n")
    cont = PhenotypeSpec(
        name="metric", kind="continuous",
        values={t: float(i) for i, t in enumerate(tips)},
    )
    cont2 = PhenotypeSpec(
        name="metric2", kind="continuous",
        values={t: float((i * 7) % 12) for i, t in enumerate(tips)},
    )
    cat = PhenotypeSpec(
        name="binary", kind="categorical",
        labels={**{t: "fg" for t in tips[:4]}, **{t: "bg" for t in tips[4:8]}},
        pairs=[(tips[0], tips[4]), (tips[1], tips[5]), (tips[2], tips[6]), (tips[3], tips[7])],
    )
    return DatasetSpec(name="toy", tree=tp, traits=[cont, cont2, cat], provenance="toy")


def test_emit_trait_matrix(tmp_path: Path) -> None:
    ds = _toy_dataset(tmp_path)
    out = tmp_path / "matrix"
    manifest = emit_trait_matrix(ds, out, quantile=0.25)

    # categorical direct: single cfg for the one cat trait
    cfg = (out / "categorical/direct/single/binary.cfg").read_text().strip().splitlines()
    assert len(cfg) == 8  # 4 pairs x 2 rows
    assert all(len(ln.split("\t")) == 3 for ln in cfg)
    assert {ln.split("\t")[1] for ln in cfg} == {"0", "1"}

    # categorical percentilized: one cfg per continuous trait
    assert (out / "categorical/percentilized/single/metric.q25.cfg").exists()
    assert (out / "categorical/percentilized/single/metric2.q25.cfg").exists()

    # continuous direct: headerless single + annotated multi
    single = (out / "continuous/direct/single/metric.traitvalues.tsv").read_text().strip().splitlines()
    assert all(len(ln.split("\t")) == 2 for ln in single)
    multi_hdr = (out / "continuous/direct/multi/traitvalues.tsv").read_text().splitlines()[0]
    assert multi_hdr.split("\t") == ["species", "metric", "metric2"]

    # continuous percentilized: value + pctile column
    pc_hdr = (out / "continuous/percentilized/single/metric.traitvalues.tsv").read_text().splitlines()[0]
    assert pc_hdr.split("\t") == ["species", "metric", "metric__pctile"]

    # multi cfg dir exists for n_trait>1 continuous -> percentilized
    assert (out / "categorical/percentilized/multi").is_dir()

    m = json.loads((out / "manifest.json").read_text())
    assert m["spec_sha256"]
    assert set(m["emitted"]) == {
        "categorical/direct", "categorical/percentilized",
        "continuous/direct", "continuous/percentilized",
    }


def test_percentilize_tails() -> None:
    spec = PhenotypeSpec(name="x", kind="continuous",
                         values={f"s{i}": float(i) for i in range(20)})
    fg, bg, mid = percentilize(spec, 0.25)
    assert len(fg) == 5 and len(bg) == 5 and len(mid) == 10
    assert "s19" in fg and "s0" in bg


def test_score_metrics_recovers_top() -> None:
    scores = {f"g{i}": float(i) for i in range(100)}
    truth = {"g99", "g98", "g97", "g5"}  # 3 easy, 1 buried
    rep = score_metrics(scores, truth, threshold=95.0, higher_is_hit=True)
    assert rep.tp == 3 and rep.fn == 1
    assert rep.positive_ranks["g99"] == 1
    assert rep.positive_ranks["g5"] == 95
    assert 0.0 <= rep.roc_auc <= 1.0 and rep.pr_auc > 0.5


def test_null_calibration_flags_uniform_vs_skewed() -> None:
    rng = np.random.default_rng(0)
    good = null_calibration(rng.uniform(size=2000))
    assert good.verdict("0.05") == "OK"
    bad = null_calibration(rng.uniform(size=2000) ** 3)  # inflated small p-values
    assert bad.verdict("0.05") != "OK"


def _fake_run(root: Path, *, planted: bool) -> Path:
    """A minimal Tier 0 replicate tree the adapter can read."""
    rep = root / ("power_x" if planted else "null_x") / "rep000"
    res = rep / "out" / "scoring"
    perm = rep / "out" / "caas_permulation"
    ct = rep / "out" / "data_exploration" / "2.CT" / "1.Traitfiles"
    for d in (res, perm, ct):
        d.mkdir(parents=True, exist_ok=True)

    genes = {}
    for i in range(1, 6):
        g = f"g{i:04d}"
        is_p = planted and i <= 2
        genes[g] = {
            "planted": is_p,
            "planted_sites": ({"3": "identical_aa", "7": "profile_shift"} if is_p else {}),
            "n_planted": 2 if is_p else 0,
        }
    (rep / "truth.json").write_text(json.dumps({
        "planted": planted, "archetype": "echo", "genes": genes,
        "phenotype": {"foreground_tips": ["spA", "spB", "spC"]},
    }))

    gh = "Gene\tn_positions\tgene_caas_score\n"
    grows = ("g0001\t2\t0.9\ng0002\t2\t0.7\n" if planted else "")
    grows += "g0003\t2\t0.10\ng0004\t2\t0.05\ng0005\t2\t0.02\n"
    (res / "gene_scores.tsv").write_text(gh + grows)

    ph = "Gene\tPosition\tpvalue\tpvalue_boot\tCAAS_score\n"
    prows = ("g0001\t3\t0.001\t0.01\t0.8\ng0002\t3\t0.002\t0.02\t0.7\n" if planted else "")
    prows += "g0003\t44\t0.2\t0.6\t0.1\n"
    (res / "position_scores.tsv").write_text(ph + prows)

    perm_h = "Gene\tPosition\tcaap_group\tn_detected\tn_cycles\tnull_pvalue_boot\n"
    perm_rows = "".join(
        f"g000{3 + j}\t{10 + j}\tUS\t1\t50\t{v}\n"
        for j, v in enumerate([0.3, 0.6])
    )
    (perm / "perm_pos_pval.tsv").write_text(perm_h + perm_rows)

    (ct / "traitfile.tab").write_text("spA\t1\t1\nspX\t0\t1\nspB\t1\t2\nspY\t0\t2\n")
    return root


def test_tier0_adapter_score(tmp_path: Path) -> None:
    from validation.harness import tier0_adapter as t0

    _fake_run(tmp_path, planted=True)
    res = t0.score(tmp_path)
    assert res.n_power_replicates == 1
    rs = res.per_replicate[0]
    assert rs.gene_planted_ranks == {"g0001": 1, "g0002": 2}
    assert rs.gene_precision_at_k == 1.0
    assert rs.site_recall_by_mechanism["identical_aa"] == 1.0
    assert rs.site_recall_by_mechanism["profile_shift"] == 0.0
    assert rs.contrast_pairs_recovered == 2
    assert rs.contrast_fg_precision == 1.0


def test_tier0_adapter_calibrate(tmp_path: Path) -> None:
    from validation.harness import tier0_adapter as t0

    _fake_run(tmp_path, planted=False)
    res = t0.calibrate(tmp_path, "perm_pos_boot")
    assert res.n_pvalues == 2
    assert res.report.n == 2


def test_map_site_to_alignment() -> None:
    # ref row: 2 leading gaps, residue 1 at col 3, gap, residue 2 at col 5
    assert truthset.map_site_to_alignment(1, "--M-K-R") == 3
    assert truthset.map_site_to_alignment(2, "--M-K-R") == 5
    assert truthset.map_site_to_alignment(3, "--M-K-R") == 7


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
