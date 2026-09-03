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


# Tier 0 adapter tests removed — Tier 0 is demoted (validation/.demoted/, D-DIR-01).


def test_map_site_to_alignment() -> None:
    # ref row: 2 leading gaps, residue 1 at col 3, gap, residue 2 at col 5
    assert truthset.map_site_to_alignment(1, "--M-K-R") == 3
    assert truthset.map_site_to_alignment(2, "--M-K-R") == 5
    assert truthset.map_site_to_alignment(3, "--M-K-R") == 7


def _fake_run(tmp: Path) -> Path:
    """A minimal <trait>_complete tree: one recovered site, one missed, one extra."""
    rc = tmp / "c4_complete"
    (rc / "caastools").mkdir(parents=True)
    (rc / "scoring").mkdir()
    (rc / "ct_disambiguation").mkdir()
    # truth site at ref 780 -> pipeline col 779 (0-based), residues A(bg)->S(fg)
    # truth site at ref 665 -> col 664, H->N   (NOT discovered here = missed)
    # an extra CAAS at col 620 not in the truth set
    (rc / "caastools" / "discovery.tab").write_text(
        "gene\tmode\tcaap_group\ttrait\tposition\tcaas\tamino_encoded\tpvalue\tpattern\n"
        "PEPC\tCAAP\tUS\tt\t779\tSSSS/AAAA\tSSSS/AAAA\t0.003\t1\n"
        "PEPC\tCAAP\tGS2\tt\t779\tSSSS/AAAA\txxxx/ssss\t0.003\t1\n"
        "PEPC\tCAAP\tUS\tt\t620\tCCCA/SSSS\tCCCA/SSSS\t0.002\t3\n"
    )
    (rc / "scoring" / "position_scores.tsv").write_text(
        "Gene\tPosition\tpvalue\trecovery_boot\tasr_score\tphen_score\tn_schemes\t"
        "scheme_set\tCAAS_score\tgate_sig\tchange_side\n"
        "PEPC\t779\t0.003\t0\t0.75\t1\t2\tGS2+US\t0.75\tTRUE\ttop\n"
        "PEPC\t620\t0.002\t0\t0.56\t1\t1\tUS\t0.50\tTRUE\ttop\n"
    )
    (rc / "ct_disambiguation" / "caas_convergence_master.csv").write_text(
        "gene,msa_pos,tag,caas,is_significant,convergence_type,caap_group,asr_path_score\n"
        "PEPC,779,CAAS_x,SSSS/AAAA,True,convergent_top,US,0.75\n"
    )
    return rc


def test_tier1_adapter_recovery(tmp_path: Path) -> None:
    from validation.harness.tier1_adapter import score_tier1

    truth = tmp_path / "t.sites.tsv"
    truth.write_text(
        "gene\tref_seq\tposition\tref_aa\talt_aa\ttier\tsource\n"
        "PEPC\tmaize\t780\tA\tS\tmutagenesis\tX\n"
        "PEPC\tmaize\t665\tH\tN\tmutagenesis\tX\n"
    )
    rep = score_tier1(truth, _fake_run(tmp_path), dataset="pepc_c4")

    got = {s.ref_pos: s for s in rep.per_site}
    assert got[780].recovered and got[780].residue_match
    assert got[780].pipeline_pos == 779 and got[780].caas_score == 0.75
    assert got[780].convergence_types == ["convergent_top"]
    assert not got[665].recovered and got[665].caas_score != got[665].caas_score  # nan

    assert rep.recall_by_tier["mutagenesis"]["recovered"] == 1
    assert rep.recall_by_tier["mutagenesis"]["n"] == 2
    lo, hi = rep.recall_by_tier["mutagenesis"]["wilson_95ci"]
    assert 0.0 <= lo <= 0.5 <= hi <= 1.0

    # 780 is the only scored recovered site -> top of 2 -> percentile 1.0
    assert rep.median_rank_percentile["mutagenesis"] == 1.0
    # the col-620 CAAS is not a truth site
    assert [e["ref_pos_est"] for e in rep.extra_calls] == [621]
    # ConDor comparator: 780 & 665 are in its call set; we recovered only 780
    assert rep.comparator["ConDor"]["they_called_we_missed"] == [665]


def test_tier1_adapter_slop_needs_residue_match(tmp_path: Path) -> None:
    from validation.harness.tier1_adapter import score_tier1

    truth = tmp_path / "t.sites.tsv"
    # truth ref 782 -> base col 781; the CAAS is at 779 (within slop 2) but its
    # residues are S/A, and here we ask for V/K -> must NOT be called recovered
    truth.write_text(
        "gene\tref_seq\tposition\tref_aa\talt_aa\ttier\tsource\n"
        "PEPC\tmaize\t782\tK\tV\tselection\tX\n"
    )
    rep = score_tier1(truth, _fake_run(tmp_path), slop=2)
    assert not rep.per_site[0].recovered


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
