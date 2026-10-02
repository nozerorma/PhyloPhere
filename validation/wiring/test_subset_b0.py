"""SUBSET_RESAMPLE_PERMS: the real labeling b_0 is always in the subset, and N = 0 keeps b_0 only.

Nextflow runs the real process on a small synthetic design, in plain mode (a trait file, resample_*.tab) and in FOP
mode (a directory of hypotheses, fop_labelings.tab and fop_pairs.tsv). The null cycles are the first N in cycle order.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_wiring as tw  # noqa: E402

PAIRS_HEADER = "hypothesis_id\tpair\tspecies1\tspecies2\tpss_score\n"
FOP_HEADER = "cycle\thypothesis_id\tpair\tspecies1\tspecies2\tpss_score\n"


def _plain(d):
    cfg = d / "cfg.tab"
    cfg.write_text("a\t1\t1\nb\t0\t1\nc\t1\t2\nd\t0\t2\n")
    res = d / "resamples"
    res.mkdir()
    (res / "resample_1.tab").write_text("".join(f"b_{k}\ta,c\tb,d\n" for k in (2, 1)) + "b_0\ta,c\tb,d\n")
    (res / "resample_2.tab").write_text("b_3\ta,d\tb,c\n")
    return cfg, res


def _fop(d):
    cfg = d / "cfg"
    cfg.mkdir()
    pairs = [("H1", "1", "a", "b", "0.5"), ("H1", "2", "c", "d", "0.25"), ("H2", "1", "a", "d", "0.75"), ("H2", "2", "c", "b", "0.125")]
    (cfg / "contrast_hypotheses_pairs.tsv").write_text(PAIRS_HEADER + "".join("\t".join(p) + "\n" for p in pairs))
    for h in ("H1", "H2"):
        rows = [p for p in pairs if p[0] == h]
        (cfg / f"traitfile_{h}.tab").write_text("".join(f"{p[2]}\t1\t{p[1]}\n{p[3]}\t0\t{p[1]}\n" for p in rows))
    res = d / "resamples"
    res.mkdir()
    lab = "".join(f"b_{k}~{h}\ta,c\tb,d\n" for k in (3, 1, 2) for h in ("H1", "H2"))
    (res / "fop_labelings.tab").write_text(lab + "b_0~H1\ta,c\tb,d\nb_0~H2\ta,c\tb,d\n")
    (res / "fop_pairs.tsv").write_text(FOP_HEADER + "".join(f"b_{k}\t" + "\t".join(p) + "\n" for k in (1, 2, 3) for p in pairs))
    return cfg, res


def _subset(tmp_path, cfg, res, n):
    r = tw._mini(tmp_path, "mini_subset.nf", "--mini_cfg", str(cfg), "--mini_resample", str(res), "--outdir", str(tmp_path / "out"),
                 "--caas_full_perms", str(n))
    listing = tmp_path / "out/subset_paths.txt"
    assert listing.exists(), r.stdout[-1500:] + r.stderr[-1500:]
    paths = [Path(p) for p in listing.read_text().split()]
    subset = next(p for p in paths if p.name == "resample_perms.tab")
    pairs = next((p for p in paths if p.name == "fop_pairs.tsv"), None)
    return [l.split("\t")[0] for l in subset.read_text().splitlines()], pairs


@tw.needs_nextflow
@pytest.mark.parametrize("n,expected", [(0, ["b_0"]), (2, ["b_0", "b_1", "b_2"])])
def test_plain_subset_holds_b0_and_the_first_n_cycles(tmp_path, n, expected):
    cfg, res = _plain(tmp_path)
    tags, pairs = _subset(tmp_path, cfg, res, n)
    assert tags == expected and pairs is None


@tw.needs_nextflow
@pytest.mark.parametrize("n,cycles", [(0, []), (2, ["b_1", "b_2"])])
def test_fop_subset_holds_b0_and_the_first_n_cycles_with_their_pairs(tmp_path, n, cycles):
    cfg, res = _fop(tmp_path)
    tags, pairs = _subset(tmp_path, cfg, res, n)
    assert sorted(tags) == sorted(f"{c}~{h}" for c in ["b_0"] + cycles for h in ("H1", "H2"))
    rows = [l.split("\t") for l in pairs.read_text().splitlines()[1:]]
    assert sorted({r[0] for r in rows}) == ["b_0"] + cycles
    b0 = {(r[1], r[2]): r[5] for r in rows if r[0] == "b_0"}
    assert b0 == {("H1", "1"): "0.5", ("H1", "2"): "0.25", ("H2", "1"): "0.75", ("H2", "2"): "0.125"}  # PSS strings verbatim
