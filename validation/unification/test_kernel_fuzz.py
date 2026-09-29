"""Seeded differential test: scalar CAAStools discovery vs the vectorized kernel on random inputs.

Each seed draws a random alignment (gaps, ambiguity codes, species absent from the alignment),
a few random hypotheses over a species pool, random gap/missing thresholds (often equal fg/bg,
which is what activates miss_pair), max_conserved, and caap_mode / miss_pair on or off. The scalar
`ct discovery` runs on the hypotheses as traitfile_H<m>.tab; the kernel replays them as the
labelings b_0~H<m>. Rows and the tested-position list must be identical.

Unlike the fixture goldens this checks both directions (a CAAS the kernel misses, one it
invents) and the branches real fixtures do not reach (missing species, asymmetric pairs).
"""
import os
import random
import subprocess
from pathlib import Path

import pandas as pd
import pytest

# CT_BIN lets the mutation check below point the test at a modified copy of `ct`
CT = Path(os.environ.get("CT_BIN", Path(__file__).resolve().parents[2] / "subworkflows/CT/local/ct"))
KEY = ["caap_group", "hyp", "position", "caas", "amino_encoded"]


def _case(rng, tmp):
    n_pairs = rng.choice([2, 3, 4])
    pool = [f"s{i}" for i in range(2 * n_pairs + rng.randint(0, 3))]
    hyps = {}
    for h in range(1, rng.randint(1, 4) + 1):
        sp = rng.sample(pool, 2 * n_pairs)
        hyps[f"H{h}"] = (sp[:n_pairs], sp[n_pairs:])  # fg[k] <-> bg[k] is pair k+1

    absent = set(rng.sample(pool, rng.choice([0, 0, 1, 2])))
    present = [s for s in pool if s not in absent]
    seqs = {s: [] for s in present}
    for _col in range(40):
        r0, r1 = rng.sample("ADKEG", 2)
        for i, s in enumerate(pool):
            if s in absent:
                continue
            x = rng.random()
            aa = "-" if x < 0.10 else "X" if x < 0.13 else rng.choice("ADKEG") if x < 0.30 else (r0 if i % 2 == 0 else r1)
            seqs[s].append(aa)
    (tmp / "G1.fasta").write_text("".join(f">{s}\n{''.join(q)}\n" for s, q in seqs.items()))

    cfg = tmp / "cfg"
    cfg.mkdir()
    for h, (fg, bg) in hyps.items():
        (cfg / f"traitfile_{h}.tab").write_text(
            "".join(f"{f}\t1\t{k + 1}\n{b}\t0\t{k + 1}\n" for k, (f, b) in enumerate(zip(fg, bg))))
    (tmp / "b0.tab").write_text("".join(f"b_0~{h}\t{','.join(fg)}\t{','.join(bg)}\n" for h, (fg, bg) in hyps.items()))

    def pair_caps(name):
        a = rng.choice([0, 1, 2])
        b = a if rng.random() < 0.6 else rng.choice([0, 1, 2])
        return [f"--max_fg_{name}", str(a), f"--max_bg_{name}", str(b), f"--max_{name}",
                str(rng.choice([a + b, a + b + 1, 4]))]

    thr = ["--max_conserved", str(rng.choice([0, 1, 2])), *pair_caps("gaps"), *pair_caps("miss")]
    flags = ["--patterns", "1,2,3"] + (["--miss_pair"] if rng.random() < 0.7 else []) + (["--caap_mode"] if rng.random() < 0.7 else [])
    return thr, flags


def _run(cmd, cwd):
    p = subprocess.run([str(c) for c in cmd], cwd=cwd, capture_output=True, text=True)
    assert p.returncode == 0, p.stdout[-800:] + p.stderr[-800:]


def _rows(path, col):
    if not path.is_file():
        return set()
    d = pd.read_csv(path, sep="\t")
    d["hyp"] = d[col].str.extract(r"(H\d+)")[0]
    return set(map(tuple, d[KEY].astype(str).itertuples(index=False, name=None)))


def _background(path):
    _g, tested = path.read_text().rstrip("\n").split("\t")
    return set() if tested == "NULL" else set(tested.split(","))


@pytest.mark.parametrize("seed", range(60))
def test_kernel_equals_scalar_on_random_input(tmp_path, seed):
    rng = random.Random(seed)
    thr, flags = _case(rng, tmp_path)
    _run([CT, "discovery", "-a", "G1.fasta", "-t", "cfg", "-o", "s.out", "--background_output", "s.bg", "--fmt", "fasta", *flags, *thr], tmp_path)
    _run([CT, "perm-replay", "-a", "G1.fasta", "-t", "cfg", "-s", "b0.tab", "-o", "k.out", "--fmt", "fasta", *flags, *thr,
          "--export_perm_discovery", "k.disc", "--export_b0_background", "k.bg"], tmp_path)
    scalar, kernel = _rows(tmp_path / "s.out", "trait"), _rows(tmp_path / "k.disc", "cycle")
    ctx = f"seed={seed} flags={flags} thr={thr}"
    assert scalar - kernel == set(), f"rows only in the scalar path: {sorted(scalar - kernel)[:3]} {ctx}"
    assert kernel - scalar == set(), f"rows only in the kernel: {sorted(kernel - scalar)[:3]} {ctx}"
    assert _background(tmp_path / "s.bg") == _background(tmp_path / "k.bg"), f"background differs {ctx}"
