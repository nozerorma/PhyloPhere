"""CT and the standalone null list alignments in file-name order and pick the toy sample from that order."""
import json
import os
import random
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_wiring as tw  # noqa: E402

NAMES = [f"GENE{i:02d}.Homo_sapiens.fa" for i in range(30)]


def _java_shuffle(items, seed):
    """Collections.shuffle(items, new java.util.Random(seed)), re-implemented as an independent oracle."""
    mask = (1 << 48) - 1
    state = [(seed ^ 0x5DEECE66D) & mask]

    def next_bits(bits):
        state[0] = (state[0] * 0x5DEECE66D + 0xB) & mask
        v = state[0] >> (48 - bits)
        return v - (1 << bits) if v >= (1 << (bits - 1)) and bits == 32 else v

    def next_int(bound):
        if bound & (bound - 1) == 0:
            return (bound * next_bits(31)) >> 31
        while True:
            bits = next_bits(31)
            val = bits % bound
            if bits - val + (bound - 1) < (1 << 31):
                return val

    out = list(items)
    for i in range(len(out), 1, -1):
        j = next_int(i)
        out[i - 1], out[j] = out[j], out[i - 1]
    return out


def test_the_oracle_matches_the_shuffle_nextflow_runs():
    # Output of `Collections.shuffle(list, new Random(seed))` in the Groovy of Nextflow 25.10.3 (Java), zero-based items
    assert _java_shuffle(list(range(10)), 42) == [4, 6, 2, 1, 7, 9, 8, 5, 3, 0]
    assert _java_shuffle(list(range(10)), 7) == [0, 1, 9, 3, 7, 4, 8, 5, 2, 6]
    assert _java_shuffle(list(range(30)), 1998) == [12, 24, 27, 9, 11, 25, 7, 15, 22, 6, 10, 28, 29, 1, 23, 16, 13, 26, 8, 21, 3, 2, 5, 20, 19, 17, 0, 4, 18, 14]


def _run(tmp_path, orders, n=7, seed=1998):
    ali = tmp_path / "ali"
    ali.mkdir()
    for name in NAMES:
        (ali / name).write_text(">a\nAC\n")
    for junk in ("genes.tsv", "notes.txt", "x.csv", "run.log", "G.map"):
        (ali / junk).write_text("x\n")
    (ali / "subdir").mkdir()
    out = tmp_path / "out.json"
    r = tw._mini(tmp_path, "mini_alignment_files.nf", "--mini_dir", str(ali), "--mini_orders", ";".join(",".join(o) for o in orders),
                 "--mini_n", str(n), "--mini_seed", str(seed), "--out", str(out))
    assert out.exists(), r.stdout[-800:] + r.stderr[-800:]
    return json.loads(out.read_text())


@tw.needs_nextflow
def test_alignments_are_listed_by_name_and_the_sample_does_not_depend_on_the_input_order(tmp_path):
    shuffled = list(NAMES)
    random.Random(5).shuffle(shuffled)
    orders = [NAMES, NAMES[::-1], shuffled]
    res = _run(tmp_path, orders)
    assert res["listed"] == NAMES  # sorted by name; tables, logs and sub-directories are left out
    expected = _java_shuffle(NAMES, 1998)[:7]
    for o in orders:
        assert res[",".join(o)] == expected


@tw.needs_nextflow
def test_another_seed_picks_another_sample_from_the_same_sorted_list(tmp_path):
    res = _run(tmp_path, [NAMES], n=5, seed=7)
    assert res[",".join(NAMES)] == _java_shuffle(NAMES, 7)[:5] != _java_shuffle(NAMES, 1998)[:5]


def test_ct_and_the_standalone_null_use_the_shared_functions_and_no_other_listing_of_alignments():
    """Source-level check: the call sites are inline workflow code that no test can enter."""
    for rel in ("workflows/ct.nf", "main.nf"):
        text = (tw.ROOT / rel).read_text()
        assert "listAlignmentFiles(" in text and "sampleAlignmentFiles(" in text, rel
        assert not re.search(r"alignParam\)\.listFiles|align_dir_f\.listFiles|Collections\.shuffle\((allFiles|all_ali_files)", text), rel
