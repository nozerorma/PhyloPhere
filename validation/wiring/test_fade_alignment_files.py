"""The FADE selection lists alignments in file-name order and draws its toy sample from that order, as CT does."""
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_wiring as tw  # noqa: E402
from test_alignment_files import NAMES, _java_shuffle  # noqa: E402

ROOT = tw.ROOT
GENES = [n.split(".")[0] for n in NAMES]
EXTRA_ELIGIBLE = ["PLAIN", "UPPER.FA", "alt.PHY", "other.fasta", "x.aln", "y.phylip"]          # FADE takes these too
EXTRA_IGNORED = ["notes.txt", "table.tsv", "tree.nwk", "run.log", "a.map", "data.json", "b.fas"]  # and these are not alignments for it


def _run(tmp_path, names, n=7, seed=1998):
    ali = tmp_path / "ali"
    ali.mkdir(parents=True)
    for name in names:
        (ali / name).write_text(">a\nAC\n")
    (ali / "subdir").mkdir()
    out = tmp_path / "out.json"
    r = tw._mini(tmp_path, "mini_fade_files.nf", "--mini_dir", str(ali), "--mini_n", str(n), "--mini_seed", str(seed), "--out", str(out))
    assert out.exists(), r.stdout[-1200:] + r.stderr[-1200:]
    return json.loads(out.read_text())


@tw.needs_nextflow
def test_the_tuples_come_in_file_name_order_whatever_the_creation_order(tmp_path):
    shuffled = list(NAMES)
    random.Random(11).shuffle(shuffled)
    res = _run(tmp_path, shuffled)
    assert res["listed"] == NAMES


@tw.needs_nextflow
def test_the_toy_genes_are_the_seeded_shuffle_of_the_sorted_list_and_equal_the_ct_sample(tmp_path):
    shuffled = list(NAMES)
    random.Random(5).shuffle(shuffled)
    res = _run(tmp_path, shuffled, n=7, seed=1998)
    expected = [n.split(".")[0] for n in _java_shuffle(NAMES, 1998)[:7]]
    assert res["toy"] == expected and res["ct_toy"] == expected
    other = _run(tmp_path / "seed7", NAMES, n=5, seed=7)
    assert other["toy"] == [n.split(".")[0] for n in _java_shuffle(NAMES, 7)[:5]] != expected[:5]


@tw.needs_nextflow
def test_the_eligible_files_are_still_the_ones_of_the_extension_filter(tmp_path):
    res = _run(tmp_path, NAMES + EXTRA_ELIGIBLE + EXTRA_IGNORED)
    assert sorted(res["listed"]) == sorted(NAMES + EXTRA_ELIGIBLE)
    assert res["listed"] == sorted(res["listed"])


@tw.needs_nextflow
def test_wanted_genes_keep_only_their_files_in_name_order(tmp_path):
    assert _run(tmp_path, NAMES)["wanted"] == ["GENE01.Homo_sapiens.fa", "GENE03.Homo_sapiens.fa"]


def test_the_selection_uses_the_shared_listing_and_sample_and_no_longer_shuffles_itself():
    text = (ROOT / "subworkflows/SELECTION/selection_prep.nf").read_text()
    assert "sampleAlignmentFiles(" in text and "listAlignmentFiles(" in text
    assert "ct_alignment_files" in text
    assert "Collections.shuffle" not in text and "new Random" not in text and "listFiles()" not in text
