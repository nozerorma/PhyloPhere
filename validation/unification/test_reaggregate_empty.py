"""reaggregate_perm_scores.py with a null that has no permuted labeling (N = 0): an explicit empty null, never a silent one."""
import csv
import gzip
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "subworkflows/CT_DISAMBIGUATION/local"
R_SCRIPT = ROOT / "subworkflows/SCORING/local/src/scoring_caas_perms.R"
COLUMNS = ["Gene", "cycle", "Position", "caap_group", "asr_path_score", "n_detected", "clust", "side"]
TABLES = ["gene_cycle_scores.tsv", "perm_pos_cycle_caas.tsv.gz", "perm_pos_quantiles.tsv", "perm_pos_sample.tsv"]


def _shard(directory, gene="G1", cycle="b_1"):
    directory.mkdir(parents=True, exist_ok=True)
    with gzip.open(directory / f"{gene}.tsv.gz", "wt", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(COLUMNS)
        w.writerow([gene, cycle, 10, "US", 0.5, 1, 0, "top"])
        w.writerow([gene, cycle, 11, "US", 0.25, 1, 0, "bottom"])


def _run(detail, out, *extra):
    return subprocess.run([sys.executable, "./reaggregate_perm_scores.py", "--detail", str(detail), "--output-dir", str(out), *extra],
                          cwd=LOCAL, capture_output=True, text=True)


def _header(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        return f.readline().rstrip("\n").split("\t")


def _rows(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        return f.read().splitlines()[1:]


def test_an_empty_null_is_written_with_the_headers_of_a_normal_run(tmp_path):
    _shard(tmp_path / "full")
    assert _run(tmp_path / "full", tmp_path / "ref").returncode == 0
    (tmp_path / "empty/b0").mkdir(parents=True)      # N = 0: only the b_0 subdirectory exists
    r = _run(tmp_path / "empty", tmp_path / "out", "--empty-null")
    assert r.returncode == 0, r.stderr[-800:]
    for name in TABLES:
        assert _header(tmp_path / "out" / name) == _header(tmp_path / "ref" / name), name
        assert _rows(tmp_path / "out" / name) == [], name
    assert not any((tmp_path / "out").glob(".*")) and not list((tmp_path / "out").glob("*empty*"))


def test_a_null_with_no_shard_still_fails_without_the_explicit_flag(tmp_path):
    (tmp_path / "empty").mkdir()
    r = _run(tmp_path / "empty", tmp_path / "out")
    assert r.returncode != 0 and "no *.tsv.gz shards" in r.stderr


def test_the_empty_null_flag_refuses_a_null_that_has_shards(tmp_path):
    _shard(tmp_path / "full")
    r = _run(tmp_path / "full", tmp_path / "out", "--empty-null")
    assert r.returncode != 0 and "unrecognized" not in r.stderr and "--empty-null" in r.stderr and "shard" in r.stderr
    assert not (tmp_path / "out/gene_cycle_scores.tsv").exists()


def test_the_empty_null_flag_does_not_excuse_a_missing_detail_path(tmp_path):
    r = _run(tmp_path / "does_not_exist", tmp_path / "out", "--empty-null")
    assert r.returncode != 0 and "unrecognized" not in r.stderr and "not found" in r.stderr


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript not installed")
def test_the_empty_tables_give_an_empty_caas_perms_rds(tmp_path):
    (tmp_path / "empty").mkdir()
    assert _run(tmp_path / "empty", tmp_path / "out", "--empty-null").returncode == 0
    r = subprocess.run(["Rscript", str(R_SCRIPT), "--gene-cycle-scores", str(tmp_path / "out/gene_cycle_scores.tsv"),
                        "--output", str(tmp_path / "caas_perms.rds")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-500:]
    out = subprocess.run(["Rscript", "-e", f'x <- readRDS("{tmp_path / "caas_perms.rds"}"); cat(length(x$caas_corStat_byrank), x$gene_stat)'],
                         capture_output=True, text=True).stdout
    assert out.strip() == "0 size_adj_max"


def test_the_accumulation_reader_names_the_empty_null_and_the_way_out(tmp_path):
    import ast
    import gzip as _gzip
    from pathlib import Path as _Path
    source = (ROOT / "subworkflows/CT_ACCUMULATION/local/src/randomization/randomize.py").read_text()
    func = next(n for n in ast.parse(source).body if isinstance(n, ast.FunctionDef) and n.name == "_iter_perm_detail_rows")
    scope = {"gzip": _gzip, "Path": _Path}      # the module's heavy imports are not needed by this reader
    exec(compile(ast.Module([func], []), "randomize.py", "exec"), scope)
    (tmp_path / "perm_pos_detail/b0").mkdir(parents=True)
    with pytest.raises(FileNotFoundError, match="--caas_full_perms 0.*naive or cons_decile"):
        list(scope["_iter_perm_detail_rows"](tmp_path / "perm_pos_detail"))
