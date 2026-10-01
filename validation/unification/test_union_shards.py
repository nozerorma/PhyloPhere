"""union_shards.sh: the union of per-gene shard directories, with hard links or copies."""
import os
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("PHYLOPHERE_ROOT", HERE.parents[1]))
SCRIPT = ROOT / "subworkflows/CT/local/scripts/union_shards.sh"


def _batches(tmp, n_batches, n_genes, b0=True):
    dirs = []
    for b in range(n_batches):
        d = tmp / f"src{b}"
        (d / "b0").mkdir(parents=True)
        for g in range(n_genes):
            (d / f"G{b}_{g}.tsv.gz").write_text(f"main {b} {g}\n")
            if b0:
                (d / "b0" / f"G{b}_{g}.tsv.gz").write_text(f"b0 {b} {g}\n")
        dirs.append(d)
    return dirs


def _union(tmp, sources, env=None, links=False):
    staged = []
    for i, s in enumerate(sources):  # Nextflow stages each input as a symbolic link
        link = tmp / f"batch_{i + 1}"
        link.symlink_to(s)
        staged.append(link.name)
    return subprocess.run(["bash", "-euo", "pipefail", str(SCRIPT), "perm_pos_detail", *staged], cwd=tmp, capture_output=True, text=True,
                          env=dict(os.environ, **(env or {})))


def test_many_shards_do_not_kill_the_probe_with_sigpipe(tmp_path):
    sources = _batches(tmp_path, 50, 400)  # 40000 shards, more than a pipe holds
    r = _union(tmp_path, sources)
    assert r.returncode == 0, r.stderr
    merged = tmp_path / "perm_pos_detail"
    assert sum(1 for _ in merged.glob("*.tsv.gz")) == 50 * 400 and sum(1 for _ in (merged / "b0").glob("*.tsv.gz")) == 50 * 400


def test_the_union_holds_the_sources_own_files_and_merges_their_b0_directories(tmp_path):
    sources = _batches(tmp_path, 3, 4)
    assert _union(tmp_path, sources).returncode == 0
    merged = tmp_path / "perm_pos_detail"
    for s in sources:
        for f in list(s.glob("*.tsv.gz")) + list((s / "b0").glob("*.tsv.gz")):
            target = merged / f.relative_to(s)
            assert target.read_text() == f.read_text() and os.stat(target).st_ino == os.stat(f).st_ino


def test_the_shards_are_copied_when_hard_links_are_refused(tmp_path):
    fake = tmp_path / "bin"
    fake.mkdir()
    (fake / "ln").write_text("#!/usr/bin/env bash\nexit 1\n")
    (fake / "ln").chmod(0o755)
    sources = _batches(tmp_path, 2, 3)
    r = _union(tmp_path, sources, env={"PATH": f"{fake}:{os.environ['PATH']}"})
    assert r.returncode == 0, r.stderr
    merged = tmp_path / "perm_pos_detail"
    for s in sources:
        for f in s.glob("*.tsv.gz"):
            assert (merged / f.name).read_text() == f.read_text() and os.stat(merged / f.name).st_ino != os.stat(f).st_ino


def test_sources_without_shards_give_an_empty_directory_and_leave_no_probe_file(tmp_path):
    sources = _batches(tmp_path, 2, 0, b0=False)
    r = _union(tmp_path, sources)
    assert r.returncode == 0, r.stderr
    assert (tmp_path / "perm_pos_detail").is_dir() and not any((tmp_path / "perm_pos_detail").glob("*.tsv.gz"))
    assert not list(tmp_path.glob("*.link_probe"))
