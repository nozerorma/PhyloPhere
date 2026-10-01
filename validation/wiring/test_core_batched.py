"""CAAS_CORE_BATCHED on PEPC with precomputed ASR: the shards do not depend on the batch layout.

Three genes: PEPC, PEPD (the same alignment and cached ASR under another name) and NOHIT (identical sequences, no
CAAS in any cycle, so `ct perm-replay` exports nothing). Nextflow runs the real process; the reference is the same two
commands run by hand.
"""
import csv
import gzip
import os
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_wiring as tw  # noqa: E402

GOLD = tw.ROOT / "validation/unification/golden/pepc_c4_complete"
CT = tw.ROOT / "subworkflows/CT/local/ct"
LOCAL = tw.ROOT / "subworkflows/CT_DISAMBIGUATION/local"
ARGS = ["--patterns", "1,2,3", "--miss_pair", "--caap_mode", "--max_conserved", "2", "--max_fg_gaps", "0", "--max_bg_gaps", "0",
        "--max_gaps", "0", "--max_fg_miss", "0", "--max_bg_miss", "0", "--max_miss", "0"]


@pytest.fixture(scope="module")
def inp(tmp_path_factory):
    d = tmp_path_factory.mktemp("core")
    with tarfile.open(GOLD / "observed_inputs.tar.gz") as t:
        t.extractall(d)
    i = d / "observed_inputs"
    for g in ("PEPD",):
        shutil.copytree(i / "asr_cache/asr_PEPC", i / "asr_cache" / f"asr_{g}")
    (i / "gene_ensembl.tsv").write_text("gene\tchr\tstart\tend\tstrand\tlength\thuman_protein_id\n" + "".join(
        f"{g}\tchr1\t100000\t102910\t+\t970\tP04711\n" for g in ("PEPC", "PEPD", "NOHIT")))
    ali = d / "ali"
    ali.mkdir()
    fa = (GOLD / "PEPC.fasta").read_text()
    for g in ("PEPC", "PEPD"):
        (ali / f"{g}.fa").write_text(fa)
    # identical sequences: no position can be convergent
    recs = [r for r in fa.split(">") if r.strip()]
    first_seq = "".join(recs[0].split("\n")[1:])
    (ali / "NOHIT.fa").write_text("".join(f">{r.split(chr(10))[0]}\n{first_seq}\n" for r in recs))
    hyps = {}
    for r in csv.DictReader(open(i / "traitfiles/contrast_hypotheses_pairs.tsv"), delimiter="\t"):
        hyps.setdefault(r["hypothesis_id"], []).append(r)
    cfg = d / "cfg"
    cfg.mkdir()
    for h, rows in hyps.items():
        (cfg / f"traitfile_{h}.tab").write_text(
            "".join(f"{r['species1']}\t1\t{r['pair']}\n{r['species2']}\t0\t{r['pair']}\n" for r in rows))
    (cfg / "contrast_hypotheses_pairs.tsv").write_text((i / "traitfiles/contrast_hypotheses_pairs.tsv").read_text())
    p = subprocess.run(["python3", str(tw.ROOT / "subworkflows/CT/local/scripts/build_b0_labelings.py"), "--config", str(cfg),
                        "--fop", "--labelings-out", "b0.tab", "--pairs-out", "/dev/null"], cwd=d, capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    (d / "resample_perms.tab").write_text((d / "b0.tab").read_text().replace("b_0~", "b_1~"))
    return d


def _nf(tmp_path, inp, size, alignments=None, reuse=None):
    tmp_path.mkdir(parents=True, exist_ok=True)
    i = inp / "observed_inputs"
    extra = []
    if alignments:
        extra += ["--mini_alignments", ",".join(str(inp / "ali" / f"{g}.fa") for g in alignments)]
    if reuse:
        extra += ["--mini_reuse", ",".join(str(p) for p in reuse)]
    r = tw._mini(tmp_path, "mini_core.nf", "--mini_cfg", str(inp / "cfg"), "--mini_resample", str(inp / "resample_perms.tab"),
                 "--mini_tree", str(i / "pruned_tree_file.nwk"), "--outdir", str(tmp_path / "out"), "--alignment", str(inp / "ali"),
                 "--ct_core_batch_size", str(size), "--ct_disambig_asr_cache_dir", str(i / "asr_cache"),
                 "--tax_id", str(i / "taxid.tsv"), "--gene_ensembl_file", str(i / "gene_ensembl.tsv"),
                 "--ct_disambig_posterior_threshold", "0.1", "--ct_disambig_asr_model", "lg", "--ali_format", "fasta",
                 "--patterns", "1,2,3", "--miss_pair", "true", "--caap_mode", "true", "--min_divergent_fraction", "0.5", *extra)
    listing = tmp_path / "out/pos_detail_dirs.txt"
    assert listing.exists(), r.stdout[-1200:] + r.stderr[-1200:]
    return [Path(l) for l in listing.read_text().split()]


def _shards(dirs):
    out = {}
    for d in dirs:
        for f in Path(d).glob("*.tsv.gz"):
            assert f.name not in out
            out[f.name] = gzip.open(f, "rt").read()
    return out


def _direct(tmp_path, inp, genes):
    """The same replay and detail-only ASR replay, run by hand."""
    i = inp / "observed_inputs"
    (tmp_path / "alignments").mkdir(parents=True)
    for g in genes:
        shutil.copy(inp / "ali" / f"{g}.fa", tmp_path / "alignments")
    (tmp_path / "m.tsv").write_text("".join(f"{g}\t{g}.fa\n" for g in genes))
    (tmp_path / "args.txt").write_text(" ".join(ARGS))
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    subprocess.run(["bash", str(tw.ROOT / "subworkflows/CT/local/scripts/run_ct_perm_replay_batch.sh"), "--batch-id", "b", "--manifest", "m.tsv",
                    "--caas-config", str(inp / "cfg"), "--resampled-path", str(inp / "resample_perms.tab"), "--workers", "2",
                    "--ali-format", "fasta", "--ct-bin", str(CT), "--export-perm-discovery", "1", "--extra-args-file", "args.txt"],
                   cwd=tmp_path, env=env, check=True, capture_output=True)
    (tmp_path / "perm_disc").mkdir()
    for f in tmp_path.glob("*.perm_replay.discovery.output"):
        f.rename(tmp_path / "perm_disc" / f.name)
    shutil.copy(inp / "resample_perms.tab", tmp_path)
    shutil.copytree(LOCAL, tmp_path / "code", dirs_exist_ok=True)
    (tmp_path / "out").mkdir()
    subprocess.run([sys.executable, str(tmp_path / "code/disambiguation_perms_main.py"), "--alignment-dir", str(inp / "ali"),
                    "--tree", str(i / "pruned_tree_file.nwk"), "--perm-discovery", "perm_disc", "--resample-dir", ".",
                    "--output-dir", "out", "--detail-only", "--asr-model", "lg", "--posterior-threshold", "0.1", "--workers", "2",
                    "--max-tasks-per-child", "50", "--asr-cache-dir", str(i / "asr_cache"), "--seed", "1998",
                    "--taxid-mapping", str(i / "taxid.tsv"), "--ensembl-genes-file", str(i / "gene_ensembl.tsv")],
                   cwd=tmp_path, check=True, capture_output=True)
    return _shards([tmp_path / "out/perm_pos_detail"])


@tw.needs_nextflow
def test_the_shards_do_not_depend_on_the_batch_size_and_equal_a_direct_run(tmp_path, inp):
    genes = ["PEPC", "PEPD", "NOHIT"]
    dirs_one = _nf(tmp_path / "b1", inp, 1, genes)
    dirs_two = _nf(tmp_path / "b2", inp, 2, genes)
    assert (len(dirs_one), len(dirs_two)) == (3, 2)  # batches of 1 and of 2 genes, in gene-name order
    for d in dirs_one + dirs_two:  # pass B scores genome-wide pools; a batch only writes its shards
        assert (d.parent / "caas_perms_out/perm_pos_detail").is_dir()
        assert not (d.parent / "caas_perms_out/gene_cycle_scores.tsv").exists()
    one, two = _shards(dirs_one), _shards(dirs_two)
    ref = _direct(tmp_path / "ref", inp, genes)
    assert set(ref) == {"PEPC.tsv.gz", "PEPD.tsv.gz"} and len(ref["PEPC.tsv.gz"]) > 1000  # NOHIT contributes nothing
    assert one == ref and two == ref


@tw.needs_nextflow
def test_a_batch_whose_genes_export_nothing_gives_an_empty_shard_directory(tmp_path, inp):
    dirs = _nf(tmp_path, inp, 1, ["NOHIT"])
    assert len(dirs) == 1 and _shards(dirs) == {}


@tw.needs_nextflow
def test_reused_exports_give_the_same_shards_as_the_replay(tmp_path, inp):
    live_dir = tmp_path / "live"
    live = _shards(_nf(live_dir, inp, 2, ["PEPC", "PEPD"]))
    published = sorted((live_dir / "out/caas_permulation/perm_disc").glob("*.perm_replay.discovery.output"))
    assert [p.name for p in published] == ["PEPC.perm_replay.discovery.output", "PEPD.perm_replay.discovery.output"]
    reused_dir = tmp_path / "reuse"
    reused = _shards(_nf(reused_dir, inp, 2, reuse=published))
    assert reused == live and live
    assert list((reused_dir / "out/caas_permulation/perm_disc").glob("*")) == []  # reused exports are not published again
