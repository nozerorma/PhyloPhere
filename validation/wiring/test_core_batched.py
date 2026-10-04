"""CAAS_CORE_BATCHED and CAAS_CORE_MERGE on PEPC with precomputed ASR: the shards do not depend on the batch layout, and the merge equals a copy union followed by the same pass B.

Three genes: PEPC, PEPD (the same alignment and cached ASR under another name) and NOHIT (identical sequences, no
CAAS in any cycle, so `ct perm-replay` exports nothing). Nextflow runs the real process; the reference is the same two
commands run by hand.
"""
import csv
import gzip
import os
import re
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
    # the real labeling (b_0) and three cycles with the labelings of the real one: the plumbing does not depend on the labels
    b0 = (d / "b0.tab").read_text()
    (d / "resample_perms.tab").write_text(b0 + "".join(b0.replace("b_0~", f"b_{k}~") for k in (1, 2, 3)))
    return d


def _nf(tmp_path, inp, size, alignments=None, reuse=None, fop_pairs=None):
    tmp_path.mkdir(parents=True, exist_ok=True)
    i = inp / "observed_inputs"
    extra = []
    if fop_pairs:
        extra += ["--mini_fop_pairs", str(fop_pairs)]
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
    batches = sorted([l for l in (d.parent / ".command.sh").read_text().split("<<'EOF'\n")[1].split("EOF")[0].splitlines()]
                     for d in dirs_two)
    assert batches == [["NOHIT\tNOHIT.fa", "PEPC\tPEPC.fa"], ["PEPD\tPEPD.fa"]]  # two-column manifests, gene-name order
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


# ── CAAS_CORE_MERGE ──────────────────────────────────────────────────────────

_FILTER = ["--gene_filter_mode", "extreme", "--iqr_multiplier", "3.0", "--extreme_threshold", "0.99", "--remove_caas_clusters", "true"]
_TABLES = ["gene_cycle_scores.tsv", "perm_pos_sample.tsv", "perm_pos_quantiles.tsv", "removed_units.tsv"]


def _merge_nf(tmp_path, inp, details, *extra):
    tmp_path.mkdir(parents=True, exist_ok=True)
    r = tw._mini(tmp_path, "mini_core_merge.nf", "--mini_details", ",".join(str(d) for d in details),
                 "--mini_lengths", str(inp / "observed_inputs/gene_ensembl.tsv"), "--outdir", str(tmp_path / "out"),
                 "--seed", "1998", *_FILTER, *extra)
    out = tmp_path / "out/caas_permulation"
    assert (out / "caas_perms.rds").exists(), r.stdout[-1200:] + r.stderr[-1200:]
    listing = tmp_path / "out/pos_detail_dir.txt"
    return out, (Path(listing.read_text().strip()) if listing.exists() else None)


def _reference_merge(tmp_path, inp, details):
    """A plain copy union of the shard directories, then the same pass B, b_0 rebuild and R step."""
    i = inp / "observed_inputs"
    tmp_path.mkdir(parents=True)
    union = tmp_path / "union"
    union.mkdir()
    for d in details:
        subprocess.run(["cp", "-RL", f"{d}/.", str(union)], check=True)
    shutil.copytree(LOCAL, tmp_path / "code", dirs_exist_ok=True)
    rem = ["--gene-lengths", str(i / "gene_ensembl.tsv"), "--gene-filter-mode", "extreme", "--iqr-multiplier", "3.0",
           "--extreme-percentile", "0.99"]
    re_ = str(tmp_path / "code/reaggregate_perm_scores.py")
    out = tmp_path / "out"
    out.mkdir()
    subprocess.run([sys.executable, re_, "--detail", str(union), "--output-dir", str(out), "--seed", "1998", *rem], check=True, capture_output=True)
    (out / "b0").mkdir()
    subprocess.run([sys.executable, re_, "--detail", str(union / "b0"), "--output-dir", str(out / "b0"), "--seed", "1998", *rem], check=True, capture_output=True)
    shutil.copytree(union / "b0", out / "b0/perm_pos_detail")  # the b_0 shards sit next to its scores
    subprocess.run(["Rscript", str(tw.ROOT / "subworkflows/SCORING/local/src/scoring_caas_perms.R"), "--gene-cycle-scores",
                    str(out / "gene_cycle_scores.tsv"), "--output", str(out / "caas_perms.rds")], check=True, capture_output=True)
    return out, union


def _same_rds(a, b):
    r = subprocess.run(["Rscript", "-e", f'q(status = !identical(readRDS("{a}"), readRDS("{b}")))'], capture_output=True)
    return r.returncode == 0


def _assert_same_outputs(got, ref, with_b0=True):
    for name in _TABLES + ["perm_pos_cycle_caas.tsv.gz"]:
        opener = gzip.open if name.endswith(".gz") else open
        assert opener(got / name, "rt").read() == opener(ref / name, "rt").read(), name
    assert _same_rds(got / "caas_perms.rds", ref / "caas_perms.rds")
    if with_b0:
        for name in _TABLES[:3] + ["removed_units.tsv", "perm_pos_cycle_caas.tsv.gz"]:
            opener = gzip.open if name.endswith(".gz") else open
            assert opener(got / "b0" / name, "rt").read() == opener(ref / "b0" / name, "rt").read(), f"b0/{name}"
        assert _shards([got / "b0/perm_pos_detail"]) == _shards([ref / "b0/perm_pos_detail"]) and _shards([ref / "b0/perm_pos_detail"])
    else:
        assert not (got / "b0").exists()


@pytest.fixture(scope="module")
def shard_batches(tmp_path_factory, inp):
    """The shard directories of two one-gene batches (PEPC, PEPD), as CAAS_CORE_BATCHED writes them, with their b_0 shards."""
    dirs = _nf(tmp_path_factory.mktemp("batches"), inp, 1, ["PEPC", "PEPD"])
    assert len(dirs) == 2 and all((d / "b0").is_dir() for d in dirs)
    return dirs


@tw.needs_nextflow
def test_the_merge_equals_a_copy_union_followed_by_the_same_pass_b_and_links_the_shards(tmp_path, inp, shard_batches):
    ref, _ = _reference_merge(tmp_path / "ref", inp, shard_batches)
    out, merged = _merge_nf(tmp_path / "nf", inp, shard_batches)
    _assert_same_outputs(out, ref)
    assert sorted(p.name for p in merged.glob("*.tsv.gz")) == ["PEPC.tsv.gz", "PEPD.tsv.gz"] and (merged / "b0").is_dir()
    for d in shard_batches:  # the union holds the batches' own files, not copies
        for f in d.glob("*.tsv.gz"):
            assert os.stat(merged / f.name).st_ino == os.stat(f).st_ino


@tw.needs_nextflow
def test_the_merge_copies_the_shards_when_hard_links_are_refused(tmp_path, inp, shard_batches):
    fake = tmp_path / "bin"
    fake.mkdir()
    # symbolic links (how Nextflow stages inputs) work; hard links fail, as across filesystems
    (fake / "ln").write_text('#!/usr/bin/env bash\nfor a in "$@"; do case "$a" in -*s*) exec /bin/ln "$@";; esac; done\nexit 1\n')
    (fake / "ln").chmod(0o755)
    (tmp_path / "path.config").write_text(f"env {{ PATH = '{fake}:' + System.getenv('PATH') }}\n")
    ref, _ = _reference_merge(tmp_path / "ref", inp, shard_batches)
    out, merged = _merge_nf(tmp_path / "nf", inp, shard_batches, "-c", str(tmp_path / "path.config"))
    _assert_same_outputs(out, ref)
    for d in shard_batches:
        for f in d.glob("*.tsv.gz"):
            assert os.stat(merged / f.name).st_ino != os.stat(f).st_ino


@tw.needs_nextflow
def test_the_merge_takes_one_shard_directory_as_the_standalone_route_gives_it(tmp_path, inp, shard_batches):
    ref, union = _reference_merge(tmp_path / "ref", inp, shard_batches)
    out, _ = _merge_nf(tmp_path / "nf", inp, [union])
    _assert_same_outputs(out, ref)


@tw.needs_nextflow
def test_the_merge_reads_a_legacy_concatenated_detail_file(tmp_path, inp, shard_batches):
    rows, header = [], None
    for f in sorted((f for d in shard_batches for f in d.glob("*.tsv.gz")), key=lambda f: f.name):  # shard order, as a directory is read
        lines = gzip.open(f, "rt").read().splitlines()
        header = header or lines[0]
        rows += lines[1:]
    legacy = tmp_path / "perm_pos_detail.tsv.gz"
    with gzip.open(legacy, "wt") as fh:
        fh.write("\n".join([header] + rows) + "\n")
    ref, _ = _reference_merge(tmp_path / "ref", inp, shard_batches)
    out, merged = _merge_nf(tmp_path / "nf", inp, [legacy])
    _assert_same_outputs(out, ref, with_b0=False)  # the legacy file holds the null cycles only
    assert merged is None


# ── the observed labeling (b_0) in the same task ─────────────────────────────

@pytest.fixture(scope="module")
def fop_pairs(inp):
    """fop_pairs.tsv with the b_0 PSS weights of the design (the null cycles carry none: equal weights)."""
    f = inp / "fop_pairs.tsv"
    with open(f, "w") as fh:
        fh.write("cycle\thypothesis_id\tpair\tspecies1\tspecies2\tpss_score\n")
        for r in csv.DictReader(open(inp / "cfg/contrast_hypotheses_pairs.tsv"), delimiter="\t"):
            fh.write("\t".join(["b_0", r["hypothesis_id"], r["pair"], r["species1"], r["species2"], r["pss_score"]]) + "\n")
    return f


def _b0_dirs(tmp_path):
    return [Path(l) for l in (tmp_path / "out/b0_observed_dirs.txt").read_text().split()]


@tw.needs_nextflow
def test_the_b0_master_does_not_depend_on_the_batch_size_and_equals_the_frozen_one(tmp_path, inp, fop_pairs):
    import pandas as pd
    genes = ["PEPC", "PEPD", "NOHIT"]
    _nf(tmp_path / "b1", inp, 1, genes, fop_pairs=fop_pairs)
    _nf(tmp_path / "b2", inp, 2, genes, fop_pairs=fop_pairs)
    one, two = _b0_dirs(tmp_path / "b1"), _b0_dirs(tmp_path / "b2")
    assert (len(one), len(two)) == (3, 2)

    def collect(dirs):
        out = {}
        for d in dirs:
            for f in d.glob("*"):
                assert f.name not in out
                out[f.name] = gzip.open(f, "rt").read() if f.suffix == ".gz" else f.read_text()
        return out
    a, b = collect(one), collect(two)
    # NOHIT has no hit under any labeling: it has a tested-positions line, no discovery rows and no master shard
    assert sorted(a) == ["NOHIT.b0.background", "PEPC.b0.background", "PEPC.b0.discovery.tsv", "PEPC.master.csv.gz",
                         "PEPD.b0.background", "PEPD.b0.discovery.tsv", "PEPD.master.csv.gz"]
    assert a == b
    gold = pd.read_csv(GOLD / "caas_convergence_master.csv", keep_default_na=False)
    key = ["msa_pos", "caap_group", "side"]
    gold = gold.sort_values(key, kind="stable").reset_index(drop=True)
    modal = [c for c in gold.columns if re.fullmatch(r"domain_\d+_(anc|top|bot)_aa", c)]
    for g in ("PEPC", "PEPD"):
        shard = next(d for d in one if (d / f"{g}.master.csv.gz").exists()) / f"{g}.master.csv.gz"
        got = pd.read_csv(gzip.open(shard, "rt"), keep_default_na=False).sort_values(key, kind="stable").reset_index(drop=True)
        assert list(got.columns) == list(gold.columns) and len(got) == len(gold) == 217
        for c in gold.columns:
            if c in ("tag_support", "gene") or c in modal:
                continue  # ids are content hashes of the gene; the modal residues are checked below
            if gold[c].dtype.kind == "f":
                assert float((got[c] - gold[c]).abs().max(skipna=True) or 0.0) <= 1e-12, c
            else:
                assert got[c].equals(gold[c]), c
        # The export lists a position's entries in a fixed order, the frozen discovery.tab in the order of the
        # file system: the rows are the same, and a modal residue differs only where two residues tie for the maximum.
        for c in modal:
            for i in got.index[got[c] != gold[c]]:
                counts = dict((r, int(n)) for r, n in (p.split(":") for p in got.loc[i, c + "_support"].split(",")))
                assert counts[got.loc[i, c]] == counts[gold.loc[i, c]] == max(counts.values()), (c, i)


@tw.needs_nextflow
def test_a_batch_that_reuses_exports_has_no_b0_slice(tmp_path, inp, fop_pairs):
    live_dir = tmp_path / "live"
    _nf(live_dir, inp, 2, ["PEPC", "PEPD"], fop_pairs=fop_pairs)
    published = sorted((live_dir / "out/caas_permulation/perm_disc").glob("*.perm_replay.discovery.output"))
    _nf(tmp_path / "reuse", inp, 2, reuse=published, fop_pairs=fop_pairs)
    dirs = _b0_dirs(tmp_path / "reuse")
    assert len(dirs) == 1 and list(dirs[0].iterdir()) == []


# ── CAAS_CORE_OBSERVED: the observed contract files ─────────────────────────

@pytest.fixture(scope="module")
def b0_batches(tmp_path_factory, inp, fop_pairs):
    """b_0 directories of two one-gene batches (PEPC, PEPD), as CAAS_CORE_BATCHED leaves them."""
    root = tmp_path_factory.mktemp("b0batches")
    _nf(root, inp, 1, ["PEPC", "PEPD"], fop_pairs=fop_pairs)
    return _b0_dirs(root)


def _observed(tmp_path, inp, b0_dirs):
    tmp_path.mkdir(parents=True, exist_ok=True)
    tw._mini(tmp_path, "mini_core_observed.nf", "--mini_b0", ",".join(str(d) for d in b0_dirs), "--mini_design", str(inp / "cfg"),
             "--outdir", str(tmp_path / "out"))
    listing = tmp_path / "out/contract_paths.txt"
    if not listing.exists():
        return {}
    got = {}
    for p in map(Path, listing.read_text().split()):
        got[p.name] = p.read_text()
        if p.name == "global_meta_caas.tsv":
            got.update({f"meta_caas/{f.name}": f.read_text() for f in p.parent.glob("*.tsv")})
        if p.name == "caas_convergence_master.csv":
            assert p.parent.name == "ct_disambiguation"
    return got


@tw.needs_nextflow
def test_the_observed_files_come_from_the_b0_slices_in_any_batch_order(tmp_path, inp, b0_batches):
    assert len(b0_batches) == 2
    forward = _observed(tmp_path / "fwd", inp, b0_batches)
    assert forward == _observed(tmp_path / "bwd", inp, b0_batches[::-1])
    # the same files as the command line gives on the same directories
    direct = tmp_path / "direct"
    subprocess.run([sys.executable, str(LOCAL / "contract_main.py"), "--b0-dirs", *map(str, b0_batches), "--design", str(inp / "cfg"),
                    "--output-dir", str(direct)], check=True, capture_output=True)
    expect = {f.name: f.read_text() for f in direct.iterdir() if f.is_file()}
    expect["caas_convergence_master.csv"] = (direct / "ct_disambiguation/caas_convergence_master.csv").read_text()
    expect.update({f"meta_caas/{f.name}": f.read_text() for f in (direct / "meta_caas").glob("*.tsv")})
    assert {k: v for k, v in forward.items() if k in expect} == expect and set(expect) <= set(forward)
    genes = [l.split("\t")[0] for l in forward["discovery.tab"].splitlines()[1:]]
    assert genes == sorted(genes) and set(genes) == {"PEPC", "PEPD"}
    assert forward["background_genes.output"] == "PEPC\nPEPD\n" and forward["background.output"].count("\n") == 2  # one line per gene
    assert len(forward["caas_convergence_master.csv"].splitlines()) == 1 + 2 * 217


@tw.needs_nextflow
def test_the_observed_files_land_in_outdir_where_the_rest_of_the_pipeline_reads_them(tmp_path, inp, b0_batches):
    """publishDir patterns that match nothing publish nothing while the run still succeeds: look in outdir."""
    _observed(tmp_path / "run", inp, b0_batches)
    out = tmp_path / "run/out"
    direct = tmp_path / "direct"
    subprocess.run([sys.executable, str(LOCAL / "contract_main.py"), "--b0-dirs", *map(str, b0_batches), "--design", str(inp / "cfg"),
                    "--output-dir", str(direct)], check=True, capture_output=True)
    for name in ("discovery.tab", "background.output", "background_genes.output"):
        assert (out / "caastools" / name).read_text() == (direct / name).read_text(), name
    assert (out / "ct_disambiguation/caas_convergence_master.csv").read_text() == \
        (direct / "ct_disambiguation/caas_convergence_master.csv").read_text()
    tables = sorted(f.name for f in (direct / "meta_caas").glob("*_meta_caas.tsv"))
    assert len(tables) == 6 and sorted(f.name for f in (out / "meta_caas/meta_caas").iterdir()) == tables
    for t in tables:
        assert (out / "meta_caas/meta_caas" / t).read_text() == (direct / "meta_caas" / t).read_text(), t


@tw.needs_nextflow
def test_batches_without_a_b0_slice_give_no_observed_file(tmp_path, inp, shard_batches):
    empty = tmp_path / "empty"
    empty.mkdir()
    assert _observed(tmp_path / "run", inp, [empty]) == {}


@tw.needs_nextflow
def test_the_core_hands_the_b0_slices_of_its_batches_to_the_observed_step_and_the_null_merge_does_not_wait_for_them(tmp_path, inp, fop_pairs):
    i = inp / "observed_inputs"
    r = tw._mini(tmp_path, "mini_core_chain.nf", "--mini_cfg", str(inp / "cfg"), "--mini_resample", str(inp / "resample_perms.tab"),
                 "--mini_tree", str(i / "pruned_tree_file.nwk"), "--mini_fop_pairs", str(fop_pairs), "--mini_lengths", str(i / "gene_ensembl.tsv"),
                 "--mini_alignments", ",".join(str(inp / "ali" / f"{g}.fa") for g in ("PEPC", "PEPD")), "--outdir", str(tmp_path / "out"),
                 "--alignment", str(inp / "ali"), "--ct_core_batch_size", "1", "--ct_disambig_asr_cache_dir", str(i / "asr_cache"),
                 "--tax_id", str(i / "taxid.tsv"), "--gene_ensembl_file", str(i / "gene_ensembl.tsv"), "--ct_disambig_posterior_threshold", "0.1",
                 "--ct_disambig_asr_model", "lg", "--ali_format", "fasta", "--patterns", "1,2,3", "--miss_pair", "true", "--caap_mode", "true",
                 "--min_divergent_fraction", "0.5", "--seed", "1998")
    listing = tmp_path / "out/chain_paths.txt"
    assert listing.exists(), r.stdout[-1500:] + r.stderr[-1500:]
    files = {Path(p).name: Path(p).read_text() if Path(p).suffix != ".rds" else "" for p in listing.read_text().split()}
    assert sorted(files) == ["background.output", "background_genes.output", "caas_convergence_master.csv", "caas_perms.rds", "discovery.tab",
                             "global_meta_caas.tsv"]
    assert len(files["caas_convergence_master.csv"].splitlines()) == 1 + 2 * 217
    assert {l.split("\t")[0] for l in files["discovery.tab"].splitlines()[1:]} == {"PEPC", "PEPD"}


@tw.needs_nextflow
def test_a_run_without_permuted_labelings_gives_the_observed_files_and_an_empty_null(tmp_path, inp, fop_pairs):
    """N = 0: the labelings file holds b_0 alone, so the null has no shard at all."""
    i = inp / "observed_inputs"
    (inp / "resample_n0.tab").write_text((inp / "b0.tab").read_text())
    r = tw._mini(tmp_path, "mini_core_chain.nf", "--mini_cfg", str(inp / "cfg"), "--mini_resample", str(inp / "resample_n0.tab"),
                 "--mini_tree", str(i / "pruned_tree_file.nwk"), "--mini_fop_pairs", str(fop_pairs), "--mini_lengths", str(i / "gene_ensembl.tsv"),
                 "--mini_alignments", ",".join(str(inp / "ali" / f"{g}.fa") for g in ("PEPC", "PEPD")), "--outdir", str(tmp_path / "out"),
                 "--alignment", str(inp / "ali"), "--ct_core_batch_size", "1", "--ct_disambig_asr_cache_dir", str(i / "asr_cache"),
                 "--tax_id", str(i / "taxid.tsv"), "--gene_ensembl_file", str(i / "gene_ensembl.tsv"), "--ct_disambig_posterior_threshold", "0.1",
                 "--ct_disambig_asr_model", "lg", "--ali_format", "fasta", "--patterns", "1,2,3", "--miss_pair", "true", "--caap_mode", "true",
                 "--min_divergent_fraction", "0.5", "--seed", "1998", "--caas_full_perms", "0")
    listing = tmp_path / "out/chain_paths.txt"
    assert listing.exists(), r.stdout[-2500:] + r.stderr[-2500:]
    paths = {Path(p).name: Path(p) for p in listing.read_text().split()}
    assert {l.split("\t")[0] for l in paths["discovery.tab"].read_text().splitlines()[1:]} == {"PEPC", "PEPD"}
    assert len(paths["caas_convergence_master.csv"].read_text().splitlines()) == 1 + 2 * 217
    assert "caas_perms.rds" in paths


@tw.needs_nextflow
def test_a_discovery_tab_that_exists_is_scored_to_the_frozen_master_and_the_meta_tables(tmp_path, inp):
    i = inp / "observed_inputs"
    disc = tmp_path / "discovery.tab"
    disc.write_text(gzip.open(GOLD / "discovery.tab.gz", "rt").read())
    (tmp_path / "ali").mkdir()
    shutil.copy(GOLD / "PEPC.fasta", tmp_path / "ali/PEPC.fa")
    r = tw._mini(tmp_path, "mini_caas_observed.nf", "--mini_discovery", str(disc), "--mini_design", str(i / "traitfiles"),
                 "--mini_tree", str(i / "pruned_tree_file.nwk"), "--mini_hyp_pairs", str(i / "traitfiles/contrast_hypotheses_pairs.tsv"),
                 "--outdir", str(tmp_path / "out"), "--alignment", str(tmp_path / "ali"), "--ct_disambig_asr_cache_dir", str(i / "asr_cache"),
                 "--tax_id", str(i / "taxid.tsv"), "--gene_ensembl_file", str(i / "gene_ensembl.tsv"), "--ct_disambig_posterior_threshold", "0.1",
                 "--ct_disambig_asr_model", "lg")
    listing = tmp_path / "out/observed_paths.txt"
    assert listing.exists(), r.stdout[-1500:] + r.stderr[-1500:]
    master, meta = sorted(map(Path, listing.read_text().split()), key=lambda p: p.suffix)  # .csv, .tsv
    import pandas as pd
    got = pd.read_csv(master, keep_default_na=False)
    gold = pd.read_csv(GOLD / "caas_convergence_master.csv", keep_default_na=False)
    assert list(got.columns) == list(gold.columns) and len(got) == len(gold) == 217
    for c in gold.columns:  # the ids of tag_support are content hashes; everything else is the frozen master
        if c == "tag_support":
            continue
        if gold[c].dtype.kind == "f":
            assert float((got[c] - gold[c]).abs().max(skipna=True) or 0.0) <= 1e-12, c
        else:
            assert got[c].equals(gold[c]), c
    assert sorted(f.name for f in meta.parent.iterdir() if f.suffix == ".tsv") == [
        "GS1_meta_caas.tsv", "GS2_meta_caas.tsv", "GS3_meta_caas.tsv", "GS4_meta_caas.tsv", "US_meta_caas.tsv", "global_meta_caas.tsv"]
    assert len(meta.read_text().splitlines()) == 1 + len(disc.read_text().splitlines()) - 1
    # published where the rest of the pipeline reads them (a pattern that matches nothing publishes nothing and the run succeeds)
    out = tmp_path / "out"
    assert (out / "ct_disambiguation/caas_convergence_master.csv").read_text() == master.read_text()
    assert (out / "meta_caas/meta_caas/global_meta_caas.tsv").read_text() == meta.read_text()
