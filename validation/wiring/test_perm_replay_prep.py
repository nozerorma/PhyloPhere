"""CAAS_PERMS_PREP on PEPC: the batch manifest has two columns and the exported perm-discovery files are published."""
import csv
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_wiring as tw  # noqa: E402

GOLD = tw.ROOT / "validation/unification/golden/pepc_c4_complete"


def _inputs(tmp_path):
    cfg = tmp_path / "cfg"
    cfg.mkdir()
    (cfg / "contrast_hypotheses_pairs.tsv").write_text((GOLD / "contrast_hypotheses_pairs.tsv").read_text())
    hyps = {}
    for r in csv.DictReader(open(cfg / "contrast_hypotheses_pairs.tsv"), delimiter="\t"):
        hyps.setdefault(r["hypothesis_id"], []).append(r)
    for h, rows in hyps.items():
        (cfg / f"traitfile_{h}.tab").write_text(
            "".join(f"{r['species1']}\t1\t{r['pair']}\n{r['species2']}\t0\t{r['pair']}\n" for r in rows))
    p = subprocess.run(["python3", str(tw.ROOT / "subworkflows/CT/local/scripts/build_b0_labelings.py"), "--config", str(cfg),
                        "--fop", "--labelings-out", "b0.tab", "--pairs-out", "/dev/null"], cwd=tmp_path, capture_output=True, text=True)
    assert p.returncode == 0, p.stdout + p.stderr
    # a permuted cycle with the labelings of the real one: the plumbing does not depend on the labels
    res = tmp_path / "resample"
    res.mkdir()
    (res / "fop_labelings.tab").write_text((tmp_path / "b0.tab").read_text().replace("b_0~", "b_1~"))
    ali = tmp_path / "ali"
    ali.mkdir()
    for g in ("PEPC", "PEPD"):
        (ali / f"{g}.fa").write_text((GOLD / "PEPC.fasta").read_text())
    return cfg, res, ",".join(str(ali / f"{g}.fa") for g in ("PEPC", "PEPD"))


@tw.needs_nextflow
def test_the_batch_manifest_has_two_columns_and_the_discovery_exports_are_published(tmp_path):
    cfg, res, alis = _inputs(tmp_path)
    r = tw._mini(tmp_path, "mini_perms_prep.nf", "--mini_alignments", alis, "--mini_cfg", str(cfg), "--mini_resample", str(res),
                 "--outdir", str(tmp_path / "out"), "--ct_perm_replay_batch_size", "2", "--caas_full_perms", "1",
                 "--ali_format", "fasta", "--patterns", "1,2,3", "--miss_pair", "true", "--caap_mode", "true",
                 "--min_divergent_fraction", "0.5")
    pub = sorted(p.name for p in (tmp_path / "out/caas_permulation/perm_disc").glob("*.perm_replay.discovery.output"))
    assert pub == ["PEPC.perm_replay.discovery.output", "PEPD.perm_replay.discovery.output"], r.stdout[-800:] + r.stderr[-800:]
    scripts = list((tmp_path / "work").glob("*/*/.command.sh"))
    batch = [s for s in scripts if "PERM_REPLAY_BATCHED" in s.read_text() or "perm_replay_batch" in s.read_text()]
    assert len(batch) == 1
    manifest = batch[0].read_text().split("<<'EOF'\n")[1].split("EOF")[0].splitlines()
    assert [len(l.split("\t")) for l in manifest] == [2, 2], manifest
