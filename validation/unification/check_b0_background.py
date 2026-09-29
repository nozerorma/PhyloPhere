#!/usr/bin/env python3
"""Per-gene check of the kernel's b_0 background against the observed background.output.

For every gene of `<run>/caastools/background.output`, replays the observed hypotheses as the
labelings `b_0~H<m>` through `ct perm-replay --export_b0_background` (same thresholds as the
run) and compares the tested-position list with the observed line. The discovery rows are
compared too, against `<run>/caastools/discovery.tab`. Exit status 1 on any difference.

Thresholds are the run's fractions, resolved as the pipeline does: int(n_pairs * fraction).
Run it on a compute node for large fixtures (one `ct perm-replay` call per gene). It can also be
streamed to a cluster without copying it: `ssh host srun ... python3 - --repo <checkout> ... < this_file`:

  srun -p high-cpu -c 4 python validation/unification/check_b0_background.py \\
       --run <results_dir> --align-dir <alignments> --tmp-dir <scratch>/.tmp
"""
import argparse
import csv
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd



def n_pairs_of(traitfile):
    with open(traitfile) as fh:
        return len({r[2] for r in csv.reader(fh, delimiter="\t") if len(r) >= 3 and r[2].strip().isdigit()})


def find_alignment(align_dir, gene):
    hits = sorted(p for p in Path(align_dir).iterdir() if p.name.split(".", 1)[0] == gene)
    return hits[0] if hits else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True, help="pipeline results dir (caastools/, data_exploration/2.CT/1.Traitfiles/)")
    ap.add_argument("--align-dir", required=True)
    ap.add_argument("--fmt", default="fasta")
    ap.add_argument("--patterns", default="1,2,3")
    ap.add_argument("--min-divergent-fraction", type=float, default=0.5)
    for k in ("bg_gaps", "fg_gaps", "gaps", "bg_miss", "fg_miss", "miss"):
        ap.add_argument(f"--max-{k.replace('_', '-')}-fraction", type=float, default=0.0)
    ap.add_argument("--no-miss-pair", action="store_true")
    ap.add_argument("--repo", help="PhyloPhere checkout providing `ct` (default: the one containing this script)")
    ap.add_argument("--tmp-dir", help="parent for the work dir (never /tmp on a cluster)")
    ap.add_argument("--limit", type=int, help="first N genes only (debugging)")
    a = ap.parse_args()

    root = Path(a.repo) if a.repo else Path(__file__).resolve().parents[2]
    ct = root / "subworkflows/CT/local/ct"
    b0_script = root / "subworkflows/CT/local/scripts/build_b0_labelings.py"
    run = Path(a.run)
    cfg = run / "data_exploration/2.CT/1.Traitfiles"
    n_pairs = n_pairs_of(cfg / "traitfile_H1.tab")
    cap = lambda f: int(n_pairs * f)
    thr = ["--max_conserved", str(int(n_pairs * (1 - a.min_divergent_fraction))),
           "--max_bg_gaps", str(cap(a.max_bg_gaps_fraction)), "--max_fg_gaps", str(cap(a.max_fg_gaps_fraction)),
           "--max_gaps", str(cap(a.max_gaps_fraction)), "--max_bg_miss", str(cap(a.max_bg_miss_fraction)),
           "--max_fg_miss", str(cap(a.max_fg_miss_fraction)), "--max_miss", str(cap(a.max_miss_fraction))]
    flags = ["--caap_mode"] + ([] if a.no_miss_pair else ["--miss_pair"])

    obs_bg = {}
    for line in open(run / "caastools/background.output"):
        g, _, pos = line.rstrip("\n").partition("\t")
        if g:
            obs_bg[g] = set() if pos in ("", "NULL") else set(pos.split(","))
    genes = sorted(obs_bg)[: a.limit]
    disc = pd.read_csv(run / "caastools/discovery.tab", sep="\t", usecols=["gene", "caap_group", "trait", "position", "caas", "amino_encoded"])
    disc["hyp"] = disc["trait"].str.extract(r"(H\d+)")[0]
    key = ["gene", "caap_group", "hyp", "position", "caas", "amino_encoded"]
    obs_rows = set(map(tuple, disc[key].astype(str).itertuples(index=False, name=None)))

    work = Path(tempfile.mkdtemp(prefix="b0bg_", dir=a.tmp_dir))
    subprocess.run([sys.executable, str(b0_script), "--config", str(cfg), "--fop", "--labelings-out", str(work / "b0.tab"),
                    "--pairs-out", "/dev/null"], check=True, capture_output=True)
    bad_bg, k_rows, missing = [], set(), []
    for g in genes:
        aln = find_alignment(a.align_dir, g)
        if aln is None:
            missing.append(g)
            continue
        p = subprocess.run([str(ct), "perm-replay", "-a", str(aln), "-t", str(cfg), "-s", str(work / "b0.tab"),
                            "-o", str(work / f"{g}.out"), "--fmt", a.fmt, "--patterns", a.patterns, *flags, *thr,
                            "--export_perm_discovery", str(work / f"{g}.disc"), "--export_b0_background", str(work / f"{g}.bg")],
                           capture_output=True, text=True)
        if p.returncode != 0:
            print(f"[{g}] perm-replay failed:\n{p.stdout[-400:]}{p.stderr[-400:]}")
            bad_bg.append(g)
            continue
        _, _, pos = (work / f"{g}.bg").read_text().rstrip("\n").partition("\t")
        kern = set() if pos in ("", "NULL") else set(pos.split(","))
        if kern != obs_bg[g]:
            bad_bg.append(g)
            print(f"[{g}] background differs: only observed {sorted(obs_bg[g] - kern, key=int)[:8]}, only kernel {sorted(kern - obs_bg[g], key=int)[:8]}")
        d = pd.read_csv(work / f"{g}.disc", sep="\t")
        if len(d):
            d["hyp"] = d["cycle"].str.extract(r"(H\d+)")[0]
            k_rows |= set(map(tuple, d[key].astype(str).itertuples(index=False, name=None)))

    checked = [g for g in genes if g not in missing]
    obs_checked = {r for r in obs_rows if r[0] in set(checked)}
    print(f"genes: {len(genes)} in observed background, {len(checked)} checked, {len(missing)} without alignment {missing[:5]}")
    print(f"background: {len(checked) - len(bad_bg)}/{len(checked)} genes identical")
    print(f"discovery rows: observed {len(obs_checked)}, kernel {len(k_rows)}, only observed {len(obs_checked - k_rows)}, only kernel {len(k_rows - obs_checked)}")
    ok = not bad_bg and not missing and obs_checked == k_rows
    shutil.rmtree(work, ignore_errors=True)
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
