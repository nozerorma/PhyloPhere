"""B0 gene removal of the null: the streaming implementation gives the removal set of the in-memory reference.

Run: python3 -I validation/unification/test_gene_removal_streaming.py
"""
import gzip
import random
import sys
import tempfile
from pathlib import Path

LOCAL = Path(__file__).resolve().parents[2] / "subworkflows" / "CT_DISAMBIGUATION" / "local"
sys.path.insert(0, str(LOCAL))

from src.core.postproc import GeneUnit, gene_unit_stats  # noqa: E402
from src.utils.gene_wrapper import _cycle_gene_removal_from_detail, iter_detail_rows  # noqa: E402

HEADER = ["Gene", "cycle", "Position", "asr_path_score", "n_detected", "clust", "side", "n_hyp", "n_sp"]


def make_detail(root: Path, rng: random.Random, n_genes=60, n_cycles=25):
    """One shard per gene; a few genes have many positions with a train so that dubious units exist."""
    lengths = {}
    for g in range(n_genes):
        gene = f"G{g:03d}"
        lengths[gene] = float(rng.randint(200, 3000))
        heavy = g % 7 == 0
        rows = []
        for c in range(1, n_cycles + 1):
            n = rng.randint(40, 90) if heavy and rng.random() < 0.5 else rng.randint(0, 6)
            for pos in rng.sample(range(1, 800), n):
                clust = 1 if (heavy and rng.random() < 0.3) else 0
                rows.append([gene, f"b_{c}", pos, 0.5, 5, clust, "none", 1, ""])
        with gzip.open(root / f"{gene}.tsv.gz", "wt", newline="") as fh:
            fh.write("\t".join(HEADER) + "\n")
            fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)
    return lengths


def reference(detail: Path, lengths, mode, iqr, ext):
    """The original B0: every detected position held in memory, stats dict per unit."""
    seen, clust = {}, {}
    for row in iter_detail_rows(detail):
        key = (row["cycle"], row["Gene"])
        seen.setdefault(key, set()).add(int(row["Position"]))
        if int(row["clust"]):
            clust[key] = True
    units = [GeneUnit(c, ge, len(p), clust.get((c, ge), False)) for (c, ge), p in seen.items()]
    out = set()
    for s in gene_unit_stats(units, lengths, iqr, ext):
        dub = mode in ("dubious", "both") and s["dubious"]
        exm = mode in ("extreme", "both") and s["extreme"]
        if dub or exm:
            out.add((s["labeling"], s["Gene"]))
    return out, sorted({u.labeling for u in units}), sum(1 for _ in iter_detail_rows(detail))


def main():
    rng = random.Random(7)
    with tempfile.TemporaryDirectory(dir=Path.home() / ".cc") as tmp:
        detail = Path(tmp) / "perm_pos_detail"
        detail.mkdir()
        lengths = make_detail(detail, rng)
        for mode in ("dubious", "extreme", "both"):
            ref, cycles, n_rows = reference(detail, lengths, mode, 3.0, 0.99)
            info = {}
            got = _cycle_gene_removal_from_detail(detail, lengths, mode, 3.0, 0.99, info=info)
            assert got == ref, f"{mode}: {len(got)} vs {len(ref)}"
            assert sorted(info["cycles"]) == cycles and info["n_rows"] == n_rows, mode
            print(f"{mode}: {len(got)} removed units identical, {n_rows} rows, {len(cycles)} cycles")
            assert mode == "extreme" or got, f"{mode}: the fixture should remove something"

        # A gene split over two shards is refused: its position sets would be counted apart.
        with gzip.open(detail / "A_split.tsv.gz", "wt", newline="") as fh:
            fh.write("\t".join(HEADER) + "\n")
            fh.write("\t".join(map(str, ["G003", "b_1", 5, 0.5, 5, 0, "none", 1, ""])) + "\n")
        try:
            _cycle_gene_removal_from_detail(detail, lengths, "dubious", 3.0, 0.99)
        except ValueError as e:
            print("non-contiguous gene refused:", e)
        else:
            raise AssertionError("non-contiguous gene accepted")
    print("OK")


if __name__ == "__main__":
    main()
