"""The merge of a block keeps the table of units the gene filter decided on (units.parquet) and the thresholds of every
cycle (cycle_thresholds.tsv), and the decisions are the ones removed_units.tsv has always listed."""
import subprocess
import sys

import pyarrow.parquet as pq

from pathlib import Path

PIPELINE = Path(__file__).resolve().parents[1]
LOCAL = PIPELINE / "subworkflows/CT_DISAMBIGUATION/local"
sys.path.insert(0, str(LOCAL))

from src.utils.null_tables import CycleReducer, concat_tables, gene_tables, write_tables  # noqa: E402
from src.utils.perm_stream import SlimRecord  # noqa: E402


def _block(tmp_path):
    units, lengths = {}, {}
    for i in range(60):                                   # ordinary genes: 1 to 3 scattered positions per cycle
        gene = f"G{i:02d}"
        lengths[gene] = 300 + 10 * i
        for cycle in ("b_1", "b_2"):
            units[(gene, cycle)] = [(5 + 40 * j, 0.3, "top") for j in range(1 + (i + (cycle == "b_2")) % 3)]
    lengths.update({"BIGCL": 900, "SHORTX": 40})
    units[("BIGCL", "b_1")] = [(100 + j, 0.5, "top") for j in range(30)]     # many positions in a train: dubious
    units[("SHORTX", "b_1")] = [(2 + j, 0.5, "top") for j in range(6)]       # few positions in 40 aa: extreme
    parts = {}
    for (gene, cycle), recs in units.items():
        parts.setdefault(gene, []).append(CycleReducer(True, True, 3, 0.7, None)(cycle, [SlimRecord(p, s, side, 1, 1) for p, s, side in recs]))
    write_tables(tmp_path / "tile", *concat_tables([gene_tables(g, r) for g, r in parts.items()]), meta={"by_labeling": False})
    (tmp_path / "lengths.tsv").write_text("gene\tlength\n" + "".join(f"{g}\t{n}\n" for g, n in lengths.items()))
    (tmp_path / "lab.tab").write_text("b_1\tsp1,sp2\tsp3,sp4\nb_2\tsp1,sp2\tsp3,sp4\n")
    r = subprocess.run([sys.executable, str(LOCAL / "merge_null_block.py"), "--tiles", str(tmp_path / "tile"), "--labelings", str(tmp_path / "lab.tab"),
                        "--output-dir", str(tmp_path / "out"), "--gene-lengths", str(tmp_path / "lengths.tsv"), "--gene-filter-mode", "both"],
                       cwd=LOCAL, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-800:]
    return tmp_path / "out"


def test_units_table_keeps_every_unit_with_its_category(tmp_path):
    out = _block(tmp_path)
    t = pq.read_table(out / "units.parquet").to_pandas()
    assert list(t.columns) == ["gene", "cycle", "n_pos", "clust", "category"]
    assert len(t) == 60 * 2 + 2                                                # 60 genes x 2 cycles, BIGCL and SHORTX in b_1
    row = t[(t.gene == "BIGCL") & (t.cycle == "b_1")].iloc[0]
    assert (row.n_pos, bool(row.clust), row.category) == (30, True, "Dubious")
    row = t[t.gene == "SHORTX"].iloc[0]
    assert (row.n_pos, bool(row.clust), row.category) == (6, True, "Extreme")
    # the categories are exactly the units of removed_units.tsv
    listed = {tuple(line.split("\t")) for line in (out / "removed_units.tsv").read_text().splitlines()[1:]}
    flagged = {(c, g) for g, c, k in zip(t.gene, t.cycle, t.category) if isinstance(k, str)}
    assert flagged == listed


def test_cycle_thresholds_list_each_cycle_with_its_units(tmp_path):
    out = _block(tmp_path)
    rows = [line.split("\t") for line in (out / "cycle_thresholds.tsv").read_text().splitlines()]
    assert rows[0] == ["cycle", "n_units", "thr_count", "thr_density"]
    assert [r[0] for r in rows[1:]] == ["b_1", "b_2"]
    assert [int(r[1]) for r in rows[1:]] == [62, 60]
