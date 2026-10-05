"""explain_positions.py: the evidence table of the N best positions of a run.

Same frozen PEPC inputs as test_observed_b0_cli (its `inp` fixture); the positions come from the frozen
position_scores.tsv and the reference for the numbers is the frozen caas_convergence_master.csv.
"""
import csv
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("PHYLOPHERE_ROOT", HERE.parents[1]))
LOCAL = ROOT / "subworkflows/CT_DISAMBIGUATION/local"
MAIN = LOCAL / "explain_positions.py"
GOLD = HERE / "golden/pepc_c4_complete"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(LOCAL))
from test_observed_b0_cli import inp  # noqa: E402,F401  (the PEPC inputs)
from src.core.evidence import EVIDENCE_COLUMNS, select_top_positions  # noqa: E402


def _scores(path=GOLD / "position_scores.tsv"):
    return list(csv.DictReader(open(path), delimiter="\t"))


def _run(inp, out, *extra, scores=GOLD / "position_scores.tsv", discovery=None, workers="2", ensembl=None):
    i = inp / "observed_inputs"
    cmd = [sys.executable, str(MAIN), "--alignment-dir", str(inp / "align"), "--tree", str(i / "pruned_tree_file.nwk"),
           "--discovery", str(discovery or inp / "b0/PEPC.b0.discovery.tsv"), "--position-scores", str(scores),
           "--design", str(i / "traitfiles"), "--output-dir", str(out),
           "--asr-model", "lg", "--posterior-threshold", "0.1", "--workers", workers, "--asr-cache-dir", str(i / "asr_cache"),
           "--taxid-mapping", str(i / "taxid.tsv"), "--ensembl-genes-file", str(ensembl or i / "gene_ensembl.tsv"), *extra]
    return subprocess.run(cmd, capture_output=True, text=True)


def _table(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def test_the_table_explains_the_chosen_positions_and_says_how_they_were_chosen(inp, tmp_path):
    p = _run(inp, tmp_path / "out", "--top", "2")
    assert p.returncode == 0, p.stdout[-1500:] + p.stderr[-1500:]
    assert sorted(f.name for f in (tmp_path / "out").iterdir()) == ["evidence_top2.tsv", "top_positions.tsv"]
    chosen = select_top_positions(_scores(), 2)
    top = _table(tmp_path / "out/top_positions.tsv")
    assert list(top.columns) == ["gene", "position", "CAAS_score", "p.emp"]
    assert [(r.gene, r.position) for r in top.itertuples()] == [(g, pos) for g, pos, _ in chosen]
    assert [float(r.CAAS_score) for r in top.itertuples()] == [info["score"] for _, _, info in chosen]
    assert top["p.emp"].tolist() == [repr(info["p_emp"]) if info["p_emp"] is not None else "NA" for _, _, info in chosen]
    assert any(v != "NA" for v in top["p.emp"])
    ev = _table(tmp_path / "out/evidence_top2.tsv")
    assert list(ev.columns) == EVIDENCE_COLUMNS
    assert {(g, m) for g, m in zip(ev.gene, ev.msa_pos)} == {(g, pos) for g, pos, _ in chosen}
    assert set(ev.domain) == {"1", "2", "3", "4"}


def test_the_domain_scores_of_the_table_pool_to_the_frozen_master(inp, tmp_path):
    """The table must not invent numbers: its rows, averaged over the hypotheses of a position, give the master's."""
    p = _run(inp, tmp_path / "out", "--top", "2")
    assert p.returncode == 0, p.stdout[-1500:] + p.stderr[-1500:]
    ev = _table(tmp_path / "out/evidence_top2.tsv")
    master = pd.read_csv(GOLD / "caas_convergence_master.csv", keep_default_na=False)
    checked = 0
    for m in master[master.msa_pos.astype(str).isin(set(ev.msa_pos)) & master.side.isin(["top", "bottom"])].itertuples():
        sel = ev[(ev.msa_pos == str(m.msa_pos)) & (ev.caap_group == m.caap_group)]
        for d in range(1, 5):
            vals = [float(v or 0.0) for v in sel[sel.domain == str(d)][f"{m.side}_domain_score"]]
            got = getattr(m, f"domain_{d}_score")
            assert (got if got != "" else 0.0) == pytest.approx(sum(vals) / len(vals), abs=1e-12), (m.msa_pos, m.caap_group, m.side, d)
            checked += 1
    assert checked >= 8


def test_top_zero_writes_the_two_headers_and_nothing_else(inp, tmp_path):
    p = _run(inp, tmp_path / "out", "--top", "0")
    assert p.returncode == 0, p.stderr[-800:]
    assert (tmp_path / "out/evidence_top0.tsv").read_text() == "\t".join(EVIDENCE_COLUMNS) + "\n"
    assert (tmp_path / "out/top_positions.tsv").read_text() == "gene\tposition\tCAAS_score\tp.emp\n"


def test_the_output_does_not_depend_on_the_number_of_workers_and_follows_the_ranking(inp, tmp_path):
    ev = []
    for w in ("1", "2"):
        p = _run(inp, tmp_path / w, "--top", "5", workers=w)
        assert p.returncode == 0, p.stderr[-800:]
        ev.append((tmp_path / w / "evidence_top5.tsv").read_bytes())
    assert ev[0] == ev[1] and ev[0].count(b"\n") > 12
    chosen = [pos for _, pos, _ in select_top_positions(_scores(), 5)]
    assert chosen != sorted(chosen, key=int), "the fixture must rank positions out of position order"
    shown = _table(tmp_path / "1/evidence_top5.tsv").msa_pos.tolist()
    assert list(dict.fromkeys(shown)) == chosen


def test_a_chosen_gene_without_alignment_fails_the_run_but_the_rest_is_written(inp, tmp_path):
    disc = (inp / "b0/PEPC.b0.discovery.tsv").read_text().splitlines()
    (tmp_path / "disc.tab").write_text("\n".join(disc + [r.replace("PEPC\t", "GHOST\t", 1) for r in disc[1:20]]) + "\n")
    rows = _scores()
    ghost = {**rows[0], "Gene": "GHOST", "CAAS_score": "9.0"}
    with open(tmp_path / "scores.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows([ghost] + rows)
    i = inp / "observed_inputs"
    (tmp_path / "ens.tsv").write_text((i / "gene_ensembl.tsv").read_text() + "GHOST\tchr1\t1\t2\t+\t970\tP04711\n")
    p = _run(inp, tmp_path / "out", "--top", "2", scores=tmp_path / "scores.tsv", discovery=tmp_path / "disc.tab", ensembl=tmp_path / "ens.tsv")
    assert p.returncode == 1, p.stdout[-1500:] + p.stderr[-1500:]
    assert "1 left out" in p.stderr + p.stdout and "GHOST" in p.stderr + p.stdout
    assert _table(tmp_path / "out/top_positions.tsv").gene.tolist()[0] == "GHOST"
    assert "incomplete" in p.stderr + p.stdout
    assert set(_table(tmp_path / "out/evidence_top2.tsv").gene) == {"PEPC"}


def test_a_chosen_position_without_discovery_rows_fails_the_run_and_is_reported(inp, tmp_path):
    rows = _scores()
    with open(tmp_path / "scores.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows([{**rows[0], "Position": "99999", "CAAS_score": "9.0"}] + rows)
    p = _run(inp, tmp_path / "out", "--top", "2", scores=tmp_path / "scores.tsv")
    assert p.returncode == 1, p.stdout[-1500:] + p.stderr[-1500:]
    assert "99999" in p.stderr + p.stdout
    assert "99999" not in set(_table(tmp_path / "out/evidence_top2.tsv").msa_pos)


def test_a_chosen_gene_with_no_discovery_rows_fails_the_run(inp, tmp_path):
    rows = _scores()
    with open(tmp_path / "scores.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows([{**rows[0], "Gene": "NODISC", "CAAS_score": "9.0"}] + rows)
    p = _run(inp, tmp_path / "out", "--top", "2", scores=tmp_path / "scores.tsv")
    assert p.returncode == 1, p.stdout[-1500:] + p.stderr[-1500:]
    assert "1 left out" in p.stderr + p.stdout and "NODISC" in p.stderr + p.stdout


def test_a_gene_outside_the_ensembl_list_fails_the_run(inp, tmp_path):
    (tmp_path / "ens.tsv").write_text("gene\tchr\tstart\tend\tstrand\tlength\thuman_protein_id\nOTHER\tchr1\t1\t2\t+\t970\tP04711\n")
    p = _run(inp, tmp_path / "out", "--top", "2", ensembl=tmp_path / "ens.tsv")
    assert p.returncode == 1, p.stdout[-1500:] + p.stderr[-1500:]
    assert "1 left out" in p.stderr + p.stdout and "PEPC" in p.stderr + p.stdout
    assert len(_table(tmp_path / "out/evidence_top2.tsv")) == 0


def test_only_the_schemes_that_were_scored_are_explained(inp, tmp_path):
    """position_scores.tsv lists the schemes a position was scored with (scheme_set, per side); discovery.tab can hold
    more of them (a removed unit), and the table must not explain what the score did not use."""
    rows, disc = _scores(), _table(inp / "b0/PEPC.b0.discovery.tsv")
    scored = {}
    for r in rows:
        scored.setdefault(r["Position"], set()).update(r["scheme_set"].split("+"))
    beyond = sorted(pos for pos in scored if set(disc[disc.position == pos].caap_group) > scored[pos])
    assert beyond, "the fixture needs a position whose discovery rows go beyond the scored schemes"
    with open(tmp_path / "scores.tsv", "w") as fh:        # put one such position and the best ones on top
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows([{**r, "CAAS_score": "9.0"} if r["Position"] == beyond[0] else r for r in rows])
    p = _run(inp, tmp_path / "out", "--top", "3", scores=tmp_path / "scores.tsv")
    assert p.returncode == 0, p.stderr[-800:]
    ev = _table(tmp_path / "out/evidence_top3.tsv")
    assert beyond[0] in set(ev.msa_pos)
    for pos, g in ev.groupby("msa_pos"):
        assert set(g.caap_group) == scored[pos], pos


def test_a_missing_scheme_set_does_not_filter(inp, tmp_path):
    rows = _scores()
    with open(tmp_path / "scores.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows([{**r, "scheme_set": "NA"} for r in rows])
    p = _run(inp, tmp_path / "out", "--top", "2", scores=tmp_path / "scores.tsv")
    assert p.returncode == 0, p.stderr[-800:]
    disc = _table(inp / "b0/PEPC.b0.discovery.tsv")
    ev = _table(tmp_path / "out/evidence_top2.tsv")
    for pos, g in ev.groupby("msa_pos"):
        assert set(g.caap_group) == set(disc[disc.position == pos].caap_group)
