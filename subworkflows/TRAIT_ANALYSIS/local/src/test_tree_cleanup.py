#!/usr/bin/env python3
# test_tree_cleanup.py — Unit tests of the species curation in tree_cleanup.py.
# PhyloPhere | subworkflows/TRAIT_ANALYSIS/local/src/
#
# Run: python3 -m pytest subworkflows/TRAIT_ANALYSIS/local/src/test_tree_cleanup.py

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from argparse import Namespace  # noqa: E402

from tree_cleanup import curate_traits, resolve_tip_taxids, write_species_outputs  # noqa: E402


def test_shared_tax_id_first_alphabetical_keeps_it():
    # Two canonical species share the real tax_id 10; 11 is taken by another name of the map.
    name_to_taxid = {"B_sp": "10", "A_sp": "10", "C_sp": "12", "X_syn": "11"}
    out = resolve_tip_taxids(["B_sp", "A_sp", "C_sp"], name_to_taxid, reserved={"10", "11", "12"})
    assert out["A_sp"][:2] == ("10", "10")
    # 10 + 1 = 11 and 12 are reserved, so the synthetic id probes forward to 13
    assert out["B_sp"][:2] == ("10", "13")
    assert "synthetic" in out["B_sp"][2]
    assert out["C_sp"][:2] == ("12", "12")


def test_species_without_tax_id_keeps_empty_id():
    out = resolve_tip_taxids(["A_sp"], {}, reserved=set())
    assert out["A_sp"] == ("", "", "no tax_id in the map")


def _decisions(rows, tips, name_to_taxid):
    header = ["species", "tax_id", "trait"]
    resolved = resolve_tip_taxids(sorted(tips), name_to_taxid, {v for v in name_to_taxid.values()})
    kept, dec = curate_traits(header, rows, "species", set(tips), name_to_taxid, resolved)
    return kept, {d[0]: d for d in dec}


def test_name_match_synonym_removal_and_duplicate():
    name_to_taxid = {"Tip_a": "1", "Syn_a": "1", "Tip_b": "2", "Orphan": "9"}
    rows = [
        ["Syn_a", "1", "0.1"],    # synonym of Tip_a, listed before the tip itself
        ["Tip_a", "1", "0.2"],    # named like the tip: wins over the synonym
        ["Tip_b", "2", "0.3"],    # plain match
        ["Orphan", "9", "0.4"],   # tax_id with no tip
        ["Unknown", "", "0.5"],   # no tax_id and no tip
    ]
    kept, dec = _decisions(rows, {"Tip_a", "Tip_b"}, name_to_taxid)
    assert [r[0] for r in kept] == ["Tip_a", "Tip_b"]
    assert dec["Tip_a"][2] == "maintained"
    assert dec["Tip_b"][2] == "maintained"
    assert dec["Syn_a"][2] == "removed" and "duplicate" in dec["Syn_a"][3]
    assert dec["Orphan"][2] == "removed" and "tax_id 9" in dec["Orphan"][3]
    assert dec["Unknown"][2] == "removed"


def test_synonym_is_renamed_when_no_row_has_the_tip_name():
    name_to_taxid = {"Tip_a": "1", "Syn_a": "1"}
    kept, dec = _decisions([["Syn_a", "1", "0.1"]], {"Tip_a"}, name_to_taxid)
    assert [r[0] for r in kept] == ["Tip_a"]
    assert dec["Syn_a"][2] == "changed" and dec["Syn_a"][1] == "Tip_a"


def test_tax_id_column_takes_the_resolved_id():
    # Two tips share the real id 1; the trait rows carry the real id and must end up distinct.
    name_to_taxid = {"A_sp": "1", "B_sp": "1"}
    kept, _ = _decisions([["A_sp", "1", "0"], ["B_sp", "1", "0"]], {"A_sp", "B_sp"}, name_to_taxid)
    ids = {r[0]: r[1] for r in kept}
    assert ids["A_sp"] == "1" and ids["B_sp"] != "1"


def test_curated_taxid_map_has_one_unique_id_per_species(tmp_path):
    name_to_taxid = {"A_sp": "1", "B_sp": "1", "C_sp": "5"}
    tips = ["A_sp", "B_sp", "C_sp"]
    resolved = resolve_tip_taxids(tips, name_to_taxid, set(name_to_taxid.values()))
    args = Namespace(traits_out=str(tmp_path / "t.csv"), species_table=str(tmp_path / "s.tsv"),
                     species_report=str(tmp_path / "r.txt"), taxid_map=str(tmp_path / "m.tsv"))
    write_species_outputs(args, ["species"], [["A_sp"]], [("A_sp", "A_sp", "maintained", "x")],
                          [("A_sp", "A_sp", "maintained", "x")], resolved, ",", {"A_sp": "Fam_a"})
    rows = [line.split("\t") for line in (tmp_path / "m.tsv").read_text().splitlines()]
    # the layout of the taxonomy file: the clade variability reads the family (column 3) and the name class (column 5)
    assert rows[0] == ["tax_id", "species", "family", "rank", "name_class"]
    assert rows[1][2] == "Fam_a" and all(r[4] == "scientific name" for r in rows[1:])
    ids = [r[0] for r in rows[1:]]
    assert sorted(r[1] for r in rows[1:]) == tips and len(set(ids)) == 3
    assert "synthetic tax_ids: 1" in (tmp_path / "r.txt").read_text()


def test_tree_only_curation_writes_the_species_outputs(tmp_path):
    # Without --traits the species table, the map and the report are still written.
    (tmp_path / "t.nwk").write_text("((A_sp:1,B_syn:1):1,C_sp:1);\n")
    (tmp_path / "ali.txt").write_text("A_sp\nC_sp\n")
    (tmp_path / "tax.tsv").write_text("tax_id\tspecies\tfamily\n1\tA_sp\tFam\n1\tB_syn\tFam\n3\tC_sp\tFam2\n")
    import subprocess
    out = subprocess.run(
        [sys.executable, str(Path(__file__).parent / "tree_cleanup.py"),
         "--tree", str(tmp_path / "t.nwk"), "--ali-sp-names", str(tmp_path / "ali.txt"),
         "--tax-id", str(tmp_path / "tax.tsv"), "--output", str(tmp_path / "o.nwk"),
         "--report", str(tmp_path / "rep.tsv"), "--species-table", str(tmp_path / "s.tsv"),
         "--species-report", str(tmp_path / "r.txt"), "--taxid-map", str(tmp_path / "m.tsv")],
        capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert not (tmp_path / "curated_traits.csv").exists()
    # B_syn shares tax_id 1 with the alignment species A_sp: it is pruned as a duplicate
    assert "B_syn" not in (tmp_path / "o.nwk").read_text()
    assert [l.split("\t")[1] for l in (tmp_path / "m.tsv").read_text().splitlines()[1:]] == ["A_sp", "C_sp"]
