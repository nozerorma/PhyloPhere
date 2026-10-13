#!/usr/bin/env python3
# test_vep_side.py: side propagation through build_vep_hgvs.py and join_vep_output.py.
# Run: python3 tests/vep_side/test_vep_side.py   (scratch files go to tests/vep_side/_work, removed afterwards)

import csv
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE.parents[1] / "subworkflows" / "VEP" / "local" / "src"
WORK = HERE / "_work"


def write(path, rows, header):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def read(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    shutil.rmtree(WORK, ignore_errors=True)
    (WORK / "maps").mkdir(parents=True)
    # One gene, alignment column 11 (Position 10, 0-based) maps to protein position 24.
    write(WORK / "maps" / "G1.map.tsv",
          [["7:100", "23", "10"], ["7:103", "24", "11"], ["7:106", "25", "12"]],
          ["hg38_nt_coord", "hg38_aa_pos", "prot_ali_col"])
    write(WORK / "genes.tsv", [["G1", "ENSP1.1"]], ["gene", "human_protein_id"])
    # The same position on both sides with different derived residues, plus a second position whose two sides
    # share the same ancestral and derived residue (same HGVS identifier for both sides).
    write(WORK / "ps.tsv", [
        ["G1", "10", "top", "W/C", "W", "C", "US"],
        ["G1", "10", "bottom", "W/C", "W", "L", "US"],
        ["G1", "11", "top", "A/G", "A", "G", "US"],
        ["G1", "11", "bottom", "A/G", "A", "G", "US"],
    ], ["Gene", "Position", "side", "caas", "ancestral_aa", "derived_aa", "caap_group"])

    subprocess.run([sys.executable, str(SRC / "build_vep_hgvs.py"), str(WORK / "ps.tsv"), str(WORK / "maps"),
                    str(WORK / "genes.tsv"), str(WORK / "ids.txt"), str(WORK / "idmap.tsv")],
                   check=True, capture_output=True)
    ids = (WORK / "ids.txt").read_text().split()
    idmap = read(WORK / "idmap.tsv")
    assert list(idmap[0].keys()) == ["hgvs_id", "Gene", "Position", "side", "caap_group"], idmap[0].keys()
    # Position 10: Trp24Cys on top, Trp24Leu on bottom; position 11 maps to column 12 (protein position 25), Ala25Gly.
    assert sorted(ids) == sorted(set(ids)), "identifiers must be unique in the VEP input"
    assert {"ENSP1.1:p.Trp24Cys", "ENSP1.1:p.Trp24Leu", "ENSP1.1:p.Ala25Gly"} == set(ids), ids
    sides = {(r["hgvs_id"], r["Position"], r["side"]) for r in idmap}
    assert ("ENSP1.1:p.Trp24Cys", "10", "top") in sides and ("ENSP1.1:p.Trp24Leu", "10", "bottom") in sides
    assert ("ENSP1.1:p.Trp24Cys", "10", "bottom") not in sides
    # The shared identifier keeps one id-map row per side.
    assert {r["side"] for r in idmap if r["hgvs_id"] == "ENSP1.1:p.Ala25Gly"} == {"top", "bottom"}

    # Fake VEP output: two transcripts for the first identifier, one row for the shared identifier.
    (WORK / "vep.txt").write_text(
        "## ENSEMBL VARIANT EFFECT PREDICTOR\n"
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tExtra\n"
        "ENSP1.1:p.Trp24Cys\t7:103\tG\tENSG1\tENST1\tTranscript\tmissense_variant\tIMPACT=MODERATE;CANONICAL=YES\n"
        "ENSP1.1:p.Trp24Cys\t7:103\tG\tENSG1\tENST2\tTranscript\tsplice_region_variant\tIMPACT=LOW\n"
        "ENSP1.1:p.Ala25Gly\t7:106\tG\tENSG1\tENST1\tTranscript\tmissense_variant\tIMPACT=MODERATE\n")
    subprocess.run([sys.executable, str(SRC / "join_vep_output.py"), str(WORK / "vep.txt"), str(WORK / "idmap.tsv"),
                    str(WORK / "out.tsv")], check=True, capture_output=True)
    out = read(WORK / "out.tsv")
    assert list(out[0].keys())[:5] == ["Gene", "Position", "side", "caap_group", "Uploaded_variation"], list(out[0].keys())
    assert "Ensembl_Gene" in out[0]
    got = sorted((r["Position"], r["side"], r["Feature"]) for r in out)
    # Trp24Cys rows belong to the top side only; the shared Ala25Gly row is emitted once per side.
    assert got == [("10", "top", "ENST1"), ("10", "top", "ENST2"), ("11", "bottom", "ENST1"), ("11", "top", "ENST1")], got
    shutil.rmtree(WORK)
    print("test_vep_side: OK")


if __name__ == "__main__":
    main()
