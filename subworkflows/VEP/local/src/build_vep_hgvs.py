#!/usr/bin/env python3
# build_vep_hgvs.py — Protein-level HGVS identifiers for Ensembl VEP, one per ancestral→derived change.
# PhyloPhere | subworkflows/VEP/local/src/

"""
BuildVepHgvs: writes the HGVS protein identifier (VEP's --format hgvs input) of the
ancestral→derived amino-acid change at every position of a scored table, and the table
that rejoins VEP's output to the rows.

The ancestral and derived residues are read from the species residue tallies of the table
(top_species_residues / bottom_species_residues, written by
CT_POSTPROC/local/src/residue_descriptors.py) and from the side, with the caas pattern as
fallback (vep_common.anc_der_from_descriptor). No reference proteome is needed.
Translating an alignment column to a protein position needs the vep_map_dir MAP files,
which this pipeline cannot derive. human_protein_id must be an Ensembl protein ID (ENSP...)
for VEP's cache lookup; a gene with another ID system is skipped without a message.

Called by:  ENSEMBL_VEP_ANNOTATE Nextflow process (ensembl_vep.nf → build_vep_hgvs.py)
Inputs:     caas_file          position_scores.tsv, with Gene, Position and the residue tally columns
            vep_map_dir        directory of the per-gene MAP files (alignment column → hg38 position)
            gene_ensembl_file  TSV with gene and human_protein_id (ENSP...)
Outputs:    output_hgvs.txt    one identifier per line, e.g. ENSP00000234875.4:p.Trp24Cys
            output_id_map.tsv  hgvs_id, Gene, Position, caap_group
"""

# ── Standard library ──────────────────────────────────────────────────────────
import csv
import sys
from pathlib import Path

# ── Package-internal ──────────────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent))
from vep_common import anc_der_from_descriptor, load_map_file  # noqa: E402


# ── Constants ─────────────────────────────────────────────────────────────────

_AA_3LETTER = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
}


# ── Inputs ────────────────────────────────────────────────────────────────────


def load_gene_protein_ids(gene_ensembl_file: str) -> dict:
    """Gene → Ensembl protein ID, skipping genes whose ID is empty or NA."""
    mapping = {}
    with open(gene_ensembl_file, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        reader.fieldnames = [c.lower() for c in reader.fieldnames]
        for row in reader:
            protein_id = row.get("human_protein_id", "").strip()
            gene = row.get("gene", "").strip()
            if gene and protein_id and protein_id.lower() not in ("", "na", "none"):
                mapping[gene] = protein_id
    return mapping


def load_caas_targets(caas_file: str) -> dict:
    """(gene, position) → list of {anc_aas, der_aas, caap_group}, for the US scheme rows only."""
    targets = {}
    with open(caas_file, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            return targets
        reader.fieldnames = [c.strip() for c in reader.fieldnames]
        col_lc = {c.lower(): c for c in reader.fieldnames}
        gene_col = col_lc.get("gene")
        pos_col = col_lc.get("position")
        if not gene_col or not pos_col:
            sys.exit(
                f"Error: {caas_file} is missing required columns 'Gene' and 'Position'."
            )
        top_col = col_lc.get("top_species_residues") or col_lc.get("top_residue_support")
        bot_col = col_lc.get("bottom_species_residues") or col_lc.get("bottom_residue_support")
        side_col = col_lc.get("side")
        caas_col = col_lc.get("caas")
        caap_col = col_lc.get("caap_group")

        for row in reader:
            caap_group = (row.get(caap_col, "US") if caap_col else "US") or "US"
            if caap_group != "US":
                continue
            gene = row[gene_col]
            try:
                position = int(row[pos_col])
            except ValueError:
                continue
            top_res = row.get(top_col, "") if top_col else ""
            bot_res = row.get(bot_col, "") if bot_col else ""
            side = row.get(side_col, "") if side_col else ""
            caas = row.get(caas_col, "") if caas_col else ""
            anc_aas, der_aas = anc_der_from_descriptor(top_res, bot_res, side, caas=caas)
            if not anc_aas or not der_aas:
                continue
            targets.setdefault((gene, position), []).append({
                "anc_aas": anc_aas, "der_aas": der_aas, "caap_group": caap_group,
            })
    return targets


# ── CLI ───────────────────────────────────────────────────────────────────────


def main():
    """Write the HGVS list and the id map for every CAAS position with a protein ID and a MAP entry."""
    if len(sys.argv) != 6:
        sys.exit(f"Usage: {sys.argv[0]} <caas_file> <vep_map_dir> <gene_ensembl_file> "
                  "<output_hgvs.txt> <output_id_map.tsv>")
    caas_file, vep_map_dir, gene_ensembl_file, output_hgvs, output_id_map = sys.argv[1:]

    protein_ids = load_gene_protein_ids(gene_ensembl_file)
    targets = load_caas_targets(caas_file)

    loaded_maps = {}
    n_written = 0
    with open(output_hgvs, "w") as hgvs_out, open(output_id_map, "w", newline="") as map_out:
        map_writer = csv.writer(map_out, delimiter="\t", lineterminator="\n")
        map_writer.writerow(["hgvs_id", "Gene", "Position", "caap_group"])

        for (gene, position), entries in targets.items():
            protein_id = protein_ids.get(gene)
            if not protein_id:
                continue
            if gene not in loaded_maps:
                loaded_maps[gene] = load_map_file(gene, vep_map_dir)
            pos_map, _strand = loaded_maps[gene]
            if not pos_map:
                continue
            prot_ali_idx = position + 1  # CAAS Position is 0-based; MAP file is 1-based
            hit = pos_map.get(prot_ali_idx)
            if not hit:
                continue
            hg38_aa_pos, _chrom, _coord = hit

            for entry in entries:
                for anc_aa in sorted(entry["anc_aas"]):
                    anc3 = _AA_3LETTER.get(anc_aa)
                    if not anc3:
                        continue
                    for der_aa in sorted(entry["der_aas"] - {anc_aa}):
                        der3 = _AA_3LETTER.get(der_aa)
                        if not der3:
                            continue
                        hgvs_id = f"{protein_id}:p.{anc3}{hg38_aa_pos}{der3}"
                        hgvs_out.write(hgvs_id + "\n")
                        map_writer.writerow([hgvs_id, gene, position, entry["caap_group"]])
                        n_written += 1

    print(f"Wrote {n_written} HGVS identifiers to {output_hgvs}", file=sys.stderr)


if __name__ == "__main__":
    main()
