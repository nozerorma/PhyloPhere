#!/usr/bin/env python3
"""
build_vep_hgvs.py  —  Build protein-level HGVS identifiers for Ensembl VEP's
--format hgvs input, for the ancestral->derived amino-acid change at every
CAAS position.

Unlike map_to_primateai.py / map_to_cosmic.py, this does not consult any
pathogenicity database for the reference amino acid: the ancestral state is
already known from this pipeline's own ASR output (the same
top_residue_support/bottom_residue_support descriptor columns
CT_POSTPROC/local/src/residue_descriptors.py writes), so no external
reference-proteome dependency is introduced by adding this annotation source.

Coordinate translation (alignment column -> real hg38 protein position) still
requires the upstream vep_map_dir MAP files (see map_to_primateai.py's
docstring) — that mapping is not something this pipeline can derive on its
own, which is why vep_map_dir stays a required external input.

human_protein_id must be an Ensembl protein stable ID (ENSP...) for VEP's
--format hgvs cache lookup to resolve it — bin/generate_ensembl_mapping.py's
BioMart query (ensembl_peptide_id attribute) already produces IDs in that
form, so this composes correctly with an auto-generated gene_ensembl_file
(see §1). A user-supplied gene_ensembl_file using a different protein-ID
system (e.g. UniProt from an older cds2prot-based run) will not resolve and
that gene is silently skipped.

Reuses map_to_primateai.py's anc_der_from_descriptor()/load_map_file() as a
library rather than re-implementing the same position-matching logic twice.

Usage
-----
    build_vep_hgvs.py <caas_file> <vep_map_dir> <gene_ensembl_file> \
        <output_hgvs.txt> <output_id_map.tsv>

Output
------
  output_hgvs.txt   — one HGVS protein identifier per line, e.g.
                       ENSP00000234875.4:p.Trp24Cys
  output_id_map.tsv — HGVS identifier -> Gene, Position, caap_group, so VEP's
                       tab output (which echoes the input identifier back
                       under Uploaded_variation) can be rejoined afterwards.
"""

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from vep_common import anc_der_from_descriptor, load_map_file  # noqa: E402

_AA_3LETTER = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
}


def load_gene_protein_ids(gene_ensembl_file: str) -> dict:
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
    """(gene, position) -> list of {anc_aas, der_aas, caap_group}."""
    targets = {}
    with open(caas_file, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        reader.fieldnames = [c.strip() for c in reader.fieldnames]
        required = {"Gene", "Position", "derived_residues"}
        if not required.issubset(reader.fieldnames):
            sys.exit(
                f"Error: {caas_file} is missing required columns {sorted(required)}; "
                "the Ensembl VEP annotation source needs the modern CT_POSTPROC "
                "descriptor output (derived_residues/top_residue_support/"
                "bottom_residue_support/caap_group/caas_side)."
            )
        for row in reader:
            caap_group = row.get("caap_group", "US") or "US"
            if caap_group != "US":
                continue
            gene = row["Gene"]
            try:
                position = int(row["Position"])
            except ValueError:
                continue
            anc_aas, der_aas = anc_der_from_descriptor(
                row.get("derived_residues", ""),
                row.get("top_residue_support", ""),
                row.get("bottom_residue_support", ""),
                row.get("side", ""),
            )
            if not anc_aas or not der_aas:
                continue
            targets.setdefault((gene, position), []).append({
                "anc_aas": anc_aas, "der_aas": der_aas, "caap_group": caap_group,
            })
    return targets


def main():
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
