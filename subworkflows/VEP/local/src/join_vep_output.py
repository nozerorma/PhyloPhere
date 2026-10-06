#!/usr/bin/env python3
# join_vep_output.py — Attach Ensembl VEP consequences to the CAAS positions they were computed for.
# PhyloPhere | subworkflows/VEP/local/src/

"""
JoinVepOutput: rejoins the --tab output of Ensembl VEP to the CAAS position that
produced each input identifier. VEP echoes the input HGVS identifier verbatim in the
first column (Uploaded_variation), which is the key of the id map written by
build_vep_hgvs.py. VEP rows whose identifier is not in the map are dropped.

Called by:  ENSEMBL_VEP_ANNOTATE Nextflow process (ensembl_vep.nf → join_vep_output.py)
Inputs:     vep_tab_output  VEP --tab file (## metadata lines, a #-prefixed header, rows)
            id_map.tsv      hgvs_id, Gene, Position, caap_group
Outputs:    output_tsv      Gene, Position, caap_group, then the VEP columns
"""

# ── Standard library ──────────────────────────────────────────────────────────
import csv
import sys


# ── Functions ─────────────────────────────────────────────────────────────────


def load_id_map(path: str) -> dict:
    """Return {hgvs_id: (Gene, Position, caap_group)} from the id-map TSV."""
    mapping = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            mapping[row["hgvs_id"]] = (row["Gene"], row["Position"], row["caap_group"])
    return mapping


# ── CLI ───────────────────────────────────────────────────────────────────────


def main():
    if len(sys.argv) != 4:
        sys.exit(f"Usage: {sys.argv[0]} <vep_tab_output> <id_map.tsv> <output_tsv>")
    vep_tab, id_map_path, output_tsv = sys.argv[1:]

    id_map = load_id_map(id_map_path)

    n_written = 0
    with open(vep_tab) as fh, open(output_tsv, "w", newline="") as out:
        header = None
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line[1:].rstrip("\n").split("\t")
                out.write("Gene\tPosition\tcaap_group\t" + "\t".join(header) + "\n")
                continue
            if header is None:
                continue
            fields = line.rstrip("\n").split("\t")
            hgvs_id = fields[0]
            hit = id_map.get(hgvs_id)
            if not hit:
                continue
            gene, position, caap_group = hit
            out.write(f"{gene}\t{position}\t{caap_group}\t" + "\t".join(fields) + "\n")
            n_written += 1

    print(f"Wrote {n_written} annotated rows to {output_tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()
