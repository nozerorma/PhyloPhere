#!/usr/bin/env python3
"""
join_vep_output.py  —  Rejoin Ensembl VEP's --tab consequence output to the
originating CAAS (Gene, Position, caap_group) via the id-map TSV
build_vep_hgvs.py wrote (VEP echoes the input HGVS identifier back verbatim
under Uploaded_variation).

Usage
-----
    join_vep_output.py <vep_tab_output> <id_map.tsv> <output_tsv>
"""

import csv
import sys


def load_id_map(path: str) -> dict:
    mapping = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            mapping[row["hgvs_id"]] = (row["Gene"], row["Position"], row["caap_group"])
    return mapping


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
