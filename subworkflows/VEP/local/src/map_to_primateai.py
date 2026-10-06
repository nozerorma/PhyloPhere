#!/usr/bin/python3
# map_to_primateai.py — Map CAAS positions to PrimateAI-3D pathogenicity scores.
# PhyloPhere | subworkflows/VEP/local/src/

"""
MapToPrimateai: For every CAAS position (Gene, Position), finds the PrimateAI-3D
missense variants of its hg38 codon that reproduce the ancestral→derived change of
the CAAS, in either orientation, and writes them with their scores.

Strategy:
  1. MAP file lookup. The gene's MAP file gives, for the 1-based alignment column
     `prot_ali_col` (= Position + 1, since Position is 0-based), the protein position
     `hg38_aa_pos` and the genomic coordinate `hg38_nt_coord`.
  2. Strand inference. Coordinates that increase down the MAP file mean the plus
     strand (codon nucleotides C, C+1, C+2); decreasing coordinates mean the minus
     strand (C, C-1, C-2).
  3. Database scan. The PrimateAI-3D table (gzip) is streamed once and the variants at
     any of the three codon nucleotides are kept when `alt_aa` is a residue of the CAAS
     that differs from the human one.

Orientation: the human reference residue of PrimateAI-3D (`ref_aa`) decides which
residues are valid alternatives. If it is an ancestral residue of the CAAS the
alternatives are the derived ones; if it is a derived residue they are the ancestral
ones; if it is neither, any CAAS residue other than the human one. When that leaves no
residue, every CAAS residue is accepted. Only rows of the `US` caap_group are used
(rows without a caap_group column count as `US`).

Called by:  PRIMATEAI_MAP Nextflow process (primateai.nf → map_to_primateai.py)
Inputs:     caas_file       position_scores.tsv, with Gene, Position and, when present, tag,
                            caas, side, amino_encoded, caap_group, ancestral_aa and derived_aa; the ancestral and
                            derived residues are those of the ASR, or the caas pattern read by side when absent
            vep_map_dir     directory of the per-gene MAP files
            primateai_gz    PrimateAI-3D hg38 table (gzip TSV, with chr, pos, ref_aa, alt_aa)
            output_tsv      path of the output table
            [position_scores_tsv]  optional fifth argument, accepted and ignored
Outputs:    output_tsv      Gene, Position, hg38_ref_aa, caas_alt_aas, caas_change,
                            caap_group, scheme_weight, then every PrimateAI-3D column;
                            header only when no CAAS row can be mapped
"""

import sys
import os
import gzip
import glob
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from vep_common import (  # noqa: E402
    anc_der_from_row,
    load_convergence_skip,
    load_map_file,
)

# ── Helpers ───────────────────────────────────────────────────────────────────

# Weight written for each caap_group (only US is processed).
SCHEME_WEIGHTS = {
    "US":  1.0,
}


def write_header_only(primateai_gz, output_tsv):
    """Write the output header (and no rows), report it on stderr and exit with status 0."""
    with gzip.open(primateai_gz, 'rt') as gz_in, open(output_tsv, 'w') as out:
        pai_header = gz_in.readline().rstrip('\n')
        out.write(
            "Gene\tPosition\t"
            "hg38_ref_aa\tcaas_alt_aas\tcaas_change\t"
            "caap_group\tscheme_weight\t"
            + pai_header + "\n"
        )
    print("No CAAS rows available for PrimateAI mapping; wrote header-only output.",
          file=sys.stderr)
    sys.exit(0)


# ── Arguments ─────────────────────────────────────────────────────────────────

if len(sys.argv) not in (5, 6):
    sys.exit(
        "Usage: map_to_primateai.py "
        "<caas_file> <vep_map_dir> <primateai_gz> <output_tsv> [position_scores_tsv]"
    )

caas_file, vep_map_dir, primateai_gz, output_tsv = sys.argv[1:5]
position_scores_tsv = sys.argv[5] if len(sys.argv) == 6 else None
skip_positions = load_convergence_skip(position_scores_tsv)


# ── Step 1: CAAS rows grouped by (Gene, Position) ─────────────────────────────

print("Loading CAAS file ...", file=sys.stderr)

caas_targets = {}  # (Gene, Position) -> list of tag dicts

with open(caas_file) as fh:
    header_line = fh.readline().rstrip('\n')
    if not header_line:
        write_header_only(primateai_gz, output_tsv)

    header = header_line.split('\t')
    col = {name.strip(): idx for idx, name in enumerate(header)}
    col_lc = {name.strip().lower(): idx for idx, name in enumerate(header)}

    tag_col = col_lc.get('tag')
    caas_col = col_lc.get('caas')
    gene_col = col_lc.get('gene')
    pos_col = col_lc.get('position')
    cside_col = col.get("side") or col_lc.get("side")
    anc_col = col_lc.get('ancestral_aa')
    der_col = col_lc.get('derived_aa')
    amino_col = col_lc.get('amino_encoded')
    caap_col = col.get('caap_group') or col_lc.get('caap_group')
    dres_col = col_lc.get('derived_residues')

    if any(c is None for c in [gene_col, pos_col]):
        print("Missing required columns (Gene, Position) in CAAS / position_scores file header.", file=sys.stderr)
        write_header_only(primateai_gz, output_tsv)

    for line in fh:
        fields = line.rstrip('\n').split('\t')
        if len(fields) < len(header):
            continue

        caap_grp = fields[caap_col] if caap_col is not None else 'US'
        if caap_grp != 'US':
            continue

        gene = fields[gene_col]
        try:
            position = int(fields[pos_col])
        except ValueError:
            continue

        if skip_positions is not None and (gene, position) in skip_positions:
            continue  # position excluded by load_convergence_skip()

        tag = fields[tag_col] if tag_col is not None else f"{gene}_{position}"
        caas_pat = fields[caas_col] if caas_col is not None else ''
        cside = fields[cside_col] if cside_col is not None else ''
        anc_val = fields[anc_col] if anc_col is not None else ''
        der_val = fields[der_col] if der_col is not None else ''
        amino_enc = fields[amino_col] if amino_col is not None else ''
        caas_change = amino_enc if amino_enc else caas_pat
        weight = SCHEME_WEIGHTS.get(caap_grp, 1.0)


        anc_aas, der_aas = anc_der_from_row(anc_val, der_val, caas_pat, cside)

        key = (gene, position)
        if key not in caas_targets:
            caas_targets[key] = []
        caas_targets[key].append({
            'tag': tag,
            'anc_aas': anc_aas,
            'der_aas': der_aas,
            'caap_group': caap_grp,
            'weight': weight,
            'caas_change': caas_change,
        })

print(f"  {len(caas_targets)} unique (Gene, Position) targets loaded.", file=sys.stderr)

if not caas_targets:
    write_header_only(primateai_gz, output_tsv)


# ── Step 2: genomic position lookup from the MAP files ────────────────────────

print("Loading MAP files and building genomic lookup ...", file=sys.stderr)

pos_lookup = {}  # (chrom, pos) -> list of (gene, caas_pos, hg38_aa_pos, tag_entries)
loaded_maps = {}  # gene -> (pos_map, strand)

for (gene, caas_pos) in caas_targets:
    if gene not in loaded_maps:
        pos_map, strand = load_map_file(gene, vep_map_dir)
        loaded_maps[gene] = (pos_map, strand)
    else:
        pos_map, strand = loaded_maps[gene]

    if not pos_map:
        continue

    # prot_ali_col of the MAP file is a 1-based alignment column; Position is 0-based.
    prot_ali_idx = caas_pos + 1
    if prot_ali_idx not in pos_map:
        continue

    hg38_aa_pos, chrom, coord = pos_map[prot_ali_idx]

    # The three nucleotides of the codon, in the direction of the strand.
    if strand == '+':
        codon_positions = [coord, coord + 1, coord + 2]
    else:
        codon_positions = [coord, coord - 1, coord - 2]

    tag_entries = caas_targets[(gene, caas_pos)]

    for nt_pos in codon_positions:
        key = (chrom, nt_pos)
        if key not in pos_lookup:
            pos_lookup[key] = []
        pos_lookup[key].append((gene, caas_pos, hg38_aa_pos, tag_entries))

print(f"  {len(pos_lookup)} genomic positions to scan.", file=sys.stderr)

if not pos_lookup:
    write_header_only(primateai_gz, output_tsv)


# ── Step 3: stream PrimateAI-3D and write the matching rows ───────────────────

print("Streaming PrimateAI-3D database ...", file=sys.stderr)

matched = 0
scanned = 0

with gzip.open(primateai_gz, 'rt') as gz_in, open(output_tsv, 'w') as out:
    pai_header = gz_in.readline().rstrip('\n')
    pai_cols = pai_header.split('\t')

    try:
        ref_aa_col = pai_cols.index('ref_aa')
        alt_aa_col = pai_cols.index('alt_aa')
        chr_col = pai_cols.index('chr')
        pos_col = pai_cols.index('pos')
    except ValueError as exc:
        sys.exit(f"Missing expected column in PrimateAI header: {exc}")

    out.write(
        "Gene\tPosition\t"
        "hg38_ref_aa\tcaas_alt_aas\tcaas_change\t"
        "caap_group\tscheme_weight\t"
        + pai_header + "\n"
    )

    for line in gz_in:
        scanned += 1
        if scanned % 5_000_000 == 0:
            print(f"  ... scanned {scanned:,} PrimateAI rows, {matched} matched",
                  file=sys.stderr)

        line = line.rstrip('\n')
        if not line:
            continue
        fields = line.split('\t')
        if len(fields) <= max(ref_aa_col, alt_aa_col):
            continue

        chrom = fields[chr_col]
        try:
            pos = int(fields[pos_col])
        except ValueError:
            continue

        key = (chrom, pos)
        if key not in pos_lookup:
            continue

        ref_aa = fields[ref_aa_col]
        alt_aa = fields[alt_aa_col]

        for gene, caas_pos, hg38_aa_pos, tag_entries in pos_lookup[key]:
            for entry in tag_entries:
                der_aas = entry['der_aas']
                anc_aas = entry['anc_aas']
                caap_grp = entry['caap_group']
                weight = entry['weight']
                caas_change = entry['caas_change']

                # Determine valid alternative amino acids to match in PrimateAI:
                # If human ref_aa is ancestral, target the derived state(s).
                # If human ref_aa is derived, target the ancestral state(s) (reverse transition in hg38).
                # Otherwise, target any CAAS state.
                if anc_aas and ref_aa in anc_aas:
                    alt_aas = der_aas - {ref_aa}
                elif der_aas and ref_aa in der_aas:
                    alt_aas = anc_aas - {ref_aa}
                else:
                    alt_aas = (der_aas | anc_aas) - {ref_aa}

                if not alt_aas:
                    alt_aas = der_aas | anc_aas

                if alt_aa not in alt_aas:
                    continue

                out.write(
                    f"{gene}\t{caas_pos}\t"
                    f"{ref_aa}\t{''.join(sorted(alt_aas))}\t{caas_change}\t"
                    f"{caap_grp}\t{weight}\t"
                    f"{line}\n"
                )
                matched += 1

print(
    f"Done. Scanned {scanned:,} PrimateAI rows → {matched} matched entries written.",
    file=sys.stderr,
)
