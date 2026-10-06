#!/usr/bin/env python3
# map_to_cosmic.py — Map CAAS positions to COSMIC somatic mutations.
# PhyloPhere | subworkflows/VEP/local/src/

"""
MapToCosmic: For every CAAS position (Gene, Position), finds the COSMIC Mutant Census
missense mutations at its hg38 codon whose amino-acid change is the ancestral→derived
change of the CAAS, and writes them with their COSMIC annotation.

Strategy:
  1. Group the CAAS rows by (Gene, Position); only the `US` caap_group is used (rows
     without a caap_group column count as `US`).
  2. Translate each position to its hg38 codon with the per-gene MAP files (alignment
     column = Position + 1) and index the three nucleotides of the codon by
     (chromosome, coordinate); the strand is inferred from the MAP file.
  3. Stream the COSMIC table (gzip TSV), look up each mutation by (chromosome,
     GENOME_START) and keep it when MUTATION_AA is a plain missense change (p.X123Y),
     X is an ancestral residue of the CAAS (when that set is known) and Y a derived one.

Called by:  COSMIC_MAP Nextflow process (cosmic.nf → map_to_cosmic.py)
Inputs:     caas_file     position_scores.tsv, with Gene, Position and, when present, tag,
                          caas, side, amino_encoded, caap_group and the residue tally
                          columns (top_species_residues / bottom_species_residues)
            vep_map_dir   directory of the per-gene MAP files
            cosmic_gz     COSMIC Mutant Census GRCh38 table (gzip TSV)
            output_tsv    path of the output table
            [position_scores_tsv]  optional fifth argument, accepted and ignored
Outputs:    output_tsv    Gene, Position, hg38_ref_aa, caas_alt_aas, caas_change,
                          caap_group, scheme_weight, CHROMOSOME, GENOME_START,
                          GENOMIC_WT_ALLELE, GENOMIC_MUT_ALLELE, MUTATION_AA,
                          MUTATION_DESCRIPTION, MUTATION_SOMATIC_STATUS; header only
                          when no CAAS row can be mapped
"""

# ── Standard library ──────────────────────────────────────────────────────────
import sys
import os
import gzip
import re
import glob
from pathlib import Path

# ── Package-internal ──────────────────────────────────────────────────────────
sys.path.insert(0, str(Path(__file__).parent))
from vep_common import (  # noqa: E402
    _support_letters,
    anc_der_from_descriptor,
    load_convergence_skip,
    load_map_file,
)


# ── Helpers ───────────────────────────────────────────────────────────────────

# Weight written for each caap_group (only US is processed).
SCHEME_WEIGHTS = {
    "US":  1.0,
}


def write_header_only(output_tsv):
    """Write the output header (and no rows) and exit with status 0."""
    with open(output_tsv, 'w') as out:
        out.write(
            "Gene\tPosition\thg38_ref_aa\tcaas_alt_aas\tcaas_change\tcaap_group\tscheme_weight\t"
            "CHROMOSOME\tGENOME_START\tGENOMIC_WT_ALLELE\tGENOMIC_MUT_ALLELE\t"
            "MUTATION_AA\tMUTATION_DESCRIPTION\tMUTATION_SOMATIC_STATUS\n"
        )
    sys.exit(0)


# ── CLI ───────────────────────────────────────────────────────────────────────


def main():
    """Map the CAAS positions of argv[1] to COSMIC mutations (see the module docstring)."""
    if len(sys.argv) not in (5, 6):
        sys.exit("Usage: map_to_cosmic.py <caas_file> <vep_map_dir> <cosmic_gz> "
                 "<output_tsv> [position_scores_tsv]")

    caas_file, vep_map_dir, cosmic_gz, output_tsv = sys.argv[1:5]
    position_scores_tsv = sys.argv[5] if len(sys.argv) == 6 else None
    skip_positions = load_convergence_skip(position_scores_tsv)

    print("Loading CAAS file ...", file=sys.stderr)
    caas_targets = {}
    
    with open(caas_file) as fh:
        header_line = fh.readline().rstrip('\n')
        if not header_line:
            write_header_only(output_tsv)

        header = header_line.split('\t')
        col = {name.strip(): idx for idx, name in enumerate(header)}
        col_lc = {name.strip().lower(): idx for idx, name in enumerate(header)}

        tag_col = col_lc.get('tag')
        caas_col = col_lc.get('caas')
        gene_col = col_lc.get('gene')
        pos_col = col_lc.get('position')
        cside_col = col.get("side") or col_lc.get("side")
        amino_col = col_lc.get('amino_encoded')
        caap_col = col.get('caap_group') or col_lc.get('caap_group')
        dres_col = col_lc.get('derived_residues')
        top_sup_col = col_lc.get('top_species_residues') or col_lc.get('top_residue_support')
        bot_sup_col = col_lc.get('bottom_species_residues') or col_lc.get('bottom_residue_support')

        if any(c is None for c in [gene_col, pos_col]):
            print("Missing required columns (Gene, Position) in CAAS / position_scores file header.", file=sys.stderr)
            write_header_only(output_tsv)

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
            amino_enc = fields[amino_col] if amino_col is not None else ''
            caas_change = amino_enc if amino_enc else caas_pat
            weight = SCHEME_WEIGHTS.get(caap_grp, 1.0)

            top_sup = fields[top_sup_col] if top_sup_col is not None else ''
            bot_sup = fields[bot_sup_col] if bot_sup_col is not None else ''

            anc_aas, der_aas = anc_der_from_descriptor(
                top_sup,
                bot_sup,
                cside,
                caas=caas_pat,
            )

            key = (gene, position)
            if key not in caas_targets:
                box = []
                caas_targets[key] = box
            caas_targets[key].append({
                'tag': tag,
                'anc_aas': anc_aas,
                'der_aas': der_aas,
                'caap_group': caap_grp,
                'weight': weight,
                'caas_change': caas_change,
            })

    if not caas_targets:
        write_header_only(output_tsv)

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

        prot_ali_idx = caas_pos + 1  # Position is 0-based, prot_ali_col 1-based
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
        write_header_only(output_tsv)

    print("Streaming COSMIC database ...", file=sys.stderr)
    matched = 0
    scanned = 0

    with gzip.open(cosmic_gz, 'rt') as gz_in, open(output_tsv, 'w') as out:
        header_line = gz_in.readline().rstrip('\n')
        header_cols = header_line.split('\t')
        
        try:
            chr_col = header_cols.index('CHROMOSOME')
            start_col = header_cols.index('GENOME_START')
            wt_col = header_cols.index('GENOMIC_WT_ALLELE')
            mut_col = header_cols.index('GENOMIC_MUT_ALLELE')
            mut_aa_col = header_cols.index('MUTATION_AA')
            desc_col = header_cols.index('MUTATION_DESCRIPTION')
            status_col = header_cols.index('MUTATION_SOMATIC_STATUS')
        except ValueError as exc:
            sys.exit(f"Missing expected column in COSMIC header: {exc}")

        out.write(
            "Gene\tPosition\thg38_ref_aa\tcaas_alt_aas\tcaas_change\tcaap_group\tscheme_weight\t"
            "CHROMOSOME\tGENOME_START\tGENOMIC_WT_ALLELE\tGENOMIC_MUT_ALLELE\t"
            "MUTATION_AA\tMUTATION_DESCRIPTION\tMUTATION_SOMATIC_STATUS\n"
        )

        for line in gz_in:
            scanned += 1
            if scanned % 2_000_000 == 0:
                print(f"  ... scanned {scanned:,} COSMIC rows, {matched} matched", file=sys.stderr)

            line = line.rstrip('\n')
            if not line:
                continue
            fields = line.split('\t')
            if len(fields) <= max(chr_col, start_col, status_col):
                continue

            chrom = fields[chr_col]
            if not chrom.startswith('chr'):
                chrom = f"chr{chrom}"
            try:
                pos = int(fields[start_col])
            except ValueError:
                continue

            key = (chrom, pos)
            if key not in pos_lookup:
                continue

            wt_allele = fields[wt_col]
            mut_allele = fields[mut_col]
            mut_aa = fields[mut_aa_col]
            mut_desc = fields[desc_col]
            somatic_status = fields[status_col]

            for gene, caas_pos, hg38_aa_pos, tag_entries in pos_lookup[key]:
                for entry in tag_entries:
                    der_aas = entry['der_aas']
                    anc_aas = entry['anc_aas']
                    caap_grp = entry['caap_group']
                    weight = entry['weight']
                    caas_change = entry['caas_change']

                    # Only a plain missense change (p.X123Y) can match the CAAS transition.
                    match_parts = re.match(r'^p\.([A-Z])\d+([A-Z])$', mut_aa)
                    if not match_parts:
                        continue

                    ref_aa = match_parts.group(1)
                    alt_aa = match_parts.group(2)

                    is_fwd = (not anc_aas or ref_aa in anc_aas) and (alt_aa in der_aas)
                    is_rev = (ref_aa in der_aas) and (alt_aa in anc_aas)
                    if not (is_fwd or is_rev):
                        continue

                    out.write(
                        f"{gene}\t{caas_pos}\t"
                        f"{ref_aa}\t{''.join(sorted(der_aas))}\t{caas_change}\t"
                        f"{caap_grp}\t{weight}\t"
                        f"{chrom}\t{pos}\t{wt_allele}\t{mut_allele}\t"
                        f"{mut_aa}\t{mut_desc}\t{somatic_status}\n"
                    )
                    matched += 1

    print(f"Done. Scanned {scanned:,} COSMIC rows → {matched} matched entries written.", file=sys.stderr)

if __name__ == "__main__":
    main()
