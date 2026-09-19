#!/usr/bin/env python3
"""Map CAAS protein positions to COSMIC somatic mutations using upstream MAP files.

Strategy:
1. Load CAAS file and group targets by (Gene, Position).
2. Load MAP files and build genomic position lookup (chrom, nt_pos) -> (gene, caas_pos).
3. Stream Cosmic_MutantCensus_v104_GRCh38.tsv.gz, match coordinates, and output matching rows.
"""

import sys
import os
import gzip
import re
import glob
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from vep_common import (  # noqa: E402
    _support_letters,
    anc_der_from_descriptor,
    load_convergence_skip,
    load_map_file,
)

SCHEME_WEIGHTS = {
    "US":  1.0,
}


def write_header_only(output_tsv):
    with open(output_tsv, 'w') as out:
        out.write(
            "Gene\tPosition\thg38_ref_aa\tcaas_alt_aas\tcaas_change\tcaap_group\tscheme_weight\t"
            "CHROMOSOME\tGENOME_START\tGENOMIC_WT_ALLELE\tGENOMIC_MUT_ALLELE\t"
            "MUTATION_AA\tMUTATION_DESCRIPTION\tMUTATION_SOMATIC_STATUS\n"
        )
    sys.exit(0)

def main():
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
        top_sup_col = col_lc.get('top_residue_support')
        bot_sup_col = col_lc.get('bottom_residue_support')

        if any(c is None for c in [tag_col, caas_col, gene_col, pos_col]):
            print("Missing required columns in CAAS file header.", file=sys.stderr)
            write_header_only(output_tsv)

        if dres_col is None:
            print(
                "WARN: `derived_residues` column absent from CAAS file — regenerate "
                "filtered_discovery.tsv through CT_POSTPROC. Falling back to raw "
                "`caas` pattern letters for this run.",
                file=sys.stderr,
            )

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
                continue  # fractional FOP rule: derived residues genuinely disagree

            tag = fields[tag_col]
            caas_pat = fields[caas_col]
            cside = fields[cside_col] if cside_col is not None else ''
            amino_enc = fields[amino_col] if amino_col is not None else ''
            caas_change = amino_enc
            weight = SCHEME_WEIGHTS.get(caap_grp, 1.0)

            if dres_col is not None:
                anc_aas, der_aas = anc_der_from_descriptor(
                    fields[dres_col],
                    fields[top_sup_col] if top_sup_col is not None else '',
                    fields[bot_sup_col] if bot_sup_col is not None else '',
                    cside,
                )
            else:
                raw_top, _, raw_bot = caas_pat.partition('/')
                anc_aas = {c for c in raw_bot.upper() if c.isalpha()} if cside == 'top' else \
                          ({c for c in raw_top.upper() if c.isalpha()} if cside == 'bottom' else set())
                der_aas = {c for c in raw_top.upper() if c.isalpha()} if cside == 'top' else \
                          ({c for c in raw_bot.upper() if c.isalpha()} if cside == 'bottom' else
                           {c for c in (raw_top + raw_bot).upper() if c.isalpha()})
                der_aas = (der_aas - anc_aas) or der_aas  # drop unchanged residues

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
    pos_lookup = {}
    loaded_maps = {}

    for (gene, caas_pos) in caas_targets:
        if gene not in loaded_maps:
            pos_map, strand = load_map_file(gene, vep_map_dir)
            loaded_maps[gene] = (pos_map, strand)
        else:
            pos_map, strand = loaded_maps[gene]

        if not pos_map:
            continue

        prot_ali_idx = caas_pos + 1
        if prot_ali_idx not in pos_map:
            continue

        hg38_aa_pos, chrom, coord = pos_map[prot_ali_idx]

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

        # Output header
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

                    # Parse amino acid change if possible
                    # Somatic mutation must be a clean missense mutation matching the CAAS transition
                    match_parts = re.match(r'^p\.([A-Z])\d+([A-Z])$', mut_aa)
                    if not match_parts:
                        continue

                    ref_aa = match_parts.group(1)
                    alt_aa = match_parts.group(2)

                    if anc_aas and ref_aa not in anc_aas:
                        continue
                    if alt_aa not in der_aas:
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
