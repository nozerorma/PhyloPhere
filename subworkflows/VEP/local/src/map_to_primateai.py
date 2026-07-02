#!/usr/bin/python3
"""Map CAAS protein positions to PrimateAI-3D scores using upstream MAP files.

Strategy
--------
For each CAAS position (Gene + Position):

1. **MAP file lookup** — read the gene's MAP file from the upstream directory.
   Match the 0-based CAAS position to `prot_ali_col - 1`.
   Extract the 1-based protein position (`hg38_aa_pos`) and genomic coordinate (`hg38_nt_coord`).

2. **Strand inference** — compare adjacent valid codon coordinates in the MAP file.
   - If coordinates increase: plus strand (+). Nucleotides are C, C+1, C+2.
   - If coordinates decrease: minus strand (-). Nucleotides are C, C-1, C-2.

3. **PrimateAI Lookup**: stream through the PrimateAI gz database and keep variants at
   any of the 3 nucleotide positions where the alternative amino acid matches the CAAS pattern.

Directional Pathogenicity Interpretation & Filtering
----------------------------------------------------
The script enforces a directional constraint to align PrimateAI lookups with evolutionary medicine context:
  `if anc_aas and ref_aa not in anc_aas: continue`
- **TOP Direction (Ancestral Human / Derived in Top)**: We look at the pathogenicity of mutating the ancestral human allele into the derived, "high-cancer" susceptibility allele.
- **BOTTOM Direction (Ancestral Human / Derived in Bottom)**: We look at the tolerance or pathogenicity of mutating the ancestral human allele (associated with the "high-cancer" susceptibility side) into the derived, "low-cancer" resistance allele.
- **Filtering**: If the human reference allele is not in the ancestral set, the position is filtered out because the human reference background does not reflect the ancestral state under the scheme. Similarly, if the derived state matches the human reference, it is skipped as no missense variant can represent a mutation to it.

Usage
-----
    map_to_primateai.py <caas_file> <vep_map_dir> <primateai_gz> <output_tsv>
"""

import sys
import os
import gzip
import re
import glob

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _uc_letters(text):
    """Return the set of uppercase A–Z letters (1-letter AA codes) in *text*."""
    return set(re.findall(r'[A-Z]', text))


AMBIGSYMS = {'-', 'X', 'B', 'Z', 'J', 'U', 'O'}

SCHEME_WEIGHTS = {
    "US":  1.0,
}


def pair_aware_caas_letters_gs(caas_pattern, amino_encoded, change_side, caap_group):
    """Return (ancestral_aas, derived_aas) as raw AA sets, with GS-aware conserved-pair exclusion.

    Uses group-aware logic while keeping raw AAs for PrimateAI lookup.
    """
    if '/' not in caas_pattern:
        return set(), _uc_letters(caas_pattern)

    raw_top, raw_bot = caas_pattern.split('/', 1)
    top_nc = set()
    bot_nc = set()

    use_enc = (caap_group != 'US'
               and amino_encoded
               and '/' in str(amino_encoded))

    if use_enc:
        enc_top, enc_bot = str(amino_encoded).split('/', 1)
        # Conserved at group level when enc_top[i] == enc_bot[i]; collect raw letters otherwise
        for enc_la, enc_ra, raw_la, raw_ra in zip(enc_top, enc_bot, raw_top, raw_bot):
            if enc_la == enc_ra:
                continue   # conserved at group level — skip
            if raw_la and raw_la.upper() not in AMBIGSYMS:
                top_nc.add(raw_la.upper())
            if raw_ra and raw_ra.upper() not in AMBIGSYMS:
                bot_nc.add(raw_ra.upper())
    else:
        # US or no amino_encoded: compare raw letters directly
        for la, ra in zip(raw_top, raw_bot):
            if la == ra:
                continue   # conserved at AA level — skip
            if la and la.upper() not in AMBIGSYMS:
                top_nc.add(la.upper())
            if ra and ra.upper() not in AMBIGSYMS:
                bot_nc.add(ra.upper())

    if change_side == 'top':
        return bot_nc, top_nc          # ancestral=bottom, derived=top
    elif change_side == 'bottom':
        return top_nc, bot_nc          # ancestral=top, derived=bottom
    else:
        # Ambiguous / none: no defined ancestral; return all non-conserved letters
        return set(), top_nc | bot_nc


def load_map_file(gene, vep_map_dir):
    """Parse Gene.*.map.tsv to retrieve codon position-to-coordinate mapping and strand."""
    pattern = os.path.join(vep_map_dir, f"{gene}.*.map.tsv")
    matching = glob.glob(pattern)
    if not matching:
        pattern = os.path.join(vep_map_dir, f"{gene}.map.tsv")
        if os.path.exists(pattern):
            matching = [pattern]
    if not matching:
        print(f"WARN: MAP file for gene {gene} not found in {vep_map_dir}", file=sys.stderr)
        return None, None

    map_file = matching[0]
    coords = []
    rows = []

    with open(map_file, 'r') as fh:
        header_line = fh.readline().rstrip('\n')
        if not header_line:
            return None, None
        header = header_line.split('\t')
        col = {name.strip(): idx for idx, name in enumerate(header)}

        # Verify required columns exist
        required = ['hg38_nt_coord', 'hg38_aa_pos', 'prot_ali_col']
        if not all(c in col for c in required):
            print(f"WARN: MAP file {map_file} is missing required columns", file=sys.stderr)
            return None, None

        for line in fh:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < len(header):
                continue

            nt_coord_val = fields[col['hg38_nt_coord']]
            aa_pos_val = fields[col['hg38_aa_pos']]
            prot_ali_val = fields[col['prot_ali_col']]

            if nt_coord_val != 'NA' and ':' in nt_coord_val:
                try:
                    pos = int(nt_coord_val.split(':')[1])
                    coords.append(pos)
                except ValueError:
                    pass

            rows.append({
                'prot_ali_col': prot_ali_val,
                'hg38_aa_pos': aa_pos_val,
                'hg38_nt_coord': nt_coord_val
            })

    if not coords:
        return None, None

    # Determine strand by coordinate trend
    is_minus = False
    if len(coords) > 1:
        diffs = [coords[i+1] - coords[i] for i in range(len(coords)-1)]
        sign_sum = sum(1 if d > 0 else -1 for d in diffs if d != 0)
        is_minus = (sign_sum < 0)
    strand = '-' if is_minus else '+'

    # Build mapping from 1-based prot_ali_col to (hg38_aa_pos, chrom, coord)
    pos_map = {}
    for r in rows:
        prot_ali = r['prot_ali_col']
        if prot_ali == 'NA' or not prot_ali:
            continue
        try:
            prot_ali_idx = int(prot_ali)
        except ValueError:
            continue

        nt_coord = r['hg38_nt_coord']
        aa_pos_str = r['hg38_aa_pos']
        if nt_coord == 'NA' or aa_pos_str == 'NA' or ':' not in nt_coord:
            continue

        chrom, pos_str = nt_coord.split(':')
        try:
            coord = int(pos_str)
            aa_pos = int(aa_pos_str)
            pos_map[prot_ali_idx] = (aa_pos, chrom, coord)
        except ValueError:
            continue

    return pos_map, strand


def write_header_only(primateai_gz, output_tsv):
    """Emit a header-only output file and exit successfully."""
    with gzip.open(primateai_gz, 'rt') as gz_in, open(output_tsv, 'w') as out:
        pai_header = gz_in.readline().rstrip('\n')
        out.write(
            "Gene\tPosition\t"
            "hg38_ref_aa\tcaas_alt_aas\t"
            "caap_group\tscheme_weight\t"
            + pai_header + "\n"
        )
    print("No CAAS rows available for PrimateAI mapping; wrote header-only output.",
          file=sys.stderr)
    sys.exit(0)

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

if len(sys.argv) != 5:
    sys.exit(
        "Usage: map_to_primateai.py "
        "<caas_file> <vep_map_dir> <primateai_gz> <output_tsv>"
    )

caas_file, vep_map_dir, primateai_gz, output_tsv = sys.argv[1:]

# ---------------------------------------------------------------------------
# Step 1: Load CAAS file → group targets by (Gene, Position)
# ---------------------------------------------------------------------------
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
    cside_col = col.get('change_side')
    amino_col = col_lc.get('amino_encoded')
    caap_col = col.get('caap_group') or col_lc.get('caap_group')

    if any(c is None for c in [tag_col, caas_col, gene_col, pos_col]):
        print("Missing required columns in CAAS file header.", file=sys.stderr)
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

        tag = fields[tag_col]
        caas_pat = fields[caas_col]
        cside = fields[cside_col] if cside_col is not None else ''
        amino_enc = fields[amino_col] if amino_col is not None else ''
        caas_change = amino_enc
        weight = SCHEME_WEIGHTS.get(caap_grp, 1.0)

        anc_aas, der_aas = pair_aware_caas_letters_gs(
            caas_pat, amino_enc, cside, caap_grp
        )

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

# ---------------------------------------------------------------------------
# Step 2: Load MAP files and build genomic position lookup
# ---------------------------------------------------------------------------
print("Loading MAP files and building genomic lookup ...", file=sys.stderr)

pos_lookup = {}  # (chrom, pos) -> list of (gene, caas_pos, hg38_aa_pos, tag_entries)
loaded_maps = {}  # gene -> pos_map

for (gene, caas_pos) in caas_targets:
    if gene not in loaded_maps:
        pos_map, strand = load_map_file(gene, vep_map_dir)
        loaded_maps[gene] = (pos_map, strand)
    else:
        pos_map, strand = loaded_maps[gene]

    if not pos_map:
        continue

    # prot_ali_col in MAP file is 1-based alignment column.
    # CAAS Position is 0-based.
    prot_ali_idx = caas_pos + 1
    if prot_ali_idx not in pos_map:
        continue

    hg38_aa_pos, chrom, coord = pos_map[prot_ali_idx]

    # Generate the 3 nucleotide coordinates for this codon
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

# ---------------------------------------------------------------------------
# Step 3: Stream PrimateAI-3D and emit matching rows
# ---------------------------------------------------------------------------
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

    # Output header
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

                alt_aas = der_aas - {ref_aa}
                if not alt_aas:
                    alt_aas = der_aas   # safety: keep all if ref subtraction empties the set

                if alt_aa not in alt_aas:
                    continue

                # Enforce ancestral constraints if ancestral state is known
                if anc_aas and ref_aa not in anc_aas:
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
