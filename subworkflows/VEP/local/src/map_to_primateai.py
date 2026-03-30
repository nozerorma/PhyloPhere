#!/usr/bin/python3
"""Map CAAS protein positions to PrimateAI-3D scores.

Strategy
--------
For each protein position query (GENE:p.N):

1. **hg38 reference amino acid** — read from the TransVar output.
   TransVar encodes this in the protein coordinates field as "p.NX" where
   X is the single-letter reference AA (e.g. "p.26E" → ref = E).

2. **Derived (changing) amino acids** — determined from the CAAS pattern
   and ``change_side``:
   * change_side = top    → the top group is the changing one
   * change_side = bottom → the bottom group is the changing one
   * change_side = none / ambiguous → use all unique AAs from both groups

   ``alt_aas = changing_group_aas − {ref_aa_hg38}``

3. **PrimateAI lookup**: stream through the gz file and keep variants where
   ``ref_aa == ref_aa_hg38  AND  alt_aa ∈ alt_aas``

   * Convergent CAAS (one changing AA)  → single PrimateAI entry per codon pos.
   * Divergent CAAS  (multiple alt AAs) → multiple entries (one per alt AA).
   * Ambiguous CAAS  (change_side=none) → all missense variants at those positions.

Usage
-----
    map_to_primateai.py <transvar_tsv> <aa2prot_csv> <caas_file> <primateai_gz> <output_tsv>

Arguments
---------
transvar_tsv   Stage 3 output: transvar panno --ccds --longest --noheader
               Columns: input | transcript | gene | strand |
                         coordinates(gDNA/cDNA/protein) | region | info
aa2prot_csv    Stage 2 output (AA2prot global.csv)
               Columns: gene | caas_pos | tag | GENE:p.N
caas_file      Original CAAS discovery TSV
               Columns: Gene | Position | tag | caas | ... | change_side | ...
primateai_gz   PrimateAI-3D.hg38.txt.gz
output_tsv     Output path for the merged results table
"""

import sys
import gzip
import re

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _uc_letters(text):
    """Return the set of uppercase A–Z letters (1-letter AA codes) in *text*."""
    return set(re.findall(r'[A-Z]', text))


def changing_side_aas(caas_pattern, change_side):
    """Return the set of amino acids on the *changing* side of the CAAS pattern.

    Parameters
    ----------
    caas_pattern : str  e.g. 'DED/EEE', 'TTT/TIA', 'QSR/RRR'
    change_side  : str  'top', 'bottom', or 'none'/'ambiguous'/''

    Returns
    -------
    set of str  – 1-letter AA codes observed on the changing side
                  (or union of both sides if direction is ambiguous)
    """
    if '/' not in caas_pattern:
        return _uc_letters(caas_pattern)

    top_str, bot_str = caas_pattern.split('/', 1)
    top_aas = _uc_letters(top_str)
    bot_aas = _uc_letters(bot_str)

    if change_side == 'top':
        return top_aas
    elif change_side == 'bottom':
        return bot_aas
    else:
        # Ambiguous / none: include all unique AAs from both groups
        return top_aas | bot_aas


def parse_codon_range(coord_str):
    """Parse TransVar coordinate string → (chrom, start, end).

    Handles:
      'chr17:g.7675087_7675089/c.523_525/p.175R'  → range
      'chr17:g.7675087/c.523/p.175R'              → single position
    Returns (None, None, None) if not parseable.
    """
    m = re.match(r'(chr[\w]+):g\.(\d+)_(\d+)(?:[/\s]|$)', coord_str)
    if m:
        return m.group(1), int(m.group(2)), int(m.group(3))
    m = re.match(r'(chr[\w]+):g\.(\d+)(?:[/\s]|$)', coord_str)
    if m:
        pos = int(m.group(2))
        return m.group(1), pos, pos
    return None, None, None


def parse_ref_aa_hg38(coord_str):
    """Extract the hg38 reference amino acid from a TransVar protein field.

    TransVar encodes it as '/p.NX' where N is the position and X is the
    single-letter AA (e.g. '/p.26E' → 'E').  Returns '' if not found.
    """
    m = re.search(r'/p\.\d+([A-Z])', coord_str)
    return m.group(1) if m else ''


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

if len(sys.argv) != 6:
    sys.exit(
        "Usage: map_to_primateai.py "
        "<transvar_tsv> <aa2prot_csv> <caas_file> <primateai_gz> <output_tsv>"
    )

transvar_tsv, aa2prot_csv, caas_file, primateai_gz, output_tsv = sys.argv[1:]

# ---------------------------------------------------------------------------
# Step 1: Load CAAS file → tag → (caas_pattern, change_side)
# ---------------------------------------------------------------------------
print("Loading CAAS file ...", file=sys.stderr)

tag_to_caas  = {}   # tag → uppercase AA pattern string
tag_to_cside = {}   # tag → change_side string

with open(caas_file) as fh:
    header = fh.readline().rstrip().split('\t')
    col = {name: idx for idx, name in enumerate(header)}

    tag_col   = col['tag']
    caas_col  = col['caas']
    cside_col = col.get('change_side')

    for line in fh:
        fields = line.rstrip('\n').split('\t')
        tag = fields[tag_col]
        if tag in tag_to_caas:
            continue
        tag_to_caas[tag]  = fields[caas_col]
        tag_to_cside[tag] = (fields[cside_col] if cside_col is not None else '')

print(f"  {len(tag_to_caas)} unique tags loaded.", file=sys.stderr)

# ---------------------------------------------------------------------------
# Step 2: Load AA2prot output → query → changing-side AA set
#
# A single GENE:p.N query can appear with multiple tags (multiple CAAS per
# position).  We union the changing-side AA sets across all associated tags.
# ---------------------------------------------------------------------------
print("Loading AA2prot output ...", file=sys.stderr)

query_changing_aas: dict[str, set] = {}   # GENE:p.N → union of changing-side AAs

with open(aa2prot_csv) as fh:
    for line in fh:
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 4:
            continue
        tag   = fields[2]
        query = fields[3]   # GENE:p.N

        caas_pat = tag_to_caas.get(tag, '')
        cside    = tag_to_cside.get(tag, '')
        if not caas_pat:
            continue

        ch_aas = changing_side_aas(caas_pat, cside)

        if query not in query_changing_aas:
            query_changing_aas[query] = set()
        query_changing_aas[query].update(ch_aas)

print(f"  {len(query_changing_aas)} unique protein-position queries.", file=sys.stderr)

# ---------------------------------------------------------------------------
# Step 3: Load TransVar output → query → (chrom, start, end, transcript, ref_aa_hg38)
# ---------------------------------------------------------------------------
print("Loading TransVar output ...", file=sys.stderr)

query_coords: dict[str, tuple] = {}   # query → (chrom, start, end, transcript, ref_aa_hg38)

with open(transvar_tsv) as fh:
    for line in fh:
        line = line.rstrip('\n')
        if not line or line.startswith('input\t'):
            continue
        fields = line.split('\t')
        if len(fields) < 5:
            continue
        query      = fields[0]
        transcript = fields[1]
        coord_str  = fields[4]

        chrom, start, end = parse_codon_range(coord_str)
        ref_aa_hg38       = parse_ref_aa_hg38(coord_str)

        if chrom and query not in query_coords:
            query_coords[query] = (chrom, start, end, transcript, ref_aa_hg38)

print(f"  {len(query_coords)} queries with genomic coordinates.", file=sys.stderr)

# ---------------------------------------------------------------------------
# Step 4: Build position → query lookup for PrimateAI streaming
#
# For each query we now have:
#   ref_aa_hg38   – ground-truth hg38 reference AA (from TransVar)
#   changing_aas  – AA set on the changing side (from CAAS + change_side)
#   alt_aas       – changing_aas − {ref_aa_hg38}  (what we look up in PAI3D)
# ---------------------------------------------------------------------------
print("Building position lookup ...", file=sys.stderr)

# pos_lookup[(chrom, pos)] = list of (query, ref_aa_hg38, alt_aas, transcript)
pos_lookup: dict[tuple, list] = {}

for query, (chrom, start, end, transcript, ref_aa_hg38) in query_coords.items():
    ch_aas  = query_changing_aas.get(query, set())
    alt_aas = ch_aas - {ref_aa_hg38}   # exclude the reference itself
    if not alt_aas:
        alt_aas = ch_aas               # safety: if empty (e.g. all same AA), keep all

    for pos in range(start, end + 1):
        key = (chrom, pos)
        if key not in pos_lookup:
            pos_lookup[key] = []
        pos_lookup[key].append((query, ref_aa_hg38, alt_aas, transcript))

print(f"  {len(pos_lookup)} genomic positions to scan.", file=sys.stderr)

# ---------------------------------------------------------------------------
# Step 5: Stream PrimateAI-3D and emit matching rows
# ---------------------------------------------------------------------------
print("Streaming PrimateAI-3D database ...", file=sys.stderr)

matched = 0
scanned = 0

with gzip.open(primateai_gz, 'rt') as gz_in, open(output_tsv, 'w') as out:
    pai_header = gz_in.readline().rstrip('\n')
    pai_cols   = pai_header.split('\t')

    try:
        ref_aa_col = pai_cols.index('ref_aa')
        alt_aa_col = pai_cols.index('alt_aa')
        chr_col    = pai_cols.index('chr')
        pos_col    = pai_cols.index('pos')
    except ValueError as exc:
        sys.exit(f"Missing expected column in PrimateAI header: {exc}")

    # Output header
    out.write(
        "transvar_query\ttransvar_transcript\t"
        "hg38_ref_aa\tcaas_alt_aas\t"
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

        for query, ref_aa_hg38, alt_aas, transcript in pos_lookup[key]:
            # Accept when:
            #   • PrimateAI ref matches the hg38 reference AA (from TransVar)
            #   • PrimateAI alt is one of the CAAS-derived changing amino acids
            #   • ref ≠ alt (missense; PrimateAI already guarantees this)
            if ref_aa == ref_aa_hg38 and alt_aa in alt_aas:
                out.write(
                    f"{query}\t{transcript}\t"
                    f"{ref_aa_hg38}\t{''.join(sorted(alt_aas))}\t"
                    f"{line}\n"
                )
                matched += 1

print(
    f"Done. Scanned {scanned:,} PrimateAI rows → {matched} matched entries written.",
    file=sys.stderr,
)
