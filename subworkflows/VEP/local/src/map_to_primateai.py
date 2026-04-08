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


AMBIGSYMS = {'-', 'X', 'B', 'Z', 'J', 'U', 'O'}

SCHEME_WEIGHTS = {
    "US":  0.5,
    "GS4": 0.2,
    "GS3": 0.1,
    "GS2": 0.1,
    "GS1": 0.1,
    "GS0": 0.0,
}


def pair_aware_caas_letters_gs(caas_pattern, amino_encoded, change_side, caap_group):
    """Return (ancestral_aas, derived_aas) as raw AA sets, with GS-aware conserved-pair exclusion.

    When amino_encoded is available (non-US grouped string from CT_DISAMBIGUATION),
    uses it positionally to detect positions that are conserved *at the group level*,
    then collects raw letters from caas at non-conserved positions — matching the
    PGLS pair-aware logic while keeping raw AAs for PrimateAI lookup.

    When amino_encoded is absent or the scheme is US, falls back to raw character
    comparison (pair-aware at the amino-acid level).

    Parameters
    ----------
    caas_pattern  : str  e.g. 'LAP/PPP', 'TTT/TNI'
    amino_encoded : str  e.g. 'aaa/ass' (group-label string from disambiguation), or ''
    change_side   : str  'top', 'bottom', or 'none'/'ambiguous'/''
    caap_group    : str  e.g. 'US', 'GS2'

    Returns
    -------
    (ancestral_aas, derived_aas) : tuple of set of str  (raw 1-letter AA codes)
        ancestral_aas – raw AAs on the non-changing side at non-conserved positions.
                        Empty when change_side is ambiguous (filter will be bypassed).
        derived_aas   – raw AAs on the changing side at non-conserved positions.
    """
    if '/' not in caas_pattern:
        return set(), _uc_letters(caas_pattern)

    raw_top, raw_bot = caas_pattern.split('/', 1)
    top_nc: set = set()
    bot_nc: set = set()

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


def write_empty_output(primateai_gz, output_tsv):
    """Emit a header-only output file and exit successfully."""
    with gzip.open(primateai_gz, 'rt') as gz_in, open(output_tsv, 'w') as out:
        pai_header = gz_in.readline().rstrip('\n')
        out.write(
            "transvar_query\ttransvar_transcript\t"
            "hg38_ref_aa\tcaas_alt_aas\t"
            + pai_header + "\n"
        )
    print("No CAAS rows available for PrimateAI mapping; wrote header-only output.",
          file=sys.stderr)
    sys.exit(0)


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
# Step 1: Load CAAS file → tag → (caas_pattern, amino_encoded, change_side, caap_group)
# ---------------------------------------------------------------------------
print("Loading CAAS file ...", file=sys.stderr)

tag_to_caas      = {}   # tag → raw AA pattern string
tag_to_amino_enc = {}   # tag → group-label pattern string (from CT_DISAMBIGUATION)
tag_to_cside     = {}   # tag → change_side string
tag_to_caap      = {}   # tag → CAAP_Group string

with open(caas_file) as fh:
    header_line = fh.readline().rstrip('\n')
    if not header_line:
        write_empty_output(primateai_gz, output_tsv)

    header = header_line.split('\t')
    col    = {name.strip(): idx for idx, name in enumerate(header)}
    col_lc = {name.strip().lower(): idx for idx, name in enumerate(header)}

    tag_col      = col_lc.get('tag')
    caas_col     = col_lc.get('caas')
    cside_col    = col.get('change_side')
    amino_col    = col_lc.get('amino_encoded')
    caap_col     = col.get('CAAP_Group') or col_lc.get('caap_group')

    if tag_col is None or caas_col is None:
        write_empty_output(primateai_gz, output_tsv)

    for line in fh:
        fields = line.rstrip('\n').split('\t')
        tag = fields[tag_col]
        if tag in tag_to_caas:
            continue
        tag_to_caas[tag]      = fields[caas_col]
        tag_to_cside[tag]     = fields[cside_col]     if cside_col  is not None else ''
        tag_to_amino_enc[tag] = fields[amino_col]     if amino_col  is not None else ''
        tag_to_caap[tag]      = fields[caap_col]      if caap_col   is not None else 'US'

print(f"  {len(tag_to_caas)} unique tags loaded.", file=sys.stderr)

# ---------------------------------------------------------------------------
# Step 2: Load AA2prot output → query → per-tag entries (one per CAAP_Group scheme)
#
# A single GENE:p.N query maps to multiple tags (one per CAAP_Group scheme).
# We keep them separate so each can carry its own GS-aware ancestral/derived
# letter sets and scheme weight, rather than unioning across schemes.
# ---------------------------------------------------------------------------
print("Loading AA2prot output ...", file=sys.stderr)

# query_tag_entries[GENE:p.N] = list of dicts:
#   { caap_group, weight, anc_aas (raw), der_aas (raw) }
query_tag_entries: dict[str, list] = {}

with open(aa2prot_csv) as fh:
    for line in fh:
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 4:
            continue
        tag   = fields[2]
        query = fields[3]   # GENE:p.N

        caas_pat  = tag_to_caas.get(tag, '')
        cside     = tag_to_cside.get(tag, '')
        amino_enc = tag_to_amino_enc.get(tag, '')
        caap_grp  = tag_to_caap.get(tag, 'US')
        if not caas_pat:
            continue

        weight = SCHEME_WEIGHTS.get(caap_grp, 0.0)
        anc_aas, der_aas = pair_aware_caas_letters_gs(
            caas_pat, amino_enc, cside, caap_grp
        )

        if query not in query_tag_entries:
            query_tag_entries[query] = []
        query_tag_entries[query].append({
            'caap_group': caap_grp,
            'weight':     weight,
            'anc_aas':   anc_aas,
            'der_aas':   der_aas,
        })

print(f"  {len(query_tag_entries)} unique protein-position queries.", file=sys.stderr)

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
# Each query entry carries the list of per-tag dicts (one per CAAP_Group),
# each with its own GS-aware ancestral/derived sets and scheme weight.
# alt_aas = der_aas − {ref_aa_hg38} is computed here once per tag.
# ---------------------------------------------------------------------------
print("Building position lookup ...", file=sys.stderr)

# pos_lookup[(chrom, pos)] = list of (query, ref_aa_hg38, transcript, tag_entries)
# tag_entries = list of { caap_group, weight, anc_aas, alt_aas }
pos_lookup: dict[tuple, list] = {}

for query, (chrom, start, end, transcript, ref_aa_hg38) in query_coords.items():
    raw_tag_entries = query_tag_entries.get(query, [])
    if not raw_tag_entries:
        continue

    # Pre-compute alt_aas (ref excluded) per tag entry
    tag_entries = []
    for entry in raw_tag_entries:
        alt_aas = entry['der_aas'] - {ref_aa_hg38}
        if not alt_aas:
            alt_aas = entry['der_aas']   # safety: keep all if ref subtraction empties the set
        tag_entries.append({
            'caap_group': entry['caap_group'],
            'weight':     entry['weight'],
            'anc_aas':    entry['anc_aas'],
            'alt_aas':    alt_aas,
        })

    for pos in range(start, end + 1):
        key = (chrom, pos)
        if key not in pos_lookup:
            pos_lookup[key] = []
        pos_lookup[key].append((query, ref_aa_hg38, transcript, tag_entries))

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

    # Output header — two extra columns: caap_group, scheme_weight
    out.write(
        "transvar_query\ttransvar_transcript\t"
        "hg38_ref_aa\tcaas_alt_aas\t"
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

        for query, ref_aa_hg38, transcript, tag_entries in pos_lookup[key]:
            if ref_aa != ref_aa_hg38:
                continue   # PrimateAI ref doesn't match hg38 reference at this position

            for entry in tag_entries:
                alt_aas  = entry['alt_aas']
                anc_aas  = entry['anc_aas']
                caap_grp = entry['caap_group']
                weight   = entry['weight']

                if alt_aa not in alt_aas:
                    continue   # not a CAAS-derived substitution for this scheme
                # Accept only when human ref is ancestral (non-changing side).
                # When ancestral is unknown (change_side=none), anc_aas is empty
                # and the filter is bypassed to preserve existing behaviour.
                if anc_aas and ref_aa_hg38 not in anc_aas:
                    continue   # human ref is not ancestral under this scheme — skip
                out.write(
                    f"{query}\t{transcript}\t"
                    f"{ref_aa_hg38}\t{''.join(sorted(alt_aas))}\t"
                    f"{caap_grp}\t{weight}\t"
                    f"{line}\n"
                )
                matched += 1

print(
    f"Done. Scanned {scanned:,} PrimateAI rows → {matched} matched entries written.",
    file=sys.stderr,
)
