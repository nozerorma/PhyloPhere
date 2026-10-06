#!/usr/bin/env python3
# vep_common.py — Ancestral/derived residues and MAP-file parsing shared by the VEP annotation scripts.
# PhyloPhere | subworkflows/VEP/local/src/

"""
VepCommon: Helpers shared by the three annotation sources, so that they read a
CAAS position and its MAP file identically: the ancestral and derived residues of
a position, and the table from alignment column to hg38 codon and strand.

Imported by: map_to_primateai.py, map_to_cosmic.py, build_vep_hgvs.py (copied next
to them in the work directory by PRIMATEAI_MAP, COSMIC_MAP and ENSEMBL_VEP_ANNOTATE)
"""

# ── Standard library ──────────────────────────────────────────────────────────
import os
import sys
import glob


# ── Residue sets ──────────────────────────────────────────────────────────────


def anc_der_from_caas(caas, side):
    """(ancestral_aas, derived_aas) of a CAAS from its pattern and the side it was found on.

    The pattern is "<residues of the top group>/<residues of the bottom group>" (e.g. "KQ/E"). The side sets the
    direction: "top" makes the bottom residues ancestral and the top ones derived, "bottom" the reverse; any other
    side leaves the ancestral set empty and takes both as derived. Residues shared by both sets are removed from the
    derived set unless that would empty it.

    Only the pattern is read. The residue tallies of a position count every species of a group, so they also hold
    residues that are not part of the CAAS (the bottom residue present in half of the top species).
    """
    raw_top, _, raw_bot = str(caas or "").partition("/")
    top_set = {c for c in raw_top.upper() if c.isalpha()}
    bot_set = {c for c in raw_bot.upper() if c.isalpha()}

    cs = str(side or "").strip().lower()
    if cs == "top":
        anc, der = bot_set, top_set          # ancestral = bottom, derived = top
    elif cs == "bottom":
        anc, der = top_set, bot_set          # ancestral = top, derived = bottom
    else:                                    # "none" / unknown: both sides derived
        anc, der = set(), top_set | bot_set

    der = (der - anc) or der
    return anc, der


def load_convergence_skip(position_scores_tsv):
    """Positions to skip, as a set of (gene, position); always None (no position is skipped).

    The position_scores argument is accepted so callers can pass the optional
    command-line argument through; it is not read.
    """
    return None


# ── MAP files ─────────────────────────────────────────────────────────────────


def load_map_file(gene, vep_map_dir):
    """Parse the MAP file of `gene` (`<gene>.*.map.tsv`, else `<gene>.map.tsv`).

    Returns (pos_map, strand), or (None, None) when the file is missing, lacks
    the columns hg38_nt_coord, hg38_aa_pos and prot_ali_col, or has no usable
    coordinate. pos_map is {prot_ali_col (1-based alignment column): (hg38_aa_pos,
    chromosome, hg38 coordinate of the first codon nucleotide)}; strand is "+" or
    "-", inferred from the trend of the coordinates down the file.
    """
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

    # Strand: the sign of the summed steps between consecutive coordinates.
    is_minus = False
    if len(coords) > 1:
        diffs = [coords[i+1] - coords[i] for i in range(len(coords)-1)]
        sign_sum = sum(1 if d > 0 else -1 for d in diffs if d != 0)
        is_minus = (sign_sum < 0)
    strand = '-' if is_minus else '+'

    # Rows without a valid alignment column or coordinate are left out of the map.
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
