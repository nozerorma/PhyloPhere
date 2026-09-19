#!/usr/bin/env python3
"""
vep_common.py  —  Shared helpers for the VEP annotation-source scripts
(map_to_primateai.py, map_to_cosmic.py, build_vep_hgvs.py).

Extracted verbatim from map_to_primateai.py (map_to_cosmic.py carried an
identical copy) so a third annotation source (Ensembl VEP) can reuse the same
CAAS-descriptor and MAP-file parsing without a fourth copy-paste, and without
importing a script that runs its own database-streaming logic at module
scope.
"""

import os
import sys
import glob


def _support_letters(support_str):
    """Raw AA letters from a `{top,bottom}_residue_support` cell.

    "L:3,S:2" -> {'L', 'S'}. Empty / malformed -> set().
    """
    out = set()
    for tok in str(support_str or "").split(","):
        tok = tok.strip()
        if not tok or ":" not in tok:
            continue
        aa = tok.split(":", 1)[0].strip().upper()
        if aa:
            out.add(aa)
    return out


def anc_der_from_descriptor(derived_residues, top_residue_support,
                            bottom_residue_support, side):
    """(ancestral_aas, derived_aas) from the upstream position-level descriptor.

    The descriptor (built in CT_POSTPROC's residue_descriptors.py) uses the
    ``caas`` left/right convention: ``top_residue_support`` letters are the top
    clade's residues at the position, ``bottom_residue_support`` the bottom
    clade's. ``side`` says which clade carries the derived change; the
    other clade's residues are the ancestral state. Only US rows reach this
    (caap_group filter below), so no GS-grouping handling is needed.

    `derived_residues` itself is not parsed — the support columns already carry
    the per-side residue letters unambiguously.

    Conservation logic: a residue present on BOTH clades did not change, so it is
    dropped from the derived set (``der - anc``). The support columns should
    already exclude it — residue_descriptors.py only fills ``mrca_<i>_<side>_aa``
    for a genuine substitution — but this keeps the filter honest if one leaks.
    A keep-all safety applies if the subtraction empties the set.
    """
    top_set = _support_letters(top_residue_support)
    bot_set = _support_letters(bottom_residue_support)
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
    """RETIRED in scoring_v2 core v3 (V3-4). The old gate skipped positions whose
    ``convergence_schemes`` was "" (fractional-FOP disagreement). core v3 dropped
    ``convergence_schemes`` entirely — a domain that does not converge just scores
    0 — so there is nothing to gate on. Always a no-op; the optional CLI arg is
    accepted but ignored for backward compatibility.
    """
    return None


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
