#                      _              _
#                     | |            | |
#   ___ __ _  __ _ ___| |_ ___   ___ | |___
#  / __/ _` |/ _` / __| __/ _ \ / _ \| / __|
# | (_| (_| | (_| \__ \ || (_) | (_) | \__ \
#  \___\__,_|\__,_|___/\__\___/ \___/|_|___/

__version__ = "2.0.0-paired"

'''
A Convergent Amino Acid Substitution identification
and analysis toolbox

Author:         Fabio Barteri (fabio.barteri@upf.edu)

Contributors:   Alejandro Valenzuela (alejandro.valenzuela@upf.edu)
                Xavier Farré (xfarrer@igtp.cat),
                David de Juan (david.juan@upf.edu).

Pair-aware implementation: Miguel Ramon (miguel.ramon@upf.edu)
CAAP implementation: Miguel Ramon (miguel.ramon@upf.edu)

MODULE NAME:    caap_id
DESCRIPTION:    Identification of CAAP (Convergent Amino Acid Properties) from MSA.
                Extends classical CAAS detection by grouping amino acids based on
                physicochemical properties and detecting convergence at the property level.
DEPENDENCIES:   pindex, alimport, hyper, caas_id


TABLE OF CONTENTS
------------------------------------------

Grouping Schemes
US    Classical CAAS (ungrouped; identity mapping: each amino acid is its own group)
GS1   Coarse biochemical recoding (6 groups; Dayhoff-like variant; custom partition)
GS2   Side-chain dipole/volume-inspired (7 groups; per Yang 2010 via Shen 2007)
GS3   Polarity and volume (6 groups; per Zhang 2000)
GS4   Fine-grained biochemical (12 groups; textbook-style functional bins)

Functions:
encode_to_groups()                          Encode amino acid string to group string for a given scheme
build_group_line_dictionary()               Transform position data to group-level format for p-value calc
_prepare_reduced_group_line_dict_for_pvalue() Exclude conserved FG/BG pairs from group-level p-value dict
check_caap_pattern()                        Check if a pattern is CAAP for a given scheme
fetch_caap()                                Main function: identifies CAAP for all traits at a position

----------

CAAP grouping follows the same logic found in:

Chen, S., & Zou, Z. (2025). Detecting Convergence of Amino Acid Physicochemical
Properties Underlying the Organismal Adaptive Convergent Evolution.
*Molecular Ecology Resources*, 25(1), e70052.
https://doi.org/10.1111/1755-0998.70052


'''

from modules.pindex import *
from modules.alimport import *
from modules.caas_id import process_position
from os.path import exists
from typing import Dict, List, Tuple, Optional


def _pair_sort_key(multiconfig, sp, trait=None):
    pair_id = multiconfig.get_pair(sp, trait) if trait is not None else multiconfig.get_pair(sp)
    if pair_id:
        try:
            return (int(pair_id), sp)
        except (ValueError, TypeError):
            return (float('inf'), sp)
    return (float('inf'), sp)


#: US: Ungrouped scheme (identity mapping; each amino acid is its own group)
US: Dict[str, str] = {
    "A": "A", "C": "C", "D": "D", "E": "E", "F": "F",
    "G": "G", "H": "H", "I": "I", "K": "K", "L": "L",
    "M": "M", "N": "N", "P": "P", "Q": "Q", "R": "R",
    "S": "S", "T": "T", "V": "V", "W": "W", "Y": "Y",
}


#: GS1: Coarse biochemical recoding (6 groups; Dayhoff-like variant?)
#: Partition: CV // AGPS // NDQE // RHK // ILMFWY // T
#: Notes: I have not been able to find a literature source for this exact partition

GS1: Dict[str, str] = {
    # t: cysteine/valine bin (special-case bin in this custom partition)
    "C": "t", "V": "t",

    # n(eutral-small): small / conformationally permissive residues (A,G,P,S)
    "A": "n", "G": "n", "P": "n", "S": "n",

    # p(olar): polar/amide/acid family (N,D,Q,E)
    "N": "p", "D": "p", "Q": "p", "E": "p",

    # b(asic): positively charged / ionizable (R,H,K)
    "R": "b", "H": "b", "K": "b",

    # h(ydrophobic): bulky non-polar incl. aromatics (I,L,M,F,W,Y)
    "I": "h", "L": "h", "M": "h", "F": "h", "W": "h", "Y": "h",

    # o(xy): hydroxyl-bearing singleton (T)
    "T": "o",
}


#: GS2: Side-chain dipole/volume-inspired (7 groups)
#: Partition: C // AGV // DE // NQHW // RK // ILFP // YMTS
#: Source: Yang 2010 (based on Shen 2007; 10.2174/092986610791760306)
GS2: Dict[str, str] = {
    # c(ysteine): "special thiol" (C)
    "C": "c",

    # s(mall): "small/simple" (A,G,V)
    "A": "s", "G": "s", "V": "s",

    # a(cidic): "negatively charged" (D,E)
    "D": "a", "E": "a",

    # n(itrogen): "N-mediated interactions / H-bonding / aromatic N" (N,Q,H,W)
    "N": "n", "Q": "n", "H": "n", "W": "n",

    # b(asic): "positively charged" (R,K)
    "R": "b", "K": "b",

    # h(ydrophobic): "core-forming / non-polar" incl. Pro (I,L,F,P)
    "I": "h", "L": "h", "F": "h", "P": "h",

    # x(hetero): "heteroatom-bearing / reactive or hydroxyl/sulfur" (Y,M,T,S)
    "Y": "x", "M": "x", "T": "x", "S": "x",
}


#: GS3: Polarity and volume (6 groups)
#: Partition: C // AGPST // NDQE // RHK // ILMV // FWY
#: Source: Zhang 2000, 10.1007/s002399910007
GS3: Dict[str, str] = {
    # c(ysteine): "special" (Cys)
    "C": "c",

    # n(eutral): "neutral and small" (A,G,P,S,T)
    "A": "n", "G": "n", "P": "n", "S": "n", "T": "n",

    # s(mall): "polar and relatively small" (N,D,Q,E)
    "N": "s", "D": "s", "Q": "s", "E": "s",

    # b(ig): "polar and relatively large" (R,H,K)
    "R": "b", "H": "b", "K": "b",

    # l(ittle): "non-polar and relatively small" (I,L,M,V)
    "I": "l", "L": "l", "M": "l", "V": "l",

    # g(rand): "non-polar and relatively large" (F,W,Y)
    "F": "g", "W": "g", "Y": "g",
}


#: GS4: Fine-grained biochemical (12 groups)
#: Partition: C // AILV // ST // NQ // DE // RH // G // P // K // M // F // WY
#: Notes: Textbook-style functional bins (historical source attribution may vary by edition).
GS4: Dict[str, str] = {
    # c(ysteine): thiol (C)
    "C": "c",

    # h(ydrophobic): aliphatic non-polar (A,I,L,V)
    "A": "h", "I": "h", "L": "h", "V": "h",

    # o(xy): hydroxyl (S,T)
    "S": "o", "T": "o",

    # p(olar): amide (N,Q)
    "N": "p", "Q": "p",

    # a(cidic): acidic (D,E)
    "D": "a", "E": "a",

    # b(asic): basic (R,H)
    "R": "b", "H": "b",

    # g(ly): glycine singleton (G)
    "G": "g",

    # r(igid): proline singleton (P)
    "P": "r",

    # k(lys): lysine singleton (K)
    "K": "k",

    # m(et): methionine singleton (M)
    "M": "m",

    # f(phenyl): phenylalanine singleton (F)
    "F": "f",

    # y(aromatic): aromatics (W,Y)
    "W": "y", "Y": "y",
}

#: Scheme registry mapping scheme names to their definitions
SCHEMES: Dict[str, Dict[str, str]] = {
    "US": US,
    "GS1": GS1,
    "GS2": GS2,
    "GS3": GS3,
    "GS4": GS4,
}


def encode_to_groups(amino_acids: str, scheme_dict: Dict[str, str]) -> str:
    """
    Encode a string of amino acids to their corresponding groups.

    Args:
        amino_acids: String of amino acids (e.g., "AAA" or "SGT")
        scheme_dict: Dictionary mapping amino acids to groups (e.g., GS1, GS2)

    Returns:
        String of concatenated group codes (e.g., "AGPSAGPSAGPS" for "AAA" in GS1)

    Example:
        >>> encode_to_groups("AAA", GS1)
        "AGPSAGPSAGPS"
        >>> encode_to_groups("RHK", GS1)
        "RHKRHKRHK"
    """
    encoded = ""
    for aa in amino_acids:
        if aa in scheme_dict:
            encoded += scheme_dict[aa]
    return encoded





def check_caap_pattern(fg_aas: str, bg_aas: str, scheme_dict: Dict[str, str],
                        max_conserved: int = 0, multiconfig=None,
                        fg_species_list=None, bg_species_list=None,
                        trait=None) -> Tuple[bool, str, str, str]:
    """
    Check if a foreground/background pattern represents a CAAP for a given scheme.

    This function encodes amino acids to groups, then checks:
    1. Groups on FG side vs groups on BG side have minimal overlap
    2. Overlap is validated against max_conserved parameter
    3. At least 2 amino acids changed on at least one side

    Args:
        fg_aas: Foreground amino acids as string (e.g., "AAA" or "AAS")
        bg_aas: Background amino acids as string (e.g., "RHK")
        scheme_dict: Dictionary mapping amino acids to groups
        max_conserved: Maximum number of amino acids allowed to overlap
        multiconfig: Multiconfig object for pair tracking (optional)
        fg_species_list: List of FG species (optional, for pair tracking)
        bg_species_list: List of BG species (optional, for pair tracking)

    Returns:
        Tuple of (is_caap: bool, pattern: str, substitution: str, conserved_pairs: str)
        pattern is "1", "2", "3", "4", or "null"
        conserved_pairs is "overlap:pair1,pair2,..." or "0:"

    Example:
        >>> check_caap_pattern("AAG", "RHK", GS1, max_conserved=0)
        (True, "1", "AAG/RHK", "0:")  # Both converge: AGPS vs RHK
    """
    class caap_result:
        def __init__(self):
            self.caap = False
            self.pattern = "null"
            self.substitution = ""
            self.conserved_pairs = "0:"

    z = caap_result()
    if fg_species_list is None:
        fg_species_list = []
    if bg_species_list is None:
        bg_species_list = []
    z.substitution = f"{fg_aas}/{bg_aas}"

    fg_aas_filtered = fg_aas
    bg_aas_filtered = bg_aas

    # Encode amino acids to groups
    # Use get() with None as default and filter out None values for non-standard AAs
    fg_groups = [scheme_dict.get(aa) for aa in fg_aas_filtered if scheme_dict.get(aa) is not None]
    bg_groups = [scheme_dict.get(aa) for aa in bg_aas_filtered if scheme_dict.get(aa) is not None]

    # Get unique groups for each side
    fg_unique = set(fg_groups)
    bg_unique = set(bg_groups)

    # Pattern classification (same as caas_id.py iscaas())
    if len(fg_unique) == 1 and len(bg_unique) == 1:
        z.pattern = "1"
    elif len(fg_unique) == 1:
        z.pattern = "2"
    elif len(bg_unique) == 1:
        z.pattern = "3"
    elif len(fg_unique) > 1 and len(bg_unique) > 1:
        z.pattern = "4"

    if len(fg_unique) == 0 or len(bg_unique) == 0:
        z.pattern = "null"
        return (z.caap, z.pattern, z.substitution, z.conserved_pairs)

    # Calculate string overlap (same logic as classical CAAS)
    # Count non-overlapping changes based on membership in the opposite side.
    # Examples:
    #   AAA/ASS -> non_overlapping_bg=2 (S,S), non_overlapping_fg=0
    #   VM/VVV  -> non_overlapping_bg=0 (all V shared), non_overlapping_fg=1 (M)
    shared_types = set(fg_unique).intersection(set(bg_unique))
    shared_fg = sum(1 for g in fg_groups if g in shared_types)
    shared_bg = sum(1 for g in bg_groups if g in shared_types)
    overlap = min(shared_fg, shared_bg)
    non_overlapping_fg = sum(1 for g in fg_groups if g not in bg_unique)
    non_overlapping_bg = sum(1 for g in bg_groups if g not in fg_unique)

    # Standard CAAP check (same logic as classical CAAS)
    if max_conserved == 0:
        # Strict mode: no overlap allowed and need at least 2 changes on one side
        if overlap == 0 and (non_overlapping_fg >= 2 or non_overlapping_bg >= 2):
            z.caap = True
        else:
            z.caap = False
    else:
        # Allow overlap up to max_conserved regardless of paired mode
        if overlap <= max_conserved and (non_overlapping_fg >= 2 or non_overlapping_bg >= 2):
            z.caap = True
        else:
            z.caap = False

        # Track conserved pairs
        if z.caap and multiconfig:
            conserved_pair_indices = []
            min_len = min(len(fg_groups), len(bg_groups))

            for i in range(min_len):
                if fg_groups[i] == bg_groups[i]:
                    if i < len(fg_species_list) and i < len(bg_species_list):
                        fg_sp = fg_species_list[i]
                        pair_id = multiconfig.get_pair(fg_sp, trait)
                        if pair_id:
                            conserved_pair_indices.append(pair_id)

            pair_list = ",".join(conserved_pair_indices) if conserved_pair_indices else ""
            z.conserved_pairs = f"{overlap}:{pair_list}"

    return (z.caap, z.pattern, z.substitution, z.conserved_pairs)


def fetch_caap(genename: str, position_obj, trait_list: List[str],
               max_fg_gaps: int, max_bg_gaps: int, max_overall_gaps: int,
               max_fg_miss: int, max_bg_miss: int, max_overall_miss: int,
               output_file: Optional[str],
               miss_pair: bool = False, max_conserved: int = 0,
               species_in_alignment: Optional[List[str]] = None,
               allowed_patterns: Optional[List[str]] = None,
               multiconfig=None, return_results: bool = False):
    """
    Identify CAAP (Convergent Amino Acid Properties) for all traits at a position.

    For each position, checks all 4 grouping schemes (GS1-GS4). When a scheme detects
    convergence, writes one output row with scheme-specific p-value calculated from
    group-level diversity.

    Args:
        genename: Gene name
        position_obj: Position object from process_position()
        trait_list: List of trait names to analyze
        max_fg_gaps: Maximum gaps allowed in foreground
        max_bg_gaps: Maximum gaps allowed in background
        max_overall_gaps: Maximum gaps overall
        max_fg_miss: Maximum missing species in foreground
        max_bg_miss: Maximum missing species in background
        max_overall_miss: Maximum missing overall
        output_file: Path to output file (or None if return_results=True)
        miss_pair: Whether to allow missing pairs
        max_conserved: Maximum number of amino acids allowed to overlap with other side's groups
        species_in_alignment: List of species in alignment (for p-value calculation)
        return_results: If True, return list of result lines instead of writing to file

    Output columns:
        gene, mode, caap_group, trait, position, caas, amino_encoded, pvalue, pattern,
        ffgn, fbgn, gfg, gbg, mfg, mbg, ffg, fbg, ms [, is_conserved_meta, conserved_pair]
    """

    # Ensure species_in_alignment is provided
    if species_in_alignment is None:
        species_in_alignment = []

    # Collect results if return_results=True
    result_lines = []

    # Filter traits by gap and missing thresholds
    valid_traits = []
    for trait in trait_list:
        if trait not in position_obj.trait2aas_fg:
            continue
        if trait not in position_obj.trait2aas_bg:
            continue

        # Check gap thresholds
        if position_obj.trait2gaps_fg.get(trait, 0) > max_fg_gaps:
            continue
        if position_obj.trait2gaps_bg.get(trait, 0) > max_bg_gaps:
            continue
        if position_obj.trait2gaps_all.get(trait, 0) > max_overall_gaps:
            continue

        # Check missing species thresholds
        if position_obj.trait2miss_fg.get(trait, 0) > max_fg_miss:
            continue
        if position_obj.trait2miss_bg.get(trait, 0) > max_bg_miss:
            continue
        if position_obj.trait2miss_all.get(trait, 0) > max_overall_miss:
            continue

        # Pair-aware filtering
        if miss_pair:
            # Check if missing thresholds are equal
            miss_thresholds_equal = False
            if max_fg_miss < 999999 and max_bg_miss < 999999 and max_fg_miss == max_bg_miss:
                miss_thresholds_equal = True

            if miss_thresholds_equal:
                miss_pairs_fg = set(position_obj.trait2miss_pairs_fg.get(trait, []))
                miss_pairs_bg = set(position_obj.trait2miss_pairs_bg.get(trait, []))
                # Only check pairs when BOTH FG and BG have missing species
                # If only one side missing → always OK (don't check pairs)
                # If both missing → only reject if from different pairs
                if miss_pairs_fg and miss_pairs_bg and miss_pairs_fg != miss_pairs_bg:
                    continue

            # Check if gap thresholds are equal
            gap_thresholds_equal = False
            if max_fg_gaps < 999999 and max_bg_gaps < 999999 and max_fg_gaps == max_bg_gaps:
                gap_thresholds_equal = True

            if gap_thresholds_equal:
                gap_pairs_fg = set(position_obj.trait2gap_pairs_fg.get(trait, []))
                gap_pairs_bg = set(position_obj.trait2gap_pairs_bg.get(trait, []))
                # Only check pairs when BOTH FG and BG have gapped species
                # If only one side gapped → always OK (don't check pairs)
                # If both gapped → only reject if from different pairs
                if gap_pairs_fg and gap_pairs_bg and gap_pairs_fg != gap_pairs_bg:
                    continue

        valid_traits.append(trait)

    # Process each valid trait across all schemes
    for trait in valid_traits:
        # Get foreground and background species lists (ungapped)
        fg_species = position_obj.trait2ungapped_fg.get(trait, [])[:]
        bg_species = position_obj.trait2ungapped_bg.get(trait, [])[:]

        if not fg_species or not bg_species:
            continue

        # Always sort species by pair ID (mandatory paired mode)
        fg_species.sort(key=lambda sp: _pair_sort_key(multiconfig, sp, trait))
        bg_species.sort(key=lambda sp: _pair_sort_key(multiconfig, sp, trait))

        # Extract amino acid for each species to create pattern string
        # Same as classical CAAS: get AA from position_obj.d[species]
        fg_aas = "".join([position_obj.d[sp].split("@")[0] for sp in fg_species])
        bg_aas = "".join([position_obj.d[sp].split("@")[0] for sp in bg_species])

        if not fg_aas or not bg_aas:
            continue

        # Check each grouping scheme
        for scheme_name, scheme_dict in SCHEMES.items():
            # Check if this scheme detects a CAAP
            is_caap, pattern, substitution, conserved_pairs = check_caap_pattern(
                fg_aas, bg_aas, scheme_dict, max_conserved,
                multiconfig, fg_species, bg_species, trait
            )

            if not is_caap:
                continue

            # Filter by allowed patterns
            if allowed_patterns is not None and pattern not in allowed_patterns:
                continue

            # Encode amino acids to groups for this scheme
            encoded = encode_to_groups(fg_aas, scheme_dict) + "/" + \
                     encode_to_groups(bg_aas, scheme_dict)

            # Prepare output row
            # Preserve the same sorted pair order used for substitution/encoded strings
            fg_species_ordered = fg_species[:]
            bg_species_ordered = bg_species[:]
            miss_species = position_obj.trait2missings.get(trait, [])

            # Format species lists
            fg_species_str = ",".join(fg_species_ordered) if fg_species_ordered else "NA"
            bg_species_str = ",".join(bg_species_ordered) if bg_species_ordered else "NA"
            miss_species_str = ",".join(miss_species) if miss_species else "NA"

            # Build output line
            output_line = [
                genename,
                "CAAP",
                scheme_name,
                trait,
                str(position_obj.position),
                substitution,
                encoded,
                pattern,
                str(len(fg_species_ordered)),  # FFGN - Final foreground number
                str(len(bg_species_ordered)),  # FBGN - Final background number
                str(position_obj.trait2gaps_fg.get(trait, 0)),  # GFG
                str(position_obj.trait2gaps_bg.get(trait, 0)),  # GBG
                str(position_obj.trait2miss_fg.get(trait, 0)),  # MFG
                str(position_obj.trait2miss_bg.get(trait, 0)),  # MBG
                fg_species_str,  # FFG
                bg_species_str,  # FBG
                miss_species_str,  # MS
            ]

            # Add conserved-pair columns when overlap tolerance is enabled
            if max_conserved > 0:
                # Parse conserved_pairs to get overlap count and pair list
                parts = conserved_pairs.split(":")
                overlap_count = int(parts[0]) if parts[0] else 0
                pair_list = parts[1] if len(parts) > 1 else ""

                # is_conserved_meta column: TRUE if overlap > 0, FALSE otherwise
                has_conserved = "TRUE" if overlap_count > 0 else "FALSE"
                output_line.append(has_conserved)
                output_line.append(conserved_pairs)

            # Write to output file or collect for return
            output_line_str = "\t".join(output_line)
            if return_results:
                result_lines.append(output_line_str)
            elif output_file:
                with open(output_file, "a") as outf:
                    outf.write(output_line_str + "\n")

    # Return results if requested
    if return_results:
        return result_lines
