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

Grouping Schemes:
GS0                         Classical CAAS (identity mapping, each AA is its own group)
GS1                         Physicochemical properties (6 groups)
GS2                         Charge and size (7 groups)
GS3                         Hydrophobicity and function (6 groups)
GS4                         Fine-grained biochemical (12 groups)
SCHEMES                     Registry mapping scheme names to definitions

Functions:
encode_to_groups()          Encode amino acid string to group string for a given scheme
build_group_line_dictionary() Transform position data to group-level format for p-value calc
check_caap_pattern()        Check if a pattern is CAAP for a given scheme
fetch_caap()                Main function: identifies CAAP for all traits at a position
'''

from modules.pindex import *
from modules.alimport import *
from modules.hyper import *
from modules.caas_id import process_position
from os.path import exists
from typing import Dict, List, Tuple, Optional

#: GS0: Classical CAAS group (identity mapping)
#: Each amino acid is its own group
GS0: Dict[str, str] = {
    "A": "A", "C": "C", "D": "D", "E": "E", "F": "F",
    "G": "G", "H": "H", "I": "I", "K": "K", "L": "L",
    "M": "M", "N": "N", "P": "P", "Q": "Q", "R": "R",
    "S": "S", "T": "T", "V": "V", "W": "W", "Y": "Y",
}

#: GS1: Physicochemical properties (6 groups)
#: CV // AGPS // NDQE // RHK // ILMFWY // T
#: Encoded as single characters: 1, 2, 3, 4, 5, 6
GS1: Dict[str, str] = {
    "C": "1",
    "V": "1",
    "A": "2",
    "G": "2",
    "P": "2",
    "S": "2",
    "N": "3",
    "D": "3",
    "Q": "3",
    "E": "3",
    "R": "4",
    "H": "4",
    "K": "4",
    "I": "5",
    "L": "5",
    "M": "5",
    "F": "5",
    "W": "5",
    "Y": "5",
    "T": "6",
}

#: GS2: Charge and size (7 groups)
#: C // AGV // DE // NQHW // RK // ILFP // YMTS
#: Encoded as: A, B, C, D, E, F, G
GS2: Dict[str, str] = {
    "C": "A",
    "A": "B",
    "G": "B",
    "V": "B",
    "D": "C",
    "E": "C",
    "N": "D",
    "Q": "D",
    "H": "D",
    "W": "D",
    "R": "E",
    "K": "E",
    "I": "F",
    "L": "F",
    "F": "F",
    "P": "F",
    "Y": "G",
    "M": "G",
    "T": "G",
    "S": "G",
}

#: GS3: Hydrophobicity and function (6 groups)
#: C // AGPST // NDQE // RHK // ILMV // FWY
#: Encoded as: X, Y, Z, U, V, W
GS3: Dict[str, str] = {
    "C": "X",
    "A": "Y",
    "G": "Y",
    "P": "Y",
    "S": "Y",
    "T": "Y",
    "N": "Z",
    "D": "Z",
    "Q": "Z",
    "E": "Z",
    "R": "U",
    "H": "U",
    "K": "U",
    "I": "V",
    "L": "V",
    "M": "V",
    "V": "V",
    "F": "W",
    "W": "W",
    "Y": "W",
}

#: GS4: Fine-grained biochemical (12 groups)
#: C // AILV // ST // NQ // DE // RH // G // P // K // M // F // WY
#: Encoded as: a, b, c, d, e, f, g, h, i, j, k, l
GS4: Dict[str, str] = {
    "C": "a",
    "A": "b",
    "I": "b",
    "L": "b",
    "V": "b",
    "S": "c",
    "T": "c",
    "N": "d",
    "Q": "d",
    "D": "e",
    "E": "e",
    "R": "f",
    "H": "f",
    "G": "g",
    "P": "h",
    "K": "i",
    "M": "j",
    "F": "k",
    "W": "l",
    "Y": "l",
}

#: Scheme registry mapping scheme names to their definitions
SCHEMES: Dict[str, Dict[str, str]] = {
    "GS0": GS0,
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


def build_group_line_dictionary(position_obj, trait: str, scheme_dict: Dict[str, str], 
                                  species_in_alignment: List[str]) -> Dict[str, str]:
    """
    Transform position data to group-level format for p-value calculation.
    
    This function creates a line_dictionary compatible with calcpval_random() from hyper.py,
    but with amino acids replaced by their group codes according to the scheme.
    
    Args:
        position_obj: Position object from process_position() containing aa2spp, aa2trait, etc.
        trait: Trait name to process
        scheme_dict: Dictionary mapping amino acids to groups
        species_in_alignment: List of all species in the alignment
    
    Returns:
        Dictionary mapping species to "GROUP@POSITION" format (e.g., {"Human": "AGPS@42"})
    """
    line_dict = {}
    
    # Get position number from position_obj
    position_num = position_obj.position
    
    # For each species in the alignment, find its amino acid and encode to group
    for species in species_in_alignment:
        # Find which amino acid this species has
        aa_found = None
        for aa, species_list in position_obj.aas2species.items():
            if species in species_list:
                aa_found = aa
                break
        
        # Encode amino acid to group
        if aa_found and aa_found != "-":  # Skip gaps
            group = scheme_dict.get(aa_found, aa_found)
            line_dict[species] = f"{group}@{position_num}"
    
    return line_dict


def check_caap_pattern(fg_aas: str, bg_aas: str, scheme_dict: Dict[str, str], 
                        max_conserved: int = 0, multiconfig=None, 
                        fg_species_list=None, bg_species_list=None) -> Tuple[bool, str, str, str]:
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
        Tuple of (is_caap: bool, pattern_type: str, substitution: str, conserved_pairs: str)
        pattern_type is "1", "2", "3", "4", or "null"
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
    z.substitution = f"{fg_aas}/{bg_aas}"
    
    # Encode amino acids to groups
    fg_groups = [scheme_dict.get(aa, aa) for aa in fg_aas if aa != "-"]
    bg_groups = [scheme_dict.get(aa, aa) for aa in bg_aas if aa != "-"]
    
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
        return (z.caap, z.pattern, z.substitution)
    
    # Calculate string overlap (same logic as classical CAAS)
    # Count how many characters from fg_groups also appear in bg_groups (with multiplicity)
    # For example: AAA/ASS has overlap=1 (one A from left matches one A in right)
    #              AAA/AAS has overlap=2 (two A's from left match two A's in right)
    
    # Count frequency of each group on both sides
    from collections import Counter
    fg_counter = Counter(fg_groups)
    bg_counter = Counter(bg_groups)
    
    # Overlap is the sum of minimum counts for each group
    overlap = sum(min(fg_counter[group], bg_counter[group]) for group in fg_counter)
    
    # Calculate non-overlapping changes
    total_fg = len(fg_groups)
    total_bg = len(bg_groups)
    non_overlapping_fg = total_fg - overlap
    non_overlapping_bg = total_bg - overlap
    
    # Standard CAAP check (same logic as classical CAAS)
    if max_conserved == 0 or not multiconfig or not multiconfig.paired_mode:
        # Strict mode: no overlap allowed and need at least 2 changes on one side
        if overlap == 0 and (non_overlapping_fg >= 2 or non_overlapping_bg >= 2):
            z.caap = True
        else:
            z.caap = False
    else:
        # Pair-aware mode: allow some overlap based on max_conserved
        # Track conserved pairs for reporting
        conserved_pair_indices = []
        min_len = min(len(fg_groups), len(bg_groups))
        
        for i in range(min_len):
            if fg_groups[i] == bg_groups[i]:
                # Track which pair this position corresponds to
                if i < len(fg_species_list) and i < len(bg_species_list):
                    fg_sp = fg_species_list[i]
                    pair_id = multiconfig.get_pair(fg_sp)
                    if pair_id:
                        conserved_pair_indices.append(pair_id)
        
        # Check if overlap is within threshold and we have enough changes
        if overlap <= max_conserved and (non_overlapping_fg >= 2 or non_overlapping_bg >= 2):
            z.caap = True
            pair_list = ",".join(conserved_pair_indices) if conserved_pair_indices else ""
            z.conserved_pairs = f"{overlap}:{pair_list}"
        else:
            z.caap = False
    
    return (z.caap, z.pattern, z.substitution, z.conserved_pairs)


def fetch_caap(genename: str, position_obj, trait_list: List[str], 
               max_fg_gaps: int, max_bg_gaps: int, max_overall_gaps: int,
               max_fg_miss: int, max_bg_miss: int, max_overall_miss: int,
               output_file: str, paired_mode: bool = False, 
               miss_pair: bool = False, max_conserved: int = 0,
               species_in_alignment: Optional[List[str]] = None,
               allowed_patterns: Optional[List[str]] = None,
               multiconfig=None):
    """
    Identify CAAP (Convergent Amino Acid Properties) for all traits at a position.
    
    For each position, checks all 5 grouping schemes (GS0-GS4). When a scheme detects
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
        output_file: Path to output file
        paired_mode: Whether to use pair-aware mode
        miss_pair: Whether to allow missing pairs
        max_conserved: Maximum number of amino acids allowed to overlap with other side's groups
        species_in_alignment: List of species in alignment (for p-value calculation)
    
    Output columns:
        Gene, Mode, CAAP_Group, Trait, Position, Substitution, Encoded, Pvalue, Pattern,
        FFGN, FBGN, GFG, GBG, MFG, MBG, FFG, FBG, MS [, ConservedPair, ConservedPairs]
    """
    
    # Ensure species_in_alignment is provided
    if species_in_alignment is None:
        species_in_alignment = []
    
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
        
        # Pair-aware filtering (only in paired mode with miss_pair enabled)
        if paired_mode and miss_pair:
            # Check if missing thresholds are equal
            miss_thresholds_equal = False
            if max_fg_miss < 999999 and max_bg_miss < 999999 and max_fg_miss == max_bg_miss:
                miss_thresholds_equal = True
            elif max_fg_miss == 999999 and max_bg_miss == 999999 and max_overall_miss < 999999:
                miss_thresholds_equal = True
            
            if miss_thresholds_equal:
                miss_pairs_fg = set(position_obj.trait2miss_pairs_fg.get(trait, []))
                miss_pairs_bg = set(position_obj.trait2miss_pairs_bg.get(trait, []))
                if miss_pairs_fg != miss_pairs_bg:
                    continue
            
            # Check if gap thresholds are equal
            gap_thresholds_equal = False
            if max_fg_gaps < 999999 and max_bg_gaps < 999999 and max_fg_gaps == max_bg_gaps:
                gap_thresholds_equal = True
            elif max_fg_gaps == 999999 and max_bg_gaps == 999999 and max_overall_gaps < 999999:
                gap_thresholds_equal = True
            
            if gap_thresholds_equal:
                gap_pairs_fg = set(position_obj.trait2gap_pairs_fg.get(trait, []))
                gap_pairs_bg = set(position_obj.trait2gap_pairs_bg.get(trait, []))
                if gap_pairs_fg != gap_pairs_bg:
                    continue
        
        valid_traits.append(trait)
    
    # Process each valid trait across all schemes
    for trait in valid_traits:
        # Get foreground and background species lists (ungapped)
        fg_species = position_obj.trait2ungapped_fg.get(trait, [])[:]
        bg_species = position_obj.trait2ungapped_bg.get(trait, [])[:]
        
        if not fg_species or not bg_species:
            continue
        
        # Sort species by pair ID if in paired mode, otherwise alphabetically
        # This ensures position indices correspond to correct pairs for conserved tracking
        if paired_mode and multiconfig:
            def pair_sort_key(sp):
                pair_id = multiconfig.get_pair(sp)
                if pair_id:
                    try:
                        return (int(pair_id), sp)
                    except (ValueError, TypeError):
                        return (float('inf'), sp)  # Non-numeric pairs go to end
                return (float('inf'), sp)  # No pair goes to end
            
            fg_species.sort(key=pair_sort_key)
            bg_species.sort(key=pair_sort_key)
        else:
            fg_species.sort()
            bg_species.sort()
        
        # Extract amino acid for each species to create pattern string
        # Same as classical CAAS: get AA from position_obj.d[species]
        fg_aas = "".join([position_obj.d[sp].split("@")[0] for sp in fg_species])
        bg_aas = "".join([position_obj.d[sp].split("@")[0] for sp in bg_species])
        
        if not fg_aas or not bg_aas:
            continue
        
        print(f"[CAAP] Position {position_obj.position}, Trait {trait}: FG_species={fg_species}, BG_species={bg_species}, FG={fg_aas}, BG={bg_aas}")
        
        # Check each grouping scheme
        for scheme_name, scheme_dict in SCHEMES.items():
            # Check if this scheme detects a CAAP
            is_caap, pattern, substitution, conserved_pairs = check_caap_pattern(
                fg_aas, bg_aas, scheme_dict, max_conserved,
                multiconfig, fg_species, bg_species
            )
            
            if not is_caap:
                continue
            
            # Filter by allowed patterns
            if allowed_patterns is not None and pattern not in allowed_patterns:
                continue
            
            print(f"CAAP found in alignment {genename} on position {position_obj.position} for scheme {scheme_name} with pattern {pattern}")
            
            # Encode amino acids to groups for this scheme
            encoded = encode_to_groups(fg_aas, scheme_dict) + "/" + \
                     encode_to_groups(bg_aas, scheme_dict)
            
            # Calculate p-value using group-level data
            # Build line_dictionary with single-character group codes for ALL species
            # The hypergeometric calculation uses all species in alignment as the population
            line_dict = build_group_line_dictionary(
                position_obj, trait, scheme_dict, species_in_alignment
            )
            
            # Calculate group-based p-value
            fg_size = len(fg_species)
            bg_size = len(bg_species)
            
            pvalue = calcpval_random(line_dict, genename, fg_size, bg_size)
            
            # Prepare output row
            # Get species lists
            fg_species = position_obj.trait2ungapped_fg.get(trait, [])
            bg_species = position_obj.trait2ungapped_bg.get(trait, [])
            miss_species = position_obj.trait2missings.get(trait, [])
            
            # Format species lists
            fg_species_str = ",".join(fg_species) if fg_species else "NA"
            bg_species_str = ",".join(bg_species) if bg_species else "NA"
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
                str(pvalue),
                pattern,
                str(len(fg_species)),  # FFGN - Final foreground number
                str(len(bg_species)),  # FBGN - Final background number
                str(position_obj.trait2gaps_fg.get(trait, 0)),  # GFG
                str(position_obj.trait2gaps_bg.get(trait, 0)),  # GBG
                str(position_obj.trait2miss_fg.get(trait, 0)),  # MFG
                str(position_obj.trait2miss_bg.get(trait, 0)),  # MBG
                fg_species_str,  # FFG
                bg_species_str,  # FBG
                miss_species_str,  # MS
            ]
            
            # Add pair-aware columns if in paired mode
            if paired_mode and max_conserved > 0:
                # Parse conserved_pairs to get overlap count and pair list
                parts = conserved_pairs.split(":")
                overlap_count = int(parts[0]) if parts[0] else 0
                pair_list = parts[1] if len(parts) > 1 else ""
                
                # ConservedPair column: TRUE if overlap > 0, FALSE otherwise
                has_conserved = "TRUE" if overlap_count > 0 else "FALSE"
                output_line.append(has_conserved)
                output_line.append(conserved_pairs)
            
            # Write to output file
            with open(output_file, "a") as outf:
                outf.write("\t".join(output_line) + "\n")
