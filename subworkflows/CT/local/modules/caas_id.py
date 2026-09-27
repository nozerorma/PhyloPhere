#                      _              _     
#                     | |            | |    
#   ___ __ _  __ _ ___| |_ ___   ___ | |___ 
#  / __/ _` |/ _` / __| __/ _ \ / _ \| / __|
# | (_| (_| | (_| \__ \ || (_) | (_) | \__ \
#  \___\__,_|\__,_|___/\__\___/ \___/|_|___/

__version__ = "2.1.0-unified"

'''
A Convergent Amino Acid Substitution and Property (CAAS / CAAP)
identification and analysis toolbox.

Consolidated module providing unified convergence detection across amino acid
identities (classical CAAS / scheme US) and physicochemical property grouping
schemes (GS1 - GS4).

Author:         Fabio Barteri (fabio.barteri@upf.edu)

Contributors:   Alejandro Valenzuela (alejandro.valenzuela@upf.edu)
                Xavier Farré (xfarrer@igtp.cat),
                David de Juan (david.juan@upf.edu).

Pair-aware implementation: Miguel Ramon (miguel.ramon@upf.edu)
CAAP implementation:       Miguel Ramon (miguel.ramon@upf.edu)

MODULE NAME:    caas_id (consolidated CAAS & CAAP)
DESCRIPTION:    Identification of convergent mutations and property shifts from MSA.
DEPENDENCIES:   pindex, alimport

TABLE OF CONTENTS
------------------------------------------
Grouping Schemes:
US    Classical CAAS (ungrouped; identity mapping: each amino acid is its own group)
GS1   Coarse biochemical recoding (6 groups; Dayhoff-like variant)
GS2   Side-chain dipole/volume-inspired (7 groups; Yang 2010 via Shen 2007)
GS3   Polarity and volume (6 groups; Zhang 2000)
GS4   Fine-grained biochemical (12 groups; textbook functional bins)

Functions:
_pair_sort_key()            Key function for sorting species by pair ID
process_position()          Processes a position from an imported alignment
encode_to_groups()          Encodes an amino acid string to group codes
check_pattern()             Core convergence & pattern check for any scheme (US or GS1-4)
check_caap_pattern()        Alias for check_pattern()
iscaas()                    Legacy wrapper for check_pattern() using scheme US
fetch_caas()                Main discovery function across specified or all schemes
fetch_caap()                Alias for fetch_caas()
'''

from os.path import exists
from typing import Dict, List, Tuple, Optional, Any

from modules.pindex import *
from modules.alimport import *





def _pair_sort_key(multiconfig, sp, trait=None):
    pair_id = multiconfig.get_pair(sp, trait) if trait is not None else multiconfig.get_pair(sp)
    if pair_id:
        try:
            return (int(pair_id), sp)
        except (ValueError, TypeError):
            return (float('inf'), sp)
    return (float('inf'), sp)

# Function process_position()
# processes a position from an imported alignment. The output will be
# a dictionary that points each aminoacid (gaps included) to the
# species sharing it.

def process_position(position, multiconfig, species_in_alignment):

    class processed_position():
        def __init__(self):
            self.position = ""
            self.aas2species = {}
            self.aas2traits = {}
            self.trait2aas_fg = {}
            self.trait2aas_bg = {}

            self.trait2ungapped_fg = {}
            self.trait2ungapped_bg = {}

            self.trait2gaps_fg = {}
            self.trait2gaps_bg = {}
            self.trait2gaps_all = {}

            self.trait2miss_fg = {}
            self.trait2miss_bg = {}
            self.trait2miss_all = {}

            self.trait2missings = {}

            self.gapped = []
            self.missing = []

            self.d = {}
            
            # Pair-aware attributes
            self.trait2miss_pairs_fg = {}
            self.trait2miss_pairs_bg = {}
            self.trait2gap_pairs_fg = {}
            self.trait2gap_pairs_bg = {}
        
    z = processed_position()
    z.d = position

    # Load missing species

    confirmed_species = set(multiconfig.s2t.keys()).intersection(set(species_in_alignment))
    z.missing = list(set(multiconfig.s2t.keys()) -  confirmed_species)

    # Load aas2species
    for x in position.keys():
        z.position = position[x].split("@")[1]
        aa = position[x].split("@")[0].upper()
        
        try:
            z.aas2species[aa].append(x)
        except:
            z.aas2species[aa] = [x]

    # Load aas2traits    
    for key in z.aas2species.keys():
        traits = []
        for species in z.aas2species[key]:
            try:
                for v in multiconfig.s2t[species]:
                    if v not in traits:
                        traits.append(v)
            except:
                pass
        
        z.aas2traits[key] = traits

        for t in traits:
            if t[-2:] == "_1" and key != "-":
                try:
                    z.trait2aas_fg[t[:-2]].append(key)
                except:
                    z.trait2aas_fg[t[:-2]] = [key]
                    pass
            if t[-2:] == "_0" and key != "-":
                try:
                    z.trait2aas_bg[t[:-2]].append(key)
                except:
                    z.trait2aas_bg[t[:-2]] = [key]
                    pass

    try:
        z.gapped = z.aas2species["-"]
    except:
        pass

    # Treat IUPAC ambiguity codes as gaps (no resolved amino acid).
    # X = any, B = Asp/Asn, Z = Glu/Gln, J = Ile/Leu, U = selenocysteine.
    # Including these in z.gapped ensures they are excluded from ungapped_fg/bg
    # and counted towards the gap quota, preventing spurious convergence calls
    # where X would be silently dropped from group encoding while the species
    # remained in the foreground/background counts and species lists.
    _AMBIGUOUS_AAS = {'X', 'B', 'Z', 'J', 'U'}
    for _ambig in _AMBIGUOUS_AAS:
        if _ambig in z.aas2species:
            z.gapped = list(set(z.gapped) | set(z.aas2species[_ambig]))

    # Determine Ungapped Species

    for trait in z.trait2aas_bg.keys():
        
        # Present species (ungapped or missing)

        pfg = list(set(multiconfig.trait2fg[trait]) - set(z.gapped + z.missing))
        z.trait2ungapped_fg[trait] = pfg

        pbg = list(set(multiconfig.trait2bg[trait]) - set(z.gapped + z.missing))
        z.trait2ungapped_bg[trait] = pbg

        # Missing in alignment

        miss_fg = list(set(multiconfig.trait2fg[trait]).intersection(set(z.missing)))
        miss_bg = list(set(multiconfig.trait2bg[trait]).intersection(set(z.missing)))

        # Number of gaps
        gfg = len(set(multiconfig.trait2fg[trait]).intersection(set(z.gapped)))
        gbg = len(set(multiconfig.trait2bg[trait]).intersection(set(z.gapped)))

        # Number of missings
        mfg = len(miss_fg)
        mbg = len(miss_bg)

        gall = gfg + gbg
        mall = mfg + mbg

        z.trait2gaps_fg[trait] = gfg
        z.trait2gaps_bg[trait] = gbg
        z.trait2gaps_all[trait] = gall

        z.trait2miss_fg[trait] = mfg
        z.trait2miss_bg[trait] = mbg
        z.trait2miss_all[trait] = mall

        try:
            z.trait2missings[trait] = (miss_fg + miss_bg)
        except:
            z.trait2missings[trait] = "none"
        
        # Track pairs for missing and gapped species
        miss_pairs_fg = set()
        miss_pairs_bg = set()
        for sp in miss_fg:
            pair = multiconfig.get_pair(sp, trait)
            if pair:
                miss_pairs_fg.add(pair)
        for sp in miss_bg:
            pair = multiconfig.get_pair(sp, trait)
            if pair:
                miss_pairs_bg.add(pair)

        z.trait2miss_pairs_fg[trait] = list(miss_pairs_fg)
        z.trait2miss_pairs_bg[trait] = list(miss_pairs_bg)

        # Get pairs for gapped species
        gap_pairs_fg = set()
        gap_pairs_bg = set()
        gapped_fg = set(multiconfig.trait2fg[trait]).intersection(set(z.gapped))
        gapped_bg = set(multiconfig.trait2bg[trait]).intersection(set(z.gapped))

        for sp in gapped_fg:
            pair = multiconfig.get_pair(sp, trait)
            if pair:
                gap_pairs_fg.add(pair)
        for sp in gapped_bg:
            pair = multiconfig.get_pair(sp, trait)
            if pair:
                gap_pairs_bg.add(pair)

        z.trait2gap_pairs_fg[trait] = list(gap_pairs_fg)
        z.trait2gap_pairs_bg[trait] = list(gap_pairs_bg)

    return z


# ---------------------------------------------------------------------------
# Physicochemical Property Grouping Schemes
# ---------------------------------------------------------------------------

#: US: Ungrouped scheme (identity mapping; classical CAAS)
US: Dict[str, str] = {
    "A": "A", "C": "C", "D": "D", "E": "E", "F": "F",
    "G": "G", "H": "H", "I": "I", "K": "K", "L": "L",
    "M": "M", "N": "N", "P": "P", "Q": "Q", "R": "R",
    "S": "S", "T": "T", "V": "V", "W": "W", "Y": "Y",
}

#: GS1: Coarse biochemical recoding (6 groups)
#: CV // AGPS // NDQE // RHK // ILMFWY // T
GS1: Dict[str, str] = {
    "C": "t", "V": "t",
    "A": "n", "G": "n", "P": "n", "S": "n",
    "N": "p", "D": "p", "Q": "p", "E": "p",
    "R": "b", "H": "b", "K": "b",
    "I": "h", "L": "h", "M": "h", "F": "h", "W": "h", "Y": "h",
    "T": "o",
}

#: GS2: Side-chain dipole/volume-inspired (7 groups; Yang 2010 via Shen 2007)
#: C // AGV // DE // NQHW // RK // ILFP // YMTS
GS2: Dict[str, str] = {
    "C": "c",
    "A": "s", "G": "s", "V": "s",
    "D": "a", "E": "a",
    "N": "n", "Q": "n", "H": "n", "W": "n",
    "R": "b", "K": "b",
    "I": "h", "L": "h", "F": "h", "P": "h",
    "Y": "x", "M": "x", "T": "x", "S": "x",
}

#: GS3: Polarity and volume (6 groups; Zhang 2000)
#: C // AGPST // NDQE // RHK // ILMV // FWY
GS3: Dict[str, str] = {
    "C": "c",
    "A": "n", "G": "n", "P": "n", "S": "n", "T": "n",
    "N": "s", "D": "s", "Q": "s", "E": "s",
    "R": "b", "H": "b", "K": "b",
    "I": "l", "L": "l", "M": "l", "V": "l",
    "F": "g", "W": "g", "Y": "g",
}

#: GS4: Fine-grained biochemical (12 groups; textbook functional bins)
#: C // AILV // ST // NQ // DE // RH // G // P // K // M // F // WY
GS4: Dict[str, str] = {
    "C": "c",
    "A": "h", "I": "h", "L": "h", "V": "h",
    "S": "o", "T": "o",
    "N": "p", "Q": "p",
    "D": "a", "E": "a",
    "R": "b", "H": "b",
    "G": "g",
    "P": "r",
    "K": "k",
    "M": "m",
    "F": "f",
    "W": "y", "Y": "y",
}

#: Scheme registry mapping scheme names to dictionaries
SCHEMES: Dict[str, Dict[str, str]] = {
    "US": US,
    "GS1": GS1,
    "GS2": GS2,
    "GS3": GS3,
    "GS4": GS4,
}


def encode_to_groups(amino_acids: str, scheme_dict: Dict[str, str]) -> str:
    """Encode an amino acid string to its group code representation.

    Args:
        amino_acids: String of amino acids (e.g., 'AAA', 'SGT').
        scheme_dict: Mapping of amino acid single-letter codes to group codes.

    Returns:
        Encoded group string.
    """
    encoded = []
    for aa in amino_acids:
        g = scheme_dict.get(aa)
        if g is not None:
            encoded.append(g)
    return "".join(encoded)


# ---------------------------------------------------------------------------
# Convergence Pattern Checking
# ---------------------------------------------------------------------------

def check_pattern(
    fg_aas: str,
    bg_aas: str,
    scheme_dict: Dict[str, str] = US,
    max_conserved: int = 0,
    multiconfig: Any = None,
    fg_species_list: Optional[List[str]] = None,
    bg_species_list: Optional[List[str]] = None,
    trait: Optional[str] = None
) -> Tuple[bool, str, str, str]:
    """Check if foreground and background residue sets represent convergence.

    Generalized pattern and convergence test. Encodes amino acids through
    scheme_dict (identity for US, coarse-grained for GS1-GS4) and tests:
    1. Overlap <= max_conserved
    2. At least 2 non-overlapping residues on foreground or background
    3. Classification into pattern 1 (strict 1-to-1), 2 (1-to-many), 3 (many-to-1), or 4 (many-to-many).

    Args:
        fg_aas: Foreground amino acids string.
        bg_aas: Background amino acids string.
        scheme_dict: Mapping from amino acid to group code (default US).
        max_conserved: Overlap tolerance (0 for strict mode).
        multiconfig: Trait configuration object for pair recovery.
        fg_species_list: Sorted foreground species names.
        bg_species_list: Sorted background species names.
        trait: Trait name.

    Returns:
        (is_convergent, pattern, substitution, conserved_pairs)
    """
    substitution = f"{fg_aas}/{bg_aas}"
    if fg_species_list is None:
        fg_species_list = []
    if bg_species_list is None:
        bg_species_list = []

    # Encode amino acids to scheme groups
    fg_groups = [scheme_dict.get(aa) for aa in fg_aas if scheme_dict.get(aa) is not None]
    bg_groups = [scheme_dict.get(aa) for aa in bg_aas if scheme_dict.get(aa) is not None]

    fg_unique = set(fg_groups)
    bg_unique = set(bg_groups)

    if len(fg_unique) == 0 or len(bg_unique) == 0:
        return (False, "null", substitution, "0:")

    # Pattern classification
    if len(fg_unique) == 1 and len(bg_unique) == 1:
        pattern = "1"
    elif len(fg_unique) == 1:
        pattern = "2"
    elif len(bg_unique) == 1:
        pattern = "3"
    else:
        pattern = "4"

    # Overlap and divergence calculations
    shared_types = fg_unique.intersection(bg_unique)
    shared_fg = sum(1 for g in fg_groups if g in shared_types)
    shared_bg = sum(1 for g in bg_groups if g in shared_types)
    overlap = min(shared_fg, shared_bg)
    non_overlapping_fg = sum(1 for g in fg_groups if g not in bg_unique)
    non_overlapping_bg = sum(1 for g in bg_groups if g not in fg_unique)

    # Strict (overlap == 0) vs relaxed (overlap <= max_conserved) check
    if max_conserved == 0:
        is_convergent = (overlap == 0 and (non_overlapping_fg >= 2 or non_overlapping_bg >= 2))
    else:
        is_convergent = (overlap <= max_conserved and (non_overlapping_fg >= 2 or non_overlapping_bg >= 2))

    conserved_pairs = "0:"
    if is_convergent and max_conserved > 0 and multiconfig:
        conserved_pair_indices = []
        min_len = min(len(fg_groups), len(bg_groups))
        for i in range(min_len):
            if fg_groups[i] == bg_groups[i]:
                if i < len(fg_species_list) and i < len(bg_species_list):
                    fg_sp = fg_species_list[i]
                    pair_id = multiconfig.get_pair(fg_sp, trait)
                    if pair_id:
                        conserved_pair_indices.append(str(pair_id))
        pair_list = ",".join(conserved_pair_indices) if conserved_pair_indices else ""
        conserved_pairs = f"{overlap}:{pair_list}"

    return (is_convergent, pattern, substitution, conserved_pairs)


# Alias for backward-compatibility with caap_id
check_caap_pattern = check_pattern


def iscaas(input_string, multiconfig=None, position_dict=None, max_conserved=0, trait=None, fg_species_list=None, bg_species_list=None):
    """Backward-compatible wrapper for classical CAAS checking (scheme US)."""
    class caaspositive:
        def __init__(self):
            self.caas = False
            self.pattern = "4"
            self.conserved_pairs = "0:"

    z = caaspositive()
    twosides = input_string.split("/")
    fg_string = twosides[0]
    bg_string = twosides[1]

    is_match, pattern, _, conserved_pairs = check_pattern(
        fg_string, bg_string, scheme_dict=US,
        max_conserved=max_conserved, multiconfig=multiconfig,
        fg_species_list=fg_species_list, bg_species_list=bg_species_list,
        trait=trait
    )
    z.caas = is_match
    z.pattern = pattern
    z.conserved_pairs = conserved_pairs
    return z


# ---------------------------------------------------------------------------
# Discovery Function
# ---------------------------------------------------------------------------

def _parse_thresh(val: Any, default: int = 999999) -> int:
    if val is None or val == "NO":
        return default
    try:
        return int(val)
    except (ValueError, TypeError):
        return default


def fetch_caas(
    genename: str,
    position_obj: Any,
    trait_list: Optional[List[str]] = None,
    output_file: Optional[str] = None,
    max_fg_gaps: Any = 999999,
    max_bg_gaps: Any = 999999,
    max_overall_gaps: Any = 999999,
    max_fg_miss: Any = 999999,
    max_bg_miss: Any = 999999,
    max_overall_miss: Any = 999999,
    multiconfig: Any = None,
    miss_pair: bool = False,
    max_conserved: int = 0,
    admitted_patterns: Optional[List[str]] = None,
    schemes: Optional[Any] = None,
    return_results: bool = False,
    **kwargs
) -> List[str]:
    """Identify convergent substitutions / properties for all valid traits at a position.

    Outputs rows formatted according to the canonical 19-column schema:
    [gene, mode, caap_group, trait, position, caas, amino_encoded, pattern,
     ffgn, fbgn, gfg, gbg, mfg, mbg, ffg, fbg, ms [, is_conserved_meta, conserved_pair]]

    Args:
        genename: Gene name.
        position_obj: Processed position object from process_position().
        trait_list: List of traits to test (defaults to multiconfig.alltraits).
        output_file: Path to append output lines (if return_results=False).
        max_fg_gaps: Maximum gaps in foreground.
        max_bg_gaps: Maximum gaps in background.
        max_overall_gaps: Maximum gaps overall.
        max_fg_miss: Maximum missing species in foreground.
        max_bg_miss: Maximum missing species in background.
        max_overall_miss: Maximum missing species overall.
        multiconfig: Multiconfig object with trait definitions.
        miss_pair: Whether to enforce pair-aware gap/missing symmetry.
        max_conserved: Conserved pair / overlap tolerance.
        admitted_patterns: List of pattern codes to accept (defaults to ["1", "2", "3"]).
        schemes: Dictionary of schemes or list of scheme names to test (defaults to all 5 schemes).
        return_results: If True, return list of result lines instead of writing to output_file.

    Returns:
        List of formatted result lines (empty list if no hits or written to file).
    """
    # Backward compatibility with keyword argument names
    if trait_list is None:
        trait_list = kwargs.get("list_of_traits")
    if trait_list is None and multiconfig:
        trait_list = getattr(multiconfig, "alltraits", [])
    if trait_list is None:
        trait_list = []

    if admitted_patterns is None:
        admitted_patterns = kwargs.get("allowed_patterns", ["1", "2", "3"])

    # Resolve gap/missing thresholds
    fg_gaps = _parse_thresh(kwargs.get("maxgaps_fg", max_fg_gaps))
    bg_gaps = _parse_thresh(kwargs.get("maxgaps_bg", max_bg_gaps))
    overall_gaps = _parse_thresh(kwargs.get("maxgaps_all", max_overall_gaps))

    fg_miss = _parse_thresh(kwargs.get("maxmiss_fg", max_fg_miss))
    bg_miss = _parse_thresh(kwargs.get("maxmiss_bg", max_bg_miss))
    overall_miss = _parse_thresh(kwargs.get("maxmiss_all", max_overall_miss))

    # Resolve schemes to test
    if schemes is None:
        active_schemes = SCHEMES
    elif isinstance(schemes, dict):
        active_schemes = schemes
    elif isinstance(schemes, (list, tuple, set)):
        active_schemes = {s: SCHEMES[s] for s in schemes if s in SCHEMES}
    else:
        active_schemes = SCHEMES

    # Filter traits based on presence, gaps, and missing species
    valid_traits = []
    for trait in trait_list:
        if trait not in position_obj.trait2aas_fg or trait not in position_obj.trait2aas_bg:
            continue
        if position_obj.trait2gaps_fg.get(trait, 0) > fg_gaps:
            continue
        if position_obj.trait2gaps_bg.get(trait, 0) > bg_gaps:
            continue
        if position_obj.trait2gaps_all.get(trait, 0) > overall_gaps:
            continue
        if position_obj.trait2miss_fg.get(trait, 0) > fg_miss:
            continue
        if position_obj.trait2miss_bg.get(trait, 0) > bg_miss:
            continue
        if position_obj.trait2miss_all.get(trait, 0) > overall_miss:
            continue

        if miss_pair:
            # Check pair symmetry when thresholds are equal
            if fg_miss < 999999 and bg_miss < 999999 and fg_miss == bg_miss:
                m_fg = set(position_obj.trait2miss_pairs_fg.get(trait, []))
                m_bg = set(position_obj.trait2miss_pairs_bg.get(trait, []))
                if m_fg and m_bg and m_fg != m_bg:
                    continue
            elif fg_miss == 999999 and bg_miss == 999999 and overall_miss < 999999:
                m_fg = set(position_obj.trait2miss_pairs_fg.get(trait, []))
                m_bg = set(position_obj.trait2miss_pairs_bg.get(trait, []))
                if m_fg and m_bg and m_fg != m_bg:
                    continue

            if fg_gaps < 999999 and bg_gaps < 999999 and fg_gaps == bg_gaps:
                g_fg = set(position_obj.trait2gap_pairs_fg.get(trait, []))
                g_bg = set(position_obj.trait2gap_pairs_bg.get(trait, []))
                if g_fg and g_bg and g_fg != g_bg:
                    continue
            elif fg_gaps == 999999 and bg_gaps == 999999 and overall_gaps < 999999:
                g_fg = set(position_obj.trait2gap_pairs_fg.get(trait, []))
                g_bg = set(position_obj.trait2gap_pairs_bg.get(trait, []))
                if g_fg and g_bg and g_fg != g_bg:
                    continue

        valid_traits.append(trait)

    result_lines = []

    # Process each trait across active schemes
    for trait in valid_traits:
        fg_species = position_obj.trait2ungapped_fg.get(trait, [])[:]
        bg_species = position_obj.trait2ungapped_bg.get(trait, [])[:]

        if not fg_species or not bg_species:
            continue

        # Sort species by pair number
        fg_species.sort(key=lambda sp: _pair_sort_key(multiconfig, sp, trait))
        bg_species.sort(key=lambda sp: _pair_sort_key(multiconfig, sp, trait))

        fg_aas = "".join([position_obj.d[sp].split("@")[0] for sp in fg_species])
        bg_aas = "".join([position_obj.d[sp].split("@")[0] for sp in bg_species])

        if not fg_aas or not bg_aas:
            continue

        for scheme_name, scheme_dict in active_schemes.items():
            is_match, pattern, substitution, conserved_pairs = check_pattern(
                fg_aas, bg_aas, scheme_dict=scheme_dict,
                max_conserved=max_conserved, multiconfig=multiconfig,
                fg_species_list=fg_species, bg_species_list=bg_species,
                trait=trait
            )

            if not is_match:
                continue
            if admitted_patterns is not None and pattern not in admitted_patterns:
                continue

            encoded = (encode_to_groups(fg_aas, scheme_dict) + "/" +
                       encode_to_groups(bg_aas, scheme_dict))

            miss_species = position_obj.trait2missings.get(trait, [])
            fg_species_str = ",".join(fg_species) if fg_species else "NA"
            bg_species_str = ",".join(bg_species) if bg_species else "NA"
            miss_species_str = ",".join(miss_species) if miss_species else "NA"

            output_line = [
                genename,
                "CAAP",
                scheme_name,
                trait,
                str(position_obj.position),
                substitution,
                encoded,
                pattern,
                str(len(fg_species)),
                str(len(bg_species)),
                str(position_obj.trait2gaps_fg.get(trait, 0)),
                str(position_obj.trait2gaps_bg.get(trait, 0)),
                str(position_obj.trait2miss_fg.get(trait, 0)),
                str(position_obj.trait2miss_bg.get(trait, 0)),
                fg_species_str,
                bg_species_str,
                miss_species_str,
            ]

            if max_conserved > 0:
                parts = conserved_pairs.split(":")
                overlap_count = int(parts[0]) if parts[0] else 0
                has_conserved = "TRUE" if overlap_count > 0 else "FALSE"
                output_line.append(has_conserved)
                output_line.append(conserved_pairs)

            line_str = "\t".join(output_line)
            result_lines.append(line_str)
            if not return_results and output_file:
                with open(output_file, "a") as outf:
                    outf.write(line_str + "\n")

    return result_lines if return_results else []


# Alias for backward-compatibility with caap_id
fetch_caap = fetch_caas

