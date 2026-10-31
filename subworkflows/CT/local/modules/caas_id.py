#                      _              _     
#                     | |            | |    
#   ___ __ _  __ _ ___| |_ ___   ___ | |___ 
#  / __/ _` |/ _` / __| __/ _ \ / _ \| / __|
# | (_| (_| | (_| \__ \ || (_) | (_) | \__ \
#  \___\__,_|\__,_|___/\__\___/ \___/|_|___/

__version__ = "2.1.0-unified"

'''
A Convergent Amino Acid Substitution (CAAS) identification and analysis toolbox.

Convergence test on amino acid identities: the residues of the foreground and background
species of a contrast are compared residue by residue (each of the 20 standard amino acids is
its own state).

Author:         Fabio Barteri (fabio.barteri@upf.edu)

Contributors:   Alejandro Valenzuela (alejandro.valenzuela@upf.edu)
                Xavier Farré (xfarrer@igtp.cat),
                David de Juan (david.juan@upf.edu).

Pair-aware implementation: Miguel Ramon (miguel.ramon@upf.edu)

MODULE NAME:    caas_id
DESCRIPTION:    Identification of convergent substitutions from MSA.
DEPENDENCIES:   pindex, alimport

TABLE OF CONTENTS
------------------------------------------
_pair_sort_key()            Key function for sorting species by pair ID
process_position()          Processes a position from an imported alignment
check_pattern()             Convergence and pattern check of a foreground/background residue pair
'''

from os.path import exists
from typing import List, Tuple, Optional, Any

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
        aa = normalize_residue(position[x].split("@")[0])
        
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
    _AMBIGUOUS_AAS = GAP_SYMBOLS - {'-'}
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


#: The 20 standard amino acids; any other symbol is not a state of the convergence test.
def residue_overlap(fg_aas, bg_aas):
    fg = [aa for aa in fg_aas.upper() if aa in STANDARD_AAS]
    bg = [aa for aa in bg_aas.upper() if aa in STANDARD_AAS]
    shared = set(fg) & set(bg)
    return min(sum(aa in shared for aa in fg), sum(aa in shared for aa in bg))


# ---------------------------------------------------------------------------
# Convergence Pattern Checking
# ---------------------------------------------------------------------------

def check_pattern(
    fg_aas: str,
    bg_aas: str,
    max_conserved: int = 0,
    multiconfig: Any = None,
    fg_species_list: Optional[List[str]] = None,
    bg_species_list: Optional[List[str]] = None,
    trait: Optional[str] = None,
    min_non_overlapping: int = 2
) -> Tuple[bool, str, str, str]:
    """Check if foreground and background residue sets represent convergence.

    Keeps the standard amino acids as states and tests:
    1. Overlap <= max_conserved
    2. At least min_non_overlapping non-overlapping residues on foreground or background
    3. Classification into pattern 1 (strict 1-to-1), 2 (1-to-many), 3 (many-to-1), or 4 (many-to-many).

    Args:
        fg_aas: Foreground amino acids string.
        bg_aas: Background amino acids string.
        max_conserved: Overlap tolerance (0 for strict mode).
        multiconfig: Trait configuration object for pair recovery.
        fg_species_list: Sorted foreground species names.
        bg_species_list: Sorted background species names.
        trait: Trait name.
        min_non_overlapping: Minimum non-overlapping residues required on FG or BG (default: 2).

    Returns:
        (is_convergent, pattern, substitution, conserved_pairs)
    """
    fg_aas = "".join(normalize_residue(aa) for aa in fg_aas)
    bg_aas = "".join(normalize_residue(aa) for aa in bg_aas)
    substitution = f"{fg_aas}/{bg_aas}"
    if fg_species_list is None:
        fg_species_list = []
    if bg_species_list is None:
        bg_species_list = []

    fg_groups = [aa for aa in fg_aas if aa in STANDARD_AAS]
    bg_groups = [aa for aa in bg_aas if aa in STANDARD_AAS]

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
    overlap = residue_overlap(fg_aas, bg_aas)
    non_overlapping_fg = sum(1 for g in fg_groups if g not in bg_unique)
    non_overlapping_bg = sum(1 for g in bg_groups if g not in fg_unique)

    # Strict (overlap == 0) vs relaxed (overlap <= max_conserved) check
    if max_conserved == 0:
        is_convergent = (overlap == 0 and (non_overlapping_fg >= min_non_overlapping or non_overlapping_bg >= min_non_overlapping))
    else:
        is_convergent = (overlap <= max_conserved and (non_overlapping_fg >= min_non_overlapping or non_overlapping_bg >= min_non_overlapping))

    conserved_pairs = "0:"
    if is_convergent and max_conserved > 0 and multiconfig:
        # The two sides can lose different members to gaps/missingness.
        # Compare by domain ID, never by independently shortened list offsets.
        fg_by_pair = {multiconfig.get_pair(sp, trait): aa
                      for sp, aa in zip(fg_species_list, fg_aas) if aa in STANDARD_AAS}
        bg_by_pair = {multiconfig.get_pair(sp, trait): aa
                      for sp, aa in zip(bg_species_list, bg_aas) if aa in STANDARD_AAS}
        conserved_pair_indices = sorted(
            (str(pid) for pid in fg_by_pair.keys() & bg_by_pair.keys()
             if pid is not None and fg_by_pair[pid] == bg_by_pair[pid]), key=int)
        conserved_pairs = f"{len(conserved_pair_indices)}:{','.join(conserved_pair_indices)}"

    return (is_convergent, pattern, substitution, conserved_pairs)
