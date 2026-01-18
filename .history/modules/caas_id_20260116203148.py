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

MODULE NAME:    caas_id
DESCRIPTION:    Identification of caas mutation from MSA.
DEPENDENCIES:   pindex, alimport


TABLE OF CONTENTS
------------------------------------------

process_position()          processes a position from an imported alignment.
                            The output will be a dictionary that points each 
                            aminoacid (gaps included) to the
                            species sharing it.

check_pattern()            checks the pattern

fetch_caas()                fetches caas per each thing
'''                                                       

from modules.pindex import *
from modules.alimport import *
from modules.hyper import *
from os.path import exists

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
            
            # Pair-aware attributes (only populated in paired mode)
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
        try:
            z.aas2species[position[x].split("@")[0]].append(x)
        except:
            z.aas2species[position[x].split("@")[0]] = [x]

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
        
        # Track pairs for missing and gapped species (only in paired mode)
        if multiconfig.paired_mode:
            # Get pairs for missing species
            miss_pairs_fg = set()
            miss_pairs_bg = set()
            for sp in miss_fg:
                pair = multiconfig.get_pair(sp)
                if pair:
                    miss_pairs_fg.add(pair)
            for sp in miss_bg:
                pair = multiconfig.get_pair(sp)
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
                pair = multiconfig.get_pair(sp)
                if pair:
                    gap_pairs_fg.add(pair)
            for sp in gapped_bg:
                pair = multiconfig.get_pair(sp)
                if pair:
                    gap_pairs_bg.add(pair)
            
            z.trait2gap_pairs_fg[trait] = list(gap_pairs_fg)
            z.trait2gap_pairs_bg[trait] = list(gap_pairs_bg)

    return z


# FUNCTION check_pattern()
# checks the pattern

def iscaas(input_string, multiconfig=None, position_dict=None, max_conserved=0, trait=None):
    
    class caaspositive():
        def __init__(self):
            self.caas = True
            self.pattern = "4"
            self.conserved_pairs = "0:"
    
    z = caaspositive()

    twosides = input_string.split("/")
    # Convert to unique amino acids for validation
    fg = list(set(twosides[0]))
    bg = list(set(twosides[1]))

    # Standard CAAS check (no overlap allowed)
    if max_conserved == 0 or not multiconfig or not multiconfig.paired_mode:
        # Original logic: strict no-overlap
        for x in fg:
            if x in bg:
                z.caas = False
                break
    else:
        # Relaxed logic: allow overlap up to max_conserved pairs
        # Need to track which species contribute each amino acid
        
        # Build species-to-AA mapping from position_dict
        species2aa = {}
        if position_dict:
            for species, aa_info in position_dict.items():
                aa = aa_info.split("@")[0]
                species2aa[species] = aa
        
        # Find overlapping amino acids
        overlap_aas = set(fg).intersection(set(bg))
        
        if not overlap_aas:
            # No overlap, this is a standard CAAS
            z.caas = True
        else:
            # When there's overlap, we need to validate it strictly
            conserved_pair_set = set()
            all_overlap_valid = True
            
            # Get all FG and BG species for this trait
            if trait and multiconfig:
                fg_species_all = multiconfig.trait2fg.get(trait, [])
                bg_species_all = multiconfig.trait2bg.get(trait, [])
                
                # For each overlapping amino acid, ALL species with it must have paired species with same AA
                for aa in overlap_aas:
                    # Find which FG and BG species have this amino acid (and are present)
                    fg_with_aa = [sp for sp in fg_species_all if species2aa.get(sp) == aa]
                    bg_with_aa = [sp for sp in bg_species_all if species2aa.get(sp) == aa]
                    
                    # Check FG species with this overlapping AA
                    for fg_sp in fg_with_aa:
                        fg_pair_id = multiconfig.get_pair(fg_sp)
                        if fg_pair_id:
                            # Find the BG species with the same pair ID
                            paired_bg_sp = None
                            for bg_sp in bg_species_all:
                                if multiconfig.get_pair(bg_sp) == fg_pair_id:
                                    paired_bg_sp = bg_sp
                                    break
                            
                            if paired_bg_sp:
                                # Check if the paired BG species has the same AA
                                if species2aa.get(paired_bg_sp) == aa:
                                    conserved_pair_set.add(fg_pair_id)
                                else:
                                    # FG species has overlapping AA but its paired BG species doesn't - reject
                                    all_overlap_valid = False
                                    break
                    
                    if not all_overlap_valid:
                        break
                    
                    # Check BG species with this overlapping AA
                    for bg_sp in bg_with_aa:
                        bg_pair_id = multiconfig.get_pair(bg_sp)
                        if bg_pair_id:
                            # Find the FG species with the same pair ID
                            paired_fg_sp = None
                            for fg_sp in fg_species_all:
                                if multiconfig.get_pair(fg_sp) == bg_pair_id:
                                    paired_fg_sp = fg_sp
                                    break
                            
                            if paired_fg_sp:
                                # Check if the paired FG species has the same AA
                                if species2aa.get(paired_fg_sp) != aa:
                                    # BG species has overlapping AA but its paired FG species doesn't - reject
                                    all_overlap_valid = False
                                    break
                    
                    if not all_overlap_valid:
                        break
            
            if not all_overlap_valid:
                # Overlap exists but not all from valid conserved pairs
                z.caas = False
            else:
                # Count how many SPECIES have non-overlapping AAs (actual changes)
                # We need at least 2 species to have changed on one side for meaningful divergence
                
                # Count FG species with non-overlapping AAs
                fg_changed_count = sum(1 for sp in fg_species_all 
                                      if species2aa.get(sp) and species2aa.get(sp) not in overlap_aas)
                
                # Count BG species with non-overlapping AAs
                bg_changed_count = sum(1 for sp in bg_species_all 
                                      if species2aa.get(sp) and species2aa.get(sp) not in overlap_aas)
                
                # At least one side must have 2+ species that changed
                if fg_changed_count < 2 and bg_changed_count < 2:
                    # Not enough changes - need at least 2 species changed on one side
                    z.caas = False
                else:
                    num_conserved_pairs = len(conserved_pair_set)
                    if num_conserved_pairs <= max_conserved:
                        z.caas = True
                        conserved_pairs_list = sorted(list(conserved_pair_set))
                        z.conserved_pairs = f"{num_conserved_pairs}:{','.join(conserved_pairs_list)}"
                    else:
                        z.caas = False

    # What is the pattern?
    if len(fg) == 1 and len(bg) == 1:
        z.pattern = "1"
    
    elif len(fg) == 1:
        z.pattern = "2"
    elif len(bg) == 1:
        z.pattern = "3"
    
    if len(fg) == 0 or len(bg) == 0:
        z.pattern = "null"
    
    return z

# FUNCTION fetch_caas():
# fetches caas per each thing

def fetch_caas(genename, processed_position, list_of_traits, output_file, maxgaps_fg, maxgaps_bg, maxgaps_all, maxmiss_fg, maxmiss_bg, maxmiss_all, multiconfig, miss_pair=False, max_conserved=0, admitted_patterns = ["1","2","3"]):

    a = set(list_of_traits)
    b = set(processed_position.trait2aas_fg.keys())
    c = set(processed_position.trait2aas_bg.keys())

    valid_traits = list(a.intersection(b).intersection(c))

    # Filter for the number of gaps and missing species

    if len(valid_traits) > 0:

        for trait in valid_traits:

            ### GAPS filtering.

            if maxgaps_fg != "NO" and processed_position.trait2gaps_fg[trait] > int(maxgaps_fg):
                if len(valid_traits) > 0:
                    valid_traits.remove(trait)

            # FIX: Changed maxgaps_fg to maxgaps_bg
            if maxgaps_bg != "NO" and processed_position.trait2gaps_bg[trait] > int(maxgaps_bg):
                if len(valid_traits) > 0:
                    valid_traits.remove(trait)

            if maxgaps_all != "NO" and processed_position.trait2gaps_fg[trait] + processed_position.trait2gaps_bg[trait] > int(maxgaps_all):
                if len(valid_traits) > 0:
                    valid_traits.remove(trait)


            ### Missing filtering.

            if maxmiss_fg != "NO" and processed_position.trait2miss_fg[trait] > int(maxmiss_fg):
                if len(valid_traits) > 0:
                    valid_traits.remove(trait)

            # FIX: Changed maxmiss_fg to maxmiss_bg
            if maxmiss_bg != "NO" and processed_position.trait2miss_bg[trait] > int(maxmiss_bg):
                if len(valid_traits) > 0:
                    valid_traits.remove(trait)

            if maxmiss_all != "NO" and processed_position.trait2miss_fg[trait] + processed_position.trait2miss_bg[trait] > int(maxmiss_all):
                if len(valid_traits) > 0:
                    valid_traits.remove(trait)
            
            ### Pair-aware filtering (only in paired mode with miss_pair enabled)
            
            if multiconfig.paired_mode and miss_pair and trait in valid_traits:
                # Check if missing thresholds are equal (either both explicitly set, or both using overall)
                miss_thresholds_equal = False
                if maxmiss_fg != "NO" and maxmiss_bg != "NO" and maxmiss_fg == maxmiss_bg:
                    miss_thresholds_equal = True
                elif maxmiss_fg == "NO" and maxmiss_bg == "NO" and maxmiss_all != "NO":
                    # Both fg and bg use the overall threshold - treat as equal
                    miss_thresholds_equal = True
                
                if miss_thresholds_equal:
                    miss_pairs_fg = set(processed_position.trait2miss_pairs_fg.get(trait, []))
                    miss_pairs_bg = set(processed_position.trait2miss_pairs_bg.get(trait, []))
                    
                    # Check if missing pairs match
                    if miss_pairs_fg != miss_pairs_bg:
                        valid_traits.remove(trait)
                        continue
                
                # Check if gap thresholds are equal (either both explicitly set, or both using overall)
                gap_thresholds_equal = False
                if maxgaps_fg != "NO" and maxgaps_bg != "NO" and maxgaps_fg == maxgaps_bg:
                    gap_thresholds_equal = True
                elif maxgaps_fg == "NO" and maxgaps_bg == "NO" and maxgaps_all != "NO":
                    # Both fg and bg use the overall threshold - treat as equal
                    gap_thresholds_equal = True
                
                if gap_thresholds_equal:
                    gap_pairs_fg = set(processed_position.trait2gap_pairs_fg.get(trait, []))
                    gap_pairs_bg = set(processed_position.trait2gap_pairs_bg.get(trait, []))
                    
                    # Check if gapped pairs match
                    if gap_pairs_fg != gap_pairs_bg:
                        if trait in valid_traits:
                            valid_traits.remove(trait)
                        continue

     
        output_traits = []

        # Filter for pattern
        for x in valid_traits:
            # Get foreground and background species lists (already filtered for ungapped)
            fg_species = processed_position.trait2ungapped_fg[x][:]
            fg_species.sort()
            bg_species = processed_position.trait2ungapped_bg[x][:]
            bg_species.sort()

            # Extract amino acid for each species to create expanded pattern
            aa_tag_fg = "".join([processed_position.d[sp].split("@")[0] for sp in fg_species])
            aa_tag_bg = "".join([processed_position.d[sp].split("@")[0] for sp in bg_species])

            tag = "/".join([aa_tag_fg, aa_tag_bg])

            check = iscaas(tag, multiconfig, processed_position.d, max_conserved, x)

            if check.caas == True and check.pattern in admitted_patterns:
                # Store pattern info including conserved pair information
                trait_info = x + "@" + tag + "@pattern" + check.pattern
                if hasattr(check, 'conserved_pairs'):
                    trait_info += "@conserved:" + str(check.conserved_pairs)
                output_traits.append(trait_info)
        

        # Print the output

        header_fields = [
            "Gene",
            "Trait",
            "Position",
            "Substitution",
            "Pvalue",
            "Pattern",
            "FFGN",
            "FBGN",
            "GFG",
            "GBG",
            "MFG",
            "MBG",
            "FFG",
            "FBG",
            "MS"
        ]
        
        # Add conserved pair columns in paired mode
        if multiconfig.paired_mode and max_conserved > 0:
            header_fields.extend(["ConservedPair", "ConservedPairs"])
        
        header = "\t".join(header_fields)

        if len(output_traits) > 0:

            if exists(output_file):
                out = open(output_file, "a")
            else:
                out = open(output_file, "w")
                print(header, file=out)

            for trait in output_traits:

                traitname = trait.split("@")[0]            
                change = trait.split("@")[1]
                thepattern = trait.split("@")[2]
                
                # Extract conserved pair info if present
                conserved_info = None
                if len(trait.split("@")) > 3 and trait.split("@")[3].startswith("conserved:"):
                    conserved_info = trait.split("@")[3].replace("conserved:", "")

                fg_species_number = str(len(processed_position.trait2ungapped_fg[traitname]))
                bg_species_number = str(len(processed_position.trait2ungapped_bg[traitname]))

                fg_ungapped = processed_position.trait2ungapped_fg[traitname]
                bg_ungapped = processed_position.trait2ungapped_bg[traitname]

                fg_ungapped.sort()
                bg_ungapped.sort()

                missings = "-"

                if len(processed_position.trait2missings[traitname]) > 0:
                    missings = ",".join(processed_position.trait2missings[traitname])
                
                # Starting the pvalue determination


                pv = calcpval_random(processed_position.d, genename, int(fg_species_number), int(bg_species_number))
                pvalue_string = str(pv)


                print("CAAS found in alignment", genename, "on position", processed_position.position, "with pvalue", pvalue_string)

                # End of the pvalue determination

                output_fields = [
                    genename,
                    traitname,
                    processed_position.position,
                    change,
                    pvalue_string,
                    thepattern,
                    fg_species_number,
                    bg_species_number,
                    str(processed_position.trait2gaps_fg[traitname]),
                    str(processed_position.trait2gaps_bg[traitname]),
                    str(processed_position.trait2miss_fg[traitname]),
                    str(processed_position.trait2miss_bg[traitname]),
                    ",".join(fg_ungapped),
                    ",".join(bg_ungapped),
                    missings
                ]
                
                # Add conserved pair columns in paired mode
                if multiconfig.paired_mode and max_conserved > 0:
                    if conserved_info:
                        conserved_count = conserved_info.split(":")[0]
                        conserved_pairs_list = conserved_info.split(":")[1] if ":" in conserved_info else ""
                        output_fields.extend([
                            "TRUE" if int(conserved_count) > 0 else "FALSE",
                            f"{conserved_count}:{conserved_pairs_list}"
                        ])
                    else:
                        output_fields.extend(["FALSE", "0:"])
                
                #pvalue_string = pvdict[genename + "@" + processed_position.position]
                print("\t".join(output_fields), file = out)
            

            out.close()