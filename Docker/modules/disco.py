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

MODULE NAME: disco.py
DESCRIPTION: runs the caas discovery on one single alignment. Returns non-validated caas candidate positions.
DEPENDENCIES: alimport.py, caas_id.py, pindex.py
CALLED BY: CT.

'''


from modules.caas_id import *
from modules.alimport import *
from modules.pindex import *
import os
from os.path import exists

### FUNCTION discovery()
### Scans one single alignment to identify the CAAS or CAAP

def discovery(input_cfg, sliced_object, max_fg_gaps, max_bg_gaps, max_overall_gaps, max_fg_miss, max_bg_miss, max_overall_miss, admitted_patterns, output_file, paired_mode=False, miss_pair=False, max_conserved=0, caap_mode=False):

    # Step 1: import the trait into a trait object (load_cfg from pindex.py)
    trait_object = load_cfg(input_cfg, paired_mode=paired_mode)

    # Step 2: import the alignment int a processed position object (slice from alimport.py)
    p = sliced_object

    # Step 3: processes the positions from imported alignment (process_position() from caas_id.py)
    processed_positions = map(functools.partial(process_position, multiconfig = trait_object, species_in_alignment = p.species), p.d)

    # Step 4: Collect results before writing to file
    results_to_write = []

    # Step 5: extract the raw caas or caap
    if caap_mode:
        # Import caap_id module for CAAP mode
        from modules.caap_id import fetch_caap
        
        # CAAP mode: detect property-based convergence
        for position in processed_positions:
            caap_results = fetch_caap( genename = p.genename,
                        position_obj = position,
                        trait_list = trait_object.alltraits,
                        
                        max_fg_gaps = int(max_fg_gaps) if max_fg_gaps != "NO" else 999999,
                        max_bg_gaps = int(max_bg_gaps) if max_bg_gaps != "NO" else 999999,
                        max_overall_gaps = int(max_overall_gaps) if max_overall_gaps != "NO" else 999999,
                        
                        max_fg_miss = int(max_fg_miss) if max_fg_miss != "NO" else 999999,
                        max_bg_miss = int(max_bg_miss) if max_bg_miss != "NO" else 999999,
                        max_overall_miss = int(max_overall_miss) if max_overall_miss != "NO" else 999999,
                        
                        output_file = None,  # Don't write yet
                        paired_mode = paired_mode,
                        miss_pair = miss_pair,
                        max_conserved = max_conserved,
                        species_in_alignment = p.species,
                        allowed_patterns = admitted_patterns,
                        multiconfig = trait_object,
                        return_results = True  # Get results instead of writing
                        )
            if caap_results:
                results_to_write.extend(caap_results)
    else:
        # CAAS mode: classical detection
        for position in processed_positions:
            caas_results = fetch_caas( p.genename,
                        position,
                        trait_object.alltraits,

                        maxgaps_bg= max_bg_gaps,
                        maxgaps_fg= max_fg_gaps,
                        maxgaps_all= max_overall_gaps,

                        maxmiss_bg= max_bg_miss,
                        maxmiss_fg= max_fg_miss,
                        maxmiss_all= max_overall_miss,
                        
                        multiconfig= trait_object,
                        miss_pair= miss_pair,
                        max_conserved= max_conserved,

                        admitted_patterns=admitted_patterns,
                        output_file = None,  # Don't write yet
                        return_results = True  # Get results instead of writing
                        )
            if caas_results:
                results_to_write.extend(caas_results)
    
    # Step 6: Only write output file if CAAS/CAAP were found
    if len(results_to_write) > 0:
        # Delete existing file if present
        if exists(output_file):
            os.system("rm -r " + output_file)
        
        # Determine header based on mode
        if caap_mode:
            header_fields = [
                "Gene",
                "Mode",
                "CAAP_Group",
                "Trait",
                "Position",
                "Substitution",
                "Encoded",
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
        else:
            header_fields = [
                "Gene",
                "Mode",
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
        if paired_mode and max_conserved > 0:
            header_fields.extend(["ConservedPair", "ConservedPairs"])
        
        header = "\t".join(header_fields)
        
        # Write header and results
        with open(output_file, "w") as outf:
            outf.write(header + "\n")
            for result_line in results_to_write:
                outf.write(result_line + "\n")
        
        print(f"Discovery complete: {len(results_to_write)} CAAS/CAAP found in {p.genename}")
    else:
        print(f"Discovery complete: No CAAS/CAAP found in {p.genename} - output file not created")