#!/usr/bin/env python3
"""The scalar CAAStools discovery of one alignment: a test-only reference, not part of the pipeline.

The pipeline finds the observed CAAS with the vectorized kernel of `ct perm-replay` (the b_0 labeling of the permulation
core). This is the position-by-position implementation that kernel is checked against (test_discovery_kernel.py,
test_kernel_fuzz.py, test_b0_discovery.py, check_b0_background.py): same options as the former `ct discovery`, same
discovery.tab columns, same background file.

    python3 reference_discovery.py -a G.fasta -t traits.tab -o discovery.tab --fmt fasta [--background_output bg] [filters]
"""
import functools
import os
import sys
from optparse import OptionParser
from os.path import exists
from pathlib import Path

CT_LOCAL = Path(os.environ.get("PHYLOPHERE_ROOT", Path(__file__).resolve().parents[2])) / "subworkflows/CT/local"
sys.path.insert(0, str(CT_LOCAL))

from modules.caas_id import *  # noqa: E402,F401,F403
from modules.alimport import *  # noqa: E402,F401,F403
from modules.pindex import *  # noqa: E402,F401,F403
from modules.runslice import runslice  # noqa: E402


### FUNCTION discovery()
### Scans one single alignment to identify the CAAS or CAAP

def discovery(input_cfg, sliced_object, max_fg_gaps, max_bg_gaps, max_overall_gaps, max_fg_miss, max_bg_miss, max_overall_miss, admitted_patterns, output_file, miss_pair=False, max_conserved=0, caap_mode=False, background_output_file=None):

    def _valid_traits_for_position(processed_position, trait_list, multiconfig,
                                   max_fg_gaps, max_bg_gaps, max_overall_gaps,
                                   max_fg_miss, max_bg_miss, max_overall_miss,
                                   miss_pair):
        valid_traits = []

        for trait in trait_list:
            if trait not in processed_position.trait2aas_fg:
                continue
            if trait not in processed_position.trait2aas_bg:
                continue

            # Gap filtering
            if max_fg_gaps != "NO" and processed_position.trait2gaps_fg.get(trait, 0) > int(max_fg_gaps):
                continue
            if max_bg_gaps != "NO" and processed_position.trait2gaps_bg.get(trait, 0) > int(max_bg_gaps):
                continue
            if max_overall_gaps != "NO" and processed_position.trait2gaps_fg.get(trait, 0) + processed_position.trait2gaps_bg.get(trait, 0) > int(max_overall_gaps):
                continue

            # Missing filtering
            if max_fg_miss != "NO" and processed_position.trait2miss_fg.get(trait, 0) > int(max_fg_miss):
                continue
            if max_bg_miss != "NO" and processed_position.trait2miss_bg.get(trait, 0) > int(max_bg_miss):
                continue
            if max_overall_miss != "NO" and processed_position.trait2miss_fg.get(trait, 0) + processed_position.trait2miss_bg.get(trait, 0) > int(max_overall_miss):
                continue

            # Pair-aware filtering
            if miss_pair:
                miss_thresholds_equal = False
                if max_fg_miss != "NO" and max_bg_miss != "NO" and max_fg_miss == max_bg_miss:
                    miss_thresholds_equal = True
                elif max_fg_miss == "NO" and max_bg_miss == "NO" and max_overall_miss != "NO":
                    miss_thresholds_equal = True

                if miss_thresholds_equal:
                    miss_pairs_fg = set(processed_position.trait2miss_pairs_fg.get(trait, []))
                    miss_pairs_bg = set(processed_position.trait2miss_pairs_bg.get(trait, []))
                    if miss_pairs_fg and miss_pairs_bg and miss_pairs_fg != miss_pairs_bg:
                        continue

                gap_thresholds_equal = False
                if max_fg_gaps != "NO" and max_bg_gaps != "NO" and max_fg_gaps == max_bg_gaps:
                    gap_thresholds_equal = True
                elif max_fg_gaps == "NO" and max_bg_gaps == "NO" and max_overall_gaps != "NO":
                    gap_thresholds_equal = True

                if gap_thresholds_equal:
                    gap_pairs_fg = set(processed_position.trait2gap_pairs_fg.get(trait, []))
                    gap_pairs_bg = set(processed_position.trait2gap_pairs_bg.get(trait, []))
                    if gap_pairs_fg and gap_pairs_bg and gap_pairs_fg != gap_pairs_bg:
                        continue

            valid_traits.append(trait)

        return valid_traits

    # Step 1: import the trait into a trait object (load_cfg from pindex.py)
    trait_object = load_cfg(input_cfg)

    # Step 2: import the alignment int a processed position object (slice from alimport.py)
    p = sliced_object

    # Step 3: processes the positions from imported alignment (process_position() from caas_id.py)
    processed_positions = map(functools.partial(process_position, multiconfig = trait_object, species_in_alignment = p.species), p.d)

    # Step 4: Collect results before writing to file
    results_to_write = []
    tested_positions = set()

    # Step 5: extract convergent mutations/properties across selected schemes
    # In caap_mode (default), all 5 schemes (US, GS1-GS4) are tested.
    # When caap_mode is disabled, only the ungrouped scheme (US) is tested.
    schemes_to_test = SCHEMES if caap_mode else {"US": US}

    for position in processed_positions:
        valid_traits = _valid_traits_for_position(
            position,
            trait_object.alltraits,
            trait_object,
            max_fg_gaps,
            max_bg_gaps,
            max_overall_gaps,
            max_fg_miss,
            max_bg_miss,
            max_overall_miss,
            miss_pair
        )
        if valid_traits:
            tested_positions.add(position.position)
        caas_results = fetch_caas(
            genename = p.genename,
            position_obj = position,
            trait_list = trait_object.alltraits,

            max_fg_gaps = int(max_fg_gaps) if max_fg_gaps != "NO" else 999999,
            max_bg_gaps = int(max_bg_gaps) if max_bg_gaps != "NO" else 999999,
            max_overall_gaps = int(max_overall_gaps) if max_overall_gaps != "NO" else 999999,

            max_fg_miss = int(max_fg_miss) if max_fg_miss != "NO" else 999999,
            max_bg_miss = int(max_bg_miss) if max_bg_miss != "NO" else 999999,
            max_overall_miss = int(max_overall_miss) if max_overall_miss != "NO" else 999999,

            output_file = None,
            miss_pair = miss_pair,
            max_conserved = max_conserved,
            species_in_alignment = p.species,
            admitted_patterns = admitted_patterns,
            multiconfig = trait_object,
            schemes = schemes_to_test,
            return_results = True
        )
        if caas_results:
            results_to_write.extend(caas_results)
    
    # Step 6: Write background coverage file (positions tested)
    if background_output_file:
        if tested_positions:
            positions_sorted = ",".join(map(str, sorted(tested_positions, key=lambda x: int(x))))
        else:
            positions_sorted = "NULL"
        with open(background_output_file, "w") as bkg_out:
            bkg_out.write(f"{p.genename}\t{positions_sorted}\n")

    # Step 7: Only write output file if CAAS/CAAP were found
    if len(results_to_write) > 0:
        # Delete existing file if present
        if exists(output_file):
            os.system("rm -r " + output_file)
        
        header_fields = [
            "gene",
            "mode",
            "caap_group",
            "trait",
            "position",
            "caas",
            "amino_encoded",
            "pattern",
            "ffgn",
            "fbgn",
            "gfg",
            "gbg",
            "mfg",
            "mbg",
            "ffg",
            "fbg",
            "ms"
        ]
        
        # Add conserved-pair columns when overlap tolerance is enabled
        if max_conserved > 0:
            header_fields.extend(["is_conserved_meta", "conserved_pair"])
        
        header = "\t".join(header_fields)
        
        # Write header and results
        with open(output_file, "w") as outf:
            outf.write(header + "\n")
            for result_line in results_to_write:
                outf.write(result_line + "\n")
        
        print(f"Discovery complete: {len(results_to_write)} CAAS/CAAP found in {p.genename}")
    else:
        print(f"Discovery complete: No CAAS/CAAP found in {p.genename} - output file not created")


def main(argv=None):
    parser = OptionParser()
    parser.add_option("-a", "--alignment", dest="single_alignment", default="none")
    parser.add_option("--fmt", dest="ali_format", default="clustal")
    parser.add_option("-t", "--traitfile", dest="config_file", default="none")
    parser.add_option("-o", "--output", dest="output_file", default="none")
    parser.add_option("--background_output", dest="background_output", default="background.output")
    parser.add_option("--patterns", dest="patterns_string", default="1,2,3")
    parser.add_option("--max_bg_gaps", dest="max_bg_gaps_string", default="NO")
    parser.add_option("--max_fg_gaps", dest="max_fg_gaps_string", default="NO")
    parser.add_option("--max_gaps", dest="max_gaps_string", default="NO")
    parser.add_option("--max_gaps_per_position", dest="max_gaps_pos_string", default="0.5")
    parser.add_option("--max_bg_miss", dest="max_bg_miss_string", default="NO")
    parser.add_option("--max_fg_miss", dest="max_fg_miss_string", default="NO")
    parser.add_option("--max_miss", dest="max_miss_string", default="NO")
    parser.add_option("--miss_pair", dest="miss_pair", action="store_true", default=False)
    parser.add_option("--max_conserved", dest="max_conserved", default="0")
    parser.add_option("--caap_mode", dest="caap_mode", action="store_true", default=False)
    options, _ = parser.parse_args(argv)
    if "none" in (options.single_alignment, options.config_file, options.output_file):
        parser.error("-a, -t and -o are required")

    discovery(
        input_cfg=options.config_file,
        sliced_object=runslice(options),
        max_fg_gaps=options.max_fg_gaps_string,
        max_bg_gaps=options.max_bg_gaps_string,
        max_overall_gaps=options.max_gaps_string,
        max_fg_miss=options.max_fg_miss_string,
        max_bg_miss=options.max_bg_miss_string,
        max_overall_miss=options.max_miss_string,
        miss_pair=options.miss_pair,
        max_conserved=int(options.max_conserved),
        caap_mode=options.caap_mode,
        admitted_patterns=options.patterns_string,
        output_file=options.output_file,
        background_output_file=options.background_output,
    )


if __name__ == "__main__":
    main()
