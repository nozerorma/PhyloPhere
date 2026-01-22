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

MODULE NAME: boot.py
DESCRIPTION: bootstrap function
DEPENDENCIES: alimport.py, caas_id.py, pindex.py
CALLED BY: ct

'''


from modules.init_bootstrap import *
from modules.disco import process_position
from modules.caas_id import iscaas
from modules.caap_id import check_caap_pattern, encode_to_groups, US, GS0, GS1, GS2, GS3, GS4
from modules.alimport import *

from os.path import exists
import functools
import time
from datetime import datetime


# UTILITY FUNCTIONS for progress tracking

def format_time(seconds):
    """Format seconds into human-readable time string.
    
    Args:
        seconds: Time in seconds
        
    Returns:
        str: Formatted time (e.g., "2h 15m", "45m 30s", "15s")
    """
    if seconds < 60:
        return f"{int(seconds)}s"
    elif seconds < 3600:
        mins = int(seconds / 60)
        secs = int(seconds % 60)
        return f"{mins}m {secs}s"
    else:
        hours = int(seconds / 3600)
        mins = int((seconds % 3600) / 60)
        return f"{hours}h {mins}m"


def calculate_eta(processed, total, elapsed):
    """Calculate estimated time to completion.
    
    Args:
        processed: Number of items processed
        total: Total number of items
        elapsed: Elapsed time in seconds
        
    Returns:
        float: Estimated seconds remaining
    """
    if processed == 0:
        return 0
    rate = processed / elapsed
    remaining = total - processed
    return remaining / rate


def log_progress(current, total, start_time, log_file=None, prefix="Progress"):
    """Log progress with timestamp and ETA.
    
    Args:
        current: Current item number (1-based)
        total: Total items
        start_time: Start time from time.time()
        log_file: Optional file path for logging (None = stdout only)
        prefix: Prefix for log message
        
    Returns:
        str: Formatted progress message
    """
    elapsed = time.time() - start_time
    pct = (current / total) * 100
    eta_seconds = calculate_eta(current, total, elapsed)
    
    timestamp = datetime.now().strftime("%H:%M:%S")
    message = f"[{timestamp}] {prefix}: {current}/{total} ({pct:.1f}%) | Elapsed: {format_time(elapsed)} | ETA: {format_time(eta_seconds)}"
    
    print(message)
    
    if log_file:
        try:
            with open(log_file, 'a') as f:
                f.write(message + "\n")
        except:
            pass
    
    return message


# FUNCTION parse_discovery_positions()
# Parses discovery output file and extracts CAAS position numbers

def parse_discovery_positions(discovery_file, genename):
    """
    Parse discovery output file and extract positions with their CAAP grouping schemes.
    Handles both classical CAAS and CAAP formats.
    
    Classical CAAS format: Gene\tMode\tTrait\tPosition\t...
    CAAP format:          Gene\tMode\tCAAP_Group\tTrait\tPosition\t...
    
    Args:
        discovery_file: Path to discovery output file
        genename: Gene name to match (e.g., "BRCA1")
    
    Returns:
        dict mapping position -> set of grouping schemes found (e.g., {"142": {"GS0", "GS1"}})
        or None if no discovery file or all positions should be tested
    """
    if not exists(discovery_file):
        print(f"Warning: Discovery file {discovery_file} not found. Processing all positions with all schemes.")
        return None
    
    # position_number -> set of grouping schemes
    position_schemes = {}
    
    try:
        with open(discovery_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("Gene"):  # Skip header
                    continue
                
                try:
                    fields = line.split('\t')
                    if len(fields) < 3:
                        continue
                    
                    gene = fields[0]
                    
                    # Match gene name first
                    if gene != genename:
                        continue
                    
                    # Detect format based on mode column
                    mode = fields[1] if len(fields) > 1 else ""
                    
                    if mode == "CAAP":
                        # CAAP format: Gene\tMode\tCAAP_Group\tTrait\tPosition\t...
                        # CAAP_Group is at index 2, Position is at index 4
                        if len(fields) >= 5:
                            scheme = fields[2]  # US, GS0, GS1, GS2, GS3, or GS4
                            position = fields[4]
                            
                            if position not in position_schemes:
                                position_schemes[position] = set()
                            position_schemes[position].add(scheme)
                    elif mode == "CAAS":
                        # Classical CAAS format: Gene\tMode\tTrait\tPosition\t...
                        # Position is at index 3, no grouping scheme (test all)
                        if len(fields) >= 4:
                            position = fields[3]
                            # For classical CAAS, mark with special 'CAAS' marker to test all schemes
                            if position not in position_schemes:
                                position_schemes[position] = set()
                            position_schemes[position].add("CAAS")
                except:
                    continue
        
        return position_schemes if len(position_schemes) > 0 else None
    
    except Exception as e:
        print(f"Warning: Error parsing discovery file: {e}. Processing all positions with all schemes.")
        return None


# FUNCTION filter_for_gaps()
# filters a trait for its gaps

def filter_for_gaps(max_bg, max_fg, max_all, gfg, gbg):

    out = True

    all_g = gfg + gbg
    
    if max_all != "NO" and all_g > int(max_all):
        out = False

    elif max_fg != "NO" and gfg > int(max_fg):
        out = False

    elif max_bg != "NO" and gbg > int(max_bg):
        out = False

    return out


# FUNCTION filter_for_missing()
# filters a trait for missing species

def filter_for_missings(max_m_bg, max_m_fg, max_m_all, mfg, mbg):

    out = True

    all_m = mfg + mbg
    
    if max_m_all != "NO" and all_m > int(max_m_all):
        out = False

    elif max_m_fg != "NO" and mfg > int(max_m_fg):
        out = False

    elif max_m_bg != "NO" and mbg > int(max_m_bg):
        out = False

    return out



def caasboot(processed_position, genename, list_of_traits, maxgaps_fg, maxgaps_bg, maxgaps_all, maxmiss_fg, maxmiss_bg, maxmiss_all, cycles, multiconfig, miss_pair=False, max_conserved=0, admitted_patterns=["1","2","3"], chunk_size=1000, caap_mode=False, discovery_schemes=None, debug_rejects=False, groups_out=None, perm_discovery_out=None):
    """Chunked bootstrap - processes traits in batches to handle large resample files
    
    Args:
        caap_mode: If True, test CAAP grouping schemes instead of classical CAAS
        discovery_schemes: Set of grouping schemes found in discovery for this position (US, GS0-GS4, or CAAS)
                          If None, test all schemes. If provided, only test those schemes.
    """
    
    a = set(list_of_traits)
    b = set(processed_position.trait2aas_fg.keys())
    c = set(processed_position.trait2aas_bg.keys())
    valid_traits = list(a.intersection(b).intersection(c))
    
    if len(valid_traits) == 0:
        if debug_rejects:
            print(f"[BOOTSTRAP DEBUG] {genename}@{processed_position.position} rejected: no valid traits after presence check (fg keys={len(processed_position.trait2aas_fg)}, bg keys={len(processed_position.trait2aas_bg)})")
        position_name = genename + "@" + str(processed_position.position)
        if caap_mode:
            # Return one line per scheme with zero counts
            schemes = [("US", US), ("GS0", GS0), ("GS1", GS1), ("GS2", GS2), ("GS3", GS3), ("GS4", GS4)]
            outlines = []
            for scheme_name, _ in schemes:
                outline = "\t".join([position_name, scheme_name, "0", str(cycles), "0.0"])
                outlines.append(outline)
            return "\n".join(outlines)
        else:
            outline = "\t".join([position_name, "0", str(cycles), "0.0"])
            return outline
    
    # Process traits in chunks to avoid memory issues with large files
    total_output_traits = []
    n_traits = len(valid_traits)
    n_chunks = (n_traits + chunk_size - 1) // chunk_size  # Ceiling division
    
    for chunk_idx in range(n_chunks):
        start_idx = chunk_idx * chunk_size
        end_idx = min(start_idx + chunk_size, n_traits)
        chunk_traits = valid_traits[start_idx:end_idx]
        
        # Filter for gaps
        filtered_traits = []
        for trait in chunk_traits:
            gfg = processed_position.trait2gaps_fg[trait]
            gbg = processed_position.trait2gaps_bg[trait]
            
            if not filter_for_gaps(maxgaps_bg, maxgaps_fg, maxgaps_all, gfg, gbg):
                if debug_rejects:
                    print(f"[BOOTSTRAP DEBUG] {genename}@{processed_position.position} trait {trait} rejected: gaps gfg={gfg} gbg={gbg} max_fg={maxgaps_fg} max_bg={maxgaps_bg} max_all={maxgaps_all}")
                continue
            
            # Filter for missings
            mfg = processed_position.trait2miss_fg[trait]
            mbg = processed_position.trait2miss_bg[trait]
            
            if not filter_for_missings(maxmiss_bg, maxmiss_fg, maxmiss_all, mfg, mbg):
                if debug_rejects:
                    print(f"[BOOTSTRAP DEBUG] {genename}@{processed_position.position} trait {trait} rejected: missings mfg={mfg} mbg={mbg} max_fg={maxmiss_fg} max_bg={maxmiss_bg} max_all={maxmiss_all}")
                continue
            
            filtered_traits.append(trait)
        
        if debug_rejects and len(filtered_traits) == 0:
            print(f"[BOOTSTRAP DEBUG] {genename}@{processed_position.position} rejected: all traits filtered out by gaps/missings")
        
        # Pattern check
        if caap_mode:
            # CAAP mode: test grouping schemes
            # If discovery_schemes provided, only test those schemes for this position
            # Otherwise test all schemes
            scheme_counts = {"US": [], "GS0": [], "GS1": [], "GS2": [], "GS3": [], "GS4": []}
            
            # Determine which schemes to test
            if discovery_schemes:
                # Only test schemes found in discovery for this position
                schemes_to_test = []
                
                if "CAAS" in discovery_schemes:
                    # Classical CAAS mode - test all schemes
                    schemes_to_test = [("US", US), ("GS0", GS0), ("GS1", GS1), ("GS2", GS2), ("GS3", GS3), ("GS4", GS4)]
                else:
                    # CAAP mode - only test discovered schemes
                    scheme_map = {"US": US, "GS0": GS0, "GS1": GS1, "GS2": GS2, "GS3": GS3, "GS4": GS4}
                    for scheme_name in discovery_schemes:
                        if scheme_name in scheme_map:
                            schemes_to_test.append((scheme_name, scheme_map[scheme_name]))
            else:
                # No discovery info - test all schemes
                schemes_to_test = [("US", US), ("GS0", GS0), ("GS1", GS1), ("GS2", GS2), ("GS3", GS3), ("GS4", GS4)]
            
            for trait in filtered_traits:
                fg_species = processed_position.trait2ungapped_fg[trait][:]
                bg_species = processed_position.trait2ungapped_bg[trait][:]
                
                # Sort species by pair number if in paired mode, otherwise alphabetically
                if multiconfig.paired_mode:
                    def pair_sort_key(sp):
                        pair_id = multiconfig.get_pair(sp)
                        if pair_id:
                            try:
                                return (int(pair_id), sp)
                            except (ValueError, TypeError):
                                return (float('inf'), sp)
                        return (float('inf'), sp)
                    
                    fg_species.sort(key=pair_sort_key)
                    bg_species.sort(key=pair_sort_key)
                else:
                    fg_species.sort()
                    bg_species.sort()
                
                # Extract amino acids per species
                fg_aas = [processed_position.d[sp].split("@")[0] for sp in fg_species]
                bg_aas = [processed_position.d[sp].split("@")[0] for sp in bg_species]

                # Test only the schemes found in discovery (or all if no discovery)
                any_match = False
                pattern_by_scheme = {}
                debug_by_scheme = {}
                for scheme_name, scheme_dict in schemes_to_test:
                    is_caap, pattern, substitution, conserved_pairs = check_caap_pattern(
                        fg_aas, bg_aas, scheme_dict, max_conserved,
                        multiconfig, fg_species, bg_species
                    )
                    pattern_by_scheme[scheme_name] = pattern
                    if is_caap and pattern in admitted_patterns:
                        any_match = True
                        scheme_counts[scheme_name].append(trait)
                        if groups_out:
                            groups_out.write(f"{trait}\t{genename}\t{processed_position.position}\tCAAP\t{scheme_name}\n")
                        if perm_discovery_out:
                            fg_species_str = ",".join(fg_species) if fg_species else "NA"
                            bg_species_str = ",".join(bg_species) if bg_species else "NA"
                            miss_species = processed_position.trait2missings.get(trait, [])
                            miss_species_str = ",".join(miss_species) if miss_species else "NA"
                            fg_aas_str = "".join(fg_aas)
                            bg_aas_str = "".join(bg_aas)
                            encoded = encode_to_groups(fg_aas_str, scheme_dict) + "/" + encode_to_groups(bg_aas_str, scheme_dict)
                            output_fields = [
                                trait,
                                genename,
                                "CAAP",
                                scheme_name,
                                trait,
                                str(processed_position.position),
                                substitution,
                                encoded,
                                "NA",
                                pattern,
                                str(len(fg_species)),
                                str(len(bg_species)),
                                str(processed_position.trait2gaps_fg.get(trait, 0)),
                                str(processed_position.trait2gaps_bg.get(trait, 0)),
                                str(processed_position.trait2miss_fg.get(trait, 0)),
                                str(processed_position.trait2miss_bg.get(trait, 0)),
                                fg_species_str,
                                bg_species_str,
                                miss_species_str,
                            ]
                            if multiconfig.paired_mode and max_conserved > 0:
                                overlap_count = conserved_pairs.split(":")[0] if conserved_pairs else "0"
                                pair_list = conserved_pairs.split(":")[1] if conserved_pairs and ":" in conserved_pairs else ""
                                output_fields.extend([
                                    "TRUE" if int(overlap_count) > 0 else "FALSE",
                                    f"{overlap_count}:{pair_list}"
                                ])
                            perm_discovery_out.write("\t".join(output_fields) + "\n")
                    if debug_rejects:
                        # Mirror CAAP overlap/change logic for debug
                        standard_aas = set("ACDEFGHIKLMNPQRSTVWY")
                        fg_filtered = [aa for aa in fg_aas if aa in standard_aas]
                        bg_filtered = [aa for aa in bg_aas if aa in standard_aas]
                        fg_groups = [scheme_dict.get(aa) for aa in fg_filtered if scheme_dict.get(aa) is not None]
                        bg_groups = [scheme_dict.get(aa) for aa in bg_filtered if scheme_dict.get(aa) is not None]
                        fg_unique = set(fg_groups)
                        bg_unique = set(bg_groups)
                        shared_types = fg_unique.intersection(bg_unique)
                        shared_fg = sum(1 for g in fg_groups if g in shared_types)
                        shared_bg = sum(1 for g in bg_groups if g in shared_types)
                        overlap = min(shared_fg, shared_bg)
                        non_overlapping_fg = sum(1 for g in fg_groups if g not in bg_unique)
                        non_overlapping_bg = sum(1 for g in bg_groups if g not in fg_unique)
                        debug_by_scheme[scheme_name] = {
                            "pattern": pattern,
                            "overlap": overlap,
                            "non_fg": non_overlapping_fg,
                            "non_bg": non_overlapping_bg,
                            "fg_groups": "".join(fg_groups),
                            "bg_groups": "".join(bg_groups),
                            "is_caap": is_caap,
                        }
                if debug_rejects and not any_match:
                    print(f"[BOOTSTRAP DEBUG] {genename}@{processed_position.position} trait {trait} rejected: no CAAP match (patterns={pattern_by_scheme})")
                    for scheme_name in sorted(debug_by_scheme.keys()):
                        info = debug_by_scheme[scheme_name]
                        print(
                            f"[BOOTSTRAP DEBUG] {genename}@{processed_position.position} {scheme_name} "
                            f"fg={info['fg_groups']} bg={info['bg_groups']} "
                            f"overlap={info['overlap']} non_fg={info['non_fg']} non_bg={info['non_bg']} "
                            f"is_caap={info['is_caap']}"
                        )
            
            # Return one line per tested scheme
            position_name = genename + "@" + str(processed_position.position)
            outlines = []
            scheme_set = {name for name, _ in schemes_to_test}
            ordered_schemes = [name for name in ["US", "GS0", "GS1", "GS2", "GS3", "GS4"] if name in scheme_set]
            for scheme_name in ordered_schemes:
                count = str(len(scheme_counts[scheme_name]))
                empval = str(len(scheme_counts[scheme_name])/cycles)
                outline = "\t".join([position_name, scheme_name, count, str(cycles), empval])
                outlines.append(outline)
            
            return "\n".join(outlines)
        else:
            # Classical CAAS mode
            for trait in filtered_traits:
                fg_species = processed_position.trait2ungapped_fg[trait][:]
                bg_species = processed_position.trait2ungapped_bg[trait][:]
                
                # Sort species by pair number if in paired mode, otherwise alphabetically
                if multiconfig.paired_mode:
                    def pair_sort_key(sp):
                        pair_id = multiconfig.get_pair(sp)
                        if pair_id:
                            try:
                                return (int(pair_id), sp)
                            except (ValueError, TypeError):
                                return (float('inf'), sp)
                        return (float('inf'), sp)
                    
                    fg_species.sort(key=pair_sort_key)
                    bg_species.sort(key=pair_sort_key)
                else:
                    fg_species.sort()
                    bg_species.sort()
                
                aa_tag_fg = "".join([processed_position.d[sp].split("@")[0] for sp in fg_species])
                aa_tag_bg = "".join([processed_position.d[sp].split("@")[0] for sp in bg_species])
                tag = "/".join([aa_tag_fg, aa_tag_bg])
                
                check = iscaas(tag, multiconfig, processed_position.d, max_conserved, trait, fg_species, bg_species)
                if check.caas == True and check.pattern in admitted_patterns:
                    total_output_traits.append(trait)
                    if groups_out:
                        groups_out.write(f"{trait}\t{genename}\t{processed_position.position}\tCAAS\tCAAS\n")
                    if perm_discovery_out:
                        fg_species_number = str(len(processed_position.trait2ungapped_fg[trait]))
                        bg_species_number = str(len(processed_position.trait2ungapped_bg[trait]))
                        fg_ungapped = processed_position.trait2ungapped_fg[trait][:]
                        bg_ungapped = processed_position.trait2ungapped_bg[trait][:]
                        if multiconfig.paired_mode:
                            def pair_sort_key(sp):
                                pair_id = multiconfig.get_pair(sp)
                                if pair_id:
                                    try:
                                        return (int(pair_id), sp)
                                    except (ValueError, TypeError):
                                        return (float('inf'), sp)
                                return (float('inf'), sp)
                            fg_ungapped.sort(key=pair_sort_key)
                            bg_ungapped.sort(key=pair_sort_key)
                        else:
                            fg_ungapped.sort()
                            bg_ungapped.sort()
                        missings = "-"
                        if len(processed_position.trait2missings.get(trait, [])) > 0:
                            missings = ",".join(processed_position.trait2missings[trait])
                        output_fields = [
                            trait,
                            genename,
                            "CAAS",
                            trait,
                            str(processed_position.position),
                            tag,
                            "NA",
                            f"pattern{check.pattern}",
                            fg_species_number,
                            bg_species_number,
                            str(processed_position.trait2gaps_fg.get(trait, 0)),
                            str(processed_position.trait2gaps_bg.get(trait, 0)),
                            str(processed_position.trait2miss_fg.get(trait, 0)),
                            str(processed_position.trait2miss_bg.get(trait, 0)),
                            ",".join(fg_ungapped),
                            ",".join(bg_ungapped),
                            missings
                        ]
                        if multiconfig.paired_mode and max_conserved > 0:
                            conserved_info = getattr(check, "conserved_pairs", "0:")
                            conserved_count = conserved_info.split(":")[0] if conserved_info else "0"
                            conserved_pairs_list = conserved_info.split(":")[1] if ":" in conserved_info else ""
                            output_fields.extend([
                                "TRUE" if int(conserved_count) > 0 else "FALSE",
                                f"{conserved_count}:{conserved_pairs_list}"
                            ])
                        perm_discovery_out.write("\t".join(output_fields) + "\n")
                elif debug_rejects:
                    print(f"[BOOTSTRAP DEBUG] {genename}@{processed_position.position} trait {trait} rejected: caas={check.caas} pattern={check.pattern} admitted={admitted_patterns}")
            
            # Return aggregated result
            position_name = genename + "@" + str(processed_position.position)
            count = str(len(total_output_traits))
            empval = str(int(count)/cycles)
            outline = "\t".join([position_name, count, str(cycles), empval])
            
            return outline

# FUNCTION boot_on_single_alignment()
# Launches the bootstrap in several lines. Returns a dictionary gene@position --> pvalue

def boot_on_single_alignment(trait_config_file, resampled_traits, sliced_object, max_fg_gaps, max_bg_gaps, max_overall_gaps, max_fg_miss, max_bg_miss, max_overall_miss, the_admitted_patterns, output_file, miss_pair=False, max_conserved=0, discovery_file=None, progress_log=None, caap_mode=False, export_groups=None, export_perm_discovery=None):
    """
    Run bootstrap analysis on a single alignment.
    
    Supports both single-file and directory-based resampled traits:
    - Single file: resampled_traits is a multicfg object loaded from one file
    - Directory: resampled_traits is the directory path (string), files loaded sequentially
    
    Args:
        resampled_traits: multicfg object OR directory path (str) containing resample_*.tab files
        progress_log: Optional file path for logging progress
        caap_mode: If True, test all CAAP grouping schemes (US, GS0-GS4) instead of classical CAAS
        ... (other parameters as before)
    """
    the_genename = sliced_object.genename
    
    groups_handle = None
    if export_groups:
        groups_handle = open(export_groups, "w")
        groups_handle.write("Cycle\tGene\tPosition\tMode\tGroup\n")
    
    perm_discovery_handle = None
    if export_perm_discovery:
        perm_discovery_handle = open(export_perm_discovery, "w")
        if caap_mode:
            header_fields = [
                "Cycle",
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
                "Cycle",
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
        if miss_pair and max_conserved > 0:
            header_fields.extend(["ConservedPair", "ConservedPairs"])
        perm_discovery_handle.write("\t".join(header_fields) + "\n")

    try:
        # Detect if resampled_traits is a directory path or a multicfg object
        if isinstance(resampled_traits, str) and os.path.isdir(resampled_traits):
            # Directory mode: sequential processing
            print(f"\n{'='*80}")
            print(f"DIRECTORY-BASED BOOTSTRAP MODE")
            print(f"{'='*80}\n")
            
            resample_dir = resampled_traits
            resample_info = get_resample_info(resample_dir)
            
            print(f"Resample directory: {resample_dir}")
            print(f"Total files: {resample_info['num_files']}")
            print(f"Total cycles: {resample_info['total_cycles']}")
            print(f"Progress log: {progress_log if progress_log else 'stdout only'}\n")
            
            # Initialize position-level result accumulators
            # position_name -> count of CAAS across all files
            position_counts = {}
            total_cycles = resample_info['total_cycles']
            
            # OPTIMIZATION: Filter to only test positions found in discovery
            positions_list = list(sliced_object.d)
            total_positions = len(positions_list)
            
            # Map position dict to (pos_dict, schemes_set) tuples
            positions_with_schemes = []
            
            if discovery_file:
                print(f"Filtering positions based on discovery results from: {discovery_file}")
                position_schemes = parse_discovery_positions(discovery_file, the_genename)
                
                if position_schemes:
                    for pos_dict in positions_list:
                        for species, aa_info in pos_dict.items():
                            pos_num = aa_info.split("@")[1]
                            if pos_num in position_schemes:
                                positions_with_schemes.append((pos_dict, position_schemes[pos_num]))
                                break
                    
                    speedup = total_positions / max(1, len(positions_with_schemes))
                    print(f"✓ Optimization: Testing only {len(positions_with_schemes)} CAAS positions (from {total_positions} total)")
                    print(f"✓ Speedup: {speedup:.1f}× fewer positions to test")
                    print(f"✓ Tests reduced: {total_positions * total_cycles:,} → {len(positions_with_schemes) * total_cycles:,}\n")
                else:
                    print("No CAAS positions found in discovery file. Processing all positions.\n")
                    positions_with_schemes = [(pos, None) for pos in positions_list]
            else:
                # No discovery file - test all positions with all schemes
                positions_with_schemes = [(pos, None) for pos in positions_list]
            
            # Process each file sequentially
            start_time = time.time()
            
            for file_idx, (file_path, file_config) in enumerate(simtrait_revive_from_dir(resample_dir), 1):
                file_start = time.time()
                
                log_progress(file_idx, resample_info['num_files'], start_time, 
                            log_file=progress_log, prefix=f"Processing file {os.path.basename(file_path)}")
                
                is_b0 = os.path.basename(file_path) == "resample_000.tab"
                
                # Process each position with its specific schemes
                for pos_dict, schemes in positions_with_schemes:
                    # Process position
                    processed_pos = process_position(pos_dict, multiconfig=file_config, species_in_alignment=sliced_object.species)
                    
                    # Run bootstrap with position-specific schemes
                    line_output = caasboot(
                        processed_pos,
                        genename=the_genename,
                        list_of_traits=file_config.alltraits,
                        maxgaps_fg=max_fg_gaps,
                        maxgaps_bg=max_bg_gaps,
                        maxgaps_all=max_overall_gaps,
                        maxmiss_fg=max_fg_miss,
                        maxmiss_bg=max_bg_miss,
                        maxmiss_all=max_overall_miss,
                        multiconfig=file_config,
                        miss_pair=miss_pair,
                        max_conserved=max_conserved,
                        admitted_patterns=the_admitted_patterns,
                        cycles=file_config.cycles,
                    caap_mode=caap_mode,
                    discovery_schemes=schemes,
                    debug_rejects=is_b0,
                    groups_out=groups_handle,
                    perm_discovery_out=perm_discovery_handle
                )
                    
                    # Accumulate counts for this position
                    # In CAAP mode, each position returns multiple lines (one per scheme)
                    if caap_mode:
                        # Multiple lines separated by newline
                        for line in line_output.split("\n"):
                            parts = line.split("\t")
                            position_name = parts[0]
                            scheme_name = parts[1]
                            count = int(parts[2])
                            
                            key = (position_name, scheme_name)
                            if key not in position_counts:
                                position_counts[key] = 0
                            position_counts[key] += count
                    else:
                        # Single line per position
                        parts = line_output.split("\t")
                        position_name = parts[0]
                        count = int(parts[1])
                        
                        if position_name not in position_counts:
                            position_counts[position_name] = 0
                        position_counts[position_name] += count
                
                file_elapsed = time.time() - file_start
                print(f"  → File completed in {format_time(file_elapsed)}\n")
            
            # Write final aggregated results
            print(f"\n{'='*80}")
            print(f"Writing aggregated results to {output_file}")
            print(f"{'='*80}\n")
            
            with open(output_file, "w") as ooout:
                if caap_mode:
                    # Sort by position then scheme
                    for key in sorted(position_counts.keys()):
                        position_name, scheme_name = key
                        count = position_counts[key]
                        empval = count / total_cycles
                        outline = "\t".join([position_name, scheme_name, str(count), str(total_cycles), str(empval)])
                        print(outline, file=ooout)
                else:
                    # Classical CAAS mode
                    for position_name in sorted(position_counts.keys()):
                        count = position_counts[position_name]
                        empval = count / total_cycles
                        outline = "\t".join([position_name, str(count), str(total_cycles), str(empval)])
                        print(outline, file=ooout)
            
            total_elapsed = time.time() - start_time
            print(f"✓ Bootstrap complete in {format_time(total_elapsed)}")
            print(f"✓ Results written to {output_file}\n")
        
        else:
            # Single file mode (backward compatibility)
            print("caastools found", resampled_traits.cycles, "resamplings")
            
            # OPTIMIZATION: Filter to only test positions found in discovery
            positions_list = list(sliced_object.d)
            total_positions = len(positions_list)
            
            # Map position dict to (pos_dict, schemes_set) tuples
            positions_with_schemes = []
            
            if discovery_file:
                print(f"\nFiltering positions based on discovery results from: {discovery_file}")
                position_schemes = parse_discovery_positions(discovery_file, the_genename)
                
                if position_schemes:
                    # Filter positions: only keep those that match discovery positions
                    for pos_dict in positions_list:
                        # Extract position number from the position dictionary
                        # Format is "AA@position_number" in values
                        for species, aa_info in pos_dict.items():
                            pos_num = aa_info.split("@")[1]
                            if pos_num in position_schemes:
                                positions_with_schemes.append((pos_dict, position_schemes[pos_num]))
                                break
                    
                    speedup = total_positions / max(1, len(positions_with_schemes))
                    print(f"✓ Optimization: Testing only {len(positions_with_schemes)} positions with discovered schemes (from {total_positions} total)")
                    print(f"✓ Speedup: {speedup:.1f}× fewer positions to test")
                    if caap_mode:
                        # Count how many scheme tests we're avoiding
                        total_scheme_tests = sum(len(schemes) if schemes and "CAAS" not in schemes else 6 for _, schemes in positions_with_schemes)
                        max_scheme_tests = len(positions_with_schemes) * 6  # 6 schemes max
                        scheme_speedup = max_scheme_tests / max(1, total_scheme_tests)
                        print(f"✓ Scheme optimization: Testing {total_scheme_tests} schemes (vs {max_scheme_tests} if testing all)")
                        print(f"✓ Scheme speedup: {scheme_speedup:.1f}× fewer scheme tests per position")
                    print(f"✓ Tests reduced: {total_positions * resampled_traits.cycles:,} → {len(positions_with_schemes) * resampled_traits.cycles:,}\n")
                else:
                    print("No positions found in discovery file. Processing all positions with all schemes.\n")
                    positions_with_schemes = [(pos, None) for pos in positions_list]
            else:
                # No discovery file - test all positions with all schemes
                positions_with_schemes = [(pos, None) for pos in positions_list]

            # Step 3 & 4: process positions with their specific schemes and run bootstrap
            output_lines = []
            for pos_dict, schemes in positions_with_schemes:
                # Process position
                processed_pos = process_position(pos_dict, multiconfig=resampled_traits, species_in_alignment=sliced_object.species)
                
                # Run bootstrap with position-specific schemes
                line_output = caasboot(
                    processed_pos,
                    list_of_traits=resampled_traits.alltraits,
                    genename=the_genename,
                    maxgaps_fg=max_fg_gaps,
                    maxgaps_bg=max_bg_gaps,
                    maxgaps_all=max_overall_gaps,
                    maxmiss_fg=max_fg_miss,
                    maxmiss_bg=max_bg_miss,
                    maxmiss_all=max_overall_miss,
                    multiconfig=resampled_traits,
                    miss_pair=miss_pair,
                    max_conserved=max_conserved,
                    admitted_patterns=the_admitted_patterns,
                    cycles=resampled_traits.cycles,
                caap_mode=caap_mode,
                discovery_schemes=schemes,
                debug_rejects=False,
                groups_out=groups_handle,
                perm_discovery_out=perm_discovery_handle
            )
                output_lines.append(line_output)

            ooout = open(output_file, "w")

            if caap_mode:
                # Each line may contain multiple scheme results separated by newlines
                for line_output in output_lines:
                    for line in line_output.split("\n"):
                        print(line, file=ooout)
            else:
                # Classical CAAS mode
                for line in output_lines:
                    print(line, file=ooout)
            
            ooout.close()
            
            print(f"Results written to {output_file}")
    finally:
        if groups_handle:
            groups_handle.close()
        if perm_discovery_handle:
            perm_discovery_handle.close()

# FUNCTION pval()
# Returns a dictionary with the pvalue

def pval(bootstrap_result):
    with open(bootstrap_result) as h:
        thelist = h.read().splitlines()
    
    d = {}

    for line in thelist:
        try:
            c = line.split("\t")
            d[c[0]] = c[2]
        except:
            pass
    
    return d
