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
    Parse discovery output file and extract positions that have CAAS.
    
    Args:
        discovery_file: Path to discovery output file
        genename: Gene name to match (e.g., "BRCA1")
    
    Returns:
        set of position numbers (as strings) that had CAAS in discovery
    """
    if not exists(discovery_file):
        print(f"Warning: Discovery file {discovery_file} not found. Processing all positions.")
        return None
    
    caas_positions = set()
    
    try:
        with open(discovery_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("Gene"):  # Skip header
                    continue
                
                try:
                    # Discovery output format: Gene\tTrait\tPosition\t...
                    fields = line.split('\t')
                    if len(fields) >= 3:
                        gene = fields[0]
                        position = fields[2]
                        
                        # Match gene name and extract position number
                        if gene == genename:
                            caas_positions.add(position)
                except:
                    continue
        
        return caas_positions if len(caas_positions) > 0 else None
    
    except Exception as e:
        print(f"Warning: Error parsing discovery file: {e}. Processing all positions.")
        return None


# FUNCTION fetch_caas():
# fetches caas per each thing

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
# filters a trait for its gaps

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



def caasboot(processed_position, genename, list_of_traits, maxgaps_fg, maxgaps_bg, maxgaps_all, maxmiss_fg, maxmiss_bg, maxmiss_all, cycles, multiconfig, miss_pair=False, max_conserved=0, admitted_patterns=["1","2","3"], chunk_size=1000):
    """Chunked bootstrap - processes traits in batches to handle large resample files"""
    
    a = set(list_of_traits)
    b = set(processed_position.trait2aas_fg.keys())
    c = set(processed_position.trait2aas_bg.keys())
    valid_traits = list(a.intersection(b).intersection(c))
    
    if len(valid_traits) == 0:
        position_name = genename + "@" + str(processed_position.position)
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
                continue
            
            # Filter for missings
            mfg = processed_position.trait2miss_fg[trait]
            mbg = processed_position.trait2miss_bg[trait]
            
            if not filter_for_missings(maxmiss_bg, maxmiss_fg, maxmiss_all, mfg, mbg):
                continue
            
            filtered_traits.append(trait)
        
        # Pair-aware filtering (if needed)
        if multiconfig.paired_mode and miss_pair and len(filtered_traits) > 0:
            pair_filtered = []
            for trait in filtered_traits:
                # Check missing pairs
                miss_thresholds_equal = False
                if maxmiss_fg != "NO" and maxmiss_bg != "NO" and maxmiss_fg == maxmiss_bg:
                    miss_thresholds_equal = True
                elif maxmiss_fg == "NO" and maxmiss_bg == "NO" and maxmiss_all != "NO":
                    miss_thresholds_equal = True
                
                if miss_thresholds_equal:
                    miss_pairs_fg = set(processed_position.trait2miss_pairs_fg.get(trait, []))
                    miss_pairs_bg = set(processed_position.trait2miss_pairs_bg.get(trait, []))
                    if miss_pairs_fg != miss_pairs_bg:
                        continue
                
                # Check gap pairs
                gap_thresholds_equal = False
                if maxgaps_fg != "NO" and maxgaps_bg != "NO" and maxgaps_fg == maxgaps_bg:
                    gap_thresholds_equal = True
                elif maxgaps_fg == "NO" and maxgaps_bg == "NO" and maxgaps_all != "NO":
                    gap_thresholds_equal = True
                
                if gap_thresholds_equal:
                    gap_pairs_fg = set(processed_position.trait2gap_pairs_fg.get(trait, []))
                    gap_pairs_bg = set(processed_position.trait2gap_pairs_bg.get(trait, []))
                    if gap_pairs_fg != gap_pairs_bg:
                        continue
                
                pair_filtered.append(trait)
            filtered_traits = pair_filtered
        
        # Pattern check
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
    
    # Return aggregated result
    position_name = genename + "@" + str(processed_position.position)
    count = str(len(total_output_traits))
    empval = str(int(count)/cycles)
    outline = "\t".join([position_name, count, str(cycles), empval])
    
    return outline

# FUNCTION boot_on_single_alignment()
# Launches the bootstrap in several lines. Returns a dictionary gene@position --> pvalue

def boot_on_single_alignment(trait_config_file, resampled_traits, sliced_object, max_fg_gaps, max_bg_gaps, max_overall_gaps, max_fg_miss, max_bg_miss, max_overall_miss, the_admitted_patterns, output_file, miss_pair=False, max_conserved=0, discovery_file=None, progress_log=None):
    """
    Run bootstrap analysis on a single alignment.
    
    Supports both single-file and directory-based resampled traits:
    - Single file: resampled_traits is a multicfg object loaded from one file
    - Directory: resampled_traits is the directory path (string), files loaded sequentially
    
    Args:
        resampled_traits: multicfg object OR directory path (str) containing resample_*.tab files
        progress_log: Optional file path for logging progress
        ... (other parameters as before)
    """
    the_genename = sliced_object.genename
    
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
        
        if discovery_file:
            print(f"Filtering positions based on discovery results from: {discovery_file}")
            caas_positions = parse_discovery_positions(discovery_file, the_genename)
            
            if caas_positions:
                filtered_positions = []
                for pos_dict in positions_list:
                    for species, aa_info in pos_dict.items():
                        pos_num = aa_info.split("@")[1]
                        if pos_num in caas_positions:
                            filtered_positions.append(pos_dict)
                            break
                
                positions_list = filtered_positions
                speedup = total_positions / max(1, len(positions_list))
                print(f"✓ Optimization: Testing only {len(positions_list)} CAAS positions (from {total_positions} total)")
                print(f"✓ Speedup: {speedup:.1f}× fewer positions to test")
                print(f"✓ Tests reduced: {total_positions * total_cycles:,} → {len(positions_list) * total_cycles:,}\n")
            else:
                print("No CAAS positions found in discovery file. Processing all positions.\n")
        
        # Process each file sequentially
        start_time = time.time()
        
        for file_idx, (file_path, file_config) in enumerate(simtrait_revive_from_dir(resample_dir), 1):
            file_start = time.time()
            
            log_progress(file_idx, resample_info['num_files'], start_time, 
                        log_file=progress_log, prefix=f"Processing file {os.path.basename(file_path)}")
            
            # Process positions with this file's resampled traits
            processed_positions = map(
                functools.partial(process_position, 
                                multiconfig=file_config, 
                                species_in_alignment=sliced_object.species),
                positions_list
            )
            
            # Run bootstrap for each position with this file's traits
            output_lines = map(
                functools.partial(
                    caasboot,
                    list_of_traits=file_config.alltraits,
                    genename=the_genename,
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
                    cycles=file_config.cycles
                ),
                processed_positions
            )
            
            # Accumulate counts for each position
            for line in output_lines:
                parts = line.split("\t")
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
        
        if discovery_file:
            print(f"\nFiltering positions based on discovery results from: {discovery_file}")
            caas_positions = parse_discovery_positions(discovery_file, the_genename)
            
            if caas_positions:
                # Filter positions: only keep those that match discovery CAAS positions
                filtered_positions = []
                for pos_dict in positions_list:
                    # Extract position number from the position dictionary
                    # Format is "AA@position_number" in values
                    for species, aa_info in pos_dict.items():
                        pos_num = aa_info.split("@")[1]
                        if pos_num in caas_positions:
                            filtered_positions.append(pos_dict)
                            break
                
                positions_list = filtered_positions
                speedup = total_positions / max(1, len(positions_list))
                print(f"✓ Optimization: Testing only {len(positions_list)} CAAS positions (from {total_positions} total)")
                print(f"✓ Speedup: {speedup:.1f}× fewer positions to test")
                print(f"✓ Tests reduced: {total_positions * resampled_traits.cycles:,} → {len(positions_list) * resampled_traits.cycles:,}\n")
            else:
                print("No CAAS positions found in discovery file. Processing all positions.\n")

        # Step 3: processes the positions from imported alignment (process_position() from caas_id.py)
        processed_positions = map(functools.partial(process_position, multiconfig = resampled_traits, species_in_alignment = sliced_object.species), positions_list)

        # Step 4: extract the raw caas

        output_lines = map(
            functools.partial(
                caasboot,
                list_of_traits = resampled_traits.alltraits,
                genename = the_genename,
                maxgaps_fg = max_fg_gaps,
                maxgaps_bg = max_bg_gaps,
                maxgaps_all = max_overall_gaps,

                maxmiss_fg = max_fg_miss,
                maxmiss_bg = max_bg_miss,
                maxmiss_all = max_overall_miss,

                multiconfig = resampled_traits,
                miss_pair = miss_pair,
                max_conserved = max_conserved,
                admitted_patterns = the_admitted_patterns,
                cycles = resampled_traits.cycles) ,processed_positions
        )

        output_lines = list(output_lines)

        ooout = open(output_file, "w")

        for line in output_lines:
            print(line, file=ooout)
        
        ooout.close()
        
        print(f"Results written to {output_file}")

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
