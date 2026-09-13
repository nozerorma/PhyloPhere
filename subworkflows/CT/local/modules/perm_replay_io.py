#                      _              _     
#                     | |            | |    
#   ___ __ _  __ _ ___| |_ ___   ___ | |___ 
#  / __/ _` |/ _` / __| __/ _ \ / _ \| / __|
# | (_| (_| | (_| \__ \ || (_) | (_) | \__ \
#  \___\__,_|\__,_|___/\__\___/ \___/|_|___/


'''
A Convergent Amino Acid Substitution identification 
and analysis toolbox

Author:         Fabio Barteri (fabio.barteri@upf.edu)

Contributors:   Alejandro Valenzuela (alejandro.valenzuela@upf.edu)
                Xavier Farré (xfarrer@igtp.cat),
                David de Juan (david.juan@upf.edu).


MODULE NAME: PERM REPLAY IO
DESCRIPTION: Reloads already-generated resample/perm-replay labeling files
DEPENDENCIES: none
CALLED BY: ct
'''


import os


# FUNCTION simtrait_revive() revive resampled trait from

def simtrait_revive(traitfile):
    """Load resampled traits from a single file (backward compatibility)"""
    
    # Class multicfg

    class multicfg():

        def __init__(self):
            self.s2t = {}
            self.alltraits = []
            self.trait2fg = {}
            self.trait2bg = {}
            self.cycles = 0
            self.paired_mode = True

            # Pair-aware attributes (for compatibility with caas_id.py)
            # Unlike the single-config case in pindex.py, ONE multicfg here holds
            # every cycle in a resample file, and a species sits in a different
            # pair from one cycle to the next — so pair lookup is keyed by
            # (cycle, species), not by species alone.
            self.trait2species2pair = {}
            self.species2pair = {}       # kept for interface compatibility; unused
            self.pair2fg_species = {}
            self.pair2bg_species = {}
            self.allpairs = []
            self._pair_cache = {}

        def get_pair(self, species, trait=None):
            """Pair id for a species within a given cycle.

            `trait` is the cycle id (b_1, b_2, ...). Without it there is no
            unambiguous answer, since the same species belongs to different pairs
            in different cycles — callers inside a per-cycle loop must pass it.
            """
            if trait is not None:
                return self.trait2species2pair.get(trait, {}).get(species)
            return self.species2pair.get(species, None)

        def update_dictionary(self, traitname, species, group, pair=None):
            try:
                self.s2t[species].append(traitname + "_" + group)
            except:
                self.s2t[species] = [traitname + "_" + group]

            if group == "1":
                try:
                    self.trait2fg[traitname].append(species)
                except:
                    self.trait2fg[traitname] = [species]

            if group == "0":
                try:
                    self.trait2bg[traitname].append(species)
                except:
                    self.trait2bg[traitname] = [species]

            if pair is not None:
                self.trait2species2pair.setdefault(traitname, {})[species] = pair
                if pair not in self.allpairs:
                    self.allpairs.append(pair)
                if group == "1":
                    self.pair2fg_species.setdefault((traitname, pair), []).append(species)
                if group == "0":
                    self.pair2bg_species.setdefault((traitname, pair), []).append(species)

            self.alltraits.append(traitname)
        
        
        def print_traits(self, outfile):
            o = open(outfile, "w")
            for x in self.trait2bg.keys():
                print("\t".join([   x,
                                    ",".join(self.trait2fg[x]),
                                    ",".join(self.trait2bg[x])
                                    ]), file = o)
            o.close()

        
    # Declare multicfg instance

    z = multicfg()

    # Open the traitfile

    with open(traitfile) as tf_handle:
        tf = tf_handle.read().splitlines()
        z.cycles = len(tf)
    
    for line in tf:
        try:
            c = line.split("\t")
            cycleid = c[0]
            fg = c[1].split(",")
            bg = c[2].split(",")

            # Pair identity is POSITIONAL: the i-th foreground species and the
            # i-th background species are the two members of pair i.
            # permulations.R writes both lists in matched pair order (see the
            # sel_fg/sel_bg accumulation in lean_contrast_selector.R), so index
            # correspondence is the contract between the two files. Recovering it
            # here is what lets conserved_pair report WHICH pairs are conserved
            # rather than only how many — the identity permulation disambiguation
            # needs. If the lists differ in length (a Tier-2 style shortfall, or a
            # hand-written resample), only the common prefix gets a pair id and
            # the rest degrade to the previous "count only" behaviour.
            paired_n = min(len(fg), len(bg))

            # Foreground update

            for i, s in enumerate(fg):
                z.update_dictionary(cycleid, s, "1", str(i + 1) if i < paired_n else None)


            # Background update

            for i, s in enumerate(bg):
                z.update_dictionary(cycleid, s, "0", str(i + 1) if i < paired_n else None)

        except:
            pass
    
    # Deduplicate alltraits list (preserves order)
    z.alltraits = list(dict.fromkeys(z.alltraits))

    return z


# FUNCTION simtrait_revive_from_dir() Load resampled traits from directory of chunked files

def simtrait_revive_from_dir(resample_dir):
    """Load resampled traits from a directory containing chunked resample files.
    
    Returns a generator that yields (file_path, multicfg_object) tuples for sequential processing.
    This allows memory-efficient processing of large resample sets.
    
    Args:
        resample_dir: Path to directory containing resample_*.tab files
        
    Yields:
        tuple: (file_path, multicfg_object) for each resample file
    """
    import glob
    import re
    
    # Find all resample files
    pattern = os.path.join(resample_dir, "resample_*.tab")
    resample_files = glob.glob(pattern)
    
    if len(resample_files) == 0:
        print(f"ERROR: No resample files found in {resample_dir}")
        print(f"Expected files matching pattern: resample_*.tab")
        exit()
    
    # Sort files by number (resample_001.tab, resample_002.tab, ...)
    def extract_number(filepath):
        match = re.search(r'resample_(\d+)\.tab$', filepath)
        return int(match.group(1)) if match else 0
    
    resample_files.sort(key=extract_number)
    
    print(f"Found {len(resample_files)} resample files in {resample_dir}")
    
    # Yield each file with its loaded multicfg object
    for file_path in resample_files:
        yield file_path, simtrait_revive(file_path)


# FUNCTION get_resample_info() Get information about resampled traits without loading all data

def get_resample_info(resample_path):
    """Get information about resampled traits (total cycles, files) without loading all data.
    
    Args:
        resample_path: Path to resample file or directory
        
    Returns:
        dict: {'total_cycles': int, 'num_files': int, 'is_directory': bool, 'files': list}
    """
    import glob
    import re
    
    if os.path.isdir(resample_path):
        # Directory mode: count files and cycles
        pattern = os.path.join(resample_path, "resample_*.tab")
        resample_files = glob.glob(pattern)
        
        def extract_number(filepath):
            match = re.search(r'resample_(\d+)\.tab$', filepath)
            return int(match.group(1)) if match else 0
        
        resample_files.sort(key=extract_number)
        
        total_cycles = 0
        for file_path in resample_files:
            with open(file_path) as f:
                total_cycles += len(f.read().splitlines())
        
        return {
            'total_cycles': total_cycles,
            'num_files': len(resample_files),
            'is_directory': True,
            'files': resample_files
        }
    else:
        # Single file mode
        with open(resample_path) as f:
            total_cycles = len(f.read().splitlines())
        
        return {
            'total_cycles': total_cycles,
            'num_files': 1,
            'is_directory': False,
            'files': [resample_path]
        }