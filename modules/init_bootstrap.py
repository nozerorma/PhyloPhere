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


MODULE NAME: INIT BOOTSTRAP
DESCRIPTION: Initialises the bootstrap through different strategy (random, phylogeny and permulations)
DEPENDENCIES: none
CALLED BY: ct
'''


import random
import os
import sys
import dendropy


# FUNCTION readtree(). Reads a tree and releases an object with some information (list of species, patristic distances matrix, bins)

def readtree(tree_file, tree_schema = "newick"):


    class topology():

        def __init__(self):
            self.treename = ""
            self.species_list = []
            self.pdd = {}
            self.farthest = ""
            self.maxd = 0

            self.dff = {}
            self.species_to_index = {}
            self.binlist = []
            self.index_dictionary = {}

    z = topology()

    # FUNCTION create_tag()
    # Two species into a single tag

    def create_tag(species1, species2):
        t = [species1, species2]
        t.sort()

        return "@".join(t)


    # Step 1: fetch the distances
    z.pdd = {}

    tree = dendropy.Tree.get(
        path=tree_file,
        schema=tree_schema)
    pdc = tree.phylogenetic_distance_matrix()

    for i, t1 in enumerate(tree.taxon_namespace[:-1]):
        for t2 in tree.taxon_namespace[i+1:]:
            thecouple = [t1.label.replace(" ", "_"), t2.label.replace(" ", "_")]
            z.species_list.append(t1.label.replace(" ", "_"))
            z.species_list.append(t2.label.replace(" ", "_"))
            thecouple.sort()
            key = "@".join(thecouple)
            z.pdd[key] = pdc(t1, t2)

    z.species_list = list(set(z.species_list))      # Eliminate duplicates

    # Step 2: fetch the farthest species

    z.pdd = dict(sorted(z.pdd.items(), key=lambda item: item[1]))
    z.farthest = list(z.pdd.items())[-1][0].split("@")[0]
    z.maxd = float(list(z.pdd.items())[-1][1])

    # Step 3: the distance from fartest (dff) dictionary

    z.dff[z.farthest] = 0

    for s in z.species_list:
        if s != z.farthest:
            t = create_tag(s, z.farthest)
            d = z.pdd[t]
            z.dff[s] = d

    # Sort the dictionary 
    z.dff = dict(sorted(z.dff.items(), key=lambda item: item[1]))

    return z


# FUNCTION simtrait() Resample trait function

def simtrait(fg_len, bg_len, template, tree_file, mode, groupfile, phenotype_values_file, cycles, simtraits_outfile, permulation_selection_strategy = "random", chunk_size = 500):
    
    # Class multicfg
    class multicfg():

        def __init__(self):
            self.s2t = {}
            self.alltraits = []
            self.trait2fg = {}
            self.trait2bg = {}
            self.paired_mode = False  # Bootstrap doesn't use paired mode
            
            # Pair-aware attributes (for compatibility with caas_id.py)
            self.species2pair = {}
            self.pair2fg_species = {}
            self.pair2bg_species = {}
            self.allpairs = []
            self._pair_cache = {}

        def get_pair(self, species):
            """Get the pair for a species (for compatibility with caas_id.py)"""
            return self.species2pair.get(species, None)

        def update_dictionary(self, traitname, species, group):
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
        
        def print_traits(self, outfile):

            o = open(outfile, "w")
            for x in self.trait2bg.keys():
                print("\t".join([   x,
                                    ",".join(self.trait2fg[x]),
                                    ",".join(self.trait2bg[x])
                                    ]), file = o)
            o.close()

    z = multicfg()

    # Step 1: import the species 

    t = readtree(tree_file)

    # WORKFLOW 1: bootstrap in random mode


    if mode == "random":

        
        # Species list.

        species = t.species_list
        allcycles = []

        # Run the bootstrap
        for i in range(1,cycles+1):
            cycle_tag = "b_" + str(i)

            z.alltraits.append(cycle_tag)

            resampled_fg = []
            resampled_bg = []

            while len(resampled_fg) < fg_len:
                k = random.choice(species)

                if k not in resampled_fg:
                    resampled_fg.append(k)
        
            while len(resampled_bg) < bg_len:
                k = random.choice(species)
                if k not in resampled_bg and k not in resampled_fg:
                    resampled_bg.append(k)
            
            allcycles.append([resampled_fg, resampled_bg, cycle_tag])
        
        for c in allcycles:

            z.alltraits.append(c[2])
            # Foreground
            for s in c[0]:
                z.update_dictionary(c[2], s, "1")

            # Background
            for s in c[1]:
                z.update_dictionary(c[2], s, "0")


        z.print_traits(simtraits_outfile)
        return z
    
    elif mode == "phylogeny-restricted-byfams":

        # Read from template

        try:

            with open(template) as tcf_handle:
                tcf = tcf_handle.read().splitlines()
        
            fg = []
            bg = []

            for l in tcf:
                try:
                    c = l.split("\t")
                    if c[1] == "1":
                        fg.append(c[0])
                    elif c[1] == "0":
                        bg.append(c[0])
                except:
                    pass
        except:
            print("ERROR: couldn't read template file. Input given:", template)
            exit()


        # Extract the groupfile

        with open(groupfile) as gf_handle:
            species_lines = gf_handle.read().splitlines()

        # Extract the species dictionaries

        s2g = {}
        g2s = {}

        for line in species_lines:
            try:
                c = line.split("\t")
                s2g[c[0]] = c[1]

                    # Update g2s
                try:
                    g2s[c[1]].append(c[0])
                except:
                    g2s[c[1]] = [c[0]]
            except:
                pass
        
        species = s2g.keys()

        allcycles = []

        # Run the bootstrap
        for i in range(1,cycles+1):
            cycle_tag = "b_" + str(i)

            z.alltraits.append(cycle_tag)

            resampled_fg = []
            resampled_bg = []

            for x in fg:
                flag = 0
                while flag != 1:
                    g = s2g[x]
                    x = random.choice(g2s[g])

                    if x not in resampled_fg:
                        resampled_fg.append(x)
                        flag = 1


            for x in bg:
                flag = 0
                while flag != 1:
                    g = s2g[x]
                    x = random.choice(g2s[g])

                    if x not in resampled_bg and x not in resampled_fg:
                        resampled_bg.append(x)
                        flag = 1


            resampled_bg.sort()
            resampled_fg.sort()
            
            allcycles.append([resampled_fg, resampled_bg, cycle_tag])

        for c in allcycles:

            z.alltraits.append(c[2])
            # Foreground
            for s in c[0]:
                z.update_dictionary(c[2], s, "1")

            # Background
            for s in c[1]:
                z.update_dictionary(c[2], s, "0")


        z.print_traits(simtraits_outfile)
        return z

    elif mode == "phylogeny-restricted-byfams":
        print("Sorry, this function is not implemented yet")
        exit()

    elif mode == "bm":

        if tree_file == "none":
            print("ERROR: for brownian motion simulation you must provide a rooted and fully dichotomic tree file in newick format.")
            print("See documentation.")
            exit()
        if phenotype_values_file == "none":
            print("ERROR: for brownian motion simulation you must provide a tsv file with phenotype values (1/0 for binary traits).")
            print("See documentation.")
            exit()
        
        # Create output directory if it doesn't exist
        if not os.path.exists(simtraits_outfile):
            os.makedirs(simtraits_outfile)
            print(f"Created output directory: {simtraits_outfile}")
        
        script_path = os.path.realpath(sys.argv[0]).split("/")
        script_path.remove(script_path[-1])
        rscript_path = "/".join(script_path) + "/permulations.r"

        os.system("ln -s " + rscript_path)
        
        r_line = " ".join([
            
            "./permulations.r",                     # Rscript
            tree_file,                              # args[1]
            template,                               # args[2]
            str(cycles),                            # args[3]
            permulation_selection_strategy,         # args[4]
            phenotype_values_file,                  # args[5]
            simtraits_outfile,                      # args[6] - now directory
            str(chunk_size)                         # args[7] - chunk size
        ])
        os.system(r_line)

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
            self.paired_mode = False  # Bootstrap doesn't use paired mode
            
            # Pair-aware attributes (for compatibility with caas_id.py)
            self.species2pair = {}
            self.pair2fg_species = {}
            self.pair2bg_species = {}
            self.allpairs = []
            self._pair_cache = {}

        def get_pair(self, species):
            """Get the pair for a species (for compatibility with caas_id.py)"""
            return self.species2pair.get(species, None)

        def update_dictionary(self, traitname, species, group):
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

            # Foreground update

            for s in fg:
                z.update_dictionary(cycleid, s, "1")


            # Background update

            for s in bg:
                z.update_dictionary(cycleid, s, "0")

        except:
            pass

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



'''
# Test
species_path = "_tests/sp2fam.210727.tab"
cfg_path = "_tests/AOB_config.txt"

st = simtrait(cfg_path, species_path, mode="phylogeny-guided")
st.print_traits("provaoutfile.tuf")

rst = simtrait_revive("provaoutfile.tuf")
rst.print_traits("provaoutfile.revived.tuf")
'''