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

MODULE NAME: pindex.py
DESCRIPTION: phenotype indexing of a single trait or multiple traits
INPUTS:      single or multiple trait binary file

TABLE OF CONTENTS
------------------------------------------
load_cfg_dictionary()       Loads the multi cfg dictionary

'''

import glob

# FUNCTION update dictionary
# A function to update a dictionary with new information. No, there is no built-in method for this.


# FUNCTION load multi cfg dictionary
# Loads the multi cfg dictionary

def load_cfg(input_path, mode = "mono"):
    import os
    if os.path.isdir(input_path):
        mode = "multi"


    class multicfg():

        def __init__(self):
            self.s2t = {}
            self.alltraits = []
            self.trait2fg = {}
            self.trait2bg = {}
            self.paired_mode = True
            
            # Pair-aware attributes
            self.species2pair = {}
            self.trait_species2pair = {}
            self.pair2fg_species = {}
            self.pair2bg_species = {}
            self.allpairs = []
            
            # Pair lookup cache for performance
            self._pair_cache = {}

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
            
            # Handle pair information (mandatory paired mode)
            if pair is not None:
                self.species2pair[species] = pair
                self.trait_species2pair[(traitname, species)] = pair
                
                if pair not in self.allpairs:
                    self.allpairs.append(pair)
                
                if group == "1":
                    try:
                        self.pair2fg_species[(traitname, pair)].append(species)
                    except:
                        self.pair2fg_species[(traitname, pair)] = [species]
                    try:
                        self.pair2fg_species[pair].append(species)
                    except:
                        self.pair2fg_species[pair] = [species]
                
                if group == "0":
                    try:
                        self.pair2bg_species[(traitname, pair)].append(species)
                    except:
                        self.pair2bg_species[(traitname, pair)] = [species]
                    try:
                        self.pair2bg_species[pair].append(species)
                    except:
                        self.pair2bg_species[pair] = [species]
        
        def get_pair(self, species, trait=None):
            """Get pair for species with per-trait caching."""
            if trait is not None:
                key = (trait, species)
                if key not in self._pair_cache:
                    p = self.trait_species2pair.get(key)
                    if p is None:
                        p = self.species2pair.get(species)
                    self._pair_cache[key] = p
                return self._pair_cache[key]
            else:
                if species not in self._pair_cache:
                    self._pair_cache[species] = self.species2pair.get(species)
                return self._pair_cache[species]


    z = multicfg()

    from modules.design_io import design_files, read_design
    files = design_files(input_path)
    for hyp, pairs in read_design(input_path).items():
        traitname = files[hyp].name
        z.alltraits.append(traitname)
        for pair, sides in pairs.items():
            for group in ('1', '0'):
                z.update_dictionary(traitname, sides[group], group, pair=str(pair))

    return z