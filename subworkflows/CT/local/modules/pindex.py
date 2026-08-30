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
update_dictionary()         a function to update a dictionary with
                            new information. No, there is no built-in method for
                            this.

load_cfg_dictionary()       Loads the multi cfg dictionary

'''

import glob

# FUNCTION update dictionary
# A function to update a dictionary with new information. No, there is no built-in method for this.

def update_dictionary(dictionary, key, value):
    try:
        dictionary[key].append(value)
    except:
        dictionary[key] = [value]

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

    if mode == "multi":
        cfg_files = [f for f in glob.glob(os.path.join(input_path, "*")) if f.endswith(".tab") or "traitfile" in os.path.basename(f)]
        for x in cfg_files:
            if not os.path.isfile(x):
                continue
            traitname = os.path.basename(x)
            z.alltraits.append(traitname)

            try:
                with open(x, encoding="utf-8", errors="ignore") as singlecfg_f:
                    singlecfg = singlecfg_f.read().splitlines()
            except Exception:
                continue
            
            for line in singlecfg:
                try:
                    c = line.split()
                    if len(c) < 3:
                        raise ValueError(
                            f"ERROR: Paired mode is mandatory but config file has only {len(c)} columns.\n"
                            f"This implementation of caastools is designed for paired contrasts.\n"
                            f"Trait config files must contain exactly: species, trait, pair.\n"
                            f"For sparse species analyses, please use the original caastools framework:\n"
                            f"Barteri et al. (https://academic.oup.com/bioinformatics/article/39/10/btad623/7319365)"
                        )
                    z.update_dictionary(traitname, c[0], c[1], pair=c[2])
                except ValueError as e:
                    raise e
                except:
                    pass
    
    elif mode == "mono":

        traitname = input_path.split("/")[-1]
        z.alltraits.append(traitname)

        with open(input_path) as singlecfg_f:
            singlecfg =  singlecfg_f.read().splitlines()
            
        for line in singlecfg:
            try:
                c = line.split()
                if len(c) < 3:
                    raise ValueError(
                        f"ERROR: Paired mode is mandatory but config file has only {len(c)} columns.\n"
                        f"This implementation of caastools is designed for paired contrasts.\n"
                        f"Trait config files must contain exactly: species, trait, pair.\n"
                        f"For sparse species analyses, please use the original caastools framework:\n"
                        f"Barteri et al. (https://academic.oup.com/bioinformatics/article/39/10/btad623/7319365)"
                    )
                z.update_dictionary(traitname, c[0], c[1], pair=c[2])
            except ValueError as e:
                raise e
            except:
                pass
    
    # Validate all species have pairs
    species_without_pairs = [s for s in z.s2t.keys() if s not in z.species2pair]
    if species_without_pairs:
        raise ValueError(
            f"ERROR: The following species lack pair assignments: {', '.join(species_without_pairs)}\n"
            f"This implementation of caastools is designed for paired contrasts.\n"
            f"For sparse species analyses, please use the original caastools framework:\n"
            f"Barteri et al. (https://academic.oup.com/bioinformatics/article/39/10/btad623/7319365)"
        )

    return z