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

def _new_multicfg():
    """An empty multicfg: the labelings container the kernel reads (see `labelings_to_multicfg`)."""

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

        
    return multicfg()


def simtrait_revive(traitfile):
    """Load resampled traits from a single file (backward compatibility)"""
    from modules.design_io import read_labelings
    return _fill(_new_multicfg(), read_labelings([traitfile]))


def labelings_to_multicfg(labelings):
    """The multicfg of in-memory labelings {tag: (foreground species, background species)}, built exactly as
    `simtrait_revive` builds it from a file (pair k is fg[k] with bg[k]); the tags are not checked against the
    cycle-tag grammar, so derived labelings (the orientations of the positional universe) can use their own."""
    return _fill(_new_multicfg(), labelings)


def _fill(z, labelings):
    z.cycles = len(labelings)
    for cycleid, (fg, bg) in labelings.items():
        for group, species in (("1", fg), ("0", bg)):
            for i, sp in enumerate(species, 1):
                z.update_dictionary(cycleid, sp, group, str(i))
    z.alltraits = list(labelings)
    return z
