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

MODULE NAME:    runslice.py
DESCRIPTION:    The slicer function.
DEPENDENCIES:   TBD
'''
from modules.alimport import *

### Function runslice (collects the slicer inputs and runs it)
def runslice(options_object):

    # Inputs (transferring parsed options_object to variables)
    the_alignment = options_object.single_alignment
    alignment_format = options_object.ali_format
    
    # Only the global gap-fraction prefilter remains. The former diversity
    # threshold was unused; scanning all design files per gene served no purpose.
    return slice(the_alignment, alignment_format, 0, float(options_object.max_gaps_pos_string))
