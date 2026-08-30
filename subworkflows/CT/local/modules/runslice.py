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
    
    # Alignment slice: 1- Calculate column treshold

    import os
    if os.path.isdir(options_object.config_file):
        import glob
        cfg_files = [f for f in glob.glob(os.path.join(options_object.config_file, "*")) if f.endswith(".tab") or "traitfile" in os.path.basename(f)]
        fg_counts = []
        bg_counts = []
        for f in cfg_files:
            if not os.path.isfile(f):
                continue
            try:
                with open(f, encoding="utf-8", errors="ignore") as cfg_handle:
                    lines = cfg_handle.read().splitlines()
                vals = [l.split("\t")[1] for l in lines if len(l.split("\t")) >= 3]
                fg_counts.append(vals.count("1"))
                bg_counts.append(vals.count("0"))
            except Exception:
                continue
        fg_species = min(fg_counts) if fg_counts else 0
        bg_species = min(bg_counts) if bg_counts else 0
    else:
        with open(options_object.config_file, encoding="utf-8", errors="ignore") as cfg_handle:
            lines = cfg_handle.read().splitlines()
        vals = [l.split("\t")[1] for l in lines if len(l.split("\t")) >= 3]
        fg_species = vals.count("1")
        bg_species = vals.count("0")



    # Alignment slicing: sum the null values (allowed_gaps + allowed_missing_species)

    sum_nulls_fg = 0

    for x in (options_object.max_fg_gaps_string,options_object.max_fg_miss_string):
        try:
            sum_nulls_fg = sum_nulls_fg + int(x)
        except:
            pass

    sum_nulls_bg = 0

    for x in (options_object.max_bg_gaps_string,options_object.max_bg_miss_string):
        try:
            sum_nulls_bg = sum_nulls_bg + int(x)
        except:
            pass



    fg_threshold = fg_species - sum_nulls_fg
    bg_threshold = bg_species - sum_nulls_bg

    c_threshold = min(fg_threshold, bg_threshold)


    # Alignment slice: 2- Filter positions (slice alignment)

    out = slice(the_alignment, alignment_format, c_threshold, float(options_object.max_gaps_pos_string))

    return out
