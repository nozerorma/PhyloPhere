# vulture_whitelist.py — Names that vulture should not flag as dead code.
# PhyloPhere | style/

# These names appear unused to static analysis but are referenced dynamically:
# argparse stores dest values as attributes on Namespace (never directly assigned),
# and several pipeline modules use getattr() or **vars(args) dispatch patterns.

# ── argparse Namespace attributes ─────────────────────────────────────────────
# Common argument destinations that are accessed via args.<name> in callers
# but not always visible in the same file (e.g., passed to sub-functions).

class Namespace:  # noqa: vulture whitelist stub
    alignment_dir = None
    genomic_info = None
    species_list = None
    metadata_caas = None
    output_prefix = None
    log_level = None
    n_permutations = None
    change_side = None
    tool = None
    traitfile = None
    tree = None
    discovery = None
    lengths = None
    output = None
    percentile = None
    gene = None
    batch_size = None
    min_sp = None
    min_pos = None
    winsorize = None
    clade = None
    perm_batches = None
    perms_per_batch = None
    input_dir = None
    output_dir = None
    trait = None
    trait_col = None
    label = None
    fg = None
    bg = None
    fasta = None
    newick = None
    phylip = None
    prune = None
    prune_data = None
    group = None
    groups = None
    caap = None
    scheme = None
    schemes = None
    nboot = None
    alpha = None
    target = None
    ref = None
    coords = None
    protein = None
    vcf = None
    annotate = None

# ── Module-level singletons accessed via import ───────────────────────────────

def logger(): pass          # module-level logger, accessed as `logger.info(...)` etc.
def get_logger(): pass      # factory used in package __init__ files

# ── Entry points called by Nextflow via subprocess ────────────────────────────

def main(): pass            # all main() functions are called by __name__ == "__main__"

# ── Conserved-pair plumbing kept latent by scoring_v2 T1 (roadmap decision G) ──
# T1 removed the conservation_gate multiplier from the ASR path score, but the
# conserved-pair machinery below stays wired so a later tier (T3) can fold it
# into the pairwise core. Guarded by test_conserved_pair_columns_survive.

class _ConservedPairLatent:  # noqa: vulture whitelist stub
    conserved_pair_scores = None
    conserved_pair_nodes = None
    conserved_pair_path_scores = None
    conserved_pair_path_nodes = None
    pair_ancestral = None
    pair_derived_top = None
    pair_derived_bot = None

def parse_conserved_ids(): pass      # path_scores.py — conserved id string parser
def _conserved_side_score(): pass    # path_scores.py — conservation-to-root walk

# ── Pytest / unittest hooks (if tests are added) ──────────────────────────────

def setUp(): pass
def tearDown(): pass
def setUpClass(): pass
def tearDownClass(): pass
