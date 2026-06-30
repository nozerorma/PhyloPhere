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

Pair-aware / BLAS implementation: Miguel Ramon (miguel.ramon@upf.edu)

MODULE NAME: boot_vec.py
DESCRIPTION: Vectorized (Level-3 BLAS) reimplementation of the per-cycle CAAS/CAAP
             counting performed scalarly by modules.boot.caasboot().

             The convergence test at a fixed alignment position is a hidden GEMM:
             the amino acid carried by each species is fixed; only the foreground /
             background labeling changes across the B bootstrap cycles. Encoding
             species -> group as a one-hot matrix G_p (shared across all cycles),
             "how many fg/bg members fall in each group, per cycle" is exactly

                 C_fg = F  @ G_p      (B x species) @ (species x groups) -> (B x groups)
                 C_bg = Bg @ G_p

             and the CAAS condition + pattern classification + gap/missing filters
             are boolean / integer reductions over the group axis of C_fg / C_bg.
             Positions are stacked column-wise into a single G_concat so one matmul
             does the whole alignment.

             This mirrors, bit-for-bit, the decision taken by:
               - caas_id.iscaas()                 (classical CAAS / scheme US)
               - caap_id.check_caap_pattern()     (CAAP schemes GS1-GS4)
               - boot.filter_for_gaps / filter_for_missings
               - boot.caasboot()'s admitted_patterns substring membership

             Parity is proven by boot_vec_equivtest.py (run before trusting it).

DEPENDENCIES: numpy, modules.caap_id (group scheme dictionaries)
CALLED BY:    modules.boot (boot_on_single_alignment, hot counting path)

NOTE ON SCOPE: this path computes the empirical-p COUNTS only (the bootstrap's sole
product in production). The debug exports (--export_groups / --export_perm_discovery)
remain on the scalar path in boot.py; the vectorized path is bypassed when they are on.
'''

import numpy as np

from modules.caap_id import US, GS1, GS2, GS3, GS4


# Ambiguity codes are treated as gaps (no resolved amino acid), exactly as
# caas_id.process_position() does. "-" is the literal gap.
_AMBIGUOUS_AAS = frozenset({"X", "B", "Z", "J", "U"})
_GAP_SYMBOLS = frozenset({"-"}) | _AMBIGUOUS_AAS

# Scheme registry, ordered as caasboot emits its per-scheme output lines.
_SCHEME_ORDER = ["US", "GS1", "GS2", "GS3", "GS4"]
_SCHEME_MAP = {"US": US, "GS1": GS1, "GS2": GS2, "GS3": GS3, "GS4": GS4}


def _threshold(value):
    """Mirror the scalar 'NO' / int-string threshold convention.

    Returns None when no filter applies, otherwise the integer cap.
    """
    if value is None:
        return None
    if isinstance(value, str):
        if value == "NO":
            return None
        return int(value)
    return int(value)


def _admitted_pattern_flags(admitted_patterns):
    """Reproduce caasboot's substring membership test for pattern admission.

    In the scalar code admitted_patterns is the raw string (e.g. "1,2,3") and the
    test is `check.pattern in admitted_patterns`, i.e. single-character substring
    membership. We precompute, for each pattern label 1..4, whether it is admitted.
    Pattern "null" is never admitted.
    """
    if admitted_patterns is None:
        container = ["1", "2", "3"]
    else:
        container = admitted_patterns  # str or list; `in` works for both
    return {k: (str(k) in container) for k in (1, 2, 3, 4)}


class VectorizedBootstrap:
    """Per-resample-file (per multiconfig) labeling masks, reused across positions.

    The foreground / background membership of each labeling (bootstrap cycle) is
    independent of alignment position, so F and Bg are built once and reused for
    every position in the gene.
    """

    def __init__(self, cfg, species_in_alignment):
        # Species universe = every species referenced by any labeling, plus every
        # species present in the alignment. Alignment-only species (never sampled
        # into any labeling) get all-zero rows in F/Bg and are harmless.
        species = set(species_in_alignment)
        for t in cfg.alltraits:
            species.update(cfg.trait2fg.get(t, ()))
            species.update(cfg.trait2bg.get(t, ()))
        self.species = sorted(species)
        self.sp_index = {sp: i for i, sp in enumerate(self.species)}
        self.n_sp = len(self.species)

        # Deduplicate trait order exactly like simtrait_revive (preserves order).
        self.alltraits = list(dict.fromkeys(cfg.alltraits))
        self.B = len(self.alltraits)

        # Labeling membership masks (B x species), float32 so the matmuls hit BLAS.
        self.F = np.zeros((self.B, self.n_sp), dtype=np.float32)
        self.Bg = np.zeros((self.B, self.n_sp), dtype=np.float32)
        for b, t in enumerate(self.alltraits):
            for sp in cfg.trait2fg.get(t, ()):
                j = self.sp_index.get(sp)
                if j is not None:
                    self.F[b, j] = 1.0
            for sp in cfg.trait2bg.get(t, ()):
                j = self.sp_index.get(sp)
                if j is not None:
                    self.Bg[b, j] = 1.0

        # Species present in this gene's alignment (others are "missing").
        self.in_alignment = np.zeros(self.n_sp, dtype=bool)
        for sp in species_in_alignment:
            j = self.sp_index.get(sp)
            if j is not None:
                self.in_alignment[j] = True

    # -- per-position encoding ------------------------------------------------

    def _position_symbols(self, pos_dict):
        """Return (symbols, gapvec, missvec) over the species universe for a position.

        symbols[i] is the uppercased amino acid for species i if it is present in
        the alignment AND not a gap/ambiguity code, else None.
        gapvec[i]  = 1.0 if species i is present in alignment but gapped/ambiguous.
        missvec[i] = 1.0 if species i is absent from the alignment (missing).
        """
        symbols = [None] * self.n_sp
        gapvec = np.zeros(self.n_sp, dtype=np.float32)
        missvec = np.zeros(self.n_sp, dtype=np.float32)

        present = np.zeros(self.n_sp, dtype=bool)
        for sp, val in pos_dict.items():
            j = self.sp_index.get(sp)
            if j is None:
                continue
            present[j] = True
            aa = val.split("@")[0].upper()
            if aa in _GAP_SYMBOLS:
                gapvec[j] = 1.0
            else:
                symbols[j] = aa

        # Missing = referenced by labelings/universe but not present in alignment.
        missvec[~present] = 1.0
        return symbols, gapvec, missvec

    @staticmethod
    def _group_block(symbols, scheme_dict):
        """Build a one-hot (species x groups) block for one position under one scheme.

        scheme_dict is None for classical CAAS (scheme US handled as identity over
        the observed symbols, exactly reproducing iscaas which treats every
        non-gap symbol as its own group, including non-standard characters).

        For an explicit scheme dict, symbols absent from the dict contribute no
        column (mirroring check_caap_pattern's scheme_dict.get(aa) is None filter),
        WITHOUT being counted as gaps.

        Returns an (n_sp x n_groups) float32 array (n_groups may be 0).
        """
        n_sp = len(symbols)
        # Map each species' symbol to a group label.
        group_labels = [None] * n_sp
        for i, aa in enumerate(symbols):
            if aa is None:
                continue
            if scheme_dict is None:
                group_labels[i] = aa  # identity: each symbol is its own group
            else:
                g = scheme_dict.get(aa)
                if g is not None:
                    group_labels[i] = g

        # Stable ordering of group columns (counting is order-invariant; ordering
        # only needs to be self-consistent within the block).
        cols = sorted({g for g in group_labels if g is not None})
        if not cols:
            return np.zeros((n_sp, 0), dtype=np.float32)
        col_index = {g: k for k, g in enumerate(cols)}
        block = np.zeros((n_sp, len(cols)), dtype=np.float32)
        for i, g in enumerate(group_labels):
            if g is not None:
                block[i, col_index[g]] = 1.0
        return block

    # -- main counting --------------------------------------------------------

    def count(self, positions_with_schemes, genename,
              maxgaps_fg, maxgaps_bg, maxgaps_all,
              maxmiss_fg, maxmiss_bg, maxmiss_all,
              max_conserved, admitted_patterns, caap_mode,
              b_chunk=20000, collect_hits=False):
        """Count CAAS/CAAP hits per (position[, scheme]) across all B labelings.

        positions_with_schemes: list of (pos_dict, schemes_set_or_None), the same
            structure boot.boot_on_single_alignment builds. In caap_mode the
            schemes_set selects which schemes to test for that position (None or a
            set containing "CAAS" -> all five; otherwise only the named schemes).

        collect_hits: when True, ALSO return, per key, the list of labeling (trait)
            names that are CAAS hits — needed to materialize the perm_discovery rows
            (BOOTSTRAP_PERMS). Hits are gathered in ascending global labeling index.

        Returns:
            collect_hits False:
                caap_mode False -> {position_name: count}
                caap_mode True  -> {(position_name, scheme_name): count}
            collect_hits True:
                (counts_dict, hits_dict) with hits_dict[key] = [trait_name, ...]
        where position_name = f"{genename}@{position_number}".
        """
        adm = _admitted_pattern_flags(admitted_patterns)
        g_fg = _threshold(maxgaps_fg)
        g_bg = _threshold(maxgaps_bg)
        g_all = _threshold(maxgaps_all)
        m_fg = _threshold(maxmiss_fg)
        m_bg = _threshold(maxmiss_bg)
        m_all = _threshold(maxmiss_all)

        # ---- Build per-position vectors and the stacked group matrix (jobs) ----
        n_pos = len(positions_with_schemes)
        gapmat = np.zeros((self.n_sp, n_pos), dtype=np.float32)
        missmat = np.zeros((self.n_sp, n_pos), dtype=np.float32)
        position_names = []

        # A "job" is a (position, scheme) unit producing one output count.
        # job = (pos_idx, scheme_name_or_None, col_start, col_end)
        jobs = []
        blocks = []
        col_cursor = 0

        for pi, (pos_dict, schemes_set) in enumerate(positions_with_schemes):
            symbols, gapvec, missvec = self._position_symbols(pos_dict)
            gapmat[:, pi] = gapvec
            missmat[:, pi] = missvec

            # Position number (matches processed_position.position).
            pos_num = None
            for v in pos_dict.values():
                parts = v.split("@")
                if len(parts) > 1:
                    pos_num = parts[1]
                    break
            position_names.append(f"{genename}@{pos_num}")

            if caap_mode:
                scheme_names = self._schemes_to_test(schemes_set)
            else:
                scheme_names = [None]  # classical -> single identity job

            for sname in scheme_names:
                scheme_dict = None if sname is None else _SCHEME_MAP[sname]
                block = self._group_block(symbols, scheme_dict)
                ng = block.shape[1]
                blocks.append(block)
                jobs.append((pi, sname, col_cursor, col_cursor + ng))
                col_cursor += ng

        # Initialize result with zero counts so positions/schemes with no hits and
        # the "no valid traits" case both report 0 (matching caasboot).
        if caap_mode:
            results = {}
            for (pi, sname, _s, _e) in jobs:
                results[(position_names[pi], sname)] = 0
        else:
            results = {position_names[pi]: 0 for pi in range(n_pos)}

        hits = {key: [] for key in results} if collect_hits else None

        if self.B == 0 or col_cursor == 0:
            # No labelings, or no groups anywhere -> all counts stay 0.
            # Gaps could still reject; counts are 0 regardless.
            return (results, hits) if collect_hits else results

        G_concat = np.hstack(blocks) if blocks else np.zeros((self.n_sp, 0), np.float32)

        # ---- Stream over cycles (B) in chunks to bound memory ----
        for start in range(0, self.B, b_chunk):
            end = min(start + b_chunk, self.B)
            Fc = self.F[start:end]      # (bc x species)
            Bc = self.Bg[start:end]

            # The big GEMMs: per-cycle group counts for every position at once.
            C_fg_all = Fc @ G_concat    # (bc x total_groups)
            C_bg_all = Bc @ G_concat

            # Per-cycle gap / missing counts per position (also GEMMs).
            gaps_fg = Fc @ gapmat       # (bc x n_pos)
            gaps_bg = Bc @ gapmat
            miss_fg = Fc @ missmat
            miss_bg = Bc @ missmat

            for (pi, sname, cs, ce) in jobs:
                C_fg = C_fg_all[:, cs:ce]
                C_bg = C_bg_all[:, cs:ce]

                # --- gap / missing labeling validity (filter_for_gaps/missings) ---
                gfg = gaps_fg[:, pi]
                gbg = gaps_bg[:, pi]
                mfg = miss_fg[:, pi]
                mbg = miss_bg[:, pi]
                valid = np.ones(C_fg.shape[0], dtype=bool)
                if g_all is not None:
                    valid &= (gfg + gbg) <= g_all
                if g_fg is not None:
                    valid &= gfg <= g_fg
                if g_bg is not None:
                    valid &= gbg <= g_bg
                if m_all is not None:
                    valid &= (mfg + mbg) <= m_all
                if m_fg is not None:
                    valid &= mfg <= m_fg
                if m_bg is not None:
                    valid &= mbg <= m_bg

                if ce == cs:
                    # No groups at this position/scheme -> no CAAS possible.
                    continue

                present_fg = C_fg > 0
                present_bg = C_bg > 0
                nfg_unique = present_fg.sum(axis=1)
                nbg_unique = present_bg.sum(axis=1)

                shared = present_fg & present_bg
                shared_fg = (C_fg * shared).sum(axis=1)
                shared_bg = (C_bg * shared).sum(axis=1)
                overlap = np.minimum(shared_fg, shared_bg)
                non_fg = (C_fg * ~present_bg).sum(axis=1)
                non_bg = (C_bg * ~present_fg).sum(axis=1)

                caas = (overlap <= max_conserved) & ((non_fg >= 2) | (non_bg >= 2))

                # --- pattern classification (mirror iscaas / check_caap_pattern) ---
                is_null = (nfg_unique == 0) | (nbg_unique == 0)
                p1 = (nfg_unique == 1) & (nbg_unique == 1)
                p2 = (nfg_unique == 1) & (nbg_unique != 1)
                p3 = (nfg_unique != 1) & (nbg_unique == 1)
                p4 = (nfg_unique != 1) & (nbg_unique != 1)
                admitted = np.zeros(C_fg.shape[0], dtype=bool)
                if adm[1]:
                    admitted |= p1
                if adm[2]:
                    admitted |= p2
                if adm[3]:
                    admitted |= p3
                if adm[4]:
                    admitted |= p4
                admitted &= ~is_null

                hit = valid & caas & admitted
                count = int(hit.sum())

                key = (position_names[pi], sname) if caap_mode else position_names[pi]
                results[key] += count

                if collect_hits and count:
                    # Map local (within-chunk) hit rows to global labeling names.
                    local_idx = np.nonzero(hit)[0]
                    hits[key].extend(self.alltraits[start + int(i)] for i in local_idx)

        return (results, hits) if collect_hits else results

    @staticmethod
    def _schemes_to_test(schemes_set):
        """Mirror caasboot's discovery_schemes -> schemes_to_test resolution."""
        if not schemes_set:
            return list(_SCHEME_ORDER)
        if "CAAS" in schemes_set:
            return list(_SCHEME_ORDER)
        return [s for s in _SCHEME_ORDER if s in schemes_set]
