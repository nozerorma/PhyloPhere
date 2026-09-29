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

MODULE NAME: perm_replay.py
DESCRIPTION: reruns and rescores resampled labelings for CAAS pattern matches
DEPENDENCIES: alimport.py, caas_id.py, pindex.py
CALLED BY: ct

'''


from modules.perm_replay_io import *
from modules.disco import process_position
from modules.caas_id import (
    US, GS1, GS2, GS3, GS4, SCHEMES,
    check_pattern, check_caap_pattern, iscaas,
    encode_to_groups, _pair_sort_key
)
from modules.alimport import *

import os
from os.path import exists
import functools
import re
import time
from datetime import datetime
import numpy as np

# ---------------------------------------------------------------------------
# FOP multi-hypothesis (Gap A) — base-cycle collapse of the fanned observed
# perm-replay. Under params.caas_perms_fop the observed perm-replay resamples over
# fop_labelings.tab, whose cycle tags are "<base>~H<m>" (one row per null cycle
# and fanned Dunn-independent alternative hypothesis). recovery_boot must be
# reported in BASE-CYCLE units: a base cycle HITS an observed (Gene@Position,
# scheme) iff ANY of its ~H<m> labelings calls a CAAS there (a cheap
# discovery-level OR — no ASR, no domain pooling; that is Gap B's job on a
# different set of cycles).
# ---------------------------------------------------------------------------
_FOP_H_SUFFIX = re.compile(r"~H\d+$")


def _fop_base_cycle(labeling_tag):
    """'b_12~H3' -> 'b_12'. A tag with no ~H<m> suffix is its own base cycle."""
    return _FOP_H_SUFFIX.sub("", labeling_tag)


def collapse_fop_hits_by_base(per_key_hits, all_labelings):
    """Collapse per-labeling CAAS hits to base-cycle units (Gap A).

    Args:
        per_key_hits: {key: [labeling_tag, ...]} — the labeling-unit hits the
            kernel produces. key is a position_name (classical)
            or (position_name, scheme_name) (caap_mode); tag is "<base>~H<m>".
        all_labelings: iterable of every labeling tag present in fop_labelings.tab
            (used for the denominator = number of distinct base cycles).

    Returns:
        (collapsed_counts, base_total) where collapsed_counts[key] is the number
        of DISTINCT base cycles with >=1 hit (the recovery_boot numerator) and
        base_total is the number of distinct base cycles (the denominator).
    """
    base_total = len({_fop_base_cycle(t) for t in all_labelings})
    collapsed = {k: len({_fop_base_cycle(t) for t in v}) for k, v in per_key_hits.items()}
    return collapsed, base_total

# ---------------------------------------------------------------------------
# Vectorized (Level-3 BLAS) counting kernel.
# ---------------------------------------------------------------------------

_CHUNK_MEM_BUDGET_MB = float(os.environ.get("CT_PERM_REPLAY_CHUNK_MEM_MB", "512"))

# Ambiguity codes count as gaps (no resolved amino acid), exactly as
# caas_id.process_position() does; "-" is the literal gap.
_AMBIGUOUS_AAS = frozenset({"X", "B", "Z", "J", "U"})
_GAP_SYMBOLS = frozenset({"-"}) | _AMBIGUOUS_AAS
_VEC_GAP_SYMBOLS = _GAP_SYMBOLS

_SCHEME_ORDER = ["US", "GS1", "GS2", "GS3", "GS4"]
_SCHEME_MAP = {"US": US, "GS1": GS1, "GS2": GS2, "GS3": GS3, "GS4": GS4}
_VEC_SCHEME_MAP = _SCHEME_MAP


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
    """Reproduce substring membership test for pattern admission.

    Precomputes, for each pattern label 1..4, whether it is admitted.
    Pattern 'null' is never admitted.
    """
    if admitted_patterns is None:
        container = ["1", "2", "3"]
    else:
        container = admitted_patterns  # str or list; `in` works for both
    return {k: (str(k) in container) for k in (1, 2, 3, 4)}


class VectorizedPermReplay:
    """Per-resample-file (per multiconfig) labeling masks, reused across positions.

    The foreground / background membership of each labeling (perm-replay cycle) is
    independent of alignment position, so F and Bg are built once and reused for
    every position in the gene.
    """

    def __init__(self, cfg, species_in_alignment):
        self._cfg = cfg
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

        # Species present in this gene's alignment (others are 'missing').
        self.in_alignment = np.zeros(self.n_sp, dtype=bool)
        for sp in species_in_alignment:
            j = self.sp_index.get(sp)
            if j is not None:
                self.in_alignment[j] = True

    # -- per-position encoding ------------------------------------------------

    def _position_symbols(self, pos_dict):
        """Return (symbols, gapvec, missvec) over the species universe for a position."""
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

        missvec[~present] = 1.0
        return symbols, gapvec, missvec

    @staticmethod
    def _group_block(symbols, scheme_dict):
        """Build a one-hot (species x groups) block for one position under one scheme."""
        n_sp = len(symbols)
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

        cols = sorted({g for g in group_labels if g is not None})
        if not cols:
            return np.zeros((n_sp, 0), dtype=np.float32)
        col_index = {g: k for k, g in enumerate(cols)}
        block = np.zeros((n_sp, len(cols)), dtype=np.float32)
        for i, g in enumerate(group_labels):
            if g is not None:
                block[i, col_index[g]] = 1.0
        return block

    def _ensure_pair_index(self):
        """Pair id of every species in every labeling (0 = not in that side), built on demand."""
        if getattr(self, "fg_pair", None) is not None:
            return
        self.fg_pair = np.zeros((self.B, self.n_sp), dtype=np.int16)
        self.bg_pair = np.zeros((self.B, self.n_sp), dtype=np.int16)
        ids = {}
        for b, t in enumerate(self.alltraits):
            for side, arr in ((self._cfg.trait2fg, self.fg_pair), (self._cfg.trait2bg, self.bg_pair)):
                for sp in side.get(t, ()):
                    j = self.sp_index.get(sp)
                    pair = self._cfg.get_pair(sp, t)
                    if j is not None and pair:  # like the scalar path, a falsy pair is ignored
                        arr[b, j] = ids.setdefault(pair, len(ids) + 1)

    def _pair_sets_differ(self, b, flag):
        """True when the fg and bg species flagged in `flag` (gap or missing) sit in different
        pair sets, both non-empty (the scalar path's `set_fg and set_bg and set_fg != set_bg`)."""
        on = flag > 0
        sfg = set(self.fg_pair[b][on & (self.fg_pair[b] > 0)].tolist())
        sbg = set(self.bg_pair[b][on & (self.bg_pair[b] > 0)].tolist())
        return bool(sfg and sbg and sfg != sbg)

    def _default_b_chunk(self, total_groups):
        bytes_per_row = max(1, 2 * (self.n_sp + total_groups) * 4)
        budget_bytes = _CHUNK_MEM_BUDGET_MB * 1024 * 1024
        return max(1, min(self.B, int(budget_bytes // bytes_per_row)))

    # -- main counting --------------------------------------------------------

    def count(self, positions_with_schemes, genename,
              maxgaps_fg, maxgaps_bg, maxgaps_all,
              maxmiss_fg, maxmiss_bg, maxmiss_all,
              max_conserved, admitted_patterns, caap_mode,
              b_chunk=None, collect_hits=False, miss_pair=False):
        """Count CAAS/CAAP hits per (position[, scheme]) across all B labelings.

        miss_pair mirrors the scalar path (caas_id.fetch_caas / disco.py): when the fg and bg
        thresholds are equal, a labeling whose fg and bg sides both have gapped (or missing)
        species, but in different pairs, is discarded.
        """
        adm = _admitted_pattern_flags(admitted_patterns)
        g_fg = _threshold(maxgaps_fg)
        g_bg = _threshold(maxgaps_bg)
        g_all = _threshold(maxgaps_all)
        m_fg = _threshold(maxmiss_fg)
        m_bg = _threshold(maxmiss_bg)
        m_all = _threshold(maxmiss_all)

        # Thresholds are 'equal' as in the scalar path: both sides capped at the same value,
        # or both uncapped with an overall cap.
        gap_pair_check = miss_pair and (
            (g_fg is not None and g_bg is not None and g_fg == g_bg)
            or (g_fg is None and g_bg is None and g_all is not None))
        miss_pair_check = miss_pair and (
            (m_fg is not None and m_bg is not None and m_fg == m_bg)
            or (m_fg is None and m_bg is None and m_all is not None))
        if gap_pair_check or miss_pair_check:
            self._ensure_pair_index()

        n_pos = len(positions_with_schemes)
        gapmat = np.zeros((self.n_sp, n_pos), dtype=np.float32)
        missmat = np.zeros((self.n_sp, n_pos), dtype=np.float32)
        position_names = []

        jobs = []
        blocks = []
        col_cursor = 0

        for pi, (pos_dict, schemes_set) in enumerate(positions_with_schemes):
            symbols, gapvec, missvec = self._position_symbols(pos_dict)
            gapmat[:, pi] = gapvec
            missmat[:, pi] = missvec

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
                scheme_names = [None]

            for sname in scheme_names:
                scheme_dict = None if sname is None else _SCHEME_MAP[sname]
                block = self._group_block(symbols, scheme_dict)
                ng = block.shape[1]
                blocks.append(block)
                jobs.append((pi, sname, col_cursor, col_cursor + ng))
                col_cursor += ng

        if caap_mode:
            results = {}
            for (pi, sname, _s, _e) in jobs:
                results[(position_names[pi], sname)] = 0
        else:
            results = {position_names[pi]: 0 for pi in range(n_pos)}

        hits = {key: [] for key in results} if collect_hits else None

        if self.B == 0 or col_cursor == 0:
            return (results, hits) if collect_hits else results

        G_concat = np.hstack(blocks) if blocks else np.zeros((self.n_sp, 0), np.float32)

        if b_chunk is None:
            b_chunk = self._default_b_chunk(col_cursor)

        for start in range(0, self.B, b_chunk):
            end = min(start + b_chunk, self.B)
            Fc = self.F[start:end]
            Bc = self.Bg[start:end]

            C_fg_all = Fc @ G_concat
            C_bg_all = Bc @ G_concat

            gaps_fg = Fc @ gapmat
            gaps_bg = Bc @ gapmat
            miss_fg = Fc @ missmat
            miss_bg = Bc @ missmat

            for (pi, sname, cs, ce) in jobs:
                C_fg = C_fg_all[:, cs:ce]
                C_bg = C_bg_all[:, cs:ce]

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
                if gap_pair_check or miss_pair_check:
                    # Sparse: only labelings that are already hits and have the tested
                    # condition on both sides can be discarded here.
                    for li in np.nonzero(hit)[0]:
                        b = start + int(li)
                        if gap_pair_check and gfg[li] > 0 and gbg[li] > 0 and \
                                self._pair_sets_differ(b, gapmat[:, pi]):
                            hit[li] = False
                        elif miss_pair_check and mfg[li] > 0 and mbg[li] > 0 and \
                                self._pair_sets_differ(b, missmat[:, pi]):
                            hit[li] = False
                count = int(hit.sum())

                key = (position_names[pi], sname) if caap_mode else position_names[pi]
                results[key] += count

                if collect_hits and count:
                    local_idx = np.nonzero(hit)[0]
                    hits[key].extend(self.alltraits[start + int(i)] for i in local_idx)

        return (results, hits) if collect_hits else results

    @staticmethod
    def _schemes_to_test(schemes_set):
        if not schemes_set:
            return list(_SCHEME_ORDER)
        if "CAAS" in schemes_set:
            return list(_SCHEME_ORDER)
        return [s for s in _SCHEME_ORDER if s in schemes_set]


def _posnum_from_posdict(pos_dict):
    for v in pos_dict.values():
        parts = v.split("@")
        if len(parts) > 1:
            return parts[1]
    return None


def _vectorized_position_counts(cfg, sliced_object, genename, positions_with_schemes,
                                max_fg_gaps, max_bg_gaps, max_overall_gaps,
                                max_fg_miss, max_bg_miss, max_overall_miss,
                                max_conserved, admitted_patterns, caap_mode,
                                collect_hits=False, miss_pair=False):
    """Run the vectorized kernel for one resample multiconfig.

    Returns the same per-(position[, scheme]) count mapping the scalar loop
    accumulates into position_counts; when collect_hits is True also returns the
    per-key list of hit labeling (trait) names for perm_discovery emission.
    """
    vb = VectorizedPermReplay(cfg, sliced_object.species)
    return vb.count(
        positions_with_schemes,
        genename,
        maxgaps_fg=max_fg_gaps, maxgaps_bg=max_bg_gaps, maxgaps_all=max_overall_gaps,
        maxmiss_fg=max_fg_miss, maxmiss_bg=max_bg_miss, maxmiss_all=max_overall_miss,
        max_conserved=max_conserved,
        admitted_patterns=admitted_patterns,
        caap_mode=caap_mode,
        collect_hits=collect_hits,
        miss_pair=miss_pair,
    )


def _emit_groups_rows(groups_out, genename, hits, caap_mode):
    """Materialize per-cycle groups debug rows from hit coordinates."""
    if not hits or not groups_out:
        return
    for key, trait_names in hits.items():
        if not trait_names:
            continue
        if caap_mode:
            position_name, scheme_name = key
        else:
            position_name, scheme_name = key, "US"
        posnum = position_name.split("@", 1)[1] if "@" in position_name else position_name
        for trait in trait_names:
            groups_out.write(f"{trait}\t{genename}\t{posnum}\tCAAP\t{scheme_name}\n")


def _emit_perm_discovery_rows(perm_discovery_out, cfg, genename, positions_with_schemes,
                              hits, caap_mode, max_conserved):
    """Materialize perm_discovery rows for the vectorized hits.

    The kernel identifies WHICH (position, scheme, labeling) are CAAS/CAAP; each hit's row
    is reconstructed with exact pair-ordered substitution strings and group encodings.
    """
    if not hits or not perm_discovery_out:
        return
    posname_to_posdict = {}
    for pos_dict, _schemes in positions_with_schemes:
        pn = genename + "@" + str(_posnum_from_posdict(pos_dict))
        posname_to_posdict[pn] = pos_dict

    def _ungapped_sorted(species_iter, pos_dict, trait=None):
        keep = [sp for sp in species_iter
                if sp in pos_dict and pos_dict[sp].split("@")[0].upper() not in _VEC_GAP_SYMBOLS]
        keep.sort(key=lambda sp: _pair_sort_key(cfg, sp, trait))
        return keep

    for key, trait_names in hits.items():
        if not trait_names:
            continue
        if caap_mode:
            position_name, scheme_name = key
            scheme_dict = _VEC_SCHEME_MAP.get(scheme_name, US)
        else:
            position_name, scheme_name, scheme_dict = key, "US", US
        pos_dict = posname_to_posdict.get(position_name)
        if not pos_dict:
            continue
        posnum = position_name.split("@", 1)[1]

        for trait in trait_names:
            fg = _ungapped_sorted(cfg.trait2fg.get(trait, ()), pos_dict, trait)
            bg = _ungapped_sorted(cfg.trait2bg.get(trait, ()), pos_dict, trait)
            fg_aas = "".join(pos_dict[sp].split("@")[0] for sp in fg)
            bg_aas = "".join(pos_dict[sp].split("@")[0] for sp in bg)

            is_match, pattern, substitution, conserved_pairs = check_pattern(
                fg_aas, bg_aas, scheme_dict=scheme_dict,
                max_conserved=max_conserved, multiconfig=cfg,
                fg_species_list=fg, bg_species_list=bg, trait=trait
            )
            encoded = encode_to_groups(fg_aas, scheme_dict) + "/" + encode_to_groups(bg_aas, scheme_dict)
            fields = [trait, genename, "CAAP", scheme_name, trait, str(posnum),
                      substitution, encoded, pattern]
            if max_conserved > 0:
                oc = conserved_pairs.split(":")[0] if conserved_pairs else "0"
                pl = conserved_pairs.split(":")[1] if conserved_pairs and ":" in conserved_pairs else ""
                fields.extend(["TRUE" if int(oc) > 0 else "FALSE", f"{oc}:{pl}"])
            perm_discovery_out.write("\t".join(fields) + "\n")


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
    Parse discovery output file and extract positions with their CAAP grouping schemes.
    Handles both classical CAAS and CAAP formats.
    
    Classical CAAS legacy format:   Gene\tMode\tTrait\tPosition\t...
    Normalized CAAS format:         Gene\tMode\tcaap_group\tTrait\tPosition\t...
    CAAP format:                    Gene\tMode\tcaap_group\tTrait\tPosition\t...
    
    Args:
        discovery_file: Path to discovery output file
        genename: Gene name to match (e.g., "BRCA1")
    
    Returns:
        dict mapping position -> set of grouping schemes found (e.g., {"142": {"GS1", "GS2"}})
        or None if no discovery file or all positions should be tested
    """
    if not exists(discovery_file):
        print(f"Warning: Discovery file {discovery_file} not found. Processing all positions with all schemes.")
        return None
    
    # position_number -> set of grouping schemes
    position_schemes = {}
    
    try:
        with open(discovery_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("gene"):  # Skip header
                    continue
                
                try:
                    fields = line.split('\t')
                    if len(fields) < 3:
                        continue
                    
                    gene = fields[0]
                    
                    # Match gene name first
                    if gene != genename:
                        continue
                    
                    # Detect format based on mode column
                    mode = fields[1] if len(fields) > 1 else ""
                    
                    if mode == "CAAP":
                        # CAAP format: Gene\tMode\tcaap_group\tTrait\tPosition\t...
                        # caap_group is at index 2, Position is at index 4
                        if len(fields) >= 5:
                            scheme = fields[2]  # US, GS1, GS2, GS3, or GS4
                            position = fields[4]
                            
                            if position not in position_schemes:
                                position_schemes[position] = set()
                            position_schemes[position].add(scheme)
                    elif mode == "CAAS":
                        # Gene mode caap_group Trait Position ...
                        if len(fields) >= 5 and fields[2] in {"US", "GS1", "GS2", "GS3", "GS4"}:
                            position = fields[4]
                        elif len(fields) >= 4:
                            position = fields[3]
                        else:
                            continue

                        # For classical CAAS, mark with special 'CAAS' marker to test all schemes
                        if position not in position_schemes:
                            position_schemes[position] = set()
                        position_schemes[position].add("CAAS")
                except:
                    continue
        
        return position_schemes if len(position_schemes) > 0 else None
    
    except Exception as e:
        print(f"Warning: Error parsing discovery file: {e}. Processing all positions with all schemes.")
        return None









# FUNCTION run_perm_replay_on_alignment()
# Launches perm-replay in several lines. Returns a dictionary gene@position --> pvalue

def run_perm_replay_on_alignment(trait_config_file, resampled_traits, sliced_object, max_fg_gaps, max_bg_gaps, max_overall_gaps, max_fg_miss, max_bg_miss, max_overall_miss, the_admitted_patterns, output_file, miss_pair=False, max_conserved=0, discovery_file=None, progress_log=None, caap_mode=False, export_groups=None, export_perm_discovery=None, fop_mode=False):
    """
    Run perm-replay on a single alignment.

    fop_mode (Gap A): resample source is the single fanned file fop_labelings.tab
    (cycle tags "<base>~H<m>"). Per-labeling CAAS hits are collapsed to base-cycle
    units before writing: a base cycle counts once iff ANY of its ~H<m> labelings
    is a CAAS at that (position, scheme). The output row format is unchanged
    (Gene@Position \\t scheme \\t hits \\t total \\t proportion) with hits/total in
    base-cycle units. Requires the single-file resample path (fop_labelings.tab is
    one file); a directory is transparently redirected to fop_labelings.tab inside
    it when present.
    
    Supports both single-file and directory-based resampled traits:
    - Single file: resampled_traits is a multicfg object loaded from one file
    - Directory: resampled_traits is the directory path (string), files loaded sequentially
    
    Args:
        resampled_traits: multicfg object OR directory path (str) containing resample_*.tab files
        progress_log: Optional file path for logging progress
        caap_mode: If True, test all CAAP grouping schemes (US, GS1-GS4) instead of classical CAAS
        ... (other parameters as before)
    """
    the_genename = sliced_object.genename

    # FOP (Gap A): fop_labelings.tab is a single file. If a directory was passed,
    # redirect to the fanned file inside it; degrade to a normal run if absent.
    if fop_mode:
        from modules.perm_replay_io import simtrait_revive
        if isinstance(resampled_traits, str) and os.path.isdir(resampled_traits):
            _fop_file = os.path.join(resampled_traits, "fop_labelings.tab")
            if os.path.exists(_fop_file):
                print(f"[FOP] base-cycle collapse over {_fop_file}")
                resampled_traits = simtrait_revive(_fop_file)
            else:
                print(f"[FOP] WARNING: {_fop_file} not found; running standard perm-replay")
                fop_mode = False
        elif isinstance(resampled_traits, str) and os.path.isfile(resampled_traits):
            print(f"[FOP] base-cycle collapse over {resampled_traits}")
            resampled_traits = simtrait_revive(resampled_traits)

    groups_handle = None
    if export_groups:
        groups_handle = open(export_groups, "w")
        groups_handle.write("Cycle\tGene\tPosition\tMode\tGroup\n")
    
    perm_discovery_handle = None
    if export_perm_discovery:
        perm_discovery_handle = open(export_perm_discovery, "w")
        header_fields = [
            "cycle",
            "gene",
            "mode",
            "caap_group",
            "trait",
            "position",
            "caas",
            "amino_encoded",
            "pattern"
        ]
        if max_conserved > 0:
            header_fields.extend(["is_conserved_meta", "conserved_pair"])
        perm_discovery_handle.write("\t".join(header_fields) + "\n")

    # Vectorized BLAS path.
    collect_hits = perm_discovery_handle is not None or groups_handle is not None or fop_mode

    try:
        # Detect if resampled_traits is a directory path or a multicfg object
        if isinstance(resampled_traits, str) and os.path.isdir(resampled_traits):
            # Directory mode: sequential processing
            print(f"\n{'='*80}")
            print(f"DIRECTORY-BASED PERM-REPLAY MODE")
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
            
            # Map position dict to (pos_dict, schemes_set) tuples
            positions_with_schemes = []
            
            if discovery_file:
                print(f"Filtering positions based on discovery results from: {discovery_file}")
                position_schemes = parse_discovery_positions(discovery_file, the_genename)
                
                if position_schemes:
                    for pos_dict in positions_list:
                        for species, aa_info in pos_dict.items():
                            pos_num = aa_info.split("@")[1]
                            if pos_num in position_schemes:
                                positions_with_schemes.append((pos_dict, position_schemes[pos_num]))
                                break
                    
                    speedup = total_positions / max(1, len(positions_with_schemes))
                    print(f"✓ Optimization: Testing only {len(positions_with_schemes)} CAAS positions (from {total_positions} total)")
                    print(f"✓ Speedup: {speedup:.1f}× fewer positions to test")
                    print(f"✓ Tests reduced: {total_positions * total_cycles:,} → {len(positions_with_schemes) * total_cycles:,}\n")
                else:
                    print("No CAAS positions found in discovery file. Processing all positions.\n")
                    positions_with_schemes = [(pos, None) for pos in positions_list]
            else:
                # No discovery file - test all positions with all schemes
                positions_with_schemes = [(pos, None) for pos in positions_list]
            
            # Process each file sequentially
            start_time = time.time()
            
            for file_idx, (file_path, file_config) in enumerate(simtrait_revive_from_dir(resample_dir), 1):
                file_start = time.time()
                
                log_progress(file_idx, resample_info['num_files'], start_time, 
                            log_file=progress_log, prefix=f"Processing file {os.path.basename(file_path)}")
                
                is_b0 = os.path.basename(file_path) == "resample_000.tab"

                # BLAS path: one set of batched matmuls for the whole file.
                if collect_hits:
                    file_counts, file_hits = _vectorized_position_counts(
                        file_config, sliced_object, the_genename, positions_with_schemes,
                        max_fg_gaps, max_bg_gaps, max_overall_gaps,
                        max_fg_miss, max_bg_miss, max_overall_miss,
                        max_conserved, the_admitted_patterns, caap_mode,
                        collect_hits=True, miss_pair=miss_pair,
                    )
                    if perm_discovery_handle is not None:
                        _emit_perm_discovery_rows(
                            perm_discovery_handle, file_config, the_genename,
                            positions_with_schemes, file_hits, caap_mode, max_conserved,
                        )
                    if groups_handle is not None:
                        _emit_groups_rows(groups_handle, the_genename, file_hits, caap_mode)
                else:
                    file_counts = _vectorized_position_counts(
                        file_config, sliced_object, the_genename, positions_with_schemes,
                        max_fg_gaps, max_bg_gaps, max_overall_gaps,
                        max_fg_miss, max_bg_miss, max_overall_miss,
                        max_conserved, the_admitted_patterns, caap_mode,
                        miss_pair=miss_pair,
                    )
                for key, count in file_counts.items():
                    position_counts[key] = position_counts.get(key, 0) + count

                file_elapsed = time.time() - file_start
                print(f"  → File completed in {format_time(file_elapsed)}\n")
            
            # Write final aggregated results
            print(f"\n{'='*80}")
            print(f"Writing aggregated results to {output_file}")
            print(f"{'='*80}\n")
            
            with open(output_file, "w") as ooout:
                if caap_mode:
                    # Sort by position then scheme
                    for key in sorted(position_counts.keys()):
                        position_name, scheme_name = key
                        count = position_counts[key]
                        empval = count / total_cycles
                        outline = "\t".join([position_name, scheme_name, str(count), str(total_cycles), str(empval)])
                        print(outline, file=ooout)
                else:
                    # Classical CAAS mode
                    for position_name in sorted(position_counts.keys()):
                        count = position_counts[position_name]
                        empval = count / total_cycles
                        outline = "\t".join([position_name, "US", str(count), str(total_cycles), str(empval)])
                        print(outline, file=ooout)
            
            total_elapsed = time.time() - start_time
            print(f"✓ Perm-replay complete in {format_time(total_elapsed)}")
            print(f"✓ Results written to {output_file}\n")
        
        else:
            # Single file mode (backward compatibility)
            resampled_traits_obj = resampled_traits
            print("caastools found", resampled_traits_obj.cycles, "resamplings")
            
            # OPTIMIZATION: Filter to only test positions found in discovery
            positions_list = list(sliced_object.d)
            total_positions = len(positions_list)
            
            # Map position dict to (pos_dict, schemes_set) tuples
            positions_with_schemes = []
            
            if discovery_file:
                print(f"\nFiltering positions based on discovery results from: {discovery_file}")
                position_schemes = parse_discovery_positions(discovery_file, the_genename)
                
                if position_schemes:
                    # Filter positions: only keep those that match discovery positions
                    for pos_dict in positions_list:
                        # Extract position number from the position dictionary
                        # Format is "AA@position_number" in values
                        for species, aa_info in pos_dict.items():
                            pos_num = aa_info.split("@")[1]
                            if pos_num in position_schemes:
                                positions_with_schemes.append((pos_dict, position_schemes[pos_num]))
                                break
                    
                    speedup = total_positions / max(1, len(positions_with_schemes))
                    print(f"✓ Optimization: Testing only {len(positions_with_schemes)} positions with discovered schemes (from {total_positions} total)")
                    print(f"✓ Speedup: {speedup:.1f}× fewer positions to test")
                    if caap_mode:
                        # Count how many scheme tests we're avoiding
                        total_scheme_tests = sum(len(schemes) if schemes and "CAAS" not in schemes else 6 for _, schemes in positions_with_schemes)
                        max_scheme_tests = len(positions_with_schemes) * 6  # 6 schemes max
                        scheme_speedup = max_scheme_tests / max(1, total_scheme_tests)
                        print(f"✓ Scheme optimization: Testing {total_scheme_tests} schemes (vs {max_scheme_tests} if testing all)")
                        print(f"✓ Scheme speedup: {scheme_speedup:.1f}× fewer scheme tests per position")
                    print(f"✓ Tests reduced: {total_positions * resampled_traits_obj.cycles:,} → {len(positions_with_schemes) * resampled_traits_obj.cycles:,}\n")
                else:
                    print("No positions found in discovery file. Processing all positions with all schemes.\n")
                    positions_with_schemes = [(pos, None) for pos in positions_list]
            else:
                # No discovery file - test all positions with all schemes
                positions_with_schemes = [(pos, None) for pos in positions_list]

            # BLAS path: batched matmuls over the whole resample object, then
            # emit the same per-(position[, scheme]) lines the scalar loop does,
            # plus perm_discovery rows and groups rows when requested.
            if collect_hits:
                counts, hits = _vectorized_position_counts(
                    resampled_traits_obj, sliced_object, the_genename, positions_with_schemes,
                    max_fg_gaps, max_bg_gaps, max_overall_gaps,
                    max_fg_miss, max_bg_miss, max_overall_miss,
                    max_conserved, the_admitted_patterns, caap_mode,
                    collect_hits=True, miss_pair=miss_pair,
                )
                if perm_discovery_handle is not None:
                    _emit_perm_discovery_rows(
                        perm_discovery_handle, resampled_traits_obj, the_genename,
                        positions_with_schemes, hits, caap_mode, max_conserved,
                    )
                if groups_handle is not None:
                    _emit_groups_rows(groups_handle, the_genename, hits, caap_mode)
            else:
                counts = _vectorized_position_counts(
                    resampled_traits_obj, sliced_object, the_genename, positions_with_schemes,
                    max_fg_gaps, max_bg_gaps, max_overall_gaps,
                    max_fg_miss, max_bg_miss, max_overall_miss,
                    max_conserved, the_admitted_patterns, caap_mode,
                    miss_pair=miss_pair,
                )
            if fop_mode:
                # Gap A: collapse per-labeling hits to base-cycle units.
                collapsed, cyc = collapse_fop_hits_by_base(hits, resampled_traits_obj.alltraits)
                counts = collapsed
                print(f"[FOP] {len(resampled_traits_obj.alltraits)} labelings -> {cyc} base cycles")
            else:
                cyc = resampled_traits_obj.cycles
            ooout = open(output_file, "w")
            if caap_mode:
                for (position_name, scheme_name), count in counts.items():
                    empval = str(count / cyc)
                    print("\t".join([position_name, scheme_name, str(count), str(cyc), empval]), file=ooout)
            else:
                for position_name, count in counts.items():
                    empval = str(count / cyc)
                    print("\t".join([position_name, "US", str(count), str(cyc), empval]), file=ooout)
            ooout.close()
            print(f"Results written to {output_file}")
    finally:
        if groups_handle:
            groups_handle.close()
        if perm_discovery_handle:
            perm_discovery_handle.close()

# FUNCTION pval()
# Returns a dictionary with the pvalue

def pval(perm_replay_result):
    with open(perm_replay_result) as h:
        thelist = h.read().splitlines()
    
    d = {}

    for line in thelist:
        try:
            c = line.split("\t")
            d[c[0]] = c[2]
        except:
            pass
    
    return d
