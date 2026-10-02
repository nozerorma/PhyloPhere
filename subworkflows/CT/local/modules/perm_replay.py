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
from modules.caas_id import (
    US, GS1, GS2, GS3, GS4, SCHEMES,
    check_pattern, check_caap_pattern, iscaas,
    encode_to_groups, _pair_sort_key, process_position
)
from modules.pindex import load_cfg
from modules.alimport import *

import os
import re
import numpy as np

# ---------------------------------------------------------------------------
# Labeling tags. Under params.multi_hypothesis the labelings file holds "<base>~H<m>" tags: one row per cycle and
# Dunn-independent hypothesis. The base cycle of a tag is the tag without its "~H<m>" suffix.
# ---------------------------------------------------------------------------
_FOP_H_SUFFIX = re.compile(r"~H\d+$")


def _fop_base_cycle(labeling_tag):
    """'b_12~H3' -> 'b_12'. A tag with no ~H<m> suffix is its own base cycle."""
    return _FOP_H_SUFFIX.sub("", labeling_tag)


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
        """Return (symbols, gapvec, missvec, ndvec) over the species universe for a position.

        ndvec flags species carrying any symbol other than '-' (ambiguity codes count): the
        scalar path registers a trait's fg/bg side (trait2aas_fg/bg) only through those.
        """
        symbols = [None] * self.n_sp
        gapvec = np.zeros(self.n_sp, dtype=np.float32)
        missvec = np.zeros(self.n_sp, dtype=np.float32)
        ndvec = np.zeros(self.n_sp, dtype=np.float32)

        present = np.zeros(self.n_sp, dtype=bool)
        for sp, val in pos_dict.items():
            j = self.sp_index.get(sp)
            if j is None:
                continue
            present[j] = True
            aa = val.split("@")[0].upper()
            if aa != "-":
                ndvec[j] = 1.0
            if aa in _GAP_SYMBOLS:
                gapvec[j] = 1.0
            else:
                symbols[j] = aa

        missvec[~present] = 1.0
        return symbols, gapvec, missvec, ndvec

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

    def _pair_ok(self, b, pi, gap_check, miss_check, gapmat, missmat, gfg, gbg, mfg, mbg, li):
        """miss_pair verdict for labeling b at position pi: False = the scalar path discards it."""
        if gap_check and gfg[li] > 0 and gbg[li] > 0 and self._pair_sets_differ(b, gapmat[:, pi]):
            return False
        if miss_check and mfg[li] > 0 and mbg[li] > 0 and self._pair_sets_differ(b, missmat[:, pi]):
            return False
        return True

    def _default_b_chunk(self, total_groups):
        bytes_per_row = max(1, 2 * (self.n_sp + total_groups) * 4)
        budget_bytes = _CHUNK_MEM_BUDGET_MB * 1024 * 1024
        return max(1, min(self.B, int(budget_bytes // bytes_per_row)))

    # -- main counting --------------------------------------------------------

    def count(self, positions_with_schemes, genename,
              maxgaps_fg, maxgaps_bg, maxgaps_all,
              maxmiss_fg, maxmiss_bg, maxmiss_all,
              max_conserved, admitted_patterns, caap_mode,
              b_chunk=None, collect_hits=False, miss_pair=False,
              background_sink=None, background_base="b_0"):
        """Count CAAS/CAAP hits per (position[, scheme]) across all B labelings.

        miss_pair mirrors the scalar path (caas_id.fetch_caas): when the fg and bg
        thresholds are equal, a labeling whose fg and bg sides both have gapped (or missing)
        species, but in different pairs, is discarded.

        background_sink (a set) collects the positions 'tested' by the labelings whose base cycle is
        `background_base`, as the scalar path's background.output does: a position is tested when at
        least one such labeling passes the gap/missing thresholds, has a symbol other than '-' on both
        sides and survives miss_pair, whether or not it yields a CAAS.
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
        ndmat = np.zeros((self.n_sp, n_pos), dtype=np.float32)
        position_names = []

        jobs = []
        blocks = []
        col_cursor = 0

        for pi, (pos_dict, schemes_set) in enumerate(positions_with_schemes):
            symbols, gapvec, missvec, ndvec = self._position_symbols(pos_dict)
            gapmat[:, pi] = gapvec
            missmat[:, pi] = missvec
            ndmat[:, pi] = ndvec

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

        track = None
        tested_pi = set()
        if background_sink is not None:
            track = np.array([_fop_base_cycle(t) == background_base for t in self.alltraits], dtype=bool)

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
            has_fg = (Fc @ ndmat) > 0
            has_bg = (Bc @ ndmat) > 0

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

                if track is not None and pi not in tested_pi:
                    for li in np.nonzero(valid & has_fg[:, pi] & has_bg[:, pi] & track[start:end])[0]:
                        if self._pair_ok(start + int(li), pi, gap_pair_check, miss_pair_check,
                                         gapmat, missmat, gfg, gbg, mfg, mbg, li):
                            tested_pi.add(pi)
                            break

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
                    # Sparse: only labelings that are already hits can be discarded here.
                    for li in np.nonzero(hit)[0]:
                        if not self._pair_ok(start + int(li), pi, gap_pair_check, miss_pair_check,
                                             gapmat, missmat, gfg, gbg, mfg, mbg, li):
                            hit[li] = False
                count = int(hit.sum())

                key = (position_names[pi], sname) if caap_mode else position_names[pi]
                results[key] += count

                if collect_hits and count:
                    local_idx = np.nonzero(hit)[0]
                    hits[key].extend(self.alltraits[start + int(i)] for i in local_idx)

        if background_sink is not None:
            background_sink.update(position_names[pi].split("@", 1)[1] for pi in tested_pi)

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
                                collect_hits=False, miss_pair=False, background_sink=None):
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
        background_sink=background_sink,
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


def _ungapped_sorted(cfg, species_iter, pos_dict, trait=None):
    """The species of `species_iter` with a resolved residue at this position, in pair order."""
    keep = [sp for sp in species_iter
            if sp in pos_dict and pos_dict[sp].split("@")[0].upper() not in _VEC_GAP_SYMBOLS]
    keep.sort(key=lambda sp: _pair_sort_key(cfg, sp, trait))
    return keep


def _reconstruct_hit(cfg, pos_dict, scheme_dict, trait, max_conserved):
    """One kernel hit with exact pair-ordered residues: the species lists and what check_pattern says.

    check_pattern is the scalar path's own verdict on the same fg/bg residues; `is_match` False means the
    kernel and the scalar CAAS test disagree.
    """
    fg = _ungapped_sorted(cfg, cfg.trait2fg.get(trait, ()), pos_dict, trait)
    bg = _ungapped_sorted(cfg, cfg.trait2bg.get(trait, ()), pos_dict, trait)
    fg_aas = "".join(pos_dict[sp].split("@")[0] for sp in fg)
    bg_aas = "".join(pos_dict[sp].split("@")[0] for sp in bg)
    is_match, pattern, substitution, conserved_pairs = check_pattern(
        fg_aas, bg_aas, scheme_dict=scheme_dict,
        max_conserved=max_conserved, multiconfig=cfg,
        fg_species_list=fg, bg_species_list=bg, trait=trait
    )
    encoded = encode_to_groups(fg_aas, scheme_dict) + "/" + encode_to_groups(bg_aas, scheme_dict)
    return dict(fg=fg, bg=bg, is_match=is_match, pattern=pattern, substitution=substitution,
                conserved_pairs=conserved_pairs, encoded=encoded)


def _conserved_fields(conserved_pairs):
    """[is_conserved_meta, conserved_pair] of a row, from check_pattern's '{count}:{pair ids}'."""
    oc = conserved_pairs.split(":")[0] if conserved_pairs else "0"
    pl = conserved_pairs.split(":")[1] if conserved_pairs and ":" in conserved_pairs else ""
    return ["TRUE" if int(oc) > 0 else "FALSE", f"{oc}:{pl}"]


def _hits_by_position(positions_with_schemes, genename):
    return {genename + "@" + str(_posnum_from_posdict(pos_dict)): pos_dict for pos_dict, _schemes in positions_with_schemes}


def _emit_perm_discovery_rows(perm_discovery_out, cfg, genename, positions_with_schemes,
                              hits, caap_mode, max_conserved):
    """Materialize perm_discovery rows for the vectorized hits.

    The kernel identifies WHICH (position, scheme, labeling) are CAAS/CAAP; each hit's row
    is reconstructed with exact pair-ordered substitution strings and group encodings.

    A hit check_pattern rejects means the kernel and the scalar CAAS test disagree; that must not pass
    silently, so it is counted and raised once every row of the call has been written.
    """
    if not hits or not perm_discovery_out:
        return
    disagreements = []
    posname_to_posdict = _hits_by_position(positions_with_schemes, genename)

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
            hit = _reconstruct_hit(cfg, pos_dict, scheme_dict, trait, max_conserved)
            if not hit["is_match"]:
                disagreements.append((trait, genename, scheme_name, posnum))
            fields = [trait, genename, "CAAP", scheme_name, trait, str(posnum),
                      hit["substitution"], hit["encoded"], hit["pattern"]]
            if max_conserved > 0:
                fields.extend(_conserved_fields(hit["conserved_pairs"]))
            perm_discovery_out.write("\t".join(fields) + "\n")

    if disagreements:
        raise RuntimeError(
            f"vectorized kernel and check_pattern disagree on {len(disagreements)} hit(s) "
            f"(labeling, gene, scheme, position), first: {disagreements[:5]}")


_B0_DISCOVERY_HEADER = ["gene", "mode", "caap_group", "trait", "position", "caas", "amino_encoded", "pattern",
                        "ffgn", "fbgn", "gfg", "gbg", "mfg", "mbg", "ffg", "fbg", "ms"]
_HYPOTHESIS = re.compile(r"H\d+")


def _hypothesis_of(name):
    """'b_0~H3' / 'traitfile_H3.tab' -> 'H3'; no hypothesis token (a single contrast) -> 'H1'."""
    m = _HYPOTHESIS.search(str(name).split("~", 1)[1] if "~" in str(name) else str(name))
    return m.group(0) if m else "H1"


def _b0_discovery_rows(design_cfg, labeling_cfg, sliced_object, genename, positions_with_schemes,
                       hits, caap_mode, max_conserved):
    """The rows discovery.tab holds for the b_0 hits of the kernel, as lists of fields.

    The kernel says which (position, hypothesis, scheme) are CAAS. caas, amino_encoded and pattern come from
    the same reconstruction as the perm-discovery export; the other columns (species counts, gaps, missing,
    species lists) are `process_position`'s, the scalar definitions, on the real design: its trait names
    (traitfile_H3.tab, not b_0~H3) and its missing species. Rows are in a fixed order: position, then trait by
    file name, then scheme (the scalar's trait order is the glob order of the trait directory, so it changes
    from one filesystem to another). `ms` lists the missing species in pair order, foreground first
    (the scalar builds it from a set, so its order is arbitrary).
    """
    if not hits:
        return []
    design_trait = {}
    for name in design_cfg.alltraits:
        design_trait.setdefault(_hypothesis_of(name), name)
    trait_rank = {name: i for i, name in enumerate(sorted(design_cfg.alltraits))}
    scheme_rank = {name: i for i, name in enumerate(SCHEMES)}
    posname_to_posdict = _hits_by_position(positions_with_schemes, genename)
    processed = {}
    rows = []
    disagreements = []

    for key, labelings in hits.items():
        tags = [t for t in labelings if _fop_base_cycle(t) == "b_0"]
        if not tags:
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
        if position_name not in processed:
            processed[position_name] = process_position(pos_dict, design_cfg, sliced_object.species)
        z = processed[position_name]

        for tag in tags:
            trait = design_trait.get(_hypothesis_of(tag))
            if trait is None:
                raise RuntimeError(f"no trait of the design matches the b_0 labeling {tag!r}: {sorted(design_trait)[:5]}")
            hit = _reconstruct_hit(labeling_cfg, pos_dict, scheme_dict, tag, max_conserved)
            order = lambda sp: _pair_sort_key(design_cfg, sp, trait)
            fg = sorted(z.trait2ungapped_fg[trait], key=order)
            bg = sorted(z.trait2ungapped_bg[trait], key=order)
            if not hit["is_match"] or set(fg) != set(hit["fg"]) or set(bg) != set(hit["bg"]):
                disagreements.append((tag, genename, scheme_name, posnum))
            missing = set(z.trait2missings[trait])
            ms = (sorted((sp for sp in design_cfg.trait2fg[trait] if sp in missing), key=order)
                  + sorted((sp for sp in design_cfg.trait2bg[trait] if sp in missing), key=order))
            fields = [genename, "CAAP", scheme_name, trait, str(posnum), hit["substitution"], hit["encoded"], hit["pattern"],
                      str(len(fg)), str(len(bg)),
                      str(z.trait2gaps_fg.get(trait, 0)), str(z.trait2gaps_bg.get(trait, 0)),
                      str(z.trait2miss_fg.get(trait, 0)), str(z.trait2miss_bg.get(trait, 0)),
                      ",".join(fg) if fg else "NA", ",".join(bg) if bg else "NA", ",".join(ms) if ms else "NA"]
            if max_conserved > 0:
                fields.extend(_conserved_fields(hit["conserved_pairs"]))
            rows.append(((int(posnum), trait_rank[trait], scheme_rank[scheme_name]), fields))

    if disagreements:
        raise RuntimeError(
            f"vectorized kernel and the scalar definitions disagree on {len(disagreements)} b_0 hit(s) "
            f"(labeling, gene, scheme, position), first: {disagreements[:5]}")
    rows.sort(key=lambda r: r[0])
    return [fields for _key, fields in rows]


def _write_b0_discovery(path, rows, max_conserved):
    """A header and the rows, and no file at all when there is no row."""
    if not rows:
        return
    header = _B0_DISCOVERY_HEADER + (["is_conserved_meta", "conserved_pair"] if max_conserved > 0 else [])
    with open(path, "w") as out:
        out.write("\t".join(header) + "\n")
        for fields in rows:
            out.write("\t".join(fields) + "\n")


def run_perm_replay_on_alignment(trait_config_file, resampled_traits, sliced_object, max_fg_gaps, max_bg_gaps, max_overall_gaps, max_fg_miss, max_bg_miss, max_overall_miss, the_admitted_patterns, miss_pair=False, max_conserved=0, caap_mode=False, export_groups=None, export_perm_discovery=None, export_b0_background=None, export_b0_discovery=None):
    """
    Replay the labelings of one file over one alignment and write the requested exports.

    `resampled_traits` is the multicfg object of the labelings file (core `simtrait_revive`): the real labeling b_0 and the
    permuted ones, any of them "<base>~H<m>" hypothesis labelings. At least one export is required:

      export_perm_discovery   the CAAS found under every labeling (the columns the permulation null reads)
      export_groups           the foreground and background groups of every labeling that has a hit
      export_b0_discovery     the discovery.tab rows of the b_0 labeling(s), written only when there is a row
      export_b0_background    '<gene>\t<positions tested by b_0 | NULL>'

    caap_mode tests the five grouping schemes (US, GS1-GS4) at every position; otherwise only the ungrouped one.
    """
    if not (export_groups or export_perm_discovery or export_b0_discovery or export_b0_background):
        raise ValueError("perm-replay needs at least one export (perm discovery, groups, b_0 discovery or b_0 background)")
    if not hasattr(resampled_traits, "cycles"):
        raise ValueError("perm-replay takes the labelings of one file (simtrait_revive), not a directory")
    the_genename = sliced_object.genename

    # positions tested by the b_0 labelings (see VectorizedPermReplay.count), one line per gene
    bg_sink = set() if export_b0_background else None

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

    collect_hits = perm_discovery_handle is not None or groups_handle is not None or bool(export_b0_discovery)

    try:
        print("caastools found", resampled_traits.cycles, "resamplings")

        # every position is tested with every scheme
        positions_with_schemes = [(pos, None) for pos in sliced_object.d]

        # Batched matmuls over the whole labelings object, then the requested exports.
        if collect_hits:
            _counts, hits = _vectorized_position_counts(
                resampled_traits, sliced_object, the_genename, positions_with_schemes,
                max_fg_gaps, max_bg_gaps, max_overall_gaps,
                max_fg_miss, max_bg_miss, max_overall_miss,
                max_conserved, the_admitted_patterns, caap_mode,
                collect_hits=True, miss_pair=miss_pair, background_sink=bg_sink,
            )
            if perm_discovery_handle is not None:
                _emit_perm_discovery_rows(
                    perm_discovery_handle, resampled_traits, the_genename,
                    positions_with_schemes, hits, caap_mode, max_conserved,
                )
            if groups_handle is not None:
                _emit_groups_rows(groups_handle, the_genename, hits, caap_mode)
            if export_b0_discovery:
                _write_b0_discovery(export_b0_discovery, _b0_discovery_rows(
                    load_cfg(trait_config_file), resampled_traits, sliced_object, the_genename,
                    positions_with_schemes, hits, caap_mode, max_conserved), max_conserved)
        else:
            _vectorized_position_counts(
                resampled_traits, sliced_object, the_genename, positions_with_schemes,
                max_fg_gaps, max_bg_gaps, max_overall_gaps,
                max_fg_miss, max_bg_miss, max_overall_miss,
                max_conserved, the_admitted_patterns, caap_mode,
                miss_pair=miss_pair, background_sink=bg_sink,
            )
    finally:
        if groups_handle:
            groups_handle.close()
        if perm_discovery_handle:
            perm_discovery_handle.close()

    if export_b0_background:
        tested = ",".join(sorted(bg_sink, key=int)) if bg_sink else "NULL"
        with open(export_b0_background, "w") as bkg:
            bkg.write(f"{the_genename}\t{tested}\n")
