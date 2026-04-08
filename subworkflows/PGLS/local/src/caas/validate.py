#!/usr/bin/env python3
# ================================ src/caas/validate.py ================================

from __future__ import annotations

import glob
import os
import sys
from typing import Dict, Set, Tuple

import numpy as np
import pandas as pd
from Bio import AlignIO

try:
    from biochem.grouping import get_grouping_scheme  # type: ignore[import-not-found]
except Exception:
    def get_grouping_scheme(aa: str, scheme: str):
        aa_u = (aa or "").strip().upper()
        return aa_u if aa_u else None

AMBIGSYMS = {"-", "X", "B", "Z", "J", "U", "O"}


def parse_gene_pos(gene_pos: str) -> Tuple[str, int]:
    try:
        g, p = gene_pos.split("_")
        return g, int(p)
    except Exception as exc:
        raise ValueError(f"Invalid GenePos value: '{gene_pos}'") from exc


def parse_multiset_conv(new_convAA: str) -> Tuple[str, str]:
    if not isinstance(new_convAA, str) or "/" not in new_convAA:
        return "", ""
    left, right = new_convAA.split("/")
    return left, right


def top_bottom_letters(new_convAA: str) -> Tuple[Set[str], Set[str]]:
    """Return (top_letters, bottom_letters) excluding conserved pairs.

    Iterates positionally over left/right characters and skips any position
    where left[i] == right[i] (conserved pair), consistent with the caastools
    hypergeometric test logic.  Works correctly for both raw amino-acid strings
    (US scheme) and pre-grouped strings produced by CT_DISAMBIGUATION.
    """
    left, right = parse_multiset_conv(new_convAA)
    top: Set[str] = set()
    bottom: Set[str] = set()
    for la, ra in zip(left, right):
        if la == ra:
            continue  # conserved pair — exclude
        if la and la not in AMBIGSYMS:
            top.add(la)
        if ra and ra not in AMBIGSYMS:
            bottom.add(ra)
    return top, bottom


def top_bottom_letters_grouped(new_convAA: str, caap_group: str) -> Tuple[Set[str], Set[str]]:
    """Group-scheme-aware variant: maps each position to its group label before
    comparing, then excludes positions conserved *at the group level*.

    Used as a fallback when raw caas strings are provided for a non-US scheme
    (i.e. amino_encoded is absent).  For pre-grouped strings (amino_encoded
    already in group labels) use top_bottom_letters() directly.
    """
    left, right = parse_multiset_conv(new_convAA)
    top: Set[str] = set()
    bottom: Set[str] = set()
    for la, ra in zip(left, right):
        gl = _map_grouped_aa(la, caap_group)
        gr = _map_grouped_aa(ra, caap_group)
        if not gl or not gr:
            continue
        if gl == gr:
            continue  # conserved at group level — exclude
        top.add(gl)
        bottom.add(gr)
    return top, bottom


def _coalesce_change_side(values) -> str:
    sides = {str(v).strip().lower() for v in values if pd.notna(v)}
    if "both" in sides or ("top" in sides and "bottom" in sides):
        return "both"
    if "top" in sides:
        return "top"
    if "bottom" in sides:
        return "bottom"
    return "none"


def _map_grouped_aa(aa: str, caap_group: str) -> str:
    aa_u = (aa or "").strip().upper()
    if not aa_u or aa_u in AMBIGSYMS:
        return ""
    grouped = get_grouping_scheme(aa_u, caap_group or "US")
    if grouped is None:
        return aa_u
    return str(grouped)


def _map_letter_set_to_grouping(letters: Set[str], caap_group: str) -> Set[str]:
    mapped: Set[str] = set()
    for aa in letters:
        g = _map_grouped_aa(aa, caap_group)
        if g:
            mapped.add(g)
    return mapped


def _derive_letters_from_extremes(aln,
                                  pos: int,
                                  high_tax: Set[str],
                                  low_tax: Set[str],
                                  sp2tax: Dict[str, str] | None,
                                  caap_group: str) -> Tuple[Set[str], Set[str]]:
    top_letters: Set[str] = set()
    bottom_letters: Set[str] = set()
    for rec in aln:
        seq_id = str(rec.id).strip()
        tax = sp2tax.get(seq_id) if sp2tax else seq_id
        if tax is None:
            continue
        grouped = _map_grouped_aa(str(rec.seq[pos]), caap_group)
        if not grouped:
            continue
        if tax in high_tax:
            top_letters.add(grouped)
        elif tax in low_tax:
            bottom_letters.add(grouped)
    return top_letters, bottom_letters


def infer_alignment_format(path: str, requested_format: str) -> str:
    if requested_format and requested_format != "auto":
        return requested_format

    suffix = os.path.splitext(path)[1].lower()
    if suffix in {".fa", ".fasta"}:
        return "fasta"

    return "phylip-relaxed"


def open_alignment_for_gene(align_dir: str, gene: str, align_fmt: str, DEBUG: bool = False):
    candidates = sorted(glob.glob(os.path.join(align_dir, f"{gene}.*")))
    if not candidates:
        raise FileNotFoundError(f"No alignment found for gene '{gene}' in {align_dir}")
    last_err = None
    for fp in candidates:
        try:
            aln = AlignIO.read(fp, infer_alignment_format(fp, align_fmt))
            return aln, fp
        except Exception as e:
            last_err = e
            if DEBUG:
                sys.stderr.write(f"[DEBUG] Failed to read alignment '{fp}' for gene '{gene}': {e}\n")
    raise ValueError(f"Could not read alignment for gene '{gene}': {last_err}")


def trait_extremes(df: pd.DataFrame, species_col: str, trait: str,
                   q_lower: float, q_upper: float,
                   extremes_method: str = "quantile",
                   use_median_separation: bool = True,
                   verbose: bool = False) -> Tuple[Set[str], Set[str], float, float]:
    species_norm = df[species_col].astype(str).str.strip()
    vals = pd.to_numeric(df[trait], errors="coerce")
    valid = vals.notna()
    if not valid.any():
        raise ValueError(f"Trait '{trait}' has no valid numeric values")

    med = float(vals[valid].median())

    if extremes_method == "quantile":
        q1 = float(vals[valid].quantile(q_lower))
        q3 = float(vals[valid].quantile(q_upper))
        if use_median_separation:
            low_mask = valid & (vals <= q1) & (vals < med)
            high_mask = valid & (vals >= q3) & (vals > med)
        else:
            low_mask = valid & (vals <= q1)
            high_mask = valid & (vals >= q3)
        if verbose:
            sys.stderr.write(
                f"[INFO] Trait '{trait}' (quantile): N_valid={valid.sum()}, "
                f"Q{int(q_lower*100)}={q1:.6g}, median={med:.6g}, "
                f"Q{int(q_upper*100)}={q3:.6g}, LOW n={low_mask.sum()}, HIGH n={high_mask.sum()}\n"
            )
        return set(species_norm[low_mask]), set(species_norm[high_mask]), q1, q3

    if extremes_method in ["1sd", "2sd"]:
        sd_factor = 1.0 if extremes_method == "1sd" else 2.0
        std_dev = float(vals[valid].std())
        low_threshold = med - sd_factor * std_dev
        high_threshold = med + sd_factor * std_dev
        low_mask = valid & (vals <= low_threshold)
        high_mask = valid & (vals >= high_threshold)
        if verbose:
            sys.stderr.write(
                f"[INFO] Trait '{trait}' ({extremes_method}): N_valid={valid.sum()}, "
                f"median={med:.6g}, SD={std_dev:.6g}, "
                f"LOW<={low_threshold:.6g}, HIGH>={high_threshold:.6g}, "
                f"LOW n={low_mask.sum()}, HIGH n={high_mask.sum()}\n"
            )
        return set(species_norm[low_mask]), set(species_norm[high_mask]), low_threshold, high_threshold

    raise ValueError(f"Unknown extremes_method: {extremes_method}")


def extremes_DEBUG_table(df: pd.DataFrame, species_col: str, trait: str,
                         q_lower: float, q_upper: float,
                         extremes_method: str = "quantile",
                         use_median_separation: bool = True) -> pd.DataFrame:
    vals = pd.to_numeric(df[trait], errors="coerce")
    valid = vals.notna()
    med = float(vals[valid].median())

    if extremes_method == "quantile":
        q1 = float(vals[valid].quantile(q_lower))
        q3 = float(vals[valid].quantile(q_upper))
        if use_median_separation:
            low_mask = valid & (vals <= q1) & (vals < med)
            high_mask = valid & (vals >= q3) & (vals > med)
        else:
            low_mask = valid & (vals <= q1)
            high_mask = valid & (vals >= q3)
        return pd.DataFrame({
            "species": df[species_col].astype(str).str.strip(),
            "trait": trait,
            "trait_value": vals,
            "is_valid": valid,
            "is_low": low_mask,
            "is_high": high_mask,
            "q1": q1,
            "median": med,
            "q3": q3,
        })

    if extremes_method in ["1sd", "2sd"]:
        sd_factor = 1.0 if extremes_method == "1sd" else 2.0
        std_dev = float(vals[valid].std())
        low_threshold = med - sd_factor * std_dev
        high_threshold = med + sd_factor * std_dev
        low_mask = valid & (vals <= low_threshold)
        high_mask = valid & (vals >= high_threshold)
        return pd.DataFrame({
            "species": df[species_col].astype(str).str.strip(),
            "trait": trait,
            "trait_value": vals,
            "is_valid": valid,
            "is_low": low_mask,
            "is_high": high_mask,
            "low_threshold": low_threshold,
            "median": med,
            "high_threshold": high_threshold,
        })

    raise ValueError(f"Unknown extremes_method: {extremes_method}")


def build_validation_table(caas_df: pd.DataFrame,
                           traits_df: pd.DataFrame,
                           species_col: str,
                           trait_name: str,
                           align_dir: str,
                           align_fmt: str,
                           q_lower: float,
                           q_upper: float,
                           keep_ambiguous: bool,
                           include_other: bool,
                           position_mode: str = "caas_only",
                           extremes_method: str = "quantile",
                           use_median_separation: bool = True,
                           DEBUG: bool = False,
                           sp2tax: Dict[str, str] | None = None,
                           tax2sp: Dict[str, str] | None = None,
                           taxid2names: Dict[str, set] | None = None,
                           caas_sp: Set[str] | None = None,
                           caas_tax_ids: Set[str] | None = None,
                           allowed_tax_ids: Set[str] | None = None) -> pd.DataFrame:
    aln_cache: Dict[str, tuple] = {}
    out_rows = []

    required_cols = ["Gene", "Position"]
    for col in required_cols:
        if col not in caas_df.columns:
            raise ValueError(f"CAAS file missing column '{col}'")
    if position_mode not in {"caas_only", "all_positions_for_caas_genes"}:
        raise ValueError("position_mode must be one of: caas_only, all_positions_for_caas_genes")

    if "caas" not in caas_df.columns:
        caas_df = caas_df.copy()
        caas_df["caas"] = ""
    if "change_side" not in caas_df.columns:
        caas_df = caas_df.copy()
        caas_df["change_side"] = "none"
    if "CAAP_Group" not in caas_df.columns:
        caas_df = caas_df.copy()
        caas_df["CAAP_Group"] = "US"

    low_sp, high_sp, q1, q3 = trait_extremes(
        traits_df,
        species_col,
        trait_name,
        q_lower,
        q_upper,
        extremes_method=extremes_method,
        use_median_separation=use_median_separation,
        verbose=DEBUG,
    )
    low_tax = {sp2tax[s] for s in low_sp if sp2tax and s in sp2tax} if sp2tax else set(low_sp)
    high_tax = {sp2tax[s] for s in high_sp if sp2tax and s in sp2tax} if sp2tax else set(high_sp)
    # Use the caller-supplied tax2sp (scientific-name-preferring) if available;
    # fall back to naive inversion only when no tax-map file was provided.
    if tax2sp is None:
        tax2sp = {tax: sp for sp, tax in (sp2tax or {}).items()}
    if allowed_tax_ids is not None:
        low_tax &= set(allowed_tax_ids)
        high_tax &= set(allowed_tax_ids)
    # Exclusion set: species used for the CAAS test should not appear in the
    # internal validation. Prefer the explicit CAAS traitfile set when provided;
    # fall back to all species in the PGLS traitfile.
    if caas_sp is not None:
        all_trait_sp_names = caas_sp
        all_trait_tax = caas_tax_ids if caas_tax_ids is not None else (
            {sp2tax[s] for s in caas_sp if sp2tax and s in sp2tax} if sp2tax else caas_sp
        )
    else:
        all_trait_sp_names = set(traits_df[species_col].astype(str).str.strip())
        all_trait_tax = (
            {sp2tax[s] for s in all_trait_sp_names if s in sp2tax}
            if sp2tax else all_trait_sp_names
        )

    s2v = dict(zip(traits_df[species_col].astype(str).str.strip(), pd.to_numeric(traits_df[trait_name], errors="coerce")))

    def _emit_internal_rows(gene: str,
                           pos: int,
                           aln,
                           top_letters: Set[str],
                           bottom_letters: Set[str],
                           substitution: str,
                           change_side: str,
                           caap_group: str,
                           is_caas_position: bool):
        for rec in aln:
            seq_id = str(rec.id).strip()
            # Resolve seq_id → tax_id: direct lookup first, then synonym reverse-lookup
            tax = None
            if sp2tax:
                tax = sp2tax.get(seq_id)
                if tax is None and taxid2names:
                    for tid, names in taxid2names.items():
                        if seq_id in names:
                            tax = tid
                            break
            else:
                tax = seq_id
            # Exclude all species in the trait file (they defined the CAAS test set)
            if tax is not None:
                if tax in all_trait_tax:
                    continue
            else:
                if seq_id in all_trait_sp_names:
                    continue
            aa_raw = str(rec.seq[pos])
            aa_obs = _map_grouped_aa(aa_raw, caap_group)
            if aa_obs in top_letters and aa_obs:
                odir = "top"
            elif aa_obs in bottom_letters and aa_obs:
                odir = "bottom"
            else:
                if not include_other:
                    continue
                odir = "other"
            # Resolve to canonical (scientific) species name using tax_id as the key
            if tax2sp and tax:
                species_name = tax2sp.get(tax, seq_id)
            elif tax2sp:
                species_name = tax2sp.get(seq_id, seq_id)
            else:
                species_name = seq_id
            out_rows.append({
                "Gene": gene,
                "Position": pos,
                "Substitution": substitution,
                "change_side": change_side,
                "CAAP_Group": caap_group,
                "is_caas_position": bool(is_caas_position),
                "species": species_name,
                "observed_AA_raw": aa_raw,
                "observed_AA": aa_obs,
                "observed_dir": odir,
                "trait": trait_name,
                "trait_value": float(s2v.get(species_name, np.nan)),
                "extreme_side": "INTERNAL",
                "q1": q1,
                "q3": q3,
            })

    caas_work = caas_df.copy()
    caas_work["Gene"] = caas_work["Gene"].astype(str).str.strip()
    caas_work["CAAP_Group"] = caas_work["CAAP_Group"].astype(str).str.strip().replace("", "US")
    caas_work["change_side"] = caas_work["change_side"].astype(str).str.strip().str.lower().replace("", "none")
    caas_work["Position"] = pd.to_numeric(caas_work["Position"], errors="coerce")
    caas_work = caas_work.dropna(subset=["Position"])
    caas_work["Position"] = caas_work["Position"].astype(int)

    if position_mode == "caas_only":
        for _, row in caas_work.iterrows():
            gene = row["Gene"]
            pos = int(row["Position"])
            caap_group = str(row.get("CAAP_Group", "US") or "US")
            site_id = f"{gene}_{pos}"
            change_side = str(row.get("change_side", "none")).strip().lower()
            if not keep_ambiguous and change_side not in {"top", "bottom"}:
                if DEBUG:
                    sys.stderr.write(f"[DEBUG] Skip ambiguous change_side '{change_side}' for {site_id}\n")
                continue

            if gene not in aln_cache:
                try:
                    aln_cache[gene] = open_alignment_for_gene(align_dir, gene, align_fmt, DEBUG=DEBUG)
                except Exception as e:
                    sys.stderr.write(f"[WARN] {e}\n")
                    aln_cache[gene] = (None, None)
            aln_entry = aln_cache.get(gene)
            if not aln_entry:
                continue
            aln = aln_entry[0]
            if aln is None:
                continue
            if pos < 0 or pos >= aln.get_alignment_length():
                sys.stderr.write(f"[WARN] Position {pos} out of bounds for gene {gene} (len={aln.get_alignment_length()})\n")
                continue

            amino_enc = str(row.get("amino_encoded", "") or "").strip()
            if caap_group != "US" and amino_enc and "/" in amino_enc:
                # amino_encoded already carries the group-label substitution string
                # as produced by CT_DISAMBIGUATION (e.g. "aaa/ass" for GS2).
                # Use it directly so the symbol space matches observed_AA exactly.
                top_letters, bottom_letters = top_bottom_letters(amino_enc)
            else:
                caas_str = str(row.get("caas", ""))
                if caap_group != "US":
                    top_letters, bottom_letters = top_bottom_letters_grouped(caas_str, caap_group)
                else:
                    top_letters, bottom_letters = top_bottom_letters(caas_str)
            if not top_letters or not bottom_letters:
                if DEBUG:
                    sys.stderr.write(f"[DEBUG] Skip {site_id} due to empty TOP/BOTTOM amino-acid sets\n")
                continue

            _emit_internal_rows(
                gene=gene,
                pos=pos,
                aln=aln,
                top_letters=top_letters,
                bottom_letters=bottom_letters,
                substitution=str(row.get("caas", "")),
                change_side=change_side,
                caap_group=caap_group,
                is_caas_position=True,
            )
    else:
        for (gene_raw, caap_group_raw), group_df in caas_work.groupby(["Gene", "CAAP_Group"], dropna=False):
            gene = str(gene_raw)
            caap_group = str(caap_group_raw or "US")
            if gene not in aln_cache:
                try:
                    aln_cache[gene] = open_alignment_for_gene(align_dir, gene, align_fmt, DEBUG=DEBUG)
                except Exception as e:
                    sys.stderr.write(f"[WARN] {e}\n")
                    aln_cache[gene] = (None, None)
            aln_entry = aln_cache.get(gene)
            if not aln_entry:
                continue
            aln = aln_entry[0]
            if aln is None:
                continue

            caas_pos_meta = {}
            for pos, g in group_df.groupby("Position"):
                caas_pos_meta[int(str(pos))] = {
                    "is_caas_position": True,
                    "change_side": _coalesce_change_side(g["change_side"].tolist()),
                    "substitution": str(g["caas"].astype(str).iloc[0]) if "caas" in g.columns else "",
                }

            aln_len = int(aln.get_alignment_length())
            for pos in range(aln_len):
                meta = caas_pos_meta.get(pos, None)
                is_caas_position = bool(meta is not None)
                change_side = (meta or {}).get("change_side", "none")
                if is_caas_position and (not keep_ambiguous) and change_side not in {"top", "bottom"}:
                    continue

                top_letters, bottom_letters = _derive_letters_from_extremes(
                    aln=aln,
                    pos=pos,
                    high_tax=high_tax,
                    low_tax=low_tax,
                    sp2tax=sp2tax,
                    caap_group=caap_group,
                )
                if not top_letters or not bottom_letters:
                    continue

                substitution = (meta or {}).get(
                    "substitution",
                    f"{''.join(sorted(top_letters))}/{''.join(sorted(bottom_letters))}",
                )
                _emit_internal_rows(
                    gene=gene,
                    pos=pos,
                    aln=aln,
                    top_letters=top_letters,
                    bottom_letters=bottom_letters,
                    substitution=substitution,
                    change_side=change_side,
                    caap_group=str(caap_group or "US"),
                    is_caas_position=is_caas_position,
                )

    cols = [
        "Gene", "Position", "Substitution", "change_side", "CAAP_Group", "is_caas_position",
        "species", "observed_AA_raw", "observed_AA", "observed_dir",
        "trait", "trait_value", "extreme_side", "q1", "q3"
    ]
    if not out_rows:
        sys.stderr.write("[WARN] No rows generated for validation table.\n")
        return pd.DataFrame(columns=cols)
    return pd.DataFrame(out_rows, columns=cols)
