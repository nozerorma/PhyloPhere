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
    left, right = parse_multiset_conv(new_convAA)
    top = {aa for aa in left if aa and aa not in AMBIGSYMS}
    bottom = {aa for aa in right if aa and aa not in AMBIGSYMS}
    return top, bottom


def open_alignment_for_gene(align_dir: str, gene: str, align_fmt: str, DEBUG: bool = False):
    candidates = sorted(glob.glob(os.path.join(align_dir, f"{gene}.*")))
    if not candidates:
        raise FileNotFoundError(f"No alignment found for gene '{gene}' in {align_dir}")
    last_err = None
    for fp in candidates:
        try:
            aln = AlignIO.read(fp, align_fmt)
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
                           extremes_method: str = "quantile",
                           use_median_separation: bool = True,
                           DEBUG: bool = False,
                           sp2tax: Dict[str, str] | None = None,
                           allowed_tax_ids: Set[str] | None = None) -> pd.DataFrame:
    aln_cache: Dict[str, tuple] = {}
    out_rows = []

    required_cols = ["Gene", "Position", "caas", "change_side"]
    for col in required_cols:
        if col not in caas_df.columns:
            raise ValueError(f"CAAS file missing column '{col}'")

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
    tax2sp = {tax: sp for sp, tax in (sp2tax or {}).items()}
    if allowed_tax_ids is not None:
        low_tax &= set(allowed_tax_ids)
        high_tax &= set(allowed_tax_ids)
    internal_sp = set(traits_df[species_col].astype(str).str.strip()) - low_sp - high_sp
    internal_tax = {sp2tax[s] for s in internal_sp if sp2tax and s in sp2tax} if sp2tax else set(internal_sp)
    if allowed_tax_ids is not None:
        internal_tax &= set(allowed_tax_ids)

    s2v = dict(zip(traits_df[species_col].astype(str), pd.to_numeric(traits_df[trait_name], errors="coerce")))

    for _, row in caas_df.iterrows():
        gene = str(row["Gene"]).strip()
        try:
            pos = int(row["Position"])
        except Exception:
            if DEBUG:
                sys.stderr.write(f"[DEBUG] Skip invalid Position '{row['Position']}' for gene '{gene}'\n")
            continue

        site_id = f"{gene}_{pos}"
        change_side = str(row["change_side"]).strip().lower()
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
        aln, _ = aln_cache.get(gene, (None, None))
        if aln is None:
            continue
        if pos < 0 or pos >= aln.get_alignment_length():
            sys.stderr.write(f"[WARN] Position {pos} out of bounds for gene {gene} (len={aln.get_alignment_length()})\n")
            continue

        top_letters, bottom_letters = top_bottom_letters(str(row["caas"]))
        if not top_letters or not bottom_letters:
            if DEBUG:
                sys.stderr.write(f"[DEBUG] Skip {site_id} due to empty TOP/BOTTOM amino-acid sets\n")
            continue

        def _emit_internal_species(tax_set: Set[str]) -> None:
            for rec in aln:
                seq_id = str(rec.id).strip()
                tax = sp2tax.get(seq_id) if sp2tax else seq_id
                if tax is None or tax not in tax_set:
                    continue
                aa = str(rec.seq[pos])
                if aa in top_letters and aa not in AMBIGSYMS:
                    odir = "top"
                elif aa in bottom_letters and aa not in AMBIGSYMS:
                    odir = "bottom"
                else:
                    if not include_other:
                        continue
                    odir = "other"
                species_name = tax2sp.get(seq_id, tax2sp.get(tax, seq_id))
                out_rows.append({
                    "Gene": gene,
                    "Position": pos,
                    "Substitution": str(row["caas"]),
                    "species": species_name,
                    "observed_AA": aa,
                    "observed_dir": odir,
                    "trait": trait_name,
                    "trait_value": float(s2v.get(species_name, np.nan)),
                    "extreme_side": "INTERNAL",
                    "q1": q1,
                    "q3": q3,
                })

        _emit_internal_species(internal_tax)

    cols = ["Gene", "Position", "Substitution", "species", "observed_AA", "observed_dir",
            "trait", "trait_value", "extreme_side", "q1", "q3"]
    if not out_rows:
        sys.stderr.write("[WARN] No rows generated for validation table.\n")
        return pd.DataFrame(columns=cols)
    return pd.DataFrame(out_rows, columns=cols)
