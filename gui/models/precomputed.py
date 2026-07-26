#!/usr/bin/env python3
# precomputed.py — Precomputed-input fallbacks, consolidated from every module tab.
# PhyloPhere | gui/models/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Every module tab used to carry its own "used when disabled" fallback fields
(discovery_from on CAAS, signification_from on Disambiguation, vep_caas_input on
VEP, ...). In practice a user configuring a standalone/resume run has to fill in
several of these across several tabs to skip straight to a later module, and two of
them (RER's rer_continuous_file/rer_perms_file, FADE's fade_json_dir_top/bottom)
were previously grayed-out-looking dead ends: those tabs had an empty fallback
group because the *per-phenotype* scoring_rer_input / scoring_fade_summary_top/
bottom fields live on PhenotypeRow instead, and nothing filled the *global*
precomputed-report-input slot documented in conf/rerconverge.config and
conf/fade.config. This dataclass gathers every global (non-per-phenotype)
precomputed/fallback path into one place so they're all visible together instead
of hidden one-per-tab.

Per-phenotype fallbacks (scoring_rer_input, scoring_fade_summary_top/bottom) stay
on PhenotypeRow (gui/models/runtime.py) — they legitimately vary per phenotype in
a batch run, so centralizing them here would be a regression, not a consolidation.
"""

# ── Standard library ──────────────────────────────────────────────────────────
from dataclasses import dataclass


@dataclass(kw_only=True)
class PrecomputedConfig:
    """Global standalone/resume-run inputs, used when the producing module is off."""

    # --- CAAS / CT (used when CAAS is disabled) ---
    discovery_from: str = ""  # --discovery_from
    resample_from: str = ""  # --resample_from
    bootstrap_from: str = ""  # --bootstrap_from

    # --- Disambiguation (used when Disambiguation is disabled) ---
    signification_from: str = ""  # --signification_from
    disambiguation_input: str = ""  # --disambiguation_input
    disambiguation_dir: str = ""  # --disambiguation_dir
    background_input: str = ""  # --background_input

    # --- Accumulation (used when Accumulation is disabled) ---
    accumulation_caas_input: str = ""  # --accumulation_caas_input
    accumulation_background_input: str = ""  # --accumulation_background_input

    # --- RERconverge (used when RER is disabled) ---
    # rer_continuous_file: precomputed RERconverge continuous-analysis RDS, renders
    # 5.RERconverge_report.html directly without running RER live (mirrors
    # fade_json_dir_top/bottom's role for FADE below).
    rer_continuous_file: str = ""  # --rer_continuous_file
    rer_perms_file: str = ""  # --rer_perms_file

    # --- FADE (used when FADE is disabled) ---
    # Directories of prior *.FADE.json results (one per gene); when set, FADE's own
    # report renders directly from them even with --fade off this invocation.
    fade_json_dir_top: str = ""  # --fade_json_dir_top
    fade_json_dir_bottom: str = ""  # --fade_json_dir_bottom

    # --- VEP (used when VEP is disabled) ---
    vep_caas_input: str = ""  # --vep_caas_input

    # --- Scoring (used when Scoring is disabled, or when a producing module is off) ---
    scoring_postproc_input: str = ""  # --scoring_postproc_input
    scoring_accum_dir: str = ""  # --scoring_accum_dir
    scoring_vep_primateai: str = ""  # --scoring_vep_primateai
    scoring_vep_cosmic: str = ""  # --scoring_vep_cosmic (COSMIC *scores* TSV fallback — distinct from VEP's own cosmic_db)
    scoring_background_input: str = ""  # --scoring_background_input
    caas_perms_file: str = ""  # --caas_perms_file
    scoring_fade_site_top: str = ""  # --scoring_fade_site_top (FADE site-level fallback, if FADE off)
    scoring_fade_site_bottom: str = ""  # --scoring_fade_site_bottom

    # --- Enrichment (used when Enrichment/Scoring is disabled) ---
    accumulation_enrichment_gene_lists_input: str = ""  # --accumulation_enrichment_gene_lists_input
    posenrich_background_file: str = ""  # --posenrich_background_file
