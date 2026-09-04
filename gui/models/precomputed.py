#!/usr/bin/env python3
# precomputed.py — Precomputed-input reuse, consolidated from every module tab.
# PhyloPhere | gui/models/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
One base_path + one checkbox per producing stage, instead of one manually-typed
global path per field. A single global path field can't work once a batch run has
more than one phenotype (TRAIT varies per row) — the same problem the old per-row
Scoring fallback overrides on PhenotypeRow (scoring_rer_input, scoring_rer_perms_input,
scoring_fade_summary_top/bottom) existed to solve, one dialog at a time. Every file
this GUI generates writes to a known, stable layout under an outdir
(see run_single.sh.j2/sbatch_array.sh.j2's RESULTS_BASE), so the actual per-phenotype
path is always base_path/<TRAIT>/<known-subpath> — computed in the generated shell
script itself (using $TRAIT), not typed in here.

Checking a box both supplies the precomputed input AND turns off the module that
would otherwise recompute it live (see gui/generation/templates/run_single.sh.j2's
PRECOMP_* wiring) — recomputing something you're also feeding a precomputed answer
for is never what you want.

Path templates (relative to base_path/<TRAIT>), verified against each process's
actual publishDir rather than guessed:
  CT/CAAS      : caastools/{discovery,resample,bootstrap}.tab, caastools/background_genes.output,
                 caastools/background.output, signification/meta_caas/global_meta_caas.tsv,
                 caas_permulation/caas_perms.rds,
                 caas_permulation/perm_pos_{pval,sample,quantiles}.tsv,
                 caas_permulation/perm_pos_detail/  (one gz shard per gene; triggers CAAS_PERMS_REBUILD:
                   SCORING re-derives the null from this rather than trusting the cached
                   caas_perms.rds, which is only valid while it holds the same gene-level
                   statistic as the observed score)
  Disambiguation: ct_disambiguation/caas_convergence_master.csv, ct_disambiguation/ (dir)
  Post-processing: postproc/gene_filtering/filtered_discovery.tsv,
                 postproc/cleaned_backgrounds/cleaned_background_main.txt
  Accumulation : accumulation/ (dir; workflows/scoring.nf recurses into
                 <dir>/{top,bottom,all}/randomization/*.csv)
  RER          : rerconverge/rer_results/<TRAIT>.continuous.{output,perms.rds},
                 rerconverge/rer_results/rerconverge_summary_<TRAIT>.tsv (glob-resolved
                 at generation time in the shell script, filename carries a variable suffix)
  FADE         : selection/fade/{top,bottom}/json (dir),
                 selection/fade/{top,bottom}/fade_summary_{top,bottom}.tsv,
                 selection/fade/{top,bottom}/fade_site_bf_{top,bottom}.tsv
  VEP          : vep/primateai_mapped.tsv, vep/cosmic_scores.tsv

Some of these are conditionally published (gene_filter_mode != 'none' for the
Post-processing pair; optional emits for RER perms, FADE summary/site, VEP scores;
whole-subworkflow opt-in for caas_perms.rds) — checking the box still wires the path,
but the file only exists if the source run actually had the matching settings on.
"""

# ── Standard library ──────────────────────────────────────────────────────────
import glob
import os
from dataclasses import dataclass


@dataclass(kw_only=True)
class PrecomputedConfig:
    """Base path + per-stage reuse toggles. See module docstring for path layout."""

    base_path: str = ""  # per-phenotype dir = base_path/<TRAIT> (no toy/postproc-mode tag)

    # --- CT / CAAS (general + 3 specific, mirroring the CAAS tab's own
    # discovery/resample/bootstrap checkboxes) ---
    use_ct: bool = False  # turns off CAAS; also wires background_input, signification_from,
    # caas_perms_file, posenrich_background_file — all derive from CT's own concat/
    # permulation outputs, not from any of the 3 specific steps individually.
    use_discovery: bool = False  # wires discovery_from
    use_resample: bool = False  # wires resample_from
    use_bootstrap: bool = False  # wires bootstrap_from

    # --- Disambiguation (live ASR/convergence compute) ---
    use_disambiguation: bool = False  # turns off Disambiguation; wires disambiguation_input/_dir

    # --- Post-processing (no separate enable toggle — always runs alongside
    # Disambiguation; see gui/generation/context.py's ct_postproc_enabled) ---
    use_postproc: bool = False  # wires the filtered-discovery/cleaned-background
    # pair that Accumulation/VEP/Scoring each fall back to, without turning off
    # Disambiguation's own live compute (check use_disambiguation for that)

    # --- Accumulation ---
    use_accumulation: bool = False  # turns off Accumulation; wires scoring_accum_dir

    # --- RERconverge ---
    use_rer: bool = False  # turns off RER; wires rer_continuous_file/rer_perms_file plus
    # the per-phenotype scoring_rer_input/scoring_rer_perms_input Scoring falls back to

    # --- FADE ---
    use_fade: bool = False  # turns off FADE; wires fade_json_dir_top/bottom plus the
    # per-phenotype scoring_fade_summary_top/bottom and scoring_fade_site_top/bottom

    # --- VEP ---
    use_vep: bool = False  # turns off VEP; wires scoring_vep_primateai/scoring_vep_cosmic


def derive_paths(config: "PrecomputedConfig", trait: str) -> list[tuple[str, str, str]]:
    """(label, path, kind) for every path implied by config's checked boxes, for one
    phenotype's TRAIT. Python-side mirror of run_single.sh.j2's own PRECOMP_OUTDIR
    construction — used by gui/generation/validate.py for real existence checks (the
    shell script builds the identical paths at generation/run time; keep both in sync
    if the layout ever changes, see this module's docstring for the source of truth).
    """
    if not config.base_path or not trait:
        return []
    outdir = os.path.join(config.base_path, trait)
    entries: list[tuple[str, str, str]] = []

    if config.use_discovery:
        entries.append(("discovery_from", os.path.join(outdir, "caastools", "discovery.tab"), "file"))
        sig = os.path.join(outdir, "signification", "meta_caas", "global_meta_caas.tsv")
        if not os.path.isfile(sig):
            sig = os.path.join(outdir, "signification", "meta_caas", "meta_caas.tsv")
        entries.append(("signification_from", sig, "file"))
        entries.append(("background_input", os.path.join(outdir, "caastools", "background_genes.output"), "file"))
        entries.append(("posenrich_background_file", os.path.join(outdir, "caastools", "background.output"), "file"))
        # The permulation outputs travel together: caas_perms.rds alone feeds the
        # FCS p.perm but leaves the report's position-level permulation section empty,
        # which is gated on the per-position files. Mirrors run_single.sh.j2's
        # PRECOMP_USE_DISCOVERY block -- keep both in sync.
        perm_dir = os.path.join(outdir, "caas_permulation")
        entries.append(("caas_perms_file", os.path.join(perm_dir, "caas_perms.rds"), "file"))
        entries.append(("caas_pos_pval_file", os.path.join(perm_dir, "perm_pos_pval.tsv"), "file"))
        entries.append(("caas_pos_sample_file", os.path.join(perm_dir, "perm_pos_sample.tsv"), "file"))
        entries.append(("caas_pos_quantiles_file", os.path.join(perm_dir, "perm_pos_quantiles.tsv"), "file"))
        # caas_pos_detail_file makes SCORING REBUILD the null (CAAS_PERMS_REBUILD)
        # instead of importing caas_perms.rds as a cached artifact. That matters
        # because a cached null is only valid while it holds the same gene-level
        # statistic as the observed score; rebuilding guarantees it by construction,
        # costs minutes (no ASR replay), and is what keeps p.perm populated --
        # fcs_enrich.R leaves it NA when it detects a stale null. Wired last so it
        # takes precedence over caas_perms_file above.
        entries.append(("caas_pos_detail_file",
                        os.path.join(perm_dir, "perm_pos_detail"), "dir"))
    if config.use_resample:
        entries.append(("resample_from", os.path.join(outdir, "caastools", "resample.tab"), "file"))
    if config.use_bootstrap:
        entries.append(("bootstrap_from", os.path.join(outdir, "caastools", "bootstrap.tab"), "file"))

    if config.use_disambiguation:
        entries.append(
            ("disambiguation_input", os.path.join(outdir, "ct_disambiguation", "caas_convergence_master.csv"), "file")
        )
        entries.append(("disambiguation_dir", os.path.join(outdir, "ct_disambiguation"), "dir"))

    if config.use_postproc:
        gene_list = os.path.join(outdir, "postproc", "gene_filtering", "filtered_discovery.tsv")
        background = os.path.join(outdir, "postproc", "cleaned_backgrounds", "cleaned_background_main.txt")
        entries.append(("accumulation_caas_input / vep_caas_input / scoring_postproc_input", gene_list, "file"))
        entries.append(("accumulation_background_input / scoring_background_input", background, "file"))

    if config.use_accumulation:
        entries.append(("scoring_accum_dir", os.path.join(outdir, "accumulation"), "dir"))

    if config.use_rer:
        rer_dir = os.path.join(outdir, "rerconverge", "rer_results")
        entries.append(("rer_continuous_file", os.path.join(rer_dir, f"{trait}.continuous.output"), "file"))
        entries.append(("rer_perms_file / scoring_rer_perms_input", os.path.join(rer_dir, f"{trait}.continuous.perms.rds"), "file"))
        matches = sorted(glob.glob(os.path.join(rer_dir, "rerconverge_summary_*.tsv")))
        if matches:
            entries.append(("scoring_rer_input", matches[0], "file"))

    if config.use_fade:
        for direction in ("top", "bottom"):
            fade_dir = os.path.join(outdir, "selection", "fade", direction)
            entries.append((f"fade_json_dir_{direction}", os.path.join(fade_dir, "json"), "dir"))
            entries.append((f"scoring_fade_summary_{direction}", os.path.join(fade_dir, f"fade_summary_{direction}.tsv"), "file"))
            entries.append((f"scoring_fade_site_{direction}", os.path.join(fade_dir, f"fade_site_bf_{direction}.tsv"), "file"))

    if config.use_vep:
        vep_dir = os.path.join(outdir, "vep")
        entries.append(("scoring_vep_primateai", os.path.join(vep_dir, "primateai_mapped.tsv"), "file"))
        entries.append(("scoring_vep_cosmic", os.path.join(vep_dir, "cosmic_scores.tsv"), "file"))

    return entries
