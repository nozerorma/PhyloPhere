#!/usr/bin/env python3
# validate.py — Pre-render validation: required fields per enabled module, CLASS rules.
# PhyloPhere | gui/generation/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Pure function, no PySide6 import. Returns a list of human-readable error strings;
an empty list means the project is renderable. Used both by tests/gui and by the
GUI's "Generate Scripts..." action (which should show these in a dialog rather than
letting render.py fail with a StrictUndefined KeyError deep inside a template).

Scope is deliberately the "essential fields" surface only (see implementation plan
§5) — this checks the fields actually exposed as GUI widgets, not the full
conf/*.config parameter space.
"""

# ── Standard library ──────────────────────────────────────────────────────────
import os

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.precomputed import derive_paths
from gui.models.project import ProjectConfig


def validate(project: ProjectConfig) -> list[str]:
    errors: list[str] = []

    def require(value: str, message: str) -> None:
        if not value.strip():
            errors.append(message)

    # --- General ---
    require(project.general.repo_dir, "General: PhyloPhere repo directory is required.")
    require(
        project.general.nextflow_plugins_dir,
        "General: Nextflow plugins directory is required (symlinked into every run's NXF_HOME).",
    )

    # --- Runtime ---
    rt = project.runtime
    pc = project.precomputed
    require(rt.alignment_dir, "Runtime: alignment directory is required.")
    require(rt.tree_file, "Runtime: species tree file is required.")
    require(rt.work_dir, "Runtime: work directory is required.")
    require(rt.results_dir, "Runtime: results directory is required.")

    if not rt.phenotype_rows:
        errors.append("Runtime: the phenotype catalogue must have at least one row.")

    has_class1 = any(row.trait_class == 1 for row in rt.phenotype_rows)
    has_class2 = any(row.trait_class == 2 for row in rt.phenotype_rows)
    if has_class1:
        require(rt.trait_file, "Runtime: trait file is required (a CLASS 1 phenotype row exists).")
        require(rt.prune_dir, "Runtime: prune directory is required (a CLASS 1 phenotype row exists).")
    if has_class2:
        require(
            rt.simple_trait_file,
            "Runtime: simple trait file is required (a CLASS 2 phenotype row exists).",
        )

    for i, row in enumerate(rt.phenotype_rows, start=1):
        require(row.trait, f"Phenotype row {i}: trait name is required.")
        if row.trait_class not in (1, 2):
            errors.append(f"Phenotype row {i}: CLASS must be 1 or 2, got {row.trait_class!r}.")
        if row.trait_class == 2:
            require(row.discrete_method, f"Phenotype row {i}: discrete method is required for CLASS 2.")
        if str(row.trait_type).strip().lower() not in ("", "auto", "ordinal", "continuous"):
            errors.append(
                f"Phenotype row {i}: TRAIT_TYPE must be blank, auto, ordinal, or continuous "
                f"(got {row.trait_type!r})."
            )

    # --- Precomputed Run tab: base_path required whenever any reuse box is checked ---
    any_precomp_checked = any(
        [
            pc.use_discovery, pc.use_resample, pc.use_bootstrap, pc.use_disambiguation,
            pc.use_postproc, pc.use_accumulation, pc.use_rer, pc.use_fade, pc.use_vep,
        ]
    )
    if any_precomp_checked:
        require(
            pc.base_path,
            "Precomputed Run: base path is required — at least one reuse checkbox is checked.",
        )

    # --- CAAS ---
    # caas_config_path is never required from this GUI. CONTRAST_SELECTION is the
    # trait-file supplier, and run_single.sh.j2 turns it on whenever CAAS **or**
    # Disambiguation is enabled — so a trait file always exists for whichever
    # consumer needs one:
    #   • CAAS on  -> CONTRAST_SELECTION's trait/tree feed CT(...) directly.
    #   • CAAS off, Disambiguation on (Precomputed Run reuse) -> ct_tool is empty so
    #     CT() never runs, but CONTRAST_SELECTION still does, and main.nf hands its
    #     trait_file_out/tree_file_out to CT_DISAMBIGUATION.
    # The --caas_config fallbacks (ct.nf, ct_disambiguation.nf) are therefore only
    # reachable by standalone non-GUI invocations, which is why the GUI can drive
    # phenotypes purely through --my_traits.
    #
    # CAAS's own output (discovery/resample/bootstrap) is only ever consumed
    # downstream by Disambiguation, via run_signification (main.nf:184-199) — so this
    # is only a real problem when Disambiguation is actually enabled. The more precise
    # version of this same check lives in the `disambig.enabled` branch below; a
    # phenotype/reporting-only run (main.nf:155's standalone --contrast_selection path,
    # or --reporting alone) has no downstream module reading CAAS output at all.
    caas = project.modules.caas
    disambig = project.modules.disambiguation
    if not caas.enabled and disambig.enabled:
        if not any([pc.use_discovery, pc.use_resample, pc.use_bootstrap]):
            errors.append(
                "CAAS is disabled but no Discovery/Resample/Bootstrap reuse box is checked on "
                "the Precomputed Run tab — downstream modules (Disambiguation, Accumulation) "
                "have no input."
            )

    # --- Disambiguation (+ Post-processing) ---
    if disambig.enabled:
        if disambig.ct_disambig_asr_mode == "precomputed":
            require(
                disambig.ct_disambig_asr_cache_dir,
                "Disambiguation: ASR cache directory is required when asr_mode=precomputed.",
            )
        if not caas.enabled and not any([pc.use_discovery, pc.use_resample, pc.use_bootstrap]):
            errors.append(
                "Disambiguation is enabled but CAAS is disabled with no reuse box checked on "
                "the Precomputed Run tab."
            )
        # Post-processing's characterization report step always runs alongside
        # Disambiguation now (no separate --ct_postproc toggle; see
        # gui/models/modules.py's DisambiguationConfig) and needs gene_ensembl_file
        # regardless of whether Scoring itself is enabled — checked here too (not
        # just under Scoring below) since exploratory-sweep runs force Scoring off
        # but still hit this requirement.
        if not project.modules.scoring.gene_ensembl_file:
            errors.append(
                "Disambiguation: Post-processing runs alongside it, so Scoring's gene-Ensembl "
                "mapping file is required for its characterization reports (and for gene "
                "filtering too, unless Gene filter mode is set to 'none') — fill it in on "
                "the Scoring tab even if Scoring itself is disabled."
            )
    # else: Disambiguation disabled means Post-processing is too (they're no longer
    # independently toggleable) — Accumulation/Scoring's own checks below already
    # cover the case where they still need Post-processing's output via pc.use_postproc.

    # --- Accumulation ---
    accum = project.modules.accumulation
    scoring = project.modules.scoring
    if accum.enabled:
        require(
            accum.accumulation_entropy_dir,
            "Accumulation: entropy directory is required when Accumulation is enabled.",
        )
        if not disambig.enabled and not pc.use_postproc:
            errors.append(
                "Accumulation is enabled but Post-processing is disabled with no "
                "'Use precomputed Post-processing output' box checked on the Precomputed Run tab."
            )
    elif scoring.enabled:
        # scoring_accum_ch (main.nf:411) is Accumulation's only downstream consumer —
        # Enrichment reads Scoring's gene lists, never Accumulation's directory
        # directly, so this only matters when Scoring is actually enabled.
        if not pc.use_accumulation:
            errors.append(
                "Accumulation is disabled but Scoring is enabled and 'Use precomputed Accumulation "
                "output' isn't checked on the Precomputed Run tab — Scoring may have no input."
            )

    # --- RERconverge ---
    rer = project.modules.rer
    if rer.enabled:
        require(rer.gene_trees, "RERconverge: gene trees file is required when RER is enabled.")

    # --- FADE ---
    # FADE parameters always have defaults; nothing strictly required.

    # --- VEP ---
    vep = project.modules.vep
    if vep.enabled:
        require(vep.vep_map_dir, "VEP: per-gene MAP directory is required when VEP is enabled.")

    # --- Scoring ---
    if scoring.enabled:
        if not disambig.enabled:
            # Already reported above (Disambiguation section) when Disambiguation
            # (and therefore Post-processing) is on, so this doesn't duplicate that
            # message for the same missing field.
            require(
                scoring.gene_ensembl_file,
                "Scoring: gene-Ensembl mapping file is required when Scoring is enabled.",
            )
        if not disambig.enabled and not pc.use_postproc:
            errors.append(
                "Scoring is enabled but Post-processing is disabled with no 'Use precomputed "
                "Post-processing output' box checked on the Precomputed Run tab."
            )
        if not rer.enabled and not pc.use_rer:
            errors.append(
                "Scoring is enabled and RERconverge is disabled — 'Use precomputed RERconverge "
                "output' isn't checked on the Precomputed Run tab."
            )
        if not project.modules.fade.enabled and not pc.use_fade:
            errors.append(
                "Scoring is enabled and FADE is disabled — 'Use precomputed FADE output' isn't "
                "checked on the Precomputed Run tab."
            )

    # --- Enrichment (+ POSENRICH) ---
    enrichment = project.modules.enrichment
    if enrichment.enabled:
        require(enrichment.gmt_dir, "Enrichment: GMT directory is required when Enrichment is enabled.")
    if enrichment.posenrich_enabled:
        # cosmic_db is NOT required: workflows/enrichment.nf and subworkflows/VEP/cosmic.nf
        # both fall back to a NO_FILE sentinel and skip the COSMIC overlap section
        # gracefully rather than failing when it's absent.
        for field_name, label in [
            ("egg_members_file", "eggNOG members file"),
            ("egg_annotations_file", "eggNOG annotations file"),
            ("domain_variability_file", "domain variability file"),
            ("ucr_positions_file", "UCR positions file"),
            ("fubar_sites_file", "FUBAR sites file"),
        ]:
            require(getattr(enrichment, field_name), f"POSENRICH: {label} is required.")

    return errors


def path_entries(project: ProjectConfig) -> list[tuple[str, str, str]]:
    """(label, path, kind) for every path-like field in the project. kind is
    "file" or "dir". Empty values are the caller's job to skip — validate() above
    already reports missing-but-required fields; this only checks paths that ARE
    filled in but don't exist on disk.

    Shared with gui/remote.py's SSH-based existence check, so local and remote
    validation always agree on exactly which fields count as paths.
    """
    general = project.general
    runtime = project.runtime
    m = project.modules
    pc = project.precomputed

    entries: list[tuple[str, str, str]] = [
        ("General: repo directory", general.repo_dir, "dir"),
        ("General: Nextflow plugins directory", general.nextflow_plugins_dir, "dir"),
        ("Runtime: work directory", runtime.work_dir, "dir"),
        ("Runtime: results directory", runtime.results_dir, "dir"),
        ("Runtime: alignment directory", runtime.alignment_dir, "dir"),
        ("Runtime: species tree", runtime.tree_file, "file"),
        ("Runtime: CLASS 1 trait file", runtime.trait_file, "file"),
        ("Runtime: CLASS 2 trait file", runtime.simple_trait_file, "file"),
        ("Runtime: prune directory", runtime.prune_dir, "dir"),
        ("Runtime: alignment species names", runtime.ali_sp_names, "file"),
        ("Runtime: taxonomy ID mapping", runtime.tax_id_file, "file"),
        ("CAAS: config file", m.caas.caas_config_path, "file"),
        ("CAAS: trait values file", m.caas.traitvalues, "file"),
        ("Disambiguation: ASR cache directory", m.disambiguation.ct_disambig_asr_cache_dir, "dir"),
        ("Accumulation: entropy directory", m.accumulation.accumulation_entropy_dir, "dir"),
        ("RERconverge: gene trees", m.rer.gene_trees, "file"),
        ("RERconverge: tested-gene universe file", m.rer.rer_universe_file, "file"),
        ("RERconverge: cross-module gene scores", m.rer.rer_gene_scores, "file"),
        ("FADE: gene-set top genes", m.fade.fade_postproc_top, "file"),
        ("FADE: gene-set bottom genes", m.fade.fade_postproc_bottom, "file"),
        ("FADE: LG substitution matrix", m.fade.lg_dat_path, "file"),
        ("FADE: tested-gene universe file", m.fade.fade_universe_file, "file"),
        ("VEP: PrimateAI-3D database", m.vep.vep_primateai_db, "file"),
        ("VEP: MAP directory", m.vep.vep_map_dir, "dir"),
        ("VEP: COSMIC database", m.vep.cosmic_db, "file"),
        ("Scoring: gene-Ensembl file", m.scoring.gene_ensembl_file, "file"),
        ("Enrichment: GMT directory", m.enrichment.gmt_dir, "dir"),
        ("Enrichment: STRING database directory", m.enrichment.string_db_dir, "dir"),
        ("Enrichment: eggNOG members file", m.enrichment.egg_members_file, "file"),
        ("Enrichment: eggNOG annotations file", m.enrichment.egg_annotations_file, "file"),
        ("Enrichment: COSMIC database", m.enrichment.cosmic_db, "file"),
        ("Enrichment: domain variability file", m.enrichment.domain_variability_file, "file"),
        ("Enrichment: UCR positions file", m.enrichment.ucr_positions_file, "file"),
        ("Enrichment: FUBAR sites file", m.enrichment.fubar_sites_file, "file"),
    ]

    # --- Precomputed Run tab: real per-phenotype derived paths, not stored strings ---
    for i, row in enumerate(runtime.phenotype_rows, start=1):
        if not row.trait:
            continue
        for name, path, kind in derive_paths(pc, row.trait):
            entries.append((f"Precomputed ({row.trait}): {name}", path, kind))

    for i, row in enumerate(runtime.phenotype_rows, start=1):
        if row.prune:
            entries.append((f"Phenotype row {i}: prune", os.path.join(runtime.prune_dir, row.prune), "file"))
        if row.prune_secondary:
            entries.append(
                (
                    f"Phenotype row {i}: prune_secondary",
                    os.path.join(runtime.prune_dir, row.prune_secondary),
                    "file",
                )
            )

    return entries


def validate_paths(project: ProjectConfig) -> list[str]:
    """Check that every filled-in path field actually exists on disk.

    Complements validate() (which only checks required-ness) — a field can be
    non-empty but point at a typo'd or not-yet-mounted path, which validate()
    can't catch since it never touches the filesystem.
    """
    problems: list[str] = []
    for label, path, kind in path_entries(project):
        if not path.strip():
            continue  # required-ness is validate()'s job, not this function's
        exists = os.path.isfile(path) if kind == "file" else os.path.isdir(path)
        if not exists:
            noun = "file" if kind == "file" else "directory"
            problems.append(f"{label}: {noun} not found — {path}")
    return problems
