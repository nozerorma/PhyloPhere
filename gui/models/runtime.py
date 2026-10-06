#!/usr/bin/env python3
# runtime.py — Runtime tab config: execution target, phenotype catalogue, dataset paths, Tower.
# PhyloPhere | gui/models/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Runtime: -resume and toy-mode toggles, local or slurm execution, the phenotype
catalogue (one row per task of the SBATCH array, rendered as the `case
$SLURM_ARRAY_TASK_ID` block of sbatch_array.sh.j2), and the dataset paths shared by
every phenotype of a batch (the "COMMON THINGS" block of that template).

The Seqera/Tower access token is deliberately not a field: it must never reach the
human-readable JSON project file. The widget keeps it as transient state and passes it
to gui/secrets_io.py, which writes it to the gitignored token.tk that the tower{} block
of conf/common.config reads.

Imported by: gui/models/project.py, gui/widgets/phenotype_table/model.py,
gui/widgets/tabs/runtime_tab.py
"""

# ── Standard library ──────────────────────────────────────────────────────────
from dataclasses import dataclass, field


@dataclass(kw_only=True)
class PhenotypeRow:
    """One arm of the generated `case $SLURM_ARRAY_TASK_ID in ... esac` block."""

    trait: str = ""  # TRAIT / --traitname
    secondary: str = ""  # SECONDARY / --secondary_trait, optional
    n_trait: str = ""  # NTRAIT / --n_trait, optional (sample size; paired with c_trait)
    c_trait: str = ""  # CTRAIT / --c_trait, optional (observed cases; paired with n_trait)
    prune: str = ""  # PRUNE filename, joined with runtime.prune_dir, optional
    prune_secondary: str = ""  # PRUNE_SEC filename, optional
    trait_type: str = ""  # TRAIT_TYPE / --trait_type, optional: "" = auto-infer,
    # "ordinal" = coded fg/bg (highest level FG, lowest BG), "continuous" = force the
    # Phylogenetic Shift Score (PSS) pair-selection path. The fg/bg partition itself is
    # produced by 4.Independent_contrasts.Rmd, not by a quantile cut.

    # The precomputed inputs of Scoring (scoring_rer_input, scoring_rer_perms_input,
    # scoring_fade_summary_top/bottom) are not row fields: gui/models/precomputed.py
    # derives them per phenotype as base_path/<trait>/...


@dataclass(kw_only=True)
class RuntimeConfig:
    # --- Top toggles ---
    resume: bool = True  # -resume
    toy_mode: bool = False  # --toy_mode
    toy_n: str = "1000"  # --toy_n
    toy_perms: str = "100"  # permulation cycles of a toy run (CAAS_FULL_PERMS; the generated script sets MAX_TRIES from it)

    # --- Execution target ---
    runtime_type: str = "slurm"  # "local" | "slurm"

    # Replaces the "phenotypes"/"phenotype" token in the names of the generated scripts
    # (run_phenotypes_local.sh, SBATCH_run_phenotypes.sh, run_phenotype_single.sh,
    # ..._exploratory.sh/..._complete.sh). Blank keeps the default names.
    script_base_name: str = ""

    # --- Directories ---
    work_dir: str = ""  # WORK_BASE: base for per-trait Nextflow work dirs
    results_dir: str = ""  # CAAS_OUTBASE: base for per-trait results

    # --- Seqera Cloud / Tower ---
    use_tower: bool = True  # -with-tower

    # --- SBATCH array-job wrapper (slurm only) ---
    sbatch_job_name: str = "phylophere"
    sbatch_partition: str = ""  # SLURM partition; "" = the cluster default
    sbatch_time: str = "144:00:00"
    sbatch_mail_user: str = ""
    sbatch_array_concurrency: str = ""  # the "%C" in --array=1-N%C; "" = uncapped

    # --- Dataset paths shared across every phenotype in the batch ---
    alignment_dir: str = ""  # --alignment (ALI_DIR)
    ali_format: str = "fasta"  # --ali_format
    tree_file: str = ""  # --tree (TREE_FILE)
    trait_file: str = ""  # --my_traits (TRAIT_FILE)
    prune_dir: str = ""  # PRUNE_DIR, joined with each row's prune/prune_secondary filename
    branch_trait: str = "LQ"  # --branch_trait
    ali_sp_names: str = ""  # --ali_sp_names (optional)
    tax_id_file: str = ""  # --tax_id (INPUT_TAX_ID)

    # --- Reporting / contrast-selection dataset shape (conf/common.config) ---
    sp_colname: str = "species"  # --sp_colname
    clade_name: str = "primates"  # --clade_name
    clade_taxid: str = "9443"  # Clade-level NCBI TaxID (e.g. for eggNOG orthogroups)
    taxon_of_interest: str = "family"  # --taxon_of_interest
    ref_species_name: str = "Homo_sapiens"  # Reference species name (in alignments / annotations)
    ref_species_taxid: str = "9606"  # Reference species NCBI TaxID (e.g. for STRING / BioMart)

    # --- Phenotype catalogue ---
    phenotype_rows: list[PhenotypeRow] = field(default_factory=list)
