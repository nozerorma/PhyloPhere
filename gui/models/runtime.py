#!/usr/bin/env python3
# runtime.py — Runtime tab config: local/slurm execution, phenotype catalogue, Tower.
# PhyloPhere | gui/models/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Covers -resume/toy-mode toggles, local-vs-slurm execution, the batched-SBATCH-array
phenotype catalogue (mirrors SBATCH_run_phenotypes_primates.sh's CASE block), and
the dataset paths shared across every phenotype in a batch ("COMMON THINGS" in the
reference script's own comment).

NOTE: the Seqera/Tower access token is intentionally NOT a field here. It must never
be persisted into the human-diffable JSON project file. The GUI widget collects it as
transient state and hands it straight to gui/secrets_io.py, which writes it to the
repo's existing gitignored token.tk (see conf/common.config's tower{} block, which
already reads that file automatically) — see implementation plan §4.
"""

# ── Standard library ──────────────────────────────────────────────────────────
from dataclasses import dataclass, field


@dataclass(kw_only=True)
class PhenotypeRow:
    """One arm of the generated `case $SLURM_ARRAY_TASK_ID in ... esac` block."""

    trait_class: int = 1  # CLASS 1 (full trait config) or CLASS 2 (simple discrete trait)
    trait: str = ""  # TRAIT / --traitname
    secondary: str = ""  # SECONDARY / --secondary_trait (CLASS 1 only)
    n_trait: str = ""  # NTRAIT / --n_trait (CLASS 1 only)
    c_trait: str = ""  # CTRAIT / --c_trait (CLASS 1 only)
    prune: str = ""  # PRUNE filename, joined with runtime.prune_dir (CLASS 1 only)
    prune_secondary: str = ""  # PRUNE_SEC filename (CLASS 1 only)
    discrete_method: str = ""  # DISCRETE / --discrete_method (CLASS 2 only)

    # Optional per-row precomputed-input overrides for Scoring (only emitted if non-empty).
    scoring_rer_input: str = ""
    scoring_rer_perms_input: str = ""  # --rer_perms_file, per-phenotype (distinct from
    # PrecomputedConfig.rer_perms_file, the global fallback used when RER itself is off)
    scoring_fade_summary_top: str = ""
    scoring_fade_summary_bottom: str = ""


@dataclass(kw_only=True)
class RuntimeConfig:
    # --- Top toggles ---
    resume: bool = True  # -resume
    toy_mode: bool = False  # --toy_mode
    toy_n: str = "1000"  # --toy_n

    # --- Execution target ---
    runtime_type: str = "slurm"  # "local" | "slurm"
    batched: bool = True  # slurm-only: generate the SBATCH array-job wrapper

    # --- Directories ---
    work_dir: str = ""  # WORK_BASE — base for per-trait Nextflow work dirs
    results_dir: str = ""  # CAAS_OUTBASE — base for per-trait results

    # --- Seqera Cloud / Tower ---
    use_tower: bool = True  # -with-tower

    # --- SBATCH array-job wrapper (slurm + batched only) ---
    sbatch_job_name: str = "phylophere"
    sbatch_partition: str = "haswell"
    sbatch_time: str = "144:00:00"
    sbatch_mail_user: str = ""
    sbatch_array_concurrency: str = ""  # the "%C" in --array=1-N%C; "" = uncapped

    # --- Dataset paths shared across every phenotype in the batch ---
    alignment_dir: str = ""  # --alignment (ALI_DIR)
    ali_format: str = "fasta"  # --ali_format
    tree_file: str = ""  # --tree (TREE_FILE)
    trait_file: str = ""  # CLASS 1's --my_traits (TRAIT_FILE)
    simple_trait_file: str = ""  # CLASS 2's --my_traits (SIMPLE_TRAIT_FILE)
    prune_dir: str = ""  # PRUNE_DIR, joined with each row's prune/prune_secondary filename
    branch_trait: str = "LQ"  # --branch_trait
    ali_sp_names: str = ""  # --ali_sp_names (optional)
    tax_id_file: str = ""  # --tax_id (INPUT_TAX_ID)

    # --- Reporting / contrast-selection dataset shape (conf/common.config) ---
    sp_colname: str = "species"  # --sp_colname
    clade_name: str = "primates"  # --clade_name
    taxon_of_interest: str = "family"  # --taxon_of_interest

    # --- Phenotype catalogue ---
    phenotype_rows: list[PhenotypeRow] = field(default_factory=list)
