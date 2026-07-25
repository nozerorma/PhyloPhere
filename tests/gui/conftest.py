#!/usr/bin/env python3
# conftest.py — Shared pytest fixtures for the GUI's generation-layer tests.
# PhyloPhere | tests/gui/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
`golden_project()` builds a ProjectConfig whose values mirror the real, hand-written
SBATCH_run_phenotypes_primates.sh / run_phenotype_single_primates.sh (same paths,
same toggles, same 2-row phenotype catalogue) — see implementation plan §8, Phase A.

The `qapp` fixture forces the offscreen Qt platform plugin before any PySide6 import
touches the platform, so widget tests run headlessly without needing QT_QPA_PLATFORM
set externally (CI, this sandbox, etc.).
"""

# ── Standard library ──────────────────────────────────────────────────────────
import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

# ── Third-party ───────────────────────────────────────────────────────────────
import pytest

# ── Local ─────────────────────────────────────────────────────────────────────
from gui.models.project import ProjectConfig
from gui.models.runtime import PhenotypeRow


@pytest.fixture(scope="session")
def qapp():
    from PySide6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])
    yield app


@pytest.fixture(autouse=True)
def _reset_remote_context():
    """gui.widgets.common.remote_context is a module-level global (see its own
    docstring for why) — reset it around every test so no test can leak a remote
    host setting into an unrelated one."""
    from gui.widgets.common import remote_context

    remote_context.set_remote_host("")
    yield
    remote_context.set_remote_host("")


def golden_project() -> ProjectConfig:
    p = ProjectConfig()

    p.general.repo_dir = "/homes/users/mramon/scratch/0.Phylophere"
    p.general.nextflow_plugins_dir = "/homes/users/mramon/.nextflow/plugins"
    p.general.seed = "1998"
    p.general.reporting = True

    p.runtime.resume = True
    p.runtime.toy_mode = False
    p.runtime.runtime_type = "slurm"
    p.runtime.batched = True
    p.runtime.work_dir = "/homes/users/mramon/scratch/3.Work_dirs2"
    p.runtime.results_dir = (
        "/homes/users/mramon/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/final"
    )
    p.runtime.use_tower = True
    p.runtime.sbatch_mail_user = "miguel.ramon@upf.edu"
    p.runtime.sbatch_array_concurrency = "2"
    p.runtime.alignment_dir = (
        "/homes/users/mramon/scratch/4.Generate_alignments_from_codons/"
        "results_ppga/20260615_200027/PROT/bmge"
    )
    p.runtime.tree_file = (
        "/homes/users/mramon/scratch/2.Primates/1.Primates_data/5.Phylogeny/"
        "science.abn7829_data_s4.nex.tree"
    )
    p.runtime.trait_file = (
        "/homes/users/mramon/scratch/2.Primates/1.Primates_data/1.Cancer_data/"
        "Neoplasia_species360/cancer_traits_processed-LQ.csv"
    )
    p.runtime.simple_trait_file = (
        "/homes/users/mramon/scratch/2.Primates/1.Primates_data/maria_caas/"
        "Datos_fenotipos/diet_traitfile_comma.csv"
    )
    p.runtime.prune_dir = (
        "/homes/users/mramon/scratch/2.Primates/1.Primates_data/1.Cancer_data/"
        "Neoplasia_species360/ZAK-CLEANUP"
    )
    p.runtime.branch_trait = "LQ"
    p.runtime.ali_sp_names = (
        "/homes/users/mramon/scratch/2.Primates/1.Primates_data/2.Alignments/ali_sp_names.txt"
    )
    p.runtime.tax_id_file = (
        "/homes/users/mramon/scratch/2.Primates/1.Primates_data/5.Phylogeny/"
        "taxid_species_family_primates.tsv"
    )
    p.runtime.phenotype_rows = [
        PhenotypeRow(
            trait_class=1,
            trait="neoplasia_prevalence",
            secondary="malignant_prevalence",
            n_trait="adult_necropsy_count",
            c_trait="neoplasia_necropsy",
            prune="neoplasia_exclude.txt",
            prune_secondary="malignant_exclude.txt",
            discrete_method="quintile",
            scoring_rer_input=(
                "/homes/users/mramon/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/"
                "cancer/neoplasia_prevalence/runtime/filter/rerconverge/rer_results/"
                "rerconverge_summary_neoplasia_prevalence.tsv"
            ),
            scoring_fade_summary_top=(
                "/homes/users/mramon/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/"
                "cancer/neoplasia_prevalence/runtime/filter/selection/fade/top/fade_summary_top.tsv"
            ),
            scoring_fade_summary_bottom=(
                "/homes/users/mramon/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/"
                "cancer/neoplasia_prevalence/runtime/filter/selection/fade/bottom/"
                "fade_summary_bottom.tsv"
            ),
        ),
        PhenotypeRow(
            trait_class=1,
            trait="malignant_prevalence",
            secondary="neoplasia_prevalence",
            n_trait="adult_necropsy_count",
            c_trait="malignant_count",
            prune="malignant_exclude.txt",
            prune_secondary="neoplasia_exclude.txt",
            discrete_method="quintile",
            scoring_rer_input=(
                "/homes/users/mramon/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/"
                "cancer/malignant_prevalence/runtime/filter/rerconverge/rer_results/"
                "rerconverge_summary_malignant_prevalence.tsv"
            ),
            scoring_fade_summary_top=(
                "/homes/users/mramon/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/"
                "cancer/malignant_prevalence/runtime/filter/selection/fade/top/fade_summary_top.tsv"
            ),
            scoring_fade_summary_bottom=(
                "/homes/users/mramon/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/"
                "cancer/malignant_prevalence/runtime/filter/selection/fade/bottom/"
                "fade_summary_bottom.tsv"
            ),
        ),
    ]

    caas = p.modules.caas
    caas.caas_config_path = p.runtime.trait_file
    caas.perms_cycles = "1000000"
    caas.caas_full_perms = "1000"
    caas.caas_permulation_enrichment = True

    p.modules.disambiguation.ct_disambig_asr_cache_dir = (
        "/homes/users/mramon/scratch/2.Primates/1.Primates_data/7.ASR_primates"
    )

    p.modules.accumulation.accumulation_n_randomizations = "1000000"
    p.modules.accumulation.accumulation_entropy_dir = (
        "/homes/users/mramon/scratch/4.Generate_alignments_from_codons/"
        "results_ppga/20260615_200027/PROT_VAR_RAW"
    )

    p.modules.rer.enabled = False  # RUN_RER=false in the reference script
    p.modules.fade.enabled = False  # RUN_FADE=false in the reference script

    p.modules.vep.vep_primateai_db = (
        "/homes/users/mramon/scratch/2.Primates/1.Primates_data/4.External_DBs/"
        "PAI3D/PrimateAI-3D.hg38.txt.gz"
    )
    p.modules.vep.vep_map_dir = (
        "/homes/users/mramon/scratch/4.Generate_alignments_from_codons/"
        "results_ppga/20260615_200027/MAP"
    )
    p.modules.scoring.gene_ensembl_file = (
        "/homes/users/mramon/scratch/2.Primates/1.Primates_data/2.Alignments/ensembl_genes.output"
    )
    enrichment = p.modules.enrichment
    enrichment.gmt_dir = (
        "/homes/users/mramon/scratch/2.Primates/1.Primates_data/4.External_DBs/GMTs"
    )
    enrichment.string_db_dir = (
        "/homes/users/mramon/scratch/2.Primates/1.Primates_data/4.External_DBs/STRING"
    )
    enrichment.egg_members_file = (
        "/homes/users/mramon/scratch/2.Primates/1.Primates_data/4.External_DBs/eggNOG/"
        "9443_members.tsv"
    )
    enrichment.egg_annotations_file = (
        "/homes/users/mramon/scratch/2.Primates/1.Primates_data/4.External_DBs/eggNOG/"
        "9443_annotations.tsv"
    )
    enrichment.cosmic_db = (
        "/homes/users/mramon/scratch/2.Primates/1.Primates_data/4.External_DBs/COSMIC/"
        "Cosmic_MutantCensus_v104_GRCh38.tsv.gz"
    )
    enrichment.domain_variability_file = (
        "/homes/users/mramon/scratch/4.Generate_alignments_from_codons/"
        "results_ppga/20260615_200027/reports/variability/domain_variability.tsv"
    )
    enrichment.ucr_positions_file = (
        "/homes/users/mramon/scratch/4.Generate_alignments_from_codons/"
        "results_ppga/20260615_200027/reports/variability/ucr_positions.tsv"
    )
    enrichment.fubar_sites_file = (
        "/homes/users/mramon/scratch/4.Generate_alignments_from_codons/"
        "results_ppga/20260615_200027/reports/positive_selection/fubar_sites.tsv"
    )

    return p
