#!/usr/bin/env python3
# about.py — Static About-tab content. No user-editable state.
# PhyloPhere | gui/models/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Modeled as a dataclass (rather than inlined constants in the widget) purely for
testability and to keep gui/widgets/ free of hardcoded prose — but it carries no
user-editable state and is not part of ProjectConfig's persisted JSON.
"""

# ── Standard library ──────────────────────────────────────────────────────────
from dataclasses import dataclass


@dataclass(frozen=True, kw_only=True)
class AboutInfo:
    name: str = "PhyloPhere"
    description: str = (
        "Nextflow DSL2 pipeline for phylogenetic comparative genome-phenome analyses, "
        "centered on CAAStools-based CAAS/CAAP discovery."
    )
    version: str = "2.0"
    repo_url: str = "https://github.com/nozerorma/PhyloPhere"
    author: str = "Miguel Ramon Alonso"
    author_email: str = "miguel.ramon@upf.edu"
    institution: str = "Institute of Evolutionary Biology (IBE), UPF-CSIC"
    license_name: str = "GNU General Public License v3.0"
    logo_relpath: str = "res/logo.png"
    attributions: tuple[str, ...] = (
        "Builds on CAAStools (github.com/linudz/caastools)",
        "CAAP grouping: Chen, S. & Zou, Z. (2025), Molecular Ecology Resources 25(1), e70052",
        "RERconverge: Kowalczyk et al. (2019), Bioinformatics 35(22)",
        "HyPhy FADE directional selection: Kosakovsky Pond et al. (2021), Mol. Biol. Evol. 37(1)",
        "PrimateAI-3D variant pathogenicity: Gao et al. (2023), Science 380(6648)",
        "COSMIC somatic mutation census: Tate et al. (2019), Nucleic Acids Res. 47(D1)",
        "DOMINO active module identification: Levi et al. (2021), Bioinformatics 37(21)",
        "STRING DB protein interaction database: Szklarczyk et al. (2023), Nucleic Acids Res. 51(D1)",
        "Ortholog Characterizator (github.com/nozerorma/ortholog_characterizator.git)",
        "Phenotype/test-scenario contributions: Maria Sanchez Bermudez",
    )
