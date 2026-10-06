#!/usr/bin/env python3
# general.py — General tab config: core nextflow.config params and the remote validation target.
# PhyloPhere | gui/models/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
GeneralConfig: seed, reporting toggle, repository and plugin paths, and the SSH host
against which paths are validated and browsed.

Environment installation (environment/install_env.sh) is not modeled here. The GUI
runs inside the `phylophere` environment (see run_gui.sh), so an "install the
environment" control would be circular on a first setup, and install_env.sh only
creates (it fails on an existing environment), so it cannot serve as a reinstall
either. First-time setup is a terminal step: `./environment/install_env.sh`.

Imported by: gui/models/project.py, gui/widgets/tabs/general_tab.py
"""

# ── Standard library ──────────────────────────────────────────────────────────
from dataclasses import dataclass


@dataclass(kw_only=True)
class GeneralConfig:
    # Display-only name, independent of the project file's path (which may be a
    # generic template filename or an autosave path). Shown in the main window's
    # title bar (see MainWindow._update_title).
    project_name: str = ""

    # --- Core params (conf/common.config) ---
    seed: str = "1998"  # --seed
    reporting: bool = True  # --reporting (RUN_REPORTING in reference script)

    # There is no prune_data toggle: pruning (--prune_data gate and --prune_list of
    # conf/common.config) is set per phenotype row and turns on whenever a row of the
    # Runtime tab's phenotype table names a PRUNE file (see PhenotypeRow.prune).

    # --- Paths needed by every generated script ---
    repo_dir: str = ""  # REPO_DIR: path to the PhyloPhere checkout containing main.nf
    nextflow_plugins_dir: str = ""  # symlinked into each run's NXF_HOME (see run_single.sh.j2)

    # --- Remote validation/browsing target ---
    # Dataset paths often live on a remote HPC cluster reached over SSH while the
    # GUI runs on a laptop. Empty means the local filesystem. Format user@host,
    # with key-based authentication (gui/remote.py runs ssh in BatchMode, so a
    # password prompt is never answered).
    remote_host: str = ""
    # Directory where the remote browse dialog (gui/widgets/common/remote_browse_dialog.py)
    # starts when a PathField is empty (e.g. "/scratch/mramon" instead of "/").
    # Only the starting point: browsing reaches anything the SSH user can read.
    remote_root_dir: str = ""
