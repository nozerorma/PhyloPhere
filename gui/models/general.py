#!/usr/bin/env python3
# general.py — General tab config: core nextflow.config knobs + remote validation target.
# PhyloPhere | gui/models/
#
# Author: Miguel Ramon (miguel.ramon@upf.edu)

"""
Environment installation (environment/install_env.sh) is deliberately NOT modeled
here — the GUI itself runs inside the `phylophere` env (see run_gui.sh), so any
"install the environment" control would be circular for a first-time setup (can't
run the GUI to install the thing the GUI needs to run), and install_env.sh is
create-only (errors on an already-existing env) so it can't double as a safe
"reinstall" action either. First-time setup stays a one-time terminal step:
`./environment/install_env.sh`.
"""

# ── Standard library ──────────────────────────────────────────────────────────
from dataclasses import dataclass


@dataclass(kw_only=True)
class GeneralConfig:
    # --- Core nextflow.config / conf/common.config knobs ---
    seed: str = "1998"  # --seed
    reporting: bool = True  # --reporting (RUN_REPORTING in reference script)

    # --- Repo-relative paths needed by every generated script ---
    repo_dir: str = ""  # REPO_DIR — path to the PhyloPhere checkout containing main.nf
    nextflow_plugins_dir: str = ""  # symlinked into each run's NXF_HOME (see run_single.sh.j2)

    # --- Remote validation/browsing target ---
    # Dataset paths often live on a remote HPC cluster reached over SSH while the
    # GUI itself runs on a laptop. Empty = validate/browse the local filesystem.
    # user@host form, relying on the user's existing SSH key-based auth (no
    # password-prompt support — see gui/remote.py, which runs ssh in BatchMode).
    remote_host: str = ""
