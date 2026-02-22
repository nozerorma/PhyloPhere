#!/usr/bin/env bash
set -euo pipefail

ENV_YML="${1:-phylophere.yml}"
ENV_NAME="phylophere"

choose_solver() {
  if command -v micromamba >/dev/null 2>&1; then
    echo "micromamba"
  elif command -v mamba >/dev/null 2>&1; then
    echo "mamba"
  elif command -v conda >/dev/null 2>&1; then
    echo "conda"
  else
    echo "none"
  fi
}

SOLVER="$(choose_solver)"
if [[ "$SOLVER" == "none" ]]; then
  echo "ERROR: Need micromamba, mamba, or conda on PATH." >&2
  exit 1
fi

echo "Using solver: $SOLVER"
echo "Creating env: $ENV_NAME from $ENV_YML"

case "$SOLVER" in
  micromamba)
    : "${MAMBA_ROOT_PREFIX:=$HOME/.micromamba}"
    export MAMBA_ROOT_PREFIX
    micromamba config set channel_priority flexible >/dev/null
    # Force env name explicitly
    micromamba env create -n "$ENV_NAME" -f "$ENV_YML" -y
    RUN=(micromamba run -n "$ENV_NAME")
    ;;
  mamba)
    mamba config --set channel_priority flexible >/dev/null
    # Force env name explicitly
    mamba env create -n "$ENV_NAME" -f "$ENV_YML" -y
    RUN=(mamba run -n "$ENV_NAME")
    ;;
  conda)
    conda config --set channel_priority flexible >/dev/null
    # Force env name explicitly
    conda env create -n "$ENV_NAME" -f "$ENV_YML" -y
    RUN=(conda run -n "$ENV_NAME")
    ;;
esac

echo "Installing R packages (CRAN + GitHub) into: $ENV_NAME"

"${RUN[@]}" Rscript -e '
options(
  repos = c(CRAN="https://cloud.r-project.org"),
  Ncpus = max(1L, parallel::detectCores() - 1L)
)

# Make sure remotes exists before GitHub installs
if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes")

# Install CRAN package without pulling/compiling deps
if (!requireNamespace("rphylopic", quietly=TRUE)) {
  install.packages("rphylopic", dependencies = FALSE)
}

# Install CRAN package without pulling/compiling deps
if (!requireNamespace("rbioapi", quietly=TRUE)) {
  install.packages("rbioapi", dependencies = TRUE)
}

# Install dependencies for RERconverge
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ggtree")
BiocManager::install("impute")

# Install GitHub package only if missing
if (!requireNamespace("RERconverge", quietly=TRUE)) {
  remotes::install_github("nclark-lab/RERconverge@251bae1", dependencies = TRUE, upgrade="never")
}

cat("OK: R deps installed\n")
'

echo "Done."
