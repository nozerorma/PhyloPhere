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
    micromamba env create -n "$ENV_NAME" -f "$ENV_YML" -y || micromamba env update -n "$ENV_NAME" -f "$ENV_YML" -y
    micromamba install -n "$ENV_NAME" -c conda-forge r-data.table compilers make pkg-config -y
    RUN=(micromamba run -n "$ENV_NAME")
    ;;
  mamba)
    mamba config --set channel_priority flexible >/dev/null
    mamba env create -n "$ENV_NAME" -f "$ENV_YML" -y || mamba env update -n "$ENV_NAME" -f "$ENV_YML" -y
    mamba install -n "$ENV_NAME" -c conda-forge r-data.table compilers make pkg-config -y
    RUN=(mamba run -n "$ENV_NAME")
    ;;
  conda)
    conda config --set channel_priority flexible >/dev/null
    conda env create -n "$ENV_NAME" -f "$ENV_YML" -y || conda env update -n "$ENV_NAME" -f "$ENV_YML" -y
    conda install -n "$ENV_NAME" -c conda-forge r-data.table compilers make pkg-config -y
    RUN=(conda run -n "$ENV_NAME")
    ;;
esac

echo "Installing R packages (CRAN + GitHub) into: $ENV_NAME"

"${RUN[@]}" Rscript -e '
options(
  repos = c(CRAN="https://cloud.r-project.org"),
  Ncpus = max(1L, parallel::detectCores() - 1L),
  buildtools.check = function(action) TRUE
)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("DT", quietly = TRUE)) install.packages("DT")

BiocManager::install(c("ggtree", "impute", "data.table", "castor"), update = FALSE, ask = FALSE)

remotes::install_github(
  "nclark-lab/RERconverge@2bd328f7530b4aca9b48c0b3997875c9b77a7026",
  dependencies = NA,
  upgrade = "never"
)

if (!requireNamespace("RERconverge", quietly = TRUE)) {
  stop("RERconverge failed to install")
}

cat("OK: R deps installed\n")
'

echo "Done."
