#!/usr/bin/env bash
set -euo pipefail

ENV_YML="${1:-phylophere.yml}"
ENV_NAME="phylophere"

find_micromamba() {
  if command -v micromamba >/dev/null 2>&1; then
    command -v micromamba
  elif [[ -n "${MAMBA_EXE:-}" && -x "$MAMBA_EXE" ]]; then
    echo "$MAMBA_EXE"
  elif [[ -x "$HOME/.local/bin/micromamba" ]]; then
    echo "$HOME/.local/bin/micromamba"
  elif [[ -x "$HOME/.micromamba/bin/micromamba" ]]; then
    echo "$HOME/.micromamba/bin/micromamba"
  else
    echo ""
  fi
}

choose_solver() {
  local mm
  mm="$(find_micromamba)"
  if [[ -n "$mm" ]]; then
    echo "$mm"
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
  echo "Please install some conda variant. Micromamba is recommended for simplicity."
  echo "https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html"
  exit 1
fi

echo "Using solver: $SOLVER"
echo "Creating env: $ENV_NAME from $ENV_YML"

if [[ "$SOLVER" == *"micromamba"* ]]; then
  : "${MAMBA_ROOT_PREFIX:=$HOME/.micromamba}"
  export MAMBA_ROOT_PREFIX
  "$SOLVER" config set channel_priority flexible >/dev/null

  echo "Resolving dependencies with micromamba (log level: info)..."
  "$SOLVER" env create -n "$ENV_NAME" -f "$ENV_YML" -y --log-level info || \
  "$SOLVER" env update -n "$ENV_NAME" -f "$ENV_YML" --log-level info

  RUN=("$SOLVER" run -n "$ENV_NAME")
elif [[ "$SOLVER" == "mamba" ]]; then
  mamba config --set channel_priority flexible >/dev/null
  echo "Resolving dependencies with mamba..."
  mamba env create -n "$ENV_NAME" -f "$ENV_YML" -y || \
  mamba env update -n "$ENV_NAME" -f "$ENV_YML"
  RUN=(mamba run -n "$ENV_NAME")
elif [[ "$SOLVER" == "conda" ]]; then
  conda config --set channel_priority flexible >/dev/null
  echo "Resolving dependencies with conda..."
  conda env create -n "$ENV_NAME" -f "$ENV_YML" -y || \
  conda env update -n "$ENV_NAME" -f "$ENV_YML"
  RUN=(conda run -n "$ENV_NAME")
fi

echo "Installing R packages (CRAN + Bioconductor + GitHub) into: $ENV_NAME"

"${RUN[@]}" Rscript -e '
# Fix Conda R bug where SHLIB_LIBADD is undefined, causing data.table source build failure
Sys.setenv(SHLIB_LIBADD = "")

options(
  repos = c(CRAN="https://cloud.r-project.org"),
  Ncpus = max(1L, parallel::detectCores() - 1L),
  buildtools.check = function(action) TRUE
)

# Make sure remotes exists before GitHub installs
install.packages("remotes")

# Install CRAN package without pulling/compiling deps
install.packages("DT", repos="https://cloud.r-project.org")

# Install dependencies for RERconverge
install.packages("BiocManager")
BiocManager::install("ggtree", ask = FALSE, update = FALSE)
BiocManager::install("impute", ask = FALSE, update = FALSE)
BiocManager::install("castor", ask = FALSE, update = FALSE)
BiocManager::install("data.table", ask = FALSE, update = FALSE)

# Prevent pkgbuild from throwing missing-toolchain warnings
options(buildtools.check = function(action) TRUE)

# RERconverge pins CXX_STD = CXX11, but modern RcppArmadillo (>= 12, and the
# 15.x build in this env) requires at least C++14 and fails compiler_check.hpp.
# Fetch the pinned commit, bump the C++ standard, and install from the local dir.
rer_ref <- "2bd328f7530b4aca9b48c0b3997875c9b77a7026"
tb <- tempfile(fileext = ".tar.gz")
utils::download.file(
  sprintf("https://github.com/nclark-lab/RERconverge/archive/%s.tar.gz", rer_ref),
  tb, quiet = TRUE
)
ex <- tempfile(); dir.create(ex)
utils::untar(tb, exdir = ex)
pkgdir <- list.files(ex, full.names = TRUE)[1]

mv <- file.path(pkgdir, "src", "Makevars")
for (f in c(mv, paste0(mv, ".win"))) {
  txt <- if (file.exists(f)) readLines(f) else character(0)
  txt <- grep("^\\s*CXX_STD\\s*=", txt, value = TRUE, invert = TRUE)
  writeLines(c(txt, "CXX_STD = CXX17"), f)
}
desc <- file.path(pkgdir, "DESCRIPTION")
writeLines(gsub("C\\+\\+11", "C++17", readLines(desc)), desc)

remotes::install_local(pkgdir, dependencies = NA, upgrade = "never")

if (!requireNamespace("RERconverge", quietly = TRUE)) {
  stop("RERconverge failed to install")
}

cat("OK: R deps installed\n")
'

echo "Done."
