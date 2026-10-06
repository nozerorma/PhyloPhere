#!/usr/bin/env bash
# install_env.sh — Create the phylophere conda environment and install the R packages that conda does not provide.
# PhyloPhere | environment/
#
# Called by:  the user, once, from the repository root (./environment/install_env.sh)
# Inputs:     $1  environment file (default: phylophere.yml in the current directory; the repository ships environment/phylophere.yml)
# Outputs:    the conda environment "phylophere": created from the file, or updated when it exists
#
# The solver is the first of micromamba, mamba and conda found on the PATH (micromamba is also
# searched at its usual install locations). The R packages installed afterwards are DT, the
# Bioconductor dependencies of RERconverge, and RERconverge itself, taken from a pinned commit
# and compiled with the C++17 standard (see the comments in the R block).

set -euo pipefail

# ── Settings ──────────────────────────────────────────────────────────────────

ENV_YML="${1:-phylophere.yml}"
ENV_NAME="phylophere"

# ── Solver selection ──────────────────────────────────────────────────────────

# Prints the path of a micromamba executable, or nothing.
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

# Prints the solver to use (micromamba path, mamba, conda) or "none".
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

# ── Conda environment ─────────────────────────────────────────────────────────

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

# ── R packages ────────────────────────────────────────────────────────────────

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
pkgdir <- tempfile("RERconverge_")
# git clone rather than download.file: codeload.github.com tarball fetches time
# out on some networks where the git protocol still works.
# system2() bypasses the shell, so pass path args raw (no quoting).
stopifnot(system2("git", c("clone", "--quiet",
  "https://github.com/nclark-lab/RERconverge.git", pkgdir)) == 0L)
stopifnot(system2("git", c("-C", pkgdir, "checkout", "--quiet", rer_ref)) == 0L)

# NOTE: no regex backslash escapes below (no "\\s", no "\\+"). This block is
# passed via `Rscript -e ${single-quoted}` and some callers (nested ssh/bash -c,
# GUI subprocess wrappers) strip one backslash level, turning "\\s" into a bare
# "\s" that R rejects as an unrecognized escape. fixed=TRUE keeps it literal.
mv <- file.path(pkgdir, "src", "Makevars")
for (f in c(mv, paste0(mv, ".win"))) {
  txt <- if (file.exists(f)) readLines(f) else character(0)
  txt <- txt[!grepl("CXX_STD", txt, fixed = TRUE)]
  writeLines(c(txt, "CXX_STD = CXX17"), f)
}
desc <- file.path(pkgdir, "DESCRIPTION")
writeLines(gsub("C++11", "C++17", readLines(desc), fixed = TRUE), desc)

remotes::install_local(pkgdir, dependencies = NA, upgrade = "never")

if (!requireNamespace("RERconverge", quietly = TRUE)) {
  stop("RERconverge failed to install")
}

cat("OK: R deps installed\n")
'

echo "Done."
