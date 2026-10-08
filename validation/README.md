# Paper validation

The two cases that support the paper: one gene with a curated truth set (Tier 1, PEPC and C4 photosynthesis in sedges) and a genome-wide case in terrestrial mammals (Tier 2). They test whether PhyloPhere recovers published signal, not whether the code is correct. The tests of the code are in the sibling directory `PhyloPhere_validation/`, outside this repository.

## Layout

| Path | Content |
|---|---|
| `tier1/input/` | fixture builders and provenance: `pepc/` (genotypic C4 annotation), `.pepc_phenotypic/` (phenotypic annotation, same alignment and tree) |
| `tier1/scripts/` | run and comparison tools: `render_tier1_scripts.py`, `compare_pepc_runs.py`, `pepc_resources_tables.py` |
| `tier1/templates/` | copy of the GUI project of the case, `tier1_pepc_c4.json` |
| `tier1/reports/` | results, genotypic against phenotypic annotation, and resource use of the PEPC runs |
| `tier1/references/` | the papers the fixture and the truth sites come from |
| `tier1/output/` | one directory per run; ignored by git |
| `tier1/previous_work/` | runs and negative controls made with the pipeline before the unified CAAS core; see its README |
| `tier2/references/` | supplementary data of the terrestrial-mammal case; the case itself has no fixture or run yet |
| `truthsets/tier1/pepc_c4.sites.tsv` | the published C4-associated PEPC sites, in maize PEPC1 numbering |

## Tier 1 procedure

All commands start from the repository root.

1. **Build the fixture.** `python3 validation/tier1/input/pepc/scripts/build.py`, then `build_cds.py` for the codon alignment. The built files are ignored by git and rebuilt by these scripts.
2. **Render the run scripts.** `python3 validation/tier1/scripts/render_tier1_scripts.py --template gui/templates/tier1_pepc_c4.json --outdir validation/tier1/output/<run> --evidence-top-n 30` writes the scripts of the GUI template into a fresh directory, with new results, work and ASR-cache directories, so a rerun never reuses an earlier run.
3. **Run them.** Copy the two generated scripts to the repository root (`/run_*.sh` is git-ignored): the batch script looks for the single-trait script there, and the latter takes its `REPO_DIR` (`main.nf`, `bin/`) from its own location. Then run the batch script locally with `bash`, or on a cluster with `sbatch`, from `output/<run>/` so that `run.log`, `started.at` and `exit.code` stay with the run.
4. **Compare runs.** `python3 validation/tier1/scripts/compare_pepc_runs.py --a <results A> --b <results B> --require-equal` reports the positions, scores and empirical p-values of two runs of the same trait, the ten truth sites side by side and the cycles of the permulation null each run holds. It exits with 1 when the runs differ beyond `--tol`.
5. **Position tables and calibration.** `python3 validation/tier1/scripts/pepc_pvalue_tables.py --results <results> [--previous <dir>]` prints the truth-set table, the tail counts of `p.emp` and `p.emp_fact`, the null of each truth site and the multiple-testing families, and compares with the position tables of another run. `Rscript validation/tier1/scripts/pepc_null_calibration.R <results>/<trait>_complete <out_prefix>` takes each null cycle as observed data against the others, with the `.fact_*` functions of `scoring_compute.R`, and writes the share of (position, cycle) pairs with p <= alpha by propensity class and the share of cycles with at least one adjusted call.
6. **Resource report.** `pepc_resources_tables.py` prints the tables of `tier1/reports/pepc_resources.md` from the trace and log files a run leaves behind.
7. **Independent check.** `tier1/input/pepc/scripts/run_ortholog_characterizator.sh` runs `ortholog_characterizator` (FUBAR and MEME) on the codon alignment, which gives a positive-selection call set that does not depend on PhyloPhere.
8. **Negative controls.** The 20 controls of `tier1/previous_work/input/pepc_negctrl/` were run with the previous pipeline and have not been repeated with the current one; `previous_work/README.md` says how to repeat them.

## What a run covers

The template runs CAAS discovery and the permulation null, the ancestral-state reconstruction of the core (`ct_disambig_asr_mode=compute`), FADE, enrichment and scoring. Accumulation, RER, VEP and POSENRICH are off: the last two need inputs the fixture does not provide (protein maps and an eggNOG annotation). The reports in `tier1/reports/` give the settings and the numbers of the runs they describe.

## Version control

Tracked: builders, specs, READMEs, templates, reports, references and truth sets. Ignored (see `.gitignore`): the built fixture files and everything under `tier1/output/`. The GUI loads its templates from `gui/templates/`, so `tier1/templates/tier1_pepc_c4.json` and `gui/templates/tier1_pepc_c4.json` must stay identical.
