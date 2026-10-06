# Paper validation

The two cases that support the paper: one gene with a curated truth set (Tier 1, PEPC and C4 photosynthesis in sedges) and a genome-wide case in terrestrial mammals (Tier 2). They test whether PhyloPhere recovers published signal, not whether the code is correct. The tests of the code are in the sibling directory `PhyloPhere_validation/`, outside this repository.

## Layout

| Path | Content |
|---|---|
| `tier1/input/` | fixture builders and provenance: `pepc/` (genotypic C4 annotation), `.pepc_phenotypic/` (phenotypic annotation, same alignment and tree), `pepc_negctrl/` (negative controls) |
| `tier1/scripts/` | run and comparison tools: `render_tier1_scripts.py`, `compare_pepc_runs.py`, `pepc_resources_tables.py` |
| `tier1/templates/` | copy of the GUI project of the case, `tier1_pepc_c4.json` |
| `tier1/reports/` | results, genotypic against phenotypic annotation, and resource use of the PEPC runs |
| `tier1/references/` | the papers the fixture and the truth sites come from |
| `tier1/output/` | one directory per run; ignored by git |
| `tier2/references/` | supplementary data of the terrestrial-mammal case; the case itself has no fixture or run yet |
| `truthsets/tier1/pepc_c4.sites.tsv` | the published C4-associated PEPC sites, in maize PEPC1 numbering |

## Tier 1 procedure

All commands start from the repository root.

1. **Build the fixture.** `python3 validation/tier1/input/pepc/scripts/build.py`, then `build_cds.py` for the codon alignment. The built files are ignored by git and rebuilt by these scripts.
2. **Render the run scripts.** `python3 validation/tier1/scripts/render_tier1_scripts.py --template gui/templates/tier1_pepc_c4.json --outdir validation/tier1/output/<run> --evidence-top-n 30` writes the scripts of the GUI template into a fresh directory, with new results, work and ASR-cache directories, so a rerun never reuses an earlier run.
3. **Run them.** Locally with `bash`, or on a cluster with `sbatch`.
4. **Compare runs.** `python3 validation/tier1/scripts/compare_pepc_runs.py --a <results A> --b <results B> --require-equal` reports the positions, scores and empirical p-values of two runs of the same trait, the ten truth sites side by side and the cycles of the permulation null each run holds. It exits with 1 when the runs differ beyond `--tol`.
5. **Resource report.** `pepc_resources_tables.py` prints the tables of `tier1/reports/pepc_resources.md` from the trace and log files a run leaves behind.
6. **Independent check.** `tier1/input/pepc/scripts/run_ortholog_characterizator.sh` runs `ortholog_characterizator` (FUBAR and MEME) on the codon alignment, which gives a positive-selection call set that does not depend on PhyloPhere.
7. **Negative controls.** `tier1/input/pepc/scripts/build_negctrl_traits.R` writes traits unrelated to the phenotype; `tier1/input/pepc_negctrl/run_negctrl_local.sh` runs them and `analyze_negctrl.py` compares the fraction of positions with `p.emp <= alpha` against alpha.

## What a run covers

The template runs CAAS discovery and the permulation null, the ancestral-state reconstruction of the core (`ct_disambig_asr_mode=compute`), FADE, enrichment and scoring. Accumulation, RER, VEP and POSENRICH are off: the last two need inputs the fixture does not provide (protein maps and an eggNOG annotation). The reports in `tier1/reports/` give the settings and the numbers of the runs they describe.

## Version control

Tracked: builders, specs, READMEs, templates, reports, references and truth sets. Ignored (see `.gitignore`): the built fixture files and everything under `tier1/output/`. The GUI loads its templates from `gui/templates/`, so `tier1/templates/tier1_pepc_c4.json` and `gui/templates/tier1_pepc_c4.json` must stay identical.
