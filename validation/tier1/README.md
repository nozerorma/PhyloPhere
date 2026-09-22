# Tier 1 — known-positive recovery on two real datasets

**Question:** on compact, well-studied datasets with an independent
positive-selection signal and a curated list of externally validated
convergent sites, how much of what we expect to see does PhyloPhere recover,
and what else does it find?

The workflow for each fixture: build it → run `ortholog_characterizator`
(FUBAR/MEME) on its CDS alignment as an independent positive-selection check →
run PhyloPhere's own pipeline via its GUI template → compare both against the
curated truth set in `../truthsets/tier1/`.

## Layout

```
tier1/
  input/        one subdirectory per dataset: alignment, tree(s), trait table,
                 species-name/gene-tree auxiliary files, build.py + build_cds.py
                 + README (provenance). Built fixture content is gitignored —
                 rebuild via each dataset's build.py/build_cds.py; only
                 build.py/build_cds.py/README/*.spec.json are tracked, per
                 validation/.gitignore.
  output/        local GUI run outputs land here (one subdir per dataset).
                 Gitignored except for a .gitkeep.
  templates/     hard copies of the GUI project files below, for version
                 tracking. The GUI itself loads from gui/templates/ (its "Load
                 template" dialog opens there) — keep both copies in sync.
  references/    papers the input data and truth sites came from.
  run_ortholog_characterizator.sh   launcher for the ortholog_characterizator
                 comparison run (see below).
```

## Datasets

| dataset | `input/` | trait(s) | genes | template |
|---------|----------|----------|-------|----------|
| PEPC / C4 photosynthesis (Cyperaceae) | `input/pepc/` | `c4` (categorical) | PEPC (1) | `gui/templates/tier1_pepc.json` |

See each `input/<dataset>/README.md` for full provenance and truth-site detail.
Truth sites live in `../truthsets/tier1/pepc_c4.sites.tsv`.

## Independent check: ortholog_characterizator

`run_ortholog_characterizator.sh {pepc}` runs
`/home/miguel/IBE-UPF/PhD/ortholog_characterizator` on the fixture's
`input/<dataset>/align_cds/` (built by `build_cds.py`) against `tree.nwk`,
producing FUBAR (and, where triggered, MEME) site-level positive-selection
calls independent of PhyloPhere's own pipeline. By default it runs
translation + positive_selection only (quality, phylogeny, and variability are
off — quality because every fixture tip is already a single curated
accession, phylogeny because `psel_species_tree` is always the fixture's own
`tree.nwk` rather than a tree built in-pipeline), with `hyphy_methods=FUBAR`
and `filter_mode=none`.

## What each GUI template runs

From `gui/templates/tier1_pepc.json`:
CAAS (discovery + resample) → CT_DISAMBIGUATION (`ct_disambig_asr_mode=compute`)
→ FADE → Enrichment (FCS gene-set enrichment; STRING/DOMINO AMI off via
`scoring_ami`/`scoring_string`) → SCORING. **Accumulation and RER are off**
(`modules.accumulation.enabled` / `modules.rer.enabled` = false). **VEP is
off** — it needs per-gene alignment-to-protein MAP files (`vep_map_dir`),
which require a codon-level CDS-to-protein mapping neither fixture provides
to PhyloPhere's own pipeline. **POSENRICH is off**
(`posenrich_enabled=false`) — it needs an eggNOG members/annotations file and
a FUBAR sites file (`egg_members_file`, `egg_annotations_file`,
`fubar_sites_file`), none of which are wired in from the
ortholog_characterizator run above.

Every other file-based input is left blank so PhyloPhere generates or fetches
it itself: `tax_id_file`, `accumulation_entropy_dir`, `gmt_dir`,
`string_db_dir`/`string_cache_dir`, `domain_variability_file`,
`ucr_positions_file`. `gene_ensembl_file` is the one exception — both
templates point it at the fixture's own `input/<dataset>/gene_ensembl.tsv`.
`taxid.tsv` stays in `input/<dataset>/` for reference but isn't wired into
either template.

Both templates run `runtime_type: "local"` with `resources.local_max_cpus=8`,
`local_max_memory=16.GB`, `local_max_time=5.day`, plus the full
per-process `resources.process_overrides` list carried in the JSON.

## Running

1. Build a fixture: `python validation/tier1/input/<dataset>/scripts/build.py`
   and `.../scripts/build_cds.py`.
2. Independent positive-selection check:
   `bash validation/tier1/run_ortholog_characterizator.sh pepc`.
3. PhyloPhere's own pipeline: open the GUI → File → Load template →
   `tier1_pepc.json` → Generate Scripts → run.
   Or point the CLI runner directly at the same JSON.
4. Compare both runs' output against `../truthsets/tier1/*.sites.tsv`.
