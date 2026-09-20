# Tier 1 — site-level truth sets (GUI-driven local runs)

**Question:** on compact, well-studied datasets with a curated list of externally
validated convergent sites, does PhyloPhere recover them, and where do they rank?

Two real datasets, each with its own input/output split, a ready-to-load GUI
project template, a narrative report, and the literature/data provenance behind
it.

`../harness/`, `../truthsets/`, `../docs/DESIGN.md` are the separate formal
metrics toolkit (precision/recall/rank-of-known-positive, null calibration) for
turning a completed run into a scored benchmark. Nothing here depends on them.

## Layout

```
tier1/
  input/        one subdirectory per dataset: alignment, tree(s), trait table,
                 species-name/gene-tree auxiliary files, build.py + README
                 (provenance). Built fixture content is gitignored — rebuild via
                 each dataset's build.py; only build.py/README/*.spec.json are
                 tracked, per validation/.gitignore.
  output/        local GUI run outputs land here (work_dir/results_dir per the
                 templates below). Gitignored.
  templates/     hard copies of the GUI project files below, for version
                 tracking. The GUI itself loads from gui/templates/ (its "Load
                 template" dialog opens there) — keep both copies in sync.
  reports/       one concise report per dataset: data source, scope, procedure,
                 results/conclusions.
  references/    papers + repositories the input data and truth sites came from.
```

## Datasets

| dataset | `input/` | trait(s) | genes | template |
|---------|----------|----------|-------|----------|
| PEPC / C4 photosynthesis (Cyperaceae) | `input/pepc/` | `c4` (categorical) | PEPC (1) | `gui/templates/tier1_pepc.json` |
| Haemoglobin / high-altitude adaptation (Sino-Himalayan tits) | `input/hb/` | `elev_mid` (continuous) + `altitude` (categorical) | HBA, HBD, HBB (3) | `gui/templates/tier1_hb_altitude.json` |

Four genes total across both datasets — see each `input/<dataset>/README.md` for
full provenance and truth-site detail.

## What each template runs

CAAS (discovery+resample) → CT_DISAMBIGUATION (`asr_mode=compute`) →
RERconverge → FADE → Accumulation → Enrichment (FCS gene-set enrichment +
STRING/DOMINO AMI) → SCORING. **VEP is off** — it needs per-gene alignment-to-protein
MAP files (`vep_map_dir`), which require a codon-level CDS alignment neither
dataset has. **POSENRICH is off** — it needs an eggNOG members/annotations file
and a FUBAR sites file, neither generated in-house yet.

Every other file-based input is left blank so PhyloPhere generates or fetches it
itself: `tax_id_file`, `gene_ensembl_file`, `accumulation_entropy_dir`, `gmt_dir`,
`string_db_dir`/`string_cache_dir`, `domain_variability_file`,
`ucr_positions_file`. The fixtures' own synthetic `taxid.tsv`/`gene_ensembl.tsv`
stay in `input/<dataset>/` for reference but aren't wired into the templates.

Both templates run `runtime_type: "local"`, resources set to the
`local_lowspec` preset (`conf/resources.config.local_lowspec`: 8 cpu / 16 GB /
5 day ceiling, per-process overrides loaded via
`gui/resource_presets.load_preset("local_lowspec")`).

## Gene trees (RER)

RER's `--gene_trees` file (`input/pepc/gene_trees.nwk`, `input/hb/gene_trees.nwk`)
is IQ-TREE output from `ortholog_characterizator`'s `PHYLOGENY` workflow
(ModelFinder `MFP` restricted to the LG family, 1000 UFBoot), run directly
against each dataset's existing protein alignment (`--prot_dir`, quality/
translation stages off — neither dataset has a codon-level CDS alignment for
those stages to consume). Raw per-gene output kept in
`input/<dataset>/gene_trees_oc/` (gitignored).

## Running

Open Phylophere's GUI → File → Load template → `tier1_pepc.json` or
`tier1_hb_altitude.json` → Generate Scripts → run. Or point the CLI runner directly at
the same JSON.
