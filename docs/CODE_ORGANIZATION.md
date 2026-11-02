# Code organization and test inventory

Inventory reviewed locally on 2026-11-02.

## Pipeline and standalone tools

`bin/` contains active pipeline helpers, their imported dependencies, and input resolvers used by runner examples:

| Files | Consumer |
|---|---|
| `build_reference_tables.py`, `compute_variability.py`, `detect_ucr.py` | `subworkflows/POSENRICH/ucr_generation.nf` |
| `compute_domain_variability.py`, `map_domain_variability.py` | `subworkflows/POSENRICH/domain_variability_generation.nf` |
| `ppga_position_references.py` | Imported by `build_reference_tables.py` for its explicit PPGA conversion command; shares domain helpers |
| `resolve_core_inputs.py`, `generate_taxid_map.py`, `generate_ensembl_mapping.py` | Runner input preparation |
| `resolve_gmts.py` | Runner GMT resolution |

Standalone reference conversion, audits, scientific validation controls and older UCR batch drivers live in the ignored local `adhoc_scripts/` directory. They are not workflow stages. The one-off def2 PP2A recovery helper and launcher live locally in `adhoc_scripts/references/def2_pp2a/` and remain ignored.

`adhoc_scripts/` contains personal/local work. Its local README records the detailed cleanup, historical dependencies and locations of preserved data. No remote files or jobs were changed during this cleanup.

## Tests remaining in workflows and subworkflows

`workflows/` and `subworkflows/` contain no remaining software test files or nf-test specifications.

The two remaining unit-test files were removed at the user's request after rerunning them successfully (**10 passed**, 2026-11-02):

- `CT_DISAMBIGUATION/local/src/phylo/test_species_mapping.py`: 3 checks of taxonomy-ID matching, unmapped species and duplicate IDs.
- `TRAIT_ANALYSIS/local/src/test_tree_cleanup.py`: 7 checks of species curation, synonyms, taxonomy outputs and tree-only CLI execution.

These tests were regression checks, not runtime dependencies. Their removal leaves the pipeline implementation unchanged and removes those automated checks for future edits. No end-to-end Nextflow or cluster run was performed.

`subworkflows/POSENRICH/accumulation_test.nf` and `local/src/accumulation_test.py` implement a statistical analysis used by POSENRICH; their names do not designate software tests.

The three former `CT_DISAMBIGUATION/local/tests/` files (`test_annotation_ids.py`, `test_asr_cache_race.py`, `test_asr_existing_cache.py`) were already absent in the working tree before this cleanup. They were not restored or counted as remaining tests. The ignored local `validation/` tree is outside this inventory.

## Relocation checks

Python syntax, shell syntax and parsing of the relocated R sources passed. The six reference/UCR Python CLIs accept `--help` when launched outside the checkout. A synthetic alignment passed through all three relocated UCR drivers and produced a populated position table. WP9 resolves the repository root correctly at its new nesting depth.

The initial cleanup recorded 144 moved files and verified the original SHA-256 checksums of 99 data, image and rendered-report files at that time. Source files retain their original-content checksums in the local manifest, before path updates.

The local Python environment lacks `pyarrow`, so the position-audit CLI and the local def2 recovery helper were checked for syntax only. Scientific audits and Slurm jobs were not executed.

The former root `scripts/` directory was subsequently consolidated into `adhoc_scripts/` (references, UCR and validation subdirectories). All standalone utilities are local and ignored by Git; pipeline code does not call them. Paths and local launchers were updated accordingly.
