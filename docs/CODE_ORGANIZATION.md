# Code organization

## Pipeline helpers

`bin/` contains the pipeline helpers, their imported dependencies and the input resolvers used by the runner scripts:

| Files | Consumer |
|---|---|
| `build_reference_tables.py`, `compute_variability.py`, `detect_ucr.py` | `subworkflows/POSENRICH/ucr_generation.nf` |
| `compute_domain_variability.py`, `map_domain_variability.py` | `subworkflows/POSENRICH/domain_variability_generation.nf` |
| `ppga_position_references.py` | Imported by `build_reference_tables.py` for its explicit PPGA conversion command; shares domain helpers |
| `resolve_core_inputs.py`, `generate_taxid_map.py`, `generate_ensembl_mapping.py` | Runner input preparation |
| `resolve_gmts.py` | Runner GMT resolution |

Standalone reference conversion, audits and one-off drivers live in the ignored local `adhoc_scripts/` directory. They are not workflow stages and the pipeline does not call them.

## Tests

The repository contains no unit-test or nf-test files. `validation/` holds the paper validation fixtures: tracked scripts, specifications and truth sets, and untracked regenerated data. `subworkflows/POSENRICH/accumulation_test.nf` and `local/src/accumulation_test.py` implement a statistical test used by POSENRICH, not a software test.
