# Wiring tests

Tests of the Nextflow wiring that run Nextflow itself, so they are skipped when `nextflow` is not on PATH.

- `test_wiring.py`: the null's gene universe (DAG edge from the cleaned background, via `-preview`), the defaults of
  `filter_maxcaas` and `caas_map_dir`, the null's post-processing arguments, the exploratory grid and the selected
  cluster file, and the real `CT_FILTER` process with and without the MAP directory.
- `test_gui_map_dir.py`: the MAP directory as a post-processing parameter in the GUI model, the templates, the
  generators and the validation, and the agreement of the single-run script's fallbacks with the config defaults.
- `mini_*.nf`: small scripts that `test_wiring.py` runs as the `main.nf` of a temporary project whose folders link
  the tree under test, so `baseDir` resolves there. They import the real modules and call the real functions and
  process instead of repeating their logic.

`PHYLOPHERE_ROOT=<checkout>` runs the tests against another tree, which shows that a test fails on the code before a
change.

```bash
python3 -m pytest -q validation/wiring/
```
