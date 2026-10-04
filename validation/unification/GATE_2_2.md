# Gate 2.2: reviewing the toy run of the b_0 core

Stage 2.2 makes the observed labeling the `b_0` slice of the permulation core. The gate checks that a run of the new code
reproduces the results of the former observed chain on the same toy design.

## What is compared

Run A is the baseline (`..._BASE_2.2_sorted_N1000`, made with the former observed chain on the sorted alignment listing:
1000 genes, N = 1000, `ct_core_batch_size` 20). Run B is the same design and settings run with `scoring_v2` at the merge of
Stage 2.2 (42d853b) or later.

| Tool | Compares | Passes when |
|---|---|---|
| `compare_b0.py --run B` | the observed results with the `b_0` slice of the same run, checkpoints A-E | A, B, C, D, E pass (|delta| <= 1e-12) |
| `compare_null.py --a A --b B --extra scoring/position_scores.tsv --extra scoring/gene_scores.tsv` | the null tables, `caas_perms.rds` and the scoring tables | every table equal (the null does not depend on the observed path: bit for bit) |
| `compare_contract.py --a A --b B` | the observed contract files: `caastools/discovery.tab`, the background files, `meta_caas/meta_caas/*`, `ct_disambiguation/caas_convergence_master.csv` | every file passes (rules in the tool's docstring) |

Differences that are expected and are not failures: the CAAS ids (content hashes now, recomputed and checked by
`compare_contract.py`), the order of the rows inside a position, and the modal-residue cells where two residues tie for the
maximum support (counted in the report as `modal_residue_tie_cells_tolerated`). Numbers of the baseline for the same genes:
checkpoint A 63 765 rows, B 9 727 rows with 25 removed units, C 9 944 rows (delta 0.0), D 3 283 rows, E 548 genes.

## Steps

1. On the cluster checkout, bring `scoring_v2` to the merge of Stage 2.2 (your own git flow) and rename the results
   directory of the baseline run before launching anything that writes to the same name.
2. Launch the toy run with the same template and settings as the baseline (toy mode, seed 1998, 1000 genes, N = 1000,
   `ct_core_batch_size` 20). Until the toy-cycles parameter below exists, the sbatch template forces `CAAS_FULL_PERMS=100`
   in toy mode: set it to 1000 (and keep `MAX_TRIES` consistent) in the generated script.
3. When it has finished, review it on a compute node:

   ```bash
   srun -p high-cpu -c 2 --mem=8G -t 20 bash validation/unification/review_gate_2_2.sh <baseline dir> <new dir>
   ```

   The script writes `compare_b0.json` and `compare_contract.json` plus one log per tool into `<new dir>/gate_2_2` and exits
   1 if any tool fails. `--dry-run` prints the three commands.
4. On a failure, read the file named in the failing line: `only_a` / `only_b` are rows present in one run only (with
   examples), `max_abs_delta` and `worst_column` locate a numeric difference, `modal_residue_cells_differing_without_a_tie`
   lists the cells that are not explained by a tie.
5. Record from the run (not part of the pass criteria): wall time, per-process peak memory and CPU
   (`sacct` TotalCPU/Elapsed and MaxRSS of the batch step, with the nodes: the same task can be several times slower on one
   node than on another), in particular `CAAS_CORE_BATCHED`, `CAAS_CORE_OBSERVED` and `CAAS_CORE_MERGE`.

## Next steps (Stage 2.3)

- **Toy-mode cycles.** A parameter to choose the number of permulation cycles of a toy run: a field next to the toy sample
  size on the Runtime tab (`gui/models/runtime.py`, `gui/widgets/tabs/runtime_tab.py`), default 100, written by
  `sbatch_array.sh.j2` in place of the fixed `CAAS_FULL_PERMS="100"`, with `MAX_TRIES` scaled with it (the two must stay
  coupled), a test, and the template and help text updated.
- **Script window.** In `MainWindow._save_generated_scripts` (`gui/widgets/main_window.py`) the "Generated scripts preview"
  window (`_show_preview`, kept in `self._preview_window`) stays open after the "Scripts saved" message is accepted: close it
  then. Test it with the Qt stub of the GUI tests if PySide6 is not installed.
- `N = 0` through the whole chain: `reaggregate_perm_scores.py` raises "no *.tsv.gz shards" when the null has no shard (no
  hit in any permuted cycle, or N = 0), and `ct_resample.nf` passes `caas_full_perms` to the resample generator unchecked.
  Then SCORING and ENRICHMENT without a null, and a discovery-only run without permuted labelings.
- Compute-mode ASR inside the batch task (time, memory, `codeml` concurrency) and the memory of `CAAS_CORE_OBSERVED` at full
  scale (the whole discovery table is held in memory once).
- Prune the legacy branches of `analyze_gene_disambiguation` (metadata read from a file, `trait_file_path`,
  `diagnostics_dir`, `asr_mode`) and the unreachable `src/asr/asr_only.py` and `node_identification.py`.
- The same name-sorted listing and seeded shuffle for the FADE selection (`subworkflows/SELECTION/selection_prep.nf`),
  after checking that its extension filter and the CT one select the same files.
- Exploratory grid on the `b_0` master; PEPC genotypic and phenotypic reruns and the negative controls (does the p.emp tail
  excess persist?); the audit of the GUI against `conf/*.config`; the document of the current null
  (`docs/CAAS_PERMULATION_EXCESS.md` is cited but does not exist).
- Delete the merged local branches (`phase2_core_wiring`, `phase21_*`, `phase22_observed_b0`) when no longer wanted.
