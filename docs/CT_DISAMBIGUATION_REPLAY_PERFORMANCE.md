# CT_DISAMBIGUATION replay performance — tiered plan

## Context

`CAAS_PERMS_DISAMBIGUATE_BATCHED` (the permulation-null replay) and `CT_DISAMBIGUATION_RUN_BATCHED`
(the observed disambiguation) are both burning far less than 1 core's worth of aggregate CPU
(~9-10% utilization on 8 allocated CPUs, per real `sacct` data from a 2026-09-16 cluster run) despite
real, non-trivial wall time (6-21 min per 20-gene batch, ~750+ batches, extrapolating to 14-30h for a
full permulation-null run).

This investigation started from the hypothesis that `gene_wrapper.py`'s `process_single_gene`/
`_perms_worker` — which walk each gene's ASR tree per-node over dict-keyed posteriors in plain Python
— needed BLAS-style vectorization, the way `subworkflows/CT/local/modules/perm_replay_vec.py` already
vectorized the CAAS discovery kernel into a GEMM.

**That hypothesis turned out to be wrong, or at least premature**, once real profiling replaced code-read
inference. The actual investigation history (documented below so a future conversation doesn't have to
re-derive it) went through three rounds of *"is this really the bottleneck?"*, each round overturning
part of the previous one. The tiers below are structured explicitly so this doesn't happen a fourth
time silently: **every tier that builds on a prior finding starts by re-checking that finding still
holds**, because everything here so far rests on a single gene (ABCB10), a single cluster (`correfoc`),
and a single day's snapshot of a live run.

### What is measured (high confidence) vs. inferred (lower confidence) — read this before trusting anything below

| Claim | Status | Evidence |
|---|---|---|
| `encoded_distribution`/`encode_aa` is ~72% of CPU time in `compute_domain_scores` | **Synthetic profile only** — my own fake tree, fake posteriors, guessed K=8 domains, 3000 fake cycles. Not real data. | cProfile run in an earlier session, no cluster access |
| The real workload replays **48,232 cycles across 20 genes** in one batch | **Measured, real** | `.command.log` on `correfoc`, job dir `.../neoplasia_prevalence_complete/20/a330bfd6392f731decb02c8021e9ca/` |
| Real tree: 229 tips, 228 internal nodes; alignment: 564 sites; posterior threshold 0.1 | **Measured, real** (one gene, ABCB10) | `rst`/`alignment_paml.phy` in `asr_ABCB10` cache dir |
| A cold-cache single-gene replay of 300 cycles: **real=36.5s, user+sys≈2.0s CPU** (~5.5% utilization) | **Measured, real**, `srun` on `correfoc` `std-cpu` | See "Tier 0 reproduction" below |
| The same 300 cycles warm-cache (re-run immediately after): **real=8.0s, user+sys≈1.5s** | **Measured, real**, same job, second invocation | — |
| `stat()` is called ~78×/cycle (23,479 calls / 300 cycles), `openat()` ~12.5×/cycle with ~19% failing | **Measured, real**, via `strace -c` | Same warm-cache run |
| The cold/warm gap is NFS latency on the per-cycle trait-file round-trip, and this is the dominant real-world cost | **CONFIRMED under real 8-way production parallelism, 2026-09-16 (Tier 0)** — see Tier 0 findings below | 6 independently-sampled live production jobs via `sstat` |
| LCA-per-domain-pair is fixed per-gene, cacheable across all cycles | **Wrong.** Domain→MRCA assignment is per-cycle (trait file encodes `pair_id = k` = position in that cycle's resampled fg/bg list) | `disambiguate_single.py:793-844`, `gene_wrapper.py:933-943` |
| `encoded_distribution` is a fixed 0/1-ish linear map, 5 schemes, cycle-independent | **Confirmed by code read** | `path_scores.py:67-105`, `grouping.py:5-136` — see correction below |
| GS2 scheme has 6 groups | **Wrong** — GS2 has **7** groups (GS1=6, GS2=7, GS3=6, GS4=12), US=20 singleton | `grouping.py:59` |
| `process_single_gene` (observed path) has the same repeated-reparse problem as `_perms_worker` | **Wrong.** It calls `parse_trait_pairs` once per gene, not once per cycle — only `_perms_worker`/the null-replay path has the N-cycles-×-reopen problem | `gene_wrapper.py:208-225` vs. `1123-1210` |
| No existing test exercises `analyze_gene_disambiguation`, `parse_trait_pairs`, or `_perms_worker` end-to-end | **Confirmed by code read** | see Tier 1 "regression net" below |

### Reproduction recipe (for re-running any of the above)

Real work dir used throughout: `/homes/users/mramon/scratch/3.Work_dirs_cancer_final_noprune_multi.fix148/neoplasia_prevalence_complete/20/a330bfd6392f731decb02c8021e9ca` on `correfoc`.
Staging area used for single-gene isolation: `~/scratch/0.Phylophere/.tmp/profile_disambig_ABCB10/` (still present on cluster; contains `run.sh`, `cycles.txt`, `perm_disc_single/`, `timing.log`, `strace_summary.log` from this session — reusable, or regenerate per the steps below).

To reproduce a single-gene, single-worker, real-data timing run:
1. `mkdir -p ~/scratch/0.Phylophere/.tmp/<name>` and symlink one gene's `perm_disc/<GENE>.*.perm_replay.discovery.output` into `<name>/perm_disc_single/` (isolates one gene; `--perm-discovery` has no per-gene CLI filter otherwise).
2. Extract a cycle-tag subset (`awk`/`sort -u`/`head` over the `cycle` column of that discovery file) into `<name>/cycles.txt`, comma-joined.
3. Run `disambiguation_perms_main.py` (same args as the real `.command.sh` in the work dir, `--workers 1`, `--cycles "$(cat cycles.txt)"`, `--output-dir <name>/out`) via `srun --partition=high-cpu --time=00:15:00 --cpus-per-task=1 --mem=6G --export=ALL`, on `correfoc`, through a login-shell wrapper (`bash -lc '...'`) since `srun`/`micromamba` are only on PATH in a login shell.
4. Wrap with bash's builtin `time` (not `/usr/bin/time` — not installed on the compute image) and optionally `strace -f -c -o strace_summary.log` for syscall breakdown.
5. **Cluster execution policy reminder** (from this repo's `CLAUDE.md`): never run compute on the login node directly — always via `srun`/`salloc` (as above) or `sbatch`. Read-only `ls`/`cat`/`head` on the login node is fine. Never write to raw `/tmp`; use `~/scratch/0.Phylophere/.tmp/` on correfoc.

---

## Tier 0 — Re-verify before touching anything (fidelity gate, ~1-2 hours)

**STATUS: DONE, 2026-09-16.** All four checks below closed without contradicting the NFS-latency
diagnosis — if anything, real production data shows the bottleneck is *worse* than the single-gene
estimate. **Go-ahead confirmed for Tier 1.**

### Findings

**Constraint that reshaped how this tier was executed**: at verification time, `mramon`'s SLURM CPU
quota on `correfoc` was already saturated — 12 `nf-CAAS_PERMULATION` jobs running (several 2-3h into
an 8-worker `CAAS_PERMS_DISAMBIGUATE_BATCHED` replay of the real cancer/neoplasia permulation-null run)
plus 12 more queued behind `QOSMaxCpuPerUserLimit`. Launching fresh isolated `srun` comparison jobs
(steps 1-2 as originally scoped) would have queued behind the user's own production run rather than
executing promptly, and would have contended for CPU with the very workload being profiled. Instead,
steps 1-2 were answered **passively, from the live production jobs themselves** — a strictly stronger
evidence source than new isolated single-gene tests would have been (real gene mix per batch, real
8-way contention, real NFS load from concurrently-running sibling jobs on the same filesystem).

**Step 1 (generalize beyond ABCB10)**: superseded by step 2's approach — each sampled production job
processes a *different* batch of 20 real genes drawn from the full ~16,100-gene set (not just ABCB10).
Six independent batches were sampled, i.e. ~120 distinct real genes, spanning alignment lengths from
564 (ABCB10, this run's smallest gene) up to 1567 residues (ABCB8/ABCC-family, this run's largest) —
confirmed via `ensembl_genes.output` for the batch's 20 genes. All show the same low-utilization
pattern; the diagnosis is not an ABCB10-specific artifact.

**Step 2 (does real 8-way parallelism hide the per-worker I/O latency?)** — **No. It does not, and
utilization is worse than the single-worker isolated estimate.** Measured via `sstat` (live per-job CPU
accounting, since `sacct`'s `TotalCPU` field does not populate for still-running jobs and `top`/`ps`/`cat`
aren't present on the compute image) against 6 independently-running `CAAS_PERMS_DISAMBIGUATE_BATCHED`
jobs, each with `AllocCPUS=8`:

| JobID | Elapsed | AveCPU (sstat) | Aggregate utilization (AveCPU / (Elapsed×8)) |
|---|---|---|---|
| 6713785 | 02:18:37 (8317s) | 00:12:39 (759s) | 1.14% |
| 6713793 | 02:15:29 (8129s) | 00:14:48 (888s) | 1.37% |
| 6713804 | 02:07:41 (7661s) | 00:18:48 (1128s) | 1.84% |
| 6713920 | 01:10:32 (4232s) | 00:20:00 (1200s) | 3.54% |
| 6713744 | 02:33:43 (9223s) | 00:14:18 (858s) | 1.16% |
| 6713688 | 02:55:26 (10526s) | 00:20:24 (1224s) | 1.45% |

Mean ≈ **1.75%** aggregate utilization across 8 allocated CPUs — i.e. real production 8-worker batches
are using on average less than 1/7th of a single core, worse than the ~9-10% figure in the original
(differently-measured) note and much worse than the ~5.5% single-worker/single-gene isolated estimate.
This is consistent with (not just "not contradicted by") the NFS-latency diagnosis: 8 workers hammering
the same per-cycle trait-file round-trip concurrently does not hide the per-worker latency behind
overlap — if anything the shared-filesystem contention across workers appears to make things slightly
worse, though this data can't cleanly separate "inherent per-worker latency" from "contention effect"
(would need single-worker vs. 8-worker on the *identical* gene batch, which the live-job approach
above cannot give — flagged as a residual gap, not worth closing before Tier 1 given the direction of
the effect is unambiguous either way).

**Step 3 (confirm file:line references haven't drifted)**: re-checked against `scoring_v2` HEAD
(`76b33f81d`, 2026-09-16) by direct read of both files. All cited references in the table above and in
Tier 1 below (`gene_wrapper.py:755-975`, `1123-1210`; `disambiguate_single.py:566-665`) still match —
no drift since the original investigation session.

**Step 4 (re-prioritize?)**: No. Tier 1 (eliminate the per-cycle trait-file round-trip in `_perms_worker`)
remains the correct first move, now with higher confidence than before — the effect is real, generalizes
across genes, and is not an artifact of single-worker measurement hiding behind parallelism.

Exit criterion met: diagnosis confirmed across ~120 real genes (6 production batches × 20 genes) and
both single-worker (original session) and real multi-worker (this session) measurements, without
consuming the user's already-saturated CPU quota.

### Tier 0 addendum — controlled single-vs-multi-worker test, 2026-09-16 (CPU quota freed)

The user stopped the production run, freeing the CPU quota, enabling the controlled apples-to-apples
test flagged above as a residual gap. Three genes spanning this batch's size range (ABCB10 564aa,
ABCC5 1439aa, ABCC8 1567aa; 300 cycles each, 900 total) were run two ways on identical inputs:

| Run | real | user+sys (CPU) | CPUs | Aggregate utilization |
|---|---|---|---|---|
| Sequential, 1 worker each (3 separate `srun --cpus-per-task=1` jobs, summed) | 15.801s | 4.821s | 1 each | ~30.5% (per-job, warm-FS) |
| Parallel, 3 workers, 1 job (`srun --cpus-per-task=3 --workers 3`, identical 3 genes + cycles) | **19.053s** | 4.122s | 3 | **7.2%** |

The parallel run was **slower in wall time than the sequential sum** despite a 3x larger CPU
allocation, and its aggregate utilization (7.2%) is ~4x worse than doing the same work one gene at a
time (30.5%). `stat()` call count scaled with the number of genes (57,784 vs. ~23,475×3), ruling out
any caching benefit from batching — the workers are not sharing filesystem work, they are contending
for it. This closes the residual gap flagged in the original Tier 0 write-up: real multi-worker
parallelism does not hide the per-cycle I/O latency, it actively degrades it, most plausibly via
concurrent NFS metadata contention (multiple workers issuing `stat`/`openat` against the same
mount concurrently creates queueing that a single worker never sees).

Note on interpreting absolute times across this addendum vs. the original single-gene ABCB10 numbers
at the top of this doc: today's filesystem is "system-wide warm" (12+ hours of the just-stopped
production run touched these same 20 genes' files extensively), so single-worker wall-times here
(~5s/300 cycles) are much faster than the original cold-cache ABCB10 measurement (36.5s/300 cycles).
The *relative* single-vs-multi-worker comparison above is unaffected by this (both runs share the same
warm state), but do not directly compare this addendum's absolute wall-clock numbers to the top-of-doc
cold-cache baseline — they measure different cache states, not different code.

---

## Tier 1 — Eliminate the per-cycle trait-file round-trip (the confirmed, highest-confidence fix)

**STATUS: implemented + local regression net green, 2026-09-16. Real-data golden-diff and before/after
timing (the two remaining exit-criterion items) are blocked on syncing this code to the cluster** — per
this repo's `CLAUDE.md`, code changes reach the cluster only through the user's own git sync process,
never by editing/copying files onto the cluster directly. Once synced, re-run the Tier-0-addendum
`srun` harness (ABCB10/ABCC5/ABCC8, single- and multi-worker) against the new code and diff
`perm_pos_detail`/`perm_pos_pval.tsv` output against a pre-change capture for the same genes/cycles.

Implementation landed in `gene_wrapper.py` (`build_cycle_inputs`, `_perms_worker`) and
`disambiguate_single.py` (`analyze_gene_disambiguation` gained a `trait_pairs=` kwarg that takes
precedence over `trait_file_path=`). `_write_cycle_trait_file` and the `cycle_traits/` work-dir were
deleted outright (dead once nothing writes to them). Regression net:
`subworkflows/CT_DISAMBIGUATION/local/src/convergence/test_perms_worker_trait_pairs.py` — 3 tests,
locally green, proving the in-memory `{1: list(zip(fg, bg))}` shape is exactly what
`parse_trait_pairs` used to derive from the file it replaces (including the FOP-mirror case where the
old trait filename embedded an "H\<n\>" token that changed the returned dict's *key* but never its
pairs — safe because `_resolve_contrast` never reads that key when there's exactly one contrast, which
is always true for a single replayed cycle).


**Target**: `_perms_worker`'s per-cycle loop currently does, for every `(gene, cycle)` pair it replays:
write once (`_write_cycle_trait_file`, `gene_wrapper.py:933-943`, called from `build_cycle_inputs`,
`gene_wrapper.py:946-975`) then re-read+re-parse via `analyze_gene_disambiguation(trait_file_path=...)`
→ `parse_trait_pairs(Path(...))` (`disambiguate_single.py:665`, itself calling
`loaders.py:parse_trait_pairs:446-503`, which does a real `.exists()` + `open()` per call, **no caching**).
The in-memory `(fg, bg)` species lists are already available from `_read_resample_labelings`
(`gene_wrapper.py:755-782`) *before* they get serialized to disk — the round-trip is pure overhead.

**Scope check (from Tier 0's signature lookup)**: this only affects `_perms_worker` (the null-replay
path). `process_single_gene` (the observed path, `gene_wrapper.py:208-225`) calls
`parse_trait_pairs` once per gene, not once per cycle — it has no repeated-reparse problem, but its call
signature into `analyze_gene_disambiguation` must keep working unchanged after this refactor.

### Implementation sketch

1. Add a new optional parameter to `analyze_gene_disambiguation`
   (`disambiguate_single.py:566-583`), e.g. `trait_pairs: Optional[Dict[int, List[Tuple[str, str]]]] = None`
   — the exact return type of `parse_trait_pairs`. At line 664-665, prefer `trait_pairs` when given,
   falling back to `parse_trait_pairs(Path(trait_file_path))` when only `trait_file_path` is supplied
   (keeps `process_single_gene`'s existing call working unchanged).
2. In `_perms_worker`'s per-cycle loop (`gene_wrapper.py:1182-1210`), replace the
   `trait_path = cycle_trait_files.get(cyc)` / `trait_file_path=Path(trait_path)` pattern
   (currently at line 1197) with a direct in-memory construction of the `trait_pairs` dict from the
   `(fg, bg)` lists already available via `_read_resample_labelings` — passed straight into
   `analyze_gene_disambiguation(trait_pairs=...)`, no file write, no file read.
3. `build_cycle_inputs` (`gene_wrapper.py:946-975`) and `_write_cycle_trait_file`
   (`gene_wrapper.py:933-943`) likely become dead code for the null-replay path once this lands — confirm
   nothing else in the FOP-mirror/`--cycles` CLI logic still depends on the on-disk trait files before
   removing them (check `disambiguation_perms_main.py`'s use of `--resample-dir`/`--cycles`, and whether
   `build_cycle_inputs`'s return values are consumed anywhere else).
4. Construct the `{pair_id: [(high_species, low_species), ...]}` shape directly from `(fg, bg)` — match
   `_write_cycle_trait_file`'s existing pairing convention (`fg[k]` ↔ `bg[k]` is pair `k+1`) exactly, since
   that convention is what `parse_trait_pairs` currently reproduces from the file; getting this pairing
   wrong would silently change scoring results.

### Regression net (must exist before merging, per Tier-0 finding that none currently does)

No test currently exercises `analyze_gene_disambiguation`, `parse_trait_pairs`, or `_perms_worker`
end-to-end (confirmed: only `test_core_v3_golden.py` and `test_pool_domains_golden.py` cover adjacent
functions; see `subworkflows/CT_DISAMBIGUATION/local/src/convergence/golden/`). Before refactoring:

1. Capture a golden fixture: run `_perms_worker` (or `analyze_gene_disambiguation` directly) for one
   small real gene and a handful of real cycles **on the current code**, save its output
   (`perm_pos_detail`/`gene_cycle_scores.tsv` rows) as a fixture, following the existing pattern in
   `subworkflows/CT_DISAMBIGUATION/local/src/convergence/golden/gen_golden.py`.
2. After the refactor, re-run the identical inputs and diff against the fixture — must be bit-for-bit
   identical (this is a pure I/O-path change, not a scoring change; any diff is a bug).
3. Also re-run the Tier-0 timing reproduction (same genes, same cycle counts, cold-cache) and compare
   wall-clock against the Tier-0 baseline numbers, to confirm the fix actually moved the needle before
   calling this tier done.

Exit criterion: golden diff is clean, and a cold-cache timing re-run on ≥1 of the Tier-0 genes shows a
measured wall-time reduction (not just a theoretical one).

---

## Tier 2 — Cut redundant `stat()` calls (secondary, contingent on Tier 1's residual profile)

The ~78 `stat()` calls/cycle and ~12.5 `openat()`/cycle (19% failing) measured in Tier 0 partly reflect
the trait-file round-trip Tier 1 removes — but not necessarily all of it (there may be other per-cycle
filesystem probing, e.g. in `find_gene_alignment`-style candidate-path checks, that Tier 1 doesn't touch).

1. **Re-run the `strace -c` profile from Tier 0 on the post-Tier-1 code.** If the stat/openat counts have
   dropped to near-zero per cycle, this tier is likely already resolved — stop here, don't chase further.
2. If a meaningful stat/openat count remains, trace it to its source (likely candidate-path probing
   somewhere in the gene/alignment lookup path, re-executed per cycle when it should be per-gene) and
   hoist it out of the per-cycle loop, same pattern as Tier 1.

Exit criterion: either "already resolved by Tier 1, no further action" or a second measured fix with its
own before/after timing.

---

## Tier 3 — `encoded_distribution` GEMM precompute (lowest confidence today, re-verify before investing)

This was the original hypothesis (per-gene tree-walk vectorization à la `perm_replay_vec.py`), and the
code-level design is sound (see below), but its measured basis is a **synthetic** cProfile run, not real
data, and Tiers 1-2 will change the shape of where real time goes. **Do not start implementation here
without re-profiling on real data first, now that Tiers 1-2 should make a real-data cProfile run fast
enough to actually do** (previously impractical at cold-cache NFS speeds).

1. **Re-profile `compute_domain_scores` with cProfile against one real gene's real cycles**, post-Tier-1/2,
   using the same `srun` harness as Tier 0 (swap `strace -c` for `python -m cProfile`). Confirm or revise
   the ~72% `encoded_distribution`/`encode_aa` figure against real posteriors, real domain counts (K),
   real scheme usage — not the synthetic guesses (fake tree, K=8, 3000 fake cycles) used originally.
2. **If confirmed as a real, meaningful fraction of remaining time**, the design (already sound
   independent of the profiling question):
   - Domain construction (`disambiguate_single.py:793-844`) confirms LCA/domain→MRCA assignment is
     per-cycle, not per-gene — so the GEMM precompute must be **per-gene, cycle-independent, on the
     posterior side only** (not an LCA cache, which would need to be per-cycle and is a much smaller win).
   - `encoded_distribution`/`encode_aa` (`path_scores.py:67-105`) is confirmed a fixed linear map:
     `group_vec = M_scheme @ aa_vec`, 5 schemes (`grouping.py:5-136`: US=20 singleton groups, GS1=6,
     **GS2=7** [corrected from earlier session's "6"], GS3=6, GS4=12 groups).
   - Posteriors are sparse/threshold-filtered (`posterior.py:263-326`) — pad-to-20 for the matmul, but
     keep the "unrecorded remainder" logic in `worst_case_any_group_probability` (`path_scores.py:130`)
     separate, since it depends on that sparsity.
   - Precompute slot: build `posterior_tensor` (`n_nodes × n_positions × 20`) once per gene inside
     `_load_gene_asr_context` (`gene_wrapper.py:1144-1157`, already gene-scoped and cycle-shared), compute
     `encoded = posterior_tensor @ M_scheme.T` once per `(gene, scheme)` (5 calls total), then index into
     the resulting dense array per `(node, position)` lookup instead of rebuilding the dict per call. What
     varies per cycle is only *which* `(node, position)` pairs get looked up (driven by that cycle's
     domain→MRCA assignment), not the tensor's contents.
3. **If the real-data profile shows this is no longer a meaningful fraction of wall time** (plausible,
   since Tiers 1-2 remove the I/O cost this profile never accounted for, and might reveal a different
   real bottleneck), stop here and re-scope based on what the real profile actually shows instead.

Exit criterion: a real-data cProfile confirming or revising the original percentage, and a go/no-go
decision on the GEMM implementation based on that — not the synthetic number.

---

## Not prioritized (flagged, not scheduled)

- `mp.Manager()` subprocess spawn/join overhead (~0.44s per job in the Tier-0 `strace` run, 3 `wait4`
  calls) — a fixed per-job cost, doesn't scale with cycle count. Worth a look only if job-count-level
  overhead becomes material (e.g. if batch count grows much further), not before.

---

## Verification checklist (applies across all tiers)

- Every tier's "exit criterion" above requires a **measured** before/after, not a code-read argument —
  this whole investigation went through one bad guess (LCA-fixed-per-gene) and one unverified synthetic
  number (72% in `encoded_distribution`) already; don't let a third one through.
- All cluster execution goes through `srun`/`sbatch`, never the login node directly, per this repo's
  `CLAUDE.md` cluster policy — reuse the reproduction recipe above.
- Any code change needs the golden-fixture regression net from Tier 1 to pass before merging, and ideally
  extended to cover Tier 2/3 changes too once they exist.
- This file's "measured vs. inferred" table should be updated in place as each tier's findings land, so
  a future conversation picking this up mid-way doesn't have to re-read the whole investigation history —
  just the current state of that table.
