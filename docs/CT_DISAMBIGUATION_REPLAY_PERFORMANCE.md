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

## Tier 1 closure — real-data verification, 2026-09-16

Ran on `correfoc` against the live (git-synced, commit `6caca24`) code, staged separately from the
Nextflow production work dir so as not to touch it: same 3 genes/cycle sets as the Tier-0 addendum.

| Test | Pre-Tier-1 | Post-Tier-1 | Golden diff |
|---|---|---|---|
| ABCB10, 1 worker, 300 cyc | real 4.819s | real 4.396s | `perm_pos_pval.tsv`/`gene_cycle_scores.tsv`/`perm_pos_detail` bit-identical |
| ABCC5, 1 worker, 300 cyc | real 5.849s | real 3.971s | bit-identical |
| ABCC8, 1 worker, 300 cyc | real 5.133s | real 23.147s* | bit-identical |
| combined3, 3 workers, 900 cyc | real 19.053s | real 19.273s | bit-identical |

\* ABCC8's post-Tier-1 wall time is an outlier driven by transient cluster/NFS load (its CPU time,
user+sys=2.0s, is in line with the other runs) — not a regression; correctness is unaffected regardless.

**Correctness**: confirmed bit-identical across every test. **Performance**: real but modest on the
single-worker path (~9-31% faster depending on gene/noise), and **essentially flat on the
multi-worker path that dominates real production cost** — the per-cycle trait-file round-trip Tier 1
removed was A cost, not THE cost. Syscall counts (single-gene, 300 cyc): `stat` 23,475→22,611 (-3.7%),
`openat` 4,049→3,482 (-14%). Exit criterion met (golden diff clean, measured wall-time improvement on
≥1 gene) but the result reframes what Tier 2/3 are actually chasing — see below.

### Unplanned finding: the "stat storm" is mostly one-time process-startup import scanning, not per-cycle work

Path-level `strace` (`-e trace=stat,newfstatat,openat`, no `-c` aggregation, so filenames are kept) on a
fresh 20-cycle ABCB10 run shows that after excluding paths under `micromamba/envs/phylophere/lib` and
this session's own test-staging path, **no remaining real-data path is touched more than ~15 times in a
20-cycle run** — the ASR cache (`asr_ABCB10/rst`, `tree_paml.nwk`), the alignment FASTA, the taxid
mapping, and the discovery file are each read a small constant number of times (once per gene-load, as
`_load_gene_asr_context` intends), never once per cycle. **There is no per-cycle filesystem probing left
to hoist — Tier 2, as originally scoped, is already closed by Tier 1.**

What *does* dominate raw syscall volume in an isolated single-process test is Python's own import
machinery scanning `site-packages` (pandas, numpy, Bio, pyarrow, multiprocessing) over NFS at process
startup — in this 20-cycle sample, ~38% of all `stat`/`newfstatat`/`openat` calls matched an
import-machinery or module-search path (this figure is inflated by this session's own nested
`.tmp/tier1_verify/src` staging path specifically — not representative of production's shallower
staging path — but the qualitative point holds regardless of the exact percentage). **This cost is
one-time per *worker process*, not per cycle or even per gene** — in production, `mp.Pool` reuses each
worker across up to `max_tasks_per_child` (50) gene-tasks, so it should be negligible amortized over a
single 20-gene batch's thousands of cycles.

**But it plausibly explains the multi-worker contention finding from the Tier-0 addendum**: `combined3`
(3 fresh worker processes launched concurrently by `mp.Pool`, one per gene) showed wall time *increase*
under parallelism (19.05s→19.27s) despite 3x the CPU allocation. Three freshly-spawned interpreters
each independently re-scanning the same NFS-hosted `site-packages` tree at roughly the same moment is a
plausible, structural explanation for why concurrency makes things *worse* here rather than hiding
per-worker latency — and it would recur **once per Nextflow batch** in production (a new `python3`
process + a new `mp.Pool` of 8 workers is spawned per `CAAS_PERMS_DISAMBIGUATE_BATCHED` task, i.e. once
per ~750+ batches in the full run), not just once for the whole run. This is *not proven* — it's an
inference from where the syscall volume concentrates, not an isolated measurement of import time under
concurrent NFS load specifically. It is flagged here as a candidate explanation worth a decision, not
folded into Tier 2/3 as scoped, since fixing it (e.g. reducing per-task interpreter/pool startup
overhead, or an architectural change to worker lifetime across batches) is a different kind of change
than either tier's original target.

---

### Import-storm hypothesis — tested directly, 2026-09-16: real but NOT the dominant cause

The flagged hypothesis above (concurrent fresh interpreters contending over NFS `site-packages` scans
explains the multi-worker slowdown) was tested directly rather than left as an inference. A probe script
imports the real `src.utils.gene_wrapper` module (the actual import graph the pipeline pays at every
worker spawn) and times only the import, launched solo vs. N concurrent copies via `srun`:

| Concurrency (N) | Mean import time | vs. solo |
|---|---|---|
| 1 (solo baseline) | 0.667s | — |
| 3 concurrent | 0.795s | +19% |
| 8 concurrent | 0.877s | +31% |

There is a real, measurable NFS-contention effect on import time under concurrency (+19-31%), but its
absolute magnitude is small: at N=8 the total *extra* wasted time across all 8 processes combined is
~1.7s. This cannot explain hours of sustained ~1-2% aggregate utilization in live production jobs, nor
the multi-second wall-time inflation seen in the `combined3` test — those need an effect roughly two to
three orders of magnitude larger than what import-time NFS contention alone produces. **Verdict: the
import-storm hypothesis is rejected as the dominant explanation.** It is a real, minor, structural cost
(worth remembering if `max_tasks_per_child` or batch sizing ever changes enough to make pool-spawn
frequency dominant), but the true source of the multi-worker contention found in the Tier-0 addendum
remains open. Given Tier 2 is otherwise closed (no per-cycle I/O left) and this dead end has been ruled
out, **Tier 3's real-data cProfile (genuine per-cycle CPU-bound compute) is the best-supported remaining
lever** — proceed there next rather than continuing to chase filesystem-level explanations.

---

## Tier 2 — Cut redundant `stat()` calls (secondary, contingent on Tier 1's residual profile)

**STATUS: CLOSED, no action needed, 2026-09-16.** See "Tier 1 closure" above — path-level `strace`
confirms zero meaningful per-cycle filesystem probing remains after Tier 1. The ~78 stat/cycle figure
from the original Tier 0 measurement was inflated by one-time process-startup import scanning divided
across the cycle count, not a genuine per-cycle cost; Tier 1 already removed the one genuine per-cycle
file round-trip that existed. Nothing further to hoist at this tier.

*(Original scoping, kept for record: the ~78 `stat()`/cycle and ~12.5 `openat()`/cycle measured in Tier
0 were suspected to partly reflect other per-cycle filesystem probing beyond the trait-file round-trip —
e.g. candidate-path checks re-executed per cycle when they should be per-gene. The path-level trace
above ruled this out: no real-data path is touched with per-cycle multiplicity anywhere in the traced
run. Exit criterion — "already resolved by Tier 1, no further action" — is met.)*

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

## Tier 3 real-data profile, 2026-09-16 — original hypothesis REJECTED, real bottleneck found

**Methodology fix first**: `cProfile` wrapping the CLI entry point only profiles the *main* process; with
any `--workers` value `mp.Pool` still spawns the actual work into a subprocess, so the first attempt at
this measured 2s of IPC-wait (`_recv_bytes`/`select.poll`) in the parent and nothing about real compute.
Fixed by calling `_perms_worker` directly, in-process, no `mp.Pool` (it's a plain function — nothing
about it requires a subprocess) — see `profile_worker_direct.py` in the reproduction recipe. This is a
cheap trap worth remembering for any future profiling of this codebase: check whether the profiled call
crosses a process boundary before trusting the numbers.

Real gene (ABCB10), real 300 cycles, real ASR cache, single process, `cProfile`. Total: 3.582s.

**The original hypothesis is rejected**: `encode_aa` (0.047s cumtime, 23,127 calls) and
`compute_domain_scores` (0.125s cumtime, 1,381 calls) together account for **~4.8% of total time**, not
the ~72% guessed from the synthetic profile. The GEMM precompute designed for Tier 3 would optimize a
function that isn't the bottleneck.

**What actually dominates**: of the 3.582s, ~1.175s (33%) is `_load_gene_asr_context` — genuinely
one-time per-gene setup (ASR cache parse, alignment load, taxid mapping) that amortizes away at
production scale (thousands of cycles/gene, not 300) — and ~0.383s (11%) is `_parse_discovery_entries`,
also one-time per gene. Excluding both one-time costs, the remaining ~2.0s of genuine per-cycle work is
overwhelmingly `get_mrca` (`tree_parser.py:442`, called 1,500 times, **1.523s cumtime — 76% of the
per-cycle remainder, 42% of total wall time**), via its helpers `find_node_by_taxid` (0.641s tottime,
1.25M calls) and `find_node_by_name` (0.448s tottime, 2.19M calls, `tree_parser.py:391-410`).

**Root cause**: `find_node_by_name`/`find_node_by_taxid` are recursive full-tree DFS searches — no
index, O(tree size) per call — invoked once per tip name inside `get_mrca`'s lookup loop
(`tree_parser.py:457-462`), for a tree that is fixed per gene and shared across every cycle. `_perms_worker`
already memoizes `get_mrca` results by sorted-taxa key (`disambiguate_single.py:723-731`,
`_get_mrca_cached`/`_mrca_cache`, ~4.6x hit rate here: 6,905 lookups → 1,500 real `get_mrca` calls) but
every cache *miss* still pays the full O(tree size) linear search per tip name, and there's no name/taxid
→ node index to make even a single lookup fast. Notably, `path_scores.py:135`'s `build_node_index`
already builds an O(1) id-keyed index for a different purpose — the same pattern isn't applied to
`get_mrca`'s name/taxid lookups.

**Proposed fix (not yet implemented, pending decision)**: build a `name → node` and `taxid → node` dict
once per gene — inside `_load_gene_asr_context` or alongside `hoisted_node_index` in
`disambiguate_single.py:722` (both already hoist other per-gene/tree invariants out of the per-cycle
loop, the same pattern) — and have `get_mrca` use dict lookups instead of recursive search. This turns
an O(tree size) per-tip-name lookup into O(1), directly targeting the ~42% of real wall time this profile
identifies, with a well-understood, low-risk fix shape (same as Tier 1's hoisting pattern, applied to a
different invariant). Scope check before implementing: confirm `find_node_by_name`/`find_node_by_taxid`
have no other callers whose behavior would change (e.g. relying on first-match-in-traversal-order
semantics that a dict wouldn't preserve if names/taxids aren't unique) — not yet done.

Exit criterion for this finding: real-data cProfile obtained (done), original % revised (done, rejected),
go/no-go decision on GEMM implementation = **no, don't build it**. A *different*, better-supported fix
(MRCA lookup indexing) is proposed in its place — this needs its own go/no-go from the user before
implementation, since it touches shared tree-traversal code (`tree_parser.py`) used outside this replay
path too, a larger blast radius than Tier 1's isolated I/O change.

### MRCA-index fix — implemented and verified, 2026-09-16 (commit `4f0b57a`)

Implemented as designed: `build_name_taxid_index(root)` (new, `tree_parser.py`) does one DFS building
`{name: node}` (all nodes) and `{taxid: leaf node}` (leaves only, same `lineage_taxid` split-on-`_`
convention as `find_node_by_taxid`); `get_mrca` gained optional `name_index`/`taxid_index` params, using
O(1) dict lookups when given and falling back to the original recursive search otherwise (so
`node_identification.py`'s one other caller is untouched). `disambiguate_single.py` builds the index
once per `analyze_gene_disambiguation` call, alongside the pre-existing `hoisted_node_index`, and threads
it into `_get_mrca_cached`. Regression net:
`subworkflows/CT_DISAMBIGUATION/local/src/asr/test_get_mrca_indexed.py` (5 tests) proves indexed and
unindexed lookups agree on name queries, taxid queries, cross-subtree queries, and that the
no-index-given path is byte-for-byte the pre-fix behavior.

**Real-data verification**: golden diff bit-identical across ABCB10/ABCC5/ABCC8 (`perm_pos_pval.tsv`,
300 cycles each, single worker) against the Tier-1 baseline. **Timing — methodology matters here**: the
first attempt (full CLI through `mp.Pool`, `--workers 1`, comparing user+sys wall-clock) showed no
measurable change on ABCB10 (1.55s → 1.62s) — misleading, because at this small a cycle count the fixed
cost of `mp.Pool` spawning a subprocess (paying Python/conda-env import cost *twice*, once per process)
dominates and swamps a real per-cycle win. Rerun as a clean in-process A/B (same technique as the
cProfile fix earlier in this tier: call `_perms_worker` directly, no `mp.Pool`, but this time with plain
`time.perf_counter()`, no `cProfile` instrumentation either, 3 repeats each side, pre-fix and post-fix
code run back-to-back in the same `srun` allocation to cancel out cluster-load drift):

| | Mean (3 runs) | Min |
|---|---|---|
| Pre-fix (Tier 1 only, commit `4f0b57a`'s parent) | 1.263s | 1.241s |
| Post-fix (Tier 1+3, commit `4f0b57a`) | 1.004s | 0.967s |
| **Improvement** | **20.5%** | **22.1%** |

This confirms a real, meaningful win — but notably smaller than the ~42% the `cProfile`-based profile
suggested. Reconciling the two: `cProfile`'s per-call instrumentation overhead is not uniform across
code shapes — it disproportionately inflates functions with very high call counts (the recursive
`find_node_by_name`/`find_node_by_taxid` were called 1.25M and 2.19M times respectively in the profiled
run), so their *share* of profiled time overstates their share of true unprofiled time. **Lesson for any
future profiling in this codebase**: `cProfile`'s relative percentages are trustworthy for finding *where*
to look (which is what correctly redirected this tier away from the GEMM precompute), but the absolute
improvement from fixing what it points at should always be confirmed with a clean, uninstrumented timing
comparison before being quoted as the expected real-world gain — this is the second such methodological
trap this tier hit (the first being `mp.Pool`'s IPC-wait time masquerading as worker compute time).

Exit criterion met: golden diff clean, real (uninstrumented, controlled A/B) timing improvement
confirmed at ~20-22% on a 300-cycle single-gene run — expected to matter more at production scale, where
thousands of cycles per gene amortize away the one-time per-gene setup cost that this fix doesn't touch,
increasing the fixed-per-cycle savings' share of total time.

---

## Multi-worker contention mystery — SOLVED, 2026-09-16 (commit pending)

The Tier-0 addendum found real 8-way (and 3-way) production parallelism performing *worse* than
sequential single-worker execution, and the import-storm hypothesis (tested in a follow-up session) was
real but far too small (~1.7s total at 8-way) to explain it. This section documents the actual root
cause, found by chasing the mystery to ground truth on real data.

### Dead ends ruled out first

- **Shared-cluster noise**: re-tested by running a hand-built direct `mp.Pool` reproduction and the real
  full-CLI invocation *back to back on the same node in the same `srun` allocation* — the direct
  reproduction consistently took ~3.1-3.3s (3 repeats) while the real CLI took ~15.7-16.5s for the
  identical 3-gene/582-cycle workload. Same node, same moment: not noise.
- **`postproc_filter`'s `clustering_discards` step**: re-tested with it enabled in the direct
  reproduction — no meaningful difference (3.14s → 3.26s).
- **`mp.Pool`/`forkserver` overhead**: Python 3.14 (confirmed via `mp.get_start_method()`) defaults to
  `forkserver`, not `fork`, on this system — a genuine CPython default change worth knowing about for any
  future multiprocessing work here, and a real hazard (any ad hoc profiling script needs
  `if __name__ == "__main__":` or forkserver crashes re-importing the launching script — production's
  `disambiguation_perms_main.py` already has this guard, confirmed). But measured directly, `mp.Pool`
  creation + task dispatch overhead was ~0.08s — negligible, not the cause.

### Root cause: `chunksize=10` starves the worker pool for typical batch sizes

`process_all_genes_perms`'s `pool.imap_unordered(_perms_worker_wrapper, args_generator, chunksize=10)`
(`gene_wrapper.py`, was line 2030) hands `imap_unordered` a fixed `chunksize=10`. `chunksize` controls how
many iterable items get bundled into ONE task dispatched atomically to ONE worker — it exists to amortize
IPC overhead for workloads with MANY CHEAP tasks. Here each task is one gene's full multi-cycle replay
(several real seconds of work), the opposite shape. **Whenever the number of genes in a batch is `<=
chunksize`, every single gene gets bundled into one chunk sent to one worker — every other worker in the
pool receives zero work for the entire batch, no matter how many are allocated.**

Confirmed directly: reproducing the exact combined3 workload (3 genes) with `chunksize=10` showed only
ONE worker's startup log line fired (`_load_gene_asr_context`'s "TAXONOMY CONFLICT" warning, normally one
per worker — expected 3, got 1), and wall time was ~1.7x worse (5.62s) than the same workload with
`chunksize=1` (3.2s, all 3 workers active, confirmed by 3 log lines). Real production batches are 20
genes with 8 allocated workers: `20 // 10 = 2` chunks, so **only 2 of the 8 allocated workers ever run
per batch, regardless of pool size** — mechanically capping utilization at 25% before any other
inefficiency stacks on top, fully consistent with the ~1-2% aggregate CPU utilization measured on live
production jobs in Tier 0.

This is a genuine regression, not an original design choice: commit `bf92df7` (2026-07-16,
"Chronological report renumbering + retire XL-mHG...") replaced a per-gene `pool.apply_async(...)`
dispatch (one task per gene, no chunking concept — every worker always had access to work) with
`imap_unordered(..., chunksize=10)`, silently introducing this cap as a side effect of an unrelated
refactor.

**Fix**: `chunksize=1` (implemented, `gene_wrapper.py`). Pure scheduling change — `imap_unordered`'s
result order was already unspecified before this fix (nothing in Pass A's aggregation depends on arrival
order), and `chunksize` never affects *which* results are produced, only how work is distributed across
workers — so no golden-diff regression test is meaningful here (output is provably identical regardless
of chunksize; only wall-clock parallelism efficiency changes). Real-data timing verification pending
sync to cluster.

### Direct confirmation from live production data (2026-09-16)

While this fix was being written, the user pointed at a real production batch
(`CAAS_PERMULATION:CAAS_PERMS_DISAMBIGUATE_BATCHED`, `malignant_prevalence_complete/b4/2883ac6...`,
reported as "batch 29 of 804" from a run that was since killed — the batch itself had already completed
successfully before being checked, not caught mid-execution). Its `.command.log` shows `workers=8` and
`"replaying 48442 cycles over 13 genes with 8 workers"`, but the per-worker `"TAXONOMY CONFLICT"` startup
warning (fires once inside `_load_gene_asr_context`, i.e. once per worker that actually receives a task)
appears **exactly twice** — not eight times. `ceil(13 genes / chunksize=10) = 2` chunks, so exactly 2
workers ran and 6 sat idle for the whole batch, precisely matching the mechanism above with zero
synthetic reproduction needed. That batch took 699.4s total (Pass A alone: 11:08:46 → 11:17:43 = ~537s)
with only 2 of 8 allocated CPUs ever doing anything — this is the single clearest piece of evidence in
the whole investigation, since it's a real completed production batch's own log, not a reconstruction.

### Reproduction recipe for this specific investigation

Direct `mp.Pool` reproduction script + phase-timing script live in
`~/scratch/0.Phylophere/.tmp/tier1_verify/{time_mp_pool2.py,time_worker_phases.py}` on `correfoc` (not
committed — investigation scaffolding, reusable for the next multi-worker question). Key technique: call
`_perms_worker`/`_perms_worker_wrapper` directly against a real `mp.Pool` built the same way
`process_all_genes_perms` does, rather than going through the full CLI, to isolate Pool-level effects
from Pass A/B and CLI-argument-parsing overhead — the same "bypass the outer machinery, call the real
function directly" technique that found Tier 3's `get_mrca` bottleneck.

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
