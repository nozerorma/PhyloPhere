# CAAS pipeline + permulation-excess — runtime & scaling analysis

> Goal: know how long each stage takes, how it scales with **genes**, **cycles**, and
> **full_perms**, and where the optimization wins are for a typical production analysis
> (**16k genes, 1M cycles, 1000 full_perms**).

## Methodology

Per-task wall-times were extracted from Nextflow work-dir timestamps
(`mtime(.exitcode) − mtime(.command.begin)`), classified by `.command.sh` content,
deduped across runs. Reference scale = the integrated test: **G₀ ≈ 1939 genes,
C₀ = 1000 cycles, P₀ = 10 full_perms**, 35 species / 3 contrast pairs, precomputed ASR
cache 3.7 GB. Numbers are **task wall-seconds** (single task); real elapsed time =
`Σ task-seconds / effective parallelism`. Extraction scripts:
`scratchpad/extract_merged.py`.

Caveat: these are local-executor runs on one workstation; absolute seconds shift with
hardware and `cpus`/`maxForks`, but the **per-unit costs and scaling exponents** are what
matter for extrapolation.

## Measured per-stage cost (test: 1939 genes / 1000 cycles / 10 full_perms)

| Stage | tasks | batch | mean/task | **per-unit cost** | scales with |
|---|---|---|---|---|---|
| **DISCOVERY** | 78 | 25 genes | 56.6 s | **2.3 s / gene** | genes × positions |
| RESAMPLE | 1 | — | ~1–4 s | ~cycles (cheap here) | **cycles** |
| **BOOTSTRAP (significance)** | 194 | 10 genes | 23.4 s | **2.34 s / gene / 1000 cycles** | **genes × cycles** |
| CT_SIGNIFICATION | 1 | — | 3 s | flat | positions |
| **CT_DISAMBIGUATION (real)** | 1 | all genes | 276 s | **0.14 s / gene** (incl. ASR load) | genes × CAAS-pos |
| CT_POSTPROC | ~6 | — | 24 s | flat-ish | positions |
| **CT_ACCUMULATION** | ~3/dir | — | 82 s | per 1M randomizations | **randomizations × genes** |
| VEP | 1 | — | 166 s | flat | positions |
| SCORING | 1 | — | 16 s | flat | genes |
| **FCS** (per report) | 4–7 | — | 58 s | ~flat to N=1000 (vectorized) | **full_perms (null)** × GMTs |
| STRING (gated) | per slice | — | 413 s | expensive | gene-list size |
| ASR_ROBUSTNESS | 1 | — | 67 s | flat | genes |
| **PERMS · SUBSET_RESAMPLE** | 1 | — | 0.2 s | trivial | — |
| **PERMS · BOOTSTRAP (full-pool)** | **1939** | **1 gene** | 2.07 s | **0.21 s / gene / full_perm** (full pool) | **genes × full_perms** |
| **PERMS · CONCAT_DISCOVERY** | 1 | — | 6 s | rows | — |
| **PERMS · DISAMBIGUATE (replay)** | 1 | all genes | *pending¹* | **load-once ASR + genes×full_perms×score** | **genes × full_perms** |
| **PERMS · AGGREGATE** | 1 | — | ~3 s | rows | full_perms |

¹ The clean validation run reaches the replay shortly; analytical estimate below.
Update this row with the measured value (`extract_merged.py` → `PERMS_DISAMBIGUATE`).

### Per-stage SERIAL task-seconds at test scale (one full run)
DISCOVERY ≈ 4 460 · BOOTSTRAP(sig) ≈ 4 540 · BOOTSTRAP_PERMS ≈ 4 010 ·
CT_DISAMBIGUATION ≈ 276 · ACCUMULATION ≈ 250–2 000 · VEP ≈ 1 000 · FCS ≈ 250 ·
STRING ≈ 4 100 (if enabled) · everything else < 100 each.
**Total ≈ 20–25 k serial task-seconds** (~6 h serial; far less wall-time under
parallelism — e.g. ÷12 cores ≈ 30–40 min, dominated by the two bootstraps).

## Scaling laws (the three knobs)

```
DISCOVERY              ∝ genes
RESAMPLE               ∝ cycles
BOOTSTRAP (sig)        ∝ genes × cycles          ← the genome-wide significance null
CT_DISAMBIGUATION      ∝ genes  (+ one-time ASR load per gene)
ACCUMULATION           ∝ randomizations × genes
SCORING / SIGNIF       ∝ genes / positions  (flat-ish)
FCS                    ∝ full_perms (null, vectorized) — near-flat to N=1000
BOOTSTRAP_PERMS        ∝ genes × full_perms      (full pool: every position, N labelings)
PERMS_DISAMBIGUATE     ∝ genes × full_perms      (ASR loaded ONCE/gene, then N replays)
```

The two terms that explode are the ones with a **product of two knobs**:
`genes × cycles` (significance bootstrap) and `genes × full_perms` (the two perms stages).

## Extrapolation to production (16k genes, 1M cycles, 1000 full_perms)

Scale factors vs test: genes **×8.25**, cycles **×1000**, full_perms **×100**.

| Stage | scale factor | serial core-seconds | ≈ core-hours |
|---|---|---|---|
| DISCOVERY | ×8.25 | 4 460 × 8.25 ≈ 3.7e4 | 10 |
| RESAMPLE | ×1000 | ~1e3–1e4 (BM permulation, 1M) | 1–3 |
| **BOOTSTRAP (significance)** | **×8250** (8.25·1000) | 4 540 × 8250 ≈ **3.7e7** | **≈ 10 400** |
| CT_DISAMBIGUATION | ×8.25 | 276 × 8.25 ≈ 2.3e3 | 0.6 |
| **CT_ACCUMULATION** | ×8.25 (1M fixed) | ~1e4–1.6e4 | 3–4 |
| **BOOTSTRAP_PERMS (full-pool)** | **×825** (8.25·100) | 4 010 × 825 ≈ **3.3e6** | **≈ 920** |
| **PERMS_DISAMBIGUATE (replay)** | **×825** scoring (load ×8.25) | ≈ **2–6e5** (est.²) | **≈ 60–170** |
| FCS | ×100 (vectorized) | ~1e2–1e3 | <1 |
| SCORING / SIGNIF / POSTPROC | ×8.25 | <1e4 total | <3 |

² Replay estimate: the real disambiguation is 276 s for 1939 genes / 1 labeling. ASR
load dominates that cost and is paid **once per gene** in replay (the load-once
invariant), so replay ≈ `load(genes) + genes × full_perms × score_per_labeling`. The
full-pool labelings score more positions than the real CAAS-positive set, so treat this
as an upper-ish band until the measured number lands.

### Headline
- **Significance BOOTSTRAP (genes × cycles) dominates the entire pipeline** at 1M cycles:
  ~10 000 core-hours — two orders of magnitude above anything else. This is **pre-existing**,
  not introduced by the permulation feature.
- **The permulation feature adds ~1000 core-hours**, almost all in **BOOTSTRAP_PERMS**
  (genes × full_perms, full pool). The replay is comparatively cheap thanks to load-once ASR.
- At realistic parallelism (e.g. a 64-core node or a SLURM array of hundreds of tasks),
  the perms feature adds **single-digit to low-tens of wall-clock hours**; the significance
  bootstrap is what sets the multi-day floor.

## Where to optimize (ranked by payoff)

1. **Significance bootstrap at 1M cycles — the real wall (pre-existing).** `genes × cycles`
   in pure-Python `caasboot` per-cycle counting (`boot.py`). Options: (a) lower `cycles`
   where the empirical-p resolution allows (1M gives p-floor 1e-6 — rarely needed);
   (b) vectorize the inner per-cycle CAAS test (NumPy over the resample matrix) — this is
   the single biggest lever; (c) it is already batched, so scheduler overhead is not the issue.

2. **BOOTSTRAP_PERMS is PER-GENE (1939 → 16k tasks).** Biggest perms-side win: **batch it**
   like `BOOTSTRAP_BATCHED` (N genes/task with an internal worker pool). At 16k one-gene
   tasks the Nextflow scheduler + container spin-up overhead (~0.2–1 s/task) becomes a large
   fraction of the 2 s/task compute. Batching 10–25 genes/task cuts that 5–25×. (File:
   `subworkflows/CT/caas_permulation.nf` → add a batched variant mirroring `ct_bootstrap.nf`.)

3. **BOOTSTRAP_PERMS tests the FULL position pool every labeling.** Most positions never
   become CAAS under any labeling. A cheap pre-pass (union of positions that are CAAS in
   *any* of the N labelings) would shrink the per-labeling work; or share the resample
   parsing across labelings (already the case). Reducing `full_perms` from 1000 is the
   blunt lever — but note p.perm floor = 1/(N+1), so N≈200–1000 is the useful range.

4. **PERMS_DISAMBIGUATE — protect the load-once invariant.** The whole feasibility rests on
   one gene = one worker task that loads ASR once and replays all N labelings in-process
   (`process_all_genes_perms`/`_perms_worker`). Ensure `ct_disambig_max_tasks_per_child`
   and the Nextflow task granularity keep it that way; if a refactor ever made it
   (gene × cycle) tasks, ASR I/O would blow up by ×N. Parallelizing the N replays *within*
   a gene (threads over cycles on the cached posteriors) is a further win.

5. **ACCUMULATION (1M randomizations × genes)** is the third-heaviest at scale (~3–4
   core-hours) — already vectorized per gene; lower-priority.

6. **STRING** is expensive per report (~400 s) and grows with gene-list size — keep it
   gated (`scoring_string`) for production-scale runs.

## Practical recommendation
For a 16k-gene / 1M-cycle / 1000-full_perm run on a cluster: the significance bootstrap is
the schedule-defining stage (batch across an array; budget ~10k core-hours). The
permulation-excess null is a **~1000 core-hour add-on**, made tractable by (a) batching
BOOTSTRAP_PERMS and (b) the load-once ASR replay. Start `full_perms` at 200–500 to confirm
the p.perm block behaves, then scale to 1000 only if the analytic-vs-permulation gap warrants it.

---

# Threading, process model & BLAS — measured + roadmap (2026-06-27)

## What "too many python forks" actually is

Investigated live during a perms run. Findings:
- The CT Python kernel (`hyper.py`) uses **scalar `scipy.stats.hypergeom`** — no matrix
  algebra. `path_scores.py` (disambiguation) imports **no numpy/scipy at all** (pure-Python
  tree walks + dict lookups). Accumulation uses numpy only for **elementwise/reduction**
  ops (`percentile`, `random`, `sum(axis=0)`) — no `dot`/`matmul`/`linalg`.
- Steady-state, a `ct bootstrap` process has **1 thread and no BLAS lib mapped**. A
  transient 12-thread burst was observed (a brief library/auto-thread spin), not a
  persistent 8×12 bloom. So the earlier "thread oversubscription" framing was overstated.
- **The real driver of the "many forks" feeling is process COUNT**: BOOTSTRAP_PERMS is
  **per-gene** (1939 → 16k one-gene tasks) at `maxForks=8`, each a fresh interpreter that
  imports scipy + loads an alignment for ~2 s of work. Thousands stream past → churn.

## Levers (two complementary, both safe; neither global)

1. **Batch BOOTSTRAP_PERMS** (primary; addresses the churn). 1939 one-gene tasks →
   ~80–200 multi-gene tasks with an internal worker pool, mirroring `BOOTSTRAP_BATCHED`
   (`ct_bootstrap.nf` + `run_ct_bootstrap_batch.sh`). Amortizes interpreter/scipy/alignment
   startup across many genes per process.
2. **Pin BLAS/OpenMP threads to 1 in CT Python stages** (cheap safety net). Applied now to
   the two perms processes (`caas_permulation.nf` BOOTSTRAP_PERMS + CAAS_PERMS_DISAMBIGUATE)
   via `export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1`.
   **Recommended (not yet applied) for the shared CT Python stages**: the significance
   `BOOTSTRAP`/`BOOTSTRAP_BATCHED` (ct_bootstrap.nf), `DISCOVERY`/`DISCOVERY_BATCHED`
   (ct_discovery.nf), `CT_DISAMBIGUATION_RUN` (ct_disambiguation.nf — `mp.Pool`, so a
   per-worker bloom would multiply by worker count), `CT_ACCUMULATION_RANDOMIZE` (ctacc_run.nf).
   **MUST stay per-stage, not in `nextflow.config` `env{}`** — RERConverge does dense
   correlation BLAS (`permpvalcor`/`getStat`) and genuinely benefits from multithreaded BLAS;
   a global pin would slow it.

## Algorithmic vectorization (BLAS) roadmap — the real win

Principle: BLAS multithreading only pays if the hot loop is reshaped into a Level-3 BLAS
(matrix multiply). The CAAS "count members per group across many labelings" operation is a
**hidden GEMM** — the same pattern already proven in `fcs_enrich.R` (`M %*% Rk` cut the FCS
null from ~5 h to ~100 s, guarded by `fcs_enrich_equivtest.R`).

**1. Bootstrap (significance) + perms bootstrap — HIGH value (the genes×cycles / genes×full_perms wall).**
Per position `p`, scheme `s`: one-hot species→group `G_p ∈ {0,1}^(species×groups)`.
For a batch of `B` labelings, masks `F,Bg ∈ {0,1}^(B×species)`:
```
C_fg = F  @ G_p     # (B×species)·(species×groups) → (B×groups)   ← GEMM
C_bg = Bg @ G_p
```
CAAS condition per labeling = boolean reductions over `C_fg`/`C_bg`. Stack positions →
batched GEMM (`einsum 'bs,psg->bpg'`), one call replacing the `for cycle in cycles` loop.
Empirical p = fraction of `B` passing. **This collapses the per-gene fork churn into a few
big BLAS-parallel matmuls** (better core utilization *and* the threads finally do real work).
Memory: chunk over cycles×positions (~0.8 GB/chunk at 10k cycles × 500 pos × 20 groups).

**2. Accumulation aggregation — MED.** Per-draw sum over the sampled gene set is
`S @ v`, `S ∈ {0,1}^(draws×genes)` membership, `v` gene values (or × categories). Replaces
`np.sum([...], axis=0)` over 1M draws with one GEMM. Sampling stays `np.random`.

**3. Discovery — LOW-MED.** Same CAAS test but ONE labeling → tall-skinny GEMM; only
position-vectorization helps. Optional.

**4. Disambiguation — DON'T.** Irregular tree walk over dict posteriors (no dense matmul)
AND I/O-bound (3.7 GB ASR load dominates the 276 s, not arithmetic). Keep load-once ASR;
BLAS-ifying the cheap part saves little.

**Hard requirement for (1)/(2):** a rewrite of the core CAAStools detection kernel must
reproduce every per-position filter (gaps, missing data, conserved pairs, admissible
patterns, the 5 grouping schemes) as **masks**, and ship an **equivalence test** proving
bit-for-bit parity with the current loop (the `fcs_enrich_equivtest.R` pattern). This is
where calibration drift hides.

> This roadmap is being spun off into its own work session/worktree (multithread / BLAS
> vectorization). This document is the hand-off record.

---

# IMPLEMENTED — BLAS vectorization of the CAAS bootstrap (2026-06-29, worktree epic-rubin-f314db)

## What landed (shared CT kernel — serves significance bootstrap AND perms bootstrap)

1. **Vectorized CAAS/CAAP counting kernel** — `subworkflows/CT/local/modules/boot_vec.py`
   (`VectorizedBootstrap`). The per-cycle `for trait in filtered_traits` loop in
   `boot.caasboot()` is replaced by Level-3 BLAS: at a position the species→group one-hot
   `G_p` is shared across cycles, so per-cycle group membership is `C_fg = F @ G_p`,
   `C_bg = Bg @ G_p` (positions stacked column-wise into one `G_concat` → a single wide
   GEMM does the whole alignment). The CAAS condition, the 5 schemes (US/GS1-GS4), the
   gap/missing filters, the conserved-pair cap (`overlap ≤ max_conserved`) and the
   pattern-admission substring test are all mask reductions over the group axis. float32
   is exact (counts ≤ n_species). Cycles are streamed in chunks (`b_chunk`) to bound memory.

2. **Equivalence test** — `subworkflows/CT/local/modules/boot_vec_equivtest.py` (mirrors
   `fcs_enrich_equivtest.R`). Drives the scalar `caasboot()` and the kernel over 18 cases
   (clean, gappy, ambiguity codes, exotic symbols, numeric vs "NO" thresholds,
   `max_conserved` 1/3, every pattern subset incl. `1,2,3,4`, all 5 schemes, discovery-
   scheme subsetting incl. `CAAS=all`, multi-chunk accumulation) → **bit-for-bit parity,
   max|Δ|=0 on every case.** Run:
   `micromamba run -n phylophere python subworkflows/CT/local/modules/boot_vec_equivtest.py`

3. **Wired into `boot.boot_on_single_alignment`** (both directory and single-file modes,
   classical + CAAP) behind a guard: the vectorized path is used for the empirical-p COUNT;
   it falls back to the scalar walk when `--export_groups` / `--export_perm_discovery` are
   requested (those debug exports need per-trait detail) or if numpy is unavailable.
   Opt-out: `CT_BOOTSTRAP_VECTORIZE=0`. End-to-end output files are **byte-identical**
   between vectorized and scalar across dir/single × classical/caap.

   **Measured:** 35 species × 300 positions × 5000 cycles → scalar 316.95 s vs
   vectorized **0.58 s = 548× faster**, counts identical. (Single-thread BLAS; the kernel
   is now BLAS-bound so it scales further with `OPENBLAS_NUM_THREADS`.)

## Thread sizing (per-stage, never global env{} — RERConverge shares this OpenBLAS)
numpy here links `libopenblasp` (pthreads build → `OPENBLAS_NUM_THREADS` is authoritative).
- `ct_bootstrap.nf` **BOOTSTRAP**: `OPENBLAS_NUM_THREADS=${task.cpus}` (now BLAS-bound), others=1.
- `ct_bootstrap.nf` **BOOTSTRAP_BATCHED**: BLAS pinned to 1 (parallelism = the worker pool over genes).
- `ct_discovery.nf` **DISCOVERY/_BATCHED**: pinned to 1 (still scalar — discovery not vectorized).
- `ct_disambiguation.nf` **CT_DISAMBIGUATION_RUN**: pinned to 1 (pure-Python tree walk, mp.Pool).
- `ctacc_run.nf` **CT_ACCUMULATION_RANDOMIZE**: pinned to 1 per worker (np.random+bincount, no BLAS yet).

## NOT done (deliberately deferred — see notes)
- **Task 2 (accumulation `randomize.py`, MED):** the inner draw loop is **already** per-draw
  vectorized (`np.bincount`, shared memory); the `np.sum([...],axis=0)` the roadmap cites
  sums per-CHUNK results (~#workers items), not 1M draws. The remaining win is batching the
  `for r in range(chunk_size)` RNG draws, which **changes per-seed bit-reproducibility** of
  this calibrated permutation null (distribution stays identical). Recommend doing it as its
  own change gated by a **statistical-equivalence** test (mean/var/p within MC tolerance),
  not bit-parity. If done, then size its BLAS to `${task.cpus}`.
- **Task 4 (batch BOOTSTRAP_PERMS):** the perms workflow files (`caas_permulation.nf`,
  `BOOTSTRAP_PERMS`, the disambiguation replay) are **not committed in this branch/worktree**
  (this runtime doc is itself untracked in the main checkout). The 548× kernel already
  collapses the per-gene compute, so BOOTSTRAP_PERMS batching is far less urgent; wire the
  same `boot_vec` path + per-stage thread sizing into those processes when they land.
