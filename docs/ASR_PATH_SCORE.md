# The ASR Path Score

*A single, continuous measure of how strong and trustworthy a convergent
amino-acid signal is, derived directly from ancestral state reconstruction (ASR).*

This document is layered on purpose:

- **§1–§3** are plain-language and need no evolutionary background — start here to
  explain the idea to a cancer-genomics audience.
- **§4–§6** are the formal model and worked examples — for supervisors and anyone
  who wants the exact definitions.
- **§7** records *why* each design choice was made (the questions reviewers ask).
- **§8** is a glossary of the few evolutionary terms used.

---

## 1. TL;DR

For each amino-acid position our pipeline flags as phenotype-associated (a **CAAS** —
Convergent Amino Acid Substitution), the ASR path score asks one question:

> *Are the changes that define this position genuinely independent, clean,
> convergent evolutionary events — or are they an artifact of shared ancestry,
> reconstruction noise, or a single non-independent origin?*

It returns one number in **[0, 1]**, built by multiplying **five factors** read
straight from the ASR posteriors. High = several lineages independently changed to
the same residue, from a clean ancestral background, with no sign that the change
predates their shared ancestors. Low = the apparent signal dissolves under any of
those checks.

It replaces three older, partly-overlapping scores (a binary ancestry gate, a
"convergence" score, and a "parallelism" score) with **one continuous signal**.

At scoring time the path score is combined with the permulation (phenotype) signal
into the final **CAAS score** (§6); the hypergeometric p-value is used only as a
significance **gate**, not as part of the score.

---

## 2. The problem, in cancer-genomics terms

Imagine you find a recurrent amino-acid change that shows up across many species
that share a trait (say, exceptional cancer resistance), and *not* in the species
that lack it. That recurrence is suggestive — but recurrence alone is weak evidence,
for exactly the reasons a recurrent mutation across tumors can be misleading:

| In tumor cohorts you worry about… | Here the analogous worry is… |
|-----------------------------------|------------------------------|
| Is the recurrent mutation a **driver** or a **passenger** that's just mutable? | Is the residue under selection, or just an easy mutation (a "hotspot")? |
| Are these **independent** tumors, or the same clone sampled twice? | Are these **independent** lineages, or one ancestral change inherited by descendants? |
| Did the samples really reach the **same** alteration? | Did the lineages converge on the **same** residue, or just each change to something? |
| How many independent samples carry it? | How many independent lineages actually changed (not merely inherited the state)? |

The ASR path score grades the **quality and independence** of the evidence. The
*strength of the phenotype association* (how surprising the partition is) is graded
separately by the permulation p-value, and the two are multiplied at the end (§6).

---

## 3. The intuition: five checks

Every CAAS position is defined by a set of **pairs** (independent phylogenetic
contrasts: one trait-positive lineage vs one trait-negative lineage sharing a
recent common ancestor). For each pair we know the present-day residues at the tips,
and ASR gives us the reconstructed posterior over residues at every ancestral node
back to the root.

The score asks five things and multiplies them:

1. **Replication of Private Isolation** (`core`) — *Did we observe at least two
   independent, clean mutations?* For each changed pair, we walk the private segment
   of its lineage (from MRCA up to its nearest LAC merge point) and compute the
   probability that it mutated independently: $\prod (1 - P(\text{derived}))$. We then
   compute the exact probability that **at least two** of these independent changes
   occurred ($P(\ge 2) = 1 - P(0) - P(1)$). Single changes score $0.0$.

2. **Independence at shared ancestors** (`independence`) — *Where these lineages
   last shared an ancestor, was the derived state already present?* We find the
   **LAC** nodes (the merge points of the changed pairs' MRCAs) and ask, at each, for
   the probability the derived state was already there. If it was, the pairs below it
   inherited the change together rather than converging independently — and the score
   craters.

3. **Parallel vs. non-parallel origins** (`mrca_diversity`) — *Did the lineages
   start from the same ancestral residue or different ones?* Changes from
   **different** starting points that land on the same residue are stronger evidence
   of selection (independent routes) than changes from an identical starting point
   (which could be a shared mutational tendency). Measured by the dot product of the
   MRCAs' posteriors, representing the exact probability of independent starts.

4. **Convergence agreement** (`derived_agreement`) — *Did the lineages reach the
   same residue?* If, on one phenotype side, the changed pairs split across several
   derived residues, that is divergence, not convergence, and the signal is divided
   across them.

5. **Conserved-pair confirmation** (`conservation_gate`) — *Did the pairs that
   did **not** acquire the derived residue hold the ancestral state deeply?* A pair
   whose tips both kept the ancestral residue confirms the contrast if that residue
   runs deep, and weakens it if the clade actually drifted. It can only
   confirm/undermine, never inflate.

---

## 4. The formal model

Notation: for a position with pairs indexed *i*, each pair has a most-recent common
ancestor (**MRCA**) node with a reconstructed posterior (its modal residue is the
*ancestral* state, ``focal_state``), and tip residues on each phenotype side. All
residues are compared in the active grouping scheme's encoded space (US = exact
residue; GS schemes = biochemical group), so "same residue" means "same group" under
a GS scheme, and ``P(group)`` sums the posterior over residues in that group.

A **changed side** is a phenotype side whose tip residue differs from the pair's
MRCA. A **changed pair** has at least one changed side. A **conserved pair**
(designated in metadata) has neither side change — both tips kept the ancestral
residue.

### 4.1 Per-pair private isolation (`pair_scores`)

The **LAC nodes** are the merge points of the changed pairs' MRCAs — the internal
nodes of the minimal subtree connecting them, found as the pairwise lowest common
ancestors. There are at most *n − 1* of them for *n* changed pairs.

For each changed pair, we identify its nearest LAC node (the deepest merge point
on its path to root). For each changed side of that pair, we walk the **private
segment** of its lineage (from the node directly above the MRCA up to, but
excluding, the nearest LAC node) and compute:

```
pair_side_score = ∏_{node ∈ private segment}  ( 1 − P(derived at node) )
```

- If the private segment is empty (e.g. sibling merge), the product evaluates to
  **1.0** (no private contamination).
- If the *modal* residue at hop 1 (directly above the MRCA) already equals the
  derived residue, a **`contaminated`** flag is raised.
- A pair's score $p_i$ is the mean of its changed sides' private scores.

### 4.2 Combinatorial convergence probability (`core`)

Rather than taking the mean of pair scores and multiplying by a heuristic count
factor, we compute the exact probability that **at least two** independent changes
occurred. This naturally incorporates count-dependence (single-change cases score
$0.0$):

$$P(0) = \prod_{i=1}^n (1 - p_i)$$

$$P(1) = \sum_{i=1}^n \left( p_i \prod_{j \neq i} (1 - p_j) \right)$$

$$\text{core} = 1 - P(0) - P(1)$$

### 4.3 Independence at shared ancestors (`independence`)

At each of the shared merge points (LAC nodes), the derived state must be absent
for the changes below it to be independent:

```
independence = ∏_{L ∈ LAC nodes}  ( 1 − Σ_{d ∈ D} P(d at L) )
```

Because it ranges over the few LAC nodes — not the full root paths — the product is
**depth-independent**. It reads ``P(derived)`` directly, with no ancestral term, so a
node sitting on a *third* residue (neither ancestral nor derived) simply lowers
``P(derived)`` instead of being misread (§7.5).

### 4.4 Parallel axis (`mrca_diversity`)

Compare the **posteriors** at the changed pairs' MRCAs in encoded space using the
**dot product** of their distributions:

```
overlap(a,b)   = Σ_g ( P_a(g) · P_b(g) )
diversity      = 1 − mean_{i<j} overlap( MRCA_i , MRCA_j )      # 0..1
diversity_mult = floor + (1 − floor) · diversity               # floor = 0.75
```

Identical MRCA posteriors (pure parallel) → `diversity 0` → `diversity_mult = floor`;
disjoint posteriors (independent starting points) → `diversity 1` → `diversity_mult = 1`.

### 4.5 Convergence axis (`derived_agreement`)

Per phenotype side, among the pairs whose tip on that side changed, count the
**distinct** encoded derived residues:

```
derived_agreement = mean_over_changed_sides ( 1 / n_distinct_derived )
```

Computed **per side**, so a pair changing on *both* sides (e.g. trait1→S, trait0→T)
is **not** penalised.

### 4.6 Conserved-pair axis (`conservation_gate`)

For each conserved pair, score how deeply its ancestral state runs to the root
(unsigned conservation = unweighted mean of ``P(ancestral)`` along the path). Then

```
conservation_gate = 0.5 + 0.5 · mean( conservation-to-root of conserved pairs )
                  = 1.0   if there are no conserved pairs (novel position)
```

A deeply conserved pair confirms the contrast (→ ~1.0); a pair that drifted weakens
it (down to 0.5) but never annihilates the real convergence. The gate can only
confirm or undermine, never boost.

### 4.7 The position score

```
asr_path_score = clip( core · independence · diversity_mult · derived_agreement · conservation_gate , 0, 1 )
```

---

## 5. Worked examples (schematic, real implementation outputs)

Three sibling pairs, MRCAs just below a clean ancestral background. `floor = 0.75`, US scheme.

**Clean convergence** — MRCAs `A/A/A`, all top tips → `S`, background solidly `A`:
`pair_scores = [1, 1, 1]`, `core = P(>= 2 changes) = 1.0`, `independence = 1`
(LAC parent 20 is `A`), `diversity_mult = 0.75` (parallel starting states),
`derived_agreement = 1`, `conservation_gate = 1` → **asr = 0.75**.

**Inherited-state trap** — top tips `A/S/S`, MRCAs `A/A/S`: pair 3's `S` was
inherited (its MRCA is already `S`), so only **one** real change (pair 2).
`pair_scores = [1]`, `core = P(>= 2 changes) = 0.0` → **asr = 0.0**. The position has
no convergence at all.

**Contaminated shared ancestor** — MRCAs `A/A/A`, top tips → `S`, but the pairs'
merge node (LAC) already shows `S` with $0.8$ posterior: `core = 1.0`, but
`independence = 1 - 0.8 = 0.2` → **asr = 0.20**. The convergence is unmasked as
likely inherited.

---

## 6. How it flows through the pipeline

1. **CT_DISAMBIGUATION** computes the score per position (per biochemical scheme):
   [`src/convergence/path_scores.py`](../subworkflows/CT_DISAMBIGUATION/local/src/convergence/path_scores.py),
   wired in
   [`disambiguate_single.py`](../subworkflows/CT_DISAMBIGUATION/local/src/convergence/disambiguate_single.py).
2. **CT_POSTPROC** passes these columns through unchanged into `filtered_discovery.tsv`.
3. **SCORING** ([`scoring_compute.R`](../subworkflows/SCORING/local/scoring_compute.R))
   builds the final **CAAS score** from **two** orthogonal axes on raw [0, 1] scales:

   ```
   phen_score = 1 − pvalue_boot          # phenotype evolution (permulation evidence)
   caas_row   = phen_score · asr_path_score   # position evolution × phenotype evolution
   CAAS_score = scheme-weighted sum of caas_row   (US = GS1..GS4 = weight 1)
   ```

   The hypergeometric `Pvalue` is **not** in the score. It is a three-way
   significance **gate** instead: `gate_all` (every scored position), `gate_sig`
   (Pvalue < 0.05), and `gate_fdr` (BH-adjusted Pvalue < 0.05).

---

## 7. Design rationale (the "why")

### 7.1 Why no hop weighting?
PAML's posterior uncertainty *already* grows toward the root, so a hop weight
penalised deep nodes twice. The walk is now an unweighted product; depth robustness
comes from ranging only over private segment nodes (for `pair_scores`) and shared LAC
nodes (for `independence`).

### 7.2 Why score only the changed sides of participating pairs?
Conservation structurally out-scores isolation, so crediting every conserved side
would let a clade that merely *held still* inflate a position above clades that
independently *made the change*. So only changes are scored — except conserved pairs,
which feed the gate (§7.4).

### 7.3 Why is the parallel axis count-independent?
Replication is scored explicitly by the combinatorial `core` probability ($P(\ge 2)$)
and by the permulation p-value. If `mrca_diversity` also rewarded count it would
triple-count; so it asks only *"how distinct are the ancestral starting points?"*,
using the dot product probability of independent starting states.

### 7.4 Why is the conserved pair a gate, not an averaged term?
A conserved pair is a contrast that did **not** change. Averaging it in would
re-inflate the problem of §7.2 and let a conserved pair *raise* a score. Framed as a
noisy-sensor gate ($0.5 + 0.5 \times P_{\text{anc}}$), the 0.5 floor represents a 50%
prior that ASR reconstruction noise is responsible if conservation drops to 0.

### 7.5 Why the LAC independence product?
The LAC product isolates exactly the merge points — the nodes that decide whether the
pairs converged *independently* — and multiplies them, so one contaminated shared
ancestor craters the score. Reading `P(derived)` directly (no ancestral term) also
resolves the "third residue at a node" edge case: such a node lowers `P(derived)`
rather than being mistaken for ancestral.

### 7.6 Why the combinatorial $P(\ge 2)$ core score?
Instead of mixing joint probabilities with an arbitrary heuristic count multiplier
like $n / (n+k)$, we treat each pair's private isolation score $p_i$ as the probability
of a true independent mutation. Computing $P(\ge 2)$ gives the exact mathematical
probability that convergent evolution (at least two independent changes) actually
occurred at this position.

---

## 8. Glossary

- **CAAS** — Convergent Amino Acid Substitution: a position where trait-positive
  species share one residue and trait-negative species share another, beyond chance.
- **ASR** — Ancestral State Reconstruction: a statistical estimate (here, PAML) of the
  amino acid at each internal tree node, with a posterior probability.
- **MRCA** — Most Recent Common Ancestor: the internal node where a pair's two
  lineages join; its modal residue is the pair's ancestral state.
- **LAC node** — a merge point of the changed pairs' MRCAs: a node where two or more
  pairs' lineages last shared an ancestor. The independence product is taken over
  these.
- **Pair / contrast** — one trait-positive lineage vs one trait-negative lineage
  sharing an MRCA; the independent unit of evidence.
- **Conserved pair** — a pair that did not acquire the derived residue (both tips kept
  the ancestral state); feeds the conservation gate, never the convergence signal.
- **Parallel vs convergent change** — *parallel* = independent lineages reach the same
  residue from the **same** ancestral residue; *convergent (non-parallel)* = from
  **different** ancestral residues. Non-parallel is stronger evidence of selection.
- **Permulation** — a phylogenetic permutation of the phenotype on the tree;
  `pvalue_boot` is its p-value, the phenotype-evolution axis of the CAAS score.
- **Grouping scheme (US / GS1–GS4)** — US compares exact residues; GS schemes compare
  biochemical groups, so conservative substitutions within a group count as "the same".
- **Posterior** — the probability ASR assigns to a residue at a node; higher = more
  confident reconstruction.
