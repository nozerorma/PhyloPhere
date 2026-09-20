# Report — Haemoglobin / high-altitude adaptation (Sino-Himalayan tits)

## a) Data source

GenBank MG772099–MG772439 (341 records: αA/αD/βA globin, 3–16 population
isolates per species/gene), per-species per-gene consensus, 16 focal taxa (14
Paridae + 2 Aegithalidae, Zhu et al. 2018 Fig S1). Source: Zhu et al. 2018,
*PNAS* 115(8):1865–1870, doi:10.1073/pnas.1720487115. Traits: `elev_mid`
(continuous, elevational-range midpoint — our encoding, Zhu analyse H/L only)
and `altitude` (categorical, 2500 m cut). Truth set: 5 curated sites
(`../../truthsets/tier1/hb_altitude.sites.tsv`) — 2 site-directed-mutagenesis
hits (Storz lab) + 3 parallel-candidate sites, ~5 independent high-altitude
origins. Comparators: Natarajan 2015/2016 (Andean birds, one site overlap);
no PCOC/ConDor run exists on this dataset. Full provenance, including two
documented "features not bugs" complications (ancestral-state retention in
*Sylviparus*, genealogical discordance in the gene trees):
[`../input/hb/README.md`](../input/hb/README.md).

## b) Scope

Continuous- and categorical-trait recovery check on a small (16-tip),
contrast-poor dataset — framed, like echolocation-type fixtures, as a hard case
and a targeted CT_DISAMBIGUATION test, **not** a high-power precision benchmark.
Fills both the `continuous · single · direct` and `categorical · single ·
direct` trait-matrix cells from one fixture.

## c) Procedure

**Legacy run** (`../../runs/tier1/hb_altitude/{elev_mid,altitude}_mdf0.5/`,
pre-GUI harness scripts, since retired): CAAS/CAAP discovery+resample →
CT_DISAMBIGUATION → FADE → SCORING. RER/accumulation/enrichment/VEP off.
`min_divergent_fraction=0.5`, `max_fop=100`, seed 1998, local; continuous run
had the PSS top-quantile gate **disabled outright** (pipeline default 0.05
left ~6/256 candidate pairs eligible on this tree).

**Current template** (`../../../gui/templates/tier1_hb_altitude.json`): same
CAAS/CT_DISAMBIGUATION/FADE chain plus **RER** (3 genes — enough for a
minimally meaningful cross-gene comparison, unlike PEPC), accumulation/VEP/
enrichment off, `local_lowspec` resource preset. Deviates from pipeline
defaults in two places, both tree-size-driven and documented in
[`../input/hb/README.md`](../input/hb/README.md) rather than tuned for
recall: `pss_top_pct=0.30` (not disabled — a milder version of the legacy
workaround) and `rer_minsp=4` (default 15 is too close to this dataset's 16
total tips). Not yet re-run under this template; §d is the legacy
CAAS+FADE-only result.

## d) Results and conclusions

**3/5 sites recovered, in both trait encodings** (categorical and continuous
agree on which sites recover):

| site | tier | rank-pct (categorical / continuous) | ASR score |
|------|------|--------------------------------------|-----------|
| αA34 A→T | mutagenesis | 0.67 / 0.78 | 0.74 |
| αA109 V→I | parallel-candidate | 0.56 / 0.44 | 0.37 |
| αD18 G→S | parallel-candidate | **1.00 / 1.00** | 0.74 |
| αA119 P→A | mutagenesis | missed | — |
| βA83 N→G | parallel-candidate | missed | — |

**Both misses are single-species autapomorphies in this fixture** (αA119A
private to *Aegithalos bonvaloti*, βA83G to *Poecile davidi*, already flagged
as "variable" in the truth set) — no convergence method, CAAS or otherwise,
can call a one-lineage change convergent. Excluding those two structurally
undetectable sites, **the pipeline recovers 3/3 of the genuinely convergent
sites**, in both encodings — a materially different and more informative
number than the raw 3/5. αD18 and αA109 only appear because a
minority-convergence pre-filter was removed pipeline-wide (documented in the
legacy walkthrough, not specific to this dataset); αA34 in the categorical
encoding only because `max_fop` was raised to the pipeline default (100, not
the harness's old hard-coded 15).

Flagging explicitly since it bears on how much weight this result should
carry: 5 truth sites on a 16-tip tree is a small-sample result. 3/3 (minus the
two structurally-unrecoverable autapomorphies) is a real, mechanistically
explained recovery — not a heuristic pattern-match — but it is not a
statistically powered precision/recall estimate, and shouldn't be reported as
one.
