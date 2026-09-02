# Tier 1 fixture — haemoglobin high-altitude adaptation (Sino-Himalayan tits)

The **continuous-trait** site-truth benchmark for Tier 1. Elevation in Paridae +
Aegithalidae of the Qinghai-Tibet Plateau; three globin genes (αA / αD / βA).
Fills `continuous · single · {direct, percentilized}` **and**
`categorical · single · direct` in one fixture.

Source: **Zhu et al. 2018**, *Divergent and parallel routes of biochemical
adaptation in high-altitude passerine birds from the Qinghai-Tibet Plateau*,
PNAS 115(8):1865–1870, doi:10.1073/pnas.1720487115.

## Build

```bash
micromamba env create -f validation/fixtures/fixture-tools.yml   # once
micromamba run -n phylophere-fixtures python3 \
    validation/fixtures/tier1/hb_altitude/build.py
```

Fetches GenBank MG772099–MG772439 (gitignored), builds per-gene AA alignments
(MAFFT L-INS-i), infers trees (IQ-TREE), roots + time-scales the species tree
(`ape::chronos`), writes the pipeline-shaped fixture:
`align/{HBA,HBD,HBB}.fasta`, `tree.nwk` (ultrametric chronogram),
`tree_substitution.nwk` (raw ML phylogram), `gene_trees/{HBA,HBD,HBB}.nwk`,
`my_traits.tsv`, `phenotype.tsv`, synthetic `ali_sp_names.txt` /
`gene_ensembl.tsv` / `taxid.tsv`. Only `build.py`, this README and
`hb_altitude.spec.json` are committed.

## Provenance

| component | source |
|-----------|--------|
| sequences | GenBank **MG772099–MG772439** (341 records: αA/αD/βA globin, 3–16 population isolates per species per gene). One **per-column majority consensus** per species per gene; initiator Met stripped from αA/βA (cleaved in the mature chain), retained for αD. |
| species set | the **16 focal taxa of Fig S1** (14 Paridae + 2 Aegithalidae). GenBank also carries extra *Aegithalos* and a pooled *Periparus ater* entry — not used. |
| tree topology | **Fig S1** (Johansson et al. 2013, MPE 69:852; Li et al. 2016, MPE 104:14). *Sylviparus* basal in Paridae; *Pseudopodoces* nested in the *Parus* group; the two Aegithalidae are the outgroup. Hard-coded in `build.py:_TOPOLOGY`. |
| branch lengths | ML (LG+G) on the concatenated αA+αD+βA alignment with the topology fixed (`iqtree -te`), then rooted on the Aegithalidae outgroup and time-scaled to an **ultrametric chronogram** (`ape::chronos`, penalised likelihood, λ=1, correlated rates; root age = 1). PhyloPhere contrast selection (modified Dunn + OU/BM PSS) assumes a time tree — on the raw phylogram *Pseudopodoces humilis*' fast-evolving long terminal branch inflates its contrast-pair diameter and the *Pseudopodoces*~*Parus_minor* truth contrast (αA A34T) scores a sub-1 Dunn. The raw phylogram is kept as `tree_substitution.nwk`. See `build.py:_datetree`. |
| gene trees | unconstrained ML per gene — reproduce the Fig S2 genealogical discordance (see below). |

ConDor / PCOC / CAAStools have **not** been run on this dataset — there is no
ready method-comparison column here (unlike PEPC). The comparators are the
Zhu 2018 analysis itself and Natarajan 2015/2016.

## Trait

Zhu et al. analyse altitude **categorically** (high / low). The fixture's
**continuous** `elev_mid` trait — per-species elevational-range midpoint — is
**our encoding**, built from:

- **7 ranges quoted verbatim in the paper text**: *P. humilis* 3100–5500,
  *P. rubidiventris* 2400–4300, *L. dichrous* 2300–4600, *P. ater aemodius*
  2100–4600, *P. minor* 0–2000, *P. palustris* 0–2100, *P. ater pekinensis*
  0–1800.
- **9 ranges read off Fig S1** (±~150 m, checked by MR 2026-08-31).

`my_traits.tsv` columns: `elev_mid` (m, continuous), `altitude` (H/L, cut at
**2500 m**), `family`. The 2500 m cut reproduces Zhu's H/L assignment for the
species they tested; *Poecile davidi* / *P. montanus* fall on the high side by
midpoint but were not functionally tested — treat their H label as soft.

`hb_altitude.spec.json` declares both traits. Categorical `pairs` are 5
phylogenetically independent high/low contrasts (each with an unambiguous H and L
member), 4 of which are Zhu's own key comparisons.

## Numbering

Mature globin chain, standard Hb convention. `align/HBA.fasta` and
`align/HBB.fasta` have the initiator Met stripped, so **alignment column N ==
mature residue N** (1:1, no gaps in the tit globin block — verified in
`build.py:_assert_landmarks`). `align/HBD.fasta` keeps the retained αD Met, same
column==residue identity. No reference row needed.

Landmarks checked on every build:
- **αA 34**: Thr in *P. humilis*, *L. dichrous* (+ *Sylviparus*, see below); Ala
  elsewhere.
- **αA 119**: Ala in *A. bonvaloti* only; Pro elsewhere.

## Truth

`validation/truthsets/tier1/hb_altitude.sites.tsv`

- **αA A34T** (mutagenesis) — parallel in *Pseudopodoces humilis* (= *Parus
  humilis*) and *Lophophanes dichrous*; SDM on the AncParidae background
  confirmed the Hb-O₂-affinity increase; the αA **gene tree** shows the two Thr34
  states are independent (CpG→CpA transitions).
- **αA P119A** (mutagenesis) — *Aegithalos bonvaloti*; SDM-confirmed
  affinity-enhancing (parallel with the bar-headed goose, outside this set).
- **parallel-candidate** (αA 109, αD 18, βA 83) — high-altitude-private in this
  fixture, ≥2 independent lineages, not individually functionally tested;
  partly alignment-derived, so **not** a headline precision number. βA83
  coincides with a Natarajan 2015 Andean site.

## Known complications (these are features, not bugs)

- **αA34 in *Sylviparus modestus*** (low altitude, basal Paridae): also Thr34.
  This is **retention of the ancestral avian state** — chicken αA is Thr34; the
  Paridae-crown ancestor (AncParidae) is the *derived* Ala34; *P. humilis* /
  *L. dichrous* reverted to Thr34 independently. So a low-altitude tip carries
  the "high-altitude" residue. A clean CAAS call at αA34 needs the contrast
  selector **not** to pick *Sylviparus* as a background partner, and ASR /
  CT_DISAMBIGUATION to score its Thr34 as retention, not convergence. This is
  precisely what those modules exist for — recovery here is a real test.
- **αA34 in *Machlolophus spilonotus*** (low altitude): Ile34 — a third state at
  the site.
- **Genealogical discordance.** The αA/αD/βA gene trees (`gene_trees/`,
  reproducing Fig S2) are strongly discordant with the species tree (ILS +
  introgression). PhyloPhere works off the **species tree**; species-tree
  "apparent parallelism" that the gene trees contradict is a false-positive
  channel the fixture can probe.
- **Small truth set + contrast-poor.** ~5 independent high-altitude origins, one
  gold parallel site (αA34) within this taxon set. Frame like echolocation: a
  hard case and a targeted CT_DISAMBIGUATION test, **not** a high-power precision
  benchmark.

## TODO

- Transcribe the AncParidae → *P. humilis* (12 αA sites) / *L. dichrous* (6 αA
  sites) substitution lists from Fig S3 for a "divergent" truth tier.
- Optional: fold in the Andean high-altitude birds (Natarajan 2015/2016) for more
  independent origins + the βA 86/94/116 truth sites (invariant here). Separate
  fixture-scope decision.
