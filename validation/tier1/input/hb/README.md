# Tier 1 fixture — haemoglobin high-altitude adaptation (Sino-Himalayan tits)

The **continuous-trait** site-truth benchmark for Tier 1. Elevation in Paridae +
Aegithalidae of the Qinghai-Tibet Plateau; three globin genes (αA / αD / βA).
Fills `continuous · single · {direct, percentilized}` **and**
`categorical · single · direct` in one fixture.

Source: **Zhu et al. 2018**, *Divergent and parallel routes of biochemical
adaptation in high-altitude passerine birds from the Qinghai-Tibet Plateau*,
PNAS 115(8):1865–1870, doi:10.1073/pnas.1720487115 (full text: PMC5828625).

## Build

```bash
micromamba env create -f validation/tier1/input/fixture-tools.yml   # once
micromamba run -n phylophere-fixtures python3 \
    validation/tier1/input/hb/build.py
```

Fetches GenBank MG772099–MG772439 (gitignored), builds per-gene AA alignments
(MAFFT L-INS-i), infers trees (IQ-TREE), roots + time-scales the species tree
(`ape::chronos`), writes the pipeline-shaped fixture:
`align/{HBA,HBD,HBB}.fasta`, `tree.nwk` (ultrametric chronogram),
`tree_substitution.nwk` (raw ML phylogram), `my_traits.tsv`, `phenotype.tsv`, synthetic `ali_sp_names.txt` /
`gene_ensembl.tsv` / `taxid.tsv`. Only `build.py`, this README and
`hb_altitude.spec.json` are committed.

## Provenance

| component | source |
|-----------|--------|
| sequences | GenBank **MG772099–MG772439**, submitted by Zhu et al. for this paper (verified: matching authors/title/journal on the individual records). αA/αD/βA globin, 3–16 population isolates per species per gene. One **per-column majority consensus** per species per gene; initiator Met stripped from αA/βA (cleaved in the mature chain), retained for αD. |
| species set | the 16 species in the paper's Fig. S1 (14 Paridae + 2 Aegithalidae), confirmed by name against the SI PDF: *Periparus ater aemodius/pekinensis*, *P. rubidiventris*, *Pardaliparus venustulus*, *Poecile palustris/davidi/montanus*, *Lophophanes dichrous*, *Parus minor/monticulus/spilonotus/humilis*, *Cyanistes cyanus*, *Sylviparus modestus*, *Aegithalos bonvaloti/fulginosus*. This fixture uses current genus assignments (e.g. *Machlolophus spilonotus* for the paper's *Parus spilonotus*, *Pseudopodoces humilis* for the paper's *Parus humilis*) from the same GenBank records' own taxonomy. |
| tree topology | Fig. S1's own caption: "The tree topology is based on data from Johansson et al. (65) and Li et al. (66)" — i.e. Johansson et al. 2013, MPE 69:852, and Li et al. 2016, MPE 104:14, as cited by the paper itself, not independently sourced by this fixture. Hard-coded in `build.py:_TOPOLOGY`; not verified node-by-node against Fig. S1's actual tree image beyond the well-established placements (*Sylviparus* as the most basal genus among those sampled here; *Pseudopodoces humilis* nested within Paridae, not a corvid). |
| branch lengths | ML (LG+G) on the concatenated αA+αD+βA alignment with the topology fixed (`iqtree -te`), then rooted on the Aegithalidae outgroup and time-scaled to an ultrametric chronogram (`ape::chronos`, penalised likelihood, λ=1, correlated rates; root age = 1) — this fixture's own processing, not from the paper. PhyloPhere contrast selection (modified Dunn + OU/BM PSS) assumes a time tree. The raw phylogram is kept as `tree_substitution.nwk`. See `build.py:_datetree`. |
| gene trees | unconstrained ML per gene, this fixture's own inference (not the paper's). The paper reports genealogical discordance among its own gene trees ("high levels of genealogical discordance among loci, likely reflecting a history of incomplete lineage sorting and/or introgressive hybridization"); whether this fixture's independently-rebuilt gene trees show the same discordance has not been checked. |

ConDor / PCOC / CAAStools have **not** been run on this dataset — there is no
ready method-comparison column here (unlike PEPC).

## Trait

Zhu et al. do not use a fixed elevation cutoff; they compare high-altitude
species against their closest lowland relatives as sister pairs. This
fixture's `elev_mid` (per-species elevational-range midpoint) and the
`altitude` H/L split (cut at 2500 m) are both **this fixture's own encoding**,
needed to give PhyloPhere continuous and categorical trait columns — neither
is the paper's own method.

Elevation ranges:
- **7 verbatim from the paper text**: *P. humilis* 3100–5500 m, *P. rubidiventris*
  2400–4300 m, *L. dichrous* 2300–4600 m, *P. ater aemodius* 2100–4600 m,
  *P. minor* 0–2000 m, *P. palustris* 0–2100 m, *P. ater pekinensis* 0–1800 m.
- **9 read off Fig. S1**'s elevation-range bar chart (±~150 m estimate). Checked
  against the actual Fig. S1 image now in hand (`validation/papers/pnas.201720487si.pdf`)
  and found consistent with the plotted bars for all 9: *Aegithalos bonvaloti*,
  *A. fulginosus*, *Cyanistes cyanus*, *Poecile davidi*, *P. montanus*,
  *Pardaliparus venustulus*, *Parus monticulus*, *Parus spilonotus*,
  *Sylviparus modestus*.

`hb_altitude.spec.json` declares both traits; categorical `pairs` are this
fixture's own contrast pairing, not from the paper.

## Numbering

Mature globin chain, standard Hb convention. `align/HBA.fasta` and
`align/HBB.fasta` have the initiator Met stripped, so alignment column N ==
mature residue N (1:1, no gaps in the tit globin block — verified in
`build.py:_assert_landmarks`). `align/HBD.fasta` keeps the retained αD Met,
same column==residue identity.

## Truth

`validation/truthsets/tier1/hb_altitude.sites.tsv`. Both sites are reported
as the paper's own site-directed-mutagenesis (SDM) results, tested on
recombinantly expressed hemoglobin:

- **αA A34T** — parallel substitution in *Pseudopodoces humilis* (= *Parus
  humilis*) and *Lophophanes dichrous*; the paper's SDM on the AncParidae
  background found it increases Hb–O₂ affinity. The paper states Ala is the
  ancestral state at this site (established via a CpG-dinucleotide decay
  argument across a broader passerine outgroup sample) and Thr34 is derived,
  arising independently in each of these two species.
- **αA P119A** — *Aegithalos bonvaloti*; the paper's own recombinant-expression
  experiments on the *Aegithalos* high/low-altitude sister pair found it
  increases Hb–O₂ affinity. The same substitution is reported elsewhere in the
  bar-headed goose (*Anser indicus*), outside this taxon set.

### AncParidae → P. humilis / L. dichrous, full substitution lists (Fig. S3)

Computed by diffing the paper's own typed AncParidae reference sequence
(Fig. S3) against this fixture's built `align/HBA.fasta` / `align/HBB.fasta`
sequences for these two species — which is also a direct check that this
fixture's per-species consensus reproduces the paper's own wild-type
sequences: the diff recovers exactly the αA34 A→T substitution reported in
the main text, with no extra or missing changes at that site.

αA-globin (141 sites), AncParidae residue → derived residue:

| position | AncParidae | *P. humilis* | *L. dichrous* |
|---|---|---|---|
| 8 | S | A | — |
| 22 | E | — | D |
| 34 | A | **T** | **T** |
| 44 | P | S | — |
| 49 | N | — | G |
| 57 | G | — | A |
| 70 | A | V | — |
| 72 | N | — | H |
| 103 | Q | H | — |
| 108 | V | A | A |
| 109 | V | I | — |
| 113 | N | S | — |

βA-globin (146 sites):

| position | AncParidae | *P. humilis* | *L. dichrous* |
|---|---|---|---|
| 25 | G | A | — |
| 43 | A | S | — |
| 44 | S | N | — |
| 51 | A | P | — |

*L. dichrous* βA-globin has no substitutions relative to AncParidae in this
fixture's sequence.

## Discrepancy found, not resolved

Checking this fixture's own alignment against Fig. S4 (the paper's own
cross-passerine survey of αA34 substitutions) surfaced an internal
inconsistency in the paper's supplementary material that this README isn't
resolving: Fig. S4's summary table marks *Sylviparus modestus* as matching
the ancestral state at αA34, but all 8 *Sylviparus modestus* isolates
deposited under this paper's own GenBank accessions (MG772099–MG772439)
translate to Thr34 — the derived state — via the same first-codon-position
G→A change reported for *P. humilis*/*L. dichrous* (`ACC`/`ACC`/`ACT`
respectively at that codon; the third-position difference in *Sylviparus* is
synonymous). This fixture's `align/HBA.fasta` reflects the GenBank sequence
data (Thr34 for *Sylviparus modestus*), not Fig. S4's table.

## TODO

- Optional: fold in the Andean high-altitude birds (Natarajan 2015 PLOS
  Genetics 11:e1005681; Natarajan 2016 PLOS Genetics 12:e1006456) for
  additional independent origins — separate fixture-scope decision, not
  started.

## GUI template wiring (`../../../../gui/templates/tier1_hb_altitude.json`)

- `taxid.tsv` / `gene_ensembl.tsv` in this directory are `build.py`'s synthetic
  placeholders (fake genomic coordinates, no real Ensembl gene ids for these
  three globin paralogs). The template leaves `tax_id_file` /
  `gene_ensembl_file` blank, so PhyloPhere generates them itself.
- `gene_trees.nwk` (gitignored) concatenates `gene_trees_oc/{HBA,HBD,HBB}.nwk`
  into one 3-line file for RER's `--gene_trees`. These are IQ-TREE trees from
  `ortholog_characterizator` (`PHYLOGENY` workflow, ModelFinder `MFP`
  restricted to the LG family + 1000 UFBoot; best-fit `LG+I` for HBA/HBB,
  `LG+R3` for HBD), run against `align/{HBA,HBD,HBB}.fasta` directly
  (`--prot_dir`, quality/translation stages off) — the sole source for this
  file; `build.py` no longer emits its own separate per-gene trees.
- `caas.pss_top_pct` is raised to **0.30** (pipeline default: 0.05) in the
  template — on this 16-tip tree the default rank-percentile PSS cut leaves
  only ~6 of 256 candidate hi>lo pairs eligible for `elev_mid`'s continuous
  contrast selection, too little headroom for the FOP harvest.
- `rer_minsp` is lowered to **4** (pipeline default: 15) — the default is close
  to this dataset's total tip count (16) and would filter too aggressively once
  RER prunes to species present in the trait file.

## Codon-level CDS and positive selection (FUBAR)

`build_cds.py` (run after `build.py`) extracts the underlying nucleotide CDS
from `hb_genbank.gb` (already-spliced mRNA records with a `/translation`
qualifier) and back-translates it into `align/{HBA,HBD,HBB}.fasta`'s existing
gap pattern, producing `align_cds/{HBA,HBD,HBB}.fasta` (codon nucleotide
alignments; gitignored, rebuild with `python3 build_cds.py`). Each species'
nucleotide consensus is checked to translate to exactly the amino-acid
sequence already in `align/`, so the two stay biologically consistent.

`align_cds/` is `ortholog_characterizator`'s `--cds_dir` input for a full run
(`quality` pass-through → BMGE codon trimming/translation → IQ-TREE gene
trees → HyPhy FUBAR/MEME), rooted the same way as the protein-only
`PHYLOGENY`-only run above. Species tree for positive selection:
`tree_substitution.nwk` (substitution units, not the time-scaled `tree.nwk`).

Results in `positive_selection/` (tracked; `fubar_sites.tsv`,
`fubar_results.tsv`, `meme_results.tsv`): per-site dN/dS posterior
probabilities and per-gene pervasive-selection counts for HBA, HBD, HBB —
this fixture's own analysis, not from the paper. Neither truth site (αA34,
αA119) is a FUBAR hit, which is expected: FUBAR tests for *pervasive*
selection across the whole tree, while both sites are lineage-specific
substitutions in 1–2 tips.
