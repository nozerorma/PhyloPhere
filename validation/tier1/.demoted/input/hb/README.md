# Tier 1 fixture — haemoglobin high-altitude adaptation (Sino-Himalayan tits)

The **continuous-trait** site-truth benchmark for Tier 1. Elevation in Paridae +
Aegithalidae of the Qinghai-Tibet Plateau; three globin genes (αA / αD / βA).
Fills `continuous · single · {direct, percentilized}` **and**
`categorical · single · direct` in one fixture.

Source: **Zhu et al. 2018**, *Divergent and parallel routes of biochemical
adaptation in high-altitude passerine birds from the Qinghai-Tibet Plateau*,
PNAS 115(8):1865-1870, doi:10.1073/pnas.1720487115 (full text: PMC5828625).

## Build

```bash
micromamba env create -f validation/tier1/input/fixture-tools.yml   # once
micromamba run -n phylophere-fixtures python3 \
    validation/tier1/input/hb/scripts/build.py
```

Fetches GenBank MG772099-MG772439 (gitignored, under `zhu2018/`), builds
per-gene AA alignments (MAFFT L-INS-i), infers trees (IQ-TREE), roots +
time-scales the species tree (`ape::chronos`), writes the pipeline-shaped
fixture: `align/{HBA,HBD,HBB}.fasta`, `tree.nwk` (ultrametric chronogram),
`tree_substitution.nwk` (raw ML phylogram), `gene_trees/{HBA,HBD,HBB}.nwk`
(unconstrained per-gene ML trees), `my_traits.tsv`, `phenotype.tsv`, synthetic
`ali_sp_names.txt` / `gene_ensembl.tsv` / `taxid.tsv`. Only `scripts/build.py`,
`scripts/build_cds.py`, this README and `hb_altitude.spec.json` are committed.

## Provenance

| component | source |
|-----------|--------|
| sequences | GenBank **MG772099-MG772439**, submitted by Zhu et al. for this paper (verified: matching authors/title/journal on the individual records). alphaA/alphaD/betaA globin, 3-16 population isolates per species per gene. One **per-column majority consensus** per species per gene; initiator Met stripped from alphaA/betaA (cleaved in the mature chain), retained for alphaD. |
| species set | the 16 species in the paper's Fig. S1 (14 Paridae + 2 Aegithalidae), confirmed by name against the SI PDF: *Periparus ater aemodius/pekinensis*, *P. rubidiventris*, *Pardaliparus venustulus*, *Poecile palustris/davidi/montanus*, *Lophophanes dichrous*, *Parus minor/monticulus/spilonotus/humilis*, *Cyanistes cyanus*, *Sylviparus modestus*, *Aegithalos bonvaloti/fulginosus*. This fixture uses current genus assignments (e.g. *Machlolophus spilonotus* for the paper's *Parus spilonotus*, *Pseudopodoces humilis* for the paper's *Parus humilis*) from the same GenBank records' own taxonomy. |
| tree topology | Fig. S1's own caption: "The tree topology is based on data from Johansson et al. (65) and Li et al. (66)" -- i.e. Johansson et al. 2013, MPE 69:852, and Li et al. 2016, MPE 104:14, as cited by the paper itself, not independently sourced by this fixture. Hard-coded in `build.py:_TOPOLOGY`, transcribed from Fig. S1 (`validation/tier1/references/hb/pnas.201720487si.pdf`, p.2): Aegithalidae outgroups all of Paridae; within Paridae, (Periparus+Pardaliparus) is sister to (Poecile+Lophophanes), and that pair is sister to a third clade of (Sylviparus,(Cyanistes,(Pseudopodoces,(Machlolophus,Parus)))), with Sylviparus the most basal member of that third clade, not of Paridae as a whole. |
| branch lengths | ML (LG+G) on the concatenated alphaA+alphaD+betaA alignment with the topology fixed (`iqtree -te`), then rooted on the Aegithalidae outgroup and time-scaled to an ultrametric chronogram (`ape::chronos`, penalised likelihood, lambda=1, correlated rates; root age = 1) -- this fixture's own processing, not from the paper. PhyloPhere contrast selection (modified Dunn + OU/BM PSS) assumes a time tree. The raw phylogram is kept as `tree_substitution.nwk`. See `build.py:_datetree`. |
| gene trees | unconstrained ML per gene (`gene_trees/{HBA,HBD,HBB}.nwk`), this fixture's own inference, meant to reproduce the paper's reported genealogical discordance among loci (Fig. S2) rather than to test against it directly. |

ConDor / PCOC / CAAStools have **not** been run on this dataset -- there is no
ready method-comparison column here (unlike PEPC).

## Trait

Zhu et al. do not use a fixed elevation cutoff; they compare high-altitude
species against their closest lowland relatives as sister pairs. This
fixture's `elev_mid` (per-species elevational-range midpoint) and the
`altitude` H/L split (cut at 2500 m) are both **this fixture's own encoding**,
needed to give PhyloPhere continuous and categorical trait columns -- neither
is the paper's own method.

Elevation ranges:
- **7 verbatim from the paper text**: *P. humilis* 3100-5500 m, *P. rubidiventris*
  2400-4300 m, *L. dichrous* 2300-4600 m, *P. ater aemodius* 2100-4600 m,
  *P. minor* 0-2000 m, *P. palustris* 0-2100 m, *P. ater pekinensis* 0-1800 m.
- **9 read off Fig. S1**'s elevation-range bar chart (+/-~150 m estimate), consistent
  with the plotted bars: *Aegithalos bonvaloti*, *A. fulginosus*, *Cyanistes cyanus*,
  *Poecile davidi*, *P. montanus*, *Pardaliparus venustulus*, *Parus monticolus*,
  *Parus spilonotus*, *Sylviparus modestus*.

`hb_altitude.spec.json` declares both traits; categorical `pairs` are this
fixture's own contrast pairing, not from the paper.

## Numbering

Mature globin chain, standard Hb convention. `align/HBA.fasta` and
`align/HBB.fasta` have the initiator Met stripped, so alignment column N ==
mature residue N (1:1, no gaps in the tit globin block -- verified in
`build.py:_assert_landmarks`). `align/HBD.fasta` keeps the retained alphaD Met,
same column==residue identity.

## Truth

`validation/truthsets/tier1/hb_altitude.sites.tsv`. Both sites are reported
as the paper's own site-directed-mutagenesis (SDM) results, tested on
recombinantly expressed hemoglobin:

- **alphaA A34T** -- parallel substitution in *Pseudopodoces humilis* (= *Parus
  humilis*) and *Lophophanes dichrous*; the paper's SDM on the AncParidae
  background found it increases Hb-O2 affinity. The paper states Ala is the
  ancestral state at this site (established via a CpG-dinucleotide decay
  argument across a broader passerine outgroup sample) and Thr34 is derived,
  arising independently in each of these two species.
- **alphaA P119A** -- *Aegithalos bonvaloti*; the paper's own recombinant-expression
  experiments on the *Aegithalos* high/low-altitude sister pair found it
  increases Hb-O2 affinity. The same substitution is reported elsewhere in the
  bar-headed goose (*Anser indicus*), outside this taxon set.

### AncParidae -> P. humilis / L. dichrous, full substitution lists (Fig. S3)

Diffing the paper's own typed AncParidae reference sequence (Fig. S3) against
this fixture's built `align/HBA.fasta` / `align/HBB.fasta` sequences for these
two species recovers exactly the alphaA34 A->T substitution reported in the main
text, with no extra or missing changes at that site -- also a check that this
fixture's per-species consensus reproduces the paper's own wild-type
sequences.

alphaA-globin (141 sites), AncParidae residue -> derived residue:

| position | AncParidae | *P. humilis* | *L. dichrous* |
|---|---|---|---|
| 8 | S | A | -- |
| 22 | E | -- | D |
| 34 | A | **T** | **T** |
| 44 | P | S | -- |
| 49 | N | -- | G |
| 57 | G | -- | A |
| 70 | A | V | -- |
| 72 | N | -- | H |
| 103 | Q | H | -- |
| 108 | V | A | A |
| 109 | V | I | -- |
| 113 | N | S | -- |

betaA-globin (146 sites):

| position | AncParidae | *P. humilis* | *L. dichrous* |
|---|---|---|---|
| 25 | G | A | -- |
| 43 | A | S | -- |
| 44 | S | N | -- |
| 51 | A | P | -- |

*L. dichrous* betaA-globin has no substitutions relative to AncParidae in this
fixture's sequence.

## Caveats

- Fig. S4 (the paper's own cross-passerine survey of alphaA34 substitutions)
  marks *Sylviparus modestus* as matching the ancestral state at alphaA34. All 8
  *Sylviparus modestus* isolates deposited under this paper's own GenBank
  accessions (MG772099-MG772439) instead translate to Thr34 -- the derived
  state -- via the same first-codon-position G->A change reported for
  *P. humilis*/*L. dichrous* (`ACC`/`ACC`/`ACT` at that codon; the
  third-position difference in *Sylviparus* is synonymous). This fixture's
  `align/HBA.fasta` reflects the GenBank sequence data (Thr34 for
  *Sylviparus modestus*), not Fig. S4's table.
- `gene_trees/{HBA,HBD,HBB}.nwk` are this fixture's own unconstrained ML gene
  trees, not the paper's. The GUI template's `gene_trees` parameter
  (`gui/templates/tier1_hb_altitude.json`) points at a single file,
  `gene_trees.nwk`, at the fixture root -- that file does not currently exist;
  the current build produces the three separate per-gene files under
  `gene_trees/` instead.

## GUI template wiring (`../../../../gui/templates/tier1_hb_altitude.json`)

- `taxid.tsv` / `gene_ensembl.tsv` in this directory are `build.py`'s synthetic
  placeholders (fake genomic coordinates, no real Ensembl gene ids for these
  three globin paralogs). The template points `gene_ensembl_file` at
  `gene_ensembl.tsv` directly; `tax_id_file` is left blank, so PhyloPhere
  generates it itself.
- `caas.pss_top_pct` is raised to **1** in the template (pipeline default:
  0.05) -- on this 16-tip tree the default rank-percentile PSS cut leaves too
  few candidate hi>lo pairs eligible for `elev_mid`'s continuous contrast
  selection.
- `rer_minsp` is lowered to **4** in the template (pipeline default: 15) -- the
  default is close to this dataset's total tip count (16) and would filter
  too aggressively once RER prunes to species present in the trait file.

## Codon-level CDS and positive selection (FUBAR)

`build_cds.py` (run after `build.py`) extracts the underlying nucleotide CDS
from `zhu2018/hb_genbank.gb` (already-spliced mRNA records with a `/translation`
qualifier) and back-translates it into `align/{HBA,HBD,HBB}.fasta`'s existing
gap pattern, producing `align_cds/{HBA,HBD,HBB}.fasta` (codon nucleotide
alignments; gitignored, rebuild with `python3 scripts/build_cds.py`). Each
species' nucleotide consensus is checked to translate to exactly the
amino-acid sequence already in `align/`, so the two stay biologically
consistent.

`align_cds/` is `ortholog_characterizator`'s `--cds_dir` input
(`validation/tier1/run_ortholog_characterizator.sh hb`: quality pass-through ->
BMGE codon trimming/translation -> FUBAR/MEME), rooted the same way as the
protein-only build above. Species tree for positive selection:
`tree_substitution.nwk` (substitution units, not the time-scaled `tree.nwk`).
Output would land under `oc_run/`, which does not currently exist in this
fixture (the pipeline has not been run against it).
