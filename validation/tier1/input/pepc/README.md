# Tier 1 fixture — PEPC C3/C4 in sedges

Single-gene, categorical-trait site-truth benchmark. C3 vs C4 photosynthesis in
Cyperaceae, PEPC (`ppc`) amino-acid alignment.

Source: **Besnard, G. et al. 2009**, *Phylogenomics of C4 photosynthesis in
sedges (Cyperaceae): multiple appearances and genetic convergence*, MBE
26(8):1909–1919, doi:10.1093/molbev/msp103.

## Build

```bash
python3 validation/tier1/input/pepc/build.py
```

Fetches the raw inputs (gitignored) and writes the pipeline-shaped fixture:
`align/PEPC.fasta`, `tree.nwk` (PhyML tree rooted on *Chrysithrix* and
time-scaled to an ultrametric chronogram via `ape::chronos`), `tree_substitution.nwk`
(raw phylogram), `my_traits.tsv`, `phenotype.tsv`, plus synthetic
`ali_sp_names.txt` / `gene_ensembl.tsv` / `taxid.tsv`. Only `build.py` and this
README are committed.

## Provenance

| file | source |
|------|--------|
| `cyp_coding.aa.coor_mays.fa` | `github.com/evolbioinfo/condor/test_data/` — 79 Cyperaceae PEPC amino-acid sequences, projected onto maize PEPC1 coordinates. Besnard et al. 2009, redistributed via PCOC (`CarineRey/pcoc/data/det/`) and ConDor. |
| `cyp_coding.phy_phyml_tree.txt` | same — PhyML tree, aLRT support on internal nodes. |
| `besnard2009_convergent_species.txt` | ConDor test data — 23 tips with the **"genotypic" C4 annotation** (presence of the A780S determinant). |
| `outgroup.txt` | ConDor test data — Chrysithrix, root; dropped from fg/bg. |

This fixture uses the amino-acid alignment as redistributed by PCOC (Rey et
al. 2018, MBE 35:2296) and ConDor (Morel et al. 2024, GBE,
doi:10.1093/gbe/evae040) — the same alignment those two tools' own published
analyses of this dataset use. Their site calls are the method-comparison
targets for this fixture's truth set.

## Numbering

Maize PEPC1 = UniProt **P04711**, 970 aa. The alignment is exactly 970 columns
on maize coordinates, so truth-set position N = alignment column N (1:1). No
reference row needed. Column 780 is S in 22/23 C4 tips and A in 54/55 C3 tips
(the A780S determinant); columns 665, 572, 540, 749 likewise match the
Besnard / ConDor substitutions.

## Trait

`my_traits.tsv`: `c4` = 1 for the 23 genotypic-C4 tips, 0 for the 55 C3 tips.
Outgroup excluded. This is the **genotypic** annotation (the tips Besnard et
al. and ConDor list as carrying the A780S determinant), as distinct from a
**phenotypic** annotation (by the plant's actual C3/C4 metabolism).

## Caveats

- Amino-acid alignment, not codon. Besnard et al. 2009 deposited nucleotide
  sequences for part of this dataset at EMBL (accessions FM208000–FM208067),
  but PCOC and ConDor both redistribute — and run their own published
  analyses on — the amino-acid alignment; this fixture does the same.
- Besnard 2009 lists 16 positive-selection codons; the truth set transcribes
  the confident subset (7 tiered + 4 weak). The rest are in the Besnard 2009
  Table 1 / Christin et al. 2007 (Genetics 177(3):1791–1808) Table 1.
- Only positions 780 (kinetics) and 665 (folding/activity) are functionally
  demonstrated (site-directed mutagenesis). The remaining truth-set positions
  are selection inferences from the source papers, not mutagenesis-confirmed.

## GUI template wiring (`../../../../gui/templates/tier1_pepc.json`)

- `taxid.tsv` / `gene_ensembl.tsv` in this directory are `build.py`'s synthetic
  placeholders (fake genomic coordinates; the PEPC `human_protein_id` is
  maize's own UniProt id). The template leaves `tax_id_file` /
  `gene_ensembl_file` blank, so PhyloPhere generates them itself.
- `gene_trees.nwk` (gitignored) is IQ-TREE output from `ortholog_characterizator`
  (`PHYLOGENY` workflow, ModelFinder `MFP` restricted to the LG family, 1000
  UFBoot), run against `align/PEPC.fasta` directly (`--prot_dir`,
  quality/translation stages off, consistent with this fixture being
  amino-acid only).
