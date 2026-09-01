# Tier 1 fixture — PEPC C3/C4 in sedges

Single-gene, categorical-trait site-truth benchmark. C3 vs C4 photosynthesis in
Cyperaceae, PEPC (`ppc`) amino-acid alignment. ≥5 independent C4 origins → the
multi-contrast structure CAAS needs.

## Build

```bash
python3 validation/fixtures/tier1/pepc/build.py
```

Fetches the raw inputs (gitignored) and writes the pipeline-shaped fixture:
`align/PEPC.fasta`, `tree.nwk`, `my_traits.tsv`, `phenotype.tsv`, plus synthetic
`ali_sp_names.txt` / `gene_ensembl.tsv` / `taxid.tsv`. Only `build.py` and this
README are committed.

## Provenance

| file | source |
|------|--------|
| `cyp_coding.aa.coor_mays.fa` | `github.com/evolbioinfo/condor/test_data/` — 79 Cyperaceae PEPC AA sequences, projected onto maize PEPC1 coordinates. Originally Besnard et al. 2009 (MBE 26:1909), provided by the authors, redistributed via PCOC (`CarineRey/pcoc/data/det/`) and ConDor. |
| `cyp_coding.phy_phyml_tree.txt` | same — PhyML tree, aLRT support on internal nodes. |
| `besnard2009_convergent_species.txt` | ConDor test data — 23 tips with the **"genotypic" C4 annotation** (presence of the A780S determinant). |
| `outgroup.txt` | ConDor test data — Chrysithrix, root; dropped from fg/bg. |

Same dataset used by PCOC (Rey et al. 2018, MBE 35:2296) and ConDor
(Morel et al. 2024, GBE, doi:10.1093/gbe/evae040) — their site calls are the
method-comparison targets.

## Numbering

Maize PEPC1 = UniProt **P04711**, 970 aa. The alignment is exactly 970 columns
on maize coordinates, so **truth-set position N = alignment column N** (1:1). No
reference row needed. Verified: column 780 is S in 22/23 C4 tips and A in 54/55
C3 tips (the A780S determinant); columns 665, 572, 540, 749 likewise match the
Besnard / ConDor substitutions.

## Trait

`my_traits.tsv`: `c4` = 1 for the 23 genotypic-C4 tips, 0 for the 55 C3 tips.
Outgroup excluded. This is the **genotypic** annotation — it matches the
convergent clades by construction (favourable to detectors). A **phenotypic**
annotation (label by the plant's metabolism, drop the ~7 C3+C4 facultative
species) is the harder, more realistic variant; ConDor showed PCOC fails on it.
Worth running both as a sensitivity analysis.

## Caveats

- AA alignment, not codon — fine for CAAS/CAAP + FADE (both AA in PhyloPhere);
  no dN/dS.
- Besnard 2009 lists **16** positive-selection codons; the truth set transcribes
  the confident subset (7 tiered + 4 weak). The rest need the Besnard 2009
  Table 1 / Christin 2007 Table 1 PDFs.
- Only positions **780** (kinetics) and **665** (folding/activity) are
  functionally demonstrated. The rest are selection inferences — recall on them
  is a soft number, report with a CI.
