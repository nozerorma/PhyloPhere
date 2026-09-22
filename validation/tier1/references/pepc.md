# References — PEPC / C4 photosynthesis

## Primary data

- Besnard, G. et al. 2009. Phylogenomics of C4 photosynthesis in sedges
  (Cyperaceae): multiple appearances and genetic convergence. *Molecular
  Biology and Evolution* 26(8):1909–1919.
  doi:[10.1093/molbev/msp103](https://doi.org/10.1093/molbev/msp103)

## Method-comparison targets (same dataset)

- Rey, C. et al. 2018. Detecting adaptive convergent amino acid evolution.
  *Molecular Biology and Evolution* 35(9):2296–2306 (PCOC method).
  doi:[10.1093/molbev/msy114](https://doi.org/10.1093/molbev/msy114)
- Morel, B. et al. 2024. ConDor: a phylogenetic approach for detecting
  convergent evolution. *Genome Biology and Evolution* 16(6):evae040.
  doi:[10.1093/gbe/evae040](https://doi.org/10.1093/gbe/evae040)

## Redistribution sources (as fetched by `../input/pepc/build.py`)

- `github.com/CarineRey/pcoc` — `data/det/` (PCOC test data, original
  Besnard 2009 sequences + Christin 2007 supplement)
- `github.com/evolbioinfo/condor` — `test_data/` (ConDor's own copy, incl.
  the genotypic C4 annotation `besnard2009_convergent_species.txt`)

## Cited but not locally available

- Christin, P.-A. et al. 2007. Evolutionary switch and genetic convergence on
  rbcS at the C3/C4 boundary. *Genetics* 177(3):1791–1808. Not present under
  `references/pepc/` and not read directly for this fixture; every claim
  attributed to it in `../../truthsets/tier1/pepc_c4.sites.tsv` traces
  through Besnard 2009's own Table 2, which cites it for the grass-side data
  — see that TSV's header for exactly what is and isn't independently
  verified.
- Bläsing OE, Westhoff P, Svensson P. 2000. J Biol Chem 275:27917–27923 (the
  position-780 kinetics paper) and Svensson P et al. 2003 (cited by Morel et
  al. 2024 for the position-665 functional claim). Neither is present under
  `references/pepc/`; both claims are corroborated only through their being
  cited by a paper that *is* available (Besnard 2009 and Morel et al. 2024
  respectively) — see `../../truthsets/tier1/pepc_c4.sites.tsv`.
