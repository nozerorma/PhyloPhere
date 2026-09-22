# References — Haemoglobin / high-altitude adaptation

## Primary data

- Zhu, X. et al. 2018. Divergent and parallel routes of biochemical adaptation
  in high-altitude passerine birds from the Qinghai-Tibet Plateau. *PNAS*
  115(8):1865–1870. doi:[10.1073/pnas.1720487115](https://doi.org/10.1073/pnas.1720487115)
  — sequence source (GenBank MG772099–MG772439), site-directed-mutagenesis
  functional validation (αA A34T, αA P119A), and the H/L altitude
  classification `../input/hb/README.md`'s `elev_mid` encoding is built from.

## Species-tree topology

- Johansson, U.S. et al. 2013. A complete multilocus species phylogeny of the
  tits and chickadees (Aves: Paridae). *Molecular Phylogenetics and Evolution*
  69(3):852–860. doi:[10.1016/j.ympev.2013.06.019](https://doi.org/10.1016/j.ympev.2013.06.019)
- Li, X. et al. 2016. Molecular phylogeny and diversification of the Old
  World Paridae. *Molecular Phylogenetics and Evolution* 104:14–20.
  doi:[10.1016/j.ympev.2016.07.021](https://doi.org/10.1016/j.ympev.2016.07.021)

## Not available, not used, not needed

Two Natarajan papers on independent high-altitude hemoglobin systems
(Andean waterfowl; bar-headed goose) were cited in an earlier draft of this
fixture's truth-set notes, but were never actually obtained, and a decision
was made not to fetch them. They would have supplied a second, independent
axis of validated positions to compare against (the way ConDor's reanalysis
does for PEPC) — HB currently has no such axis; every row in
`hb_altitude.sites.tsv` traces to Zhu et al. 2018 alone. Do not treat any
claim as if it were backed by these:

- Natarajan, C. et al. 2015. Convergent evolution of hemoglobin function in
  high-altitude Andean waterfowl involves limited parallelism at the
  molecular sequence level. *PLOS Genetics* 11(12):e1005681.
  doi:[10.1371/journal.pgen.1005681](https://doi.org/10.1371/journal.pgen.1005681)
- Natarajan, C. et al. 2016. Molecular basis of hemoglobin adaptation in
  the high-flying bar-headed goose. *PLOS Genetics* 12(12):e1006456.
  doi:[https://doi.org/10.1371/journal.pgen.1007331](https://doi.org/10.1371/journal.pgen.1007331
)

(Two *different* Natarajan papers, on deer mouse hemoglobin, were briefly
present under `references/hb/` under the mistaken assumption they were the
ones above — confirmed unrelated by reading them, and removed.)

## Partially incorporated

Zhu et al. 2018 Fig. S3 (AncParidae → *P. humilis* / *L. dichrous*
substitution lists) — the full site-by-site diff is transcribed in
`../../truthsets/tier1/hb_altitude.sites.tsv`'s own header comment (as
documentation, not as scored truth-set rows); see that file for exactly
which positions were checked and which weren't.
