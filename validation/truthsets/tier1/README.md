# Tier 1 truth sets — provenance

Curated known-positive lists. Small, version-controlled, hand-checked. Every row
traces to a citation key resolved here.

Coordinates in `*.sites.tsv` are in the `ref_seq` sequence's ungapped 1-based
numbering. `harness.truthset.map_site_to_alignment` converts to alignment columns
using the gapped `ref_seq` row from the fixture alignment.

## Status

| file | level | rows | completeness |
|------|-------|------|--------------|
| `rh1_spectral.sites.tsv` | site | 12 | core done; full ~30-site set is TODO (Hauser 2017 SI) |
| `pepc_c4.sites.tsv` | site | 3 | only pos 780 validated; ~20-codon set is TODO (Christin 2007 / Besnard 2009) |
| `echolocation.genes.tsv` | gene | 6 | core done; also serves as null negative control |

Recall should not be quoted as a headline figure until the TODO rows are
transcribed from the source tables. Until then, use the core rows for
"are the validated positives recovered and highly ranked" only.

## Citation keys

- **Yokoyama2008** — Yokoyama S. "Evolution of dim-light and color vision pigments." Annu Rev Genomics Hum Genet 9:259–282.
- **Yokoyama2014** — Yokoyama S et al. "Epistatic adaptive evolution of human color vision." PLoS Genet / related RH1 mutagenesis series.
- **Hauser2017** — Hauser FE et al. 2017. "Accelerated evolution and functional divergence of the dim light visual pigment... convergent RH1 evolution in fish." Mol Biol Evol 34(5):1155–1170.
- **Hunt2001** — Hunt DM et al. 2001. Spectral tuning and molecular evolution of rod visual pigments in deep-sea fish.
- **Blasing2000** — Bläsing OE, Westhoff P, Svensson P. 2000. "Evolution of C4 phosphoenolpyruvate carboxylase in Flaveria... Ser-780 is essential for high activity/low malate sensitivity." J Biol Chem 275:27917–27923.
- **Christin2007** — Christin P-A et al. 2007. "C4 photosynthesis evolved in grasses via parallel adaptive genetic changes." Curr Biol 17:1241–1247.
- **Besnard2009** — Besnard G et al. 2009. "Phylogenomics of C4 photosynthesis in sedges (Cyperaceae): multiple appearances and genetic convergence." Mol Biol Evol 26:1909–1919.
- **Christin2012** — Christin P-A et al. 2012. "Adaptive evolution of C4 photosynthesis through recurrent lateral gene transfer." / PEPC convergence follow-ups.
- **Liu2010** — Liu Y et al. 2010. "Convergent sequence evolution between echolocating bats and dolphins (Prestin)." Curr Biol 20:R53–R54.
- **Liu2011** — Liu Y et al. 2011. "Convergent evolution of KCNQ4 in echolocating bats and toothed whales." (hearing gene series)
- **Davies2012** — Davies KTJ et al. 2012. "Parallel signatures of sequence evolution among hearing genes in echolocating mammals." Mol Biol Evol 29:685–698.
- **Shen2012** — Shen Y-Y et al. 2012. "Parallel evolution of auditory genes for echolocation in bats and toothed whales."
- **Parker2013** — Parker J et al. 2013. "Genome-wide signatures of convergent evolution in echolocating mammals." Nature 502:228–231. (claim later shown mostly neutral)
- **Thomas2015** — Thomas GWC, Hahn MW. 2015. "Determining the null model for detecting adaptive convergence... echolocation reanalysis." Mol Biol Evol 32:1232–1236.
- **ZouZhang2015** — Zou Z, Zhang J. 2015. "No genome-wide protein sequence convergence for echolocation." Mol Biol Evol 32:1237–1241.
