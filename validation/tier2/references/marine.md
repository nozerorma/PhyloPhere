# References: marine mammals (Tier 2)

Layout: one directory per study under `marine/`, named `<first author><year>_<journal id>`, with the
supplementary files of that study in `supplementary_data/`. Original file names are kept. The primary
papers themselves are not stored here (see the last section).

Evidence level in the last column of each table: `site` (amino acid position list), `gene` (gene list or
gene-level statistic), `assay` (functional experiment), `input` (data that can feed a run).

## Input data of the RER set

| Directory | Study | Content | Level |
|---|---|---|---|
| `marine/partha2019_msz107/` | Partha et al. 2019, *Mol Biol Evol* 36(8):1817. doi:[10.1093/molbev/msz107](https://doi.org/10.1093/molbev/msz107) | Supplementary figures and Tables S1 to S12 (GO enrichment of the top accelerated genes, subterranean and marine). No trees and no alignments. | gene |
| `marine/chikina2016_msw112/` | Chikina, Robinson, Clark 2016, *Mol Biol Evol* 33(9):2182. doi:[10.1093/molbev/msw112](https://doi.org/10.1093/molbev/msw112) | `species_key_mammals.txt` (58 species, UCSC assembly code and common name), amino acid alignments of 7 genes (ALDH3B2, AQP10, CAPN8, MYH7, MYOF, NEB, PGAM2) and `Supp_Tables_Legends.xlsx` (Tables S1 to S8). | input, gene |

The marine tips of the species key are `turTru2` (dolphin), `orcOrc1` (killer whale), `odoRosDi` (Pacific
walrus), `lepWed1` (Weddell seal) and `triMan1` (manatee). The 19,149 gene trees and the full alignments
of the RER set are not in these supplements.

## Site-level evidence

| Directory | Study | Files | Level |
|---|---|---|---|
| `marine/foote2015_ng3198/` | Foote et al. 2015, *Nat Genet*. doi:[10.1038/ng.3198](https://doi.org/10.1038/ng.3198) | `MOESM21` (PDF): Supplementary Tables 1, 2, 4, 5, 6, 8, 14 and Supplementary Figure 1; Table 14 lists the positions of parallel substitutions. `MOESM22` to `MOESM28` (xlsx): Tables 3, 7, 9, 10, 11, 12, 13 (positive selection per branch, GO enrichment). | site, gene |
| `marine/zhou2015_srep16550/` | Zhou et al. 2015, *Sci Rep* 5:16550. doi:[10.1038/srep16550](https://doi.org/10.1038/srep16550) | `srep16550-s1.pdf`: S1 (133 parallel substitutions in 132 genes), S2 (26 genes), S3 (70 genes). | site, gene |
| `marine/yuan2021_pnas2106080118/` | Yuan et al. 2021, *PNAS*. doi:[10.1073/pnas.2106080118](https://doi.org/10.1073/pnas.2106080118) | SI appendix and datasets `sd01` to `sd08`; the position list is SI Table S24 (195 substitutions in 172 genes). | site |
| `marine/xu2025_xinn100798/` | Xu et al. 2025, *The Innovation*. doi:[10.1016/j.xinn.2025.100798](https://doi.org/10.1016/j.xinn.2025.100798) | `mmc2.xlsx` (Table S7, 84 sites in 82 genes). The 23 MB `mmc1.pdf` is not kept. | site, assay |
| `marine/huang2024_isci110394/` | Huang et al. 2024, *iScience*. doi:[10.1016/j.isci.2024.110394](https://doi.org/10.1016/j.isci.2024.110394) | `mmc1.pdf`, `mmc2` to `mmc10` (Table S15: 232 sites in 215 genes). Cetaceans and sirenians only. | site |
| `marine/chen2025_caap_men70052/` | Chen and Zou 2025 (CAAP), *Mol Ecol Resour*. doi:[10.1111/1755-0998.70052](https://doi.org/10.1111/1755-0998.70052) | `s001.xlsx` (Tables S11, S12: sites of about 15 genes, alignment-column numbering), `s002.pdf`. | site |

The numbering frame of each list differs and is not always stated; the lists cannot be compared before
their coordinates are remapped to one reference.

## Gene-level evidence

| Directory | Study | Files | Level |
|---|---|---|---|
| `marine/huelsmann2019_aaw6671/` | Huelsmann et al. 2019, *Sci Adv* 5:eaaw6671. doi:[10.1126/sciadv.aaw6671](https://doi.org/10.1126/sciadv.aaw6671) | `aaw6671_sm.pdf` and Tables S1 to S5 in one xlsx (genes lost in cetaceans). | gene |
| `marine/meyer2018_aap7714/` | Meyer et al. 2018, *Science*, convergent losses of *PON1*. doi:[10.1126/science.aap7714](https://doi.org/10.1126/science.aap7714) | Supplementary text (PDF), Tables S1 (per-gene loss test), S4, S5, S7 (enrichment), S9, S10 (sequences), S13 (enzyme assays). Tables S2, S3, S6, S8, S11 and S12 are not in the set. | gene, assay |
| `marine/treaster2021_msab226/` | Treaster et al. 2021 (TRACCER), *Mol Biol Evol*. doi:[10.1093/molbev/msab226](https://doi.org/10.1093/molbev/msab226) | `Treaster_TRACCER_Supplemental.xlsx`: S1 (marine genes), S2 (marine GO terms), S3 (marine control). | gene |
| `marine/saputra2021_msab068/` | Saputra et al. 2021 (permulations), *Mol Biol Evol*. doi:[10.1093/molbev/msab068](https://doi.org/10.1093/molbev/msab068) | Supplementary walkthrough and one supplementary figure; the study publishes no gene sets. | gene |

## Not stored here

The papers (PDF) are not in this directory. The studies above need their main text for the methods and
for the numbering of each list. Experimental papers cited in the search reports are not listed in any
table above because their supplements carry no position table: the rhodopsin resurrection papers
(Dungan and Chang) and the review by Clark et al. 2025,
doi:[10.1146/annurev-genet-030325-041233](https://doi.org/10.1146/annurev-genet-030325-041233).

Search reports (outside the repository): `~/.cc/tier2_search/{A_site_level,B_functional,C_gene_level}.md`.
