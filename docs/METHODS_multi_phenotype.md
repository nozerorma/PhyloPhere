# Materials and Methods

## Identification of Convergent Amino Acid Substitutions across Two Primate Phenome Projects

---

### Overview

We applied a unified comparative-genomics pipeline — PhyloPhere (Ramon et al., in prep.) — to identify convergent amino acid substitutions (CAAS) associated with phenotypic variation across two independent biological projects in primates. Both projects exploit the same genomic resource: 16,133 multiple sequence alignments (MSAs) of orthologous protein-coding genes derived from the most comprehensive primate genomic catalogue available to date, encompassing 233 species (Kuderna et al., 2023). The primate time-calibrated species phylogeny from the same study was used as the reference tree throughout.

Both analyses were executed through PhyloPhere's Nextflow-based orchestration framework (Di Tommaso et al., 2017). The complete analytical configuration is encoded in `run_phenotypes.sh`, which defines two phenotype classes with distinct pre-processing requirements but a shared downstream CAAS discovery and post-processing chain. Post-processing was run in **filter mode**, and ancestral state reconstruction (ASR) for disambiguation was run in **precomputed mode** using pre-computed ASR caches to reduce computational overhead.

---

## Project 1: Comparative Oncology in Primates

### Background

This project is an ongoing collaboration between the IBE-UPF group and the Comparative Oncology Alliance, aimed at uncovering the genomic determinants of differential cancer susceptibility across the primate clade. Neoplastic disease data were compiled from the Species360 Zoological Information Management System (ZIMS), a global veterinary database covering captive primate populations.

### Phenotypic data

Two complementary tumour-burden phenotypes were defined from necropsy records:

- **`neoplasia_prevalence`**: species-level prevalence of all neoplastic disease (benign and malignant tumours combined), expressed as the proportion of adults found with neoplasia at necropsy.
- **`malignant_prevalence`**: species-level prevalence restricted to malignant tumours only.

Both are frequency/proportion phenotypes derived from individual clinical records. Two auxiliary population-size variables were carried through the analysis:

- **`n_trait`** (`adult_necropsy_count`): total number of adult necropsies per species, serving as the population denominator for prevalence estimation and confidence-interval (CI) calculation.
- **`c_trait`** (`neoplasia_necropsy` / `malignant_necropsy`): absolute case count per species, matched to the primary phenotype (swapped accordingly when each trait is used as primary vs. secondary).

A continuous **branch-colouring trait** (`branch_trait = LQ`, Longevity Quotient) was also carried through for downstream visualisation of evolutionary shifts. LQ serves two complementary roles in the cancer analyses:

1. **Phenotypic co-visualisation.** In the phenotype exploration report (`2.Phenotype_exploration.Rmd`), LQ is reconstructed onto the phylogenetic tree via ancestral state reconstruction and rendered as a continuous branch-colour overlay. This provides a qualitative readout of how longevity co-varies with cancer prevalence across the clade, allowing identification of lineages where both traits shift concordantly.
2. **Direct quantitative correlation.** The same exploration report generates scatter plots and Pearson/Spearman correlation statistics between LQ and each cancer phenotype (neoplasia prevalence or malignancy prevalence). Where the sample is sufficient, phylogenetically generalised least squares (PGLS) regression is additionally computed to account for shared ancestry. These correlations serve as a biological validation layer — long-lived species are expected to display elevated cancer suppression mechanisms — and provide interpretive context for the CAAS results.

LQ is deliberately excluded from `secondary_trait` to avoid conflating the CAAS contrast-group definition (which already uses the cross-cancer phenotype as secondary). Its quantitative relationship with cancer phenotypes is fully captured at the reporting and exploration stage.

Both phenotypes were run twice — each as the primary trait and as a secondary (co-variate) trait in the other's run — enabling cross-phenotype characterisation of CAAS.

### Species pruning

Because prevalence estimates are unreliable when necropsy sample sizes are small, species failing per-trait quality thresholds were excluded prior to analysis using curated exclusion lists:

- `neoplasia_exclude.txt`: species excluded from the neoplasia analysis.
- `malignant_exclude.txt`: species excluded from the malignancy analysis.

Exclusion lists are independent, reflecting the fact that data completeness differs between the two phenotypes; a species may be excluded from one analysis while being retained for the other.

### Contrast selection

Contrast groups (sets of species with high vs. low phenotypic values) were formed using PhyloPhere's `CONTRAST_SELECTION` subworkflow. Because both phenotypes are derived from count data (proportions with known denominators), **confidence-interval (CI) estimation** was enabled via `n_trait` and `c_trait`. This allows the contrast selection algorithm to account for sampling uncertainty in prevalence estimates when classifying species into phenotypic groups, ensuring that only species with sufficiently precise prevalence estimates contribute to the contrast.

---

## Project 2: Genetic Basis of Dietary Diversity in Primates

### Background

This project constitutes the continuation of the Master's thesis of **María Sánchez Bermúdez** (MU Bioinformàtica i Bioestadística, Universitat Oberta de Catalunya, defended January 2026), carried out in close collaboration with the IBE-UPF group under the supervision of **Dr. David de Juan Sopeña** and **Dr. José Francisco Sánchez Herrero**. The project studies the genomic mechanisms underlying the diverse dietary strategies observed across 159 primate species, spanning 16 families, with a particular focus on the differential consumption of meat (vertebrates), invertebrates (insects), leaves, fruit, seeds, and fermented fruit (ethanol). The study also explores the genetic bases of ethanol tolerance, given that multiple primate species regularly consume fermented fruit in the wild.

### Trait file format and preparation

The PhyloPhere pipeline requires trait data in **wide CSV format**: one row per species, each phenotype as a distinct column, comma-separated, with mandatory columns `species` and `family` (or the column specified by `taxon_of_interest`). The pipeline reads files using `read.csv(file, sep=",")` and validates that the `traitname` parameter corresponds to an existing column name in the loaded dataframe.

The trait file used for CLASS 2 runs (`diet_traitfile.csv`, produced by María Sánchez Bermúdez) is already in **wide format** (one row per species, each phenotype as a distinct column), which matches the structural requirement of the pipeline. However, the file is **semicolon-separated**, whereas `commons.R` assumes comma separation for files with a `.csv` extension. A one-line conversion is therefore required before first use:

```r
write.csv(read.csv2("diet_traitfile.csv"), "diet_traitfile_comma.csv", row.names = FALSE)
```

The resulting comma-separated file was used as the input trait CSV for all CLASS 2 pipeline runs.

### Phenotypic data

María Sánchez Bermúdez compiled and curated dietary composition data from published literature and established databases (primarily EltonTraits; Wilman et al., 2014), covering 159 primate species. She defined the following phenotypic variables and indices:

**Raw dietary fractions** (percentage of diet derived from each food source, drawn from EltonTraits and primary-literature sources):


| Variable      | Description                       |
| --------------- | ----------------------------------- |
| `Diet.Inv`    | Invertebrate (insect) consumption |
| `Diet.Vert`   | Vertebrate (meat) consumption     |
| `Diet.Fruit`  | Fruit consumption                 |
| `Diet.Nect`   | Nectar consumption                |
| `Diet.Seed`   | Seed consumption                  |
| `Diet.PlantO` | Other plant material consumption  |

**Trophic binary flags**:


| Variable | Description                  |
| ---------- | ------------------------------ |
| `Herb`   | Herbivore flag (binary: 0/1) |
| `Carn`   | Carnivore flag (binary: 0/1) |

**Trophic specialisation indices** (derived by María Sánchez Bermúdez; quantify the degree of specialism towards a given dietary guild relative to phylogenetic expectation):


| Index          | Description                         |
| ---------------- | ------------------------------------- |
| `frug_idx`     | Frugivory index                     |
| `fol_idx`      | Folivory index                      |
| `ins_idx`      | Insectivory index                   |
| `herb_idx`     | Herbivory index                     |
| `omn_idx`      | Omnivory index                      |
| `omn_spec_idx` | Omnivory specialisation index       |
| `folfrug_idx`  | Folivore–frugivore composite index |

**Dietary diversity metrics** (entropy-based):


| Variable       | Description                                       |
| ---------------- | --------------------------------------------------- |
| `Shannon_norm` | Shannon diversity of diet composition, normalised |
| `Shannon_spec` | Shannon diversity at species level                |

**Ethanol exposure**:


| Variable  | Description                                                                                     |
| ----------- | ------------------------------------------------------------------------------------------------- |
| `Ethanol` | Level of ethanol consumption/exposure (0–3 ordinal scale, from field and captive observations) |

### No population-size metrics

Unlike the cancer project, no individual-level clinical records are available for the dietary phenotypes: dietary data are observational estimates derived from the literature. Accordingly, **no `n_trait` or `c_trait`** were defined, CI calculation was disabled in contrast selection, and **no species pruning** was applied. Each phenotype was run independently through the pipeline in a serial fashion.

---

## Shared Analytical Pipeline

Both projects share an identical downstream analytical chain, described below.

### 1. Data reporting and exploration

For each phenotype, an exploratory HTML report was generated using PhyloPhere's `REPORTING` subworkflow, summarising the distribution of phenotypic values across species and families, flagging potential outliers, and visualising phylogenetic patterns. For the cancer project, reporting was run after pruning; for the diet project, it was run on the unfiltered data.

### 2. Contrast selection

Species were classified into two contrast groups (high vs. low phenotypic value) using PhyloPhere's `CONTRAST_SELECTION` subworkflow. The algorithm accounts for phylogenetic structure to identify sets of species with convergently extreme phenotypic values. For the cancer phenotypes, CI support was enabled; for the diet phenotypes, CI was disabled (see above).

### 3. CAAS discovery (CAAStools)

Convergent amino acid substitutions were identified using CAAStools (Barteri et al., 2023), executed through PhyloPhere's `CT` subworkflow, which chains three tools:

- **Discovery** (`ct_tool = discovery`): identifies positions in each MSA where the amino acid state systematically differs between the high and low phenotype groups, computing a binomial p-value for the observed pattern.
- **Resample** (`ct_tool = resample`): generates *N* random permutations of species–group assignments, preserving group sizes and respecting phylogenetic and phenotypic structure, to build a null distribution of CAAS counts per position.
- **Bootstrap** (`ct_tool = bootstrap`): for each real CAAS position, counts the fraction of random permutations yielding the same substitution, yielding an empirical p-value.

The number of permutation cycles was set to **100** in toy mode (for development and validation) and **1,000,000** in full-scale runs.

### 4. CT signification

Empirical bootstrap p-values were computed and aggregated across genes using PhyloPhere's `CT_SIGNIFICATION` subworkflow. This step produces a per-gene summary of CAAS burden and a global metadata table (`global_meta_caas.tsv`) that feeds into the disambiguation step.

### 5. CT disambiguation

The `CT_DISAMBIGUATION` subworkflow filters CAAS by assessing evolutionary convergence at each candidate position using ancestral state reconstruction (ASR). ASR was run in **precomputed mode** (`ct_disambig_asr_mode = precomputed`), leveraging pre-computed ancestral probability distributions to maximise computational efficiency. This step produces the `caas_convergence_master.csv` file, which serves as input for post-processing.

### 6. CT post-processing

The `CT_POSTPROC` subworkflow was run in **filter mode** (`caas_postproc_mode = filter`), applying a gene-level filtering procedure to retain only genes with statistically significant CAAS burden (based on the combined evidence from bootstrap p-values and convergence assessment). Outputs include a filtered discovery table, cleaned background gene lists, and per-gene HTML characterisation reports.

### 7. Functional enrichment

Two complementary enrichment strategies were applied to the filtered gene set:

- **Over-Representation Analysis (ORA)** via PhyloPhere's `ORA` subworkflow, using WebGestalt (Wang et al., 2017) against Gene Ontology (GO) and KEGG databases.
- **STRING protein–protein interaction enrichment** via PhyloPhere's `STRING` subworkflow, using the rbioapi R package to query the STRING database (Szklarczyk et al., 2023).

### 8. CT accumulation

Gene-level significance of CAAS accumulation — i.e., whether individual genes carry more CAAS than expected by chance given the background gene set — was assessed using PhyloPhere's `CT_ACCUMULATION` subworkflow. This subworkflow performs a permutation test by randomly re-assigning CAAS across background genes. The number of randomisations was set to **1,000** in toy mode and **10,000** in full-scale runs.

---

## Computational Infrastructure

All analyses were run locally using PhyloPhere's `local` Nextflow profile. The pipeline was executed from the project root using `nextflow run main.nf`. Work directories and results were written to a dedicated external storage volume. Nextflow Tower monitoring was enabled for all runs.

---

## References

- Barteri, F., et al. (2023). CAAStools: a toolbox to identify and test Convergent Amino Acid Substitutions. *Bioinformatics*, 39(11), btad623.
- Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology*, 35, 316–319.
- Kuderna, L.F.K., et al. (2023). A global catalog of whole-genome diversity from 233 primate species. *Science*, 380(6648), 906–913.
- Sánchez Bermúdez, M. (2026). *Estudio de los mecanismos genéticos asociados a los distintos patrones de dieta en primates mediante genómica comparativa*. Master's Thesis (MU Bioinformàtica i Bioestadística, UOC). Supervised by Dr. D. de Juan Sopeña & Dr. J.F. Sánchez Herrero.
- Szklarczyk, D., et al. (2023). The STRING database in 2023: protein–protein association networks and functional enrichment analyses for any of 14 090 organisms. *Nucleic Acids Research*, 51(D1), D638–D646.
- Valenzuela, A., et al. (2025). Characterising the genome–phenome relationship across 224 primate species. *Nature*, in press.
- Wang, J., et al. (2017). WebGestalt 2017: a more comprehensive, powerful, flexible and interactive gene set enrichment analysis toolkit. *Nucleic Acids Research*, 45(W1), W130–W137.
- Wilman, H., et al. (2014). EltonTraits 1.0: Species-level foraging attributes of the world's birds and mammals. *Ecology*, 95(7), 2027.
