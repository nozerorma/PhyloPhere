# Tier 2 fixture: marine lifestyle in 58 placental mammals

Genome-wide, binary trait. Marine foreground: dolphin and killer whale (Cetacea), walrus and Weddell seal
(Pinnipedia), manatee (Sirenia): three independent origins with 2+2+1 tips. Input criteria and truth tiers are
in `../../CRITERIA.md`; the studies behind them are indexed in `../../references/marine.md`.

## Build

```bash
python3 validation/tier2/input/marine/scripts/build.py \
    --ucsc-dir <dir with the three UCSC files below> \
    --rds <UCSC62mammaltrees.rds> \
    --align-out <directory for the alignments>
```

`--stage trees` rebuilds only the tree and gene-tree files, `--stage align` only the alignments and tables.
The alignments (17,816 FASTA, about 630 MB) are written to `validation/phylophere_tier2/marine/align`, which git ignores.

| input | source |
|---|---|
| `UCSC62mammaltrees.rds` | RERconverge trees object of Partha et al. 2019: 19,149 gene trees, 62 tips, master tree; public folder `Trees_GitHub/Mammal_62species` of the Clark lab |
| `knownCanonical.exonAA.fa.gz` | `hgdownload.soe.ucsc.edu/goldenPath/hg19/multiz100way/alignments/`, 459,332,005 bytes; exon-wise amino acid alignments of 100 species |
| `knownCanonical.txt.gz`, `kgXref.txt.gz`, `kgAlias.txt.gz` | `hgdownload.soe.ucsc.edu/goldenPath/hg19/database/`; canonical transcripts, symbols and aliases (20,553,777 bytes for the last) |

## Files

| file | content |
|---|---|
| `tree.nwk` | master tree of the object, 58 placentals, rooted between Atlantogenata and the other placentals, dated as an ultrametric chronogram (`ape::chronos`, lambda 1, correlated rates, root age 1; the pipeline assumes a time tree) |
| `tree_substitution.nwk` | the same tree in substitution units; it is the input of the dating |
| `gene_trees.nwk` | `GENE<TAB>newick`, pruned to the 58 placentals; 19,149 genes |
| `rer_genes.tsv` | gene and number of placental species in its tree |
| `<align-out>/<GENE>.fasta` | stitched exons of the canonical transcript, 58 placentals, tips named `Genus_species`, columns without data in any of the 58 removed, final stop column removed, genes with fewer than 15 species not written |
| `my_traits.tsv` | `species`, `marine`, `control`, `order` |
| `prune/marine_tips.txt` | the five marine tips; pruned in the terrestrial-control run |
| `taxid.tsv`, `ali_sp_names.txt` | NCBI tax ids (live lookup) and the 58 species names |
| `gene_ensembl.tsv` | gene, hg19 span of the canonical transcript, strand, human protein length, UniProt accession |
| `unmapped_genes.txt`, `gene_aliases.tsv` | RER genes without an alignment; RER genes matched through an alias (RER symbol, current symbol) |
| `species_map.tsv` | tip label of the RER object and species name |

The template of the GUI project is `../../templates/tier2_marine.json`.

## Conventions

- **Gene assignment.** The identifiers of the UCSC alignments (`uc0xxxxx`) are not in the current gene tables.
  Each UCSC transcript is assigned by coordinates to the canonical transcript of `knownCanonical.txt.gz` that covers
  the largest share of its coding exon bases (smallest one on a tie); canonical transcripts of named genes, then those
  with a UniProt accession, are preferred over unnamed `ENSG` genes, which overlap them. The gene symbol is read from
  `kgXref.txt.gz`. A RER gene is matched by symbol; when the current gene set does not know the symbol, by the single
  current symbol that `kgAlias.txt.gz` gives it (`gene_aliases.tsv` lists these). A gene is written when its
  alignment has at least 15 species after cleaning. Of 19,124 RER genes with at least 15 placental species, 17,824 are
  assigned and 17,816 written; the 1,300 unassigned (`unmapped_genes.txt`) are `LOC`, `LINC` and `MIR` entries, cDNA
  clones, pseudogenes and older symbols that `kgAlias` does not carry (`NHP2L1`, `ATP5D`, `TCEB1`, `HIST1H…`), plus
  73 whose alias is ambiguous or already held by another gene.
- **Numbering.** The human row of each FASTA, ungapped, is the hg19 canonical protein; its position is the
  reference numbering of the truth sets. The alignment column of a position is the column of that residue in the
  FASTA. The transcript is the one of the 2016 UCSC alignment, which can differ from the current canonical isoform
  (NEIL1: 476 residues).
- **Residues.** UCSC writes the stop codon as `Z`; the final stop is removed and any other non-standard symbol is
  written as `X`.
- **Species names.** The walrus is `odoRosDiv1` in the alignment and `odoRosDi` in the RER object. NCBI has no
  species *Tupaia chinensis*: the Chinese tree shrew is written `Tupaia_belangeri`, the species that contains the
  subspecies *T. b. chinensis*. `canFam3` is `Canis_lupus`, `musFur1` is `Mustela_putorius`.
- **Taxonomic column.** `order` is the column named by `taxon_of_interest`. Cetaceans are in Artiodactyla.
- **Controls.** `control` marks aardvark, alpaca, Bactrian camel, little brown bat and David's myotis, the
  terrestrial control foreground of Chikina et al. 2016; the second trait row prunes the marine tips.
- **Alignment location.** `alignment_dir` of the template points to `validation/phylophere_tier2/marine/align`. A run on
  a cluster needs the alignment directory and the other fixture files copied there and the paths of the template changed.
- **Pinniped ancestral branch.** The template codes tips only. The marine or non-marine coding of the
  ancestral pinniped branch is not set in the template.

## Toy

`toy_genes.tsv` and `scripts/make_toy.py` define a 121-gene toy of the fixture for pipeline checks (`../../reports/marine_toy_results.md`). It is not the genome-wide case.

- **Genes of interest (81)**: listed in `toy_genes.tsv` with the reason of each: the experimental T1 genes (RHO, APPL1, NEIL1), 46 genes reported in at least two of the four replication lists of the T2 tier (Zhou 2015, Yuan 2021, Xu 2025, Huang 2024), 15 genes of Table 1 of Foote et al. 2015, 12 T3 genes named by Chikina et al. 2016 and 7 genes whose alignments are in the supplement of that paper. The rule that chose them beyond the reason column is not recorded.
- **Background genes (40)**: drawn from the alignments that are not genes of interest with a seeded shuffle (`--seed`, 1998 by default; `--n-background`, 40 by default).
- **Outputs**: `toy/genes.tsv` (gene and role), `toy/gene_trees.nwk` and `toy/gene_ensembl.tsv`, restricted to the 121 genes, and the toy alignment directory (`--out-align`), which git ignores. Species, tree, traits and tax ids are those of the fixture.

```bash
python3 validation/tier2/input/marine/scripts/make_toy.py \
    --align-dir validation/phylophere_tier2/marine/align --out-align validation/phylophere_tier2/toy/align
```

The template of the GUI project is `../../templates/tier2_marine_toy.json` (100 null cycles).
