# Tier 1 fixture — PEPC C3/C4 in sedges, **phenotypic** annotation

Sibling fixture to [`../pepc/`](../pepc/README.md), same gene/alignment/tree,
**different trait definition**. `../pepc/` uses ConDor's ("Besnard-derived")
**genotypic** annotation: `c4 = 1` iff the tip's ppc-1 sequence carries the
A780S substitution. This fixture uses ConDor's **phenotypic** annotation:
`c4 = 1` iff the plant species itself was independently classified as C4 by
anatomical/physiological survey — **fully decoupled from the ppc-1 sequence**,
avoiding the circularity of the genotypic annotation (see below).

## Why this fixture exists

Besnard et al. (2009)'s own gene tree colors "C4 ppc" branches by presence of
serine at maize-numbered position 780 (Fig. 3 caption) — i.e. genotypic C4
status is *defined* by that residue. Morel et al. (2024, ConDor) make this
explicit for the exact 78-tip dataset both fixtures share:

> "They then predicted the metabolism associated with the sedge PEPC
> sequences according to the presence or absence of the A780S mutation...
> this dataset, using the genotype-based annotation by Besnard et al. (2009),
> thus contains 78 sequences, 23 annotated as C4." (Morel et al. 2024, p.11)

And they show what this costs empirically: PCOC recovers 7/11 true positives
under the genotypic annotation but **0/11 under the phenotypic one**, on the
same underlying data (Morel et al. 2024, Table 2; "there is a perfect match
between the convergent clades and the mutations [under genotypic annotation],
which is advantageous for PCOC. However, on this dataset, PCOC fails with the
phenotypic annotation..."). Recovering position 780 under `../pepc/`'s
genotypic trait is close to definitional, not a real test of any method's
power. This fixture is the non-circular counterpart.

## Trait source: Bruhl & Wilson (2007)

`c4` here follows Bruhl, J.J. & Wilson, K.L. (2007), *"Towards a Comprehensive
Survey of C3 and C4 Photosynthetic Pathways in Cyperaceae,"* Aliso 23(1):
99–148 (`validation/tier1/references/pepc/Photosynthetic Pathways in
Cyperaceae.pdf`) — the same source Morel et al. (2024) cite for their own
"phenotypic" annotation. Two levels of evidence were used, per tip:

1. **Species-level individual record** in Bruhl & Wilson's Appendix 1
   checklist (3395 records), matched via the tip's **original Besnard/PCOC
   short code** (from `code_accession.py`, resolved via
   `align_cds/verify.tsv`'s accession column), not the fixture's modern
   GenBank-taxonomy tip name — several tips carry a different genus under
   Bruhl & Wilson's (2007-contemporary) nomenclature than under current
   taxonomy (e.g. `Machaerina_articulata` is listed as *Baumea articulata*;
   `Cyperus_sanguinolentus` as *Pycreus sanguinolentus*; `Cyperus_distichus`
   as *Volkiella disticha*; `Rhynchospora_colorata` as *Kyllinga colorata*;
   `Carex_capensis`/`Carex_phleoides`/`Carex_uncinata` as *Schoenoxiphium*/
   *Uncinia* species, both since folded into *Carex* by the Global Carex
   Group).
2. **Genus-level fallback** (Bruhl & Wilson's Table 1) when no individual
   record was found: used only for genera the paper classifies as
   "consistently C3" or "consistently C4" (i.e. not one of the five variable
   genera — *Abildgaardia*, *Cyperus*, *Eleocharis*, *Fimbristylis*,
   *Rhynchospora* — nor *Kyllinga*/*Pycreus*, which the paper also resolves
   at genus level as consistently C4).

## Dropped tips: *Eleocharis baldwinii* / *E. vivipara* (7 accessions)

Per Morel et al. (2024): *"we annotated each gene using the annotation of the
plant species in which it was sequenced from Bruhl and Wilson (2007), and we
removed the 7 genes from Eleocharis baldwinii and Eleocharis vivipara that
perform both C3 and C4 metabolisms... The 7 sequences from E. baldwinii and
E. vivipara were pruned from the provided phylogeny."* Both species are
documented C3/C4-intermediate at the species level (Bruhl & Wilson 2007), so
no single phenotypic label applies — dropped from alignment, both trees, and
trait files, exactly as ConDor does. This independently reproduces Morel et
al.'s own reported genotypic breakdown for these 7 accessions (5 C4 + 2 C3),
tip-for-tip, which we verified against our own data before building this
fixture (see prior session analysis, not reproduced here).

Fixture size: **71 tree tips (70 in-group + 1 outgroup)**, vs. `../pepc/`'s
78 (77 in-group + 1 outgroup) — 7 fewer, all *Eleocharis*.

## Result: 20 C4 / 50 C3 (in-group)

Morel et al. (2024) report 71 sequences / 22 C4 for their phenotypic dataset
(matches our 71-tip total exactly, including the outgroup). **Our own
tip-by-tip determination from Bruhl & Wilson (2007) gives 20 C4, not 22** —
a 2-tip discrepancy we could not resolve: ConDor's own per-tip phenotype
table (their supplementary Tables S1–S4) is not a species list, it's
mutation/method-level results, so it doesn't let us check individual calls
against theirs. We report our own sourced determination as-is rather than
force it to match Morel's aggregate count. If this matters for downstream
interpretation, the two extra ConDor C4 calls most likely live among the
`Fimbristylis`/`Eleocharis`-adjacent or *Volkiella*/*Kyllinga*/*Pycreus*
segregate-genus cases below, where genus reassignment made lookup least
certain.

## Discordant tips (genotypic ≠ phenotypic)

Besides the 7 dropped *Eleocharis* tips (genuinely intermediate, not
"discordant"), four tips carry a different label under phenotypic vs.
genotypic annotation — real, not lookup noise, each individually confirmed
against Appendix 1:

| tip | genotypic (A780S) | phenotypic (Bruhl & Wilson) | note |
|---|---|---|---|
| `Fimbristylis_dichotoma_FM208032` | C3 | **C4** | species-level call is C4; this accession's own sequence lacks A780S |
| `Fimbristylis_ferruginea_FM208034` | C3 | **C4** | same pattern |
| `Fimbristylis_littoralis_FM208036` | C3 | **C4** | same pattern |
| `Cyperus_distichus` (*Volkiella disticha*) | C4 | **C3** | opposite direction — carries A780S despite the species being phenotypically C3 |

The three `Fimbristylis` cases push the *same* direction as Morel's
documented "4 sequences flip from phenotypic-C4 to genotypic-C3" pattern.
`Cyperus_distichus`/*Volkiella disticha* flips the other way and isn't
described by that sentence at all — an unresolved oddity, reported rather
than smoothed over.

## Files

Same shape as `../pepc/`: `align/PEPC.fasta` (970 cols, amino acid),
`align_cds/PEPC.fasta` + `verify.tsv` (codon-level, for
`ortholog_characterizator`/FUBAR), `tree.nwk` (chronogram) /
`tree_substitution.nwk` (phylogram), both pruned from `../pepc/`'s trees with
`dendropy.Tree.prune_taxa_with_labels` (branch lengths of retained tips
unchanged), `my_traits.tsv` / `phenotype.tsv` (phenotypic `c4` label),
`ali_sp_names.txt` / `taxid.tsv` (71 tips, outgroup included),
`gene_ensembl.tsv` (identical placeholder, gene-level not tip-level).

Built by `/tmp/build_phenotypic.py` (not committed — one-off derivation
script; rerun against `../pepc/`'s current files if they change). Every
tip's genotypic label was checked against `../pepc/my_traits.tsv` before
writing (`assert`-verified, not just spot-checked).

## References

- Besnard G, et al. 2009. Phylogenomics of C4 photosynthesis in sedges
  (Cyperaceae): multiple appearances and genetic convergence. MBE
  26(8):1909–1919. doi:10.1093/molbev/msp103.
- Bruhl JJ, Wilson KL. 2007. Towards a comprehensive survey of C3 and C4
  photosynthetic pathways in Cyperaceae. Aliso 23(1):99–148.
  doi:10.5642/aliso.20072301.11.
- Morel B, et al. 2024. ConDor: accurate detection of convergent mutations.
  Genome Biol Evol 16(4):evae040. doi:10.1093/gbe/evae040.
