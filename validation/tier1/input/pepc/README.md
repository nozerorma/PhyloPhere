# Tier 1 fixture — PEPC C3/C4 in sedges

Single-gene, categorical-trait site-truth benchmark. C3 vs C4 photosynthesis in
Cyperaceae, PEPC (`ppc`) amino-acid alignment.

Source: **Besnard, G. et al. 2009**, *Phylogenomics of C4 photosynthesis in
sedges (Cyperaceae): multiple appearances and genetic convergence*, MBE
26(8):1909–1919, doi:10.1093/molbev/msp103.

## Build

```bash
python3 validation/tier1/input/pepc/scripts/build.py
```

Fetches the raw inputs (gitignored, under `besnard2009/`) and writes the
pipeline-shaped fixture: `align/PEPC.fasta`, `tree.nwk` (PhyML tree rooted on
*Chrysithrix* and time-scaled to an ultrametric chronogram via `ape::chronos`),
`tree_substitution.nwk` (raw phylogram), `my_traits.tsv`, `phenotype.tsv`,
plus synthetic `ali_sp_names.txt` / `gene_ensembl.tsv` / `taxid.tsv`. Only
`scripts/build.py`, `scripts/build_cds.py`, `scripts/code_accession.py` and
this README are committed.

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

## Fixture size and dropped taxon

The raw alignment has 79 tips (23 C4, 55 C3, 1 outgroup). `Abildgaar`, one of
the 55 C3 tips, has no resolvable GenBank accession (see below) and is
dropped from every output — alignment, tree, and trait files. Fixture size is
**78 tips (23 C4, 54 C3, 1 outgroup)**.

## GenBank accession resolution

`scripts/code_accession.py`'s `CODE_ACCESSION` table maps each of the
alignment's original PCOC/ConDor short codes (`Cyp_era1`, `Ele_bal2`, etc.)
to its Besnard et al. 2009 ppc-1 GenBank accession, cross-checked against the
paper's own Supplementary Table S1. 76 of 79 codes match a named species via
a single accession. Two short codes resolve to real species under a
different genus than the alignment's own labels suggest:

- **`Baumea`** = **FM207991**, whose GenBank record's `ORGANISM` field is
  ***Machaerina articulata*** (an older synonym, *Baumea articulata*, appears
  only in the record's `DEFINITION` line).
- **`Machaeri`** = **FM208044**, ***Machaerina scirpoidea***.

**`Abildgaar` has no resolvable accession**: it is not in Table S1, not in
either supplementary tree figure, and not reachable via any of the ppc-1
accession numbers not otherwise claimed by a named taxon in this dataset.

The 9 species with 2–4 ppc-1 accessions (*C. eragrostis* ×2, *C. ustulatus*
×2, *E. baldwinii* ×4, *E. gracilis* ×2, *E. limosa* ×3, *E. palustris* ×2,
*E. vivipara* ×3, *Fimbristylis* ×2 each for dichotoma/ferruginea/littoralis,
*Hellmuthia membranacea* ×2, *Rhynchospora globosa* ×2) keep every accession
as its own tip; no deduplication is performed. Which physical accession is
which existing code is determined by best-offset local protein identity
between each candidate accession's own GenBank `/translation` and the
existing tip's amino-acid sequence in `align/PEPC.fasta` (sliding the shorter
sequence across the longer one and taking the offset with the highest match
rate). Every one resolves to ≥98% identity at its best offset; *Eleocharis
limosa*'s three accessions (FM208022/23/24) are >99% identical to each other,
so which physical accession lands under which of its three existing codes is
nominal.

The paper's own ppc-1 sampling includes *Fimbristylis hygrophila* (FM207989,
C3); this alignment does not carry it (only *F. dichotoma*, *F. ferruginea*,
*F. littoralis* are present).

## Tip renaming (real species names)

`align/PEPC.fasta`, `tree.nwk`, `tree_substitution.nwk`, `my_traits.tsv`,
`phenotype.tsv`, `ali_sp_names.txt` and `taxid.tsv` all use real species
binomials as tip identifiers, not the original PCOC/ConDor short codes, so
PhyloPhere's own tax-id auto-inference has something real to look up.
`build.py`'s `_species_names()` derives the rename from
`code_accession.py`'s `CODE_ACCESSION` table by way of each accession's own
GenBank `ORGANISM` annotation, i.e. current NCBI-accepted taxonomy, not
necessarily Besnard 2009's own label (e.g. `Baumea` → `Machaerina_articulata`,
`Killinga` → `Rhynchospora_colorata`). `build_cds.py` re-derives the identical
map from the same table, so a lookup by the tip names in `align/PEPC.fasta`
always agrees with `code_accession.py`.

The 12 species with 2–4 ppc-1 accessions get an accession suffix so tips stay
unique, e.g. `Cyperus_eragrostis_FM208065` / `Cyperus_eragrostis_FM208066`;
the other 66 single-accession species use the plain binomial.

## Codon-level CDS (`build_cds.py`, `align_cds/PEPC.fasta`)

```bash
python3 validation/tier1/input/pepc/scripts/build_cds.py   # after build.py
```

`build_cds.py` back-translates `align/PEPC.fasta`'s 970-column amino-acid
alignment into a real codon alignment, for use as `ortholog_characterizator`'s
`--cds_dir` input (quality → translation/BMGE → phylogeny → positive
selection): it walks the amino-acid alignment column by column, emitting one
real codon per non-gap column (from `code_accession.py`'s per-tip GenBank
accession) and `"---"` per gap column, so the codon alignment shares
`align/PEPC.fasta`'s exact coordinate system.

All 78 tips get a real codon row; `align_cds/verify.tsv` records
`build_cds.py`'s own per-tip verification (coverage + exact-match flag
against `align/PEPC.fasta`, using `difflib.SequenceMatcher`'s coverage of the
fixture's row rather than `.ratio()`, since several accessions' real
sequences cover more of the gene than the fixture's 455-aa fragment). All 78
tips pass at 100% coverage.

## Positive selection (`ortholog_characterizator`)

`validation/tier1/run_ortholog_characterizator.sh pepc` runs
`ortholog_characterizator` (translation → BMGE → FUBAR) against
`align_cds/PEPC.fasta`, rooted on `tree.nwk` via `--psel_species_tree`.
`PHYLOGENY` is off by default in that launcher — `POSITIVE_SELECTION` always
uses this fixture's own `tree.nwk` directly, not a tree derived from
`PHYLOGENY`'s own gene-tree/species-tree inference. Output would land under
`oc_run/`, which does not currently exist in this fixture (the pipeline has
not been run against it).

## Caveats

- `align/PEPC.fasta` is amino-acid, not codon — PCOC and ConDor both
  redistribute, and run their own published analyses on, the amino-acid
  alignment; this fixture uses the same one for GUI-driven runs. The codon
  alignment (`align_cds/PEPC.fasta`) is a separate track for
  `ortholog_characterizator`'s CDS-based pipeline.
- Besnard 2009 lists 16 positive-selection codons; the truth set transcribes
  the confident subset (7 tiered + 4 weak). The rest are in the Besnard 2009
  Table 1 / Christin et al. 2007 (Genetics 177(3):1791–1808) Table 1.
- Only positions 780 (kinetics) and 665 (folding/activity) are functionally
  demonstrated (site-directed mutagenesis). The remaining truth-set positions
  are selection inferences from the source papers, not mutagenesis-confirmed.
- `gene_trees.nwk` at the fixture root uses the pre-rename PCOC/ConDor short
  codes (`Abildgaar`, `Cyp_era1`, etc.), not the species-name tips used by
  every other file in this fixture, and still includes `Abildgaardia`. The
  GUI template's `gene_trees` parameter points at this file as-is; it would
  need regenerating against the current tip names to match the rest of the
  fixture.

## GUI template wiring (`../../../../gui/templates/tier1_pepc.json`)

- `gene_ensembl.tsv` in this directory is `build.py`'s synthetic placeholder
  (fake genomic coordinates; the PEPC `human_protein_id` is maize's own
  UniProt id) — the template points `gene_ensembl_file` at it directly.
- `taxid.tsv` carries real NCBI tax_ids (`build.py`'s `_ncbi_tax_ids`, one
  lookup per base species, shared across a multi-accession group's sibling
  tips, via `bin/generate_taxid_map.py`'s `resolve_taxids`: live NCBI eutils
  first, falling back to `ete3`'s local taxonomy dump). The template still
  leaves `tax_id_file` blank, so PhyloPhere generates it itself at run time
  rather than reading `taxid.tsv` directly; `generate_taxid_map.py` resolves
  tip labels directly, including the accession-suffixed multi-accession tips,
  via a genus+species fallback.
