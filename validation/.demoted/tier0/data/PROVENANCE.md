# Tier 0 substitution-matrix data

PAML-format empirical amino-acid rate matrices, copied verbatim from the PAML
4.10.10 conda package (`paml-4.10.10-h7b50bb2_0/dat/`) on the UPF cluster.

| file | model | reference | Tier 0 use |
|------|-------|-----------|------------|
| `wag.dat` | WAG | Whelan & Goldman 2001, Mol Biol Evol 18:691 | **default generative matrix.** Combined with site-specific Dirichlet profiles it is a different model family from the pipeline's LG-single-profile ASR/FADE. |
| `dayhoff.dat` | Dayhoff | Dayhoff et al. 1978 | alternate generative matrix, further from LG; use for a stronger misspecification check. |
| `lg.dat` | LG | Le & Gascuel 2008, Mol Biol Evol 25:1307 | **do not simulate under this** — it is (close to) the inference model. Kept only for a deliberate "model matches" sanity run where the null should be its easiest. |

Format: 19 lines lower-triangular exchangeabilities (line *k* has *k* entries),
then one line of 20 equilibrium frequencies, then `//` and notes. Amino-acid
order `A R N D C Q E G H I L K M F P S T W Y V`. Parsed by `tier0.model.load_rate_matrix`.
