# Input data

Not committed (large). Pull from the run:

```
RUN=correfoc-vs.s.upf.edu:/homes/users/mramon/scratch/2.Primates/2.Primates_results/CAAS_RESULTS/cancer_no_prune_multi05_toy/malignant_prevalence_toy_complete
scp mramon@$RUN/postproc/gene_filtering/filtered_discovery.tsv .
scp mramon@$RUN/data_exploration/2.CT/1.Traitfiles/contrast_hypotheses_pairs.tsv .
scp mramon@$RUN/scoring/position_scores.tsv .
scp mramon@$RUN/scoring/gene_scores.tsv .
```

Then:

```
Rscript da_frac_verify.R ../../subworkflows/SCORING/local/src \
  filtered_discovery.tsv contrast_hypotheses_pairs.tsv \
  position_scores.tsv gene_scores.tsv out 0.8 0.75
```

`contrast_hypotheses_pairs.tsv` and `position_scores.tsv` are kept here (small);
`filtered_discovery.tsv` and `gene_scores.tsv` are not.
