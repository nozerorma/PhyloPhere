# Tier 2 input

## `species_tree_meredith_plus.nwk`

Species tree of the 58 placental mammals of the RER set ("Meredith+" topology), rooted, with UCSC
assembly codes as tip labels.

- **Source:** the consensus topology printed in the supplementary text of Meyer et al. 2018
  (*Science*, doi:[10.1126/science.aap7714](https://doi.org/10.1126/science.aap7714);
  `references/marine/meyer2018_aap7714/supplementary_data/NIHMS995232-supplement-1.pdf`, paragraph
  "Our consensus topology, with branch lengths re-estimated from 10 randomly selected genes").
  The extraction removes line breaks and spaces only.
- **Definition:** Meredith et al. 2011 topology with three local resolutions: star-nosed mole outside
  hedgehog and shrew; cow outside Tibetan antelope, sheep and goat; ursids outside mustelids and
  pinnipeds. Kowalczyk et al. 2020 (*eLife*, Figure 1 supplement 7) and Saputra et al. 2021 describe
  the RER trees with this topology.
- **Check:** the master tree of `UCSC62mammaltrees.rds`, pruned to these 58 tips and rooted at the
  armadillo, has Robinson-Foulds distance 0 to this tree. Against `hg19.100way.nh` the distance is 22
  (maximum 110).
- **Branch lengths:** those of the Meyer supplement (10 genes). They are not the branch lengths of the
  gene trees of the RER set; use the topology only.
- **Tip codes:** the walrus is `odoRosDi` here and `odoRosDiv1` in the UCSC 100-way files.
