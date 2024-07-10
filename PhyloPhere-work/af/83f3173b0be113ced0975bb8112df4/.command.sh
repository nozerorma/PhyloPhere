#!/bin/bash -euo pipefail
ct bootstrap \
    -a A3GALT2.Homo_sapiens.filter2.phy \
    -t /home/miguel/IBE-UPF/TFM/Master_projects/PhyloPhere/Data/6.CAAS_traitfiles/Ranked/df4_sp.tab \
    -s science.abn7829_data_s4.nex.pruned.tree.resampled.output \
    -o A3GALT2.Homo_sapiens.filter2.bootstraped.output \
    --fmt phylip-relaxed \
    --patterns 1,2,3         --max_bg_gaps 0         --max_fg_gaps 0         --max_gaps 0         --max_gaps_per_position 0.5         --max_bg_miss 0         --max_fg_miss 0         --max_miss 0
