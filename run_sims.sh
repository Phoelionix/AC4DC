set -x
#!/usr/bin/env

#for f in input/batch_lys_all_light/lys_all_light-{133..150}.mol; do ./ac4dc "$f"; done
# ./ac4dc_gd input/galli/lys_galli_HF_L_edge.mol
# ./ac4dc_gd input/galli/lys_galli_HF_L_edge2.mol
# for f in input/batch_lys_full/lys_all_light-{130..150}.mol; do ./ac4dc "$f"; done
# for f in input/batch_lys_full/lys_full-{130..150}.mol; do ./ac4dc "$f"; done
#./ac4dc input/galli/lys_galli_LF_L_edge.mol
#./ac4dc input/galli/lys_galli_LF_L_edge.mol
# ./ac4dc input/royle/aluminium_royle_C
# ./ac4dc input/royle/aluminium_royle_B

# Fixing sims where transition energy code stuffed up
NUM=(143 144 145 149 150)
for n in ${NUM[@]}; do ./ac4dc "input/batch_lys_all_light/lys_all_light-${n}.mol"; done
for n in ${NUM[@]}; do ./ac4dc "input/batch_lys_full/lys_full-${n}.mol"; done

