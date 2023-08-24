set -x
#!/usr/bin/env

#NUM=($(seq 25))
#NUM=(17 18 19 20 21 22 23 24 25)
NUM=(23 24 25)
for n in ${NUM[@]}; do python3.9 scripts/light_ions_plot.py "lys_full_fixed-${n}_1" "lys_all_light_fixed-${n}_2"; done

#NUM=(22)
#for n in ${NUM[@]}; do python3.9 scripts/light_ions_plot.py "lys_full_fixed-${n}_1" "lys_all_light-${n}_patched"; done