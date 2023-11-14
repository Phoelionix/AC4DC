set -x
#!/usr/bin/env

# Full batch 
#for f in input/batch_lys/lys-{1..150}.mol; do ./ac4dc "$f"; done

# Specific files in batch
# NUM=(3 5)
# for n in ${NUM[@]}; do ./ac4dc "input/_batches/batch_lys/lys_-${n}.mol"; done

####################################

#for f in input/_batches/batch_ClpP/ClpP-{1..36}.mol; do ./ac4dc "$f"; done

TAG=("Ne" "Ar" "Zn" "Se" "Kr" "Zr" "Xe")
for t in ${TAG[@]}; do ./ac4dc "input/superheavy/SH_${t}.mol"; done