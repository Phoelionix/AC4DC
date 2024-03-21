set -x
#!/usr/bin/env

TAG=("N" "S" "Fe" "Se" "Kr" "I" "Gd")
for t in ${TAG[@]}; do python generate_batch.py "input/superheavy/SH2_${t}/SH2_${t}.mol"; done
