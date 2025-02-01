set -x
#!/usr/bin/env

TAG=("N" "S" "Fe" "Se" "Kr" "I" "Gd")
for t in ${TAG[@]}; do python3.9 generate_batch.py "input/templates/doped_targets/targets/SH2_${t}.mol"; done
