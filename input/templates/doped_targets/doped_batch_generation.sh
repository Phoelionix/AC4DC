set -x
#!/usr/bin/env

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
BATCH_GENERATION_SCRIPT=$SCRIPT_DIR/../../../generate_batch.py 

### Specific targets ###
# TARGET_STEM="SH2"
# TAG=("N" "S" "Fe" "Se" "Kr" "Ag" "Xe" "Gd")
# for t in ${TAG[@]}; do python3.9 $BATCH_GENERATION_SCRIPT "$SCRIPT_DIR/targets/${TARGET_STEM}_${t}.mol"; done

# ### All targets ###
for f in $SCRIPT_DIR/targets/*; do python3.9 $BATCH_GENERATION_SCRIPT "$f"; done