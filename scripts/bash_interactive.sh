set -x

#!/usr/bin/env
FILE=$1
echo $1
echo $FILE
basename "$FILE"

TARGET="$(basename -- $FILE)"
echo $TARGET
./ac4dc_v3 $FILE
python3 ~/AC4DC/scripts/generate_interactive.py $TARGET