#!/usr/bin/env bash

# Without args: performs scripts on most recent folder in __Molecular.
# With args: Performs scripts on each handle (output data folder name) passed as an arg.

scripts=( "generate_interactive.py" "spatial_bound.py" )

if [ $# -eq 0 ]; then
handles=$(ls output/__Molecular -Art | tail -n 1)
else
handles=$@ # all args
fi

for handle in  ${handles[@]};
do
    for script in ${scripts[@]};
    do
    python3.9 ~/AC4DC/scripts/$script $handle
    done
done