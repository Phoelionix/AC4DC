set -x
#!/usr/bin/env

# Background water files for use 
#Idea: Each unit cell corresponds to a unique water drop. We use a slightly different size to randomise.
# Note that water *inside* is generated too
# Ions are not included currently. Need to make a psf file and use -ion... though it will only include salt. And probably insignificant to scattering pattern regardless.
PDB="4et8"
OUTFOLDER="lys_water"
NUM=15



mkdir -p $OUTFOLDER
for ((i=1,T=500;i<=$NUM;i++,T+=5))
do 
    t=$(printf  %.3f\\n "$(($T))e-2")
    echo $t
    ./solvate -t $t -n 2 $PDB $OUTFOLDER/sol_$PDB-$i
done


#for t in {5,5.05,5.1,5.15,5.2,5.25,5.3,5.35,5.4,5.45,5.5}; do ./solvate -t $thickness$ -n "2 4et8 lys_water/sol_4et8-$t"