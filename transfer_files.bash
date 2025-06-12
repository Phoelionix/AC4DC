TAG=("misc" "N" "S" "Gd" "Fe"   "Se"  )
for t in ${TAG[@]}; do scp -r passmores@spartan.hpc.unimelb.edu.au:/data/gpfs/projects/punim1915/AC4DC/output/__Molecular/Sprtn-SH2_${t} output/; done