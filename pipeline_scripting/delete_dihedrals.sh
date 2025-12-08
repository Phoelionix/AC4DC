set -u

working_folder=$1 # full path

topology_folder=$working_folder/topology/

for topol in $topology_folder/topol_*.itp; do
    awk '$0=="[ dihedrals ]"{print; exit} {print}' $topol > ${topol}_tmp
    mv ${topol}_tmp $topol
done