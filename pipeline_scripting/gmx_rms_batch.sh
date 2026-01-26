set -u

folder_stem=hemoglobin_2QSP
#dt=0.00005
dt=0
reuse_xvg=true # skip to plotting, using the .xvg files generated on last call. If false, old xvg files are deleted.
reuse_xtc=true # Reuse unwrapped xtc if generated.

#reference_structure=~/AC4DC/scripts/scattering/targets/2qsp_unit_Hfix.pdbv # This won't have the solvent atoms added
reference_structure=~/cr-md/2qsp_template/2qsp_unit_Hfix.gro
index_file=~/cr-md/2qsp_template/index.ndx

#for atom_selection in "system" "HEME" "Fe"; do  
for atom_selection in "system" "HEME"; do  
    for trajectory_folder in $(dirname $0)/../scripts/scattering/targets/$folder_stem-*/; do
        bash $(dirname $0)/.gmx_rms_inner_logic.sh $atom_selection $(realpath $trajectory_folder) $dt $reuse_xvg $reuse_xtc $reference_structure $index_file
done;   
done;


#rm $(dirname $0)/output/*/*.xvg

