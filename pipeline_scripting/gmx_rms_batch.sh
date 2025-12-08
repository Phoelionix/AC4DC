set -u

#dt=0.00005
dt=0
reuse_xvg=false # skip to plotting, using the .xvg files generated on last call. If false, old xvg files are deleted.
reuse_xtc=true # Reuse unwrapped xtc


for atom_selection in "HEME" "Fe"; do  
    for width_fs in 3 10; do 
        bash $(dirname $0)/.gmx_rms_inner_logic.sh $atom_selection $width_fs $dt $reuse_xvg $reuse_xtc

done;   
done;


#rm $(dirname $0)/output/*/*.xvg

