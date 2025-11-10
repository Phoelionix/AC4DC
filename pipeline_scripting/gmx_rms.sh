set -u

gmx='/home/speno/programs/bin'
# trajectory=$1 # Can be .pdb snapshots, not just .xtc or .trr
# reference_structure=$2

atom_group=13  # 0 for whole system. 2 for non-H. 13 for HEME. 

fluence_indices=$(echo {7..14})
#dt=0.0004
dt=0.00005

reuse_xvg=false # skip to plotting, using the .xvg files generated on last call. If false, old xvg files are deleted.

#TODO should average trajectories first



# for f in input/_batches/batch_SH2_Zn/SH2_Zn-{2..9}.mol; do ./ac4dc "$f"; done
reference_structure=scripts/scattering/targets/hemoglobin_solv_Hfix_frame0.pdb 
#reference_structure=scripts/scattering/targets/hemoglobin_solv_Hfix.gro # TODO!!!! Try to use again # Can't use this structure as extends beyond GROMACS unit cell

if ! $reuse_xvg; then
    rm $(dirname $0)/output/*/*.xvg
fi

xvg_folders=""
for fluence_idx in $fluence_indices; do
    xvg_out_dir=$(realpath $(dirname $0)/output//$fluence_idx/)
    xvg_folders="$xvg_folders${xvg_out_dir}/ "
    mkdir -p $xvg_out_dir
    for trajectory in scripts/scattering/targets/converted_charges_10fs_${fluence_idx}_hemoglobin_solv_Hfix/output_*.xtc; do
        [ -e "$trajectory" ] || continue

        tmp=${trajectory##*/output_}
        MD_run_idx=${tmp%.xtc}
        
        for expected_path in $trajectory $reference_structure; do
            if [ ! -f $expected_path ]; then
                echo $expected_path not found
                exit
            fi
        done

        if ! $reuse_xvg; then
            rm -f tmp.xtc
            printf 0 | "$gmx/trjconv" -f $trajectory -s $reference_structure  -o  tmp.xtc  -pbc nojump -dt $dt


            printf "$atom_group\n$atom_group\n" | "$gmx/g_rms" -s $reference_structure -f tmp.xtc -o $xvg_out_dir/$MD_run_idx.xvg 
            rm -f tmp.xtc
        fi
    done

done 


echo "Plotting"
python3.9 $(dirname $0)/plotting/plot_all_rms.py $atom_group $xvg_folders
#rm $(dirname $0)/output/*/*.xvg

