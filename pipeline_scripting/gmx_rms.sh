set -u

gmx='/home/speno/programs/bin'
# trajectory=$1 # Can be .pdb snapshots, not just .xtc or .trr
# reference_structure=$2

atom_group=13  # 2 for non-H. 13 for HEME. 

fluence_idx=3

# run_indices=($(seq 2 9))
# echo ${run_indices[@]}
run_indices=$(echo {1..10})

out_dir=$(realpath $(dirname $0)/output/$fluence_idx/)
mkdir -p $out_dir

# for f in input/_batches/batch_SH2_Zn/SH2_Zn-{2..9}.mol; do ./ac4dc "$f"; done
for run_idx in $run_indices; do
    trajectory=scripts/scattering/targets/converted_charges_10fs_${fluence_idx}_hemoglobin_solv_Hfix/hemoglobin_uppsula_moldstruct${run_idx}.pdb
    #reference_structure=scripts/scattering/targets/hemoglobin_solv_Hfix.gro # Can't use this structure as extends beyond GROMACS unit cell
    reference_structure=scripts/scattering/targets/hemoglobin_solv_Hfix_frame0.pdb 

    for expected_path in $trajectory $reference_structure; do
        if [ ! -f $expected_path ]; then
            echo $expected_path not found
            exit
        fi
    done
printf "$atom_group\n$atom_group\n" | "$gmx/g_rms" -s $reference_structure -f $trajectory -o $out_dir/$run_idx.xvg 
done 

xvg_args=""
for r in $run_indices; do
    xvg_args="$xvg_args${out_dir}/${r}.xvg "
done

echo "Plotting $xvg_args"
python3.9 $(dirname $0)/plotting/plot_xvg.py $xvg_args