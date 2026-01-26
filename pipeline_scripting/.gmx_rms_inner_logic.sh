set -u


atom_selection=$1 
trajectory_folder=$2
dt=$3
reuse_xvg=$4
reuse_xtc=$5
reference_structure=$6
index_file=$7

gmx='/home/speno/programs/bin'

cd $(dirname $0)
# trajectory=$1 # Can be .pdb snapshots, not just .xtc or .trr
# reference_structure=$2

num_to_sample=99999


#reference_structure=../scripts/scattering/targets/hemoglobin_solv_Hfix_frame0.pdb 
#reference_structure=../scripts/scattering/targets/hemoglobin_solv_Hfix.gro 
#index_file=~/cr-md/hemoglobin_uppsala_template/index.ndx




if [ $atom_selection = "system" ]; then
    atom_group_idx=0
elif [ $atom_selection = "nonH" ]; then
    atom_group_idx=2
elif [ $atom_selection = "HEME" ]; then
    atom_group_idx=13
elif [ $atom_selection = "Fe" ]; then
    atom_group_idx=18
else
    echo "Invalid atom_selection " $atom_slection 
    exit
fi
#parent_dir=output/RMSD/$atom_selection/${width_fs}fs
xvg_out_dir=output/RMSD/$atom_selection/$(basename $trajectory_folder)
nojump_dir=output/unwrapped_trj/$(basename $trajectory_folder)/

if ! $reuse_xvg; then
    rm -f $xvg_out_dir/*.xvg
fi

#xvg_out_dir="$parent_dir/$fluence_idx/"
mkdir -p $xvg_out_dir $nojump_dir

i=0
for trajectory in $trajectory_folder/output_*.xtc; do
    [ -e "$trajectory" ] || continue

    tmp=${trajectory##*/output_}
    MD_run_idx=${tmp%.xtc}
    
    for expected_path in $trajectory $reference_structure; do
        if [ ! -f $expected_path ]; then
            echo $expected_path not found
            exit
        fi
    done

    tmp="${tmp}"
    tmp=${tmp##*/}
    tmp="${tmp%.xtc}"
    trj_basename=$tmp

    unwrapped_xtc=$nojump_dir/${trj_basename}-nojump.xtc

    if ! $reuse_xvg; then
        if [ ! -f $unwrapped_xtc ]  ||  [ ! $reuse_xtc = "true" ]; then 
            rm -f $unwrapped_xtc
            printf 0 | "$gmx/trjconv" -f $trajectory -s $reference_structure  -pbc nojump -dt $dt -n $index_file -o  $unwrapped_xtc 
            # TODO delete if empty
        fi
        printf "$atom_group_idx\n$atom_group_idx\n" | "$gmx/g_rms" -s $reference_structure -f $unwrapped_xtc -n $index_file -o $xvg_out_dir/$MD_run_idx.xvg 
    fi
    i=$((i+1))
    if [ "$i" -ge $num_to_sample ]; then
        break
    fi 
done
echo "Plotting"
xvg_folders=$xvg_out_dir
python3.9 plotting/plot_all_rms.py $(basename $trajectory_folder) $atom_selection $xvg_folders 


#rm $(dirname $0)/output/*/*.xvg

