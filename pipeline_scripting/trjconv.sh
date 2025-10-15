set -u

# convert trajectories file to pdb snapshots

gmx='/home/speno/programs/bin' # path to where GROMACS is installed

working_folder=$1
gromacs_file_path=$2
out_folder=$3


dt=0.001
#out_folder="$(realpath ~/AC4DC/scripts/scattering/targets/)" 

handle=$(basename $working_folder)

num=0


#gro_file_handle=lys_example

xtc_file=$working_folder/output_${num}.xtc
struct_file=$gromacs_file_path



out_file=$working_folder/${handle}_moldstruct${num}.pdb
printf "0" | "$gmx/trjconv" -f $xtc_file -s $struct_file -dt $dt  -o  $out_file # printf "0" selects group 0 (the whole system) to output.

remove_sol='true'


# Not necessary, ignore in scattering code by default. But good to save space.
if $remove_sol; then 
    sed -i '/SOL/d' $out_file # -i option means edit in place
fi
 

# TODO in some script we need to make sure not overwriting 
mv $out_file $out_folder
