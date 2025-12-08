set -u

# convert trajectories file to pdb snapshots

gmx='/home/speno/programs/bin' # path to where GROMACS is installed

xtc_file=$1
out_folder=$2
gromacs_file_path=$3
sample_interval=$4



i=1
out_path="$out_folder/output_$i.pdb"
while [ -f $out_path ]; do
    i=$((i+1))
    out_path="$out_folder/output_$i.pdb"
    if [ "$i" -ge 9999 ]; then
        echo "Critical error: Couldn't find a filename that doesn't exist!"
        exit
    fi 
done

#dt=0.0004
#dt=0.001
#dt=0.001
#out_folder="$(realpath ~/AC4DC/scripts/scattering/targets/)" 

handle=$(basename $working_folder)



#gro_file_handle=lys_example

struct_file=$gromacs_file_path



#printf "0" | "$gmx/trjconv"  -f $xtc_file -s $struct_file -pbc nojump -dt $dt  -o  $out_path # printf "0" selects group 0 (the whole system) to output.
printf "0" | "$gmx/trjconv"  -f $xtc_file -s $struct_file -pbc nojump -skip $sample_interval  -o  $out_path # printf "0" selects group 0 (the whole system) to output.

remove_sol='true'


# Not necessary, ignore in scattering code by default. But good to save space.
if $remove_sol; then 
    sed -i '/SOL/d' $out_path # -i option means edit in place
    sed -i '/TIP/d' $out_path # -i option means edit in place
    sed -i '/PEG/d' $out_path # -i option means edit in place
    sed -i '/SOD/d' $out_path # -i option means edit in place
    sed -i '/CLA/d' $out_path # -i option means edit in place
fi
 
