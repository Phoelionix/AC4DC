set -u

# convert trajectories file to pdb snapshots

gmx='/home/speno/programs/bin' # path to where GROMACS is installed

xtc_file=$1
out_path=$2
gromacs_file_path=$3
sample_interval=$4


handle=$(basename $working_folder)


struct_file=$gromacs_file_path


rm -rf $out_path 
printf "0" | "$gmx/trjconv"  -f $xtc_file -s $struct_file -pbc nojump -skip $sample_interval  -o  $out_path # printf "0" selects group 0 (the whole system) to output.
