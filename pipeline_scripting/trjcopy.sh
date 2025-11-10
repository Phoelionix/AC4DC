set -u

gmx='/home/speno/programs/bin' # path to where GROMACS is installed

xtc_file=$1
out_folder=$2

handle=$(basename $working_folder)

#gro_file_handle=lys_example

i=1
out_path="$out_folder/output_$i.xtc"
while [ -f $out_path ]; do
    i=$((i+1))
    out_path="$out_folder/output_$i.xtc"
    if [ "$i" -ge 9999 ]; then
        echo "Critical error: Couldn't find a filename that doesn't exist!"
        exit
    fi 
done

mv $xtc_file $out_path